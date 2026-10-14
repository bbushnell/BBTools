package fileIO;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;
import stream.bam.BgzfSettings;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;

/**
 * Buffered, condensed FASTA block reader used by the FastaStreamer2 readers.
 * The primary API, nextLine(IntList), returns a newly allocated backing array and
 * records the meaningful header/sequence newline positions in the supplied list.
 * Physical sequence lines are joined and a terminal CR is removed from each line;
 * the backing-array capacity can exceed the meaningful output length.
 * The no-argument compatibility API trims that array but discards the boundaries.
 * This is not a general FASTA validator.
 * <p>Reading and reset require a single owner. Synchronized open/close do not make
 * concurrent consumption or reconfiguration safe. Configure shared SIMD and input
 * settings before opening a reader; scalar boundary scanning is also supported.
 *
 * @author Brian Bushnell
 * @contributor Shinobu (condensed block and lifecycle contracts)
 * @date November 11, 2025
 */
public final class ByteFile1Fc extends ByteFile{

	/*--------------------------------------------------------------*/
	/*----------------           Startup            ----------------*/
	/*--------------------------------------------------------------*/

	/** Prints or counts a positional block range, then closes/reopens/closes.
	 * Parses arguments before opening input. After construction, thrown processing
	 * closes input and preserves the primary failure. Normal completion checks sticky
	 * input status and, in print mode, borrowed stdout without closing it.
	 * Optional arguments enable supported SIMD or native BGZF through shared settings.
	 * Counts describe blocks and add one byte per block, not exact printed bytes. */
	public static void main(String[] args){
		long first=0, last=100;
		boolean speedtest=false;
		if(args.length>1){
			if(args[1].equalsIgnoreCase("speedtest")){
				speedtest=true;
				first=0;
				last=Long.MAX_VALUE;
			}else{
				first=Integer.parseInt(args[1]);
				last=first+100;
			}
		}
		if(args.length>2){
			if(args[2].equalsIgnoreCase("simd")){
				//[FIO-079] The optional request must respect Vector availability and supported width.
				if(Vector.simd256){Shared.SIMD=true;}
			}else{
				last=Integer.parseInt(args[2]);
			}
		}
		if(args.length>3){
			if(args[3].equalsIgnoreCase("native")){
				ReadWrite.ALLOW_NATIVE_BGZF=BgzfSettings.USE_MULTITHREADED_BGZF=true;
			}
		}
		ByteFile1Fc tf=new ByteFile1Fc(args.length>0 ? args[0] : "stdin", true);
		boolean error=false;
		Throwable failure=null;
		try{
			speedtest(tf, first, last, !speedtest);
			//[FIO-080] Keep the normal reopen diagnostic without discarding earlier errors.
			error|=tf.close();
			tf.reset();
		}catch(RuntimeException | Error t){
			failure=t;
			throw t;
		}finally{
			try{
				error|=tf.close();
			}catch(RuntimeException | Error t){
				if(failure==null){throw t;}
				if(failure!=t){failure.addSuppressed(t);}
			}
			if(failure!=null && error){
				failure.addSuppressed(new RuntimeException("ByteFile1Fc: error reading or closing "+tf.name()));
			}
		}
		if(error){throw new RuntimeException("ByteFile1Fc: error reading or closing "+tf.name());}
		if(!speedtest && System.out.checkError()){
			throw new RuntimeException("ByteFile1Fc: error writing output");
		}
	}

	/** Skips initial blocks, then prints or counts the requested range until EOF.
	 * Byte totals add one per trimmed block, rather than measuring input or printed bytes. */
	private static void speedtest(ByteFile1Fc bf, long first, long last, boolean reprint){
		Timer t=new Timer();
		long records=0;
		long bytes=0;
		for(long i=0; i<first; i++){bf.nextLine();}
		if(reprint){
			for(long i=first; i<last; i++){
				byte[] s=bf.nextLine();
				if(s==null){break;}
				records++;
				bytes+=s.length+1;
				System.out.print("*");
				System.out.print(new String(s));
				System.out.println("*");
			}
			System.err.println("\n");
			System.err.println("Records: "+records);
			System.err.println("Bytes: "+bytes);
		}else{
			for(long i=first; i<last; i++){
				byte[] s=bf.nextLine();
				if(s==null){break;}
				records++;
				bytes+=s.length+1;
			}
		}
		t.stop();

		if(!reprint){
			System.err.println(Tools.timeLinesBytesProcessed(t, records, bytes, 8));
			System.err.println("Bytes: "+bytes);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens a named input, using FASTA as the fallback format.
	 * @param fname Input filename or supported stream name
	 * @param allowSubprocess_ Whether the input factory may use subprocesses */
	public ByteFile1Fc(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTA, null, allowSubprocess_, false));
	}

	/** Borrows an input descriptor and opens its stream with current static settings.
	 * @param ff Descriptor configured for reading */
	public ByteFile1Fc(FileFormat ff){
		super(ff);
		if(verbose){System.err.println("ByteFile1F("+ff+")");}
		is=open();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Reader Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Closes/reopens input and resets buffer positions and inherited batch counters.
	 * Sticky errors survive reset. Reopening a borrowed stream does not rewind it. */
	@Override
	public final void reset(){
		close();
		is=open();
		superReset();
		firstRecord=true;
		listPos=0;
		positions.clear();
	}

	/** Finishes an open input once; normal return clears its reference and lineNum.
	 * @return Sticky caught-read or finishReading error status */
	@Override
	public synchronized final boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		if(!open){return errorState;}
		open=false;
		assert(is!=null);
		errorState|=ReadWrite.finishReading(is, name(), ff.subprocess);

		is=null;
		lineNum=-1;
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		return errorState;
	}

	/** Returns a length-correct block while discarding its record boundaries.
	 * May copy again to trim backing-array padding; block condensation already copies
	 * input bytes. Production Streamer2 readers use nextLine(IntList) instead.
	 * @return Trimmed block, or null after exhaustion/closure
	 * @deprecated Use {@link #nextLine(IntList)} to receive meaningful boundaries
	 * and avoid the additional trimming copy.
	 */
	@Override
	@Deprecated
	public final byte[] nextLine(){
		final IntList newlines=new IntList();
		final byte[] block=nextLine(newlines);
		if(block==null){return null;}
		//nextLine(IntList) returns the oversized ByteBuilder backing array; the final newline
		//marks the true end of content, so trim the trailing NUL padding for a correct length.
		final int len=newlines.size>0 ? newlines.get(newlines.size-1)+1 : 0;
		return len==block.length ? block : KillSwitch.copyOf(block, len);
	}
	
	/** Condenses the next complete-record prefix into a newly allocated backing array.
	 * Output newline indices delimit headers and joined sequences; the last index plus
	 * one gives the meaningful length. Array capacity is not a record-length oracle.
	 * Condensation clears/refills newlines; an already-closed reader returns null
	 * without changing the supplied list. EOF processing closes exhausted input.
	 * @param newlines Caller-owned list receiving output newline positions
	 * @return Backing array, possibly with unused capacity, or null when no data remains */
	public final byte[] nextLine(IntList newlines){
		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}
		int lastRecordLoc=fillBuffer();
		byte[] records=condense(buffer, positions, newlines, lastRecordLoc);
		return records;
	}
	
	/** Synchronously wraps one trimmed compatibility block in a numbered list.
	 * Creates no worker and exports no per-record boundary list.
	 * @return One-element batch, or null after exhaustion/closure */
	@Override
	public final ListNum<byte[]> nextListAsync(){
		return nextList(new IntList());
	}
	
	/** Wraps one trimmed compatibility block and advances the inherited batch ID.
	 * FIO-081 records that this overload currently ignores the supplied boundary list.
	 * @param newlines Currently unused; no boundary positions are exported
	 * @return One-element batch, or null when no block is available */
	public final ListNum<byte[]> nextList(IntList newlines){
		//TODO: Probable bug - FIO-081: this overload discards the supplied boundary list.
		//No production caller found; live Streamer2 readers use nextLine(IntList) directly.
		byte[] records=nextLine();
		if(records==null){return null;}
		ArrayList<byte[]> list=new ArrayList<byte[]>(1);
		list.add(records);
		return new ListNum<byte[]>(list, nextID++);
	}

//	/** Fill buffer and find all FASTA boundaries */
//	private int fillBuffer(){
//		// Shift remaining data to start
//		if(bstart>0 && bstart<bstop){
//			int extra=bstop-bstart;
//			System.arraycopy(buffer, bstart, buffer, 0, extra);
//			bstop=extra;
//			bstart=0;
//		}else if(bstart>=bstop){
//			bstart=0;
//			bstop=0;
//		}
//
//		// Clear position list
//		listPos=0;
//		positions.clear();
//
//		// Read data
//		int lastRecordLoc=-1;
//		while(lastRecordLoc<0){
//			if(bstop==buffer.length){
//				buffer=KillSwitch.copyOf(buffer, buffer.length*2);
//			}
//			
//			int r=0;
//			try{
//				r=is.read(buffer, bstop, buffer.length-bstop);
//			}catch(IOException e){
//				if(!Shared.anomaly){e.printStackTrace();}
//			}catch(NullPointerException e){
//				if(!Shared.anomaly){e.printStackTrace();}
//			}
//			
////			final int from=Math.max(0, bstop-1);//Not safe because buffer can come in with residual newlines
//			bstop+=r;
//			
//			if(firstRecord){
//			    assert(bstop<1 || buffer[0]==carrot) : "File does not start with '>' - "+ff.name()+" is not a valid FASTA file";
//			    firstRecord=false;
//			}
//
//			// Find all \n> boundaries in the buffer
//			positions.clear();
//			Vector.findSymbols(buffer, 0, bstop, slashn, positions);
//			for(int i=positions.size-1; i>=0 && lastRecordLoc<0; i--) {
//				int loc=positions.array[i];
//				if(loc+1<bstop && buffer[loc+1]==carrot) {lastRecordLoc=loc;}
//			}
//			
//			if(r<=0){// EOF
//				if(bstop>0 && buffer[bstop-1]!='\n') {positions.add(bstop);}
//				break;
//			}
//		}
//		return lastRecordLoc;
//	}
	
	/** Compacts residual input and scans until a record boundary or nonpositive read result.
	 * Keeps all physical newline positions, adding a synthetic end position for an
	 * unterminated tail. Growth uses capped long arithmetic. Both caught read exception types
	 * latch sticky status so subsequent EOF and close cannot hide the failure.
	 * @return LF index before the last detected next header, or -1 if none was found */
	private int fillBuffer(){
		// Shift remaining data to start
		if(bstart>0 && bstart<bstop){
			int extra=bstop-bstart;
			System.arraycopy(buffer, bstart, buffer, 0, extra);
			bstop=extra;
			bstart=0;
		}else if(bstart>=bstop){
			bstart=0;
			bstop=0;
		}

		// Clear position list
		listPos=0;
		positions.clear();

		// Read data
		int lastRecordLoc=-1;
		int scanStart=0;//First byte not yet scanned in this refill
		
		while(lastRecordLoc<0){
			if(bstop==buffer.length){
				//G11's reader fix: long arithmetic prevents overflow above 1GiB. Fail at the cap
				//instead of issuing read(...,0), which would silently truncate the record as false EOF.
				final int maxCap=Integer.MAX_VALUE-8;
				if(buffer.length>=maxCap){
					throw new RuntimeException("ByteFile1Fc: FASTA record in "+ff.name()+" exceeds reader capacity ("
						+maxCap+" bytes) with no '\\n>' boundary; records >~2GB require a streaming reader, not this whole-record buffer.");
				}
				buffer=KillSwitch.copyOf(buffer, (int)Math.min(2L*buffer.length, maxCap));
			}
			
			int r=0;
			try{
				r=is.read(buffer, bstop, buffer.length-bstop);
			}catch(IOException e){
				errorState=true;
				if(!Shared.anomaly){e.printStackTrace();}
			}catch(NullPointerException e){
				//[FIO-078] EOF processing and later close must preserve this caught read failure.
				errorState=true;
				if(!Shared.anomaly){e.printStackTrace();}
			}
			
//			final int from=Math.max(0, bstop-1);//Not safe because buffer can come in with residual newlines
			//EOF returns -1; only bytes actually read extend the buffer. Preserve a final base without newline.
			if(r>0){bstop+=r;}
			
			if(firstRecord){
				assert(bstop<1 || buffer[0]==carrot) : "File does not start with '>' - "+ff.name()+" is not a valid FASTA file";
				firstRecord=false;
			}

			// Find all \n> boundaries in the buffer
			// Scan only newly read bytes; previously found newline positions remain valid.
			Vector.findSymbols(buffer, scanStart, bstop, slashn, positions);
			scanStart=bstop; // Update scanStart for next iteration
			
			for(int i=positions.size-1; i>=0 && lastRecordLoc<0; i--){
				int loc=positions.array[i];
				if(loc+1<bstop && buffer[loc+1]==carrot){lastRecordLoc=loc;}
			}
			
			if(r<=0){// EOF
				if(bstop>0 && buffer[bstop-1]!='\n'){positions.add(bstop);}
				break;
			}
		}
		return lastRecordLoc;
	}
	
	/** Copies a complete prefix, joins sequence lines and records output boundaries.
	 * Trims one terminal CR per physical line and emits header/sequence newlines.
	 * Retains the record after the last detected boundary, or consumes all data
	 * when no next-header boundary was found. Empty input closes the reader.
	 * @param buffer Input storage, with meaningful bytes ending at bstop
	 * @param newlines1 Physical LF positions, including any synthetic EOF boundary
	 * @param newlines2 Replaced with emitted newline positions
	 * @param lastRecordLoc LF before the retained next header, or -1 at the final prefix
	 * @return Newly allocated builder backing array, or null when no prefix remains */
	private byte[] condense(byte[] buffer, IntList newlines1, IntList newlines2, int lastRecordLoc){
		//At EOF, use the final physical newline's index, or the synthetic boundary for an unterminated line.
		final int lastByte=lastRecordLoc>0 ? lastRecordLoc :
			(bstop>0 && buffer[bstop-1]==slashn ? bstop-1 : bstop);
		newlines2.clear();
		if(lastByte<1){
			close();
			return null;
		}
		ByteBuilder bb=new ByteBuilder(lastByte+1);
		int prev=-1;
		if(verbose){System.err.println("Entered condense: lastRecordLoc="+
			lastRecordLoc+", bstop="+bstop+", newlines="+newlines1);}
		for(int i=0; i<newlines1.size && prev<lastByte; i++){
			int loc=newlines1.get(i);
			int len=loc-prev-1;
			assert(len>=0) : "len="+len+", loc="+loc+", i="+i+", prev="+prev+
				", bb.length="+bb.length+", newlines1.size="+newlines1.size;
			//[FIO-082] Empty records add sequence newlines, so output can outgrow the physical prefix.
			bb.ensureExtra(len);
			System.arraycopy(buffer, prev+1, bb.array, bb.length, len);//Copy excluding newline
			if(verbose){System.err.println("Added "+len+" bytes.");}
			bb.length+=len;
			if(bb.endsWith('\r')){bb.length--;}//Trim \r
			if(buffer[prev+1]==carrot){//Record start
				newlines2.add(bb.length);
				bb.nl();
				if(verbose){System.err.println("Added header start newline.");}
			}
			if(loc>=lastByte || buffer[loc+1]==carrot){//Record end
				newlines2.add(bb.length);
				bb.nl();
				if(verbose){System.err.println("Added record stop newline.");}
			}
			prev=loc;
		}
		//Shift residual
		int residual=bstop-lastRecordLoc-1;
		if(lastRecordLoc>0 && residual>0){
			System.arraycopy(buffer, lastRecordLoc+1, buffer, 0, residual);
			bstop=residual;
		}else{
			bstop=0;
		}
		//assert(false) : "*"+bb.toString()+"*";
		return bb.array;
	}

	/** Rejects pushback for this block reader.
	 * @param record Unused argument
	 * @throws UnsupportedOperationException Always */
	@Override
	public void pushBack(byte[] record){
		throw new UnsupportedOperationException("pushBack not supported for ByteFile1F");
	}

	/** Opens configured input and clears buffered boundary state, preserving errors.
	 * @return The opened stream, also stored in is
	 * @throws RuntimeException If already open or if input creation fails */
	private final synchronized InputStream open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened ByteFile1F "+name());
		}
		open=true;
		is=ReadWrite.getInputStream(name(), BUFFERED, allowSubprocess(), true);
		bstart=0;
		bstop=0;
		firstRecord=true;
		listPos=0;
		positions.clear();
		return is;
	}

	/** Returns the local lifecycle flag; not a concurrent-completion guarantee. */
	@Override
	public boolean isOpen(){return open;}

	/** Exposes the current input without transferring ownership; null after normal close. */
	@Override
	public final InputStream is(){return is;}

	/** Returns the exposed counter field, which block consumption does not advance. */
	@Override
	public final long lineNum(){return lineNum;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Local open state, not volatile. */
	private boolean open=false;
	/** Reused input buffer; output blocks are copied into separate storage. */
	private byte[] buffer=new byte[bufferlen];
	/** Unused historical empty-result constant. */
	private static final byte[] blankLine=new byte[0];
	/** Buffered input bounds; the current condensation path keeps bstart at zero. */
	private int bstart=0, bstop=0;
	/** Current backend reference, closed according to finishReading ownership rules. */
	public InputStream is;
	/** Historical counter field; reset by close but not advanced by block reads. */
	public long lineNum=-1;
	/** Physical input LF indices and any synthetic EOF boundary. */
	private IntList positions=new IntList();
	/** Historical cursor, cleared during refill/reset but unused by current condensation. */
	private int listPos=0;
	/** Whether the initial input-header assertion is still required. */
	private boolean firstRecord=true;

	/** Compile-time diagnostic switch. */
	public static final boolean verbose=false;
	/** Input-factory buffering request, sampled when opening. */
	public static boolean BUFFERED=false;
	/** Initial input capacity, sampled by each new reader. */
	public static int bufferlen=262144;

	/** Sticky caught-read and input-finish error status. */
	private boolean errorState=false;
}
