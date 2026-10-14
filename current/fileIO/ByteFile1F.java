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
import structures.IntList;
import structures.ListNum;

/**
 * Buffered FASTA record reader using LF-before-header boundaries instead of lines.
 * Returns copied byte blocks with sequence newlines retained; no sequence parsing or
 * FASTA validation is performed. Each block retains its leading header marker;
 * the final block also retains trailing bytes.
 * Boundary detection uses the shared SIMD setting with a scalar fallback.
 * This class is directly constructible and is not selected by ByteFile factories.
 * Reading, batching and reset require a single owner; synchronized open/close alone
 * do not make concurrent reading or reconfiguration safe.
 * 
 * @author Brian Bushnell
 * @contributor Shinobu (record-reader contracts)
 * @date November 11, 2025
 */
public final class ByteFile1F extends ByteFile {

	/*--------------------------------------------------------------*/
	/*----------------           Startup            ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs a positional record-printing or speed diagnostic, then close/reset/close.
	 * Parses arguments before opening input. After construction, exceptional exits
	 * close the reader and preserve the primary failure. Ordinary completion checks
	 * accumulated input status and, in print mode, borrowed stdout.
	 * Optional arguments select a record range, capability-gated SIMD and native BGZF.
	 * Printed byte totals add one byte per returned block, not the exact input-file
	 * byte count. This is not a format validator. */
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
				//FIO-074: keep scalar scanning when the optional Vector API is unavailable.
				if(Vector.simd256){Shared.SIMD=true;}
			}else {
				last=Integer.parseInt(args[2]);
			}
		}
		if(args.length>3){
			if(args[3].equalsIgnoreCase("native")){
				ReadWrite.ALLOW_NATIVE_BGZF=BgzfSettings.USE_MULTITHREADED_BGZF=true;
			}
		}
		ByteFile1F tf=new ByteFile1F(args.length>0 ? args[0] : "stdin", true);
		boolean error=false;
		Throwable failure=null;
		try{
			speedtest(tf, first, last, !speedtest);
			// Retain the successful diagnostic reopen without losing earlier errors.
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
				failure.addSuppressed(new RuntimeException("ByteFile1F: error reading or closing "+tf.name()));
			}
		}
		if(error){throw new RuntimeException("ByteFile1F: error reading or closing "+tf.name());}
		if(!speedtest && System.out.checkError()){
			throw new RuntimeException("ByteFile1F: error writing output");
		}
	}

	/** Skips to first, then prints or counts blocks in [first,last), stopping at EOF. */
	private static void speedtest(ByteFile1F bf, long first, long last, boolean reprint){
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
				System.out.println(new String(s));
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

	/** Opens a named input with FASTA as the fallback format.
	 * @param fname Input filename or supported stream name
	 * @param allowSubprocess_ Whether decompression subprocesses are permitted */
	public ByteFile1F(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTA, null, allowSubprocess_, false));
	}

	/** Borrows the input descriptor and opens its stream using current static settings.
	 * @param ff Descriptor configured for reading */
	public ByteFile1F(FileFormat ff){
		super(ff);
		if(verbose){System.err.println("ByteFile1F("+ff+")");}
		is=open();
	}

	/*--------------------------------------------------------------*/
	/*----------------        Reader Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Closes and reopens the input, clearing buffered positions and inherited batch
	 * counters. Prior close errors remain latched. Reopening a borrowed stream does
	 * not itself rewind that stream. */
	@Override
	public final void reset(){
		close();
		is=open();
		superReset();
		firstRecord=true;
		listPos=0;
		positions.clear();
	}

	/** Finishes the input once, releases this stream reference and resets lineNum.
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

	/** Returns the next copied FASTA block, or null after exhaustion/closure.
	 * Searches for an initial header marker, then splits at LF-before-header positions.
	 * Header markers and internal newlines remain. A caught input exception is
	 * reported through close status; buffered prefix bytes may still be returned.
	 * @return Copied record bytes, the shared empty array on the empty-boundary path,
	 * or null when no bytes remain */
	@Override
	public final byte[] nextLine(){
		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}

		// For the first record, find the first '>'
		if(firstRecord){
			if(bstop==0){fillBuffer();}
			while(bstart<bstop && buffer[bstart]!=carrot){bstart++;}
			if(bstart>=bstop){
				close();
				return null;
			}
			firstRecord=false;
		}

		// Need more positions?
		if(listPos>=positions.size()){
			fillBuffer();
			if(listPos>=positions.size() && bstart>=bstop){
				// EOF
				close();
				return null;
			}
		}

		lineNum++;

		// Determine record boundaries
		final int limit;
		if(listPos<positions.size()){
			// Next boundary is at positions[listPos] (points to \n before >)
			limit=positions.get(listPos++);
		}else{
			// No more boundaries, use rest of buffer
			limit=bstop;
		}

		if(bstart>=limit){
			// Empty record or at end
			if(bstart<bstop){
				bstart++; // Skip the >
			}
			return blankLine;
		}

		byte[] record=KillSwitch.copyOfRange(buffer, bstart, limit);

		// Move past the delimiter LF, preserving the next header marker.
		if(limit<bstop){
			bstart=limit+1; //FIO-072: the boundary points to LF, not to the next record start.
		}else{
			bstart=bstop;
		}

		return record;
	}
	
	/** Collects up to 200 blocks synchronously and assigns the next batch ID.
	 * Creates no worker thread. Direct calls require the reader's single owner.
	 * @return A nonempty batch, or null when no block is available */
	@Override
	public final ListNum<byte[]> nextListAsync(){
		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}
		
		final int listSize=200;
		final ArrayList<byte[]> list=new ArrayList<byte[]>(listSize);
		for(int i=0; i<listSize; i++){
			byte[] record=nextLine();
			if(record==null){
				break;
			}
			list.add(record);
		}
		return list.isEmpty() ? null : new ListNum<byte[]>(list, nextID++);
	}

	/** Compacts unread bytes, then reads until a boundary or a nonpositive read result.
	 * Boundary positions refer to the LF byte. Grows the buffer with capped long
	 * arithmetic and rejects a full maximum-capacity buffer instead of treating it
	 * as EOF. This does not validate headers or distinguish a zero-byte read from EOF. */
	private void fillBuffer(){
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
		while(positions.isEmpty()){
			if(bstop==buffer.length){
				//Use long arithmetic above 1GiB; never treat a full capped record buffer as EOF.
				final int maxCap=Integer.MAX_VALUE-8;
				if(buffer.length>=maxCap){
					throw new RuntimeException("ByteFile1F: FASTA record in "+name()+" exceeds reader capacity ("+maxCap+" bytes)");
				}
				buffer=KillSwitch.copyOf(buffer, (int)Math.min(2L*buffer.length, maxCap));
			}
			
			int r=-1;
			try{
				r=is.read(buffer, bstop, buffer.length-bstop);
			}catch(IOException e){
				errorState=true; //FIO-073: preserve the read failure even if closing succeeds.
				if(!Shared.anomaly){e.printStackTrace();}
			}catch(NullPointerException e){
				errorState=true;
				if(!Shared.anomaly){e.printStackTrace();}
			}
			
			if(r<=0){break;} // EOF
			
			final int from=Math.max(0, bstop-1);
			bstop+=r;

			// Find all \n> boundaries in the buffer
			Vector.findFastaHeaders(buffer, from, bstop, positions);
		}
	}

	/** Rejects pushback; this record reader does not implement it.
	 * @param record Unused record
	 * @throws UnsupportedOperationException Always */
	@Override
	public void pushBack(byte[] record){
		throw new UnsupportedOperationException("pushBack not supported for ByteFile1F");
	}

	/** Opens the configured input and resets buffer/boundary state, preserving errors.
	 * @return The stream also stored in is
	 * @throws RuntimeException If already open, or if the input factory fails */
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

	/** Reports the local open flag; not a concurrent-completion guarantee. */
	@Override
	public boolean isOpen(){return open;}

	/** Exposes the current stream without transferring ownership; null after close. */
	@Override
	public final InputStream is(){return is;}

	/** Returns the zero-based block counter, or -1 before reading/after close. */
	@Override
	public final long lineNum(){return lineNum;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Local lifecycle flag, not volatile. */
	private boolean open=false;
	/** Reused input storage; each returned nonempty block is copied out. */
	private byte[] buffer=new byte[bufferlen];
	/** Shared result for an empty boundary interval. */
	private static final byte[] blankLine=new byte[0];
	/** Unread interval within buffer. */
	private int bstart=0, bstop=0;
	/** Current input reference; finishReading handles its closure or borrowed-stream policy. */
	public InputStream is;
	/** Block counter, incremented before extracting each interval. */
	public long lineNum=-1;
	/** LF indices before header markers in the current buffer. */
	private IntList positions=new IntList();
	/** Next unconsumed boundary index. */
	private int listPos=0;
	/** Whether the initial header search remains to be performed. */
	private boolean firstRecord=true;

	/** Shared diagnostic switch; configure before reading. */
	public static boolean verbose=false;
	/** Input factory buffering request, sampled when opening. */
	public static boolean BUFFERED=false;
	/** Initial buffer capacity for newly constructed readers. */
	public static int bufferlen=65536;

	/** Sticky caught-read or finishReading failure status; reset does not clear it. */
	private boolean errorState=false;
}
