package fileIO;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;
import stream.bam.BgzfSettings;
import structures.IntList;
import structures.ListNum;


/**
 * Single-consumer byte-line reader with a reusable input buffer and optional
 * cached newline scanning selected by Shared.SIMD. Opens input in the constructor.
 * Nonempty returned lines are copies; empty lines share a zero-length array.
 * LF is removed, as is one preceding CR; a final unterminated line retains its CR.
 * No text decoding is performed.
 *
 * <p>Keep Shared.SIMD stable while consuming a reader. Direct line/batch reads,
 * pushback and reset are not mutually synchronized; synchronized open/close do
 * not make concurrent reading safe. The inherited nextList method serializes
 * its own callers and supplies record origins. Scalar batches contain at most
 * 200 lines; SIMD batches use the inherited line/byte targets.
 *
 * <p>Always close and check the returned error flag, including after EOF. The
 * SIMD batch path can return null without closing the stream. Read errors latch
 * an error flag rather than guaranteeing complete data or immediate failure;
 * bytes from successful reads remain available. Reset reopens the same source,
 * clears byte/newline/batch state, and preserves that error flag. It cannot
 * promise rewinding a non-reopenable source such as standard input.
 *
 * @author Brian Bushnell, Shinobu
 */
public final class ByteFile1 extends ByteFile{


	/** Diagnostic range printer or throughput loop, followed by close/reset/close.
	 * Positional options are filename (default stdin), first line or speedtest,
	 * exclusive last line or simd, and optional native BGZF. Throws after cleanup
	 * if reading or closing reported an error. The simd option enables cached SIMD
	 * scanning only when Vector.simd256 reports support. This is not an input validator. */
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
				//FIO-009: retain scalar fallback when the optional Vector API is unavailable.
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
		ByteFile1 tf=new ByteFile1(args.length>0 ? args[0] : "stdin", true);
		speedtest(tf, first, last, !speedtest);

		//FIO-011: retain errors across the diagnostic reopen and report a failing exit.
		boolean error=tf.close();
		tf.reset();
		error|=tf.close();
		if(error){throw new RuntimeException("ByteFile1: error reading or closing "+tf.name());}
	}

	/** Skips first lines and processes the half-open range [first,last), or until EOF.
	 * Timing includes skipping; counters exclude skipped lines and estimate bytes as
	 * returned line length plus one, not original file bytes. Reprint uses the default
	 * charset and println terminators; it does not preserve original line endings.
	 * @param bf Open reader with a stable SIMD setting
	 * @param first Number of initial lines to skip
	 * @param last Exclusive line limit
	 * @param reprint Whether to print selected lines instead of only timing statistics */
	private static void speedtest(ByteFile1 bf, long first, long last, boolean reprint){
		Timer t=new Timer();
		long lines=0;
		long bytes=0;
		for(long i=0; i<first; i++){bf.nextLine();}
		if(reprint){
			for(long i=first; i<last; i++){
				byte[] s=bf.nextLine();
				if(s==null){break;}

				lines++;
				bytes+=s.length+1;
				System.out.println(new String(s));
			}

			System.err.println("\n");
			System.err.println("Lines: "+lines);
			System.err.println("Bytes: "+bytes);
		}else{
			boolean nextList=false;
			if(nextList){
				for(ListNum<byte[]> ln=bf.nextList(); ln!=null; ln=bf.nextList()){
					for(byte[] line : ln.list){
						lines++;
						bytes+=line.length+1;
					}
				}
			}else if(!Shared.SIMD){
				for(long i=first; i<last; i++){
					byte[] s=bf.nextLine();
					if(s==null){break;}
					lines++;
					bytes+=s.length+1;
				}
			}else{
				for(long i=first; i<last; i++){
					byte[] s=bf.nextLineList();
					if(s==null){break;}
					lines++;
					bytes+=s.length+1;
				}
			}
		}
		t.stop();

		if(!reprint){
			System.err.println(Tools.timeLinesBytesProcessed(t, lines, bytes, 8));
			System.err.println("Bytes: "+bytes);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Opens the named input with format detection and a TEXT fallback.
	 * @param fname Non-null input name, including ReadWrite's supported special names
	 * @param allowSubprocess_ Whether input may use a decoding subprocess */
	public ByteFile1(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess_, false));
	}

	/** Opens input described by ff; the initial buffer uses the current bufferlen.
	 * @param ff Non-null input description */
	public ByteFile1(FileFormat ff){
		super(ff);
		if(verbose){System.err.println("ByteFile1("+ff+")");}
		is=open();
	}

	/*--------------------------------------------------------------*/
	/*----------------          Lifecycle           ----------------*/
	/*--------------------------------------------------------------*/

	/** Closes and reopens the same named input, discarding buffered input and clearing
	 * newline positions, list IDs and record origins. Retains the accumulated error flag.
	 * Reopening a non-seekable source does not guarantee a return to its beginning. */
	@Override
	public final void reset(){
		close();
		//FIO-007: open() resets cached SIMD positions along with byte offsets; stale positions
		//after a partial read previously made the reopened input appear exhausted.
		is=open();
		superReset();
	}

	/** Finishes the underlying input and any associated decoder, clears the stream
	 * reference and resets lineNum to -1. Repeated calls return the retained error flag.
	 * Does not coordinate with direct reads in progress.
	 * @return True if any read or close error has been recorded on this instance */
	@Override
	public synchronized final boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		if(!open){return errorState;}
		open=false;
		assert(is!=null);
		//		assert(false) : name()+","+allowSubprocess();
		errorState|=ReadWrite.finishReading(is, name(), ff.subprocess);

		is=null;
		lineNum=-1;
		//		pushBack=null;
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		return errorState;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Line reads          ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads one line, selecting the cached-position path when Shared.SIMD is true.
	 * Removes LF and one CR immediately before it; retains an unterminated final CR.
	 * An empty EOF closes the stream. A final nonempty line is returned first.
	 * @return Copied line bytes, shared empty array, or null for exhausted/closed input */
	@Override
	public final byte[] nextLine(){
		if(Shared.SIMD){return nextLineList();}
		if(verbose){System.err.println("Reading line "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}

		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}

		int nlpos=bstart;

		while(nlpos<bstop && buffer[nlpos]!=slashn){nlpos++;}

		if(nlpos>=bstop){//At this point we are at the last character which may or may not be a newline
			nlpos=fillBuffer();
		}

		// Check if buffer is empty after fill attempt
		if(bstop<1){
			close();
			return null;
		}

		lineNum++;

		// Determine the limit for copying
		// If nlpos >= bstop, we have data but no newline (EOF case) - use bstop
		// Otherwise, nlpos points to the newline
		final int limit;
		if(nlpos>=bstop){
			// No newline found, but we have data - use everything remaining
			limit=bstop;
		}else{
			// Found newline - exclude it (and any preceding \r)
			limit=(nlpos>bstart && buffer[nlpos-1]==slashr) ? nlpos-1 : nlpos;
		}

		if(bstart==limit){//Empty line.
			bstart=(nlpos<bstop) ? nlpos+1 : bstop;
			return blankLine;
		}
		
		byte[] line=KillSwitch.copyOfRange(buffer, bstart, limit);

		assert(line.length>0) : bstart+", "+nlpos+", "+limit;

		// Advance bstart past the newline (if there was one)
		bstart=(nlpos<bstop) ? nlpos+1 : bstop;

		return line;
	}

	/** Reads one line using cached newline positions. Requires Shared.SIMD to remain
	 * true so refills populate those positions; normal callers use nextLine instead.
	 * Uses the same newline stripping and empty-EOF closure as the scalar line path.
	 * @return Copied line bytes, shared empty array, or null for exhausted/closed input */
	public final byte[] nextLineList(){
		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}

		if(listPos>=positions.size()){
			fillBuffer(); // Now handles Vector.findSymbols internally
		}

		final int nlpos;
		if(listPos>=positions.size()){//No newlines
			nlpos=bstop;
		}else{
			nlpos=positions.get(listPos++);
		}

		// Check if buffer is empty after fill attempt
		if(bstop<1){
			close();
			return null;
		}

		lineNum++;

		// Determine the limit for copying
		final int limit;
		if(nlpos>=bstop){
			limit=bstop;
		}else{
			limit=(nlpos>bstart && buffer[nlpos-1]==slashr) ? nlpos-1 : nlpos;
		}

		if(bstart==limit){//Empty line
			bstart=(nlpos<bstop) ? nlpos+1 : bstop;
			return blankLine;
		}
		
		final byte[] line;
//		if(bstart+1==limit && buffer[bstart]==plus) {//Slows it down
//			line=plusLine;//Not really safe from modification...
//			assert(plusLine[0]==plus);
//		}else {
			line=KillSwitch.copyOfRange(buffer, bstart, limit);
//		}

		assert(line.length>0) : bstart+", "+nlpos+", "+limit;

		bstart=(nlpos<bstop) ? nlpos+1 : bstop;

		return line;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Batch reads         ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds a numbered batch of at most 200 lines, ignoring inherited target limits.
	 * nextLine performs counting and empty-EOF closure. Does not assign record origin.
	 * @return Nonempty numbered batch, or null when no line is available */
	private final ListNum<byte[]> nextListScalar(){
		final int listSize=200;
		final ArrayList<byte[]> list=new ArrayList<byte[]>(listSize);
		for(int i=0; i<listSize; i++){
			byte[] line=nextLine();
			if(line==null){
				close();
				break;
			}
			list.add(line);
			//FIXED [fileIO/ByteFile1#001]: removed a redundant 'lineNum++' here. nextLine() (called above) already
			//increments lineNum once per line (historical L171), so this second increment made lineNum() report 2x the real
			//line count on the scalar nextList path (Shared.SIMD off). The SIMD nextList and nextLineList paths
			//correctly increment once; this was the lone outlier.
		}
		return list.isEmpty() ? null : new ListNum<byte[]>(list, nextID++);
	}

	/** Builds a numbered batch without locking or assigning firstRecordNum.
	 * Scalar mode delegates to the fixed 200-line path. SIMD mode uses positive
	 * TARGET_LIST_SIZE and TARGET_LIST_BYTES, counting stripped line bytes and
	 * allowing a whole line to exceed the byte target. Null on SIMD EOF need not
	 * close the stream; the caller must close and inspect the error flag.
	 * @return Nonempty numbered batch, or null when no line is available */
	@Override
	public final ListNum<byte[]> nextListAsync(){
		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}
		if(!Shared.SIMD){return nextListScalar();}

		final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
		final ArrayList<byte[]> list=new ArrayList<byte[]>(slimit);
		int bytes=0;
		while(list.size()<slimit && bytes<blimit && open){
			if(listPos>=positions.size()){
				fillBuffer();

				// If still no newlines after fill, emit any trailing line (a final line
				// lacking a terminating newline) and finish; matches nextLineList, which
				// synthesizes nlpos=bstop for this case. Without this the trailing line was
				// dropped on the SIMD path, and the leftover bstart<bstop tripped
				// fillBuffer's assert(bstart>0) on the next call.
				if(positions.size()==0){
					if(bstop>bstart){//Trailing data with no final newline; mirror nextLineList's nlpos>=bstop branch (limit=bstop, no \r strip)
						lineNum++;
						list.add(KillSwitch.copyOfRange(buffer, bstart, bstop));
						bstart=bstop;//Consume it so the next fillBuffer sees bstart==bstop and cleanly EOFs
					}
					break;
				}
			}

			// Process positions until we have enough lines OR run out of positions
			final int iters=Math.min(positions.size()-listPos, slimit-list.size());
			for(int i=0; i<iters && bytes<blimit; i++){
				final int nlpos=positions.get(listPos++);

				// Check if buffer is empty
				if(bstop<1){
					close();
					break;
				}
				lineNum++;

				// Determine the limit for copying
				final int limit;
				if(nlpos>=bstop){
					limit=bstop;
				}else{
					limit=(nlpos>bstart && buffer[nlpos-1]==slashr) ? nlpos-1 : nlpos;
				}

				if(bstart==limit){//Empty line
					bstart=(nlpos<bstop) ? nlpos+1 : bstop;
					list.add(blankLine);
				}else{
					byte[] line=KillSwitch.copyOfRange(buffer, bstart, limit);
					assert(line.length>0) : bstart+", "+nlpos+", "+limit;
					bstart=(nlpos<bstop) ? nlpos+1 : bstop;
					list.add(line);
					bytes+=line.length;
				}
			}
		}

		return list.isEmpty() ? null : new ListNum<byte[]>(list, nextID++);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Buffer helpers       ----------------*/
	/*--------------------------------------------------------------*/

	/** Prints the populated buffer region for diagnostics, making CR/LF visible. */
	private final void printBuffer(){
		for(int i=0; i<bstop; i++){
			char c=(char)buffer[i];
			if(c=='\n'){
				System.err.println("\\n");
			}else if(c==slashr){
				System.err.print("\\r");
			}else{
				System.err.print(c);
			}
		}
	}

	/** Compacts a newline-free unread suffix and reads until a newline, EOF or caught error.
	 * Grows the array for long lines up to Integer.MAX_VALUE-8, then throws. Rebuilds
	 * newline positions when SIMD is enabled. Read exceptions latch errorState;
	 * only successful read calls advance the populated byte range.
	 * @return First newline offset, or bstop if no newline was found */
	private int fillBuffer(){
		if(bstart<bstop){ //Shift end bytes to beginning
			assert(bstart>0);
			int extra=bstop-bstart;
			for(int i=0; i<extra; i++, bstart++){
				buffer[i]=buffer[bstart];
				assert(buffer[i]!=slashn);
			}
			bstop=extra;
		}else{
			bstop=0;
		}

		bstart=0;
		
		// Clear SIMD state if applicable
		if(Shared.SIMD){
			listPos=0;
			positions.clear();
		}
		
		int len=bstop;
		int r=-1;
		while(len==bstop){//hit end of input without encountering a newline
			if(bstop==buffer.length){
				//Use long arithmetic above 1GiB; a full capped buffer must not become a zero-byte read/false EOF.
				final int maxCap=Integer.MAX_VALUE-8;
				if(buffer.length>=maxCap){
					throw new RuntimeException("ByteFile1: line in "+name()+" exceeds reader capacity ("+maxCap+" bytes)");
				}
				buffer=KillSwitch.copyOf(buffer, (int)Math.min(2L*buffer.length, maxCap));
			}
			//FIO-010: a caught exception must not reuse a prior positive read count
			//and advance bstop over bytes that were never successfully read.
			r=-1;
			try{
				r=is.read(buffer, bstop, buffer.length-bstop);
			}catch(IOException e){
				if(!errorState){//Warn loudly ONCE per file; errorState doubles as the spam-guard (decoupled from Shared.anomaly)
					System.err.println("Error: truncated or corrupt input for '"+name()+"'; data may be incomplete.");
					if(verbose){e.printStackTrace(); System.err.println("open="+open);}
				}
				errorState=true;//Best-effort: still return what was read, but latch the error so close() reports it (->nonzero exit)
			}catch(NullPointerException e){
				if(!errorState){//Same latch+warn-once for a null stream mid-read
					System.err.println("Error: read failure (null stream) for '"+name()+"'; data may be incomplete.");
					if(verbose){e.printStackTrace(); System.err.println("open="+open);}
				}
				errorState=true;
			}
			if(r>0){
				bstop=bstop+r;
				while(len<bstop && buffer[len]!=slashn){len++;}
			}else{
				len=bstop;
				break;
			}
		}

		assert(r==-1 || buffer[len]==slashn);
		
		// Find all newlines with SIMD if enabled
		if(Shared.SIMD){
			Vector.findSymbols(buffer, 0, bstop, slashn, positions);
		}
		
		return len;
	}

	/**
	 * Prepends a copy of the supplied bytes and an LF to unread buffered data.
	 * Grows/shifts storage when needed and refreshes cached newline positions.
	 * The next read applies normal LF/CRLF parsing to those bytes. Does not decrement
	 * counters or reopen closed input; repeated returned lines count again.
	 * @param line Non-null line bytes to prepend, normally without a line terminator
	 */
	@Override
	public void pushBack(byte[] line){
		if(bstart>line.length){
			bstart--;
			buffer[bstart]='\n';
			for(int i=0, j=bstart-line.length; i<line.length; i++, j++){
				buffer[j]=line[i];
			}
			bstart=bstart-line.length;
			refreshLinePositions();
			return;
		}

		int bLen=bstop-bstart;
		final long required=(long)bLen+line.length+1;
		final int maxCap=Integer.MAX_VALUE-8;
		if(required>maxCap){
			throw new RuntimeException("ByteFile1: pushback in "+name()+" requires "+required+" bytes; capacity is "+maxCap);
		}
		final int newLen=(int)required;
		int rShift=line.length+1-bstart;
		assert(rShift>0) : bstop+", "+bstart+", "+line.length;
		while(newLen>buffer.length){
			//This could get big if pushback is used often,
			//unless special steps are taken to prevent it, like leaving extra space for pushbacks.
			buffer=Arrays.copyOf(buffer, (int)Math.min(2L*buffer.length, maxCap));
		}

		Tools.shiftRight(buffer, rShift);

		for(int i=0; i<line.length; i++){
			buffer[i]=line[i];
		}
		buffer[line.length]='\n';
		bstart=0;
		bstop=newLen;
		refreshLinePositions();
	}

	/** Pushback changes offsets and adds a newline; cached SIMD positions must describe the unread buffer. */
	private void refreshLinePositions(){
		assert(bstart>=0 && bstart<=bstop && bstop<=buffer.length) : "Invalid unread range after pushback: "+bstart+","+bstop+","+buffer.length;
		listPos=0;
		positions.clear();
		if(Shared.SIMD){Vector.findSymbols(buffer, bstart, bstop, slashn, positions);}
	}

	/**
	 * Opens the named input and resets byte offsets and cached newline positions.
	 * Synchronizing this helper does not protect the unsynchronized read methods.
	 * Does not clear the error flag or reset the base class's batch counters.
	 * @return Opened InputStream for the file
	 * @throws RuntimeException if file is already open
	 */
	private final synchronized InputStream open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened TextFile "+name());
		}
		open=true;
		is=ReadWrite.getInputStream(name(), BUFFERED, allowSubprocess(), true);
		bstart=-1;
		bstop=-1;
		listPos=0;
		positions.clear();
		return is;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Accessors           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the open flag; SIMD batch EOF alone may leave it true. */
	@Override
	public boolean isOpen(){return open;}

	/** Borrows the underlying stream, or null after close. Reading it directly bypasses
	 * already buffered bytes and invalidates normal reader sequencing. */
	@Override
	public final InputStream is(){return is;}

	/** Returns the zero-based returned-line counter, or -1 before reading/after close.
	 * Pushback does not rewind it; it is not a physical file offset or final total. */
	@Override
	public final long lineNum(){return lineNum;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Input lifecycle flag; not a lock or a complete EOF predicate. */
	private boolean open=false;
	/** Owned input storage, expanded when a line exceeds its current capacity. */
	private byte[] buffer=new byte[bufferlen];
	/** Shared zero-length array for empty lines. */
	private static final byte[] blankLine=new byte[0];
	/** Unread half-open byte interval; both start at -1 after opening. */
	private int bstart=0, bstop=0;
	/** Underlying stream owned by this reader; callers must not replace it. */
	public InputStream is;
	/** Returned-line counter; direct external assignment changes lineNum(). */
	public long lineNum=-1;
	/** Newline offsets in the current buffer when SIMD mode is active. */
	private IntList positions=new IntList();
	/** Index of the next unused cached newline offset. */
	private int listPos=0;

	/** Diagnostic logging switch. Configure before using readers. */
	public static boolean verbose=false;
	/** Whether ReadWrite should add stream buffering on open; this reader buffers either way. */
	public static boolean BUFFERED=false;
	/** Positive initial byte-array capacity for newly constructed readers. */
	public static int bufferlen=65536;

	/** Read/close error latch retained across reset and reported by close(). */
	private boolean errorState=false;


}
