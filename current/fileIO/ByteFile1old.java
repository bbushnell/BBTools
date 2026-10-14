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
 * Retained byte-line reader with scalar scanning and a cached newline-position path.
 * ByteFile's factory selects other implementations; this class is constructed directly.
 * LF terminates a line, and one CR immediately before LF is removed. An unterminated
 * final line is returned as data, including a trailing CR. Nonempty results are copies;
 * empty lines share a zero-length array.
 *
 * nextLine selects its path from Shared.SIMD. Batches contain at most 200 lines,
 * independently of ByteFile's target-size settings. The inherited nextList wrapper
 * assigns record origins; direct nextListAsync calls do not. Keep read-mode settings
 * stable during iteration because the two paths share a cursor but cache boundaries.
 *
 * Cursor, buffer and boundary-cache state belong to the instance. Direct iteration
 * and reset are not serialized by synchronized open/close. Reset reopens the input,
 * reuses the buffers and resets list numbering; accumulated read/closure errors remain.
 * Single-line EOF closes input, while SIMD batch EOF can return null with input open.
 * Caught read exceptions latch the status returned by close. The diagnostic checks this
 * status; cached pushBack behavior remains an unresolved finding at its implementation site.
 *
 * @author Brian Bushnell, Shinobu
 */
public final class ByteFile1old extends ByteFile {

	/** Runs the positional diagnostic and closes/reopens/closes input on normal return.
	 * Parses positions before opening input. After construction, processing failures
	 * trigger cleanup while preserving the primary exception. Normal completion checks
	 * accumulated input errors and, in print mode, borrowed stdout without closing it.
	 * @param args Input (default stdin), first line or speedtest, exclusive last line or
	 * simd, and optional native BGZF mode; numeric positions are parsed as ints
	 */
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
				// FIO-089: the optional request must respect Vector availability and supported width.
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
		ByteFile1old tf=new ByteFile1old(args.length>0 ? args[0] : "stdin", true);
		boolean error=false;
		Throwable failure=null;
		try{
			speedtest(tf, first, last, !speedtest);
			// FIO-088: retain the reopen diagnostic without losing earlier read/close errors.
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
				failure.addSuppressed(new RuntimeException("ByteFile1old: error reading or closing "+tf.name()));
			}
		}
		if(error){throw new RuntimeException("ByteFile1old: error reading or closing "+tf.name());}
		if(!speedtest && System.out.checkError()){
			throw new RuntimeException("ByteFile1old: error writing output");
		}
	}

	/** Skips a prefix, then prints lines or reports a timed count through the exclusive bound.
	 * Diagnostic byte counts add one byte per line regardless of the original terminator.
	 * The local nextList switch preserves an inactive alternative benchmark path.
	 * @param bf Reader whose current position is consumed
	 * @param first Number of lines to skip before processing
	 * @param last Exclusive upper line position
	 * @param reprint Whether to print using the platform default character encoding
	 */
	private static void speedtest(ByteFile1old bf, long first, long last, boolean reprint){
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

	/** Creates an input descriptor with TEXT as the fallback format and opens it.
	 * @param fname Input path or recognized special input name
	 * @param allowSubprocess_ Whether the descriptor permits subprocess decompression
	 */
	public ByteFile1old(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess_, false));
	}

	/** Borrows the input descriptor and opens its stream immediately.
	 * @param ff Nonnull descriptor configured for reading
	 */
	public ByteFile1old(FileFormat ff){
		super(ff);
		if(verbose){System.err.println("ByteFile1("+ff+")");}
		is=open();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Closes and reopens input, clearing cached boundaries and base-class list counters.
	 * Retains the allocated buffers and accumulated read/closure error flag.
	 */
	@Override
	public final void reset(){
		close();
		is=open();
		superReset();
	}

	/** Finishes input/subprocess cleanup and drops the stream reference.
	 * An open reader's local line counter becomes -1; repeated close retains prior status.
	 * @return Accumulated caught-read and finishReading error status
	 */
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

	/** Returns the next line, selecting cached newline scanning when Shared.SIMD is set.
	 * Removes LF and an immediately preceding CR, but preserves CR at unterminated EOF.
	 * @return Copied nonempty line, shared empty array, or null for exhausted/closed input
	 */
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
		// If nlpos>=bstop, we have data but no newline (EOF case) - use bstop
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

	/** Reads one line through cached newline offsets, refilling/scanning when exhausted.
	 * Vector.findSymbols chooses its own scalar/SIMD implementation. When a scan finds
	 * no LF before EOF, this method uses bstop as the remaining line's exclusive limit.
	 * @return Copied nonempty line, shared empty array, or null after closing exhausted input
	 */
	public final byte[] nextLineList(){

		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}

		if(listPos>=positions.size()){
			listPos=0;
			positions.clear();
			fillBuffer();
			Vector.findSymbols(buffer, 0, bstop, slashn, positions);
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
		// If nlpos>=bstop, we have data but no newline (EOF case) - use bstop
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
	
	/** Collects up to 200 nextLine results and allocates one successive batch ID.
	 * A short batch's EOF probe can close input and reset the local line counter.
	 * @return Nonempty batch, or null if the first read is exhausted
	 */
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
		}
		return list.isEmpty() ? null : new ListNum<byte[]>(list, nextID++);
	}

	/** Collects at most 200 lines with scalar nextLine calls or cached newline offsets.
	 * The cached path also emits a nonempty unterminated tail. Its empty EOF branch
	 * can leave input open; callers still close the reader. Assigns nextID, while only
	 * the inherited nextList wrapper assigns firstRecordNum. This method is not synchronized.
	 * @return Nonempty batch, or null when no next line exists
	 */
	@Override
	public final ListNum<byte[]> nextListAsync(){
		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}
		if(!Shared.SIMD){return nextListScalar();}
		
		final int listSize=200;
		final ArrayList<byte[]> list=new ArrayList<byte[]>(listSize);

		while(list.size()<listSize && open){
			if(listPos>=positions.size()){
				listPos=0;
				positions.clear();
				fillBuffer();
				Vector.findSymbols(buffer, 0, bstop, slashn, positions);

				// At EOF, a nonempty tail is one final line even without a newline position.
				if(positions.size()==0){
					if(bstop<1){break;}
					positions.add(bstop);
				}
			}

			// Process positions until we have enough lines OR run out of positions
			final int iters=Math.min(positions.size-listPos, listSize-list.size());
			for(int i=0; i<iters; i++){
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
				}
			}
		}

		return list.isEmpty() ? null : new ListNum<byte[]>(list, nextID++);
	}

	/** Prints meaningful buffered bytes to stderr, rendering LF and CR as visible markers. */
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

	/** Compacts the unread newline-free suffix and reads until a newline or nonpositive read.
	 * Doubles a full buffer when needed. Caught read exceptions latch error status
	 * even when diagnostic printing is suppressed.
	 * @return First LF offset, or bstop when no LF was found before reading ended
	 */
	private int fillBuffer(){
		if(bstart<bstop){//Shift end bytes to beginning
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
		int len=bstop;
		int r=-1;
		while(len==bstop){//hit end of input without encountering a newline
			if(bstop==buffer.length){
				buffer=KillSwitch.copyOf(buffer, buffer.length*2);
			}
			try{
				r=is.read(buffer, bstop, buffer.length-bstop);
			}catch(IOException e){//java.io.IOException: Stream Closed
				errorState=true;// FIO-087: a later successful close must not hide a failed read.
				//TODO: This should be avoided rather than caught.  It happens when a stream is shut down with e.g. "reads=100".
				if(!Shared.anomaly){
					e.printStackTrace();
					System.err.println("open="+open);
				}
			}catch(NullPointerException e){//Can be thrown by java.util.zip.Inflater.ensureOpen(Inflater.java:389)
				errorState=true;
				//TODO: This should be avoided rather than caught.  It happens when a stream is shut down with e.g. "reads=100".
				if(!Shared.anomaly){
					e.printStackTrace();
					System.err.println("open="+open);
				}
			}
			if(r>0){
				bstop=bstop+r;
				//Historical scanner also recognized old-style Mac line endings; this implementation uses LF.

				while(len<bstop && buffer[len]!=slashn){len++;}
				//Historical Vector.find alternative: x0.85 for short sequences, x1.00 for 150bp FASTQ.
			}else{
				len=bstop;
				break;
			}
		}

		assert(r==-1 || buffer[len]==slashn);
		return len;
	}

	/** Prepends the supplied bytes and LF to buffered input without changing lineNum.
	 * The cached-position path is not invalidated; that unresolved behavior is FIO-090.
	 * @param line Nonnull line bytes to copy into the input buffer
	 */
	@Override
	public void pushBack(byte[] line){
		//TODO: Probable bug - FIO-090: buffer edits leave cached SIMD positions stale; pushBack repair is deferred under IO-037.
		if(bstart>line.length){
			bstart--;
			buffer[bstart]='\n';
			for(int i=0, j=bstart-line.length; i<line.length; i++, j++){
				buffer[j]=line[i];
			}
			bstart=bstart-line.length;
			return;
		}

		int bLen=bstop-bstart;
		int newLen=bLen+line.length+1;
		int rShift=line.length+1-bstart;
		assert(rShift>0) : bstop+", "+bstart+", "+line.length;
		while(newLen>buffer.length){
			//This could get big if pushback is used often,
			//unless special steps are taken to prevent it, like leaving extra space for pushbacks.
			buffer=Arrays.copyOf(buffer, buffer.length*2);
		}

		Tools.shiftRight(buffer, rShift);

		for(int i=0; i<line.length; i++){
			buffer[i]=line[i];
		}
		buffer[line.length]='\n';
		bstart=0;
		bstop=newLen;
	}

	/** Opens the named input using current buffering/subprocess settings and clears cursors.
	 * Reuses allocated storage and preserves accumulated read/closure errors.
	 * @return Input stream retained in is
	 */
	private final synchronized InputStream open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened TextFile "+name());
		}
		open=true;
		is=ReadWrite.getInputStream(name(), BUFFERED, allowSubprocess(), true);
		bstart=-1;
		bstop=-1;
		// Cached newline offsets belong to the previous input buffer, not the reopened stream.
		positions.clear();
		listPos=0;
		return is;
	}

	/** Returns the local open flag, which can remain true after cached batch EOF. */
	@Override
	public boolean isOpen(){return open;}

	/** Returns the current input handle, or null after close; external reads share its cursor. */
	@Override
	public final InputStream is(){return is;}

	/** Returns the zero-based consumed-line counter, reset to -1 by close/open via reset. */
	@Override
	public final long lineNum(){return lineNum;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Whether this instance has requested an open input and not yet closed it. */
	private boolean open=false;
	/** Reused byte storage, initially sized from bufferlen and grown for long lines. */
	private byte[] buffer=new byte[bufferlen];
	/** Shared zero-length result for an empty line. */
	private static final byte[] blankLine=new byte[0];
	/** Unread byte interval [bstart,bstop), or -1 endpoints immediately after opening. */
	private int bstart=0, bstop=0;
	/** Active backend handle exposed for compatibility; close delegates cleanup to ReadWrite. */
	public InputStream is;
	/** Zero-based consumed-line count; close resets it to -1. */
	public long lineNum=-1;
	/** Cached LF offsets, with an optional bstop sentinel for a final unterminated batch line. */
	private IntList positions=new IntList();
	/** Next unread entry in positions. */
	private int listPos=0;

	/** Shared diagnostic verbosity. */
	public static boolean verbose=false;
	/** Whether new/opened input streams request ReadWrite buffering. */
	public static boolean BUFFERED=false;
	/** Initial byte capacity captured by newly constructed readers. */
	public static int bufferlen=65536;

	/** Sticky caught-read and finishReading error status, preserved across reset. */
	private boolean errorState=false;

}
