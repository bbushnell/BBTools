package fileIO;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ListNum;

/**
 * Legacy byte-line reader retained for I/O debugging and DiskBench's method=qf path.
 * Originally documented as a tool for investigating NERSC I/O slowdowns; its active
 * benchmark caller means this class is not an unreferenced deletion candidate.
 * Splits on LF, strips a preceding CR (also at unterminated EOF), and keeps empty
 * lines. Nonempty results copy bytes; empty results share a zero-length array.
 *
 * Intended for one reading owner. Synchronized batch/close methods do not make
 * scalar reading, reset, pushBack or public stream mutation generally thread-safe.
 * Read errors stop the current refill and latch a sticky flag;
 * callers must not treat returned partial data as verified successful input.
 * @author Brian Bushnell
 */
public class QuickFile{
	
	/*--------------------------------------------------------------*/
	/*----------------           Startup            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Prints a line range or measures remaining payload bytes, then closes/reopens/closes.
	 * Retains input-open-before-argument-parse error precedence. After construction,
	 * exceptional exits close the reader and preserve the primary failure. Ordinary
	 * completion checks accumulated input status and, in print mode, borrowed stdout.
	 * @param args Optional input name (default stdin), first index or speedtest, and last index
	 */
	public static void main(String[] args){
		QuickFile tf=new QuickFile(args.length>0 ? args[0] : "stdin", true);
		boolean speedtest=false, error=false;
		Throwable failure=null;
		try{
			long first=0, last=100;
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
			if(args.length>2){last=Integer.parseInt(args[2]);}
			speedtest(tf, first, last, !speedtest);
			// Retain the successful diagnostic reopen without losing earlier read errors.
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
				failure.addSuppressed(new RuntimeException("QuickFile: error reading or closing "+tf.name()));
			}
		}
		if(error){throw new RuntimeException("QuickFile: error reading or closing "+tf.name());}
		if(!speedtest && System.out.checkError()){
			throw new RuntimeException("QuickFile: error writing output");
		}
	}
	
	/**
	 * Skips first lines and processes through exclusive last, stopping at null.
	 * Counts returned payload bytes, excluding stripped terminators; printing uses the
	 * default charset and platform line separator. Does not close on exceptional exit.
	 * @param tf Borrowed reader
	 * @param first Number of initial nextLine calls; negative skips none
	 * @param last Exclusive logical bound, not an additional count
	 * @param reprint Print decoded lines and counters instead of throughput
	 */
	private static void speedtest(QuickFile tf, long first, long last, boolean reprint){
		Timer t=new Timer();
		long lines=0;
		long bytes=0;
		for(long i=0; i<first; i++){tf.nextLine();}
		if(reprint){
			for(long i=first; i<last; i++){
				byte[] s=tf.nextLine();
				if(s==null){break;}
				lines++;
				bytes+=s.length;
				System.out.println(new String(s));
			}
			System.err.println("\n");
			System.err.println("Lines: "+lines);
			System.err.println("Bytes: "+bytes);
		}else{
			for(long i=first; i<last; i++){
				byte[] s=tf.nextLine();
				if(s==null){break;}
				lines++;
				bytes+=s.length;
			}
		}
		t.stop();
		if(!reprint){
			System.err.println(Tools.timeLinesBytesProcessed(t, lines, bytes, 8));
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Resolves a text input descriptor and opens its stream immediately.
	 * @param fname Input name, including supported ReadWrite pseudo-inputs
	 * @param allowSubprocess_ Whether the descriptor permits subprocess decompression
	 */
	public QuickFile(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess_, false));
	}
	
	/**
	 * Retains the read descriptor and opens using the current BUFFERED setting.
	 * @param ff_ Nonnull descriptor asserted to be in read mode
	 */
	public QuickFile(FileFormat ff_){
		ff=ff_;
		assert(ff.read()) : ff;
		if(verbose){System.err.println("ByteFile1("+ff+")");}
		is=open();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Closes and reopens the same input name, clearing pushBack and batch IDs.
	 * Keeps sticky errorState; reopening stdin does not make it rewindable.
	 * Not atomic with concurrent scalar reads or public stream replacement.
	 */
	public final void reset(){
		close();
		is=open();
		pushBack=null;
		nextID=0;
	}
	
	/**
	 * Finishes the current stream through ReadWrite and accumulates its error status.
	 * A normal close clears the stream, line counter and pending pushBack reference;
	 * batch IDs are retained until reset. Already-closed calls return sticky status.
	 * @return True if this reader has recorded an error
	 */
	public synchronized final boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		if(!open){return errorState;}
		open=false;
		assert(is!=null);
		errorState|=ReadWrite.finishReading(is, name(), ff.subprocess);
		is=null;
		lineNum=-1;
		pushBack=null;
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		return errorState;
	}
	
	/**
	 * Returns the next LF-delimited payload, including empty lines, or null when closed/EOF.
	 * Strips a final CR from the payload even without LF; does not split on lone CR.
	 * Copies nonempty payloads and shares empty results. EOF automatically closes.
	 * The pushBack slot is currently not replayed, and read failure can resemble EOF
	 * or partial input; see the separate source findings rather than assuming success.
	 * @return Payload bytes, or null after exhaustion/closure
	 */
	public byte[] nextLine(){
		if(verbose){System.err.println("Reading line "+this.getClass().getName()+" for "+name()+"; open="+open+"; errorState="+errorState);}
		if(!open || is==null){
			if(Shared.WINDOWS){System.err.println("Attempting to read from a closed file: "+name());}
			return null;
		}
		int nlpos=bstart;
		while(nlpos<bstop && buffer[nlpos]!=slashn){nlpos++;}
		if(nlpos>=bstop){
			nlpos=fillBuffer();
		}
		if(nlpos<0 || bstop<1){
			close();
			return null;
		}
		lineNum++;
		//Limit is the position after the last position to copy.
		//Limit equals nlpos unless there was a \r before the \n.
		final int limit=(nlpos>bstart && buffer[nlpos-1]==slashr) ? nlpos-1 : nlpos;
		if(bstart==limit){//Empty line.
			bstart=nlpos+1;
			return blankLine;
		}
		byte[] line=KillSwitch.copyOfRange(buffer, bstart, limit);
		assert(line.length>0) : bstart+", "+nlpos+", "+limit;
		bstart=nlpos+1;
		return line;
	}
	
	/** Retained debugging allocation; unused by this reader. */
	final byte[] dummy=new byte[100];
	
	/**
	 * Compacts a trailing fragment and reads until LF or a nonpositive read result.
	 * Doubles a full buffer with unchecked integer capacity arithmetic. IOException
	 * is printed and latched; the successfully read prefix remains available.
	 * @return LF index or end of retained bytes; no terminator is included
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
			r=-1;//A failed read must not reuse the preceding byte count (FIO-059).
			try{
				r=is.read(buffer, bstop, buffer.length-bstop);
			}catch(IOException e){
				errorState=true;
				e.printStackTrace();
				System.err.println("open="+open);
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
		return len;
	}
	
	/**
	 * Opens the configured backend and initializes buffer cursors to the unopened sentinel.
	 * Sets open before acquisition; a thrown acquisition failure has no local rollback.
	 * @return Newly acquired stream, also assigned to is
	 */
	private final synchronized InputStream open(){
		if(open){
			throw new RuntimeException("Attempt to open already-opened TextFile "+name());
		}
		open=true;
		is=ReadWrite.getInputStream(name(), BUFFERED, allowSubprocess());
		bstart=-1;
		bstop=-1;
		return is;
	}
	
	/**
	 * Collects remaining nextLine results until null; normal EOF closes the input.
	 * @return New mutable list of payload arrays, with shared zero-length entries
	 */
	public final ArrayList<byte[]> toByteLines(){
		byte[] s=null;
		ArrayList<byte[]> list=new ArrayList<byte[]>(4096);
		for(s=nextLine(); s!=null; s=nextLine()){
			list.add(s);
		}
		return list;
	}
	
	/**
	 * Counts remaining lines, then resets after normal loop completion.
	 * The returned count is not necessarily the file total; reset is not stdin rewind.
	 * @return Number of nonnull nextLine results encountered
	 */
	public final long countLines(){
		byte[] s=null;
		long count=0;
		for(s=nextLine(); s!=null; s=nextLine()){count++;}
		reset();
		return count;
	}
	
	/**
	 * Collects at most 200 lines and assigns dense batch IDs starting at zero.
	 * Empty lines are entries; exhaustion before the first entry returns null.
	 * Synchronization does not cover independent calls to the scalar reader.
	 * @return Newly allocated batch wrapper/list, or null when no lines remain
	 */
	public synchronized final ListNum<byte[]> nextList(){
		byte[] line=nextLine();
		if(line==null){return null;}
		ArrayList<byte[]> list=new ArrayList<byte[]>(200);
		list.add(line);
		for(int i=1; i<200; i++){
			line=nextLine();
			if(line==null){break;}
			list.add(line);
		}
		ListNum<byte[]> ln=new ListNum<byte[]>(list, nextID);
		nextID++;
		return ln;
	}
	
	/**
	 * Accepts stdin/stdin.* and jar: names by prefix, otherwise tests filesystem existence.
	 * This is a name heuristic, not proof that an input can be opened or read.
	 * @return Whether the heuristic accepts the configured input name
	 */
	public final boolean exists(){
		return name().equals("stdin") || name().startsWith("stdin.") || name().startsWith("jar:") || new File(name()).exists();//TODO Ugly and unsafe hack for files in jars
	}
	
	/**
	 * Stores a borrowed reference in the dormant pushBack slot; nextLine does not replay it.
	 * Asserts that the slot is empty. Close/reset clear it. This is not working push-back.
	 * @param line Reference to store, possibly null
	 */
	public final void pushBack(byte[] line){
		//TODO: Probable bug - FIO-060: replay in nextLine is disabled, so stored bytes are
		//never emitted. Keep dormant until a supported caller warrants a separate repair.
		assert(pushBack==null);
		pushBack=line;
	}
	
	/** @return Configured input name from the borrowed descriptor. */
	public final String name(){return ff.name();}
	/** @return Whether the descriptor permits subprocess-backed input. */
	public final boolean allowSubprocess(){return ff.allowSubprocess();}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed read descriptor; backend flags come from it. */
	public final FileFormat ff;

	/** Legacy selector placeholder, not consulted by QuickFile. */
	public static boolean FORCE_MODE_BF1=false;//!(Shared.GENEPOOL || Shared.DENOVO || Shared.CORI || Shared.WINDOWS);
	/** Legacy selector placeholder, not consulted by QuickFile. */
	public static boolean FORCE_MODE_BF2=false;
	/** Legacy selector placeholder, not consulted by QuickFile. */
	public static boolean FORCE_MODE_BF3=false;
	
	/** LF/CR separators and retained legacy format-marker constants. */
	protected final static byte slashr='\r', slashn='\n', carrot='>', plus='+', at='@';//, tab='\t';

	/** Unused debugging scalar retained with its original initializer order. */
	long a=1;
	/** Unused debugging scalar. */
	long b=2;
	/** Unused debugging scalar. */
	long c=3;
	/** Unused debugging scalar. */
	long d=4;
	/** Unused debugging reference. */
	byte[] p0=null;
	/** Unused debugging reference. */
	byte[] p1=null;
	/** Unused debugging reference. */
	byte[] p2=null;
	/** Unused debugging reference. */
	byte[] p3=null;
	/** Dormant borrowed line slot; see FIO-060. */
	private byte[] pushBack=null;
	/** Next batch ID; reset restores zero. */
	private long nextID=0;
	
	/** Local acquisition/closure flag, set before backend acquisition succeeds. */
	private boolean open=false;
	/** Per-instance byte buffer, initially sized from bufferlen during construction. */
	private byte[] buffer=new byte[bufferlen];
	/** Shared immutable-length empty-line result. */
	private static final byte[] blankLine=new byte[0];
	/** Cursor and exclusive filled-byte bound; open assigns both -1. */
	private int bstart=0, bstop=0;
	/** Current backend; public replacement must not race reading/closure. */
	public InputStream is;
	/** Zero-based returned-line index; normal close resets it to -1. */
	public long lineNum=-1;
	
	/** Enables diagnostic prints. */
	public static boolean verbose=false;
	/** Buffering request read at each backend open. */
	public static boolean BUFFERED=false;
	/** Initial instance-buffer capacity, read during object initialization. */
	public static int bufferlen=16384;

	/** Sticky read/finishReading error flag; reset does not clear it. */
	private boolean errorState=false;
	
}
