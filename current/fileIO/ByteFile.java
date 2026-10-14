package fileIO;
import java.io.File;
import java.io.InputStream;
import java.util.ArrayList;

import shared.Shared;
import structures.ListNum;

/**
 * Byte-line reader contract, implementation selection, and collection helpers.
 * Factories select ByteFile1 through ByteFile4 using explicit choices, process-wide
 * flags, thread count and LOW_MEMORY; selection is a policy, not a per-file benchmark.
 * Concrete readers define buffering, line counters, reset support and lifecycle.
 * The inherited nextList method serializes its own callers and assigns record origins;
 * direct nextLine/nextListAsync calls, resets and overriding methods have separate
 * synchronization rules. Configure shared selection/batch settings before use.
 *
 * @author Brian Bushnell
 * @contributor Shinobu (reader-local preference and contracts)
 */
public abstract class ByteFile{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
//	public static final ByteFile makeByteFile(String fname){
//		return makeByteFile(fname, false, true);
//	}
	
	/**
	 * Opens ByteFile1 directly, regardless of selection and allow flags.
	 *
	 * @param fname The filename to read
	 * @param allowSubprocess Whether to allow subprocess decompression
	 * @return ByteFile1 instance for the specified file
	 */
	public static final ByteFile makeByteFile1(String fname, boolean allowSubprocess){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess, false);
		return new ByteFile1(ff);
	}
	
	/**
	 * Creates an input descriptor with TEXT as fallback and uses the default reader policy.
	 *
	 * @param fname The filename to read
	 * @param allowSubprocess Whether to allow subprocess decompression
	 * @return Selected ByteFile1, ByteFile2, ByteFile3 or ByteFile4 instance
	 */
	public static final ByteFile makeByteFile(String fname, boolean allowSubprocess){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TEXT, null, allowSubprocess, false);
		return makeByteFile(ff);
	}
	
	/**
	 * Opens a reader using forced-mode or automatic selection.
	 * @param ff The FileFormat describing the input file
	 * @return ByteFile instance for the specified format
	 */
	public static final ByteFile makeByteFile(FileFormat ff){
		return makeByteFile(ff, 0);
	}
	
	/**
	 * Creates a reader with explicit type selection. An allowed explicit type
	 * takes priority over global force flags. Type 0 or an unavailable explicit
	 * type uses forced/automatic selection; automatic selection can choose BF4.
	 *
	 * @param ff The FileFormat describing the input file
	 * @param type Implementation type (0=auto, 1-4=ByteFile1-4)
	 * @return ByteFile instance of the specified or selected type
	 */
	public static final ByteFile makeByteFile(FileFormat ff, int type){
		type=pickType(type);
		if(type==4){return new ByteFile4(ff);}
		if(type==3){return new ByteFile3(ff);}
		if(type==2){return new ByteFile2(ff);}
		return new ByteFile1(ff);
	}

	/**
	 * Opens a reader with a preference for this call, preserving global force-flag priority.
	 * A forced mode wins; otherwise an allowed preferred type wins over automatic selection.
	 * Type 0 or an unavailable preference falls back to the existing automatic policy.
	 * Does not change any shared selection flags. Configure shared flags before opening readers.
	 * @param ff Input format descriptor
	 * @param preferredType Preferred implementation (0=auto, 1-4=ByteFile1-4)
	 * @return Newly opened reader
	 * @throws IllegalArgumentException if preferredType is outside 0-4
	 */
	public static final ByteFile makeByteFileWithPreference(final FileFormat ff, final int preferredType){
		if(preferredType<0 || preferredType>4){
			throw new IllegalArgumentException("ByteFile preference must be 0-4: "+preferredType);
		}
		final int forced=forcedType();
		return makeByteFile(ff, forced>0 ? forced : preferredType);
	}
	
	/**
	 * Borrows the input descriptor and asserts that it is configured for reading.
	 * @param ff_ The FileFormat for this ByteFile instance
	 */
	protected ByteFile(FileFormat ff_){
		ff=ff_;
		assert(ff.read()) : ff;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Collects remaining lines from the current position using nextLine.
	 * Retains the returned arrays without copying them. Does not explicitly close or
	 * reset the reader; backend EOF handling still applies.
	 * @return Newly allocated list of remaining line arrays
	 */
	public final ArrayList<byte[]> toByteLines(){
		byte[] s=null;
		ArrayList<byte[]> list=new ArrayList<byte[]>(4096);
		for(s=nextLine(); s!=null; s=nextLine()){list.add(s);}
		return list;
	}
	
	/** Collects remaining batches through nextListAsync, retaining their line arrays.
	 * Does not use nextList's monitor or record-origin numbering, or explicitly close
	 * or reset the reader. Overrides determine batching behavior.
	 * @return Newly allocated list of remaining line arrays */
	public final ArrayList<byte[]> toByteLinesList(){
		ArrayList<byte[]> list=new ArrayList<byte[]>(TARGET_LIST_SIZE*8);
		for(ListNum<byte[]> ln=nextListAsync(); ln!=null; ln=nextListAsync()){
			list.addAll(ln.list);
		}
		return list;
	}
	
	/**
	 * Opens the selected reader, collects lines, and closes it before returning.
	 * Reported read/close errors reject the collected data. A collection exception
	 * remains primary if closing also fails; the closure failure is suppressed on it.
	 * @param ff The FileFormat to read from
	 * @return ArrayList containing all lines as byte arrays
	 * @throws RuntimeException if the reader reports an input or closure error
	 */
	public static final ArrayList<byte[]> toLines(FileFormat ff){
		final ByteFile bf=makeByteFile(ff);
		Throwable failure=null;
		try{
			return bf.toByteLines();
		}catch(RuntimeException | Error t){
			failure=t;
			throw t;
		}finally{
			//FIO-002: close on rejection and do not return data after a latched error.
			try{
				if(bf.close()){throw new RuntimeException("Error reading or closing "+bf.name());}
			}catch(RuntimeException | Error t){
				if(failure==null){throw t;}
				if(failure!=t){failure.addSuppressed(t);}
			}
		}
	}
	
	/**
	 * Creates an input descriptor with TEXT as fallback and subprocesses allowed,
	 * then delegates to toLines.
	 * Shares that helper's closure and error-reporting contract.
	 * @param fname The filename to read
	 * @return ArrayList containing all lines as byte arrays
	 */
	public static final ArrayList<byte[]> toLines(String fname){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TEXT, null, true, false);
		return toLines(ff);
	}
	
	/**
	 * Consumes and counts remaining lines, then requests reset from the backend.
	 * Requires a reset-capable reader to restore subsequent reading; not every reader
	 * selected by the factory supports that operation. Does not check a close-error flag.
	 * @return Number of lines consumed from the current position
	 */
	public final long countLines(){
		byte[] s=null;
		long count=0;
		for(s=nextLine(); s!=null; s=nextLine()){count++;}
		//TODO: Probable bug (FIO-003) - automatic selection can return ByteFile4, whose reset cannot restart reading.
		reset();
		return count;
	}

	/** Requests rewind/reopening according to the concrete reader's reset support.
	 * Reset is not supported by every implementation returned by the factory. */
	public abstract void reset();
	/** Resets this base class's list-ID and record-origin counters to zero.
	 * Does not reopen input or reset subclass state. */
	final void superReset(){
		nextID=0;
		nextRecordNum=0;
	}
	
	/**
	 * Serializes calls through this inherited implementation, delegates batching to
	 * nextListAsync, and sets firstRecordNum from the cumulative size of its returned
	 * batches. The delegate supplies the list ID and batch limits. Direct calls to
	 * other consumption methods are not coordinated by this monitor and are not
	 * included in its record-origin count. Subclasses can override this method.
	 * @return Next numbered batch, or null when the delegate has no more batches
	 */
	public ListNum<byte[]> nextList(){
		synchronized(nlSync){
			return numberList(nextListAsync());
		}
	}

	/** Assigns the zero-based record origin for lists returned by {@link #nextList()}. */
	private final ListNum<byte[]> numberList(ListNum<byte[]> ln){
		if(ln!=null){
			ln.firstRecordNum=nextRecordNum;
			nextRecordNum+=ln.size();
		}
		return ln;
	}
	
	/**
	 * Builds a batch without synchronization by repeatedly calling nextLine.
	 * With positive limits, includes at most TARGET_LIST_SIZE lines and stops after
	 * their total length reaches TARGET_LIST_BYTES. The last whole line can exceed
	 * the byte threshold; line terminators are not added to the byte count.
	 * Assigns successive list IDs from nextID, but does not set firstRecordNum.
	 * Subclass overrides can use different batching limits and mechanisms.
	 * @return Batch containing at least one line, or null if the first read returns null
	 */
	public ListNum<byte[]> nextListAsync(){
		byte[] line=nextLine();
		if(line==null){return null;}
		final int slimit=TARGET_LIST_SIZE, blimit=TARGET_LIST_BYTES;
		ArrayList<byte[]> list=new ArrayList<byte[]>(slimit);
		list.add(line);
		int bytes=line.length;
		
		for(int i=1; i<slimit && bytes<blimit; i++){
			line=nextLine();
			if(line==null){break;}
			list.add(line);
			bytes+=line.length;
		}
		ListNum<byte[]> ln=new ListNum<byte[]>(list, nextID);
		nextID++;
		return ln;
	}
	
	/**
	 * Recognizes stdin/stdin.* and jar: names without opening them; otherwise uses
	 * File.exists. A true result does not establish readability or a regular file.
	 * @return True for a recognized special name or an existing filesystem entry
	 */
	public final boolean exists(){
		return name().equals("stdin") || name().startsWith("stdin.") || name().startsWith("jar:") || new File(name()).exists(); //TODO Ugly and unsafe hack for files in jars
	}

	/**
	 * Exposes the backend's input stream, which may be null after closure.
	 * Buffered or producer-driven readers can have consumed bytes ahead of returned
	 * lines; direct access must be coordinated with that backend.
	 * @return Backend-owned input stream, or null as defined by the implementation
	 */
	public abstract InputStream is();
	/**
	 * Returns the backend's line counter. Numbering, reset/closure values and whether
	 * it observes producer progress or delivered lines depend on the implementation.
	 * @return Implementation-specific line counter; not a universal one-based index
	 */
	public abstract long lineNum();
	
	/** Requests closure according to the backend's lifecycle.
	 * @return True if the backend reports a read or close error */
	public abstract boolean close();
	
	/**
	 * Returns the next line according to the concrete reader's delimiter and buffering
	 * rules. Empty arrays can represent blank lines; null signals no next line.
	 * EOF alone does not establish an error-free read; inspect close's error result.
	 * @return Next line array, possibly empty, or null
	 */
	public abstract byte[] nextLine();
	
//	public final void pushBack(byte[] line){
//		assert(pushBack==null);
//		pushBack=line;
//	}
	
	/**
	 * Requests that a line be returned again by nextLine. Support, copying, counter
	 * effects and interaction with batching depend on the backend; some reject pushback.
	 * @param line The line to push back for re-reading
	 */
	public abstract void pushBack(byte[] line);
	
	/**
	 * Returns whether the file is currently open for reading.
	 * Implementation varies by subclass.
	 * @return true if the file is open
	 */
	public abstract boolean isOpen();
	
	/** Returns the filename from the associated FileFormat */
	public final String name(){return ff.name();}
	/** Returns whether subprocess decompression is allowed */
	public final boolean allowSubprocess(){return ff.allowSubprocess();}
	
	/** Applies explicit type, forced mode, then automatic selection in that order.
	 * Disallowed explicit modes fall back to policy; force flags bypass allow checks.
	 * Automatic low-memory/fewer-than12-thread selection always uses BF1. Otherwise
	 * it tries allowed BF4, BF3, BF2, then falls back to BF1. */
	private static final int pickType(int type){
		if(type==4 && !ALLOW_BF4){type=0;}
		else if(type==3 && !ALLOW_BF3){type=0;}
		else if(type==2 && !ALLOW_BF2){type=0;}
		else if(type==1 && !ALLOW_BF1){type=0;}
		if(type>0){return type;}
		assert(type==0) : type;
		
		final int threads=Shared.threads();
		final int forced=forcedType();
		if(forced>0){return forced;}
		
		if(Shared.LOW_MEMORY || threads<12){return 1;}
		if(ALLOW_BF4){return 4;}
		if(ALLOW_BF3){return 3;}
		if(ALLOW_BF2){return 2;}
		return 1;
	}

	/** Returns the first forced mode in BF1/BF4/BF3/BF2 priority, or zero. */
	private static final int forcedType(){
		if(FORCE_MODE_BF1){return 1;}
		if(FORCE_MODE_BF4){return 4;}
		if(FORCE_MODE_BF3){return 3;}
		if(FORCE_MODE_BF2){return 2;}
		return 0;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** The FileFormat describing this file's characteristics and location */
	public final FileFormat ff;
	/** Per-reader monitor used only by the inherited nextList implementation. */
	private final String nlSync=new String(hashCode()+"");
	
	/** Highest-priority process-wide forced mode; consulted after an explicit type. */
	public static boolean FORCE_MODE_BF1=false;
	/** Force BF2 when no higher-priority forced mode or allowed explicit type wins. */
	public static boolean FORCE_MODE_BF2=false;
	/** Force usage of ByteFile3 implementation */
	public static boolean FORCE_MODE_BF3=false;
	/** Force usage of ByteFile4 implementation */
	public static boolean FORCE_MODE_BF4=false;

	/** Allows explicit BF1 selection; automatic fallback can still choose BF1 when false. */
	public static boolean ALLOW_BF1=true;
	/** Allows explicit/automatic BF2 selection; force flags take precedence. */
	public static boolean ALLOW_BF2=true;
	/** Allows explicit/automatic BF3 selection; disabled pending the historical FASTA concern. */
	public static boolean ALLOW_BF3=false;//Hung on fasta input
	/** Allows explicit/automatic BF4 selection; callers remain responsible for closure. */
	public static boolean ALLOW_BF4=true;//Hangs if there are limited reads and it is never closed.
	
	/** Shared ASCII delimiters and common FASTA/FASTQ marker bytes. */
	protected final static byte slashr='\r', slashn='\n', carrot='>', plus='+', at='@';//, tab='\t';
	/** Shared one-byte plus marker. */
	protected static final byte[] plusLine=new byte[] {plus};
	
	/** Positive line limit used by base batching and selected backend paths. */
	public static int TARGET_LIST_SIZE=800;
	/** Positive whole-line byte threshold used by base batching and selected backends. */
	public static int TARGET_LIST_BYTES=262144;
	
//	byte[] pushBack=null;
	/** Next batch ID for implementations using this base counter. */
	protected long nextID=0;
	/** Cumulative lines returned through the inherited nextList wrapper since superReset. */
	private long nextRecordNum=0;
	
}
