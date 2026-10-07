package stream;

import java.util.ArrayList;
import java.util.Set;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.FileFormat;
import parse.Parse;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Base class for routing named read batches to buffered multi-destination outputs.
 * List submissions route nonnull Read.obj values, cast to String, through subclass buffers.
 * Nonthreaded mode routes on the caller; threaded mode transfers lists to this Thread,
 * which the caller must start before using its termination-and-wait operations.
 * Subclasses own destination writers, dumping thresholds, residual handling and reporting.
 * This outer transfer mode is separate from any threads used by destination writers.
 * @author Brian Bushnell
 * @date May 14, 2019
 */
public abstract class BufferedMultiCros extends Thread{

	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates the configured default implementation with current default threading and stream limits.
	 * Passes explicit output options unchanged to the fully specified factory.
	 * @return The constructed implementation; this method does not call its start method
	 */
	public static BufferedMultiCros make(String out1, String out2, boolean overwrite, boolean append,
			boolean allowSubprocess, boolean useSharedHeader, int defaultFormat){
		return make(out1, out2, overwrite, append, allowSubprocess, useSharedHeader,
				defaultFormat, defaultThreaded, defaultMcrosType, defaultMaxStreams);
	}

	/**
	 * Selects MultiCros2 through MultiCros6 for types 2 through 6, or MultiWriter for type 7.
	 * Forwards output options to the chosen constructor. Type 2 fixes its stream limit at one;
	 * the others receive maxStreams. Constructor side effects belong to the selected class;
	 * in particular, MultiCros4 starts retirement helpers during construction.
	 * @param out1 Primary output pattern containing a percent placeholder
	 * @param out2 Optional second output pattern
	 * @param overwrite Whether existing destination files may be replaced
	 * @param append Whether appending is permitted; concrete implementations choose reopening policy
	 * @param allowSubprocess Whether output descriptors permit subprocesses
	 * @param useSharedHeader Whether to request shared SAM headers
	 * @param defaultFormat Fallback output format
	 * @param threaded Whether the outer list-submission path uses a transfer thread
	 * @param mcrosType Implementation selector, from 2 through 7
	 * @param maxStreams Requested limit on open destinations
	 * @return The constructed implementation
	 * @throws RuntimeException If mcrosType is not supported
	 */
	public static BufferedMultiCros make(String out1, String out2, boolean overwrite, boolean append,
			boolean allowSubprocess, boolean useSharedHeader, int defaultFormat, boolean threaded,
			int mcrosType, int maxStreams){

		BufferedMultiCros mcros=null;
		if(mcrosType==2){//Slow, synchronous mcros type
			mcros=new MultiCros2(out1, out2, overwrite, append, allowSubprocess, useSharedHeader, defaultFormat, threaded);
		}else if(mcrosType==3){//Faster, asynchronous type
			mcros=new MultiCros3(out1, out2, overwrite, append, allowSubprocess, useSharedHeader, defaultFormat, threaded, maxStreams);
		}else if(mcrosType==4){//Threaded file closing
			mcros=new MultiCros4(out1, out2, overwrite, append, allowSubprocess, useSharedHeader, defaultFormat, threaded, maxStreams);
		}else if(mcrosType==5){//New retirement ordering by timer
			mcros=new MultiCros5(out1, out2, overwrite, append, allowSubprocess, useSharedHeader, defaultFormat, threaded, maxStreams);
		}else if(mcrosType==6){//New retirement ordering by timer
			mcros=new MultiCros6(out1, out2, overwrite, append, allowSubprocess, useSharedHeader, defaultFormat, threaded, maxStreams);
		}else if(mcrosType==7){//Writer-based; same structure as 6 with ZT Writers instead of CROS
			mcros=new MultiWriter(out1, out2, overwrite, append, allowSubprocess, useSharedHeader, defaultFormat, threaded, maxStreams);
		}else{
			throw new RuntimeException("Bad mcrosType: "+mcrosType);
		}
		return mcros;
	}

	/**
	 * Captures output options and computes queue, retirement and memory settings.
	 * Requires a percent placeholder in each supplied pattern. With no explicit second
	 * pattern, a hash placeholder in the first is expanded into separate 1 and 2 patterns.
	 * Threaded mode allocates a capacity-eight transfer queue; this constructor does not
	 * start the outer thread. Memory thresholds have fixed minimums plus configured
	 * fractions of currently available memory; subclasses decide how to apply them.
	 * @param pattern1_ Nonnull primary output pattern
	 * @param pattern2_ Optional second output pattern
	 * @param overwrite_ Whether existing outputs may be replaced
	 * @param append_ Stored append permission
	 * @param allowSubprocess_ Stored subprocess permission
	 * @param useSharedHeader_ Stored shared-header preference
	 * @param defaultFormat_ Fallback output format
	 * @param threaded_ Whether to allocate the outer transfer queue
	 * @param maxStreams_ Destination limit used to derive the retirement batch size
	 */
	public BufferedMultiCros(String pattern1_, String pattern2_,
			boolean overwrite_, boolean append_, boolean allowSubprocess_, boolean useSharedHeader_,
			int defaultFormat_, boolean threaded_, int maxStreams_){
		assert(pattern1_!=null && pattern1_.indexOf('%')>=0);
		assert(pattern2_==null || pattern2_.indexOf('%')>=0); //#001 FIXED (author had marked "Possible bug"): was pattern1_.indexOf, which is redundant with the assert above and never validated pattern2_. A non-null pattern2_ must contain '%' or its per-name replaceFirst is a no-op -> all R2 output collapses into one file. (Twin of the same bug in old MultiCros.java ctor.)

		//Perform # expansion for twin files
		if(pattern2_==null && pattern1_.indexOf('#')>=0){
			pattern1=pattern1_.replaceFirst("#", "1");
			pattern2=pattern1_.replaceFirst("#", "2");
		}else{
			pattern1=pattern1_;
			pattern2=pattern2_;
		}

		overwrite=overwrite_;
		append=append_;
		allowSubprocess=allowSubprocess_;
		useSharedHeader=useSharedHeader_;

		defaultFormat=defaultFormat_;

		threaded=threaded_;
		transferQueue=threaded ? new ArrayBlockingQueue<ArrayList<Read>>(8) : null;
		maxStreams=maxStreams_;

		//Significantly impacts performance.
		//Higher numbers give more retires but less time per retire.
		//Optimal seems to be around 4-6, at least for 16 streams.
		streamsToRetire=Tools.mid(2, (maxStreams+1)/3, 16);

		final long bytes=Shared.memAvailable();
		memLimitLower=Tools.max(50000000, (long)(memLimitLowerMult*bytes));
		memLimitMid=Tools.max(70000000, (long)(memLimitMidMult*bytes));
		memLimitUpper=Tools.max(90000000, (long)(memLimitUpperMult*bytes));
	}

	/*--------------------------------------------------------------*/
	/*----------------           Parsing            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses shared multi-output defaults and memory fractions, affecting later instances.
	 * Buffer-size keys accept mixed case; implementation, threading and stream-limit keys
	 * are compared literally. Unrecognized keys return false without changing settings.
	 * @param arg Original argument, currently unused
	 * @param a Parsed key
	 * @param b Parsed value
	 * @return Whether the key was recognized
	 */
	public static boolean parseStatic(String arg, String a, String b){
		if(a.equals("mcrostype")){
			defaultMcrosType=Integer.parseInt(b);
		}else if(a.equals("threaded")){
			defaultThreaded=Parse.parseBoolean(b);
		}else if(a.equals("streams")){
			defaultMaxStreams=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("readsPerBuffer")){
			defaultReadsPerBuffer=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("bytesPerBuffer")){
			defaultBytesPerBuffer=Integer.parseInt(b);
		}else if(a.equalsIgnoreCase("memLimitLowerMult") || a.equals("mllmult") || a.equals("mllm")){
			memLimitLowerMult=Float.parseFloat(b);
			assert(memLimitLowerMult>=0 && memLimitLowerMult<1);
		}else if(a.equalsIgnoreCase("memLimitMidMult") || a.equals("mlmmult") || a.equals("mlmm")){
			memLimitMidMult=Float.parseFloat(b);
			assert(memLimitMidMult>=0 && memLimitMidMult<1);
		}else if(a.equalsIgnoreCase("memLimitUpperMult") || a.equals("mlumult") || a.equals("mlum")){
			memLimitUpperMult=Float.parseFloat(b);
			assert(memLimitUpperMult>=0 && memLimitUpperMult<1);
		}else{
			return false;
		}

		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Abstract Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the concrete implementation's success indication; this query is not a completion wait. */
	public abstract boolean finishedSuccessfully();

	/**
	 * Adds a single read to the specified output buffer.
	 * Called by this base class's routing path, including on its transfer thread.
	 * External callers must not bypass that queue in threaded mode.
	 * @param r Read to add to buffer
	 * @param name Name of destination buffer
	 */
	abstract void add(Read r, String name);

	/** Asks the subclass to dump eligible buffers; return units follow that implementation. */
	abstract long dumpAll();

	/**
	 * Dumps all residual reads to the specified output stream.
	 * @param rosu Destination stream for residual reads
	 * @return Number of residual reads that were dumped
	 */
	public abstract long dumpResidual(ConcurrentReadOutputStream rosu);

	/** Runs concrete output draining and finalization after normal submissions finish.
	 * @return Implementation-specific dump count */
	abstract long closeInner();

	/** Generates a report on how many reads went to each output file.
	 * @return ByteBuilder containing the formatted report */
	public abstract ByteBuilder report();

	/**
	 * Gets timing information for shutting down output threads.
	 * Default implementation throws RuntimeException for unsupported subclasses.
	 * @return Formatted timing information
	 * @throws RuntimeException if not implemented by subclass
	 */
	public String printRetireTime(){
		throw new RuntimeException("printRetireTime not available for "+getClass().getName());
	}

	/**
	 * Gets timing information for creating output threads.
	 * Default implementation throws RuntimeException for unsupported subclasses.
	 * @return Formatted timing information
	 * @throws RuntimeException if not implemented by subclass
	 */
	public String printCreateTime(){
		throw new RuntimeException("printCreateTime not available for "+getClass().getName());//#003 fix: message said printRetireTime (copy-paste)
	}

	/** Returns destination keys under the concrete implementation's collection-view contract. */
	public abstract Set<String> getKeys();

	/*--------------------------------------------------------------*/
	/*----------------        Final Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Completes normal output after submissions have stopped.
	 * Threaded mode queues its terminal marker and waits for this started thread;
	 * nonthreaded mode invokes closeInner directly. Repeated-close behavior belongs
	 * to the concrete implementation; this base does not make close idempotent.
	 */
	public final void close(){
		if(threaded){poisonAndWait();}
		else{closeInner();}
	}

	/** Gets the primary file pattern.
	 * @return The pattern1 file pattern string */
	public final String fname(){return pattern1;}

	/** Checks if this stream has detected an error condition.
	 * @return true if an error state has been detected */
	public final boolean errorState(){return errorState;}

	/**
	 * Submits a list for named-read routing without copying the list or its reads.
	 * Keep submitted contents stable while routing and downstream writers still use them;
	 * subclasses may retain references after dequeue or a synchronous routing return.
	 * Nonthreaded mode routes immediately and then invokes the subclass load check.
	 * @param list Nonnull list containing nonnull reads; unnamed reads are ignored by routing
	 */
	public final void add(ArrayList<Read> list){
		if(threaded){//Send to the transfer queue
			addToQueue(list);
		}else{//Add the reads from this thread
			addToBuffers(list);
		}
	}

	/**
	 * Distributes individual reads to their designated buffers.
	 * Ignores reads whose obj is null; otherwise casts obj to String and delegates by name.
	 * Increments readsInTotal once per routed list entry, excluding separate mate counts,
	 * and performs the subclass load check after the complete list.
	 * @param list List of reads to add to buffers
	 */
	private final void addToBuffers(ArrayList<Read> list){
		for(Read r : list){
			if(r.obj!=null){
				String name=(String)r.obj;
				readsInTotal++;
				add(r, name);//Reads without a name in the obj field get ignored here.
			}
		}
		handleLoad0();
	}

	/** Called after adding a list of reads to handle load management.
	 * Default implementation does nothing; subclasses may override. */
	void handleLoad0(){
		//Do nothing
	}

	/*--------------------------------------------------------------*/
	/*----------------       Threaded Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Routes lists until the identity-based terminal token arrives, then invokes closeInner.
	 * Requires threaded mode. Interrupted queue takes invoke exceptionKill; other escaping
	 * Throwables also set errorState before exceptionKill. Historical rationale is retained
	 * beside the catch block. This method is the started outer thread's entry point.
	 */
	@Override
	public final void run(){
		assert(threaded) : "This should only be called in threaded mode.";
		try{
			for(ArrayList<Read> list=transferQueue.take(); list!=poisonToken; list=transferQueue.take()){
				if(verbose){System.err.println("Got list; size="+transferQueue.size());}//#002 fix: a misplaced escaped-quote made this print the LITERAL text 'size="+transferQueue.size())' instead of the value (debug-only).
				addToBuffers(list);
				if(verbose){System.err.println("Added list; size="+transferQueue.size());}
			}
			closeInner();
		}catch(InterruptedException e){
			//Terminate JVM if something goes wrong
			KillSwitch.exceptionKill(e);
		}catch(Throwable t){
			//An uncaught Throwable here previously just killed THIS thread; waitForFinish()'s join then
			//saw a TERMINATED thread and treated it as success, so a failed dump (e.g. unwritable output
			//file) reported full yield and exit 0 with NOTHING written (replicated via demuxbyname.sh,
			//blocked-output-file injection, 2026-09-05; nonthreaded mode was already loud). Crash the
			//whole process instead: never report success on lost output.
			errorState=true;
			KillSwitch.exceptionKill(t);
		}
	}

	/** Queues the terminal token after the final submission in threaded mode; does not wait for termination. */
	public final void poison(){
		assert(threaded) : "This should only be called in threaded mode.";
		addToQueue(poisonToken);
	}

	/**
	 * Puts a list into the transfer queue, retrying at most ten interrupted attempts.
	 * Each put may wait for capacity; ten attempts are not a time limit. Exhausting the
	 * attempts invokes KillSwitch.kill. Requires the threaded-mode queue to exist.
	 * @param list Read batch or the private terminal token
	 * @return Whether a queue put completed
	 */
	boolean addToQueue(ArrayList<Read> list){
		boolean success=false;
		for(int i=0; i<10 && !success; i++){
			try{
				transferQueue.put(list);
				success=true;
			}catch(InterruptedException e){
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		if(!success){
			KillSwitch.kill("Something went wrong when adding to "+getClass().getName());
		}
		return success;
	}

	/** Queues termination after normal submissions and waits for the already-started outer thread. */
	public final void poisonAndWait(){
		assert(threaded) : "This should only be called in threaded mode.";
		poison();
		waitForFinish();
	}

	/**
	 * Waits until this started object's state is TERMINATED, retrying one-second joins.
	 * Requires threaded mode; does not start the thread or enqueue a terminal marker.
	 * Interrupted joins are printed and retried. This wait alone does not report success.
	 */
	public final void waitForFinish(){
		assert(threaded);
		if(verbose){System.err.println("Waiting for finish.");}
		while(this.getState()!=Thread.State.TERMINATED){
			if(verbose){System.err.println("Attempting join.");}
			try{
				this.join(1000);
			}catch(InterruptedException e){
				e.printStackTrace();
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolved destination patterns after optional twin-file hash expansion. */
	public final String pattern1, pattern2;

	/** True if an error was encountered during processing */
	boolean errorState=false;

	/** Stored permission for concrete implementations to replace existing output. */
	final boolean overwrite;

	/** Permission to append to existing files */
	final boolean append;

	/** Permission to spawn subprocesses (e.g., for pigz compression) */
	final boolean allowSubprocess;

	/** Default output file format when unclear from file extension */
	final int defaultFormat;

	/** Writer-buffer setting consumed by concrete implementations; some factories ignore it. */
	int rswBuffers=1;

	/** Whether to print the shared header in SAM output files */
	final boolean useSharedHeader;

	/** Lower estimated-byte threshold supplied to subclass load/retirement policy. */
	final long memLimitLower;

	/** Middle estimated-byte threshold supplied to subclass load policy. */
	final long memLimitMid;

	/** Upper estimated-byte threshold supplied to subclass load policy. */
	final long memLimitUpper;

	/** Maximum number of active streams allowed for MCros3+ implementations */
	public final int maxStreams;

	/** Number of streams to retire at a time for load balancing */
	public final int streamsToRetire;

	/** Per-buffer entry setting interpreted by subclasses; some use it only for initial capacity. */
	public int readsPerBuffer=defaultReadsPerBuffer;

	/** Estimated-byte target interpreted by subclass buffering policy. */
	public int bytesPerBuffer=defaultBytesPerBuffer;

	/** Minimum destination read total for output; concrete implementations define counting units. */
	public long minReadsToDump=0;

	/** Residual totals maintained by each concrete implementation. */
	public long residualReads=0, residualBases=0;

	/** Number of named list entries routed by addToBuffers, excluding separate mate counts. */
	long readsInTotal=0;

	/** Subclass-maintained in-flight read counter; units and updates belong to the implementation. */
	long readsInFlight=0;

	/** Current number of estimated bytes held in buffers */
	long bytesInFlight=0;

	/** Queue for transferring read lists when MultiCros runs in threaded mode */
	private final ArrayBlockingQueue<ArrayList<Read>> transferQueue;

	/** Special empty list used to signal thread termination in threaded mode */
	private final ArrayList<Read> poisonToken=new ArrayList<Read>(0);

	/** Whether list submission uses the outer transfer queue and Thread. */
	public final boolean threaded;

	/** Whether to use LogLog for tracking cardinality of each output file */
	public boolean trackCardinality=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Available-memory fraction for the lower threshold; a fixed minimum also applies. */
	private static float memLimitLowerMult=0.20f;
	/** Available-memory fraction for the middle threshold; a fixed minimum also applies. */
	private static float memLimitMidMult=0.40f;
	/** Available-memory fraction for the upper threshold; a fixed minimum also applies. */
	private static float memLimitUpperMult=0.60f;
	/** Default outer transfer mode used by the short factory overload. */
	public static boolean defaultThreaded=true;
	/** Default requested destination limit used by the short factory overload. */
	public static int defaultMaxStreams=12;
	/** Default implementation selector used by the short factory overload. */
	public static int defaultMcrosType=6;
	/** Initial readsPerBuffer setting for newly constructed instances. */
	public static int defaultReadsPerBuffer=32000;
	/** Initial bytesPerBuffer setting for newly constructed instances. */
	public static int defaultBytesPerBuffer=16000000;

	/** Shared diagnostic-output switch. */
	public static boolean verbose=false;

}
