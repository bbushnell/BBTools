package stream;

import java.util.ArrayList;
import java.util.Set;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.FileFormat;
import shared.KillSwitch;
import shared.Parse;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Allows output of reads to multiple different output streams.
 * Each output stream is controlled by a buffer,
 * which stores reads until there is a sufficient quantity to dump.
 * 
 * @author Brian Bushnell
 * @date May 14, 2019
 *
 */
public abstract class BufferedMultiCros extends Thread {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructors         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static BufferedMultiCros make(String out1, String out2, boolean overwrite, boolean append, 
			boolean allowSubprocess, boolean useSharedHeader, int defaultFormat) {
		return make(out1, out2, overwrite, append, allowSubprocess, useSharedHeader, 
				defaultFormat, defaultThreaded, defaultMcrosType, defaultMaxStreams);
	}
	
	public static BufferedMultiCros make(String out1, String out2, boolean overwrite, boolean append, 
			boolean allowSubprocess, boolean useSharedHeader, int defaultFormat, boolean threaded,
			int mcrosType, int maxStreams) {

		BufferedMultiCros mcros=null;
		if(mcrosType==2){//Slow, synchronous mcros type
			mcros=new MultiCros2(out1, out2, overwrite, append, true, useSharedHeader, FileFormat.FASTQ, threaded);
		}else if(mcrosType==3){//Faster, asynchronous type
			mcros=new MultiCros3(out1, out2, overwrite, append, true, useSharedHeader, FileFormat.FASTQ, threaded, maxStreams);
		}else if(mcrosType==4){//Threaded file closing
			mcros=new MultiCros4(out1, out2, overwrite, append, true, useSharedHeader, FileFormat.FASTQ, threaded, maxStreams);
		}else if(mcrosType==5){//New retirement ordering by timer
			mcros=new MultiCros5(out1, out2, overwrite, append, true, useSharedHeader, FileFormat.FASTQ, threaded, maxStreams);
		}else if(mcrosType==6){//New retirement ordering by timer
			mcros=new MultiCros6(out1, out2, overwrite, append, true, useSharedHeader, FileFormat.FASTQ, threaded, maxStreams);
		}else{
			throw new RuntimeException("Bad mcrosType: "+mcrosType);
		}
		return mcros;
	}
	
	/**
	 * Primary constructor.
	 * @param pattern1_ Name pattern for file 1; must contain % (required)
	 * @param pattern2_ Name pattern for file 2; must contain % (optional)
	 * @param overwrite_ Permission to overwrite
	 * @param append_ Permission to append to existing files (this should generally be false)
	 * @param allowSubprocess_ Allow subprocesses such as pigz, bgzip, or samtools
	 * @param useSharedHeader_ Print the stored header (from an input sam file) in all output sam files 
	 * @param defaultFormat_ Assume files are in this format if they don't have a valid extension
	 * @param threaded_ Run this mcros in its own thread
	 * @param maxStreams_ Max allowed number of concurrent open streams
	 */
	public BufferedMultiCros(String pattern1_, String pattern2_,
			boolean overwrite_, boolean append_, boolean allowSubprocess_, boolean useSharedHeader_, 
			int defaultFormat_, boolean threaded_, int maxStreams_){
		assert(pattern1_!=null && pattern1_.indexOf('%')>=0);
		assert(pattern2_==null || pattern1_.indexOf('%')>=0);
		
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

	/** True if no errors were encountered */
	public abstract boolean finishedSuccessfully();

	/** 
	 * Add a single read.  Should not be used in threaded mode.
	 * Should only be used by this class.
	 * @param r Read to add.
	 * @param name Name of destination buffer.
	 */
	abstract void add(Read r, String name);
	
	/** 
	 * Dump all buffered reads to disk, except when minReadsToDump forbids it.
	 * @return Number of reads dumped.
	 */
	abstract long dumpAll();
	
	/** 
	 * Dump all residual reads to this stream.
	 * @param rosu Destination stream.
	 * @return Number of residual reads dumped.
	 */
	public abstract long dumpResidual(ConcurrentReadOutputStream rosu);
	
	/** Dump everything and close any open streams. */
	abstract long closeInner();
	
	/** Generate a report on how many reads went to each file */
	public abstract ByteBuilder report();
	
	/** Time for shutting down output threads */
	public String printRetireTime() {
		throw new RuntimeException("printRetireTime not available for "+getClass().getName());
	}
	
	/** Time for shutting down output threads */
	public String printCreateTime() {
		throw new RuntimeException("printRetireTime not available for "+getClass().getName());
	}
	
	public abstract Set<String> getKeys();
	
	/*--------------------------------------------------------------*/
	/*----------------        Final Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Shut this down and perform any cleanup needed. */
	public final void close(){
		if(threaded){poisonAndWait();}
		else{closeInner();}
	}
	
	/** Primary file pattern */
	public final String fname(){return pattern1;}
	
	/** Return true if this stream has detected an error */
	public final boolean errorState(){
		return errorState;
	}
	
	/** 
	 * Send a list of reads to an output buffer.
	 * The reads must have a name attached to the object field in order to be written. 
	 */
	public final void add(ArrayList<Read> list) {
		if(threaded){//Send to the transfer queue
			addToQueue(list);
		}else{//Add the reads from this thread
			addToBuffers(list);
		}
	}
	
	/** Send individual reads to their designated buffer */
	private final void addToBuffers(ArrayList<Read> list) {
		for(Read r : list){
			if(r.obj!=null){
				String name=(String)r.obj;
				readsInTotal++;
				add(r, name);//Reads without a name in the obj field get ignored here.
			}
		}
		handleLoad0();
	}
	
	/** Called to handle load after adding a list */
	void handleLoad0() {
		//Do nothing
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Threaded Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	/** For threaded mode */
	public final void run(){
		assert(threaded) : "This should only be called in threaded mode.";
		try {
			for(ArrayList<Read> list=transferQueue.take(); list!=poisonToken; list=transferQueue.take()){
				if(verbose){System.err.println("Got list; size=\"+transferQueue.size())");}
				addToBuffers(list);
				if(verbose){System.err.println("Added list; size="+transferQueue.size());}
			}
		} catch (InterruptedException e) {
			//Terminate JVM if something goes wrong
			KillSwitch.exceptionKill(e);
		}
		closeInner();
	}
	
	/** Indicate that no more reads will be sent, for threaded mode */
	public final void poison(){
		assert(threaded) : "This should only be called in threaded mode.";
		addToQueue(poisonToken);
	}
	
	boolean addToQueue(ArrayList<Read> list) {
		boolean success=false;
		for(int i=0; i<10 && !success; i++) {
			try {
				transferQueue.put(list);
				success=true;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		if(!success) {
			KillSwitch.kill("Something went wrong when adding to "+getClass().getName());
		}
		return success;
	}
	
	/** Indicate that no more reads will be sent, for threaded mode */
	public final void poisonAndWait(){
		assert(threaded) : "This should only be called in threaded mode.";
		poison();
		waitForFinish();
	}
	
	/** Wait for this object's thread to terminate */
	public final void waitForFinish(){
		assert(threaded);
		if(verbose){System.err.println("Waiting for finish.");}
		while(this.getState()!=Thread.State.TERMINATED){
			if(verbose){System.err.println("Attempting join.");}
			try {
				this.join(1000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output file patterns containing a % symbol */
	public final String pattern1, pattern2;
	
	/** True if an error was encountered */
	boolean errorState=false;
	
	/** File overwrite permission */
	final boolean overwrite;
	
	/** File append permission */
	final boolean append;
	
	/** Subprocess spawning permission (e.g., for pigz) */
	final boolean allowSubprocess;
	
	/** Output file format, if unclear from file extension */
	final int defaultFormat;
	
	/** Buffers for each ReadStreamWriter */
	int rswBuffers=1;
	
	/** Print the shared header (for sam files) */
	final boolean useSharedHeader;
	
	/** Don't retire below this limit */
	final long memLimitLower;
	
	/** Possibly take some action */
	final long memLimitMid;
	
	/** Dump everything if this limit is reached from buffered reads */
	final long memLimitUpper;

	/** Allow this many active streams, for MCros3+ */
	public final int maxStreams;
	
	/** Retire this many streams at a time */
	public final int streamsToRetire;
	
	/** Dump a buffer once it holds this many reads */
	public int readsPerBuffer=defaultReadsPerBuffer;
	
	/** Dump a buffer once it holds this many bytes (estimated) */
	public int bytesPerBuffer=defaultBytesPerBuffer;
	
	/** Never write files with fewer than this many reads */
	public long minReadsToDump=0;

	/** Number of reads encountered that were not written */
	public long residualReads=0, residualBases=0;
	
	long readsInTotal=0;
	
	/** Current number of buffered reads */
	long readsInFlight=0;
	
	/** Current number of buffered bytes (estimated) */
	long bytesInFlight=0;
	
	/** Used when MultiCros is run in threaded mode */
	private final ArrayBlockingQueue<ArrayList<Read>> transferQueue;
	
	/** Signal to terminate when in threaded mode */
	private final ArrayList<Read> poisonToken=new ArrayList<Read>(0);
	
	/** True if this object is intended to run in a separate thread */
	public final boolean threaded;
	
	/** Use a LogLog to track cardinality for each output file */
	public boolean trackCardinality=false;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private static float memLimitLowerMult=0.20f;
	private static float memLimitMidMult=0.40f;
	private static float memLimitUpperMult=0.60f;
	public static boolean defaultThreaded=true;
	public static int defaultMaxStreams=12;
	public static int defaultMcrosType=6;
	public static int defaultReadsPerBuffer=32000;
	public static int defaultBytesPerBuffer=16000000;
	
	public static boolean verbose=false;

}
