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
 * Abstract base class for buffered multi-stream output of reads to different files.
 * Each output stream is controlled by a buffer that stores reads until sufficient quantity accumulates to dump efficiently.
 * Supports threaded and non-threaded modes.
 * @author Brian Bushnell
 * @date May 14, 2019
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
	
	public BufferedMultiCros(String pattern1_, String pattern2_,
			boolean overwrite_, boolean append_, boolean allowSubprocess_, boolean useSharedHeader_, 
			int defaultFormat_, boolean threaded_, int maxStreams_){
		assert(pattern1_!=null && pattern1_.indexOf('%')>=0);
		assert(pattern2_==null || pattern1_.indexOf('%')>=0); //Possible bug: should check pattern2_.indexOf('%')>=0, not pattern1_
		
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

	public abstract boolean finishedSuccessfully();

	/**
	 * Adds a single read to the specified output buffer.
	 * Should not be used in threaded mode and should only be called by this class.
	 * @param r Read to add to buffer
	 * @param name Name of destination buffer
	 */
	abstract void add(Read r, String name);
	
	abstract long dumpAll();
	
	/**
	 * Dumps all residual reads to the specified output stream.
	 * @param rosu Destination stream for residual reads
	 * @return Number of residual reads that were dumped
	 */
	public abstract long dumpResidual(ConcurrentReadOutputStream rosu);
	
	/** Dumps everything and closes any open streams.
	 * @return Number of reads processed during shutdown */
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
	public String printRetireTime() {
		throw new RuntimeException("printRetireTime not available for "+getClass().getName());
	}
	
	/**
	 * Gets timing information for creating output threads.
	 * Default implementation throws RuntimeException for unsupported subclasses.
	 * @return Formatted timing information
	 * @throws RuntimeException if not implemented by subclass
	 */
	public String printCreateTime() {
		throw new RuntimeException("printRetireTime not available for "+getClass().getName());
	}
	
	public abstract Set<String> getKeys();
	
	/*--------------------------------------------------------------*/
	/*----------------        Final Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Shuts down the MultiCros and performs any needed cleanup.
	 * Handles both threaded and non-threaded shutdown procedures. */
	public final void close(){
		if(threaded){poisonAndWait();}
		else{closeInner();}
	}
	
	/** Gets the primary file pattern.
	 * @return The pattern1 file pattern string */
	public final String fname(){return pattern1;}
	
	/** Checks if this stream has detected an error condition.
	 * @return true if an error state has been detected */
	public final boolean errorState(){
		return errorState;
	}
	
	public final void add(ArrayList<Read> list) {
		if(threaded){//Send to the transfer queue
			addToQueue(list);
		}else{//Add the reads from this thread
			addToBuffers(list);
		}
	}
	
	/**
	 * Distributes individual reads to their designated buffers.
	 * Only processes reads that have a name in their obj field.
	 * @param list List of reads to add to buffers
	 */
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
	
	/** Called after adding a list of reads to handle load management.
	 * Default implementation does nothing; subclasses may override. */
	void handleLoad0() {
		//Do nothing
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Threaded Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Main execution method for threaded mode.
	 * Continuously processes read lists from the transfer queue until poisoned.
	 * Terminates JVM if interrupted during operation.
	 */
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
	
	/** Signals that no more reads will be sent in threaded mode.
	 * Adds the poison token to the transfer queue to terminate the thread. */
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
	
	/** Poisons the transfer queue and waits for thread termination.
	 * Combines poison() and waitForFinish() for clean shutdown. */
	public final void poisonAndWait(){
		assert(threaded) : "This should only be called in threaded mode.";
		poison();
		waitForFinish();
	}
	
	/** Waits for this object's thread to terminate.
	 * Repeatedly attempts to join with 1 second timeout until thread terminates. */
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
	
	public final String pattern1, pattern2;
	
	/** True if an error was encountered during processing */
	boolean errorState=false;
	
	final boolean overwrite;
	
	/** Permission to append to existing files */
	final boolean append;
	
	/** Permission to spawn subprocesses (e.g., for pigz compression) */
	final boolean allowSubprocess;
	
	/** Default output file format when unclear from file extension */
	final int defaultFormat;
	
	/** Number of buffers for each ReadStreamWriter */
	int rswBuffers=1;
	
	/** Whether to print the shared header in SAM output files */
	final boolean useSharedHeader;
	
	/** Lower memory threshold below which streams are not retired */
	final long memLimitLower;
	
	/** Middle memory threshold where some action may be taken */
	final long memLimitMid;
	
	/** Upper memory threshold where everything is dumped if reached */
	final long memLimitUpper;

	/** Maximum number of active streams allowed for MCros3+ implementations */
	public final int maxStreams;
	
	/** Number of streams to retire at a time for load balancing */
	public final int streamsToRetire;
	
	/** Trigger buffer dump when it contains this many reads */
	public int readsPerBuffer=defaultReadsPerBuffer;
	
	/** Trigger buffer dump when it contains this many estimated bytes */
	public int bytesPerBuffer=defaultBytesPerBuffer;
	
	public long minReadsToDump=0;

	public long residualReads=0, residualBases=0;
	
	long readsInTotal=0;
	
	/** Current number of reads held in buffers */
	long readsInFlight=0;
	
	/** Current number of estimated bytes held in buffers */
	long bytesInFlight=0;
	
	/** Queue for transferring read lists when MultiCros runs in threaded mode */
	private final ArrayBlockingQueue<ArrayList<Read>> transferQueue;
	
	/** Special empty list used to signal thread termination in threaded mode */
	private final ArrayList<Read> poisonToken=new ArrayList<Read>(0);
	
	public final boolean threaded;
	
	/** Whether to use LogLog for tracking cardinality of each output file */
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
