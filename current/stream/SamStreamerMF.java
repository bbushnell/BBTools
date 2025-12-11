package stream;

import java.io.PrintStream;
import java.util.ArrayDeque;

import fileIO.FileFormat;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ListNum;

/**
 * Multi-file SAM stream processor that loads multiple SAM files rapidly using multiple threads.
 * Manages concurrent reading from multiple SAM files with dynamic load balancing,
 * allowing efficient processing of large numbers of SAM files simultaneously.
 *
 * @author Brian Bushnell
 * @date March 6, 2019
 */
public class SamStreamerMF {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Program entry point for command-line execution.
	 * Creates a SamStreamerMF instance and runs a test processing loop.
	 * @param args Command-line arguments: [comma-separated file paths] [optional thread count]
	 */
	public static final void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		int threads=Shared.threads();
		if(args.length>1){threads=Integer.parseInt(args[1]);}
		SamStreamerMF x=new SamStreamerMF(args[0].split(","), threads, false, -1);
		
		//Run the object
		x.start();
		x.test();
		
		t.stop("Time: ");
	}
	
	/**
	 * Constructs a SamStreamerMF from file names.
	 * Converts file names to FileFormat objects and delegates to primary constructor.
	 *
	 * @param fnames_ Array of SAM file names to process
	 * @param threads_ Number of threads to use for processing
	 * @param saveHeader_ Whether to save and propagate SAM headers
	 * @param maxReads_ Maximum number of reads to process, or -1 for no limit
	 */
	public SamStreamerMF(String[] fnames_, int threads_, boolean saveHeader_, long maxReads_){
		this(FileFormat.testInput(fnames_, FileFormat.SAM, null, true, false), threads_, saveHeader_, maxReads_);
	}
	
	/**
	 * Primary constructor for SamStreamerMF.
	 * Initializes the multi-file SAM streaming system with specified parameters.
	 *
	 * @param ffin_ Array of FileFormat objects representing SAM files to process
	 * @param threads_ Number of threads to use for processing
	 * @param saveHeader_ Whether to save and propagate SAM headers
	 * @param maxReads_ Maximum number of reads to process, or -1 for no limit
	 */
	public SamStreamerMF(FileFormat[] ffin_, int threads_, boolean saveHeader_, long maxReads_){
		fname=ffin_[0].name();
		threads=threads_;
		ffin=ffin_;
		saveHeader=saveHeader_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	
	final void test(){
		for(ListNum<Read> list=nextReads(); list!=null; list=nextReads()){
			if(verbose){outstream.println("Got list of size "+list.size());}
		}
	}
	
	
	public final void start(){
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the reads in separate threads
		spawnThreads();
		
		if(verbose){outstream.println("Finished; closing streams.");}
	}

	public final ListNum<Read> nextList(){return nextReads();}
	public final ListNum<Read> nextReads(){
		ListNum<Read> list=null;
		assert(activeStreamers!=null);
		synchronized(activeStreamers){
			if(activeStreamers.isEmpty()){return null;}
			while(list==null && !activeStreamers.isEmpty()){
				Streamer srs=activeStreamers.poll();
				list=srs.nextList();
				if(list!=null){activeStreamers.add(srs);}
				else{
					readsProcessed+=srs.readsProcessed();
					basesProcessed+=srs.basesProcessed();
//					if(srs.header!=null){//Should be automatic now
//						SamReadInputStream.setSharedHeader(srs.header);
//					}
					
					if(!streamerSource.isEmpty()){
						srs=streamerSource.poll();
						srs.start();
						activeStreamers.add(srs);
					}
				}
			}
		}
		return list;
	}
//	public abstract ListNum<SamLine> nextLines();
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Initializes and starts the thread management system for concurrent SAM processing.
	 * Creates streamer queues, calculates optimal active streamer count, and starts initial threads. */
	void spawnThreads(){
		final int maxActive=Tools.max(2, Tools.min((Shared.threads()+4)/5, ffin.length, MAX_FILES));
		streamerSource=new ArrayDeque<Streamer>(ffin.length);
		activeStreamers=new ArrayDeque<Streamer>(maxActive);
		for(int i=0; i<ffin.length; i++){
			Streamer srs=StreamerFactory.makeSamOrBamStreamer(ffin[i], threads, saveHeader && i==0, false, maxReads, true);
			streamerSource.add(srs);
		}
		while(activeStreamers.size()<maxActive && !streamerSource.isEmpty()){
			Streamer srs=streamerSource.poll();
			srs.start();
			activeStreamers.add(srs);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	protected String fname;
	
	/*--------------------------------------------------------------*/

	protected long readsProcessed=0;
	protected long basesProcessed=0;

	/** Maximum reads to process before stopping; -1 means no limit */
	protected long maxReads=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	final boolean saveHeader;

	/** Array of input file formats representing SAM files to process */
	final FileFormat[] ffin;
	
//	final Streamer[] streamers;
	/** Queue of unstarted streamers waiting to be activated */
	private ArrayDeque<Streamer> streamerSource;
	private ArrayDeque<Streamer> activeStreamers;
	
	final int threads;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int DEFAULT_THREADS=6;
	public static int MAX_FILES=8;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages */
	protected PrintStream outstream=System.err;
	public static final boolean verbose=false;
	public static final boolean verbose2=false;
	/** Indicates whether an error was encountered during processing */
	public boolean errorState=false;
	
}
