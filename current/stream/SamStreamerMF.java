package stream;

import java.io.PrintStream;
import java.util.ArrayDeque;

import fileIO.FileFormat;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ListNum;

/**
 * Multiplexes batches from several SAM/BAM readers in blocking round-robin order.
 * A selected slow child can delay the next call; batches and their IDs pass through
 * unchanged, without global renumbering or a merged file-order guarantee.
 * Construct with nonempty inputs, start once, drain to null, then inspect errorState.
 * There is no close or restart API. Child selection and header publication depend
 * on shared settings; configure them before starting and avoid concurrent changes.
 *
 * @author Brian Bushnell, Shinobu
 * @date March 6, 2019
 */
public class SamStreamerMF{

	/**
	 * Program entry point for command-line execution.
	 * Creates a SamStreamerMF instance and runs a test processing loop.
	 * Aborts the JVM on a reported input error before printing completion time.
	 * @param args Command-line arguments: [comma-separated file paths] [optional thread count]
	 */
	public static final void main(final String[] args){
		//Start a timer immediately upon code entrance.
		final Timer t=new Timer();

		//Create an instance of this class
		int threads=Shared.threads();
		if(args.length>1){threads=Integer.parseInt(args[1]);}
		final SamStreamerMF x=new SamStreamerMF(args[0].split(","), threads, false, -1);

		//Run the object
		x.start();
		x.test();

		//STR-029: nextReads folds child failures into the aggregate flag; reject it before reporting completion.
		if(x.errorState){KillSwitch.kill("Error reading SAM/BAM files: "+args[0]);}
		t.stop("Time: ");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Resolves descriptors with SAM as the fallback format, without starting readers.
	 *
	 * @param fnames_ Nonempty array of input paths; subprocess input is permitted
	 * @param threads_ Per-child reader-selection hint, not a total thread budget
	 * @param saveHeader_ Request shared-header publication from the first file only
	 * @param maxReads_ Limit forwarded separately to each selected reader, in its units; negative means unlimited
	 */
	public SamStreamerMF(final String[] fnames_, final int threads_, final boolean saveHeader_, final long maxReads_){
		this(FileFormat.testInput(fnames_, FileFormat.SAM, null, true, false), threads_, saveHeader_, maxReads_);
	}

	/**
	 * Retains the input descriptors and configuration without starting readers.
	 * The array is not copied; do not modify it while this instance is in use.
	 *
	 * @param ffin_ Nonempty array of nonnull descriptors accepted by StreamerFactory
	 * @param threads_ Per-child reader-selection hint, not a total thread budget
	 * @param saveHeader_ Request shared-header publication from the first file only
	 * @param maxReads_ Limit forwarded separately to each selected reader, in its units; negative means unlimited
	 */
	public SamStreamerMF(final FileFormat[] ffin_, final int threads_, final boolean saveHeader_, final long maxReads_){
		fname=ffin_[0].name();
		threads=threads_;
		ffin=ffin_;
		saveHeader=saveHeader_;
		maxReads=maxReads_;//STR-028: retain the requested per-child limit for spawnThreads.
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Drains all batches for the demo; main checks the accumulated error flag. */
	final void test(){
		for(ListNum<Read> list=nextReads(); list!=null; list=nextReads()){
			if(verbose){outstream.println("Got list of size "+list.size());}
		}
	}

	/** Resets counters and starts the initial children. Call once before consumption;
	 * the error flag is not reset, and this method does not wait for completion.
	 */
	public final void start(){
		//Reset counters
		readsProcessed=0;
		basesProcessed=0;

		//Process the reads in separate threads
		spawnThreads();

		if(verbose){outstream.println("Finished; closing streams.");}
	}

	/** Alias for nextReads with the same blocking and terminal semantics.
	 * @return Unmodified child batch, or null after all children have been retired
	 */
	public final ListNum<Read> nextList(){return nextReads();}

	/** Reads the next active child while holding the active-queue monitor.
	 * Nonnull batches, including empty ones, return unchanged and requeue the child.
	 * A null child result is treated as exhaustion: fold its currently reported counters
	 * and error flag, then activate the next queued reader. Retirement does not close
	 * or join the child. Library callers must inspect errorState after draining; this
	 * does not guarantee worker completion after interruption. Batch and Read numeric
	 * IDs remain child-local.
	 * @return Child batch, or null when no active children remain
	 */
	public final ListNum<Read> nextReads(){
		ListNum<Read> list=null;
		assert(activeStreamers!=null);
		synchronized(activeStreamers){
			if(activeStreamers.isEmpty()){return null;}
			while(list==null && !activeStreamers.isEmpty()){
				Streamer srs=activeStreamers.poll();
				list=srs.nextList();
				if(list!=null){activeStreamers.add(srs);}else{
					readsProcessed+=srs.readsProcessed();
					basesProcessed+=srs.basesProcessed();
					errorState|=srs.errorState();//[stream/SamStreamerMF#002] fold the exhausted sub-streamer's
					//error state (truncated/corrupt SAM among many) so it isn't silently dropped - the 2b reader-fold.

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

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs all children and starts up to the activation target.
	 * The target uses Shared.threads(), the file count and MAX_FILES, with a floor
	 * of two; actual starts also require available files. Each child receives threads
	 * as its factory hint, maxReads unchanged, an ordered=false request and Read
	 * conversion enabled. The selected backend determines ordering behavior.
	 * Only file zero receives the shared-header request.
	 */
	void spawnThreads(){
		final int maxActive=Tools.max(2, Tools.min((Shared.threads()+4)/5, ffin.length, MAX_FILES));
		streamerSource=new ArrayDeque<Streamer>(ffin.length);
		activeStreamers=new ArrayDeque<Streamer>(maxActive);
		for(int i=0; i<ffin.length; i++){
			final Streamer srs=StreamerFactory.makeSamOrBamStreamer(ffin[i], threads, saveHeader && i==0, false, maxReads, true);
			streamerSource.add(srs);
		}
		while(activeStreamers.size()<maxActive && !streamerSource.isEmpty()){
			final Streamer srs=streamerSource.poll();
			srs.start();
			activeStreamers.add(srs);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Name of the first input, retained for subclasses. */
	protected String fname;

	/** Sum of child read counts observed at retirement after a null batch result. */
	protected long readsProcessed=0;
	/** Sum of child base counts observed at retirement after a null batch result. */
	protected long basesProcessed=0;

	/** Limit passed separately to each child in its own units; negative means unlimited. */
	protected long maxReads=-1;

	/** Whether to request shared-header publication from file zero. */
	final boolean saveHeader;

	/** Caller-owned descriptor array; retained without copying. */
	final FileFormat[] ffin;

	/** Queue of unstarted children in input-file order. */
	private ArrayDeque<Streamer> streamerSource;
	/** Round-robin queue; its monitor serializes batch retrieval and child retirement. */
	private ArrayDeque<Streamer> activeStreamers;

	/** Factory thread hint supplied separately to every child. */
	final int threads;

	/** Destination for optional status messages. */
	protected PrintStream outstream=System.err;
	/** Aggregate errors observed at child retirement; inspect after draining. */
	public boolean errorState=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Legacy setting; this class does not use it to select a default thread count. */
	public static int DEFAULT_THREADS=6;
	/** Input to the activation target, which has a floor of two rather than a strict cap. */
	public static int MAX_FILES=8;
	/** Enables optional status messages when recompiled with true. */
	public static final boolean verbose=false;
	/** Reserved verbosity switch; currently unused. */
	public static final boolean verbose2=false;

}
