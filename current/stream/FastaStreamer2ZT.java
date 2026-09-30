package stream;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import fileIO.ByteFile1Fc;
import fileIO.FileFormat;
import shared.KillSwitch;
import shared.Shared;
import structures.IntList;
import structures.ListNum;

/**
 * Read names honor Shared.TRIM_READ_DESCRIPTION using legacy byte-whitespace semantics.
 * Synchronous FASTA reader using ByteFile1Fc record blocks; no worker thread.
 * Sequence lines are joined by the backend. Returned reads have copied bases,
 * US-ASCII names, null qualities and no mates. Interleaved input is unsupported.
 * Limits, counters and numeric IDs count retained records after sampling.
 * Start once before reading, and close in a caller finally block if iteration
 * exits early or throws. Synchronized operations serialize on this instance;
 * unsynchronized status getters are not cross-thread completion barriers.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date November 12, 2025
 */
public class FastaStreamer2ZT implements Streamer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Configures a FASTA input path without opening it.
	 * @param fname_ Input path; subprocess input is allowed by the descriptor
	 * @param pairnum_ Pair marker (0 or 1), without linking mates
	 * @param maxReads_ Maximum retained records; negative means unlimited
	 */
	public FastaStreamer2ZT(String fname_, int pairnum_, long maxReads_){
		this(FileFormat.testInput(fname_, FileFormat.FASTA, null, true, false), pairnum_, maxReads_);
	}

	/**
	 * Configures a noninterleaved input descriptor without opening it.
	 * The initial read flags reflect descriptor/global amino settings.
	 * @param ffin_ Nonnull input descriptor
	 * @param pairnum_ Pair marker (0 or 1), checked by assertion
	 * @param maxReads_ Maximum retained records; negative means unlimited
	 */
	public FastaStreamer2ZT(FileFormat ffin_, int pairnum_, long maxReads_){
		ffin=ffin_;
		fname=ffin_.name();
		pairnum=pairnum_;
		flag=(ffin.amino() || Shared.AMINO_IN ? Read.AAMASK : 0);
		assert(pairnum==0 || pairnum==1) : pairnum;
		interleaved=(ffin.interleaved());
		assert(!interleaved) : "FastaStreamer2ZT does not support interleaved files";
		maxReads=(maxReads_<0 ? Long.MAX_VALUE : maxReads_);

		if(verbose){outstream.println("Made FastaStreamer2ZT");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Opens ByteFile1Fc and resets counters. Call once: this is not a restart
	 * operation and does not reset terminal/list/error state or close an old input.
	 */
	@Override
	public synchronized void start(){
		if(verbose){outstream.println("FastaStreamer2ZT.start() called.");}

		//Reset counters
		readsProcessed=0;
		basesProcessed=0;

		//Open the file
		bf=new ByteFile1Fc(ffin);

		if(verbose){outstream.println("FastaStreamer2ZT started.");}
	}

	/**
	 * Closes the current backend and folds its reported error flag.
	 * Clears the backend reference only after normal return; exceptions propagate.
	 * Does not itself mark iteration finished. Repeated successful close is safe.
	 */
	@Override
	public synchronized void close(){
		if(bf!=null){
			errorState|=bf.close();//Preserve the backend's reported read/close errors.
			bf=null;
		}
	}

	/** Returns the configured input path. */
	@Override
	public String fname(){return fname;}

	/** Returns an unsynchronized hint; false only after nextList observes termination. */
	@Override
	public boolean hasMore(){return !finished;}

	/** Returns the cached backend error flag; thrown exceptions are not universally latched. */
	@Override
	public boolean errorState(){return errorState;}

	/** Returns false; this reader never links mates. */
	@Override
	public boolean paired(){return false;}

	/** Returns the pair marker assigned to each retained read. */
	@Override
	public int pairnum(){return pairnum;}

	/**
	 * Returns the number of retained reads constructed and added to batches.
	 * A later parsing exception can prevent return of an already-counted batch.
	 */
	@Override
	public synchronized long readsProcessed(){return readsProcessed;}

	/** Returns the summed lengths of counted retained reads; see readsProcessed(). */
	@Override
	public synchronized long basesProcessed(){return basesProcessed;}

	/**
	 * Sets the retention rate and resets the generator for subsequent records.
	 * Header/base slices are copied before sampling; only retained records become Reads.
	 * @param rate Retention probability; at least 1 disables sampling, 0 drops all
	 * @param seed Nonnegative for reproducibility; negative selects a time-based seed
	 */
	@Override
	public synchronized void setSampleRate(float rate, long seed){
		samplerate=rate;
		randy=(rate>=1f ? null : shared.Shared.random(seed));
	}

	/**
	 * Returns the next nonempty retained batch, or null after EOF/limit/closed input.
	 * Fully discarded blocks are skipped. Batch sizes follow backend blocks, with
	 * contiguous batch IDs and firstRecordNum equal to the preceding retained count.
	 * Normal termination closes input. Parsing/Read-construction exceptions propagate.
	 * Returned lists and copied record data remain owned by the caller.
	 */
	@Override
	public synchronized ListNum<Read> nextList(){
		if(finished){return null;}

		ListNum<Read> list=nextListSingle();

		if(list==null || list.size()==0){
			finished=true;
			close();
			return null;
		}

		listNum++;
		return list;
	}

	/** Always throws UnsupportedOperationException; FASTA has no SamLine representation. */
	@Override
	public ListNum<SamLine> nextLines(){
		throw new UnsupportedOperationException("FASTA does not support SamLine");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses blocks until a retained batch is available or input/limit is exhausted.
	 * Uses recorded newline boundaries, not the capacity of the backend byte array.
	 * Constructs Reads with current global constructor validation settings; no
	 * additional explicit validation is performed here. Caller holds this monitor.
	 */
	private ListNum<Read> nextListSingle(){
		if(bf==null){return null;}

		ArrayList<Read> readList=new ArrayList<Read>();
		ListNum<Read> reads=new ListNum<Read>(readList, listNum);
		reads.firstRecordNum=readsProcessed;

		//[stream/FastaStreamer2ZT#001] Skip fully discarded blocks until data, EOF or the retained limit.
		//A block that is fully subsampled-out must NOT be returned empty: nextList() treats an empty
		//list as EOF, so a single unsampled block (1-record blocks for large contigs!) used to silently
		//truncate the rest of the file under samplerate<1. Limit exhaustion also ends iteration.
		while(readList.isEmpty() && readsProcessed<maxReads){
			// Get block of records with newline positions
			byte[] block=bf.nextLine(newlines);
			if(block==null || block.length==0){break;}

			for(int i=0, nl0=-1; i<newlines.size() && readsProcessed<maxReads; i++){
				int nl1=newlines.get(i);
				int nl2=(newlines.size()>i+1 ? newlines.get(i+1) : nl1);
				assert(block[nl0+1]=='>') : nl0+", "+(char)block[nl0+1];
				final byte[] header=KillSwitch.copyOfRange(block, nl0+2, nl1);
				final byte[] bases=(nl2>nl1 ? KillSwitch.copyOfRange(block, nl1+1, nl2) : null);
				if(samplerate>=1f || randy.nextFloat()<samplerate){
					Read r=new Read(bases, null, new String(header, 0, ReadHeader.end(header, 0, Shared.TRIM_READ_DESCRIPTION), StandardCharsets.US_ASCII), readsProcessed, flag);
					r.setPairnum(pairnum);
					readList.add(r);
					readsProcessed++;
					basesProcessed+=r.length();
				}
				if(bases!=null){
					i++;
					nl0=nl2;
				}else{nl0=nl1;}
			}
		}

		return readList.isEmpty() ? null : reads;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path */
	public final String fname;

	/** Primary input file */
	final FileFormat ffin;

	/** Record-block backend; null before start and after successful close. */
	private ByteFile1Fc bf;

	/** Pair marker applied to retained reads. */
	final int pairnum;
	/** Descriptor setting, asserted false during construction. */
	final boolean interleaved;
	/** Reused boundaries for the meaningful region of each backend block. */
	private final IntList newlines=new IntList(256);

	/** Retained reads added to batches, including a batch whose later parsing throws. */
	protected long readsProcessed=0;
	/** Total length of counted retained reads. */
	protected long basesProcessed=0;

	/** Maximum retained records, normalized to Long.MAX_VALUE for a negative limit. */
	final long maxReads;
	/** Read-constructor flags; configure before consuming input. */
	public int flag;

	/** ID for the next returned nonempty list. */
	private long listNum=0;

	/** True after nextList observes normal termination; close alone does not set it. */
	private boolean finished=false;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Print status messages to this output stream */
	protected PrintStream outstream=System.err;
	/** Print verbose messages */
	public static final boolean verbose=false;
	/** Cached backend error reports, folded on close; not a catch-all exception flag. */
	public boolean errorState=false;
	/** Retention threshold; at least 1 bypasses the generator. */
	private float samplerate=1f;
	/** Sampling generator, null when sampling is disabled. */
	private shared.Random randy=null;

}