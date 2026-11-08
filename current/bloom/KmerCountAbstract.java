package bloom;

/**
 * Shared counters, configuration and histogram helpers for the counting implementations.
 * Configuration is process-wide, not isolated between instances. Set it before starting
 * workers; callers coordinating simultaneous counting operations must also coordinate
 * these fields. Individual implementations support different subsets of the options.
 *
 * @author Brian Bushnell
 * @date Dec 2, 2014
 */
public abstract class KmerCountAbstract{

	/*--------------------------------------------------------------*/
	/*----------------       Static Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Counts how many nonnegative entries have each value, without changing the input.
	 * Values at least 1999 share the final bin.
	 * @param count Non-null array of nonnegative counts
	 * @return Newly allocated 2000-bin frequency histogram
	 */
	protected static final long[] transformToFrequency(int[] count){
		long[] freq=new long[2000];
		int max=freq.length-1;
		for(int i=0; i<count.length; i++){
			int x=count[i];
			x=min(x, max);
			freq[x]++;
		}
		return freq;
	}

	/** Sums a non-null integer array using a long accumulator. */
	protected static final long sum(int[] array){
		long x=0;
		for(int y : array){x+=y;}
		return x;
	}

	/** Sums a non-null long array with ordinary signed-long overflow semantics. */
	protected static final long sum(long[] array){
		long x=0;
		for(long y : array){x+=y;}
		return x;
	}

	/** Returns the smaller value. */
	protected static final int min(int x, int y){return x<y ? x : y;}
	/** Returns the larger value. */
	protected static final int max(int x, int y){return x>y ? x : y;}

	/*--------------------------------------------------------------*/
	/*----------------        Shared State          ----------------*/
	/*--------------------------------------------------------------*/

	/** Minimum Phred score used by quality-aware counting paths. */
	public static byte minQuality=6;
	/** Accumulated reads: ReadCounter counts mates separately; older counters count roots. */
	public static long readsProcessed=0;
	/** Shared reader-dependent input limit: FASTQ roots/pairs or SAM records; negative is unlimited. */
	public static long maxReads=-1;

	/** Minimum k-mer correctness probability used by probability-aware counting paths. */
	public static float minProb=0.5f;

	/** Accepted occurrences submitted for counting, including buffered multiplicities. */
	public static long keysCounted=0;
	/** Array-update calls issued; not the number of cells whose values changed. */
	public static long increments=0;

	/** Compile-time diagnostic switch for subclasses using the inherited flag. */
	public static final boolean verbose=false;
	/** Join reads with a known positive insert before supported short-k counting paths. */
	public static boolean PREJOIN=false;
	/** Shared canonical-mode contract checked by supported short-k counting paths. */
	public static boolean CANONICAL=false;
	/** Retain repeated keys within a read pair in ReadCounter's short-k path. */
	public static boolean KEEP_DUPLICATE_KMERS=true;
	/** Positive values select up to this many sampled windows per read in that path. */
	public static int KMERS_PER_READ=0;
	/** Mask for index-based subsampling relative to read numericID; zero disables it. */
	public static int IDMASK=0;
	/** Upper worker limit used by ReadCounter's ordinary input-counting overload. */
	public static int MAX_COUNT_THREADS=1024;
	/** Apply the sketch hash cutoff before compatible single-key array updates. */
	public static boolean SKETCH_MODE=false;
	/** Store the sketch hash instead of the original key when SKETCH_MODE is enabled. */
	public static boolean STORE_HASHED=false;
	/** Sort and combine compatible single-key updates in per-worker buffers. */
	public static boolean BUFFERED=false;
	/** Buffer flush threshold; independent two-hash updates use their direct path. */
	public static int BUFFERLEN=3000; //Optimal is in the range of 2000-8000 for Clumpified 2x150bp data.
	//SORT_SERIAL is obsolete; use the shared parallel-sort setting.

	/** Largest k routed through ReadCounter's primitive-long nucleotide path. */
	public static int maxShortKmerLength=31;

}
