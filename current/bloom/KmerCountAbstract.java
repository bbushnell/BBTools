package bloom;

/**
 * @author Brian Bushnell
 * @date Dec 2, 2014
 *
 */
public abstract class KmerCountAbstract {

	/**
	 * Transforms a count array into a frequency histogram.
	 * Creates a histogram showing how many k-mers occur at each count level,
	 * with counts capped at maximum histogram size.
	 *
	 * @param count Array of k-mer counts
	 * @return Frequency histogram where index represents count value and value represents frequency
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

	/**
	 * Calculates the sum of all elements in an integer array.
	 * @param array Array to sum
	 * @return Total sum of all array elements
	 */
	protected static final long sum(int[] array){
		long x=0;
		for(int y : array){x+=y;}
		return x;
	}

	/**
	 * Calculates the sum of all elements in a long array.
	 * @param array Array to sum
	 * @return Total sum of all array elements
	 */
	protected static final long sum(long[] array){
		long x=0;
		for(long y : array){x+=y;}
		return x;
	}

	protected static final int min(int x, int y){return x<y ? x : y;}
	protected static final int max(int x, int y){return x>y ? x : y;}
	
	/** Minimum quality score for considering bases during k-mer counting */
	public static byte minQuality=6;
	/** Total number of reads processed during k-mer counting operations */
	public static long readsProcessed=0;
	/** Maximum number of reads to process, or -1 for unlimited */
	public static long maxReads=-1;
	
	/** Minimum probability threshold for k-mer operations */
	public static float minProb=0.5f;
	
	/** Total number of k-mers counted during processing */
	public static long keysCounted=0;
	/** Total number of increment operations performed on k-mer counts */
	public static long increments=0;
	
	/** Flag for verbose output during k-mer counting operations */
	public static final boolean verbose=false;
	/** Whether to pre-join overlapping read pairs before k-mer counting */
	public static boolean PREJOIN=false;
	/**
	 * Whether to use canonical k-mer representation for strand-independent counting
	 */
	public static boolean CANONICAL=false;
	/** Whether to count duplicate k-mers within a single read multiple times */
	public static boolean KEEP_DUPLICATE_KMERS=true;
	/** Maximum number of k-mers to count per read, or 0 for unlimited */
	public static int KMERS_PER_READ=0;
	/** Bit mask for read ID processing during k-mer counting */
	public static int IDMASK=0;
	/** Maximum number of threads to use for k-mer counting operations */
	public static int MAX_COUNT_THREADS=1024;
	/** Whether to operate in sketch mode for approximate k-mer counting */
	public static boolean SKETCH_MODE=false;
	/** Whether to store hashed k-mer values instead of original k-mers */
	public static boolean STORE_HASHED=false;
	/** Whether to use buffered I/O for improved performance during counting */
	public static boolean BUFFERED=false;
	/** Buffer length for buffered I/O operations in k-mer counting */
	public static int BUFFERLEN=3000; //Optimal is in the range of 2000-8000 for Clumpified 2x150bp data.
//	public static boolean SORT_SERIAL=false; //Not needed, see parallel sort flag
	
	/** Maximum k-mer length for direct encoding without hashing */
	public static int maxShortKmerLength=31;
	
}
