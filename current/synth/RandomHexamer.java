package synth;

import shared.Random;

import dna.AminoAcid;

/**
 * Assigns reusable, sequence-dependent keep probabilities to the first k bases of a read.
 * The default k is six, modeling unequal primer efficiencies during random hexamer priming.
 * Each encoded k-mer receives one probability when the shared table is built; subsequent
 * reads with that prefix reuse it and draw separately to decide whether to keep the read.
 *
 * Configure and initialize this shared static state before starting worker threads.
 * Do not change settings or rebuild the table while readers are using it. The construction
 * monitor does not make the unsynchronized setters, fast paths or keep calls safe for
 * arbitrary concurrent configuration. RandomReadsMG initializes before starting its workers.
 *
 * @author Brian Bushnell
 * @author Isla
 */
public class RandomHexamer{

	/**
	 * Samples the shared probability for the read's first k bases.
	 * Reads shorter than k and prefixes containing an unencoded ASCII base are kept without
	 * consuming a random draw. Otherwise, the encoded prefix selects a table entry and the
	 * read is kept when that value exceeds the next random float. Lowercase bases are accepted;
	 * U is encoded as T by AminoAcid.baseToNumber. Bases after the prefix are not examined.
	 *
	 * @param bases Non-null sequence whose examined prefix contains ASCII base symbols
	 * @param randy Caller-owned generator, used only for a fully encoded prefix of length k
	 * @return true for a bypassed read or a successful probability draw
	 * @throws AssertionError if the table is not initialized and assertions are enabled
	 */
	public static boolean keep(byte[] bases, Random randy){
		assert(initialized) : "RandomHexamer must be initialized prior to use.";
		if(bases.length<k){return true;}
		//Generate a 2-bit number from the first k bases
		int kmer=0;
		for(int i=0; i<k; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];
			if(x<0){return true;} // Skip ambiguous sequences
			kmer=(kmer<<2) | AminoAcid.baseToNumber[b];
		}
//		assert(false) : bases.length+", "+k+", "+probs[kmer]+", "+Arrays.toString(probs);
		return probs[kmer]>randy.nextFloat(); // Higher prob = better priming = keep read
	}
	
	/**
	 * Applies the k, exponent and minimum-probability settings in that order.
	 * A changed value invalidates the current table; call initialize before using it again.
	 * This method does not build a table or provide atomic configuration for concurrent callers.
	 *
	 * @param k_ Prefix length in [1,15], default 6
	 * @param power_ Positive exponent applied to a uniform draw, default 0.5
	 * @param minProb_ Lower endpoint of the keep-probability transform in [0,1), default 0.1
	 */
	public static void set(int k_, float power_, float minProb_){
		setK(k_);
		setPower(power_);
		setMinProb(minProb_);
	}
	
	/**
	 * Sets the lower endpoint of the probability transform, invalidating the table if changed.
	 * The transformed draw is scaled and shifted into the requested interval, not clamped
	 * after sampling. A positive minimum prevents a zero keep probability for any encoded k-mer.
	 * @param minProb_ Lower endpoint in [0,1), default 0.1
	 */
	public static void setMinProb(float minProb_){
		initialized=initialized && (minProb==minProb_);
		minProb=minProb_;
		assert(minProb>=0 && minProb<1) : minProb_;
	}
	
	/**
	 * Sets the exponent of the probability transform, invalidating the table if changed.
	 * For a fixed draw strictly between zero and one, a larger exponent lowers its keep
	 * probability toward minProb; a smaller exponent raises it toward one.
	 * @param power_ Positive exponent, default 0.5
	 */
	public static void setPower(float power_){
		initialized=initialized && (power==power_);
		power=power_;
		assert(power>0) : power;
	}
	
	/**
	 * Sets the encoded prefix length, invalidating the table if changed.
	 * Initialization allocates one float per possible k-mer; the encoding range does not
	 * guarantee that the corresponding table fits in the available heap.
	 * @param k_ Prefix length in [1,15], default 6
	 */
	public static void setK(int k_){
		initialized=initialized && (k==k_);
		k=k_;
		assert(k>0 && k<=15) : k;
	}
	
	/** Returns the configured prefix length, whether or not the table is initialized. */
	public static int getK(){return k;}
	
	/**
	 * Initializes the probability table using the configured prefix length.
	 * @param randy Generator for table construction; unused when the current table is reused
	 * @return true after construction or reuse; failures propagate to the caller
	 */
	public static boolean initialize(Random randy){
		return initialize(randy, k);
	}
	
	/**
	 * Builds or reuses the shared table for the specified prefix length.
	 * If already initialized for this k, the table is reused without consulting randy;
	 * passing a new generator or seed alone does not regenerate probabilities. Otherwise,
	 * construction consumes one float draw per table entry and marks the table initialized
	 * after all entries are assigned. The inner check avoids rebuilding after another
	 * initialization has completed under the monitor; the class-level threading restrictions apply.
	 *
	 * This overload validates the same supported range as setK, but checks it before
	 * changing shared state. Callers must ensure adequate heap for the resulting table.
	 *
	 * @param randy Generator used for construction, or unused if the table is reused
	 * @param k_ Encoded prefix length in [1,15]
	 * @return true after construction or reuse; failures propagate to the caller
	 */
	public static boolean initialize(Random randy, int k_){
		assert(k_>0 && k_<=15) :
			"Prefix length must be in [1,15] for positive length and 2-bit int indexing: "+k_;
		if(initialized && k==k_){return true;}
		
		synchronized(RandomHexamer.class){
			if(initialized && k==k_){return true;} // Double-checked locking
			
			k=k_;
			final int cardinality=1<<(2*k);
			probs=new float[cardinality];
			for(int i=0; i<probs.length; i++){
				probs[i]=function(randy);
			}
			initialized=true;
			return true;
		}
	}
	
	/**
	 * Transforms one uniform float draw into a table-entry probability.
	 * The formula is {@code minProb+(1-minProb)*random^power}; it scales and shifts
	 * the powered draw rather than imposing a cutoff. Float rounding may produce one.
	 *
	 * @param randy Random number generator
	 * @return Keep probability between minProb and one for valid settings and draws
	 */
	private static float function(Random randy){
		return minProb+(1-minProb)*(float)Math.pow(randy.nextFloat(), power);
	}
	
	/** K-mer length for hexamer analysis (default 6) */
	private static int k=6;
	/** Table validity flag; changed settings clear it. Not a concurrent publication mechanism. */
	private static boolean initialized=false;
	/** Shared probability table indexed by the encoded prefix, replaced during initialization. */
	private static float[] probs;
	/** Positive exponent; larger values lower probabilities for draws strictly between zero and one. */
	private static float power=0.5f;
	/** Lower endpoint of the affine probability transform, default 0.1. */
	private static float minProb=0.1f;
}
