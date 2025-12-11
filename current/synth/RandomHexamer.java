package synth;

import java.util.Random;

import dna.AminoAcid;

/**
 * Models hexamer priming bias for realistic coverage patterns in synthetic read generation.
 * Creates sequence-dependent coverage variation by assigning different priming efficiencies
 * to different k-mer sequences, simulating the non-random nature of "random" hexamer priming
 * used in library preparation.
 *
 * @author Brian Bushnell
 * @author Isla Winglet
 */
public class RandomHexamer {

	/**
	 * Determines whether to keep a read based on hexamer priming bias at the start position.
	 * Evaluates the first k bases to compute a k-mer value and compares against the
	 * probability table to simulate sequence-dependent priming efficiency.
	 *
	 * @param bases The sequence bases to check
	 * @param randy Random number generator for probabilistic decision
	 * @return true if the read should be kept based on priming efficiency, false otherwise
	 */
	public static boolean keep(byte[] bases, Random randy) {
		assert(initialized) : "RandomHexamer must be initialized prior to use.";
		if(bases.length<k) {return true;}
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
	 * Sets all main parameters for hexamer bias modeling.
	 * Convenience method to configure k-mer length, power distribution, and minimum probability.
	 *
	 * @param k_ The k-mer length (typically 6 for hexamers)
	 * @param power_ The exponent for power-law distribution (typically 0.5)
	 * @param minProb_ Minimal relative occurrence frequency (typically 0.1)
	 */
	public static void set(int k_, float power_, float minProb_) {
		setK(k_);
		setPower(power_);
		setMinProb(minProb_);
	}
	
	/**
	 * Sets the minimum probability threshold for any k-mer.
	 * Values below this threshold will be set to this minimum, preventing
	 * complete exclusion of any k-mer sequence.
	 * @param minProb_ Minimal relative occurrence frequency (typically 0.1)
	 */
	public static void setMinProb(float minProb_) {
		initialized=initialized && (minProb==minProb_);
		minProb=minProb_;
		assert(minProb>=0 && minProb<1) : minProb_;
	}
	
	/**
	 * Sets the power parameter for the probability distribution function.
	 * Lower values create more bias (greater variation between k-mers),
	 * while higher values create more uniform distribution.
	 * @param power_ The exponent for power-law distribution (typically 0.5)
	 */
	public static void setPower(float power_) {
		initialized=initialized && (power==power_);
		power=power_;
		assert(power>0) : power;
	}
	
	/**
	 * Sets the k-mer length for hexamer bias analysis.
	 * Determines the number of bases used to compute priming probability.
	 * @param k_ The k-mer length (typically 6 for hexamers, max 15)
	 */
	public static void setK(int k_) {
		initialized=initialized && (k==k_);
		k=k_;
		assert(k>0 && k<=15) : k;
	}
	
	/** Gets the current k-mer length setting */
	public static int getK() {return k;}
	
	/**
	 * Initializes the hexamer probability table using the current k value.
	 * Convenience method that calls initialize(randy, k).
	 * @param randy Random number generator for probability assignment
	 * @return true if initialization was successful
	 */
	public static boolean initialize(Random randy) {
		return initialize(randy, k);
	}
	
	/**
	 * Initializes the hexamer probability table with specified k-mer length.
	 * Creates probability array for all possible k-mers using power-law distribution
	 * to simulate realistic priming bias patterns. Uses double-checked locking
	 * for thread-safe initialization.
	 *
	 * @param randy Random number generator for probability assignment
	 * @param k_ The k-mer length to use for initialization
	 * @return true if initialization was successful
	 */
	public static boolean initialize(Random randy, int k_) {
		if(initialized && k==k_){return true;}
		
		synchronized(RandomHexamer.class){
			if(initialized && k==k_){return true;} // Double-checked locking
			
			k=k_;
			final int cardinality=1<<(2*k);
			probs=new float[cardinality];
			for(int i=0; i<probs.length; i++) {
				probs[i]=function(randy);
			}
			initialized=true;
			return true;
		}
	}
	
	/**
	 * Generates a biased probability value using power-law distribution.
	 * Applies the formula: minProb + (1-minProb) * random^power to create
	 * non-uniform priming probabilities across k-mers.
	 *
	 * @param randy Random number generator
	 * @return Probability value between minProb and 1.0
	 */
	private static float function(Random randy){
		return minProb+(1-minProb)*(float)Math.pow(randy.nextFloat(), power);
	}
	
	/** K-mer length for hexamer analysis (default 6) */
	private static int k=6;
	/** Flag indicating whether the probability table has been initialized */
	private static boolean initialized=false;
	/** Probability array indexed by k-mer numeric value */
	private static float[] probs;
	/** Power parameter for probability distribution (lower values = more bias) */
	private static float power=0.5f;
	/** Minimum probability threshold for any k-mer sequence */
	private static float minProb=0.1f;
}