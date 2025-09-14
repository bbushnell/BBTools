package synth;

import java.util.Arrays;
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
 * @date June 30, 2025
 */
public class RandomHexamer {

	/**
	 * Determines whether to keep a read based on hexamer priming bias at the start position.
	 * @param bases The sequence bases to check
	 * @param randy Random number generator
	 * @return True if the read should be kept based on priming efficiency
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
	 * Sets minimum probability for hexamers.
	 * @param k_ The k-mer length (typically 6 for hexamers)
	 * @param power_ The exponent (typically 0.5 for hexamers)
	 * @param minProb_ Minimal relative occurrence frequency (typically 0.1 for hexamers)
	 */
	public static void set(int k_, float power_, float minProb_) {
		setK(k_);
		setPower(power_);
		setMinProb(minProb_);
	}
	
	/**
	 * Sets minimum probability for hexamers.
	 * @param minProb_ Minimal relative occurrence frequency (typically 0.1 for hexamers)
	 */
	public static void setMinProb(float minProb_) {
		initialized=initialized && (minProb==minProb_);
		minProb=minProb_;
		assert(minProb>=0 && minProb<1) : minProb_;
	}
	
	/**
	 * Sets the power for random kmer probability density function.
	 * @param power_ The exponent (typically 0.5 for hexamers)
	 */
	public static void setPower(float power_) {
		initialized=initialized && (power==power_);
		power=power_;
		assert(power>0) : power;
	}
	
	/**
	 * Sets the k-mer length for hexamer analysis.
	 * @param k_ The k-mer length (typically 6 for hexamers)
	 */
	public static void setK(int k_) {
		initialized=initialized && (k==k_);
		k=k_;
		assert(k>0 && k<=15) : k;
	}
	
	/**
	 * Gets the current k-mer length.
	 * @return Current k value
	 */
	public static int getK() {return k;}
	
	/**
	 * Initializes the hexamer probability table with default k value.
	 * @param randy Random number generator for probability assignment
	 * @return True if initialization was successful
	 */
	public static boolean initialize(Random randy) {
		return initialize(randy, k);
	}
	
	/**
	 * Initializes the hexamer probability table with specified k value.
	 * Uses power-law distribution to create realistic priming bias.
	 * @param randy Random number generator for probability assignment
	 * @param k_ The k-mer length to use
	 * @return True if initialization was successful
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
	 * Generates a biased probability using power-law distribution.
	 * @param randy Random number generator
	 * @return Probability value between minProb and 1.0
	 */
	private static float function(Random randy){
		return minProb+(1-minProb)*(float)Math.pow(randy.nextFloat(), power);
	}
	
	/** K-mer length for hexamer analysis */
	private static int k=6;
	/** Whether the probability table has been initialized */
	private static boolean initialized=false;
	/** Probability array indexed by k-mer value */
	private static float[] probs;
	/** Power parameter for probability distribution (lower = more bias) */
	private static float power=0.5f;
	/** Minimum probability for any k-mer */
	private static float minProb=0.1f;
}