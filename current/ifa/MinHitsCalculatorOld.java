package ifa;

import java.util.BitSet;

import map.IntHashMap;
import shared.Random;
import shared.Shared;

/**
 * Legacy cached estimator of a seed-hit threshold for contiguous query windows.
 * Models substitutions and central wildcards, without clipping or sampling stride.
 * Interior probabilities use simulated error-position draws with replacement, so
 * a trial can contain fewer distinct errors than the substitution allowance.
 * The result describes this model; it is not a general alignment-detection guarantee.
 * Use each instance serially: its cache and random state mutate without synchronization.
 * Configure iterations before use; changing it does not invalidate cached thresholds.
 * @author Brian Bushnell
 */
public class MinHitsCalculatorOld{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Stores model parameters and creates an owned central-wildcard pattern.
	 * Thresholds are calculated lazily on cache misses, not during construction.
	 *
	 * @param k_ K-mer length
	 * @param maxSubs_ Substitution allowance; number of random position draws in each simulated trial
	 * @param midMaskLen_ Number of wildcard bases in middle of k-mer
	 * @param minProb_ Upper-tail fraction for simulation; 0 requests all windows,
	 * at least 1 uses a loss bound, and negative values request one hit
	 */
	public MinHitsCalculatorOld(int k_, int maxSubs_, int midMaskLen_, float minProb_){
		k=k_;
		maxSubs=maxSubs_;
		midMaskLen=midMaskLen_;
		minProb=minProb_;

		// Pre-compute wildcard pattern for efficient simulation
		wildcards=makeWildcardPattern(k, midMaskLen);

		// Retained legacy masks; the simulation uses wildcards instead.
		kMask=~((-1)<<(2*k));

		// Calculate the stored middle mask without changing initialization order.
		int bitsPerBase=2;
		int bits=midMaskLen*bitsPerBase;
		int shift=((k-midMaskLen)/2)*bitsPerBase;
		midMask=~((~((-1)<<bits))<<shift);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates an owned wildcard pattern; true positions do not need to match.
	 * With unequal flanks, the central block starts toward the lower-index side.
	 *
	 * @param k K-mer length
	 * @param midMaskLen Number of consecutive wildcard bases in middle
	 * @return Boolean array where true indicates wildcard position
	 */
	private boolean[] makeWildcardPattern(int k, int midMaskLen){
		boolean[] wildcards=new boolean[k];
		// Default false: non-wildcard positions must match exactly

		// Round the left flank down when the two flanks cannot have equal length.
		int start=(k-midMaskLen)/2;
		for(int i=0; i<midMaskLen; i++){
			wildcards[start+i]=true;
		}
		return wildcards;
	}

	/**
	 * Counts k-mers that would still match despite errors, accounting for wildcards.
	 * A k-mer is considered error-free if no errors fall on non-wildcard positions.
	 *
	 * @param errors Error positions; read only
	 * @param wildcards Wildcard positions within a window; read only
	 * @param queryLen Logical query length, independent of the BitSet storage capacity
	 * @return Number of surviving windows, checking every start at stride one
	 */
	private int countErrorFreeKmers(BitSet errors, boolean[] wildcards, int queryLen){
		int count=0;
		// Use the logical query extent, not errors.size(), which is storage capacity.
		for(int i=0; i<=queryLen-k; i++){
			boolean errorFree=true;
			for(int j=0; j<k && errorFree; j++){
				errorFree=wildcards[j] || (!errors.get(i+j));
			}
			count+=errorFree ? 1 : 0;
		}
		return count;
	}

	/**
	 * Calculates a raw threshold, using shortcuts outside the open probability interval.
	 * For probabilities between 0 and 1, uses a contiguous query of validKmers+k-1
	 * bases and returns the largest hit count whose histogram upper tail reaches
	 * floor(iterations*minProb). Repeated error draws collapse in the BitSet.
	 * For minProb at least 1, subtracts (k-midMaskLen)*maxSubs from validKmers;
	 * the result may be negative until minHits clamps it.
	 *
	 * @param validKmers Number of windows represented by the contiguous simulation
	 * @return Raw threshold; callers clamp negative values to zero
	 */
	private int simulate(int validKmers){
		// Deterministic loss bound; the zero-probability branch requests all windows.
		if(minProb>=1){
			return validKmers-(k-midMaskLen)*maxSubs;
		}else if(minProb==0){
			return validKmers;
		}else if(minProb<0){return 1;}

		// Build histogram of surviving k-mer counts
		int[] histogram=new int[validKmers+1];
		int queryLen=validKmers+k-1;// Length needed to generate validKmers windows
		BitSet errors=new BitSet(queryLen);// Reused for all trials

		// Run Monte Carlo simulation
		for(int iter=0; iter<iterations; iter++){
			errors.clear();

			// Draw with replacement: duplicate positions do not add another distinct error.
			for(int i=0; i<maxSubs; i++){
				int pos=randy.nextInt(queryLen);
				errors.set(pos);
			}

			// Count k-mers that survive the errors
			int errorFreeKmers=countErrorFreeKmers(errors, wildcards, queryLen);
			histogram[errorFreeKmers]++;
		}

		// The integer target truncates the requested fraction of trials.
		int targetCount=(int)(iterations*minProb);
		int cumulative=0;

		// Walk down from highest hit count to find percentile threshold
		for(int hits=validKmers; hits>=0; hits--){
			cumulative+=histogram[hits];
			if(cumulative>=targetCount){
				return hits;
			}
		}

		return 0;// Fallback when the histogram walk does not reach the target
	}

	/**
	 * Returns a cached model threshold, simulating on a miss and clamping it to zero.
	 * The key is only the window count; changing iterations leaves prior entries intact.
	 * Cache-miss order consumes the instance's random stream and can affect later estimates.
	 * @param validKmers Number of valid windows, represented as contiguous during simulation
	 * @return Nonnegative seed threshold for this model, not a general detection guarantee
	 */
	public int minHits(int validKmers){
		if(!validKmerToMinHits.contains(validKmers)){
			int minHits=Math.max(0, simulate(validKmers));
			validKmerToMinHits.put(validKmers, minHits);
		}
		return validKmerToMinHits.get(validKmers);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Bases per seed window. */
	private final int k;

	/** Substitution allowance and number of random position draws per trial. */
	private final int maxSubs;

	/** Number of central wildcard positions in each window. */
	private final int midMaskLen;

	/** Stored legacy two-bit window mask; not read by this implementation. */
	private final int kMask;

	/** Stored legacy central-bit mask; simulation reads wildcards instead. */
	private final int midMask;

	/** Requested histogram upper-tail fraction, with separate endpoint/negative shortcuts. */
	private final float minProb;

	/** Owned immutable-in-use pattern; true positions ignore errors. */
	private final boolean[] wildcards;

	/** Instance cache keyed by window count; accessed without synchronization. */
	private final IntHashMap validKmerToMinHits=new IntHashMap();

	/** Instance generator initialized with seed 1; the current Shared factory creates a new generator. */
	private final Random randy=Shared.threadLocalRandom(1);

	/** Trial count for subsequent cache misses; configure before use, no cache invalidation. */
	public static int iterations=100000;
}
