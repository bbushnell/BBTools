package aligner;

import java.util.Random;

import shared.Shared;
import shared.Tools;
import structures.IntHashMap;

/**
 * Calculates the minimum seed hits required to detect indel-free alignments at a target probability.
 * Uses Monte Carlo simulation to model wildcards, error patterns, and clipping limits.
 * Optimized version with branchless operations and direct error storage.
 *
 * @author Brian Bushnell
 * @contributor Noire
 * @date December 30, 2024
 */
public class MinHitsCalculator2 {

	/**
	 * Constructs the calculator and precomputes wildcard masks.
	 *
	 * @param k_ K-mer length
	 * @param maxSubs_ Maximum allowed substitutions
	 * @param minid_ Minimum identity allowed
	 * @param midMaskLen_ Number of wildcard bases in the middle of the k-mer
	 * @param minProb_ Minimum detection probability (0.0-1.0)
	 * @param maxClip_ Maximum clipping allowed (fraction <1 or absolute ≥1)
	 * @param kStep_ K-mer step size (1 for all kmers, 2 for every other kmer, etc.)
	 */
	public MinHitsCalculator2(int k_, int maxSubs_, float minid_, int midMaskLen_, float minProb_, float maxClip_, int kStep_){
		k=k_;
		maxSubs0=maxSubs_;
		minid=minid_;
		midMaskLen=midMaskLen_;
		minProb=minProb_;
		maxClipFraction=maxClip_;
		kStep=Math.max(1, kStep_);

		// Build wildcard mask for error checking (1-bit per position)
		// Start with all k bits set
		int wildcardMask_=(1<<k)-1;
		// Clear the middle bits (wildcard positions)
		int midStart=(k-midMaskLen)/2;
		for(int i=0; i<midMaskLen; i++){
			wildcardMask_&=(~(1<<(midStart+i)));
		}
		wildcardMask=wildcardMask_;
	}

	/**
	 * Counts k-mers unaffected by errors, honoring wildcard positions.
	 * Branchless implementation using bit operations.
	 *
	 * @param errors Array of error positions (0 or 1 for each position)
	 * @param queryLen Query length
	 * @param step Step size for sampling k-mer positions (models reference indexing step)
	 * @return Number of error-free k-mers
	 */
	private int countErrorFreeKmers(int[] errors, int queryLen, int step){
		int count=0;
		int errorPattern=0;
		int len=0;
		final int stepMask=step-1;
		final int stepTarget=(k-1)&stepMask;

		for(int i=0; i<queryLen; i++){
			// Roll error pattern: shift right and add new bit
			errorPattern=(errorPattern>>1)|(errors[i]<<(k-1));

			// Branchless len tracking: reset to 0 on error, else increment
			len=(len+1)*(1-errors[i]);

			// Branchless check and count
			boolean valid=(len>=k) && (i&stepMask)==stepTarget && ((errorPattern&wildcardMask)==0);
			count+=valid ? 1 : 0;
		}
		return count;
	}

	/**
	 * Runs Monte Carlo simulation to find the minimum hits satisfying the probability target.
	 * @param validKmers Number of valid k-mers in the query
	 * @return Minimum hits needed
	 */
	private int simulate(int validKmers){
		// Calculate effective clipping limit for this query length
		int queryLen=validKmers+k-1;
		final int maxSubs=Math.min(maxSubs0, (int)(queryLen*(1-minid)));
		int maxClips=(maxClipFraction<1 ? (int)(maxClipFraction*queryLen) : (int)maxClipFraction);

		// Deterministic case: require all possible hits
		if(minProb>=1){
			int unmasked=(Tools.max(2, k-midMaskLen));// Number of kmers impacted by a sub
			return Math.max(1, validKmers-(unmasked*maxSubs)-maxClips);
		}else if(minProb==0){
			return validKmers;
		}else if(minProb<0){
			return 1;
		}

		// Build histogram of surviving k-mer counts
		int[] histogram=new int[validKmers+1];
		int[] errors=new int[queryLen]; // Direct error storage (0 or 1)

		// Run Monte Carlo simulation
		for(int iter=0; iter<iterations; iter++){
			// Clear errors
			for(int i=0; i<queryLen; i++){errors[i]=0;}

			// Place maxSubs random errors
			for(int i=0; i<maxSubs; i++){
				int pos=randy.nextInt(queryLen);
				errors[pos]=1;
			}

			// Count k-mers that survive the errors
			int errorFreeKmers=countErrorFreeKmers(errors, queryLen, kStep);
			histogram[errorFreeKmers]++;
		}

		// Find threshold that captures minProb fraction of cases
		int targetCount=(int)(iterations*minProb);
		int cumulative=0;

		// Walk down from highest hit count to find percentile threshold
		for(int hits=validKmers; hits>=0; hits--){
			cumulative+=histogram[hits];
			if(cumulative>=targetCount){
				// Don't exceed theoretical maximum after clipping
				return Math.min(hits, validKmers-maxSubs-maxClips);
			}
		}

		return Math.max(1, validKmers-maxSubs-maxClips); // Fallback
	}

	/**
	 * Returns the minimum seed hits for the given valid k-mer count, caching results.
	 * @param validKmers Valid k-mer count
	 * @return Minimum hits needed
	 */
	public int minHits(int validKmers){
		int minHits=validKmerToMinHits.get(validKmers);
		if(minHits<0 && !validKmerToMinHits.contains(validKmers)){
			synchronized(validKmerToMinHits) {
				if(!validKmerToMinHits.contains(validKmers)){
					minHits=Math.max(0, simulate(validKmers));
					validKmerToMinHits.put(validKmers, minHits);
				}else{
					minHits=validKmerToMinHits.get(validKmers);
				}
			}
		}
		return minHits;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int k;
	private final int maxSubs0;
	private final float minid;
	private final int midMaskLen;
	private final float maxClipFraction;
	private final float minProb;
	public final int kStep;
	private final int wildcardMask;
	private final IntHashMap validKmerToMinHits=new IntHashMap();
	private final Random randy=Shared.threadLocalRandom(1);
	public static int iterations=100000;
}
