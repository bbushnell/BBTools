package aligner;

import java.util.BitSet;
import java.util.Random;

import shared.Shared;
import shared.Tools;
import structures.IntHashMap;

/**
 * Calculates the minimum seed hits required to detect indel-free alignments at a target probability.
 * Uses Monte Carlo simulation to model wildcards, error patterns, and clipping limits.
 *
 * @author Brian Bushnell
 * @contributor Isla SOS
 * @date June 4, 2025
 */
public class MinHitsCalculator {

	/**
	 * Constructs the calculator and precomputes wildcard masks.
	 *
	 * @param k_ K-mer length
	 * @param maxSubs_ Maximum allowed substitutions
	 * @param minid_ Minimum identity allowed
	 * @param midMaskLen_ Number of wildcard bases in the middle of the k-mer
	 * @param minProb_ Minimum detection probability (0.0-1.0)
	 * @param maxClip_ Maximum clipping allowed (fraction <1 or absolute ≥1)
	 */
	public MinHitsCalculator(int k_, int maxSubs_, float minid_, int midMaskLen_, float minProb_, float maxClip_){
		k=k_;
		maxSubs0=maxSubs_;
		minid=minid_;
		midMaskLen=midMaskLen_;
		minProb=minProb_;
		maxClipFraction=maxClip_;

		// Pre-compute wildcard pattern for efficient simulation
		wildcards=makeWildcardPattern(k, midMaskLen);

		// Calculate bit mask for k-mer (may not be needed)
		kMask=~((-1)<<(2*k));

		// Calculate middle mask for wildcards (may not be needed)
		int bitsPerBase=2;
		int bits=midMaskLen*bitsPerBase;
		int shift=((k-midMaskLen)/2)*bitsPerBase;
		midMask=~((~((-1)<<bits))<<shift);
	}

	/**
	 * Builds a boolean array marking wildcard positions within a k-mer.
	 * @param k K-mer length
	 * @param midMaskLen Count of wildcard bases
	 * @return Boolean array with true for wildcard positions
	 */
	private boolean[] makeWildcardPattern(int k, int midMaskLen){
		boolean[] wildcards=new boolean[k];
		// Default false: non-wildcard positions must match exactly

		// Set wildcard positions to true (middle positions, right-shifted for even k)
		int start=(k-midMaskLen)/2;
		for(int i=0; i<midMaskLen; i++){
			wildcards[start+i]=true;
		}
		return wildcards;
	}

	/**
	 * Counts k-mers unaffected by errors, honoring wildcard positions.
	 *
	 * @param errors BitSet of error positions
	 * @param wildcards Wildcard position map
	 * @param queryLen Query length
	 * @return Number of error-free k-mers
	 */
	private int countErrorFreeKmers(BitSet errors, boolean[] wildcards, int queryLen){
		int count=0;

		// Check each possible k-mer position in query
		for(int i=0; i<=queryLen-k; i++){
			boolean errorFree=true;

			// Check each position within this k-mer
			for(int j=0; j<k && errorFree; j++){
				errorFree=wildcards[j]||(!errors.get(i+j));
			}
			if(errorFree){count++;}
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
			return Math.max(1, validKmers-(unmasked*maxSubs0)-maxClips);
		}else if(minProb==0){
			return validKmers;
		}else if(minProb<0){
			return 1;
		}
		
		// Build histogram of surviving k-mer counts
		int[] histogram=new int[validKmers+1];
		BitSet errors=new BitSet(queryLen); // Reuse BitSet for efficiency

		// Run Monte Carlo simulation
		for(int iter=0; iter<iterations; iter++){
			errors.clear();

			// Place maxSubs random errors in query
			for(int i=0; i<maxSubs0; i++){
				int pos=randy.nextInt(queryLen);
				errors.set(pos);
			}

			// Count k-mers that survive the errors
			int errorFreeKmers=countErrorFreeKmers(errors, wildcards, queryLen);
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
				return Math.min(hits, validKmers-maxSubs0-maxClips);
			}
		}

		return Math.max(1, validKmers-maxSubs0-maxClips); // Fallback
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
	private final int kMask;
	private final int midMask;
	private final float minProb;
	private final boolean[] wildcards;
	private final IntHashMap validKmerToMinHits=new IntHashMap();
	private final Random randy=Shared.threadLocalRandom(1);
	public static int iterations=100000;
}