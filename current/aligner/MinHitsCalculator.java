package aligner;

import java.util.BitSet;
import java.util.Random;

import shared.Shared;
import shared.Tools;
import structures.IntHashMap;

/** 
 * Calculates the minimum number of seed hits needed
 * to ensure all valid indel-free alignments are found
 * with a specified probability threshold.
 * Uses Monte Carlo simulation to account for wildcard
 * masking, error patterns, and clipping allowances.
 * 
 * @author Brian Bushnell
 * @contributor Isla SOS
 * @date June 4, 2025
 */
public class MinHitsCalculator {

	/**
	 * Constructor for minimum hits calculator.
	 * @param k_ K-mer length
	 * @param maxSubs_ Maximum allowed substitutions in alignment
	 * @param midMaskLen_ Number of wildcard bases in middle of k-mer
	 * @param minProb_ Minimum probability of detecting valid alignments (0.0-1.0)
	 * @param maxClip_ Maximum clipping allowed (fraction <1 or absolute ≥1)
	 */
	public MinHitsCalculator(int k_, int maxSubs_, int midMaskLen_, float minProb_, float maxClip_){
		k=k_;
		maxSubs=maxSubs_;
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
	 * Creates a boolean array indicating which k-mer positions are wildcarded.
	 * Wildcard positions are set to true for efficient short-circuiting.
	 * @param k K-mer length
	 * @param midMaskLen Number of consecutive wildcard bases in middle
	 * @return Boolean array where true indicates wildcard position
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
	 * Counts k-mers that would still match despite errors, accounting for wildcards.
	 * A k-mer is considered error-free if no errors fall on non-wildcard positions.
	 * @param errors BitSet indicating error positions in query sequence
	 * @param wildcards Boolean array indicating wildcard positions in k-mer
	 * @param queryLen Length of query sequence
	 * @return Number of k-mers that remain matchable
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
	 * Simulates random error patterns to determine minimum seed hits needed.
	 * Uses Monte Carlo simulation to find the threshold that ensures
	 * the specified probability of detecting valid alignments.
	 * @param validKmers Number of valid k-mers in query sequence
	 * @return Minimum number of seed hits needed
	 */
	private int simulate(int validKmers){
		// Calculate effective clipping limit for this query length
		int queryLen=validKmers+k-1;
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
		BitSet errors=new BitSet(queryLen); // Reuse BitSet for efficiency

		// Run Monte Carlo simulation
		for(int iter=0; iter<iterations; iter++){
			errors.clear();

			// Place maxSubs random errors in query
			for(int i=0; i<maxSubs; i++){
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
				return Math.min(hits, validKmers-maxSubs-maxClips);
			}
		}

		return Math.max(1, validKmers-maxSubs-maxClips); // Fallback
	}

	/**
	 * Gets the minimum number of seed hits required for a query with
	 * the specified number of valid k-mers. Results are cached.
	 * @param validKmers Number of valid k-mers in query sequence
	 * @return Minimum seed hits needed to ensure detection probability
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

	/** K-mer length */
	private final int k;
	/** Maximum allowed substitutions in alignment */
	private final int maxSubs;
	/** Number of wildcard bases in middle of k-mer */
	private final int midMaskLen;
	/** Maximum clipping allowed as fraction or absolute value */
	private final float maxClipFraction;
	/** Bit mask for k-mer (may not be needed) */
	private final int kMask;
	/** Bit mask for middle wildcard positions (may not be needed) */
	private final int midMask;
	/** Minimum probability of detecting valid alignments (0.0-1.0) */
	private final float minProb;
	/** Boolean array indicating wildcard positions in k-mer */
	private final boolean[] wildcards;
	/** Cache mapping valid k-mer count to minimum required hits */
	private final IntHashMap validKmerToMinHits=new IntHashMap();
	/** Random number generator for simulation */
	private final Random randy=Shared.threadLocalRandom(1);
	/** Number of Monte Carlo iterations for simulation */
	public static int iterations=100000;
}