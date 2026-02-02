package aligner;

import shared.Random;

import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;
import map.IntHashMap2;

/**
 * Calculates the minimum seed hits required to detect indel-free alignments at a target probability.
 * Optimized version with "Fail Fast" failure budget and optional (currently disabled) SIMD path.
 *
 * @author Brian Bushnell
 * @contributor Amber
 * @date February 2, 2026
 */
public class MinHitsCalculator3 {

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
	public MinHitsCalculator3(int k_, int maxSubs_, float minid_, int midMaskLen_, float minProb_, float maxClip_, int kStep_){
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
		final int stepMask=step-1;
		final int stepTarget=(k-1)&stepMask;

		for(int i=0; i<queryLen; i++){
			// Roll error pattern: shift right and add new bit
			errorPattern=(errorPattern>>1)|(errors[i]<<(k-1));

			// Check if we have a full k-mer (i>=k-1) at a step position with no errors in non-wildcard positions
			boolean valid=(i>=k-1) && (i&stepMask)==stepTarget && ((errorPattern&wildcardMask)==0);
			count+=valid ? 1 : 0;
		}
		return count;
	}

	private static int upperBoundValid(int validKmers, int subs, int k){
		return Math.max(0, validKmers-subs);
	}

	private static int lowerBoundValid(int validKmers, int subs, int k, int mm){
		int kEff=k-mm;
		return Math.max(0, validKmers-subs*kEff);
	}

	private static int expectedUpperBoundValid(int validKmers, int subs, int k, int mm){
		int kEff=k-mm;
		return (int)Math.ceil(Math.max(0, validKmers-subs*kEff*0.45f));
	}

	/** * If this returns a value of at least minHits, it is at least feasible
	 * that simulation could render this a usable kmer length.
	 * @param validKmers Number of valid k-mers in the query
	 * @return A safe upper bound of remaining valid kmers; 0 means failure
	 */
	private int simulateFast(int validKmers) {
		int queryLen=validKmers+k-1;
		final int maxSubs=Math.min(maxSubs0, (int)(queryLen*(1-minid)));
		return expectedUpperBoundValid(validKmers, maxSubs, k, midMaskLen);
	}

	/**
	 * Runs Monte Carlo simulation with "Fail Fast" failure budget.
	 * @param validKmers Number of valid k-mers in the query
	 * @return Minimum hits needed; 0 means failure
	 */
	private int simulate(int validKmers, int iters){
		if(simulateFast(validKmers)<1) {return 0;}

		int queryLen=validKmers+k-1;
		final int maxSubs=Math.min(maxSubs0, (int)(queryLen*(1-minid)));
		int maxClips=(maxClipFraction<1 ? (int)(maxClipFraction*queryLen) : (int)maxClipFraction);

		// Deterministic case: require all possible hits
		if(minProb>=1){
			int unmasked=(Tools.max(2, k-midMaskLen));// Number of kmers impacted by a sub
			return Math.max(1, validKmers-(unmasked*maxSubs)-maxClips);
		}else if(minProb<=0){
			return validKmers;
		}

		// Calculate the Failure Budget
		// If we exceed this many zero-hit results, it is statistically impossible to meet minProb.
		final int maxFailures = (int)(iters * (1.0f - minProb));

		// Heuristic: If we are in the early game (iter < Limit) and have already blown a fraction of the budget, abort.
		// "if (histogram[0]*4 > maxZeros && iter*16 < iters)"
		final int heuristicIterLimit=iters/16;
		final int heuristicFailureLimit=(maxFailures+3)/4;

		int currentFailures = 0;

		int[] histogram=new int[validKmers+1];

		// Check SIMD capability
		final int lanes = Vector.INT_LANES;
		final boolean useSIMD = (USE_SIMD && lanes > 1 && Shared.SIMD);

		if(useSIMD) {
			// Vectorized Batched Mode
			final int batchSize = lanes;
			final int[] errorBuffer = new int[queryLen * batchSize];
			final int[] results = new int[batchSize];
			final int[] touchedIndices = new int[maxSubs * batchSize];

			for(int iter=0; iter<iters; iter+=batchSize){
				// Fill errors for 'batchSize' simulations
				int touchedCount = 0;
				for(int l=0; l<batchSize; l++){
					for(int s=0; s<maxSubs; s++){
						int pos = randy.nextInt(queryLen);
						int idx = pos * batchSize + l; // Interleaved layout
						if(errorBuffer[idx] == 0) {
							errorBuffer[idx] = 1;
							touchedIndices[touchedCount++] = idx;
						}
					}
				}

				// Run Vector Kernel
				Vector.countErrorFreeKmersBatch(errorBuffer, results, k, queryLen, kStep, wildcardMask);

				// Process Results
				for(int l=0; l<batchSize; l++){
					int count = results[l];
					if(count == 0) {
						currentFailures++;
						if(currentFailures > maxFailures) {
							return 0; // Fail Fast (Absolute)
						}
						// Note: Heuristic check skipped for SIMD for simplicity, but could be added
					}
					histogram[count]++;
				}

				// Cleanup error buffer
				for(int i=0; i<touchedCount; i++){
					errorBuffer[touchedIndices[i]] = 0;
				}
			}
		} else {
			// Scalar Fallback (Original Logic)
			int[] errors=new int[queryLen];
			for(int iter=0; iter<iters; iter++){
				// Clear errors
				for(int i=0; i<queryLen; i++){errors[i]=0;}

				// Place maxSubs random errors
				for(int i=0; i<maxSubs; i++){
					int pos=randy.nextInt(queryLen);
					errors[pos]=1;
				}

				// Count k-mers that survive the errors
				int errorFreeKmers=countErrorFreeKmers(errors, queryLen, kStep);

				if(errorFreeKmers == 0) {
					currentFailures++;
					// Absolute budget check
					if(currentFailures > maxFailures) {
						return 0;
					}
					// Heuristic check: Abort if we fail too much too early (trash detection)
					if(iter < heuristicIterLimit && currentFailures > heuristicFailureLimit) {
						return 0;
					}
				}

				histogram[errorFreeKmers]++;
			}
		}

		// Print histogram if verbose
		if(verbose){
			printHistogram(histogram);
		}

		// Find threshold that captures minProb fraction of cases
		int targetCount=(int)(iters*minProb);
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
					minHits=Math.max(0, simulate(validKmers, iterations));
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

	final int k;
	private final int maxSubs0;
	private final float minid;
	final int midMaskLen;
	final float maxClipFraction;
	private final float minProb;
	public final int kStep;
	private final int wildcardMask;
	private final IntHashMap2 validKmerToMinHits=new IntHashMap2();
	private final Random randy=Shared.threadLocalRandom(1);
	public static int iterations=200000;
	private static final boolean verbose=false; // Set to true for debugging

	/** Controls whether SIMD path is attempted. Currently disabled due to RNG setup overhead. */
	private static final boolean USE_SIMD=false;

	/*--------------------------------------------------------------*/
	/*----------------        Debug Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Main method for standalone testing and debugging.
	 * Usage: java MinHitsCalculator3 verbose=true k=13 validKmers=50 maxsubs=5 minid=0.9 midmask=1 minprob=0.99 maxclip=0.25 kstep=1 iterations=10000
	 */
	public static void main(String[] args){
		int k=13, validKmers=50, maxSubs=5, midMaskLen=1, kStep=1, iters=10000;
		float minid=0.9f, minProb=0.99f, maxClip=0.25f;

		for(String arg : args){
			String[] split=arg.split("=");
			if(split.length<2){continue;}
			String a=split[0].toLowerCase(), b=split[1];

			if(a.equals("verbose")){/*verbose=Boolean.parseBoolean(b);*/assert(false) : "Verbose is final.";}
			else if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("validkmers")){validKmers=Integer.parseInt(b);}
			else if(a.equals("maxsubs")){maxSubs=Integer.parseInt(b);}
			else if(a.equals("minid")){minid=Float.parseFloat(b);}
			else if(a.equals("midmask") || a.equals("midmasklen")){midMaskLen=Integer.parseInt(b);}
			else if(a.equals("minprob")){minProb=Float.parseFloat(b);}
			else if(a.equals("maxclip")){maxClip=Float.parseFloat(b);}
			else if(a.equals("kstep") || a.equals("step")){kStep=Integer.parseInt(b);}
			else if(a.equals("iterations")){iters=Integer.parseInt(b);}
		}
		iterations=iters;

		System.err.println("MinHitsCalculator3 testing:");
		System.err.println("  k="+k+" validKmers="+validKmers+" maxSubs="+maxSubs+" minid="+minid);
		System.err.println("  midMaskLen="+midMaskLen+" minProb="+minProb+" maxClip="+maxClip+" kStep="+kStep);
		System.err.println("  iterations="+iterations+" verbose="+verbose);
		System.err.println("  SIMD available: "+Shared.SIMD+", lanes="+Vector.INT_LANES+", useSimd="+USE_SIMD);

		Timer t=new Timer();
		MinHitsCalculator3 mhc=new MinHitsCalculator3(k, maxSubs, minid, midMaskLen, minProb, maxClip, kStep);
		t.stopAndPrint();

		t.start();
		int minHits=mhc.minHits(validKmers);
		t.stopAndPrint();

		System.err.println("\nResult: minHits="+minHits);
	}

	/**
	 * Prints histogram of error-free k-mer counts.
	 * @param histogram Array where histogram[i] is count of iterations with i error-free kmers
	 */
	static void printHistogram(int[] histogram){
		if(!verbose){return;}
		System.err.println("Histogram (errorFreeKmers -> count):");
		for(int i=0; i<histogram.length; i++){
			if(histogram[i]>0){
				System.err.println("  "+i+": "+histogram[i]);
			}
		}
	}
}