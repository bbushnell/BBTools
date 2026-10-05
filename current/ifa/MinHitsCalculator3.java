package ifa;

import shared.Random;

import shared.Shared;
import shared.Timer;
import shared.Tools;
import simd.Vector;
import map.IntHashMap2;

/**
 * Experimental cached seed-threshold estimator with a pre-screen and zero-hit failure budgets.
 * Uses a contiguous-window model, sampled starts, central wildcards and replacement error draws.
 * Unlike version 2, a zero fast-screen estimate returns zero before probability shortcuts;
 * a negative probability otherwise returns the unsampled window count. Preserve those distinctions.
 * The SIMD batch branch is disabled by a compile-time flag; lane information is still read.
 * Model heuristics and the one-hit endpoint floor are not universal detection guarantees.
 * Use instances serially: initial cache reads are outside the miss-path lock (IFA-005).
 * Configure iterations before use; changing it does not invalidate cached thresholds.
 *
 * @author Brian Bushnell
 * @contributor Amber
 * @date February 2, 2026
 */
public class MinHitsCalculator3{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Stores model parameters and creates a one-bit-per-position matching mask.
	 * The stride is clamped to at least one; simulation and caching are lazy.
	 *
	 * @param k_ Seed window length; current aligners use 1 through 15
	 * @param maxSubs_ Absolute substitution allowance, further limited by query length and identity
	 * @param minid_ Identity fraction used to derive the effective substitution allowance
	 * @param midMaskLen_ Number of wildcard bases in the middle of the k-mer
	 * @param minProb_ Model probability setting; see simulate for screening and shortcut behavior
	 * @param maxClip_ Maximum clipping allowed (fraction {@code <1} or absolute ≥1)
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

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Counts sampled windows without errors at unmasked positions.
	 * Rolls an error bitmask and samples complete windows at the requested stride.
	 *
	 * @param errors Per-position 0/1 error flags; read only
	 * @param queryLen Logical query length, within errors
	 * @param step Positive sampling stride; the first complete window starts at zero
	 * @return Number of surviving sampled windows
	 */
	private int countErrorFreeKmers(int[] errors, int queryLen, int step){
		int count=0;
		int errorPattern=0;
		//FIXED IFA-004: query sampling accepts any positive stride; step-1 masks only work for powers of two.
		int nextSample=k-1;

		for(int i=0; i<queryLen; i++){
			// Roll error pattern: shift right and add new bit
			errorPattern=(errorPattern>>1)|(errors[i]<<(k-1));

			// Check if we have a full k-mer (i>=k-1) at a step position with no errors in non-wildcard positions
			final boolean sampled=(i==nextSample);
			if(sampled){nextSample+=step;}
			boolean valid=sampled && ((errorPattern&wildcardMask)==0);
			count+=valid ? 1 : 0;
		}
		return count;
	}

	/** Unused legacy window-count estimate; k is retained but not used. */
	private static int upperBoundValid(int validKmers, int subs, int k){
		return Math.max(0, validKmers-subs);
	}

	/** Unused legacy loss estimate in unsampled-window units, clamped to zero. */
	private static int lowerBoundValid(int validKmers, int subs, int k, int mm){
		int kEff=k-mm;
		return Math.max(0, validKmers-subs*kEff);
	}

	/** Screening heuristic using a 0.45 loss factor; not a certified mathematical upper bound. */
	private static int expectedUpperBoundValid(int validKmers, int subs, int k, int mm){
		int kEff=k-mm;
		return (int)Math.ceil(Math.max(0, validKmers-subs*kEff*0.45f));
	}

	/**
	 * Estimates retained unsampled windows using the 0.45 loss heuristic and effective substitutions.
	 * A result below one causes simulate to return zero before any probability shortcut.
	 * This is a screening heuristic, not a certified mathematical upper bound.
	 * @param validKmers Valid-window count for the contiguous model
	 * @return Nonnegative screening estimate
	 */
	private int simulateFast(int validKmers){
		int queryLen=validKmers+k-1;
		final int maxSubs=Math.min(maxSubs0, (int)(queryLen*(1-minid)));
		return expectedUpperBoundValid(validKmers, maxSubs, k, midMaskLen);
	}

	/**
	 * Calculates a raw threshold, rejecting a zero fast-screen estimate before all shortcuts.
	 * If screening passes, probability at least one uses a sampled-window loss bound with
	 * a one-hit floor, zero requests all sampled starts, and a negative value returns validKmers.
	 * Interior probabilities draw error positions with replacement. A truncated zero-hit
	 * budget and a scalar early-trial heuristic can return zero before all trials finish.
	 * Otherwise histogram selection uses the truncated iters*minProb upper-tail target,
	 * with a post-selection validKmers-maxSubs-maxClips cap; trial positions are not clipped.
	 * The dormant batch path can process a full final batch beyond iters and skips the
	 * early-trial heuristic. This documentation does not establish SIMD/scalar equivalence.
	 * @param validKmers Windows represented by a contiguous query of validKmers+k-1 bases
	 * @param iters Requested trial count, also used to derive failure and selection budgets
	 * @return Raw threshold, or zero on a screening/budget rejection; caller clamps negatives
	 */
	private int simulate(int validKmers, int iters){
		if(simulateFast(validKmers)<1){return 0;}

		int queryLen=validKmers+k-1;
		final int maxSubs=Math.min(maxSubs0, (int)(queryLen*(1-minid)));
		int maxClips=(maxClipFraction<1 ? (int)(maxClipFraction*queryLen) : (int)maxClipFraction);

		// Probability shortcuts, after the fast screening policy above.
		//FIXED IFA-008: thresholds count sampled windows, as do the seed consumers and simulation counter.
		if(minProb>=1){
			final int sampled=(int)(((long)validKmers+kStep-1)/kStep);
			int unmasked=(Tools.max(2, k-midMaskLen));// Number of kmers impacted by a sub
			return Math.max(1, sampled-(unmasked*maxSubs)-maxClips);
		}else if(minProb<=0){
			// Preserve this experimental class's existing negative-probability behavior separately.
			return minProb==0 ? (int)(((long)validKmers+kStep-1)/kStep) : validKmers;
		}

		// Calculate the Failure Budget
		// Absolute cutoff for observed zero-hit trials; not a general probability guarantee.
		final int maxFailures=(int)(iters*(1.0f-minProb));

		// Heuristic: If we are in the early game (iter < Limit) and have already blown a fraction of the budget, abort.
		// "if (histogram[0]*4 > maxZeros && iter*16 < iters)"
		final int heuristicIterLimit=iters/16;
		final int heuristicFailureLimit=(maxFailures+3)/4;

		int currentFailures=0;

		int[] histogram=new int[validKmers+1];

		// Check SIMD capability
		final int lanes=Vector.INT_LANES;
		final boolean useSIMD=(USE_SIMD && lanes>1 && Shared.SIMD);

		if(useSIMD){
			// Vectorized Batched Mode
			final int batchSize=lanes;
			final int[] errorBuffer=new int[queryLen*batchSize];
			final int[] results=new int[batchSize];
			final int[] touchedIndices=new int[maxSubs*batchSize];

			//TODO: Probable bug IFA-012 - a partial final batch still counts every lane, but selection uses iters.
			//USE_SIMD is false, so this is dormant; review trial/budget accounting before enabling the branch.
			for(int iter=0; iter<iters; iter+=batchSize){
				// Fill errors for 'batchSize' simulations
				int touchedCount=0;
				for(int l=0; l<batchSize; l++){
					for(int s=0; s<maxSubs; s++){
						int pos=randy.nextInt(queryLen);
						int idx=pos*batchSize+l;// Interleaved layout
						if(errorBuffer[idx]==0){
							errorBuffer[idx]=1;
							touchedIndices[touchedCount++]=idx;
						}
					}
				}

				// Run Vector Kernel
				Vector.countErrorFreeKmersBatch(errorBuffer, results, k, queryLen, kStep, wildcardMask);

				// Process Results
				for(int l=0; l<batchSize; l++){
					int count=results[l];
					if(count==0){
						currentFailures++;
						if(currentFailures>maxFailures){
							return 0;// Absolute failure-budget cutoff
						}
						// Note: Heuristic check skipped for SIMD for simplicity, but could be added
					}
					histogram[count]++;
				}

				// Cleanup error buffer
				for(int i=0; i<touchedCount; i++){
					errorBuffer[touchedIndices[i]]=0;
				}
			}
		}else{
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

				if(errorFreeKmers==0){
					currentFailures++;
					// Absolute budget check
					if(currentFailures>maxFailures){
						return 0;
					}
					// Heuristic check: Abort if we fail too much too early (trash detection)
					if(iter<heuristicIterLimit && currentFailures>heuristicFailureLimit){
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
				// Preserve the model's post-selection substitution/clipping cap.
				return Math.min(hits, validKmers-maxSubs-maxClips);
			}
		}

		return Math.max(1, validKmers-maxSubs-maxClips); // Fallback
	}

	/**
	 * Returns a cached model threshold, calculating and clamping to zero on a miss.
	 * The key contains only window count; changing iterations leaves old entries intact.
	 * Miss order consumes the instance random stream and can affect later estimates.
	 * Initial cache reads are unlocked despite the synchronized miss path (IFA-005).
	 * @param validKmers Valid-window count, represented as contiguous in the model
	 * @return Nonnegative model threshold, not a universal detection guarantee
	 */
	public int minHits(int validKmers){
		//TODO: IFA-005 - concurrent use needs review: these reads race with cache writes; current Query setup is serial.
		int minHits=validKmerToMinHits.get(validKmers);
		if(minHits<0 && !validKmerToMinHits.contains(validKmers)){
			synchronized(validKmerToMinHits){
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

	/** Bases per seed window. */
	final int k;
	/** Absolute substitution allowance before the identity-based limit. */
	private final int maxSubs0;
	/** Identity fraction used with the modeled query length. */
	private final float minid;
	/** Number of central wildcard positions. */
	final int midMaskLen;
	/** Clipping setting: fraction below 1, otherwise absolute bases. */
	final float maxClipFraction;
	/** Probability setting for screening, endpoint shortcuts and histogram selection. */
	private final float minProb;
	/** Positive sampling stride, clamped during construction. */
	public final int kStep;
	/** One bit per window position: set bits require an error-free position, cleared bits are wildcards. */
	private final int wildcardMask;
	/** Instance window-count cache; initial reads are outside the miss-path lock. */
	private final IntHashMap2 validKmerToMinHits=new IntHashMap2();
	/** Instance generator seeded with 1; the current Shared factory creates a new generator. */
	private final Random randy=Shared.threadLocalRandom(1);
	/** Requested trial count for future cache misses; changes do not invalidate cached entries. */
	public static int iterations=200000;
	/** Compile-time diagnostic switch, not a mutable command-line option. */
	private static final boolean verbose=false;

	/** Compile-time batch switch, currently false; the historical rationale cites RNG setup overhead. */
	private static final boolean USE_SIMD=false;

	/*--------------------------------------------------------------*/
	/*----------------        Debug Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Runs a standalone model diagnostic and writes settings, timings and threshold to stderr.
	 * Example arguments: k=13 validkmers=50 maxsubs=5 minid=0.9 midmask=1
	 * minprob=0.99 maxclip=0.25 kstep=1 iterations=10000.
	 * Passing verbose asserts under -ea because the diagnostic switch is final.
	 * Sets the static iterations field before constructing the calculator.
	 * @param args Recognized key=value settings; keys are case insensitive
	 */
	public static void main(String[] args){
		int k=13, validKmers=50, maxSubs=5, midMaskLen=1, kStep=1, iters=10000;
		float minid=0.9f, minProb=0.99f, maxClip=0.25f;

		for(String arg : args){
			String[] split=arg.split("=");
			if(split.length<2){continue;}
			String a=split[0].toLowerCase(), b=split[1];

			if(a.equals("verbose")){/*verbose=Boolean.parseBoolean(b);*/assert(false) : "Verbose is final.";
			}else if(a.equals("k")){k=Integer.parseInt(b);
			}else if(a.equals("validkmers")){validKmers=Integer.parseInt(b);
			}else if(a.equals("maxsubs")){maxSubs=Integer.parseInt(b);
			}else if(a.equals("minid")){minid=Float.parseFloat(b);
			}else if(a.equals("midmask") || a.equals("midmasklen")){midMaskLen=Integer.parseInt(b);
			}else if(a.equals("minprob")){minProb=Float.parseFloat(b);
			}else if(a.equals("maxclip")){maxClip=Float.parseFloat(b);
			}else if(a.equals("kstep") || a.equals("step")){kStep=Integer.parseInt(b);
			}else if(a.equals("iterations")){iters=Integer.parseInt(b);}
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
	 * Prints nonzero histogram bins when compiled verbose; currently does nothing.
	 * @param histogram Trial counts indexed by surviving sampled-window count; read only
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
