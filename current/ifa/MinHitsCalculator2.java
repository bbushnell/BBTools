package ifa;

import shared.Random;

import shared.Shared;
import shared.Timer;
import shared.Tools;
import map.IntHashMap2;

/**
 * Cached seed-threshold estimator used by Query, with a rolling one-bit error mask.
 * Models a contiguous query, sampled windows, central wildcards and replacement error draws.
 * A fast heuristic can reduce the trial count; zero-hit budgets can return zero early.
 * Clipping adjusts threshold expressions rather than trial positions. Neither the
 * histogram nor the one-hit endpoint floor is a universal detection guarantee.
 * Use instances serially: initial cache reads are outside the miss-path lock (IFA-005).
 * Configure iterations before use; existing cached thresholds are not invalidated.
 *
 * @author Brian Bushnell
 * @contributor Noire
 * @date December 30, 2025
 */
public class MinHitsCalculator2{

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
			//TODO: The expression could return 1 or 0 using clever Math.max subexpressions
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
	 * A result below one makes simulate reduce the requested trial count tenfold;
	 * this estimate is not a rigorous upper bound or an independent acceptance decision.
	 * @param validKmers Valid-window count for the contiguous model
	 * @return Nonnegative screening estimate
	 */
	private int simulateFast(int validKmers){
		int queryLen=validKmers+k-1;
		final int maxSubs=Math.min(maxSubs0, (int)(queryLen*(1-minid)));
		return expectedUpperBoundValid(validKmers, maxSubs, k, midMaskLen);
	}

	/**
	 * Calculates a raw model threshold, after optionally reducing iters by integer division by ten.
	 * For probability at least one, uses a sampled-window loss bound with a one-hit floor;
	 * zero requests every sampled window and negative values return one.
	 * Interior probabilities draw error positions with replacement. A ceil-based zero-hit
	 * budget and an early-trial heuristic can return zero before all trials finish.
	 * Otherwise the histogram is walked from high counts to low until its upper tail
	 * reaches the truncated iters*minProb target, then capped by validKmers-maxSubs-maxClips.
	 * Clipping does not remove trial positions. The caller clamps negative results to zero.
	 * @param validKmers Number of windows represented by a contiguous query of validKmers+k-1 bases
	 * @param iters Requested trial count, possibly reduced by the fast heuristic
	 * @return Raw threshold, or zero on an early rejection
	 */
	private int simulate(int validKmers, int iters){
		if(simulateFast(validKmers)<1){iters/=10;}
		// Calculate effective clipping limit for this query length
		int queryLen=validKmers+k-1;
		final int maxSubs=Math.min(maxSubs0, (int)(queryLen*(1-minid)));
		int maxClips=(maxClipFraction<1 ? (int)(maxClipFraction*queryLen) : (int)maxClipFraction);

		// Probability shortcuts, after the fast screening policy above.
		//FIXED IFA-008: thresholds count sampled windows, as do the seed consumers and simulation counter.
		if(minProb>=1){
			final int sampled=(int)(((long)validKmers+kStep-1)/kStep);
			int unmasked=(Tools.max(2, k-midMaskLen));// Number of kmers impacted by a sub
			return Math.max(1, sampled-(unmasked*maxSubs)-maxClips);
		}else if(minProb==0){
			return (int)(((long)validKmers+kStep-1)/kStep);
		}else if(minProb<0){
			return 1;
		}

		// Build histogram of surviving k-mer counts
		int[] histogram=new int[validKmers+1];
		int[] errors=new int[queryLen];// Owned 0/1 flags, reused for all trials

		// Run Monte Carlo simulation
		final int maxZeros=(int)(Math.ceil((1-minProb)*iters));
		final int earlyIters=iters/16, earlyZeros=(maxZeros+3)/4;
		for(int iter=0; iter<iters; iter++){
			// Clear errors
			for(int i=0; i<queryLen; i++){errors[i]=0;}

			// Duplicate random draws set the same flag, rather than adding a distinct error.
			for(int i=0; i<maxSubs; i++){
				int pos=randy.nextInt(queryLen);
				errors[pos]=1;
			}

			// Count k-mers that survive the errors
			int errorFreeKmers=countErrorFreeKmers(errors, queryLen, kStep);
			histogram[errorFreeKmers]++;
			if(histogram[0]>maxZeros || (histogram[0]>earlyZeros && iter<earlyIters)){
				return 0;
			}//Early exit, success is unlikely

			// Print iterations if verbose
			if(verbose){
				System.err.println("\nIteration "+(iter+1)+" (validKmers="+validKmers+"):");
				printSequence(errors, queryLen);
				System.err.println("Error-free kmers: "+errorFreeKmers);
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

		System.err.println("MinHitsCalculator2 testing:");
		System.err.println("  k="+k+" validKmers="+validKmers+" maxSubs="+maxSubs+" minid="+minid);
		System.err.println("  midMaskLen="+midMaskLen+" minProb="+minProb+" maxClip="+maxClip+" kStep="+kStep);
		System.err.println("  iterations="+iterations+" verbose="+verbose);
		Timer t=new Timer();
		MinHitsCalculator2 mhc=new MinHitsCalculator2(k, maxSubs, minid, midMaskLen, minProb, maxClip, kStep);
		t.stopAndPrint();
		System.err.println("Wildcard mask (int bits): "+Integer.toBinaryString(mhc.wildcardMask));
		System.err.println("Wildcard mask visual:");
		for(int i=k-1; i>=0; i--){
			System.err.print((mhc.wildcardMask&(1<<i))!=0 ? "m" : "W");
		}
		System.err.println();

		t.start();
		int minHits=mhc.minHits(validKmers);
		t.stopAndPrint();

		System.err.println("\nResult: minHits="+minHits);
	}

	/**
	 * Prints m/S error flags when compiled verbose; currently does nothing.
	 * @param errors Per-position 0/1 flags; read only
	 * @param queryLen Logical query length
	 */
	static void printSequence(int[] errors, int queryLen){
		if(!verbose){return;}
		StringBuilder sb=new StringBuilder(queryLen);
		for(int i=0; i<queryLen; i++){
			sb.append(errors[i]==1 ? 'S' : 'm');
		}
		System.err.println(sb.toString());
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
