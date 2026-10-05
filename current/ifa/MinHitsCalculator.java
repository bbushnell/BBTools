package ifa;

import java.util.BitSet;

import map.IntHashMap2;
import shared.Random;
import shared.Shared;
import shared.Timer;
import shared.Tools;

/**
 * Cached seed-threshold estimator for a contiguous-window query model.
 * Samples starts at the configured stride and models central wildcards and random
 * substitution positions. Error draws use replacement; clipping adjusts the selected
 * threshold rather than removing trial positions. Endpoint shortcuts use separate bounds.
 * These model results are not a general alignment-detection guarantee, and the endpoint
 * one-hit floor can remain positive when no sampled seed survives.
 * Use instances serially: cache misses are synchronized, but cache reads occur outside
 * that lock (IFA-005). Configure the static iteration count before use; changing it
 * neither invalidates cached thresholds nor synchronizes existing calculations.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date June 4, 2025
 */
public class MinHitsCalculator{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Stores parameters and creates an owned wildcard pattern; simulation is lazy.
	 * The stride is clamped to at least one. Other parameters are stored directly.
	 *
	 * @param k_ K-mer length
	 * @param maxSubs_ Absolute substitution allowance, further limited by query length and minid_
	 * @param minid_ Identity fraction used in the query-length substitution limit
	 * @param midMaskLen_ Number of wildcard bases in the middle of the k-mer
	 * @param minProb_ Histogram upper-tail fraction; 0 requests all sampled starts,
	 * at least 1 uses a loss bound with a one-hit floor, and negative values request one hit
	 * @param maxClip_ Maximum clipping allowed (fraction &lt;1 or absolute ≥1)
	 * @param kStep_ Kmer step size (1 for all kmers, 2 for every other kmer, etc.)
	 */
	public MinHitsCalculator(int k_, int maxSubs_, float minid_, int midMaskLen_, float minProb_, float maxClip_, int kStep_){
		k=k_;
		maxSubs0=maxSubs_;
		minid=minid_;
		midMaskLen=midMaskLen_;
		minProb=minProb_;
		maxClipFraction=maxClip_;
		kStep=Math.max(1, kStep_);

		// Pre-compute wildcard pattern for efficient simulation
		wildcards=makeWildcardPattern(k, midMaskLen);

		// Retained legacy masks; the model uses wildcards instead.
		kMask=~((-1)<<(2*k));

		// Store the central-bit mask in its original initialization position.
		int bitsPerBase=2;
		int bits=midMaskLen*bitsPerBase;
		int shift=((k-midMaskLen)/2)*bitsPerBase;
		midMask=~((~((-1)<<bits))<<shift);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Builds an owned central-wildcard pattern; true positions ignore substitutions.
	 * Unequal flanks place the block toward the lower-index side.
	 * @param k K-mer length
	 * @param midMaskLen Count of wildcard bases
	 * @return Boolean array with true for wildcard positions
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
	 * Counts k-mers unaffected by errors, honoring wildcard positions.
	 *
	 * @param errors Error positions; read only
	 * @param wildcards Wildcard positions within each window; read only
	 * @param queryLen Logical query length, independent of BitSet capacity
	 * @param step Positive sampling stride; the first sampled window starts at zero
	 * @return Number of surviving sampled windows
	 */
	private int countErrorFreeKmers(BitSet errors, boolean[] wildcards, int queryLen, int step){
		int count=0;

		// Check every step-th k-mer position in query
		for(int i=0; i<=queryLen-k; i+=step){
			boolean errorFree=true;

			// Check each position within this k-mer
			for(int j=0; j<k && errorFree; j++){
				errorFree=wildcards[j] || (!errors.get(i+j));
			}
			if(errorFree){count++;}
		}
		return count;
	}

	/**
	 * Calculates a raw threshold for a contiguous query of validKmers+k-1 bases.
	 * Effective substitutions are the lesser of maxSubs0 and the truncated identity
	 * allowance; clipping is a truncated fraction of that query length or an absolute count.
	 * For an interior probability, draws error positions with replacement, counts surviving
	 * sampled windows, and walks the histogram from high counts to low until the upper tail
	 * reaches the truncated iterations*minProb target. The selected count is capped by
	 * validKmers-maxSubs-maxClips; trial masks themselves are not clipped.
	 * Endpoint and fallback branches retain their own floors and loss expressions.
	 * @param validKmers Number of valid windows, represented as contiguous in the simulation
	 * @return Raw model threshold; minHits clamps negative results to zero
	 */
	private int simulate(int validKmers){
		// Calculate effective clipping limit for this query length
		int queryLen=validKmers+k-1;
		final int maxSubs=Math.min(maxSubs0, (int)(queryLen*(1-minid)));
		int maxClips=(maxClipFraction<1 ? (int)(maxClipFraction*queryLen) : (int)maxClipFraction);
		
		//FIXED IFA-011: endpoint thresholds use the sampled-window units of the simulation.
		if(minProb>=1){
			int unmasked=(Tools.max(2, k-midMaskLen));// Number of kmers impacted by a sub
			final int sampled=(int)((validKmers+(long)kStep-1)/kStep);
			return Math.max(1, sampled-(unmasked*maxSubs)-maxClips);
		}else if(minProb==0){
			return (int)((validKmers+(long)kStep-1)/kStep);
		}else if(minProb<0){
			return 1;
		}
		
		// Build histogram of surviving k-mer counts
		int[] histogram=new int[validKmers+1];
		BitSet errors=new BitSet(queryLen);// Reused for all trials

		// Run Monte Carlo simulation
		for(int iter=0; iter<iterations; iter++){
			errors.clear();

			// Duplicate draws collapse to one error position in the BitSet.
			for(int i=0; i<maxSubs; i++){
				int pos=randy.nextInt(queryLen);
				errors.set(pos);
			}

			// Count k-mers that survive the errors
			int errorFreeKmers=countErrorFreeKmers(errors, wildcards, queryLen, kStep);
			histogram[errorFreeKmers]++;

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

		// Truncate the requested upper-tail fraction to a trial count.
		int targetCount=(int)(iterations*minProb);
		int cumulative=0;

		// Walk down from highest hit count to find percentile threshold
		for(int hits=validKmers; hits>=0; hits--){
			cumulative+=histogram[hits];
			if(cumulative>=targetCount){
				// Preserve the model's post-selection substitution/clipping cap.
				return Math.min(hits, validKmers-maxSubs-maxClips);
			}
		}

		return Math.max(1, validKmers-maxSubs0-maxClips);// Legacy fallback uses the absolute substitution allowance
	}

	/**
	 * Returns a cached threshold, computing and storing a nonnegative result on a miss.
	 * Cache keys contain only the window count; changing iterations leaves old entries intact.
	 * Miss order consumes the instance random stream and can affect later estimates.
	 * The synchronized miss path does not protect the initial unlocked cache reads.
	 * @param validKmers Valid-window count for the contiguous model
	 * @return Cached nonnegative model threshold, not a universal detection guarantee
	 */
	public int minHits(int validKmers){
		//TODO: IFA-005 - cache reads here are outside the write lock; concurrent callers could race map mutation.
		//The standalone main is serial; no production concurrent caller has been established for this version.
		int minHits=validKmerToMinHits.get(validKmers);
		if(minHits<0 && !validKmerToMinHits.contains(validKmers)){
			synchronized(validKmerToMinHits){
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

	/** Bases per seed window. */
	private final int k;
	/** Absolute substitution allowance before applying the query-length identity limit. */
	private final int maxSubs0;
	/** Identity fraction used to derive the effective substitution allowance. */
	private final float minid;
	/** Number of central wildcard positions. */
	private final int midMaskLen;
	/** Clipping setting: fraction below 1, otherwise absolute bases. */
	private final float maxClipFraction;
	/** Stored legacy window mask; not read by the simulation. */
	private final int kMask;
	/** Stored legacy central-bit mask; the simulation reads wildcards instead. */
	private final int midMask;
	/** Histogram upper-tail fraction, with separate endpoint and negative shortcuts. */
	private final float minProb;
	/** Positive sampling stride, clamped during construction. */
	final int kStep;
	/** Owned wildcard pattern, read without modification after construction. */
	private final boolean[] wildcards;
	/** Window-count cache; misses lock this map, but initial reads remain unlocked. */
	private final IntHashMap2 validKmerToMinHits=new IntHashMap2();
	/** Instance generator seeded with 1; the current Shared factory constructs a new generator. */
	private final Random randy=Shared.threadLocalRandom(1);
	/** Trial count for future cache misses; changes do not invalidate existing cache entries. */
	public static int iterations=100000;
	/** Compile-time switch for per-trial diagnostics, not a mutable command-line option. */
	private static final boolean verbose=false;

	/*--------------------------------------------------------------*/
	/*----------------        Debug Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Runs a standalone model diagnostic, writing parameters, mask, timings and threshold to stderr.
	 * Example arguments: k=13 validkmers=50 maxsubs=5 minid=0.9 midmask=1
	 * minprob=0.99 maxclip=0.25 kstep=1 iterations=10000.
	 * Passing verbose triggers an assertion when assertions are enabled; the debug switch is final.
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

		System.err.println("MinHitsCalculator testing:");
		System.err.println("  k="+k+" validKmers="+validKmers+" maxSubs="+maxSubs+" minid="+minid);
		System.err.println("  midMaskLen="+midMaskLen+" minProb="+minProb+" maxClip="+maxClip+" kStep="+kStep);
		System.err.println("  iterations="+iterations+" verbose="+verbose);

		Timer t=new Timer();
		MinHitsCalculator mhc=new MinHitsCalculator(k, maxSubs, minid, midMaskLen, minProb, maxClip, kStep);
		t.stopAndPrint();
		System.err.println("Wildcard mask (boolean[]):");
		for(int i=0; i<mhc.wildcards.length; i++){
			System.err.print(mhc.wildcards[i] ? "W" : "m");
		}
		System.err.println();

		t.start();
		int minHits=mhc.minHits(validKmers);
		t.stopAndPrint();

		System.err.println("\nResult: minHits="+minHits);
	}

	/**
	 * Prints logical query positions as m (match) or S (substitution) when compiled verbose.
	 * Does nothing under the current false verbose constant.
	 * @param errors Error positions; read only
	 * @param queryLen Logical query length
	 */
	static void printSequence(BitSet errors, int queryLen){
		if(!verbose){return;}
		StringBuilder sb=new StringBuilder(queryLen);
		for(int i=0; i<queryLen; i++){
			sb.append(errors.get(i) ? 'S' : 'm');
		}
		System.err.println(sb.toString());
	}

	/**
	 * Prints nonzero histogram bins to stderr when compiled verbose; currently does nothing.
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
