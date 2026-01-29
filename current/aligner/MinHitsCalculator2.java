package aligner;

import shared.Random;

import shared.Shared;
import shared.Timer;
import shared.Tools;
import map.IntHashMap2;

/**
 * Calculates the minimum seed hits required to detect indel-free alignments at a target probability.
 * Uses Monte Carlo simulation to model wildcards, error patterns, and clipping limits.
 * Optimized version with branchless operations and direct error storage.
 *
 * @author Brian Bushnell
 * @contributor Noire
 * @date December 30, 2025
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

	/*--------------------------------------------------------------*/
	/*----------------        Debug Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Main method for standalone testing and debugging.
	 * Usage: java MinHitsCalculator2 verbose=true k=13 validKmers=50 maxsubs=5 minid=0.9 midmask=1 minprob=0.99 maxclip=0.25 kstep=1 iterations=10000
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
	 * Prints a sequence with errors marked. Format: mmmmmSmmmmSmmmm where m=match, S=substitution.
	 * @param errors Array of error positions (0 or 1)
	 * @param queryLen Length of query sequence
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
