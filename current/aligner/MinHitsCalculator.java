package aligner;

import java.util.BitSet;
import shared.Random;

import shared.Shared;
import shared.Timer;
import shared.Tools;
import map.IntHashMap2;

/**
 * Calculates the minimum seed hits required to detect indel-free alignments at a target probability.
 * Uses Monte Carlo simulation to model wildcards, error patterns, and clipping limits.
 *
 * @author Brian Bushnell
 * @contributor Isla
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
	 * @param step Step size for sampling k-mer positions (models reference indexing step)
	 * @return Number of error-free k-mers
	 */
	private int countErrorFreeKmers(BitSet errors, boolean[] wildcards, int queryLen, int step){
		int count=0;
		
		// Check every step-th k-mer position in query
		for(int i=0; i<=queryLen-k; i+=step){
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
	final int kStep;
	private final boolean[] wildcards;
	private final IntHashMap2 validKmerToMinHits=new IntHashMap2();
	private final Random randy=Shared.threadLocalRandom(1);
	public static int iterations=100000;
	private static final boolean verbose=false; // Set to true for debugging

	/*--------------------------------------------------------------*/
	/*----------------        Debug Methods         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Main method for standalone testing and debugging.
	 * Usage: java MinHitsCalculator verbose=true k=13 validKmers=50 maxsubs=5 minid=0.9 midmask=1 minprob=0.99 maxclip=0.25 kstep=1 iterations=10000
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
	 * Prints a sequence with errors marked. Format: mmmmmSmmmmSmmmm where m=match, S=substitution.
	 * @param errors BitSet of error positions
	 * @param queryLen Length of query sequence
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