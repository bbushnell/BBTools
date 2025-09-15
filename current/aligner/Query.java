package aligner;

import dna.AminoAcid;
import shared.Tools;

/**
 * Represents a query sequence for alignment with k-mer indexing.
 * Optimized for scenarios with small query sets and large reference databases.
 * Pre-computes forward and reverse k-mer indices plus metadata for efficient alignment.
 * 
 * @author Brian Bushnell
 * @contributor Isla  
 * @date June 3, 2025
 */
public class Query {

	/**
	 * Constructs a Query with sequence data and pre-computed indices.
	 * @param name_ Query sequence identifier
	 * @param nid Numeric ID for the query
	 * @param bases_ Forward sequence bases
	 * @param quals_ Quality scores (may be null)
	 */
	public Query(String name_, long nid, byte[] bases_, byte[] quals_) {
		name=name_; numericID=nid; bases=bases_; quals=quals_;
		rbases=AminoAcid.reverseComplementBases(bases); // Pre-compute reverse complement
		int[][] index=makeIndex(bases);
		kmers=index[0]; rkmers=index[1]; // Forward and reverse k-mer arrays
		validKmers=(kmers==null ? 0 : Tools.countGreaterThan(kmers, -1)); // Count non-masked kmers
		minHits=(mhc==null ? 0 : mhc.minHits(validKmers)); // Calculate minimum hits needed
		maxClips=(maxClip<1 ? (int)(maxClip*bases.length) : (int)maxClip); // Max clipping allowed
	}
	
	/** Returns the length of the query sequence */
	public int length() {return bases.length;}

	/**
	 * Creates k-mer indices for forward and reverse orientations.
	 * Applies middle-masking and homopolymer filtering.
	 * @param sequence The sequence to index
	 * @return Array containing [forward_kmers, reverse_kmers]
	 */
	private static int[][] makeIndex(byte[] sequence){
		if(sequence.length<k || !indexKmers){return blankIndex;}
		assert(indexKmers);

		final int shift=2*k, shift2=shift-2, mask=~((-1)<<shift); // Bit manipulation constants
		int kmer=0, rkmer=0, len=0; // Rolling k-mer state
		int[][] ret=new int[2][sequence.length-k+1]; // Forward and reverse arrays
		int kmerCount=0;

		for(int i=0, idx=-k+1; i<sequence.length; i++, idx++){
			final byte b=sequence[i];
			final int x=AminoAcid.baseToNumber[b], x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask; // Roll forward k-mer
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask; // Roll reverse k-mer

			if(x<0){len=0; rkmer=0;}else{len++;} // Reset on ambiguous base
			if(idx>=0) {
				if(len>=k && !AminoAcid.isHomopolymer(kmer, k, blacklistRepeatLength)){
					ret[0][idx]=(kmer&midMask); ret[1][idx]=(rkmer&midMask); // Apply middle mask
					kmerCount++;
				}else {ret[0][idx]=ret[1][idx]=-1;} // Mark invalid k-mers
			}
		}
		Tools.reverseInPlace(ret[1]); // Reverse complement array needs reversal
		return kmerCount<1 ? blankIndex : ret;
	}

	/**
	 * Creates a mask that zeros out middle bases of k-mers for fuzzy matching.
	 * @param k K-mer length
	 * @param maskLen Number of middle bases to mask
	 * @return Bit mask for k-mer filtering
	 */
	public static int makeMidMask(int k, int maskLen) {
		if(maskLen<1) {return -1;}
		int bitsPerBase=2;
		assert(k>maskLen+1);
		int bits=maskLen*bitsPerBase, shift=((k-maskLen)/2)*bitsPerBase;
		int middleMask=~((~((-1)<<bits))<<shift); // Create mask with middle bits zeroed
		return middleMask;
	}
	
	/** Sets all indexing parameters for Query creation */
	public static void setMode(int k_, int midMaskLen_, boolean indexKmers_) {
		setK(k_); setMidMaskLen(midMaskLen_); indexKmers=indexKmers_;
		assert(!indexKmers || (k>0 && k<16 && midMaskLen<k+1));
	}
	
	/** Sets k-mer length for indexing */
	private static void setK(int x) {
		k=Math.max(x, 0); assert(k<16); indexKmers=k>0;
	}
	
	/** Sets middle masking length for fuzzy k-mer matching */
	private static void setMidMaskLen(int x) {
		midMaskLen=Math.max(x, 0);
		assert(midMaskLen>=0 && (midMaskLen<1 || midMaskLen<k-1));
		midMask=makeMidMask(k, midMaskLen); // Recompute mask
	}

	/** Query sequence name/identifier */
	public final String name;
	/** Numeric identifier for the query */
	public final long numericID;
	/** Forward sequence bases */
	public final byte[] bases;
	/** Reverse complement sequence bases */
	public final byte[] rbases;
	/** Quality scores (may be null) */
	public final byte[] quals;
	/** Forward k-mer index array */
	public final int[] kmers;
	/** Reverse k-mer index array */
	public final int[] rkmers;
	/** Number of valid (non-masked) k-mers */
	public final int validKmers;
	/** Minimum hits required for alignment consideration */
	public final int minHits;
	/** Maximum bases that can be clipped */
	public final int maxClips;

	/** K-mer length for indexing */
	private static int k=11;
	/** Length of middle region to mask in k-mers */
	private static int midMaskLen=1;
	/** Bit mask for middle-masking k-mers */
	public static int midMask=makeMidMask(k, midMaskLen);
	/** Whether to create k-mer indices */
	public static boolean indexKmers=true;
	/** Maximum homopolymer length before blacklisting */
	public static int blacklistRepeatLength=2;
	/** Empty index for sequences too short to index */
	private static final int[][] blankIndex=new int[2][];
	/** Calculator for minimum hits based on query length */
	public static MinHitsCalculator mhc;
	/** Maximum fraction/count of bases that can be clipped */
	public static float maxClip=0.25f;
}