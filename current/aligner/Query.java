package aligner;

import java.util.concurrent.atomic.AtomicLong;

import dna.AminoAcid;
import shared.Tools;

/**
 * Represents a query sequence for alignment with k-mer indexing.
 * Uses static configuration arrays to self-select optimal k-mer length.
 * Pre-computes forward and reverse k-mer indices plus metadata for efficient alignment.
 * @author Brian Bushnell
 * @contributor Isla, Amber
 * @date June 3, 2025
 */
public class Query {

	/**
	 * Constructs a Query object.
	 * Automatically selects the best k-mer length from the static calculators array.
	 *
	 * @param name_ Sequence name
	 * @param nid Numeric ID
	 * @param bases_ Sequence bases
	 * @param quals_ Sequence qualities
	 */
	public Query(String name_, long nid, byte[] bases_, byte[] quals_) {
		name=name_; 
		numericID=nid; 
		bases=bases_; 
		quals=quals_;
		rbases=AminoAcid.reverseComplementBases(bases);

		if(calculators!=null){
			// Self-select the best calculator
			int bestIndex=calculators.length-1; // Default to shortest K / last calculator
			for(int i=0; i<calculators.length; i++){
				MinHitsCalculator2 mhc=calculators[i];
				int vk=bases.length-mhc.k+1;
				if(vk>0 && mhc.minHits(vk)>=minSeedHits){
					bestIndex=i;
					break;
				}
			}

			// Apply configuration
			calculatorIndex=bestIndex;
			MinHitsCalculator2 mhc=calculators[bestIndex];
			k=mhc.k;
			midMaskLen=mhc.midMaskLen;
			midMask=makeMidMask(k, midMaskLen);
			maxClipFraction=mhc.maxClipFraction;

			// Build Index
			int[][] index=makeIndex(bases, k, midMask);
			kmers=index[0]; 
			rkmers=index[1];

			// Calculate Hit Requirements
			validKmers=(kmers==null ? 0 : Tools.countGreaterThan(kmers, -1));
			if(validKmers>0){
				minHits=mhc.minHits(validKmers);
				maxMisses=(validKmers/mhc.kStep)-minHits;
			}else{
				minHits=0;
				maxMisses=Integer.MAX_VALUE;
			}
		}else{
			// Fallback for non-indexed queries or uninitialized statics
			calculatorIndex=-1;
			k=0; midMaskLen=0; midMask=0; maxClipFraction=maxClip;
			kmers=rkmers=null;
			validKmers=0; minHits=0; maxMisses=Integer.MAX_VALUE;
		}

		// Calculate derived constants
		maxClips=(maxClipFraction<1 ? (int)(maxClipFraction*bases.length) : (int)maxClipFraction);
	}

	public int length(){return bases.length;}

	/**
	 * Creates k-mer indices for forward and reverse orientations.
	 */
	private static int[][] makeIndex(byte[] sequence, int k, int mask){
		if(sequence.length<k || k<1){return blankIndex;}

		final int shift=2*k, shift2=shift-2, bitMask=~((-1)<<shift);
		int kmer=0, rkmer=0, len=0;
		int[][] ret=new int[2][sequence.length-k+1];
		int kmerCount=0;

		for(int i=0, idx=-k+1; i<sequence.length; i++, idx++){
			final byte b=sequence[i];
			final int x=AminoAcid.baseToNumber[b], x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&bitMask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&bitMask;

			if(x<0){len=0; rkmer=0;}else{len++;}
			if(idx>=0){
				if(len>=k && !AminoAcid.isHomopolymer(kmer, k, blacklistRepeatLength)){
					ret[0][idx]=(kmer&mask);
					ret[1][idx]=(rkmer&mask); 
					kmerCount++;
				}else{ret[0][idx]=ret[1][idx]=-1;}
			}
		}
		Tools.reverseInPlace(ret[1]);
		return kmerCount<1 ? blankIndex : ret;
	}

	public static int makeMidMask(int k, int maskLen){
		if(maskLen<1){return -1;}
		int bitsPerBase=2;
		assert(k>maskLen+1);
		int bits=maskLen*bitsPerBase, shift=((k-maskLen)/2)*bitsPerBase;
		int middleMask=~((~((-1)<<bits))<<shift);
		return middleMask;
	}

	/** Must be called before creating Query objects */
	public static void setCalculators(MinHitsCalculator2[] calcs, int minSeedHits_){
		calculators=calcs;
		minSeedHits=minSeedHits_;
		if(calcs!=null && calcs.length>0){
			kArray=new int[calcs.length];
			midMaskArray=new int[calcs.length];
			for(int i=0; i<calcs.length; i++){
				kArray[i]=calcs[i].k;
				midMaskArray[i]=calcs[i].midMaskLen;
			}
		}else {
			kArray=midMaskArray=null;
		}
	}

	// Fields
	public final String name;
	public final long numericID;
	public final byte[] bases;
	public final byte[] rbases;
	public final byte[] quals;

	// Indexing Fields
	public final int calculatorIndex; // Points to the index in the static calculators array
	public final int k;
	public final int midMaskLen;
	public final int midMask;
	public final int[] kmers;
	public final int[] rkmers;

	// Alignment Parameters
	public final int validKmers;
	public final int minHits;
	public final int maxMisses;
	public final int maxClips;
	public final float maxClipFraction;

	// Statics
	public static int blacklistRepeatLength=2;
	private static final int[][] blankIndex=new int[2][];

	// Static Configuration Arrays
	public static MinHitsCalculator2[] calculators;
	public static int[] kArray;
	public static int[] midMaskArray;
	public static float maxClip=0.25f; // Default fallback
	public static int minSeedHits=1;
	public AtomicLong alignments=new AtomicLong(0);
}