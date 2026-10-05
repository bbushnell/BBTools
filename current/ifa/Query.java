package ifa;

import java.util.concurrent.atomic.AtomicLong;

import dna.AminoAcid;
import shared.Tools;

/**
 * Query bases, reverse complement, optional k-mer indexes, and alignment counters.
 * Borrows the supplied bases and qualities; owns the reverse complement and nonempty
 * indexes. Do not mutate these arrays after construction: derived values are not rebuilt.
 * Calculator selection follows the supplied order, choosing the first eligible entry
 * or falling back to the last; it is not a search for a globally optimal k.
 * Configure the mutable statics before construction, as the aligners do before starting
 * workers. Configuration changes are not synchronized. The atomic alignment counter
 * alone does not make array or configuration mutation safe during shared use.
 * @author Brian Bushnell
 * @contributor Isla, Amber
 * @date June 3, 2025
 */
public class Query{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Borrows input arrays and builds a reverse complement and optional seed indexes.
	 * With calculators configured, selects the first whose positive potential-window
	 * count produces at least minSeedHits; otherwise uses the last calculator. The
	 * final seed threshold uses the actual valid-window count after indexing.
	 * A null calculator array selects the nonindexed fallback. A nonnull calculator
	 * array must be nonempty and contain only nonnull calculators. Empty-array
	 * rejection remains IFA-003.
	 *
	 * @param name_ Sequence name
	 * @param nid Numeric ID
	 * @param bases_ Nonnull ASCII bases, retained without copying
	 * @param quals_ Optional qualities, retained without copying
	 */
	public Query(String name_, long nid, byte[] bases_, byte[] quals_){
		name=name_;
		numericID=nid;
		bases=bases_;
		quals=quals_;
		rbases=AminoAcid.reverseComplementBases(bases);

		if(calculators!=null){
			// Caller order determines preference; use the last entry if none qualifies.
			int bestIndex=calculators.length-1;
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
				//FIXED IFA-009: each prescan visits a subset of these valid windows, on either strand.
				//Dividing by kStep can undercount query samples and wrongly count reference sampling as query skipping.
				//This conservative budget may scan farther, but final seed counts still enforce minHits.
				maxMisses=validKmers-minHits;
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

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of bases in the borrowed forward sequence. */
	public int length(){return bases.length;}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Builds masked, two-bit windows for each strand, ordered by that strand's start.
	 * Undefined bases and windows rejected by AminoAcid.isHomopolymer yield -1 in
	 * both orientations. Reverse slots are reversed after the forward scan.
	 * @param sequence ASCII bases; read only
	 * @param k Window length; callers use 1 through 15, and values below 1 disable indexing
	 * @param mask Bit mask applied to both orientations after repeat filtering
	 * @return Two owned arrays, or the shared pair of null entries if no valid window exists
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

	/**
	 * Makes a mask clearing maskLen central bases in a two-bit k-mer.
	 * When centering is uneven, the low-bit offset is rounded down.
	 * @param k K-mer length; current callers use 1 through 15
	 * @param maskLen Number of bases to mask; positive values require k greater than maskLen+1
	 * @return All bits set when maskLen is below 1; otherwise the central bits are cleared
	 */
	public static int makeMidMask(int k, int maskLen){
		if(maskLen<1){return -1;}
		int bitsPerBase=2;
		assert(k>maskLen+1);
		int bits=maskLen*bitsPerBase, shift=((k-maskLen)/2)*bitsPerBase;
		int middleMask=~((~((-1)<<bits))<<shift);
		return middleMask;
	}

	/**
	 * Sets global construction parameters and copies calculator k/mask lengths into arrays.
	 * Call before constructing queries, with no concurrent configuration changes.
	 * The calculator array itself is borrowed; existing Query instances are not rebuilt.
	 * @param calcs Null for nonindexed queries, otherwise a nonempty array of nonnull calculators
	 * @param minSeedHits_ Threshold used to select the first eligible calculator
	 */
	public static void setCalculators(MinHitsCalculator2[] calcs, int minSeedHits_){
		//TODO: IFA-003 - empty arrays are retained here but construction selects index -1; no rejection guard yet.
		calculators=calcs;
		minSeedHits=minSeedHits_;
		if(calcs!=null && calcs.length>0){
			kArray=new int[calcs.length];
			midMaskArray=new int[calcs.length];
			for(int i=0; i<calcs.length; i++){
				kArray[i]=calcs[i].k;
				midMaskArray[i]=calcs[i].midMaskLen;
			}
		}else{
			kArray=midMaskArray=null;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Query name used in emitted alignments. */
	public final String name;
	/** Caller-supplied numeric identifier. */
	public final long numericID;
	/** Borrowed forward ASCII bases; mutation invalidates derived indexes and reverse bases. */
	public final byte[] bases;
	/** Owned reverse complement of bases at construction time. */
	public final byte[] rbases;
	/** Borrowed qualities, or null. */
	public final byte[] quals;

	/** Selected construction-time calculator slot; -1 when calculators was null. */
	public final int calculatorIndex;
	/** Selected k-mer length; 0 for nonindexed queries. */
	public final int k;
	/** Selected number of masked central bases; 0 for nonindexed queries. */
	public final int midMaskLen;
	/** Two-bit seed mask; 0 for nonindexed queries. */
	public final int midMask;
	/** Forward-start seed keys, with -1 for invalid windows; null when no usable seeds exist. */
	public final int[] kmers;
	/** Reverse-start seed keys, with -1 for invalid windows; null when no usable seeds exist. */
	public final int[] rkmers;

	/** Number of valid windows before query/reference sampling. */
	public final int validKmers;
	/** Calculator seed threshold; 0 for nonindexed or seedless queries. */
	public final int minHits;
	/** Conservative prescan miss budget from all valid windows; MAX_VALUE for nonindexed/seedless queries. */
	public final int maxMisses;
	/** Whole-base clipping allowance: fractional length product below 1, otherwise absolute count, truncated to int. */
	public final int maxClips;
	/** Selected clipping setting: fraction below 1, otherwise absolute bases despite the field name. */
	public final float maxClipFraction;

	/*--------------------------------------------------------------*/
	/*----------------      Static Configuration    ----------------*/
	/*--------------------------------------------------------------*/

	/** Maximum repeat period passed to AminoAcid.isHomopolymer; below 1 disables that filter. */
	public static int blacklistRepeatLength=2;
	/** Shared no-index sentinel; both entries are null and must remain unchanged. */
	private static final int[][] blankIndex=new int[2][];

	/** Borrowed ordered calculators used by subsequent constructors; null disables indexing. */
	public static MinHitsCalculator2[] calculators;
	/** Snapshot of calculator k values at setCalculators, or null for null/empty input. */
	public static int[] kArray;
	/** Snapshot of calculator mask lengths at setCalculators, or null for null/empty input. */
	public static int[] midMaskArray;
	/** Clipping setting used only when no calculator is selected; fraction below 1 or absolute bases. */
	public static float maxClip=0.25f;
	/** Global construction-time calculator eligibility threshold. */
	public static int minSeedHits=1;
	/** Mutable per-query count; workers increment it to choose the primary alignment. */
	public AtomicLong alignments=new AtomicLong(0);
}
