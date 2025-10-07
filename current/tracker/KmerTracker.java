package tracker;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Tools;
import structures.IntRingBufferCond;

/**
 * Counts short k-mers in sequences with optional fixed-window tracking.
 * Provides compositional metrics using GC-binned normalization for strand-independent
 * dinucleotide analysis. For k=2, calculates frequencies within compositional bins
 * (e.g., 0-GC dimers: AA/AT/TA/TT, 2-GC dimers: CC/CG/GC/GG) for composition-independent
 * sequence signature analysis.
 *
 * @author Brian Bushnell
 * @date October 2, 2025
 */
public class KmerTracker{
	
	public KmerTracker(int k_) {this(k_, 0);}

	/**
	 * Constructs a k-mer tracker with optional windowed counting.
	 * @param k_ K-mer length (must be 1-15)
	 * @param window_ Fixed window size for rolling counts (0 for unlimited)
	 */
	public KmerTracker(int k_, int window_) {
		assert(k_>0 && k_<16);
		k=k_;
		bits=2*k;
		mask=~((-1)<<bits);
		window=window_;
		counts=new long[mask+1];
		buffer=(window>0 ? new IntRingBufferCond(window) : null);
	}

	/**
	 * Counts all k-mers in a sequence.
	 * Resets on ambiguous bases (non-ACGT).
	 * @param bases Sequence to process
	 */
	public void add(final byte[] bases){
		if(bases==null || bases.length<k){return;}

		int kmer=0;
		int len=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumber[b];// Element i is: 0 for 'A', 1 for 'C', 2 for 'G', 3 for 'T', -1 otherwise
			kmer=(((kmer<<2)|x)&mask);
			if(x>=0){
				len++;
				count++;
				if(len>=k) {counts[kmer]++;}
			}else{len=kmer=0;}
		}
	}

	/**
	 * Adds a single base to rolling k-mer counter.
	 * Updates count when k valid bases have been seen.
	 * @param b Base to add
	 */
	public void add(byte b) {
		int x=AminoAcid.baseToNumber[b];// Element i is: 0 for 'A', 1 for 'C', 2 for 'G', 3 for 'T', -1 otherwise
		kmer=(((kmer<<2)|x)&mask);
		if(x>=0){
			len++;
			count++;
			if(len>=k) {counts[kmer]++;}
		}else{len=kmer=0;}
	}

	/**
	 * Adds a base with windowed counting.
	 * Maintains fixed-size window by incrementing new k-mer and decrementing evicted k-mer.
	 * @param b Base to add
	 * @return True if this completed a new valid window
	 */
	public boolean addWindowed(byte b) {
		int x=AminoAcid.baseToNumber[b];
		kmer=(((kmer<<2)|x)&mask);
		if(x>=0){
			len++;
			count++;
			if(len>=k) {
				counts[kmer]++;
				int old=buffer.add(kmer);
				if(old>=0) {counts[old]--;}
				return buffer.isFull();
			}
		}else{len=kmer=0;}
		return false;
	}

	/** Returns GC content from instance counts */
	public float GC() {return GC(counts);}
	/** Returns strand bias metric from instance counts */
	public float strandedness() {return strandedness(counts, k);}
	/** Returns AA+TT fraction within 0-GC dimers from instance counts */
	public float AAAT() {return AAAT(counts);}
	/** Returns CC+GG fraction within 2-GC dimers from instance counts */
	public float CCCG() {return CCCG(counts);}
	/** Returns homopolymer dimer fraction from instance counts */
	public float HH() {return HH(counts);}
	/** Returns purine/pyrimidine dimer fraction from instance counts */
	public float PP() {return PP(counts);}
	/** Returns hydrophobic metric from instance counts */
	public float HMH() {return HMH(counts);}
	/** Returns combined homopolymer and purine/pyrimidine metric from instance counts */
	public float HHPP() {return HHPP(counts);}
	public float ACTG() {return ACTG(counts);}
	public float ACAG() {return ACAG(counts);}
	public float CAGA() {return CAGA(counts);}
	public float CCMCG() {return CCMCG(counts);}
	public float ATMTA() {return ATMTA(counts);}
	public float AT() {return AT(counts);}

	/**
	 * Calculates GC content from k-mer counts.
	 * Works for any k by examining terminal base.
	 * @param counts K-mer count array
	 * @return GC fraction (0.0-1.0)
	 */
	public static float GC(long[] counts) {//Works for any k
		final int mask=0b11;
		long[] acgt=new long[4];
		for(int kmer=0; kmer<counts.length; kmer++) {
			final long count=counts[kmer];
			final int masked=kmer&mask;
			acgt[masked]+=count;
		}
		long gc=acgt[1]+acgt[2];
		long at=acgt[0]+acgt[3];
		return gc/(float)(at+gc);
	}

	/**
	 * Calculates strand bias as deviation from strand symmetry.
	 * Compares each k-mer to its reverse complement.
	 * @param counts Dimer count array (must be length 16)
	 * @param k K-mer length
	 * @return Strand bias metric (0.0=perfect symmetry, 1.0=maximum bias)
	 */
	public static float strandedness(long[] counts, int k) {//I assume k must be 2, need to check
		assert(counts.length==16);
		final int mask=~((-1)<<(2*k));
		assert(mask==counts.length-1);
		long lower=0, upper=0;
		for(int kmer=0, limit=counts.length/2; kmer<limit; kmer++) {
			long a=counts[kmer];
			long b=counts[mask&(~kmer)];
			lower+=Math.min(a, b);
			upper+=Math.max(a, b);
		}
		//        return lower/(float)(Long.max(1, upper));//Old scale
		return (2*upper/(float)(upper+lower))-1;//Nice 0-1 scale.
	}

	/**
	 * Calculates AA+TT fraction within 0-GC dinucleotides.
	 * Uses GC-binned normalization for composition-independent metric.
	 * @param counts Dimer count array (must be length 16)
	 * @return (AA+TT)/(AA+AT+TA+TT)
	 */
	public static float AAAT(long[] counts) {
		assert(counts.length==16);
		long AA=counts[0b0000], TT=counts[0b1111];
		long AT=counts[0b0011], TA=counts[0b1100];
		return (AA+TT)/(float)(AA+TT+AT+TA);
	}

	/**
	 * Calculates AT-TA signal within 0-GC dinucleotides.
	 * Uses GC-binned normalization for composition-independent metric.
	 * @param counts Dimer count array (must be length 16)
	 * @return 0.5f*(1+((AT-TA)/(float)(AA+TT+AT+TA)));
	 */
	public static float ATMTA(long[] counts) {
		assert(counts.length==16);
		long AA=counts[0b0000], TT=counts[0b1111];
		long AT=counts[0b0011], TA=counts[0b1100];
		return 0.5f*(1+((AT-TA)/(float)(AA+TT+AT+TA)));
	}

	/**
	 * Calculates AT ratio within 0-GC dinucleotides.
	 * Uses GC-binned normalization for composition-independent metric.
	 * @param counts Dimer count array (must be length 16)
	 * @return (AT)/(AA+AT+TA+TT)
	 */
	public static float AT(long[] counts) {
		assert(counts.length==16);
		long AA=counts[0b0000], TT=counts[0b1111];
		long AT=counts[0b0011], TA=counts[0b1100];
		return AT/(float)(AA+TT+AT+TA);
	}
	
	/**
	 * Calculates CC+GG fraction within 2-GC dinucleotides.
	 * Uses GC-binned normalization for composition-independent metric.
	 * @param counts Dimer count array (must be length 16)
	 * @return (CC+GG)/(CC+CG+GC+GG)
	 */
	public static float CCCG(long[] counts) {
		assert(counts.length==16);
		long CC=counts[0b0101], GG=counts[0b1010];
		long CG=counts[0b0110], GC=counts[0b1001];
		return (CC+GG)/(float)(CC+GG+CG+GC);
	}
	
	public static float CCMCG(long[] counts) {
		assert(counts.length==16);
		long CC=counts[0b0101], GG=counts[0b1010];
		long CG=counts[0b0110], GC=counts[0b1001];
		return 0.5f*(1+(CC+GG-CG)/(float)(CC+GG+CG+GC));
	}

	/**
	 * Calculates homopolymer dimer fraction (AA, CC, GG, TT).
	 * @param counts Dimer count array (must be length 16)
	 * @return Homopolymer fraction of all dimers
	 */
	public static float HH(long[] counts) {
		assert(counts.length==16);
		long AA=counts[0b0000], TT=counts[0b1111];
		long AT=counts[0b0011], TA=counts[0b1100];
		long CC=counts[0b0101], GG=counts[0b1010];
		long CG=counts[0b0110], GC=counts[0b1001];
		return (AA+CC+GG+TT)/(float)(AA+TT+AT+TA+CC+GG+CG+GC);
	}

	/**
	 * Calculates homopolymer dimer fraction from int counts.
	 * @param counts Dimer count array (must be length 16)
	 * @return Homopolymer fraction of all dimers
	 */
	public static float HH(int[] counts) {
		assert(counts.length==16);
		long AA=counts[0b0000], TT=counts[0b1111];
		long AT=counts[0b0011], TA=counts[0b1100];
		long CC=counts[0b0101], GG=counts[0b1010];
		long CG=counts[0b0110], GC=counts[0b1001];
		return (AA+CC+GG+TT)/(float)(AA+TT+AT+TA+CC+GG+CG+GC);
	}

	/**
	 * Calculates purine-purine/pyrimidine-pyrimidine dimer fraction.
	 * Purines: A (00), G (10); Pyrimidines: C (01), T (11).
	 * Tests second bit of each base (0=purine, 1=pyrimidine).
	 * @param counts Dimer count array (must be length 16)
	 * @return Fraction of dimers with matching purine/pyrimidine type
	 */
	public static float PP(long[] counts) {
		assert(counts.length==16);
		//Purine: A=00, G=10
		//Pyramidine: C=01, T=11
		final int mask=0b0101;
		final int purine=0;
		final int pyramidine=0b0101;
		long purineCount=0;
		long pyramidineCount=0;
		long deltaCount=0;
		for(int kmer=0; kmer<counts.length; kmer++) {
			final long count=counts[kmer];
			final int masked=kmer&mask;
			if(masked==purine) {purineCount+=count;}
			else if(masked==pyramidine) {pyramidineCount+=count;}
			else {deltaCount+=count;}
		}
		long pp=purineCount+pyramidineCount;
		return (pp)/(float)(pp+deltaCount);
	}
	
	public static float ACTG(long[] counts) {
		assert(counts.length==16);
		long AC=counts[0b0001], TG=counts[0b1110];
		long AG=counts[0b0010], CT=counts[0b0111];
		long TC=counts[0b1101], GA=counts[0b1000];
		long GT=counts[0b1011], CA=counts[0b0100];
		return (AC+TG+GT+CA)/(float)(AC+AG+CA+GA+TC+TG+CT+GT);
	}
	
	public static float ACAG(long[] counts) {
		assert(counts.length==16);
		long AC=counts[0b0001], TG=counts[0b1110];
		long AG=counts[0b0010], CT=counts[0b0111];
		long TC=counts[0b1101], GA=counts[0b1000];
		long GT=counts[0b1011], CA=counts[0b0100];
		return 0.5f*(1+(AC+GT-AG-CT)/(float)(AC+AG+CA+GA+TC+TG+CT+GT));
	}
	
	public static float CAGA(long[] counts) {
		assert(counts.length==16);
		long AC=counts[0b0001], TG=counts[0b1110];
		long AG=counts[0b0010], CT=counts[0b0111];
		long TC=counts[0b1101], GA=counts[0b1000];
		long GT=counts[0b1011], CA=counts[0b0100];
		return 0.5f*(1+(CA+TG-GA-TC)/(float)(AC+AG+CA+GA+TC+TG+CT+GT));
	}
	
	public static float CAGA(int[] counts) {
		assert(counts.length==16);
		long AC=counts[0b0001], TG=counts[0b1110];
		long AG=counts[0b0010], CT=counts[0b0111];
		long TC=counts[0b1101], GA=counts[0b1000];
		long GT=counts[0b1011], CA=counts[0b0100];
		return 0.5f*(1+(CA+TG-GA-TC)/(float)(AC+AG+CA+GA+TC+TG+CT+GT));
	}

	/**
	 * Calculates hydrophobic metric from AT-rich vs GC-rich homopolymers.
	 * Normalized to 0.0-1.0 range.
	 * @param counts Dimer count array (must be length 16)
	 * @return Hydrophobic metric
	 */
	public static float HMH(long[] counts) {
		return Math.max(0, 0.5f*(AAAT(counts)-CCCG(counts)+1));
	}

	/**
	 * Calculates combined homopolymer and purine/pyrimidine metric.
	 * @param counts Dimer count array (must be length 16)
	 * @return Average of HH and PP metrics
	 */
	public static float HHPP(long[] counts) {
		return 0.5f*(HH(counts)+PP(counts));
	}

	/**
	 * Adds counts from another count array.
	 * @param counts2 Count array to add
	 */
	public void add(final long[] counts2){Tools.add(counts, counts2);}

	/**
	 * Adds counts from another KmerTracker.
	 * @param tracker Tracker to merge
	 */
	public void add(KmerTracker tracker){add(tracker.counts);}

	/** Resets rolling k-mer state without clearing accumulated counts */
	public void reset() {len=kmer=0;}

	/** Clears all state and accumulated counts */
	public void clearAll() {
		count=len=kmer=0;
		Arrays.fill(counts, 0);
		if(buffer!=null) {buffer.clear();}
	}
	
	public long count() {return count;}
	public void resetCount() {count=0;}

	/** Current rolling k-mer value */
	private int kmer=0;
	/** Current run length of valid bases */
	private long len=0;
	/** Monotonic counter, not reset by Ns */
	private long count=0;

	/** K-mer length */
	public final int k;
	/** Bits per k-mer (2*k) */
	public final int bits;
	/** Bitmask for k-mer extraction */
	public final int mask;
	/** Window size for rolling counts (0 for unlimited) */
	public final int window;

	/** K-mer counts array (length = 4^k) */
	public final long[] counts;
	/** Ring buffer for windowed counting (null if unlimited) */
	public final IntRingBufferCond buffer;

}
