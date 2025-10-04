package tracker;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Tools;
import structures.IntRingBufferCond;

/** 
 * Counts short kmers.
 * Optionally uses a ringbuffer to efficiently maintain a fixed-length window.
 * Provides some statistical methods for calculating metrics on dimers.
 * 
 * @author Brian Bushnell
 * @date October 2, 2025
 */
public class KmerTracker{
	
	public KmerTracker(int k_, int window_) {
		assert(k_>0 && k_<16);
		k=k_;
		bits=2*k;
		mask=~((-1)<<bits);
		window=window_;
		counts=new long[mask+1];
		buffer=(window>0 ? new IntRingBufferCond(window) : null);
	}
	
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
				if(len>=k) {counts[kmer]++;}
			}else{len=kmer=0;}
		}
	}
	
	public void add(byte b) {
		int x=AminoAcid.baseToNumber[b];// Element i is: 0 for 'A', 1 for 'C', 2 for 'G', 3 for 'T', -1 otherwise
		kmer=(((kmer<<2)|x)&mask);
		if(x>=0){
			len++;
			if(len>=k) {counts[kmer]++;}
		}else{len=kmer=0;}
	}
	
	/** @return True if this made a new valid window */
	public boolean addWindowed(byte b) {
		int x=AminoAcid.baseToNumber[b];
		kmer=(((kmer<<2)|x)&mask);
		if(x>=0){
			len++;
			if(len>=k) {
				counts[kmer]++;
				int old=buffer.add(kmer);
				if(old>=0) {counts[old]--;}
				return buffer.isFull();
			}
		}else{len=kmer=0;}
		return false;
	}
	
	public float GC() {return GC(counts);}
	public float strandedness() {return strandedness(counts, k);}
	public float AAAT() {return AAAT(counts);}
	public float CCCG() {return CCCG(counts);}
	public float HH() {return HH(counts);}
	public float PP() {return PP(counts);}
	public float HMH() {return HMH(counts);}
	public float HHPP() {return HHPP(counts);}
	
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
//		return lower/(float)(Long.max(1, upper));//Old scale
		return (2*upper/(float)(upper+lower))-1;//Nice 0-1 scale.
	}
	
	public static float AAAT(long[] counts) {
		assert(counts.length==16);
		long AA=counts[0b0000], TT=counts[0b1111];
		long AT=counts[0b0011], TA=counts[0b1100];
		return (AA+TT)/(float)(AA+TT+AT+TA);
	}
	
	public static float CCCG(long[] counts) {
		assert(counts.length==16);
		long CC=counts[0b0101], GG=counts[0b1010];
		long CG=counts[0b0110], GC=counts[0b1001];
		return (CC+GG)/(float)(CC+GG+CG+GC);
	}
	
	public static float HH(long[] counts) {
		assert(counts.length==16);
		long AA=counts[0b0000], TT=counts[0b1111];
		long AT=counts[0b0011], TA=counts[0b1100];
		long CC=counts[0b0101], GG=counts[0b1010];
		long CG=counts[0b0110], GC=counts[0b1001];
		return (AA+CC+GG+TT)/(float)(AA+TT+AT+TA+CC+GG+CG+GC);
	}
	
	public static float HH(int[] counts) {
		assert(counts.length==16);
		long AA=counts[0b0000], TT=counts[0b1111];
		long AT=counts[0b0011], TA=counts[0b1100];
		long CC=counts[0b0101], GG=counts[0b1010];
		long CG=counts[0b0110], GC=counts[0b1001];
		return (AA+CC+GG+TT)/(float)(AA+TT+AT+TA+CC+GG+CG+GC);
	}
	
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
	
	public static float HMH(long[] counts) {
		return 0.5f*(AAAT(counts)-CCCG(counts)+1);
	}
	
	public static float HHPP(long[] counts) {
		return 0.5f*(HH(counts)+PP(counts));
	}
	
	public void add(final long[] counts2){Tools.add(counts, counts2);}
		
	public void add(KmerTracker tracker){add(tracker.counts);}

	void reset() {kmer=len=0;}

	void clearAll() {
		kmer=len=0;
		Arrays.fill(counts, 0);
	}
	
	private int kmer=0;
	private int len=0;
	
	public final int k;
	public final int bits;
	public final int mask;
	public final int window;
	
	public final long[] counts;
	public final IntRingBufferCond buffer;
	
}
