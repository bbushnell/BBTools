package sort;

import dna.AminoAcid;
import shared.Tools;
import stream.Read;

/**
 * Sorts similarly to Clumpify
 * @author Brian Bushnell
 * @date Oct 5, 2022
 *
 */
public final class ReadComparatorClump extends ReadComparator {
	
	private ReadComparatorClump(){}
	
	@Override
	public int compare(Read a, Read b) {
		int x=compareInner(a, b);
		if(x==0){x=compareInner(a.mate, b.mate);}
		if(x==0){x=a.id.compareTo(b.id);}
		return ascending*x;
	}

	private static int compareInner(Read a, Read b) {
		if(a==b){return 0;}
		if(a==null){return 1;}
		if(b==null){return -1;}
		if(a.numericID!=b.numericID){return a.numericID>b.numericID ? 1 : -1;}
		if(a.strand()!=b.strand()){return a.strand()-b.strand();}
		if(a.start!=b.start){return a.start-b.start;}
		return 0;
	}
	
	/** Finds the global maximum, forward and reverse */
	public static final long set(Read r){
		if(r.length()<k){return setShort(r);}
		
		final byte[] bases=r.bases;
		long kmer=0;
		long rkmer=0;
		int len=0;
		
//		if(bases==null || bases.length<k){return -1;}
		
		long topCode=Long.MIN_VALUE;
		long topKmer=Long.MIN_VALUE;
		int topStrand=0;
		int topStop=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber[b];
			long x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
			
			if(x<0){
				len=0;
			}else{len++;}
			
			if(len>=k){
				final long kmax=Tools.max(kmer, rkmer);
				final long code=Tools.hash64shift(kmax);
				
				if(code>topCode){
					topKmer=kmax;
					topCode=code;
					topStrand=(kmax==kmer ? 0 : 1);
					topStop=i;
				}
			}
		}
		if(topCode==Long.MIN_VALUE){
			return setShort(r);
		}
		r.numericID=topKmer;
		r.setStrand(topStrand);
		r.start=topStop;
		return topKmer;
	}
	
	/** Generates a key when the read is shorter than k */
	public static final long setShort(Read r){
		final byte[] bases=r.bases;
		final int max=Tools.min(bases.length, k);
		long kmer=0;
		long rkmer=0;
		
		for(int i=0; i<max; i++){
			byte b=bases[i];
			long x=AminoAcid.baseToNumber0[b];
			long x2=AminoAcid.baseToComplementNumber0[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
		}

		final long kmax=Tools.max(kmer, rkmer);
		r.numericID=kmax;
		r.setStrand((kmax==kmer) ? 0 : 1);
		r.start=max-1;
		return kmax;
	}
	
	@Override
	public void setAscending(boolean asc){
		ascending=(asc ? 1 : -1);
	}
	
//	public void setK(int k_){
//		k=k_;
//		assert(k>0 && k<=32) : k;
//		
//		shift=2*k;
//		shift2=shift-2;
//		mask=(shift>63 ? -1L : ~((-1L)<<shift));
//	}
	
	public static final ReadComparatorClump comparator=new ReadComparatorClump();
	
	private int ascending=-1;
	
	private static final int k=31;
	private static final int shift=2*k;
	private static final int shift2=shift-2;
	private static final long mask=(shift>63 ? -1L : ~((-1L)<<shift));;
	
}
