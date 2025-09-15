package sort;

import dna.AminoAcid;
import shared.Tools;
import stream.Read;

/**
 * @author Brian Bushnell
 * @date Oct 27, 2014
 *
 */


public class ReadComparatorTopological5Bit extends ReadComparator{
	
	private ReadComparatorTopological5Bit(){}
	
	@Override
	public int compare(Read r1, Read r2) {
		return ascending*compare(r1, r2, true);
	}
	
	public int compare(Read r1, Read r2, boolean compareMates) {
		
		if(r1.numericID!=r2.numericID){return r1.numericID>r2.numericID ? 1 : -1;}
		
		int x=compareVectors(r1.bases, r2.bases, 12);
		if(x!=0){return x;}
		
		if(r1.mate!=null && r2.mate!=null){
			x=compareVectors(r1.mate.bases, r2.mate.bases, 0);
		}
		if(x!=0){return x;}

		if(r1.bases!=null && r2.bases!=null && r1.length()!=r2.length()){return r1.length()-r2.length();}
		if(r1.mate!=null && r2.mate!=null && r1.mate.bases!=null && r2.mate.bases!=null
				&& r1.mateLength()!=r2.mateLength()){return r1.mateLength()-r2.mateLength();}
		
		x=compareVectors(r1.quality, r2.quality, 0);
		if(x!=0){return 0-x;}
		
		if(r1.mate!=null && r2.mate!=null){
			x=compareVectors(r1.mate.quality, r2.mate.quality, 0);
		}
		if(x!=0){return 0-x;}
		
		return r1.id.compareTo(r2.id);
	}
	
	public int compareVectors(final byte[] a, final byte[] b, final int start){
		if(a==null || b==null){
			if(a==null && b!=null){return 1;}
			if(a!=null && b==null){return -1;}
			return 0;
		}
		final int lim=Tools.min(a.length, b.length);
		for(int i=start; i<lim; i++){
			if(a[i]<b[i]){return -1;}
			if(a[i]>b[i]){return 1;}
		}
		return 0;
	}
	
	public int compareVectorsN(final byte[] a, final byte[] b, final int start){
		if(a==null || b==null){
			if(a==null && b!=null){return 1;}
			if(a!=null && b==null){return -1;}
			return 0;
		}
		final int lim=Tools.min(a.length, b.length);
		for(int i=start; i<lim; i++){
			if(a[i]=='N' && b[i]!='N'){return 1;}
			if(a[i]!='N' && b[i]=='N'){return -1;}
			if(a[i]<b[i]){return -1;}
			if(a[i]>b[i]){return 1;}
		}
		return 0;
	}
	
	public static long genKmer(Read r) {
		long kmer=genKmer(r.bases);
		r.numericID=kmer;
		return kmer;
	}
	
	public static long genKmer(byte[] bases){
		final byte[] lookup=AminoAcid.symbolTo5Bit;
		final int k=12;
		final int max=Tools.min(bases.length, 12);
		long kmer=0;
		
		for(int i=0; i<max; i++){
			byte b=bases[i];
			long x=lookup[b];
			kmer=((kmer<<5)|x);
		}
		if(max<k){kmer<<=(5*(k-max));}
		assert(kmer>=0);
		return kmer;
	}

	@Override
	public void setAscending(boolean asc) {
		ascending=(asc ? 1 : -1);
	}
	
	public static final ReadComparatorTopological5Bit comparator=new ReadComparatorTopological5Bit();
	
	int ascending=1;
}
