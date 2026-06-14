package sort;

import dna.AminoAcid;
import shared.Tools;
import stream.Read;



/**
 * Specialized comparator for sorting reads using a 5-bit topological encoding.
 * Compares reads based on k-mer, sequence content, length, and quality.
 * Generates 12-mer using 5-bit nucleotide encoding for efficient comparison.
 *
 * @author Brian Bushnell
 * @date Oct 27, 2014
 */
public class ReadComparatorTopological5Bit extends ReadComparator{

	/** Private constructor; use the ascending/descending singletons. */
	private ReadComparatorTopological5Bit(int mult_){mult=mult_;}

	/**
	 * Compares two reads using topological sorting criteria.
	 * Delegates to the three-parameter compare method with mate comparison enabled.
	 *
	 * @param r1 First read to compare
	 * @param r2 Second read to compare
	 * @return Negative if r1 < r2, positive if r1 > r2, zero if equal
	 */
	@Override
	public int compare(Read r1, Read r2) {
		return mult*compare(r1, r2, true);
	}

	/** Compares two reads by their precomputed 12-mer key (stored in numericID by {@link #genKmer}),
	 * then remaining bases, mate bases, lengths, quality (higher first), and finally string id. */
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

	/** Lexicographically compares two byte arrays from index {@code start} to the shorter length;
	 * null sorts after non-null. */
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

	/** Like {@link #compareVectors} but treats 'N' as greater than any specific base. */
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

	/** Computes the read's leading 12-mer key and stores it in the read's numericID for sorting. */
	public static long genKmer(Read r) {
		long kmer=genKmer(r.bases);
		r.numericID=kmer;
		return kmer;
	}

	/** Builds a 12-symbol, 5-bits-per-symbol key from the leading bases (left-padded if shorter than 12). */
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
	public ReadComparator getComparator(boolean asc){return asc ? ascending : descending;}

	@Override
	public final int ascendingMult() {return mult;}
	@Override
	public final String name() {return "Topology5Bit";}

	/** Sort direction multiplier: +1 for ascending, -1 for descending (immutable). */
	private final int mult;

	/** Ascending-order singleton. */
	public static final ReadComparatorTopological5Bit ascending=new ReadComparatorTopological5Bit(1);
	/** Descending-order singleton. */
	public static final ReadComparatorTopological5Bit descending=new ReadComparatorTopological5Bit(-1);
	/** Default (ascending) singleton; alias retained for selection/identity call sites. */
	public static final ReadComparatorTopological5Bit comparator=ascending;
}
