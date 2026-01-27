package sort;

import stream.Read;

/**
 * Sorts reads based on quality and length, prioritizing high-quality and longer reads.
 * Compares Read objects using a multi-stage sorting algorithm that evaluates error rates,
 * total length, and read IDs for consistent ordering.
 *
 * @author Brian Bushnell
 * @date Nov 9, 2016
 */
public final class ReadQualityComparator extends ReadComparator {
	
	private ReadQualityComparator(){}
	
	/**
	 * Compares two reads based on quality and length criteria.
	 * Applies ascending/descending order multiplier to the comparison result.
	 *
	 * @param a First read to compare
	 * @param b Second read to compare
	 * @return Negative if a < b, positive if a > b, zero if equal
	 */
	@Override
	public int compare(Read a, Read b) {
		int x=compareInner(a, b);
		return ascending*x;
	}

	private static int compareInner(Read a, Read b) {
		if(a==b){return 0;}
		if(a==null){return 1;}
		if(b==null){return -1;}
		Read a2=a.mate, b2=b.mate;
		int alen=a.length()+a.mateLength();
		int blen=b.length()+b.mateLength();

		double aerrors=Read.expectedErrors(a.bases, a.quality, true, a.length());
		double berrors=Read.expectedErrors(b.bases, b.quality, true, b.length());
		if(a2!=null){aerrors+=Read.expectedErrors(a2.bases, a2.quality, true, a2.length());}
		if(b2!=null){berrors+=Read.expectedErrors(b2.bases, b2.quality, true, b2.length());}
		
		double arate=aerrors/alen;
		double brate=berrors/blen;
		
		if(arate!=brate){return arate>brate ? 1 : -1;}
		if(alen!=blen){return blen-alen;}
		int x=a.id.compareTo(b.id);
		if(x!=0){return x;}
		if(a.numericID>b.numericID){return 1;}
		if(b.numericID>a.numericID){return -1;}
		return 0;
	}
	
	public static final ReadQualityComparator comparator=new ReadQualityComparator();
	
	private int ascending=1;
	
	/** Sets the sorting order direction.
	 * @param asc true for ascending order, false for descending */
	@Override
	public void setAscending(boolean asc){
		ascending=(asc ? 1 : -1);
	}
	
	@Override
	public final int ascendingMult() {return ascending;}
	@Override
	public final String name() {return "Quality";}
	
}
