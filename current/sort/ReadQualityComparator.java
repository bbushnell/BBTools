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

	/** Private constructor; use the ascending/descending singletons. */
	private ReadQualityComparator(int mult_){mult=mult_;}

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
		return mult*x;
	}

	/** Compares two reads (plus mates) by combined expected-error rate ascending, then total
	 * length descending, then string id, then numericID. Lower error rate sorts first. */
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

		double arate=aerrors/Math.max(alen, 1);
		double brate=berrors/Math.max(blen, 1);

		if(arate!=brate){return arate>brate ? 1 : -1;}
		if(alen!=blen){return blen-alen;}
		int x=a.id.compareTo(b.id);
		if(x!=0){return x;}
		if(a.numericID>b.numericID){return 1;}
		if(b.numericID>a.numericID){return -1;}
		return 0;
	}

	@Override
	public ReadComparator getComparator(boolean asc){return asc ? ascending : descending;}

	@Override
	public final int ascendingMult() {return mult;}
	@Override
	public final String name() {return "Quality";}

	/** Sort direction multiplier: +1 for ascending, -1 for descending (immutable). */
	private final int mult;

	/** Ascending-order singleton. */
	public static final ReadQualityComparator ascending=new ReadQualityComparator(1);
	/** Descending-order singleton. */
	public static final ReadQualityComparator descending=new ReadQualityComparator(-1);
	/** Default (ascending) singleton; alias retained for selection/identity call sites. */
	public static final ReadQualityComparator comparator=ascending;

}
