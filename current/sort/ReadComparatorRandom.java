package sort;

import stream.Read;


/**
 * Comparator for sorting reads based on their random values.
 * Compares reads using their internal rand field for randomized ordering.
 * @author Brian Bushnell
 * @date Mar 6, 2017
 */
public final class ReadComparatorRandom extends ReadComparator{

	/** Private constructor; use the ascending/descending singletons. */
	private ReadComparatorRandom(int mult_){mult=mult_;}

	/**
	 * Compares two reads based on their random values.
	 * Applies the sort direction multiplier to the comparison result.
	 *
	 * @param r1 First read to compare
	 * @param r2 Second read to compare
	 * @return Negative if r1 < r2, positive if r1 > r2, zero if equal
	 */
	@Override
	public int compare(Read r1, Read r2) {
		return compareInner(r1, r2)*mult;
	}

	/** Compares two reads by their assigned random value, producing a shuffled order. */
	public static int compareInner(Read r1, Read r2) {
		if(r1.rand<r2.rand){return -1;}
		if(r1.rand>r2.rand){return 1;}
		return 0;
	}

	@Override
	public ReadComparator getComparator(boolean asc){return asc ? ascending : descending;}

	@Override
	public final int ascendingMult() {return mult;}
	@Override
	public final String name() {return "Random";}

	/** Sort direction multiplier: +1 for ascending, -1 for descending (immutable). */
	private final int mult;

	/** Ascending-order singleton. */
	public static final ReadComparatorRandom ascending=new ReadComparatorRandom(1);
	/** Descending-order singleton. */
	public static final ReadComparatorRandom descending=new ReadComparatorRandom(-1);
	/** Default (ascending) singleton; alias retained for selection/identity call sites. */
	public static final ReadComparatorRandom comparator=ascending;

}
