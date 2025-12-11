package sort;

import stream.Read;


/**
 * Comparator for sorting reads based on their random values.
 * Compares reads using their internal rand field for randomized ordering.
 * @author Brian Bushnell
 * @date Mar 6, 2017
 */
public final class ReadComparatorRandom extends ReadComparator{
	
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
	
	public static int compareInner(Read r1, Read r2) {
		if(r1.rand<r2.rand){return -1;}
		if(r1.rand>r2.rand){return 1;}
		return 0;
	}
	
	public static final ReadComparatorRandom comparator=new ReadComparatorRandom();

	@Override
	public void setAscending(boolean asc) {
		mult=asc ? 1 : -1;
	}
	
	private int mult=1;
	
}
