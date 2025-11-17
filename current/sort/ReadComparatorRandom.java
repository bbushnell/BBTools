package sort;

import stream.Read;

/**
 * @author Brian Bushnell
 * @date Mar 6, 2017
 *
 */

public final class ReadComparatorRandom extends ReadComparator{
	
	@Override
	public int compare(Read r1, Read r2) {
		return compareInner(r1, r2)*mult;
	}
	
	/**
	 * Performs the core random value comparison between two reads.
	 * Compares the rand field of each read without applying sort direction.
	 *
	 * @param r1 First read to compare
	 * @param r2 Second read to compare
	 * @return -1 if r1.rand < r2.rand, 1 if r1.rand > r2.rand, 0 if equal
	 */
	public static int compareInner(Read r1, Read r2) {
		if(r1.rand<r2.rand){return -1;}
		if(r1.rand>r2.rand){return 1;}
		return 0;
	}
	
	/** Singleton instance of the random read comparator */
	public static final ReadComparatorRandom comparator=new ReadComparatorRandom();

	@Override
	public void setAscending(boolean asc) {
		mult=asc ? 1 : -1;
	}
	
	/** Sort direction multiplier: 1 for ascending, -1 for descending */
	private int mult=1;
	
}
