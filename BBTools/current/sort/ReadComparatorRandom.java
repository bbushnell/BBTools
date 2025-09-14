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
