package sort;

import stream.Read;

/**
 * @author Brian Bushnell
 * @date Oct 27, 2014
 *
 */

public final class ReadComparatorName extends ReadComparator {
	
	/**
	 * Private constructor to prevent direct instantiation. Use static comparator instance.
	 */
	private ReadComparatorName(){}
	
	@Override
	public int compare(Read r1, Read r2) {
		int x=compareInner(r1, r2);
		return ascending*x;
	}
	
	/**
	 * Core comparison logic for sorting reads by name.
	 * Null IDs are considered less than non-null IDs.
	 * For reads with identical names, uses pair number as tiebreaker.
	 *
	 * @param r1 First read to compare
	 * @param r2 Second read to compare
	 * @return Negative if r1 < r2, zero if equal, positive if r1 > r2
	 */
	public static int compareInner(Read r1, Read r2) {
		
		if(r1.id==null && r2.id==null){return r1.pairnum()-r2.pairnum();}
		if(r1.id==null){return -1;}
		if(r2.id==null){return 1;}
		int x=r1.id.compareTo(r2.id);
		if(x==0){return r1.pairnum()-r2.pairnum();}
		return x;
	}
	
	/** Sort direction multiplier: 1 for ascending, -1 for descending */
	private int ascending=1;
	
	@Override
	public void setAscending(boolean asc){
		ascending=(asc ? 1 : -1);
	}

	/** Singleton instance for name-based read comparison */
	public static final ReadComparatorName comparator=new ReadComparatorName();
	
}
