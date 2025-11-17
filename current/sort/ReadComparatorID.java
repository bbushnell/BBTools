package sort;

import stream.Read;

/**
 * @author Brian Bushnell
 * @date Oct 27, 2014
 *
 */

public final class ReadComparatorID extends ReadComparator{
	
	@Override
	public int compare(Read r1, Read r2) {
		return compareInner(r1, r2)*mult;
	}
	
	/**
	 * Core comparison logic implementing three-level hierarchical sorting.
	 * First compares by numeric ID, then by pair number, finally by string ID lexicographically.
	 * Always returns ascending order result regardless of sort direction setting.
	 *
	 * @param r1 First read to compare
	 * @param r2 Second read to compare
	 * @return Negative if r1 < r2, positive if r1 > r2, zero if equal
	 */
	public static int compareInner(Read r1, Read r2) {
		if(r1.numericID<r2.numericID){return -1;}
		else if(r1.numericID>r2.numericID){return 1;}
		
		int p1=r1.pairnum(), p2=r2.pairnum();
		if(p1<p2){return -1;}
		else if(p1>p2){return 1;}
		
		return r1.id.compareTo(r2.id);
	}

	/** Singleton instance for convenient access to the comparator */
	public static final ReadComparatorID comparator=new ReadComparatorID();

	@Override
	public void setAscending(boolean asc) {
		mult=asc ? 1 : -1;
	}
	
	/**
	 * Multiplier applied to comparison result: 1 for ascending, -1 for descending
	 */
	private int mult=1;
	
}
