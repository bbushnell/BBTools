package sort;

import stream.Read;


/**
 * Hierarchical comparator for sorting reads by numeric and lexicographic identifiers.
 * Implements a three-level sorting strategy: numeric ID (primary), pair number (secondary),
 * and string ID lexicographically (tertiary). Supports both ascending and descending sort orders.
 *
 * @author Brian Bushnell
 * @date Oct 27, 2014
 */
public final class ReadComparatorID extends ReadComparator{
	
	/**
	 * Compares two reads using hierarchical sorting criteria with configurable order.
	 * Applies the comparison multiplier to support ascending/descending sort.
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
		if(r1.numericID<r2.numericID){return -1;}
		else if(r1.numericID>r2.numericID){return 1;}
		
		int p1=r1.pairnum(), p2=r2.pairnum();
		if(p1<p2){return -1;}
		else if(p1>p2){return 1;}
		
		return r1.id.compareTo(r2.id);
	}

	public static final ReadComparatorID comparator=new ReadComparatorID();

	/** Configures the sort order direction.
	 * @param asc true for ascending order, false for descending order */
	@Override
	public void setAscending(boolean asc) {
		mult=asc ? 1 : -1;
	}
	
	private int mult=1;
	
}
