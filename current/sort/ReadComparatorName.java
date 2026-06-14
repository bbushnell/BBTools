package sort;

import stream.Read;


/**
 * Comparator for sorting reads by sequence name/identifier.
 * Compares reads primarily by their ID string, with pair number as tiebreaker.
 * Handles null IDs by treating them as less than non-null IDs.
 * @author Brian Bushnell
 */
public final class ReadComparatorName extends ReadComparator {
	
	/** Private constructor; use the ascending/descending singletons. */
	private ReadComparatorName(int mult_){mult=mult_;}
	
	/**
	 * Compares two reads for sorting by name/ID.
	 * Applies ascending/descending order based on current setting.
	 *
	 * @param r1 First read to compare
	 * @param r2 Second read to compare
	 * @return Negative if r1 < r2, zero if equal, positive if r1 > r2
	 */
	@Override
	public int compare(Read r1, Read r2) {
		int x=compareInner(r1, r2);
		return mult*x;
	}
	
	/**
	 * Core comparison logic for sorting reads by name, treating null IDs as less than non-null and using pair number as a tiebreaker.
	 * Always returns ascending-order result; caller applies sort direction.
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
	
	/** Sort direction multiplier: +1 for ascending, -1 for descending (immutable). */
	private final int mult;

	@Override
	public ReadComparator getComparator(boolean asc){return asc ? ascending : descending;}

	@Override
	public final int ascendingMult() {return mult;}
	@Override
	public final String name() {return "Name";}

	/** Ascending-order singleton. */
	public static final ReadComparatorName ascending=new ReadComparatorName(1);
	/** Descending-order singleton. */
	public static final ReadComparatorName descending=new ReadComparatorName(-1);
	/** Default (ascending) singleton; alias retained for selection/identity call sites. */
	public static final ReadComparatorName comparator=ascending;
	
}
