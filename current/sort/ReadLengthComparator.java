package sort;

import stream.Read;

/**
 * Comparator for sorting reads by length, with longest reads first by default.
 * Uses multi-level comparison: read length, mate length, string ID, then numeric ID.
 * @author Brian Bushnell
 * @date Jul 19, 2013
 */
public final class ReadLengthComparator extends ReadComparator {
	
	/**
	 * Private constructor to enforce singleton usage via the static comparator instance.
	 */
	private ReadLengthComparator(){}
	
	/**
	 * Compares two reads for sorting by length.
	 * Uses hierarchical comparison: primary read length, mate read length,
	 * string ID, then numeric ID as tiebreakers.
	 *
	 * @param a First read to compare
	 * @param b Second read to compare
	 * @return Negative if a < b, positive if a > b, zero if equal
	 */
	@Override
	public int compare(Read a, Read b) {
		int x=compareInner(a, b);
		if(x==0){x=compareInner(a.mate, b.mate);}
		if(x==0){x=a.id.compareTo(b.id);}
		if(x==0){x=a.numericID>b.numericID ? 1 : a.numericID<b.numericID ? -1 : 0;}
		return ascending*x;
	}

	/**
	 * Compares two individual reads by length only, treating nulls as greater than non-nulls.
	 * @param a First read to compare (may be null)
	 * @param b Second read to compare (may be null)
	 * @return Length difference (a.length - b.length), or null ordering
	 */
	private static int compareInner(Read a, Read b) {
		if(a==b){return 0;}
		if(a==null){return 1;}
		if(b==null){return -1;}
		int x=a.length()-b.length();
		return x;
	}
	
	/** Singleton instance for length-based read comparison. */
	public static final ReadLengthComparator comparator=new ReadLengthComparator();
	
	/**
	 * Sort direction multiplier: -1 for descending (default longest first), 1 for ascending.
	 */
	private int ascending=-1;
	
	/** Sets the sort order for length comparison.
	 * @param asc true for ascending (shortest first), false for descending (longest first) */
	@Override
	public void setAscending(boolean asc){
		ascending=(asc ? 1 : -1);
	}
	
	@Override
	public final int ascendingMult() {return ascending;}
	@Override
	public final String name() {return "Length";}
	
}
