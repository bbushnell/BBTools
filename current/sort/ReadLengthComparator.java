package sort;

import stream.Read;

/**
 * Comparator for sorting reads by length, with longest reads first by default.
 * Uses multi-level comparison: read length, mate length, string ID, then numeric ID.
 * @author Brian Bushnell
 * @date Jul 19, 2013
 */
public final class ReadLengthComparator extends ReadComparator {

	/** Private constructor; use the ascending/descending singletons (or the default longest-first comparator). */
	private ReadLengthComparator(int mult_){mult=mult_;}

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
		return mult*x;
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

	@Override
	public ReadComparator getComparator(boolean asc){return asc ? ascending : descending;}

	@Override
	public final int ascendingMult() {return mult;}
	@Override
	public final String name() {return "Length";}

	/** Sort direction multiplier: +1 for ascending (shortest first), -1 for descending (longest first). */
	private final int mult;

	/** Ascending (shortest-first) singleton. */
	public static final ReadLengthComparator ascending=new ReadLengthComparator(1);
	/** Descending (longest-first) singleton. */
	public static final ReadLengthComparator descending=new ReadLengthComparator(-1);
	/** Default singleton; longest-first (descending) by historical default. Alias for selection/identity call sites. */
	public static final ReadLengthComparator comparator=descending;

}
