package sort;

import java.util.Comparator;

import assemble.Contig;
import shared.Tools;

/**
 * Comparator for sorting contigs by length with secondary tie-breaking criteria.
 * Sorts longest contig first by default (descending order), with options for ascending order.
 * Uses multiple comparison levels: length, coverage, sequence content, and ID.
 *
 * @author Brian Bushnell
 * @date Jul 12, 2018
 */
public final class ContigLengthComparator implements Comparator<Contig> {

	/** Private constructor; use the ascending/descending singletons. */
	private ContigLengthComparator(int mult_){mult=mult_;}

	/**
	 * Compares two contigs using multi-level sorting criteria.
	 * Primary: length difference, Secondary: coverage, Tertiary: sequence content, Quaternary: ID.
	 * Result is multiplied by the direction multiplier to control sort direction.
	 *
	 * @param a First contig to compare
	 * @param b Second contig to compare
	 * @return Negative/zero/positive as first contig is less/equal/greater than second
	 */
	@Override
	public int compare(Contig a, Contig b) {
		int x=compareInner(a, b);
		if(x==0){x=a.coverage>b.coverage ? 1 : a.coverage<b.coverage ? -1 : 0;}
		if(x==0){x=compareVectors(a.bases, b.bases);}
		if(x==0){x=a.id>b.id ? 1 : a.id<b.id ? -1 : 0;}
		return mult*x;
	}

	/**
	 * Primary length-based comparison with null safety; treats null contigs as greater than non-null.
	 * @param a First contig (may be null)
	 * @param b Second contig (may be null)
	 * @return Length difference (a.length - b.length) or null ordering
	 */
	private static int compareInner(Contig a, Contig b) {
		if(a==b){return 0;}
		if(a==null){return 1;}
		if(b==null){return -1;}
		int x=a.length()-b.length();
		return x;
	}

	/**
	 * Lexicographic comparison of sequence byte arrays; treats null as greater than non-null.
	 * @param a First byte array (may be null)
	 * @param b Second byte array (may be null)
	 * @return Negative, zero, or positive for ordering
	 */
	private static int compareVectors(final byte[] a, final byte[] b){
		if(a==null || b==null){
			if(a==null && b!=null){return 1;}
			if(a!=null && b==null){return -1;}
			return 0;
		}
		final int lim=Tools.min(a.length, b.length);
		for(int i=0; i<lim; i++){
			if(a[i]<b[i]){return -1;}
			if(a[i]>b[i]){return 1;}
		}
		return 0;
	}

	/** @return the directional singleton (ascending/shortest-first if true, else descending/longest-first). */
	public static ContigLengthComparator getComparator(boolean asc){return asc ? ascending : descending;}

	/** Sort direction multiplier: +1 for ascending (shortest first), -1 for descending (longest first). */
	private final int mult;

	/** Ascending (shortest-first) singleton. */
	public static final ContigLengthComparator ascending=new ContigLengthComparator(1);
	/** Descending (longest-first) singleton. */
	public static final ContigLengthComparator descending=new ContigLengthComparator(-1);
	/** Default singleton; longest-first (descending) by historical default. Alias for call sites. */
	public static final ContigLengthComparator comparator=descending;

}
