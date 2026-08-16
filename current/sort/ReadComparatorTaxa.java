package sort;

import stream.Read;
import tax.TaxNode;
import tax.TaxTree;


/**
 * Comparator for sorting genomic reads based on their taxonomic classification.
 * Compares reads by hierarchically traversing taxonomic nodes from species to family levels.
 * Handles unclassified reads with name-based fallback sorting.
 *
 * @author Brian Bushnell
 * @date Oct 27, 2014
 */
public final class ReadComparatorTaxa extends ReadComparator {

	/** Private constructor; use the ascending/descending singletons. */
	private ReadComparatorTaxa(int mult_){mult=mult_;}

	/**
	 * Compares two reads based on their taxonomic classification.
	 * Delegates to internal comparison method and applies sort direction.
	 *
	 * @param r1 First read to compare
	 * @param r2 Second read to compare
	 * @return Negative, zero, or positive integer for less than, equal, or greater than
	 */
	@Override
	public int compare(Read r1, Read r2) {
		int x=compareInner(r1, r2);
		return mult*x;
	}

	/** Compares two reads by taxonomy parsed from their headers, climbing to the family (then genus,
	 * then species) level to find a distinguishing rank; falls back to name order for equal/unknown taxa. */
	private static int compareInner(Read r1, Read r2) {
		final TaxNode a0=tree.parseNodeFromHeader(r1.id, true);
		final TaxNode b0=tree.parseNodeFromHeader(r2.id, true);
		final int x=compareNodes(a0, b0, tree);
		//x==0 means both null or the identical taxon; reads then fall back to name order (unchanged behavior).
		if(x!=0){return x;}
		return ReadComparatorName.compareInner(r1, r2);
	}

	/**
	 * Orders two taxonomic nodes exactly as this comparator orders reads: climb both to family, then
	 * (if still equal) genus, then (from the originals) species, breaking the tie at the first
	 * distinguishing rank by extended level and then taxID.  Returns 0 when the two resolve to the
	 * same taxon (both null, or the identical node) so the caller can apply its own tiebreak; an
	 * unresolvable (null) node sorts last.  This is the reusable, Read-independent core of the taxa
	 * sort, so spectra records (clade.Clade) and DDL sketch records (ddl.DDLRecord) can sort by the
	 * same ordering bbsort.sh taxa uses -- clustering taxonomically-similar records adjacently, which
	 * is what lets a bgzipped per-taxID database deflate tightly.
	 * @param a0 First node, or null if its taxID was unresolvable.
	 * @param b0 Second node, or null.
	 * @param tree Loaded taxonomy tree (ancestor lookups use it).
	 * @return Negative, zero, or positive for a0 sorting before, equal to, or after b0.
	 */
	public static int compareNodes(final TaxNode a0, final TaxNode b0, final TaxTree tree){
		if(a0==null || b0==null){
			if(a0==null && b0==null){return 0;}
			return (a0==null ? 1 : -1);
		}
		if(a0==b0){return 0;}
		TaxNode a=a0, b=b0;

		while(a.id!=a.pid && a.levelExtended<TaxTree.FAMILY_E){
			TaxNode x=tree.getNode(a.pid);
			if(x.levelExtended>TaxTree.FAMILY_E){break;}
			a=x;
		}
		while(b.id!=b.pid && b.levelExtended<TaxTree.FAMILY_E){
			TaxNode x=tree.getNode(b.pid);
			if(x.levelExtended>TaxTree.FAMILY_E){break;}
			b=x;
		}

		if(a==b){
			while(a.id!=a.pid && a.levelExtended<TaxTree.GENUS_E){
				TaxNode x=tree.getNode(a.pid);
				if(x.levelExtended>TaxTree.GENUS_E){break;}
				a=x;
			}
			while(b.id!=b.pid && b.levelExtended<TaxTree.GENUS_E){
				TaxNode x=tree.getNode(b.pid);
				if(x.levelExtended>TaxTree.GENUS_E){break;}
				b=x;
			}

			if(a==b){
				a=a0;
				b=b0;
				while(a.id!=a.pid && a.levelExtended<TaxTree.SPECIES_E){
					TaxNode x=tree.getNode(a.pid);
					if(x.levelExtended>TaxTree.SPECIES_E){break;}
					a=x;
				}
				while(b.id!=b.pid && b.levelExtended<TaxTree.SPECIES_E){
					TaxNode x=tree.getNode(b.pid);
					if(x.levelExtended>TaxTree.SPECIES_E){break;}
					b=x;
				}

				if(a==b){
					return compareSimple(a0, b0);
				}
			}
		}

		return compareSimple(a, b);
	}

	/** Orders two taxonomic nodes by extended level, then by taxonomic id. */
	private static int compareSimple(final TaxNode a, final TaxNode b){
		if(a.levelExtended!=b.levelExtended){return a.levelExtended-b.levelExtended;}
		return a.id-b.id;
	}

	@Override
	public ReadComparator getComparator(boolean asc){return asc ? ascending : descending;}

	@Override
	public final int ascendingMult() {return mult;}
	@Override
	public final String name() {return "Taxa";}

	/** Sort direction multiplier: +1 for ascending, -1 for descending (immutable). */
	private final int mult;

	/** Taxonomy tree used to resolve and climb taxonomic nodes; must be loaded before sorting. */
	public static TaxTree tree;

	/** Ascending-order singleton. */
	public static final ReadComparatorTaxa ascending=new ReadComparatorTaxa(1);
	/** Descending-order singleton. */
	public static final ReadComparatorTaxa descending=new ReadComparatorTaxa(-1);
	/** Default (ascending) singleton; alias retained for selection/identity call sites. */
	public static final ReadComparatorTaxa comparator=ascending;

}
