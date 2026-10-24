package sort;

import shared.Tools;
import stream.Read;



/**
 * Comparator for Read objects that performs hierarchical topological ordering of
 * primary bases and length, optional mates, qualities, and identifiers. Supports
 * ascending or descending order and exposes a singleton comparator instance.
 *
 * @author Brian Bushnell
 * @date Oct 27, 2014
 */
public class ReadComparatorTopological extends ReadComparator{

	/** Private constructor; use the ascending/descending singletons. */
	private ReadComparatorTopological(int mult_){mult=mult_;}

	/**
	 * Overrides the base comparator to apply the current ascending/descending
	 * direction while comparing two reads.
	 *
	 * @param r1 First read to compare
	 * @param r2 Second read to compare
	 * @return Negative if {@code r1 < r2}, zero if equal, positive if {@code r1 > r2}
	 */
	@Override
	public int compare(Read r1, Read r2){
		return mult*compare(r1, r2, true);
	}

	/**
	 * Performs the full hierarchical comparison: primary bases and length, optional
	 * mate bases and length, quality scores (inverted), optional mate qualities
	 * (inverted), numeric ID, then string ID.
	 *
	 * @param r1 First read to compare
	 * @param r2 Second read to compare
	 * @param compareMates Whether to include mate presence, bases, length, and qualities
	 * @return Negative if {@code r1 < r2}, zero if equal, positive if {@code r1 > r2}
	 */
	public int compare(Read r1, Read r2, boolean compareMates){

		int x=compareVectors(r1.bases, r2.bases);
		if(x!=0){return x;}

		if(r1.bases!=null && r2.bases!=null && r1.length()!=r2.length()){return r1.length()-r2.length();}

		//FIXED sort/ReadComparatorTopological#001: primary length precedes
		//optional mate keys; mixed mate presence gets a total order.
		if(compareMates){
			if(r1.mate==null || r2.mate==null){
				if(r1.mate==null && r2.mate!=null){return 1;}
				if(r1.mate!=null && r2.mate==null){return -1;}
			}else{
				x=compareVectors(r1.mate.bases, r2.mate.bases);
				if(x!=0){return x;}
				if(r1.mate.bases!=null && r2.mate.bases!=null
						&& r1.mateLength()!=r2.mateLength()){return r1.mateLength()-r2.mateLength();}
			}
		}

		x=compareVectors(r1.quality, r2.quality);
		if(x!=0){return 0-x;}

		if(compareMates && r1.mate!=null && r2.mate!=null){
			x=compareVectors(r1.mate.quality, r2.mate.quality);
		}
		if(x!=0){return 0-x;}

		if(r1.numericID!=r2.numericID){return r1.numericID>r2.numericID ? 1 : -1;}

		return r1.id.compareTo(r2.id);
	}

	/**
	 * Lexicographically compares two byte arrays representing read bases or qualities.
	 * Handles null arrays by treating null as greater than non-null.
	 *
	 * @param a First byte array to compare
	 * @param b Second byte array to compare
	 * @return Negative if {@code a < b}, zero if equal up to minimum length,
	 * positive if {@code a > b}
	 */
	public int compareVectors(final byte[] a, final byte[] b){
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

	/**
	 * Lexicographically compares two byte arrays, treating 'N' as greater than other
	 * bases before falling back to standard comparison.
	 *
	 * @param a First byte array to compare
	 * @param b Second byte array to compare
	 * @return Negative if {@code a < b}, zero if equal up to minimum length,
	 * positive if {@code a > b}
	 */
	public int compareVectorsN(final byte[] a, final byte[] b){
		if(a==null || b==null){
			if(a==null && b!=null){return 1;}
			if(a!=null && b==null){return -1;}
			return 0;
		}
		final int lim=Tools.min(a.length, b.length);
		for(int i=0; i<lim; i++){
			if(a[i]=='N' && b[i]!='N'){return 1;}
			if(a[i]!='N' && b[i]=='N'){return -1;}
			if(a[i]<b[i]){return -1;}
			if(a[i]>b[i]){return 1;}
		}
		return 0;
	}

	@Override
	public ReadComparator getComparator(boolean asc){return asc ? ascending : descending;}

	@Override
	public final int ascendingMult(){return mult;}
	@Override
	public final String name(){return "Topology";}

	/** Sort direction multiplier: +1 for ascending, -1 for descending (immutable). */
	private final int mult;

	/** Ascending-order singleton. */
	public static final ReadComparatorTopological ascending=new ReadComparatorTopological(1);
	/** Descending-order singleton. */
	public static final ReadComparatorTopological descending=new ReadComparatorTopological(-1);
	/** Default (ascending) singleton; alias retained for selection/identity call sites. */
	public static final ReadComparatorTopological comparator=ascending;
}
