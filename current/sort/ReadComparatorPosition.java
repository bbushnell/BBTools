package sort;

import stream.Read;
import stream.SamLine;
import var2.ScafMap;

/**
 * @author Brian Bushnell
 * @date November 20, 2016
 *
 */

public final class ReadComparatorPosition extends ReadComparator {
	
	/** Private constructor to enforce singleton pattern */
	private ReadComparatorPosition(){}
	
	@Override
	public int compare(Read r1, Read r2) {
		int x=compareInner(r1, r2);
		return ascending*x;
	}
	
	/**
	 * Core comparison logic for reads using SAM alignment data and read IDs.
	 * Delegates to SAM-based comparison, then falls back to read ID comparison.
	 *
	 * @param r1 First read to compare
	 * @param r2 Second read to compare
	 * @return Negative if r1 < r2, positive if r1 > r2, zero if equal
	 */
	public static int compareInner(Read r1, Read r2) {
		int x=compareInner(r1.samline, r2.samline);
		if(x!=0){return x;}
		if(r1.id==null && r2.id==null){return 0;}
		if(r1.id==null){return -1;}
		if(r2.id==null){return 1;}
		return r1.id.compareTo(r2.id);
	}
	
	/**
	 * Compares SAM alignment records by genomic position and metadata.
	 * Hierarchical comparison: scaffold number, position, strand, mate position, pair number.
	 * Sets scaffold numbers from ScafMap if not already assigned.
	 *
	 * @param a First SAM alignment record
	 * @param b Second SAM alignment record
	 * @return Negative if a < b, positive if a > b, zero if equal
	 */
	public static int compareInner(SamLine a, SamLine b) {
		if(a.scafnum<0){a.setScafnum(scafMap);}
		if(b.scafnum<0){b.setScafnum(scafMap);}
		if(a.scafnum!=b.scafnum){return a.scafnum-b.scafnum;}
		if(a.pos!=b.pos){return a.pos-b.pos;}
		if(a.strand()!=b.strand()){return a.strand()-b.strand();}
		if(a.pnext!=b.pnext){return a.pnext-b.pnext;}
		if(a.pairnum()!=b.pairnum()){return a.pairnum()-b.pairnum();}
		return 0;
	}
	
	/** Multiplier for sort direction: 1 for ascending, -1 for descending */
	private int ascending=1;
	
	@Override
	public void setAscending(boolean asc){
		ascending=(asc ? 1 : -1);
	}

	/** Singleton instance for position-based read comparison */
	public static final ReadComparatorPosition comparator=new ReadComparatorPosition();
	/** Scaffold mapping used to resolve scaffold numbers from names */
	public static ScafMap scafMap=null;
	
}
