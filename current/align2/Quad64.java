package align2;

/** Mutable hit-list cursor ordered by a long site value, then seed column.
 * BBIndex5 assigns packed int coordinates widened with &0xFFFFFFFFL; thus bit31
 * does not make a coordinate negative. compareTo itself uses signed long order.
 * Row indexes the caller-owned hit list, not a dynamic-programming matrix.
 * Poll before changing site, then reinsert. Mutable, unsynchronized, and not
 * intended as a hash/set key.
 * @author Brian Bushnell
 * @date December 21, 2010
 */
public class Quad64 implements Comparable<Quad64>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** The int site is sign-extended here; BBIndex5 starts at zero, then assigns
	 * its unsigned coordinate separately. Seed columns are nonnegative indices. */
	public Quad64(final int col_, final int row_, final int val_){
		column=col_;
		row=row_;
		site=val_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Unsupported for heap cursors: the legacy assertion intentionally rejects calls.
	 * Quad64Heap uses compareTo and reference identity, not this method.
	 * With assertions disabled the historical site-only comparison still executes;
	 * it ignores column and is not a supported value-equality contract.
	 * @param other Unused; equals must never be called on Quad64
	 * @return never returns normally under assertions
	 */
	@Override
	public boolean equals(final Object other){
		assert(false) : "Quad64.equals is unsupported: Quad64Heap uses compareTo and reference identity";
		return site==((Quad64)other).site;
	}

	/** Returns the low 32 site bits; this does not make the mutable cursor a supported key. */
	@Override
	public int hashCode(){return (int)site;}

	/**
	 * Compares Quad64 objects by site value, then by column for tie-breaking.
	 * @param other Quad64 to compare against
	 * @return Negative if this < other, positive if this > other, 0 if equal
	 */
	@Override
	public int compareTo(final Quad64 other){
		//Comparison-based, not (site-other.site): site is a long, so subtract-then-cast-to-int could
		//mis-order (the commented-out version below is exactly that trap). Column tie-break is a
		//non-negative bounded index, so its subtraction is overflow-safe.
		return site>other.site ? 1 : site<other.site ? -1 : column-other.column;
//		int x=site-other.site;
//		return(x>0 ? 1 : x<0 ? -1 : column-other.column);
	}

	/** Returns string representation showing column, row, and site values in the format "(column,row,site)".
	 * @return Formatted string representation */
	@Override
	public String toString(){
		return("("+column+","+row+","+site+")");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Immutable seed/hit-list column used for tie-breaking. */
	public final int column;
	/** Current index into list. */
	public int row;
	/** Adjusted packed coordinate; BBIndex5 assigns an unsigned int widened to long. */
	public long site;
	/** Caller-owned hit positions; the cursor does not copy their contents. */
	public int[] list;

}
