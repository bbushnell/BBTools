package align2;

/**
 * Represents a coordinate position in alignment matrices with column, row, and site values.
 * Used for storing alignment positions in dynamic programming matrices with natural ordering based on site values.
 * @author Brian Bushnell
 * @date December 21, 2010
 */
public class Quad64 implements Comparable<Quad64>{
	
	public Quad64(int col_, int row_, int val_){
		column=col_;
		row=row_;
		site=val_;
	}
	
	/**
	 * Intentionally forbidden — asserts false. Quad64 lives only in Quad64Heap (compareTo) and
	 * raw arrays, never as a hash/set key, so equals must never be invoked; with -ea (always on)
	 * any call crashes loudly as a tripwire. The site comparison below is unreachable.
	 * @param other Unused; equals must never be called on Quad64
	 * @return never returns normally under assertions
	 */
	@Override
	public boolean equals(Object other){
		assert(false);
		return site==((Quad64)other).site;
	}
	
	/** Returns the site value as the hash code.
	 * @return Hash code derived from site */
	@Override
	public int hashCode(){return (int)site;}
	
	/**
	 * Compares Quad64 objects by site value, then by column for tie-breaking.
	 * @param other Quad64 to compare against
	 * @return Negative if this < other, positive if this > other, 0 if equal
	 */
	@Override
	public int compareTo(Quad64 other) {
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
	
	/** Column position in alignment matrix. */
	public final int column;
	/** Row position in alignment matrix. */
	public int row;
	/** Site value used for position scoring and comparison. */
	public long site;
	/** Array for storing additional position-related data. */
	public int list[];
	
}
