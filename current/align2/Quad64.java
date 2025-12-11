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
	 * Compares this Quad64 with another object for equality based on site value.
	 * An assertion is present in the implementation; primarily used for debugging.
	 * @param other Object to compare
	 * @return true if site values match; may assert in debug scenarios
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
