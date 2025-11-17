package align2;

/**
 * Represents a coordinate position in alignment matrices with column, row, and site values.
 * Used for storing alignment positions in dynamic programming matrices with natural ordering
 * based on site values. Provides basic coordinate operations and comparison functionality.
 *
 * @author Brian Bushnell
 * @date December 21, 2010
 */
public class Quad64 implements Comparable<Quad64>{
	
	/**
	 * Constructs a Quad64 with specified column, row, and site values.
	 * @param col_ Column position in alignment matrix
	 * @param row_ Row position in alignment matrix
	 * @param val_ Site value for position scoring and comparison
	 */
	public Quad64(int col_, int row_, int val_){
		column=col_;
		row=row_;
		site=val_;
	}
	
	@Override
	public boolean equals(Object other){
		assert(false);
		return site==((Quad64)other).site;
	}
	
	@Override
	public int hashCode(){return (int)site;}
	
	@Override
	public int compareTo(Quad64 other) {
		return site>other.site ? 1 : site<other.site ? -1 : column-other.column;
//		int x=site-other.site;
//		return(x>0 ? 1 : x<0 ? -1 : column-other.column);
	}
	
	@Override
	public String toString(){
		return("("+column+","+row+","+site+")");
	}
	
	/** Column position in alignment matrix */
	public final int column;
	/** Row position in alignment matrix */
	public int row;
	/** Site value used for position scoring and comparison */
	public long site;
	/** Array for storing additional position-related data */
	public int list[];
	
}
