package align2;

/**
 * Represents a coordinate point with column, row, and site position values.
 * Used in alignment algorithms for tracking positions and values within
 * dynamic programming matrices or similar data structures.
 * Implements Comparable for sorting based on site position and column.
 *
 * @author Brian Bushnell
 */
public class Quad implements Comparable<Quad>{
	
	/**
	 * Constructs a Quad with specified column, row, and site values.
	 * @param col_ The column position
	 * @param row_ The row position
	 * @param val_ The site value
	 */
	public Quad(int col_, int row_, int val_){
		column=col_;
		row=row_;
		site=val_;
	}
	
	@Override
	public boolean equals(Object other){
		return site==((Quad)other).site;
	}
	
	@Override
	public int hashCode(){return site;}
	
	@Override
	public int compareTo(Quad other) {
		int x=site-other.site;
		return(x==0 ? column-other.column : x);
	}
	
	@Override
	public String toString(){
		return("("+column+","+row+","+site+")");
	}
	
	/** The column position (immutable) */
	public final int column;
	/** The row position */
	public int row;
	/** The site value used for equality and primary sorting */
	public int site;
	/**
	 * Array for storing additional integer values associated with this coordinate
	 */
	public int list[];
	
}
