package align2;

/**
 * Represents a coordinate point with column, row, and site position values.
 * Used in alignment algorithms for tracking positions within dynamic programming matrices, with Comparable ordering by site then column.
 * @author Brian Bushnell
 */
public class Quad implements Comparable<Quad>{
	
	public Quad(int col_, int row_, int val_){
		column=col_;
		row=row_;
		site=val_;
	}
	
	/**
	 * Compares this Quad with another object for equality based on site value.
	 * @param other The object to compare with
	 * @return true if both Quads have the same site value
	 */
	@Override
	public boolean equals(Object other){
		//align2/Quad#001 FIXED: equals was site-only, contradicting compareTo (site,column) — a latent
		//landmine. equals/hashCode are never invoked today (Quad lives only in QuadHeap and raw Quad[]
		//arrays, never a hash/tree key), so aligning equals to compareTo's identity (site AND column)
		//changes no active path; hashCode=site stays contract-valid since equal Quads always share site.
		return site==((Quad)other).site && column==((Quad)other).column;
	}
	
	/** Returns the site value as the hash code. */
	@Override
	public int hashCode(){return site;}
	
	/**
	 * Compares this Quad with another for ordering: primary by site value, secondary by column value.
	 * @param other The Quad to compare with
	 * @return Negative if this < other, positive if this > other, zero if equal
	 */
	@Override
	public int compareTo(Quad other) {
		int x=site-other.site;
		return(x==0 ? column-other.column : x);
	}
	
	/** Returns a string representation of this Quad in the format "(column,row,site)".
	 * @return String representation of this Quad */
	@Override
	public String toString(){
		return("("+column+","+row+","+site+")");
	}
	
	/** The column position (immutable). */
	public final int column;
	/** The row position. */
	public int row;
	/** The site value used for equality and primary sorting. */
	public int site;
	/**
	 * Array for storing additional integer values associated with this coordinate.
	 */
	public int list[];
	
}
