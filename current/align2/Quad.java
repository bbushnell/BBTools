package align2;

/** Mutable cursor for merging index hit lists, ordered by site then column.
 * BBIndex's object walk assigns one cursor per seed column; row selects an entry
 * in list, and site is that entry adjusted by the seed's query offset. These are
 * hit-list coordinates, not dynamic-programming matrix cells.
 * Equality uses site and column; row and list do not participate. Do not mutate
 * site while this object is a hash key or while its heap position depends on it.
 * The index walker polls, updates, and reinserts the cursor. Not thread-safe.
 * @author Brian Bushnell
 */
public class Quad implements Comparable<Quad>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a cursor; the caller assigns its hit list separately. */
	public Quad(final int col_, final int row_, final int val_){
		column=col_;
		row=row_;
		site=val_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns true for a Quad with the same site and column, regardless of row/list.
	 * Null and objects of unrelated types are unequal. */
	@Override
	public boolean equals(final Object other){
		if(!(other instanceof Quad)){return false;}
		final Quad q=(Quad)other;
		return site==q.site && column==q.column;
	}

	/** Returns site; equal cursors share this value, though different columns collide. */
	@Override
	public int hashCode(){return site;}

	/**
	 * Compares this Quad with another for ordering: primary by site value, secondary by column value.
	 * @param other The Quad to compare with
	 * @return Negative if this < other, positive if this > other, zero if equal
	 */
	@Override
	public int compareTo(final Quad other){
		//TODO: Probable bug - subtraction can overflow for unrestricted int values
		// (MAX_VALUE compared with -1 returns negative). Index-walk reachability
		// needs separate validation before changing its ordering or primitive analogues.
		final int x=site-other.site;
		return(x==0 ? column-other.column : x);
	}

	/** Returns a string representation of this Quad in the format "(column,row,site)".
	 * @return String representation of this Quad */
	@Override
	public String toString(){
		return("("+column+","+row+","+site+")");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Immutable seed/hit-list column, used to break ties at the same site. */
	public final int column;
	/** Current entry index within list. */
	public int row;
	/** The site value used for equality and primary sorting. */
	public int site;
	/** Hit positions supplied by the index; ownership remains with the caller. */
	public int[] list;

}
