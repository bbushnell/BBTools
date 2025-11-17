package structures;

/**
 * Interface for objects that have a settable location and are comparable.
 * Provides a contract for data structures that maintain positional information
 * and can be ordered based on their location values.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public interface SetLoc<T> extends Comparable<T> {
	
	/** Sets the location of this object.
	 * @param newLoc The new location value to assign */
	public void setLoc(int newLoc);
	/** Returns the current location of this object */
	public int loc();
	
}
