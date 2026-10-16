package structures;

/**
 * Comparable element that can store its current {@link HeapLoc} array index.
 * {@code HeapLoc} uses this callback to repair an entry in place after the
 * sort value changes.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public interface SetLoc<T> extends Comparable<T>{

	/** Records the one-based heap index, or a negative value when absent.
	 * @param newLoc New heap location */
	public void setLoc(int newLoc);
	/** Returns the last heap location assigned by {@link #setLoc(int)}. */
	public int loc();

}
