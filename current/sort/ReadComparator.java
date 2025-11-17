package sort;

import java.util.Comparator;

import stream.Read;

/**
 * @author Brian Bushnell
 * @date Nov 9, 2016
 *
 */
public abstract class ReadComparator implements Comparator<Read> {
	
	/** Sets the sort order for this comparator.
	 * @param asc true for ascending order, false for descending order */
	public abstract void setAscending(boolean asc);
	
}
