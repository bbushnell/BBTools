package sort;

import java.util.Comparator;

import stream.Read;

/**
 * @author Brian Bushnell
 * @date Nov 9, 2016
 *
 */
public abstract class ReadComparator implements Comparator<Read> {
	
	public abstract void setAscending(boolean asc);
	
}
