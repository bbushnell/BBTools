package sort;

import java.util.Comparator;

import stream.Read;

/**
 * Abstract base class for implementing comparators to sort Read objects.
 * Provides template method pattern for configurable Read sorting strategies.
 * Used as foundation for specialized Read sorting implementations in BBTools.
 *
 * @author Brian Bushnell
 * @date Nov 9, 2016
 */
public abstract class ReadComparator implements Comparator<Read> {
	
	public abstract void setAscending(boolean asc);
	
}
