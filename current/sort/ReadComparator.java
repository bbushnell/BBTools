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
	
	/** @return the directional singleton for the requested order (ascending if true, else descending). */
	public abstract ReadComparator getComparator(boolean ascending);

	/** @return true if this comparator currently sorts in ascending order. */
	public final boolean ascending() {return ascendingMult()>0;}

	/** @return the direction multiplier: +1 for ascending, -1 for descending. */
	public abstract int ascendingMult();

	/** @return a short human-readable name for this sort order. */
	public abstract String name();
	
}
