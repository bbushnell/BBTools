package sort;

import java.util.Comparator;

import stream.Read;

/**
 * Base API for named read orders with ascending and descending variants.
 * Implementations define their keys, null handling and shared-state requirements.
 *
 * @author Brian Bushnell
 * @date Nov 9, 2016
 */
public abstract class ReadComparator implements Comparator<Read>{
	
	/** Selects an order; implementations may return a singleton or share immutable configuration.
	 * @param ascending True for ascending, false for descending
	 * @return Comparator configured for the requested direction */
	public abstract ReadComparator getComparator(boolean ascending);

	/** @return true if this comparator currently sorts in ascending order. */
	public final boolean ascending(){return ascendingMult()>0;}

	/** @return the direction multiplier: +1 for ascending, -1 for descending. */
	public abstract int ascendingMult();

	/** @return a short human-readable name for this sort order. */
	public abstract String name();
	
}
