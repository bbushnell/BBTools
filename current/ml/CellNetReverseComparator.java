package ml;

import java.util.Comparator;

/**
 * Reverses {@link CellNet#compareTo(CellNet)} for best-first sorting.
 * Uses that method's selected metric, tie breakers and invalid-metric behavior;
 * this wrapper does not impose a separate ordering policy.
 *
 * @author Brian Bushnell
 * @date October 25, 2013
 */
public class CellNetReverseComparator implements Comparator<CellNet>{

	/**
	 * Negates the first network's natural comparison result.
	 * @param a Non-null first network
	 * @param b Second network, passed unchanged to the natural comparator
	 * @return Negative when a ranks better, positive when b ranks better, or zero for a tie
	 */
	@Override
	public int compare(CellNet a, CellNet b){
		return -a.compareTo(b);
	}

}
