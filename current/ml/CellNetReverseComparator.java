package ml;

import java.util.Comparator;

/**
 * Reverses the natural ordering of CellNet neural networks for sorting purposes.
 * When used with sorting algorithms, this comparator arranges networks so that
 * the best performing networks appear first (in descending order of quality).
 * Since CellNet.compareTo() returns positive values for better networks,
 * this reversal creates standard best-first ordering.
 *
 * @author Brian Bushnell
 * @date October 25, 2013
 */
public class CellNetReverseComparator implements Comparator<CellNet>{

	@Override
	public int compare(CellNet a, CellNet b) {
		return -a.compareTo(b);
	}

}
