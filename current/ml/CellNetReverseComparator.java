package ml;

import java.util.Comparator;

public class CellNetReverseComparator implements Comparator<CellNet>{

	@Override
	public int compare(CellNet a, CellNet b) {
		return -a.compareTo(b);
	}

}
