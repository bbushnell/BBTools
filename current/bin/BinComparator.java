package bin;

import java.util.Comparator;

/**
 * Comparator for Bin objects that prioritizes bins by contamination level, size, and ID.
 * Establishes sorting order for bin processing and display in downstream analysis.
 * Lower contamination bins are prioritized first, followed by larger bins, then lower IDs.
 *
 * @author Brian Bushnell
 * @date February 2025
 */
class BinComparator implements Comparator<Bin> {
	
	@Override
	public int compare(Bin a, Bin b) {
		if(a.contam!=b.contam) {
			return a.contam<b.contam ? -1 : 1;
		}
		if(a.size()!=b.size()) {
			return a.size()<b.size() ? 1 : -1;
		}
		return a.id()-b.id();
	}
	
}