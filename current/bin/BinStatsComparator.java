package bin;

import java.util.Comparator;

/**
 * Comparator for BinStats objects that provides multi-level sorting criteria.
 * Sorts bins by contamination level (ascending), then by size (descending),
 * then by ID (ascending) as a tie-breaker.
 * @author Brian Bushnell
 */
class BinStatsComparator implements Comparator<BinStats> {
	
	@Override
	public int compare(BinStats a, BinStats b) {
		if(a.contam!=b.contam) {
			return a.contam<b.contam ? -1 : 1;
		}
		if(a.size!=b.size) {
			return a.size<b.size ? 1 : -1;
		}
		return a.id-b.id;
	}
	
}