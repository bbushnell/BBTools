package bin;

import java.util.Comparator;

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