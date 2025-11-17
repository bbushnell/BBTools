package bin;

import java.util.Comparator;

/**
 * Implements a custom Comparator for Bin objects with multi-level sorting logic.
 * Provides a consistent comparison method for sorting Bin objects based on score, size, and ID.
 * Sorts primarily by bin score (ascending order), with secondary sorting by bin size,
 * and final tiebreaker using bin ID for deterministic ordering.
 *
 * @author Brian Bushnell
 * @date Feb 4, 2025
 */
class ScoreComparator implements Comparator<Bin> {
	
	/** Private constructor prevents external instantiation of multiple instances */
	private ScoreComparator() {}
	
	@Override
	public int compare(Bin a, Bin b) {
		if(a.score!=b.score) {return a.score<b.score ? -1 : 1;}
		if(a.size()!=b.size()) {return a.size()<b.size() ? -1 : 1;}
		return a.id()-b.id();
	}
	
	/** Singleton instance of ScoreComparator for reuse across the application */
	public static final ScoreComparator comparator=new ScoreComparator();
	
}