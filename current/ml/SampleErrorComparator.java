package ml;

import java.util.Comparator;

/**
 * Dormant Sample comparator placeholder for error-based ordering.
 * No Java callers outside this class were found in the reviewed tree.
 * Its comparison method unconditionally asserts under {@code -ea}; with
 * assertions disabled, it delegates to Sample's natural ordering.
 *
 * @author Brian Bushnell
 * @deprecated Dormant and assert-failing; retained only as legacy scaffolding.
 */
@Deprecated
public final class SampleErrorComparator implements Comparator<Sample>{//Seems to be unused...

	/** Private constructor to enforce singleton use. */
	private SampleErrorComparator(){}
	
	/**
	 * Always trips assertions when they are enabled; otherwise returns
	 * {@link Sample#compareTo(Sample)}.
	 */
	@Override
	public int compare(Sample a, Sample b){
//		if(a.errorValue*b.errorValue<0) {
//			return a.errorValue<0 ? -1 : 1;
//		}
		assert(false); // Possible bug: This always fails, making comparator non-functional
		return a.compareTo(b);
	}
	
	/** Dormant singleton comparator instance. */
	public static SampleErrorComparator COMPARATOR=new SampleErrorComparator();
	
}
