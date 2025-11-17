package ml;

import java.util.Comparator;

/**
 * Deprecated comparator for Sample objects based on error values.
 * This class appears to be unused and contains non-functional comparison logic.
 * The compare method always asserts false and delegates to Sample's natural ordering.
 *
 * @author Brian Bushnell
 * @deprecated This class is marked as deprecated and seems to be unused
 */
@Deprecated
public final class SampleErrorComparator implements Comparator<Sample>{//Seems to be unused...

	/** Private constructor prevents external instantiation */
	private SampleErrorComparator() {}
	
	@Override
	public int compare(Sample a, Sample b) {
//		if(a.errorValue*b.errorValue<0) {
//			return a.errorValue<0 ? -1 : 1;
//		}
		assert(false); // Possible bug: This always fails, making comparator non-functional
		return a.compareTo(b);
	}
	
	/** Singleton instance of the SampleErrorComparator */
	public static SampleErrorComparator COMPARATOR=new SampleErrorComparator();
	
}
