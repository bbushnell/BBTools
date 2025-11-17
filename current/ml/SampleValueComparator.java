package ml;

import java.util.Comparator;

/**
 * Comparator for sorting samples by their neural network output values.
 * Implements a three-way comparison mechanism for Sample objects, sorting based
 * on their primary result value with an ID-based tie-breaker for deterministic ordering.
 *
 * @author Brian Bushnell
 * @date 2014
 */
public final class SampleValueComparator implements Comparator<Sample>{

	/** Private constructor to enforce singleton pattern */
	private SampleValueComparator() {}
	
	@Override
	public int compare(Sample a, Sample b) {
		float ar=a.result[0], br=b.result[0];
		return ar>br ? 1 : ar<br ? -1 : a.id-b.id;
	}
	
	/** Singleton instance for consistent access across the application */
	public static SampleValueComparator COMPARATOR=new SampleValueComparator();
	
}
