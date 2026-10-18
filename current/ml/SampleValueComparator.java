package ml;

import java.util.Comparator;

/**
 * Sorts Samples by first prediction value in ascending order, with NaNs after
 * all non-NaN values including infinities.
 * <p>
 * Equal values, signed-zero ties and both-NaN ties fall through to
 * {@link Sample#id}.
 * @author Brian Bushnell
 * @date 2014
 */
public final class SampleValueComparator implements Comparator<Sample>{

	/** Private constructor to enforce singleton use. */
	private SampleValueComparator(){}
	
	/**
	 * Compares first prediction values, then falls back to sample id for
	 * equal values, signed-zero ties and both-NaN ties.
	 */
	@Override
	public int compare(Sample a, Sample b){
		float ar=a.result[0], br=b.result[0];
		if(ar>br){return 1;}
		if(ar<br){return -1;}
		final boolean aNan=Float.isNaN(ar), bNan=Float.isNaN(br);
		if(aNan!=bNan){return aNan ? 1 : -1;}
		return a.id-b.id;
	}
	
	/** Singleton comparator instance. */
	public static SampleValueComparator COMPARATOR=new SampleValueComparator();
	
}
