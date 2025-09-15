package ml;

import java.util.Comparator;

public final class SampleValueComparator implements Comparator<Sample>{

	private SampleValueComparator() {}
	
	@Override
	public int compare(Sample a, Sample b) {
		float ar=a.result[0], br=b.result[0];
		return ar>br ? 1 : ar<br ? -1 : a.id-b.id;
	}
	
	public static SampleValueComparator COMPARATOR=new SampleValueComparator();
	
}
