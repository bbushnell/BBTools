package ml;

import java.util.Comparator;

@Deprecated
public final class SampleErrorComparator implements Comparator<Sample>{//Seems to be unused...

	private SampleErrorComparator() {}
	
	@Override
	public int compare(Sample a, Sample b) {
//		if(a.errorValue*b.errorValue<0) {
//			return a.errorValue<0 ? -1 : 1;
//		}
		assert(false);
		return a.compareTo(b);
	}
	
	public static SampleErrorComparator COMPARATOR=new SampleErrorComparator();
	
}
