package ml;

public class Seed implements Comparable<Seed>{

	Seed(long netSeed_, float pivot_){
		netSeed=netSeed_;
		pivot=pivot_;
	}

//	Seed(long netSeed_, long annealSeed_, float pivot_){
//		netSeed=netSeed_;
//		annealSeed=annealSeed_;
//		pivot=pivot_;
//	}
	
	@Override
	public int compareTo(Seed s) {
		if(pivot!=s.pivot) {
			return pivot>s.pivot ? 1 : -1;
		}
		if(netSeed!=s.netSeed) {
			return netSeed>s.netSeed ? 1 : -1;
		}
//		if(annealSeed!=s.annealSeed) {
//			return annealSeed>s.annealSeed ? 1 : -1;
//		}
		return 0;
	}

	@Override
	public boolean equals(Object o) {
		return equals((Seed)o);
	}
	
	public boolean equals(Seed s) {
		return s.netSeed==netSeed;// && s.annealSeed==annealSeed;
	}
	
	@Override
	public String toString() {
		return netSeed+/*", "+annealSeed+*/", "+pivot;
	}
	
	final long netSeed;
//	final long annealSeed;
	final float pivot;
	
}
