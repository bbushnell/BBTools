package ml;

/**
 * Represents a deterministic seed configuration for neural network training.
 * Provides a comparable representation of network initialization states and performance metrics.
 * Used to ensure reproducible random states and facilitate seed comparison during model training.
 * @author Brian Bushnell
 */
public class Seed implements Comparable<Seed>{

	/**
	 * Constructs a Seed with network seed and performance pivot.
	 * @param netSeed_ The network initialization seed value
	 * @param pivot_ The performance pivot value for comparison
	 */
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
	
	/**
	 * Checks equality with another Seed based on network seed only.
	 * @param s The Seed to compare against
	 * @return true if network seeds are equal, false otherwise
	 */
	public boolean equals(Seed s) {
		return s.netSeed==netSeed;// && s.annealSeed==annealSeed;
	}
	
	@Override
	public String toString() {
		return netSeed+/*", "+annealSeed+*/", "+pivot;
	}
	
	/** The network initialization seed value used for reproducible random states */
	final long netSeed;
//	final long annealSeed;
	/** The performance pivot value used for seed comparison and sorting */
	final float pivot;
	
}
