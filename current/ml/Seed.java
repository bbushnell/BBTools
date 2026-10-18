package ml;

/**
 * Seed-scan result carrying a network seed and validation pivot.
 * ScannerThread returns these to Trainer, which sorts them by pivot and netSeed
 * before training the winning seeds. Equality and hash codes are seed-only for
 * any hash owner; natural ordering intentionally differs and includes pivot.
 * @author Brian Bushnell
 */
public class Seed implements Comparable<Seed>{

	/**
	 * Constructs a seed-scan result.
	 * @param netSeed_ network initialization seed
	 * @param pivot_ validation pivot used for ordering
	 */
	public Seed(long netSeed_, float pivot_){
		netSeed=netSeed_;
		pivot=pivot_;
	}

//	Seed(long netSeed_, long annealSeed_, float pivot_){
//		netSeed=netSeed_;
//		annealSeed=annealSeed_;
//		pivot=pivot_;
//	}
	
	/**
	 * Orders finite pivots ascending, with NaN pivots after all non-NaN values
	 * including infinities; signed-zero ties and both-NaN ties fall through to
	 * ascending network seed order.
	 */
	@Override
	public int compareTo(Seed s){
		if(pivot!=s.pivot){
			final boolean nan=Float.isNaN(pivot), sNan=Float.isNaN(s.pivot);
			if(nan || sNan){
				if(!nan){return -1;}
				if(!sNan){return 1;}
			}else{
				return pivot>s.pivot ? 1 : -1;
			}
		}
		if(netSeed!=s.netSeed){
			return netSeed>s.netSeed ? 1 : -1;
		}
//		if(annealSeed!=s.annealSeed){
//			return annealSeed>s.annealSeed ? 1 : -1;
//		}
		return 0;
	}

	@Override
	public boolean equals(Object o){
		return o instanceof Seed && equals((Seed)o);
	}

	/** Returns a seed-only hash matching {@link #equals(Object)}. */
	@Override
	public int hashCode(){return Long.hashCode(netSeed);}
	
	/**
	 * Checks seed-only equality; this intentionally differs from
	 * {@link #compareTo(Seed)}, which orders by pivot before seed.
	 */
	public boolean equals(Seed s){
		return s!=null && s.netSeed==netSeed;
	}
	
	/** Returns the network seed and pivot separated by a comma. */
	@Override
	public String toString(){
		return netSeed+", "+pivot;
	}
	
	/** Network initialization seed returned from ScannerThread. */
	final long netSeed;
//	final long annealSeed;
	/** Validation pivot used for seed comparison and sorting. */
	final float pivot;
	
}
