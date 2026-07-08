package cardinality;

/**
 * Pluggable stream-key source for referee harnesses (Brian's design,
 * 2026-07-05): each complexity pattern is a subclass, so high-complexity and
 * low-complexity gauntlets are the same program drinking from different
 * bottles.  Construction preserves the historical RNG call order
 * (Random(seed) -> seedBase=nextLong() -> SplittableRandom(nextLong())), so
 * "fresh" and "minrand" produce bit-identical streams to every crown table
 * ever published by WidthWt53/WidthWt53C.
 *
 * Contract: next() returns the next key and maintains distinct(); the caller
 * owns draw budgeting via maxDraws(maxTrue).
 *
 * @author Amber (design: Brian)
 * @date July 2026
 */
public abstract class StreamSource {

	static final long GOLD=0x9E3779B97F4A7C15L;

	protected final long seedBase;
	protected final java.util.SplittableRandom rng;
	protected long distinct=0;

	protected StreamSource(long seed){
		final java.util.Random seedr=new java.util.Random(seed);
		seedBase=seedr.nextLong();
		rng=new java.util.SplittableRandom(seedr.nextLong());
	}

	/** Returns the next stream key; updates distinct(). */
	public abstract long next();

	/** True distinct count so far — the referee's ground truth. */
	public final long distinct(){return distinct;}

	/** Draw budget for a target of maxTrue distincts. */
	public long maxDraws(long maxTrue){return 4L*maxTrue;}

	public static StreamSource make(String name, long maxTrue, long seed){
		switch(name){
			case "fresh": case "HC": return new Fresh(seed);
			case "minrand": return new MinRand(seed, (int)maxTrue);
			case "unifpool": return new UnifPool(seed, (int)maxTrue);
			case "zipfish": return new Zipfish(seed, maxTrue);
			case "onehot": return new Onehot(seed);
			default: throw new IllegalArgumentException("Unknown pattern: "+name
				+" (know: fresh, minrand, unifpool, zipfish, onehot)");
		}
	}

	/** Every key new: the classic high-complexity stream. */
	public static final class Fresh extends StreamSource {
		Fresh(long seed){super(seed);}
		@Override public long next(){distinct++; return seedBase+distinct*GOLD;}
		@Override public long maxDraws(long maxTrue){return maxTrue;}
	}

	/** min(2 uniforms) over a pool: the historical LC crown pattern —
	 * near-constant complexity, small indices hot. */
	public static final class MinRand extends StreamSource {
		private final boolean[] seen;
		private final int pool;
		MinRand(long seed, int pool_){super(seed); pool=pool_; seen=new boolean[pool];}
		@Override public long next(){
			final int pos=Math.min(rng.nextInt(pool), rng.nextInt(pool));
			if(!seen[pos]){seen[pos]=true; distinct++;}
			return seedBase+pos*GOLD;
		}
	}

	/** Uniform over a pool: coupon-collector — complexity FALLS over time. */
	public static final class UnifPool extends StreamSource {
		private final boolean[] seen;
		private final int pool;
		UnifPool(long seed, int pool_){super(seed); pool=pool_; seen=new boolean[pool];}
		@Override public long next(){
			final int pos=rng.nextInt(pool);
			if(!seen[pos]){seen[pos]=true; distinct++;}
			return seedBase+pos*GOLD;
		}
	}

	/** Power-law pool via inverse transform (density skewed toward small
	 * indices) — Zipf-flavored duplication with O(1) memory beyond seen[]. */
	public static final class Zipfish extends StreamSource {
		private final boolean[] seen;
		private final long pool;
		Zipfish(long seed, long pool_){super(seed); pool=pool_; seen=new boolean[(int)pool];}
		@Override public long next(){
			final int pos=(int)(pool*Math.pow(rng.nextDouble(), 3));
			if(!seen[pos]){seen[pos]=true; distinct++;}
			return seedBase+pos*GOLD;
		}
		@Override public long maxDraws(long maxTrue){return 8L*maxTrue;}
	}

	/** 75% re-hits on a 16-element hot set, else fresh: heavy duplication
	 * with near-constant complexity. */
	public static final class Onehot extends StreamSource {
		Onehot(long seed){super(seed);}
		@Override public long next(){
			if(rng.nextDouble()<0.75){return seedBase+rng.nextInt(16)*GOLD-GOLD*20;}
			distinct++;
			return seedBase+distinct*GOLD;
		}
		@Override public long maxDraws(long maxTrue){return 8L*maxTrue;}
	}
}
