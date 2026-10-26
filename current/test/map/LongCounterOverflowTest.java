package test.map;

import java.math.BigInteger;
import java.util.Random;
import map.IntLongHashMap2;
import map.IntLongHashMap3;
import map.LongLongHashMap2;

/** Checks positive saturation independently of the maps' array positions.
 * @author Shinobu
 */
public final class LongCounterOverflowTest{

	public static void main(String[] args){
		check(args.length==1, "Expected expect=baseline|fixed");
		final boolean fixed=args[0].equals("expect=fixed");
		check(fixed || args[0].equals("expect=baseline"), "Unknown expectation");
		final long expected=fixed ? Long.MAX_VALUE : Long.MIN_VALUE+5;
		final IntLongHashMap2 a=new IntLongHashMap2(64);
		final IntLongHashMap3 b=new IntLongHashMap3(64);
		a.put(17, Long.MAX_VALUE-1); b.put(17, Long.MAX_VALUE-1);
		check(a.increment(17, 7)==expected && a.get(17)==expected, "IntLongHashMap2 positive overflow");
		check(b.increment(17, 7)==expected && b.get(17)==expected, "IntLongHashMap3 positive overflow");
		a.put(18, -10); b.put(18, -10);
		check(a.increment(18, 3)==-7 && b.increment(18, 3)==-7, "Negative counts are not positive overflow");
		a.put(19, Long.MAX_VALUE); b.put(19, Long.MAX_VALUE);
		check(a.getAndIncrement(19)==Long.MAX_VALUE && a.get(19)==Long.MIN_VALUE, "Map2 getAndIncrement retains its unclamped contract");
		check(b.getAndIncrement(19)==Long.MAX_VALUE && b.get(19)==Long.MIN_VALUE, "Map3 getAndIncrement retains its unclamped contract");
		if(fixed){checkBoundaryGrid(a, b);}
		System.out.println("PASS\tpositiveOverflow="+expected+"\tfixed="+fixed);
	}

	/** Uses unbounded addition as an oracle for upper saturation and retained lower wraparound. */
	private static void checkBoundaryGrid(final IntLongHashMap2 a, final IntLongHashMap3 b){
		final LongLongHashMap2 c=new LongLongHashMap2(64);
		final long[] values={Long.MIN_VALUE, Long.MIN_VALUE+1, -Long.MAX_VALUE, -1, 0, 1,
			Long.MAX_VALUE-1, Long.MAX_VALUE};
		final BigInteger maximum=BigInteger.valueOf(Long.MAX_VALUE);
		for(long old : values){
			for(long delta : values){
				checkPair(a, b, c, maximum, old, delta);
			}
		}
		final Random random=new Random(7231);
		for(int i=0; i<1024; i++){checkPair(a, b, c, maximum, random.nextLong(), random.nextLong());}
		System.out.println("PASS\tboundaryPairs=64\trandomPairs=1024\tmaps=3");
	}

	/** Compares actual map operations with unbounded arithmetic, independently of the saturation formula. */
	private static void checkPair(final IntLongHashMap2 a, final IntLongHashMap3 b,
		final LongLongHashMap2 c, final BigInteger maximum, final long old, final long delta){
		final long expected=BigInteger.valueOf(old).add(BigInteger.valueOf(delta)).min(maximum).longValue();
		a.put(31, old); b.put(31, old); c.put(31, old);
		check(a.increment(31, delta)==expected && a.get(31)==expected, "Map2 boundary old="+old+", delta="+delta);
		check(b.increment(31, delta)==expected && b.get(31)==expected, "Map3 boundary old="+old+", delta="+delta);
		check(c.increment(31, delta)==expected && c.get(31)==expected, "Long map boundary old="+old+", delta="+delta);
		final long unit=BigInteger.valueOf(old).add(BigInteger.ONE).min(maximum).longValue();
		a.put(31, old); b.put(31, old); c.put(31, old);
		check(a.increment(31)==unit && a.get(31)==unit, "Map2 unit increment old="+old);
		check(b.increment(31)==unit && b.get(31)==unit, "Map3 unit increment old="+old);
		check(c.increment(31)==unit && c.get(31)==unit, "Long map unit increment old="+old);
	}

	private static void check(boolean ok, String message){
		if(!ok){throw new AssertionError(message);}
	}
}
