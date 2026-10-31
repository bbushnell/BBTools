package test.rand;

import java.lang.reflect.Constructor;
import java.lang.reflect.Method;

import rand.FastRandom;
import rand.FastRandomAES;
import rand.FastRandomXoshiro;
import rand.JavaRandom;
import shared.Random;

/** Checks rand package contracts with finite deterministic cases.
 * @author Sucrose */
public final class RandContractTest{

	/** Prevents instantiation. */
	private RandContractTest(){}

	/** Runs all rand contract checks.
	 * @param args No arguments are expected
	 * @throws Exception If FastRandomSIMD should load but cannot */
	public static void main(final String[] args) throws Exception{
		if(args.length>0){throw new IllegalArgumentException("No arguments expected");}
		int metricGroups=0, skipped=0, failures=0;
		failures+=checkGaussianReset("FastRandom", new RandFactory(){
			@Override
			public Random make(final long seed){return new FastRandom(seed);}
		});
		metricGroups++;
		failures+=checkGaussianReset("FastRandomAES", new RandFactory(){
			@Override
			public Random make(final long seed){return new FastRandomAES(seed);}
		});
		metricGroups++;
		failures+=checkGaussianReset("FastRandomXoshiro", new RandFactory(){
			@Override
			public Random make(final long seed){return new FastRandomXoshiro(seed);}
		});
		metricGroups++;
		final RandFactory simdFactory=simdFactory();
		if(simdFactory==null){
			skipped++;
		}else{
			failures+=checkGaussianReset("FastRandomSIMD", simdFactory);
			metricGroups++;
		}
		failures+=checkJavaRandomSequences(987654321L);
		metricGroups++;
		failures+=checkJavaRandomSequences(-987654321L);
		metricGroups++;
		failures+=checkWideNextIntSequences("FastRandom", new RandFactory(){
			@Override
			public Random make(final long seed){return new FastRandom(seed);}
		});
		metricGroups++;
		failures+=checkWideNextIntSequences("FastRandomAES", new RandFactory(){
			@Override
			public Random make(final long seed){return new FastRandomAES(seed);}
		});
		metricGroups++;
		failures+=checkWideNextIntSequences("FastRandomXoshiro", new RandFactory(){
			@Override
			public Random make(final long seed){return new FastRandomXoshiro(seed);}
		});
		metricGroups++;
		if(simdFactory!=null){
			failures+=checkWideNextIntSequences("FastRandomSIMD", simdFactory);
			metricGroups++;
		}
		failures+=checkWideNextIntSequences("JavaRandom", new RandFactory(){
			@Override
			public Random make(final long seed){
				final JavaRandom randy=new JavaRandom();
				randy.setSeed(seed);
				return randy;
			}
		});
		metricGroups++;
		System.out.println("RAND_CONTRACT_RESULT metricGroups="+metricGroups+
				" skipped="+skipped+" failures="+failures);
		if(failures>0){throw new AssertionError("rand contract failures="+failures);}
	}

	/** Verifies setSeed clears the cached second Box-Muller Gaussian.
	 * @param label Diagnostic generator label
	 * @param factory Generator factory
	 * @return One on failure, otherwise zero */
	private static int checkGaussianReset(final String label, final RandFactory factory){
		try{
			final Random randy=factory.make(123456789L);
			final double first=randy.nextGaussian();
			randy.setSeed(123456789L);
			final double afterReset=randy.nextGaussian();
			if(Double.doubleToLongBits(first)==Double.doubleToLongBits(afterReset)){return 0;}
			System.out.println("FAIL gaussianReset label="+label+" first="+
					Double.toHexString(first)+" afterReset="+
					Double.toHexString(afterReset));
		}catch(Throwable t){
			System.out.println("FAIL gaussianReset label="+label+" exception="+t);
			t.printStackTrace(System.out);
		}
		return 1;
	}

	/** Verifies JavaRandom keeps Java's seeded nextInt/nextLong sequence.
	 * @param seed Seed to compare
	 * @return One on failure, otherwise zero */
	private static int checkJavaRandomSequences(final long seed){
		try{
			int failures=checkJavaRandomRawSequences(seed);
			final Method nextIntRange, nextLongBound;
			try{
				nextIntRange=java.util.Random.class.getMethod("nextInt",
						int.class, int.class);
				nextLongBound=java.util.Random.class.getMethod("nextLong",
						long.class);
			}catch(NoSuchMethodException nsme){
				System.out.println("SKIP javaRandom JDK bounded-overload oracle "+
						"seed="+seed+" exception="+nsme);
				return failures;
			}
			failures+=checkJavaRandomNextIntRangeSequence(seed, nextIntRange,
					17, 1000020);
			failures+=checkJavaRandomNextIntRangeSequence(seed, nextIntRange,
					17, 1048593);
			failures+=checkJavaRandomNextIntRangeSequence(seed, nextIntRange,
					Integer.MIN_VALUE, Integer.MAX_VALUE);
			failures+=checkJavaRandomNextIntRangeSequence(seed, nextIntRange,
					Integer.MIN_VALUE, Integer.MAX_VALUE-1);
			failures+=checkJavaRandomNextLongBoundSequence(seed, nextLongBound,
					1024L);
			failures+=checkJavaRandomNextLongBoundSequence(seed, nextLongBound,
					1000000007L);
			return failures;
		}catch(Throwable t){
			System.out.println("FAIL javaRandom sequence seed="+seed+" exception="+t);
			t.printStackTrace(System.out);
		}
		return 1;
	}

	/** Verifies JavaRandom keeps Java's raw seeded nextInt/nextLong sequence.
	 * @param seed Seed to compare
	 * @return One on failure, otherwise zero */
	private static int checkJavaRandomRawSequences(final long seed){
		final JavaRandom randy=new JavaRandom();
		final java.util.Random expected=new java.util.Random(seed);
		randy.setSeed(seed);
		for(int i=0; i<16; i++){
			final int a=randy.nextInt(), b=expected.nextInt();
			if(a!=b){
				System.out.println("FAIL javaRandom nextInt seed="+seed+
						" index="+i+" observed="+a+" expected="+b);
				return 1;
			}
			final long c=randy.nextLong(), d=expected.nextLong();
			if(c!=d){
				System.out.println("FAIL javaRandom nextLong seed="+seed+
						" index="+i+" observed="+c+" expected="+d);
				return 1;
			}
		}
		return 0;
	}

	/** Verifies one JavaRandom nextInt(origin,bound) sequence by reflection.
	 * @param seed Seed to compare
	 * @param method java.util.Random.nextInt(int,int) method from the running JDK
	 * @param origin Lower bound, inclusive
	 * @param bound Upper bound, exclusive
	 * @return One on failure, otherwise zero
	 * @throws Exception If reflective method invocation fails */
	private static int checkJavaRandomNextIntRangeSequence(final long seed,
			final Method method, final int origin, final int bound) throws Exception{
		final JavaRandom randy=new JavaRandom();
		final java.util.Random expected=new java.util.Random(seed);
		randy.setSeed(seed);
		for(int i=0; i<16; i++){
			final int observed=randy.nextInt(origin, bound);
			final int value=((Integer)method.invoke(expected,
					Integer.valueOf(origin), Integer.valueOf(bound))).intValue();
			if(observed!=value){
				System.out.println("FAIL javaRandom nextIntRange seed="+seed+
						" origin="+origin+" bound="+bound+" index="+i+
						" observed="+observed+" expected="+value);
				return 1;
			}
		}
		return 0;
	}

	/** Verifies one JavaRandom nextLong(bound) sequence by reflection.
	 * @param seed Seed to compare
	 * @param method java.util.Random.nextLong(long) method from the running JDK
	 * @param bound Upper bound, exclusive
	 * @return One on failure, otherwise zero
	 * @throws Exception If reflective method invocation fails */
	private static int checkJavaRandomNextLongBoundSequence(final long seed,
			final Method method, final long bound) throws Exception{
		final JavaRandom randy=new JavaRandom();
		final java.util.Random expected=new java.util.Random(seed);
		randy.setSeed(seed);
		for(int i=0; i<16; i++){
			final long observed=randy.nextLong(bound);
			final long value=((Long)method.invoke(expected,
					Long.valueOf(bound))).longValue();
			if(observed!=value){
				System.out.println("FAIL javaRandom nextLongBound seed="+seed+
						" bound="+bound+" index="+i+" observed="+observed+
						" expected="+value);
				return 1;
			}
		}
		return 0;
	}

	/** Verifies wide nextInt(origin,bound) sequences against raw rejection oracles.
	 * @param label Diagnostic generator label
	 * @param factory Generator factory
	 * @return One on failure, otherwise zero */
	private static int checkWideNextIntSequences(final String label, final RandFactory factory){
		int failures=0;
		failures+=checkNextIntRangeSequence(label, factory, 314159265L,
				Integer.MIN_VALUE, Integer.MAX_VALUE);
		failures+=checkNextIntRangeSequence(label, factory, 271828182L,
				Integer.MIN_VALUE, Integer.MAX_VALUE-1);
		return failures;
	}

	/** Verifies one nextInt(origin,bound) sequence against raw rejection sampling.
	 * @param label Diagnostic generator label
	 * @param factory Generator factory
	 * @param seed Seed for the actual and oracle generators
	 * @param origin Lower bound, inclusive
	 * @param bound Upper bound, exclusive
	 * @return One on failure, otherwise zero */
	private static int checkNextIntRangeSequence(final String label,
			final RandFactory factory, final long seed, final int origin, final int bound){
		try{
			final Random actual=factory.make(seed);
			final Random raw=factory.make(seed);
			for(int i=0; i<16; i++){
				final int observed=actual.nextInt(origin, bound);
				final int expected=nextIntByRejection(raw, origin, bound);
				if(observed!=expected){
					System.out.println("FAIL rangeSequence label="+label+" seed="+seed+
							" origin="+origin+" bound="+bound+" index="+i+
							" observed="+observed+" expected="+expected);
					return 1;
				}
			}
			return 0;
		}catch(Throwable t){
			System.out.println("FAIL rangeSequence label="+label+" seed="+seed+
					" origin="+origin+" bound="+bound+" exception="+t);
			t.printStackTrace(System.out);
		}
		return 1;
	}

	/** Computes nextInt(origin,bound) with rejection sampling from raw nextInt().
	 * @param randy Generator for raw nextInt() values
	 * @param origin Lower bound, inclusive
	 * @param bound Upper bound, exclusive
	 * @return Next accepted value */
	private static int nextIntByRejection(final Random randy, final int origin,
			final int bound){
		int r=randy.nextInt();
		while(r<origin || r>=bound){
			r=randy.nextInt();
		}
		return r;
	}

	/** Builds a FastRandomSIMD factory when the Vector module is present.
	 * @return Factory, or null when the Vector module is absent
	 * @throws Exception If FastRandomSIMD cannot load with the Vector module present */
	private static RandFactory simdFactory() throws Exception{
		try{
			Class.forName("jdk.incubator.vector.LongVector");
		}catch(ClassNotFoundException | LinkageError e){
			return null;
		}
		final Class<?> clazz=Class.forName("rand.FastRandomSIMD");
		final Constructor<?> constructor=clazz.getConstructor(long.class);
		return new RandFactory(){
			@Override
			public Random make(final long seed) throws Exception{
				return (Random)constructor.newInstance(seed);
			}
		};
	}

	/** Constructs seeded shared.Random instances. */
	private interface RandFactory{
		/** Returns a generator seeded with the requested value.
		 * @param seed Seed value
		 * @return Seeded generator
		 * @throws Exception If reflective construction fails */
		Random make(long seed) throws Exception;
	}

}
