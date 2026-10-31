package test.rand;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;

import parse.Parse;
import rand.FastRandom;
import rand.FastRandomAES;
import rand.FastRandomXoshiro;
import rand.JavaRandom;
import shared.Random;

/** Runs single-threaded rand hot-loop benchmarks with retained checksums.
 * @author Sucrose */
public final class RandBenchmark{

	/** Ordinary bounded-int upper bound. */
	private static final int INT_BOUND=1_000_003;
	/** Ordinary ranged-int lower bound. */
	private static final int INT_RANGE_ORIGIN=17;
	/** Ordinary ranged-int upper bound. */
	private static final int INT_RANGE_BOUND=1_000_020;
	/** Ordinary bounded-long upper bound. */
	private static final long LONG_BOUND=1_000_000_007L;
	/** Default number of random values to draw. */
	private static final long DEFAULT_OPERATIONS=100_000_000L;
	/** Default positive deterministic seed. */
	private static final long DEFAULT_SEED=123456789L;

	/** Prevents instantiation. */
	private RandBenchmark(){}

	/** Runs the selected benchmark.
	 * @param args key=value options: kind, mode, operations, seed
	 * @throws Exception If reflective SIMD construction fails */
	public static void main(final String[] args) throws Exception{
		String kind="fast", mode="intBounded";
		long operations=DEFAULT_OPERATIONS, seed=DEFAULT_SEED;
		for(final String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<0){throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String key=arg.substring(0, equals).toLowerCase();
			final String value=arg.substring(equals+1);
			if("kind".equals(key)){kind=value;}
			else if("mode".equals(key)){mode=value;}
			else if("operations".equals(key)){operations=Parse.parseKMG(value);}
			else if("seed".equals(key)){seed=Parse.parseKMG(value);}
			else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		if(operations<=0){throw new IllegalArgumentException("operations must be positive");}
		if(seed<0){throw new IllegalArgumentException("seed must be nonnegative");}

		final Random randy=make(kind, seed);
		final ThreadMXBean bean=ManagementFactory.getThreadMXBean();
		final boolean cpu=bean.isCurrentThreadCpuTimeSupported() &&
				bean.isThreadCpuTimeEnabled();
		final long startCpu=cpu ? bean.getCurrentThreadCpuTime() : -1;
		final long start=System.nanoTime();
		final long checksum=run(randy, mode, operations);
		final long wallNanos=System.nanoTime()-start;
		final long cpuNanos=cpu ? bean.getCurrentThreadCpuTime()-startCpu : -1;
		System.out.println("RAND_BENCHMARK_RESULT kind="+kind+" mode="+mode+
				" seed="+seed+" operations="+operations+" checksum="+checksum+
				" wallNanos="+wallNanos+" cpuNanos="+cpuNanos);
	}

	/** Runs the benchmark loop.
	 * @param randy Generator to exercise
	 * @param mode Benchmark mode name
	 * @param operations Iteration count
	 * @return Retained checksum */
	private static long run(final Random randy, final String mode, final long operations){
		long checksum=0;
		if("intBounded".equals(mode)){
			for(long i=0; i<operations; i++){
				checksum+=randy.nextInt(INT_BOUND);
			}
		}else if("intRange".equals(mode)){
			for(long i=0; i<operations; i++){
				checksum+=randy.nextInt(INT_RANGE_ORIGIN, INT_RANGE_BOUND);
			}
		}else if("longBounded".equals(mode)){
			for(long i=0; i<operations; i++){
				checksum^=randy.nextLong(LONG_BOUND);
				checksum=Long.rotateLeft(checksum, 1);
			}
		}else{
			throw new IllegalArgumentException("Unknown mode: "+mode);
		}
		return checksum;
	}

	/** Constructs a generator.
	 * @param kind fast, xoshiro, aes, java, or simd
	 * @param seed Positive seed
	 * @return Seeded generator
	 * @throws Exception If reflective SIMD construction fails */
	private static Random make(final String kind, final long seed) throws Exception{
		if("fast".equals(kind)){return new FastRandom(seed);}
		if("xoshiro".equals(kind)){return new FastRandomXoshiro(seed);}
		if("aes".equals(kind)){return new FastRandomAES(seed);}
		if("java".equals(kind)){
			final JavaRandom randy=new JavaRandom();
			randy.setSeed(seed);
			return randy;
		}
		if("simd".equals(kind)){
			return (Random)Class.forName("rand.FastRandomSIMD").
					getConstructor(long.class).newInstance(seed);
		}
		throw new IllegalArgumentException("Unknown kind: "+kind);
	}

}
