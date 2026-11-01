package test.shared;

import rand.FastRandom;

/** Matched timing and seeded-output checks for valid shared.Random int defaults.
 * A wrapper exposes only FastRandom's primitive methods, leaving int defaults
 * inherited. The untimed oracle uses the independently maintained FastRandom
 * bounded methods; its unbounded low-word output supplies the power-of-two mask.
 * Oracle work is outside the measured loop but included in process wall/CPU time.
 * @author Sucrose
 */
public final class RandomDefaultTimingTest{

	/** Runs mode=power2|odd|shifted with positive rounds/warmup and a fixed positive seed. */
	public static void main(String[] args){
		long rounds=10000000, warmup=1000000;
		String name="power2";
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Expected key=value: "+arg);}
			if(split[0].equals("rounds")){rounds=Long.parseLong(split[1]);}
			else if(split[0].equals("warmup")){warmup=Long.parseLong(split[1]);}
			else if(split[0].equals("mode")){name=split[1];}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		final int mode=name.equals("power2") ? 0 : name.equals("odd") ? 1 : name.equals("shifted") ? 2 : -1;
		if(mode<0 || rounds<1 || warmup<1 || rounds>Long.MAX_VALUE/1024 || warmup>Long.MAX_VALUE/1024){
			throw new IllegalArgumentException("Valid mode and positive checksum-safe rounds/warmup required");
		}
		require(loop(new Defaults(SEED), warmup, mode)==oracle(new FastRandom(SEED), warmup, mode), "warmup seeded checksum");
		final long expected=oracle(new FastRandom(SEED), rounds, mode);
		final Defaults random=new Defaults(SEED);
		final long start=System.nanoTime();
		final long sum=loop(random, rounds, mode);
		final long elapsed=System.nanoTime()-start;
		require(sum==expected, "measured seeded checksum: "+sum+" != "+expected);
		System.out.println("RANDOM_DEFAULT_TIMING mode="+name+" rounds="+rounds+" warmup="+warmup+
				" seed="+SEED+" nanos="+elapsed+" checksum="+sum);
	}

	/** Times only inherited positive-bound calls; no allocation or mode dispatch per draw. */
	private static long loop(Defaults random, long rounds, int mode){
		long sum=0;
		if(mode==0){
			for(long i=0; i<rounds; i++){sum+=random.nextInt(1024);}
		}else if(mode==1){
			for(long i=0; i<rounds; i++){sum+=random.nextInt(1001);}
		}else{
			for(long i=0; i<rounds; i++){sum+=random.nextInt(-100, 900);}
		}
		return sum;
	}

	/** Uses existing concrete-generator methods to check the preserved valid-input stream.
	 * FastRandom's own power-of-two bounded method uses different bits, so that case
	 * explicitly masks its unbounded low-word result as required by the interface default.
	 */
	private static long oracle(FastRandom random, long rounds, int mode){
		long sum=0;
		if(mode==0){
			for(long i=0; i<rounds; i++){sum+=random.nextInt()&1023;}
		}else if(mode==1){
			for(long i=0; i<rounds; i++){sum+=random.nextInt(1001);}
		}else{
			for(long i=0; i<rounds; i++){sum+=random.nextInt(-100, 900);}
		}
		return sum;
	}

	/** Checks expected output independently of JVM assertion settings. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Delegates required primitives while deliberately inheriting every int method. */
	private static final class Defaults implements shared.Random{
		/** Creates the same deterministic core used by the independent reference path. */
		Defaults(long seed){source=new FastRandom(seed);}
		/** Supplies raw words without delegating to any concrete int method. */
		@Override
		public long nextLong(){return source.nextLong();}
		/** Required bounded-long primitive, unused by positive int paths. */
		@Override
		public long nextLong(long bound){return source.nextLong(bound);}
		/** Required seed operation. */
		@Override
		public void setSeed(long seed){source.setSeed(seed);}
		/** Required Gaussian primitive. */
		@Override
		public double nextGaussian(){return source.nextGaussian();}
		/** Private core; no instance is shared between measured/reference/warmup paths. */
		private final FastRandom source;
	}

	/** Positive seed avoids the core's time-seeding convention for negative values. */
	private static final long SEED=20261031L;
}
