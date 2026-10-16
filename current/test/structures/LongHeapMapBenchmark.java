package test.structures;

import parse.Parse;
import structures.LongHeapMap;

/** Single-thread benchmark for retained LongHeapMap increment paths.
 * @author Jean
 */
public final class LongHeapMapBenchmark{

	/** Parses controls, runs warmups, and prints measured compute times. */
	public static void main(final String[] args){
		final long startupNanos=System.nanoTime();
		final Config config=new Config(args);
		final long parsedNanos=System.nanoTime();

		if(config.roundRobin){
			runProfile("roundrobin", config);
		}
		if(config.replaceMinimum){
			runProfile("replace", config);
		}
		System.out.println("SETUP\tsetupNanos="+(parsedNanos-startupNanos));
	}

	/** Runs warmup and measured repetitions for one profile. */
	private static void runProfile(final String profile, final Config config){
		for(int i=0; i<config.warmups; i++){
			runOnce(profile, config.iterations, false, i);
		}
		for(int i=0; i<config.reps; i++){
			final Result result=runOnce(profile, config.iterations, true, i);
			System.out.println("TIME\tprofile="+profile+"\trep="+i+
					"\titerations="+config.iterations+"\tnanos="+result.nanos+
					"\tchecksum="+result.checksum);
		}
	}

	/** Dispatches one benchmark run with a fresh map. */
	private static Result runOnce(final String profile, final long iterations,
			final boolean measured, final int rep){
		if("roundrobin".equals(profile)){
			return roundRobin(iterations, measured, rep);
		}
		if("replace".equals(profile)){
			return replaceMinimum(iterations, measured, rep);
		}
		throw new AssertionError("Unknown profile "+profile);
	}

	/** Repeatedly increments 64 retained keys in round-robin order. */
	private static Result roundRobin(final long iterations, final boolean measured,
			final int rep){
		final LongHeapMap map=new LongHeapMap(64);
		for(int i=0; i<64; i++){
			checkEquals("roundrobin seed "+i, 1, map.increment(i, 1));
		}

		long checksum=0;
		final long start=System.nanoTime();
		for(long i=0; i<iterations; i++){
			checksum+=map.increment(i&63, 1);
		}
		final long nanos=System.nanoTime()-start;

		checkRoundRobin(map, iterations, checksum, measured, rep);
		return new Result(nanos, checksum);
	}

	/** Replaces the minimum with successive unique increasing keys at capacity 64. */
	private static Result replaceMinimum(final long iterations, final boolean measured,
			final int rep){
		final LongHeapMap map=new LongHeapMap(64);
		for(int i=0; i<64; i++){
			checkEquals("replace seed "+i, 1, map.increment(i, 1));
		}

		long checksum=0;
		final long start=System.nanoTime();
		for(long i=0; i<iterations; i++){
			checksum+=map.increment(i+64, 1);
		}
		final long nanos=System.nanoTime()-start;

		checkReplace(map, iterations, checksum, measured, rep);
		return new Result(nanos, checksum);
	}

	/** Checks the exact final state after the round-robin profile. */
	private static void checkRoundRobin(final LongHeapMap map, final long iterations,
			final long checksum, final boolean measured, final int rep){
		checkEquals("roundrobin size", 64, map.size());
		checkEquals("roundrobin peek", 0, map.peek());
		final long fullCycles=iterations/64;
		final int extra=(int)(iterations&63);
		long expectedChecksum=0;
		for(int i=0; i<64; i++){
			final int count=1+(int)fullCycles+(i<extra ? 1 : 0);
			check(map.contains(i), "roundrobin missing retained key "+i);
			checkEquals("roundrobin count "+i, count, map.map.get(i));
			final long increments=fullCycles+(i<extra ? 1 : 0);
			expectedChecksum+=increments*(increments+3)/2;
		}
		checkEquals("roundrobin checksum measured="+measured+" rep="+rep,
				expectedChecksum, checksum);
	}

	/** Checks the exact final state after the replacement profile. */
	private static void checkReplace(final LongHeapMap map, final long iterations,
			final long checksum, final boolean measured, final int rep){
		checkEquals("replace size", 64, map.size());
		checkEquals("replace peek", iterations, map.peek());
		for(long key=iterations; key<iterations+64; key++){
			check(map.contains(key), "replace missing retained key "+key);
			checkEquals("replace count "+key, 1, map.map.get(key));
		}
		checkEquals("replace checksum measured="+measured+" rep="+rep,
				iterations, checksum);
	}

	/** Checks an int oracle. */
	private static void checkEquals(final String name, final int expected, final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Checks a long oracle. */
	private static void checkEquals(final String name, final long expected,
			final long observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Immutable parsed benchmark controls. */
	private static final class Config{

		/** Parses flag=value benchmark arguments. */
		Config(final String[] args){
			long iterations_=1000000;
			int warmups_=3, reps_=5;
			boolean roundRobin_=true, replaceMinimum_=true;
			for(String arg : args){
				final String[] split=arg.split("=", 2);
				check(split.length==2, "Expected flag=value, got "+arg);
				final String flag=split[0], value=split[1];
				if("iterations".equals(flag)){
					iterations_=parseLong(flag, value, 0);
				}else if("warmups".equals(flag)){
					warmups_=parseInt(flag, value, 0);
				}else if("reps".equals(flag)){
					reps_=parseInt(flag, value, 1);
				}else if("profile".equals(flag)){
					roundRobin_="both".equals(value) || "roundrobin".equals(value);
					replaceMinimum_="both".equals(value) || "replace".equals(value);
					check(roundRobin_ || replaceMinimum_,
							"profile must be both, roundrobin, or replace");
				}else{
					throw new AssertionError("Unknown flag "+flag);
				}
			}
			check(iterations_<=MAX_ITERATIONS,
					"iterations="+iterations_+" exceeds "+MAX_ITERATIONS);
			iterations=iterations_;
			warmups=warmups_;
			reps=reps_;
			roundRobin=roundRobin_;
			replaceMinimum=replaceMinimum_;
		}

		/** Parses a nonnegative long with a lower bound. */
		private static long parseLong(final String flag, final String value,
				final long min){
			try{
				final long parsed=Parse.parseKMG(value);
				check(parsed>=min, flag+" must be at least "+min);
				return parsed;
			}catch(NumberFormatException e){
				throw new AssertionError(flag+" requires a long: "+value, e);
			}
		}

		/** Parses a nonnegative int with a lower bound. */
		private static int parseInt(final String flag, final String value,
				final int min){
			final long parsed=parseLong(flag, value, min);
			check(parsed<=Integer.MAX_VALUE, flag+" is too large: "+parsed);
			return (int)parsed;
		}

		/** Upper bound keeps retained counts positive and checksum exact. */
		private static final long MAX_ITERATIONS=10000000000L;

		final long iterations;
		final int warmups;
		final int reps;
		final boolean roundRobin;
		final boolean replaceMinimum;
	}

	/** One timed computation and its loop-carried checksum. */
	private static final class Result{

		/** Stores one measured result. */
		Result(final long nanos_, final long checksum_){
			nanos=nanos_;
			checksum=checksum_;
		}

		final long nanos;
		final long checksum;
	}
}
