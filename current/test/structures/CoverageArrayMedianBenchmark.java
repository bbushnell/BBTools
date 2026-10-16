package test.structures;

import parse.Parse;
import shared.Shared;
import structures.CoverageArray;
import structures.CoverageArray2;
import structures.CoverageArray2A;
import structures.CoverageArray3;
import structures.CoverageArray3A;

/** Single-backend benchmark for CoverageArray.median() on deterministic depths.
 * @author Sucrose
 */
public final class CoverageArrayMedianBenchmark{

	/** Runs warmups and one timed block for the requested backend. */
	public static void main(final String[] args){
		final Config config=new Config(args);
		Shared.setThreads(1);
		final CoverageArray coverage=CoverageArray.makeArray(0, LENGTH, config.backend.type);

		long warmupChecksum=0;
		for(int i=0; i<config.warmup; i++){
			warmupChecksum=accumulate(warmupChecksum,
					runOnce(coverage, config.backend, config.expect, false).checksum);
		}

		long checksum=0;
		long nanos=0;
		for(int i=0; i<config.reps; i++){
			final Result result=runOnce(coverage, config.backend, config.expect, true);
			checksum=accumulate(checksum, result.checksum);
			nanos+=result.nanos;
		}

		System.out.println("TIME"
				+"\tbackend="+config.backend.name
				+"\texpect="+config.expect
				+"\tlength="+LENGTH
				+"\twarmup="+config.warmup
				+"\treps="+config.reps
				+"\tnanos="+nanos
				+"\tchecksum="+checksum
				+"\twarmupChecksum="+warmupChecksum);
	}

	/** Restores original depths, times median() alone when requested, then checks. */
	private static Result runOnce(final CoverageArray coverage,
			final Backend backend, final String expect, final boolean measured){
		fillOriginal(coverage);
		coverage.setMedian(-1);
		final long start=(measured ? System.nanoTime() : 0);
		final int median=coverage.median();
		final long nanos=(measured ? System.nanoTime()-start : 0);
		return new Result(nanos, checkExpected(coverage, backend, expect, median));
	}

	/** Fills every coordinate with its deterministic original depth. */
	private static void fillOriginal(final CoverageArray coverage){
		for(int i=0; i<LENGTH; i++){
			coverage.set(i, originalDepth(i));
		}
	}

	/** Checks median(), selected coordinates, and touched bounds. */
	private static long checkExpected(final CoverageArray coverage,
			final Backend backend, final String expect, final int median){
		checkEquals("median", MEDIAN, median);
		final boolean sorted="baseline".equals(expect) && backend.baselineReorders;
		long checksum=median;
		for(int i=0; i<PROBES.length; i++){
			final int probe=PROBES[i];
			final int expected=(sorted ? SORTED_PROBE_DEPTHS : ORIGINAL_PROBE_DEPTHS)[i];
			final int observed=coverage.get(probe);
			checkEquals("depth "+probe, expected, observed);
			checksum=accumulate(checksum, observed);
		}
		checkEquals("minIndex", 0, coverage.minIndex);
		checkEquals("maxIndex", LENGTH-1, coverage.maxIndex);
		checksum=accumulate(checksum, coverage.minIndex);
		checksum=accumulate(checksum, coverage.maxIndex);
		return checksum;
	}

	/** Returns a reproducible depth in the range 0..1023. */
	private static int originalDepth(final int index){
		return (index*1103515245+12345)&DEPTH_MASK;
	}

	/** Mixes one exact check value into a signed long checksum. */
	private static long accumulate(final long checksum, final long value){
		return checksum*1000003L+value;
	}

	/** Checks an exact integer oracle. */
	private static void checkEquals(final String name, final int expected,
			final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Immutable parsed benchmark controls. */
	private static final class Config{

		/** Parses required flag=value benchmark arguments. */
		Config(final String[] args){
			Backend backend_=null;
			String expect_=null;
			int reps_=-1, warmup_=-1;
			for(String arg : args){
				final String[] split=arg.split("=", 2);
				check(split.length==2, "Expected flag=value, got "+arg);
				final String flag=split[0], value=split[1];
				if("backend".equals(flag)){
					backend_=Backend.parse(value);
				}else if("expect".equals(flag)){
					check("baseline".equals(value) || "fixed".equals(value),
							"expect must be baseline or fixed");
					expect_=value;
				}else if("reps".equals(flag)){
					reps_=parseInt(flag, value, 1);
				}else if("warmup".equals(flag)){
					warmup_=parseInt(flag, value, 0);
				}else{
					throw new AssertionError("Unknown flag "+flag);
				}
			}
			check(backend_!=null, "Missing backend=ca2|ca2a|ca3|ca3a");
			check(expect_!=null, "Missing expect=baseline|fixed");
			check(reps_>=1, "Missing reps=");
			check(warmup_>=0, "Missing warmup=");
			backend=backend_;
			expect=expect_;
			reps=reps_;
			warmup=warmup_;
		}

		/** Parses a nonnegative int with a lower bound. */
		private static int parseInt(final String flag, final String value,
				final int min){
			try{
				final long parsed=Parse.parseKMG(value);
				check(parsed>=min, flag+" must be at least "+min);
				check(parsed<=Integer.MAX_VALUE, flag+" is too large: "+parsed);
				return (int)parsed;
			}catch(NumberFormatException e){
				throw new AssertionError(flag+" requires an integer: "+value, e);
			}
		}

		final Backend backend;
		final String expect;
		final int reps;
		final int warmup;
	}

	/** Supported coverage backends and their baseline live-array behavior. */
	private enum Backend{
		CA2("ca2", CoverageArray2.class, true),
		CA2A("ca2a", CoverageArray2A.class, false),
		CA3("ca3", CoverageArray3.class, true),
		CA3A("ca3a", CoverageArray3A.class, false);

		/** Stores one backend descriptor. */
		Backend(final String name_, final Class<? extends CoverageArray> type_,
				final boolean baselineReorders_){
			name=name_;
			type=type_;
			baselineReorders=baselineReorders_;
		}

		/** Parses a command-line backend name. */
		static Backend parse(final String name){
			for(Backend backend : values()){
				if(backend.name.equals(name)){return backend;}
			}
			throw new AssertionError("backend must be ca2, ca2a, ca3, or ca3a");
		}

		final String name;
		final Class<? extends CoverageArray> type;
		final boolean baselineReorders;
	}

	/** One median timing and its checked coordinate fingerprint. */
	private static final class Result{

		/** Stores one measurement. */
		Result(final long nanos_, final long checksum_){
			nanos=nanos_;
			checksum=checksum_;
		}

		final long nanos;
		final long checksum;
	}

	private static final int LENGTH=1048576;
	private static final int DEPTH_MASK=1023;
	private static final int MEDIAN=511;
	private static final int[] PROBES={0, 1, LENGTH/2, LENGTH-2, LENGTH-1};
	private static final int[] ORIGINAL_PROBE_DEPTHS={57, 678, 57, 863, 460};
	private static final int[] SORTED_PROBE_DEPTHS={1023, 1023, 511, 0, 0};

}
