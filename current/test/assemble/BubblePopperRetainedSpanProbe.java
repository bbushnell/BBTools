package test.assemble;

import assemble.BubblePopper;
import parse.PreParser;
import shared.Shared;

/**
 * Reproduces the retained-span arithmetic expression extracted from BubblePopper cross-K merges.
 * This scalar probe exercises only the denominator expression and coverage division; it does not
 * construct Contigs or run an end-to-end graph merge.
 * @author Brian Bushnell, Jean
 */
public final class BubblePopperRetainedSpanProbe{

	/** Static fixture entry point only. */
	private BubblePopperRetainedSpanProbe(){}

	/** Runs retained-span denominator checks, or a requested scalar timing mode. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, BubblePopperRetainedSpanProbe.class, false);
		try{
			if(pp.args.length==0){runFunctionalChecks();}
			else{runPerformance(pp.args);}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Runs the three retained-span denominator cases without constructing graph objects. */
	private static void runFunctionalChecks(){
		checkCase("normal", 3, 5, 7.0f, 11.0f, 8, 9.5f, 9.5f);
		checkCase("integerMax", 1073741823, 1073741824, 5.0f, 5.0f,
				Integer.MAX_VALUE, 5.0f, 5.0f);
		checkCase("overflow", 1073741824, 1073741824, 5.0f, 5.0f,
				Integer.MIN_VALUE, -5.0f, 5.0f);
		checkRejects("zero", 0, 0);
		System.out.println("RETAINED_SPAN_ARITHMETIC_PASS");
	}

	/** Runs one no-allocation scalar arithmetic loop after validating every input. */
	private static void runPerformance(final String[] args){
		String mode=null;
		long iterations=-1;
		long warmup=-1;
		for(String arg : args){
			final int equals=arg.indexOf('=');
			require(equals>0 && equals<arg.length()-1,
					"Expected key=value option, observed: "+arg);
			final String key=arg.substring(0, equals);
			final String value=arg.substring(equals+1);
			if("mode".equalsIgnoreCase(key)){
				require("legacy".equals(value) || "fixed".equals(value),
						"mode must be legacy or fixed: "+value);
				mode=value;
			}else if("iterations".equalsIgnoreCase(key)){
				iterations=parseDecimal(key, value);
			}else if("warmup".equalsIgnoreCase(key)){
				warmup=parseDecimal(key, value);
			}else{
				throw new IllegalArgumentException("Unknown option: "+key);
			}
		}
		require(mode!=null, "Performance mode requires mode=legacy or mode=fixed.");
		require(iterations>0, "Performance mode requires positive iterations=.");
		require(warmup>=0, "Performance mode requires nonnegative warmup=.");

		final PerfInputs inputs=makePerfInputs();
		checkPerfInputs(inputs);
		final long warmupChecksum="legacy".equals(mode) ? loopLegacy(inputs, warmup) :
				loopFixed(inputs, warmup);
		final long start=System.nanoTime();
		final long checksum="legacy".equals(mode) ? loopLegacy(inputs, iterations) :
				loopFixed(inputs, iterations);
		final long nanos=System.nanoTime()-start;
		System.out.println("RETAINED_SPAN_PERF_META"+
				"\tmode="+mode+
				"\titerations="+iterations+
				"\twarmup="+warmup+
				"\tseed="+PERF_SEED+
				"\tcases="+PERF_CASES);
		System.out.println("RETAINED_SPAN_PERF_RESULT"+
				"\tnanos="+nanos+
				"\tchecksum="+checksum+
				"\twarmupChecksum="+warmupChecksum);
	}

	/** Checks the legacy expression and the fixed production helper on one scalar case. */
	private static void checkCase(final String name, final int retainedLeft,
			final int retainedRight, final float leftCoverage, final float rightCoverage,
			final int expectedDenominator, final float expectedLegacyCoverage,
			final float expectedFixedCoverage){
		final int denominator=legacyDenominator(retainedLeft, retainedRight);
		require(denominator==expectedDenominator,
				name+" denominator="+denominator+" expected="+expectedDenominator);
		final float coverage=legacyCoverage(retainedLeft, retainedRight, leftCoverage,
				rightCoverage);
		require(Float.floatToIntBits(coverage)==Float.floatToIntBits(expectedLegacyCoverage),
				name+" coverage="+coverage+" expected="+expectedLegacyCoverage+
				" coverageSum="+legacyCoverageSum(retainedLeft, retainedRight,
						leftCoverage, rightCoverage));
		final float fixed=BubblePopper.retainedSpanCoverage(retainedLeft, retainedRight,
				leftCoverage, rightCoverage);
		require(Float.floatToIntBits(fixed)==Float.floatToIntBits(expectedFixedCoverage),
				name+" fixedCoverage="+fixed+" expected="+expectedFixedCoverage);
		System.out.println("RETAINED_SPAN_ARITHMETIC_CASE\t"+name+
				"\tdenominator="+denominator+"\tcoverage="+coverage+
				"\tfixedCoverage="+fixed+"\tcoverageSum="+
				legacyCoverageSum(retainedLeft, retainedRight, leftCoverage, rightCoverage));
	}

	/** Confirms invalid retained-span pairs are rejected instead of yielding NaN. */
	private static void checkRejects(final String name, final int retainedLeft,
			final int retainedRight){
		try{
			BubblePopper.retainedSpanCoverage(retainedLeft, retainedRight, 0, 0);
			throw new IllegalStateException(name+" invalid spans unexpectedly accepted.");
		}catch(IllegalArgumentException expected){
			System.out.println("RETAINED_SPAN_REJECTED\t"+name);
		}
	}

	/** Parses a decimal long without accepting a negative or nonnumeric option. */
	private static long parseDecimal(final String key, final String value){
		try{
			final long parsed=Long.parseLong(value);
			require(parsed>=0, key+" must be nonnegative: "+value);
			return parsed;
		}catch(NumberFormatException nfe){
			throw new IllegalArgumentException(key+" must be a decimal long: "+value, nfe);
		}
	}

	/** Builds deterministic, int-safe retained spans and finite coverages for timing. */
	private static PerfInputs makePerfInputs(){
		final int[] retainedLeft=new int[PERF_CASES];
		final int[] retainedRight=new int[PERF_CASES];
		final float[] leftCoverage=new float[PERF_CASES];
		final float[] rightCoverage=new float[PERF_CASES];
		long state=PERF_SEED;
		for(int i=0; i<PERF_CASES; i++){
			state=nextState(state);
			retainedLeft[i]=1+(int)((state>>>20)&0xFFFFF);
			state=nextState(state);
			retainedRight[i]=1+(int)((state>>>20)&0xFFFFF);
			state=nextState(state);
			leftCoverage[i]=1.0f+(int)((state>>>16)&0xFFFF)/128.0f;
			state=nextState(state);
			rightCoverage[i]=1.0f+(int)((state>>>16)&0xFFFF)/128.0f;
		}
		return new PerfInputs(retainedLeft, retainedRight, leftCoverage, rightCoverage);
	}

	/** Advances the deterministic setup generator. */
	private static long nextState(long state){
		state^=(state<<13);
		state^=(state>>>7);
		state^=(state<<17);
		return state;
	}

	/** Verifies every performance input before any timing loop starts. */
	private static void checkPerfInputs(final PerfInputs inputs){
		for(int i=0; i<PERF_CASES; i++){
			final int retainedLeft=inputs.retainedLeft[i];
			final int retainedRight=inputs.retainedRight[i];
			require(retainedLeft>0 && retainedRight>0,
					"nonpositive timing span at "+i);
			require(legacyDenominator(retainedLeft, retainedRight)>0,
					"nonpositive timing denominator at "+i);
			final float leftCoverage=inputs.leftCoverage[i];
			final float rightCoverage=inputs.rightCoverage[i];
			require(Float.isFinite(leftCoverage) && Float.isFinite(rightCoverage),
					"nonfinite timing coverage at "+i);
			final float legacy=legacyCoverage(retainedLeft, retainedRight,
					leftCoverage, rightCoverage);
			final float fixed=BubblePopper.retainedSpanCoverage(retainedLeft,
					retainedRight, leftCoverage, rightCoverage);
			require(Float.floatToIntBits(legacy)==Float.floatToIntBits(fixed),
					"legacy/fixed mismatch at "+i+": "+legacy+" != "+fixed);
		}
	}

	/** Runs the historical int-denominator expression without allocating. */
	private static long loopLegacy(final PerfInputs inputs, final long iterations){
		long checksum=0;
		for(long i=0; i<iterations; i++){
			final int index=(int)i&PERF_CASE_MASK;
			final float coverage=legacyCoverage(inputs.retainedLeft[index],
					inputs.retainedRight[index], inputs.leftCoverage[index],
					inputs.rightCoverage[index]);
			checksum=updateChecksum(checksum, coverage);
		}
		return checksum;
	}

	/** Runs the widened production helper without allocating. */
	private static long loopFixed(final PerfInputs inputs, final long iterations){
		long checksum=0;
		for(long i=0; i<iterations; i++){
			final int index=(int)i&PERF_CASE_MASK;
			final float coverage=BubblePopper.retainedSpanCoverage(
					inputs.retainedLeft[index], inputs.retainedRight[index],
					inputs.leftCoverage[index], inputs.rightCoverage[index]);
			checksum=updateChecksum(checksum, coverage);
		}
		return checksum;
	}

	/** Mixes one primitive float result into the loop checksum. */
	private static long updateChecksum(final long checksum, final float coverage){
		return checksum*0x9E3779B97F4A7C15L+Float.floatToIntBits(coverage);
	}

	/** Preserves BubblePopper's current float product/addition before widening to double. */
	private static double legacyCoverageSum(final int retainedLeft, final int retainedRight,
			final float leftCoverage, final float rightCoverage){
		return leftCoverage*retainedLeft+rightCoverage*retainedRight;
	}

	/** Returns the historical signed int denominator. */
	private static int legacyDenominator(final int retainedLeft, final int retainedRight){
		return retainedLeft+retainedRight;
	}

	/** Evaluates the historical retained-span coverage expression. */
	private static float legacyCoverage(final int retainedLeft, final int retainedRight,
			final float leftCoverage, final float rightCoverage){
		final double coverageSum=legacyCoverageSum(retainedLeft, retainedRight,
				leftCoverage, rightCoverage);
		return (float)(coverageSum/legacyDenominator(retainedLeft, retainedRight));
	}

	/** Rejects unexpected arithmetic without depending on JVM assertion settings. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Number of precomputed timing cases; a power of two for mask indexing. */
	private static final int PERF_CASES=1024;

	/** Bitmask for cycling through {@link #PERF_CASES} timing cases. */
	private static final int PERF_CASE_MASK=PERF_CASES-1;

	/** Fixed generator seed recorded by performance-mode metadata. */
	private static final long PERF_SEED=0x5eed5eed1234abcdL;

	/** Precomputed retained-span and coverage arrays for scalar timing. */
	private static final class PerfInputs{

		/** Stores deterministic spans and coverages. */
		PerfInputs(final int[] retainedLeft_, final int[] retainedRight_,
				final float[] leftCoverage_, final float[] rightCoverage_){
			retainedLeft=retainedLeft_;
			retainedRight=retainedRight_;
			leftCoverage=leftCoverage_;
			rightCoverage=rightCoverage_;
		}

		/** Positive retained left spans with int-safe pairwise sums. */
		final int[] retainedLeft;

		/** Positive retained right spans with int-safe pairwise sums. */
		final int[] retainedRight;

		/** Finite left coverage values. */
		final float[] leftCoverage;

		/** Finite right coverage values. */
		final float[] rightCoverage;
	}
}
