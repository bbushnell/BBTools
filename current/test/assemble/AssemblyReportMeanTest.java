package test.assemble;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import assemble.SameKAssemblyReport;
import parse.PreParser;
import shared.Shared;

/** Checks the real timing-mean helper against literal raw-bit results and explicit
 * invalid-input errors. Arrays must remain unchanged, including NaN payload bits.
 * Reflection preserves the package-private API. This does not test report parsing,
 * decimal formatting, timing exclusions or the complete reporter CLI.
 * @author Brian Bushnell, Sucrose
 */
public final class AssemblyReportMeanTest{

	/** Accepts only standard launcher flags consumed by PreParser. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, AssemblyReportMeanTest.class, false);
		try{
			if(pp.args.length!=0){throw new IllegalArgumentException("Mean fixture takes no program arguments.");}
			run();
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Captures every returned result before reporting aggregate failure; setup errors propagate. */
	private static void run() throws Exception{
		final Method mean=SameKAssemblyReport.class.getDeclaredMethod("timingMean", double[].class, double.class);
		mean.setAccessible(true);
		check(INPUTS.length==9 && EXPECTED.length==9 && NAMES.length==9,
				"Nine literal valid inputs and independent result oracles are required.");
		check(INVALID.length==9 && INVALID_SUMS.length==9 && ERRORS.length==9 && INVALID_NAMES.length==9,
				"Nine explicit invalid input contracts are required.");
		final StringBuilder report=new StringBuilder("case\tinput_bits\tsum_bits\texpected\tobserved\tunchanged\tpassed\n");
		int failed=0;
		for(int c=0; c<INPUTS.length; c++){
			final double[] values=INPUTS[c].clone();
			double sum=0;
			for(double value : values){
				check(Double.isFinite(value) && value>=0, "Valid mean case has invalid input: "+NAMES[c]);
				sum+=value;
			}
			final double observed;
			try{observed=((Double)mean.invoke(null, values, sum)).doubleValue();}
			catch(InvocationTargetException e){throw new IllegalStateException("Mean unexpectedly rejected "+NAMES[c], e.getCause());}
			final boolean unchanged=equalBits(values, INPUTS[c]);
			final boolean passed=unchanged && Double.doubleToRawLongBits(observed)==EXPECTED[c];
			if(!passed){failed++;}
			report.append(NAMES[c]).append('\t').append(bits(INPUTS[c])).append('\t').append(bits(sum)).append('\t')
				.append(Long.toHexString(EXPECTED[c])).append('\t').append(bits(observed)).append('\t')
				.append(unchanged).append('\t').append(passed).append('\n');
		}
		for(int c=0; c<INVALID.length; c++){
			final double[] values=INVALID[c]==null ? null : INVALID[c].clone();
			String observed="returned";
			try{mean.invoke(null, values, INVALID_SUMS[c]);}
			catch(InvocationTargetException e){
				final Throwable cause=e.getCause();
				if(!(cause instanceof IllegalArgumentException)){
					throw new IllegalStateException("Wrong mean failure type: "+INVALID_NAMES[c], cause);
				}
				observed=cause.getMessage();
			}
			final boolean unchanged=equalBits(values, INVALID[c]);
			final boolean passed=unchanged && ERRORS[c].equals(observed);
			if(!passed){failed++;}
			report.append(INVALID_NAMES[c]).append('\t').append(bits(INVALID[c])).append('\t').append(bits(INVALID_SUMS[c]))
				.append('\t').append(ERRORS[c]).append('\t').append(observed).append('\t')
				.append(unchanged).append('\t').append(passed).append('\n');
		}
		System.out.print(report);
		System.out.println("ASSEMBLY_REPORT_MEAN assertions="+AssemblyReportMeanTest.class.desiredAssertionStatus()+
				" cases=18 failed="+failed);
		System.out.flush();
		if(failed!=0){throw new AssertionError("Assembly report mean contract failed: "+failed+" of 18");}
	}

	/** Compares borrowed inputs without normalizing signed zero or NaN payloads. */
	private static boolean equalBits(final double[] actual, final double[] expected){
		if(actual==null || expected==null){return actual==expected;}
		if(actual.length!=expected.length){return false;}
		for(int i=0; i<actual.length; i++){
			if(Double.doubleToRawLongBits(actual[i])!=Double.doubleToRawLongBits(expected[i])){return false;}
		}
		return true;
	}

	/** Raw IEEE-754 representation avoids the report's decimal formatter. */
	private static String bits(final double value){return Long.toHexString(Double.doubleToRawLongBits(value));}

	/** Distinguishes null, empty and each exact input bit pattern. */
	private static String bits(final double[] values){
		if(values==null){return "null";}
		if(values.length==0){return "empty";}
		final StringBuilder text=new StringBuilder();
		for(int i=0; i<values.length; i++){
			if(i>0){text.append(',');}
			text.append(bits(values[i]));
		}
		return text.toString();
	}

	/** Fixture setup checks stay active when Java assertions are disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException("Mean fixture setup: "+message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private static final String[] NAMES={"ordinary-unequal", "fractional", "all-zero", "zero-large",
			"unequal-finite", "equal-overflow", "unequal-overflow", "maximum-pair", "overflow-with-zeros"};
	/** Templates are cloned; the fixture accumulates their sum, never their expected mean. */
	private static final double[][] INPUTS={{4, 2}, {1, 2, 4}, {0, 0}, {0x1p1023, 0},
			{0x1p1023, 0x1p1022}, {0x1p1023, 0x1p1023}, {0x1p1023, 0x1.8p1023},
			{Double.MAX_VALUE, Double.MAX_VALUE}, {0x1p1023, 0x1p1023, 0, 0}};
	/** Literal means:3,7/3,0,H/2,3H/4,H,5H/4,MAX,H/2, where H=2^1023. */
	private static final long[] EXPECTED={0x4008000000000000L, 0x4002aaaaaaaaaaabL, 0L,
			0x7fd0000000000000L, 0x7fd8000000000000L, 0x7fe0000000000000L,
			0x7fe4000000000000L, 0x7fefffffffffffffL, 0x7fd0000000000000L};
	private static final String[] INVALID_NAMES={"null-array", "empty-array", "negative-value", "nan-value",
			"infinite-value", "nan-sum", "negative-sum", "negative-infinite-sum", "zero-overflow"};
	/** NaN payload is deliberately noncanonical, so the preservation check must use raw bits. */
	private static final double[][] INVALID={null, {}, {-1, 2}, {1, Double.longBitsToDouble(0x7ff8000000000042L)},
			{1, Double.POSITIVE_INFINITY}, {1}, {1}, {1}, {0, 0}};
	private static final double[] INVALID_SUMS={0, 0, 1, 1, Double.POSITIVE_INFINITY,
			Double.NaN, -1, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY};
	/** Exact contract messages distinguish intended rejection from unrelated setup failures. */
	private static final String[] ERRORS={"Timing mean requires a nonempty array.", "Timing mean requires a nonempty array.",
			"Timing mean requires finite nonnegative values: index=0 value=-1.0",
			"Timing mean requires finite nonnegative values: index=1 value=NaN",
			"Timing mean requires finite nonnegative values: index=1 value=Infinity",
			"Timing mean requires a nonnegative sum or positive overflow: NaN",
			"Timing mean requires a nonnegative sum or positive overflow: -1.0",
			"Timing mean requires a nonnegative sum or positive overflow: -Infinity",
			"Timing mean overflow requires a positive maximum."};
}
