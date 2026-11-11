package test.assemble;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import assemble.SameKAssemblyReport;
import parse.PreParser;
import shared.Shared;

/** Checks the real report median with literal finite inputs and raw-bit expectations.
 * Reflection preserves the package-private production API. Each case owns its array;
 * the sorted-array oracle also verifies the helper's documented in-place mutation.
 * This isolates arithmetic from decimal parsing, output formatting and other summary
 * calculations. It is not a whole-report or performance test.
 * @author Brian Bushnell, Sucrose
 */
public final class AssemblyReportMedianTest{

	/** Runs all seven cases; accepts only standard launcher arguments consumed by PreParser. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, AssemblyReportMedianTest.class, false);
		try{
			if(pp.args.length!=0){throw new IllegalArgumentException("Median fixture takes no program arguments.");}
			run();
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Saves all returned observations before failing the aggregate contract; setup errors are distinct. */
	private static void run() throws Exception{
		final Method median=SameKAssemblyReport.class.getDeclaredMethod("median", double[].class);
		median.setAccessible(true);
		check(INPUTS.length==7 && SORTED.length==7 && EXPECTED_BITS.length==7 && NAMES.length==7,
				"Seven independent input, sorted-array and scalar oracles are required.");
		final StringBuilder report=new StringBuilder("case\tinput_bits\texpected_bits\tobserved_bits\tsorted_bits\tsorted_ok\tpassed\n");
		int failed=0;
		for(int c=0; c<NAMES.length; c++){
			final double[] values=INPUTS[c].clone();
			check(values.length>0, "Report callers supply nonempty arrays: "+NAMES[c]);
			for(double value : values){
				check(Double.isFinite(value) && value>=0, "Report caller domain is finite/nonnegative: "+NAMES[c]);
			}
			final double observed;
			try{observed=((Double)median.invoke(null, (Object)values)).doubleValue();}
			catch(InvocationTargetException e){throw new IllegalStateException("Median helper threw unexpectedly: "+NAMES[c], e.getCause());}
			final long observedBits=Double.doubleToLongBits(observed);
			final boolean sorted=equalBits(values, SORTED[c]);
			final boolean passed=sorted && observedBits==EXPECTED_BITS[c];
			if(!passed){failed++;}
			report.append(NAMES[c]).append('\t').append(bits(INPUTS[c])).append('\t')
				.append(Long.toHexString(EXPECTED_BITS[c])).append('\t').append(Long.toHexString(observedBits)).append('\t')
				.append(bits(values)).append('\t').append(sorted).append('\t').append(passed).append('\n');
		}
		System.out.print(report);
		System.out.println("ASSEMBLY_REPORT_MEDIAN assertions="+AssemblyReportMedianTest.class.desiredAssertionStatus()+
				" cases=7 failed="+failed);
		System.out.flush();
		if(failed!=0){throw new AssertionError("Assembly report median contract failed: "+failed+" of 7");}
	}

	/** Compares to the separately listed literal order, including signed-zero bits. */
	private static boolean equalBits(final double[] actual, final double[] expected){
		if(actual.length!=expected.length){return false;}
		for(int i=0; i<actual.length; i++){
			if(Double.doubleToLongBits(actual[i])!=Double.doubleToLongBits(expected[i])){return false;}
		}
		return true;
	}

	/** Hexadecimal IEEE-754 bits avoid decimal parser and ByteBuilder double formatting. */
	private static String bits(final double[] values){
		final StringBuilder text=new StringBuilder();
		for(int i=0; i<values.length; i++){
			if(i>0){text.append(',');}
			text.append(Long.toHexString(Double.doubleToLongBits(values[i])));
		}
		return text.toString();
	}

	/** Fixture/setup checks remain active with Java assertions disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException("Median fixture setup: "+message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private static final String[] NAMES={"ordinary-even", "ordinary-odd", "zero-pair", "zero-large",
			"equal-large", "equal-maximum", "unequal-large"};
	/** Fresh clones are passed to production; these literal input templates remain unchanged. */
	private static final double[][] INPUTS={{4, 2}, {4, 0, 2}, {0, 0}, {0x1.0p1023, 0},
			{0x1.0p1023, 0x1.0p1023}, {Double.MAX_VALUE, Double.MAX_VALUE}, {0x1.8p1023, 0x1.0p1023}};
	/** Independent expected order, never computed by sorting or averaging in the fixture. */
	private static final double[][] SORTED={{2, 4}, {0, 2, 4}, {0, 0}, {0, 0x1.0p1023},
			{0x1.0p1023, 0x1.0p1023}, {Double.MAX_VALUE, Double.MAX_VALUE}, {0x1.0p1023, 0x1.8p1023}};
	/** Exact results: 3, 2, 0, 2^1022, 2^1023, maximum finite double, and 1.25*2^1023. */
	private static final long[] EXPECTED_BITS={0x4008000000000000L, 0x4000000000000000L, 0L,
			0x7fd0000000000000L, 0x7fe0000000000000L, 0x7fefffffffffffffL, 0x7fe4000000000000L};
}
