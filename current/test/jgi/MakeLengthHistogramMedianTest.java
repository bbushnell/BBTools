package test.jgi;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import jgi.MakeLengthHistogram;

/**
 * Finite direct checks for the production MakeLengthHistogram median helper.
 * Cases use independent expected bin indexes and no CLI-scale generated reads.
 * @author Jean
 */
public final class MakeLengthHistogramMedianTest{

	/** Runs every median-bin fixture and reports each observed and expected bin. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=0){throw new IllegalArgumentException("Expected no arguments");}
		final MakeLengthHistogramMedianTest test=new MakeLengthHistogramMedianTest(medianMethod());
		test.check("empty", 0);
		test.check("empty bins", 0, 0, 0);
		test.check("small lower even", 0, 1, 1);
		test.check("small odd", 1, 1, 1, 1);
		test.check("near Long.MAX exact odd", 1, (1L<<62)-1, 1, (1L<<62)-1);
		test.check("beyond Long.MAX delegated", 1, 0, Long.MAX_VALUE, Long.MAX_VALUE);
		test.check("large even round-down boundary", 1, 1L<<53, 2, 1L<<53);
		test.check("large even round-up boundary", 0, (1L<<53)+3, 1, (1L<<53)+2);
		System.out.println("MakeLengthHistogramMedianTest cases="+test.cases+
			" failures="+test.failures);
		if(test.failures!=0){
			throw new AssertionError("MakeLengthHistogram median failures="+test.failures);
		}
	}

	/** Stores the reflected helper used for exactly one production call per case. */
	private MakeLengthHistogramMedianTest(final Method medianBin_){
		medianBin=medianBin_;
	}

	/** Invokes the production median helper once and compares the observed bin. */
	private void check(final String name, final int expected, final long... histogram)
			throws Exception{
		cases++;
		final int observed=median(histogram);
		System.out.println(name+"\tobserved="+observed+"\texpected="+expected);
		if(observed!=expected){failures++;}
	}

	/** Calls the private static helper, preserving its underlying failure cause. */
	private int median(final long[] histogram) throws Exception{
		try{
			return ((Integer)medianBin.invoke(null, new Object[]{histogram})).intValue();
		}catch(final InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof Exception){throw (Exception)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw new RuntimeException(cause);
		}
	}

	/** Finds the private production helper without changing production accessibility. */
	private static Method medianMethod() throws NoSuchMethodException{
		final Method method=MakeLengthHistogram.class.getDeclaredMethod(
			"medianBin", long[].class);
		method.setAccessible(true);
		return method;
	}

	/** Production helper under test. */
	private final Method medianBin;
	/** Number of finite fixtures executed. */
	private int cases=0;
	/** Number of fixture results that missed their independent oracle. */
	private int failures=0;

}
