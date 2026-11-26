package test.jgi;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import jgi.CountDuplicates;

/**
 * Finite direct checks for the production CountDuplicates pair-length helper.
 * Cases use literal boundary lengths and no giant input reads.
 * @author Jean
 */
public final class CountDuplicatesPairLengthTest{

	/** Runs every pair-length fixture and reports each observed and expected sum. */
	public static void main(final String[] args) throws Exception{
		if(args.length!=0){throw new IllegalArgumentException("Expected no arguments");}
		final CountDuplicatesPairLengthTest test=
				new CountDuplicatesPairLengthTest(pairLengthMethod());
		test.check("zero zero", 0L, 0, 0);
		test.check("small", 5L, 2, 3);
		test.check("max zero", 2147483647L, Integer.MAX_VALUE, 0);
		test.check("zero max", 2147483647L, 0, Integer.MAX_VALUE);
		test.check("max one", 2147483648L, Integer.MAX_VALUE, 1);
		test.check("one max", 2147483648L, 1, Integer.MAX_VALUE);
		test.check("half overflow", 2147483648L, 1073741824, 1073741824);
		test.check("max max", 4294967294L, Integer.MAX_VALUE, Integer.MAX_VALUE);
		System.out.println("CountDuplicatesPairLengthTest cases="+test.cases+
				" failures="+test.failures);
		if(test.failures!=0){
			throw new AssertionError("CountDuplicates pair-length failures="+test.failures);
		}
	}

	/** Stores the reflected helper used for exactly one production call per case. */
	private CountDuplicatesPairLengthTest(final Method pairLength_){
		pairLength=pairLength_;
	}

	/** Invokes the production pair-length helper once and compares the observed sum. */
	private void check(final String name, final long expected, final int length1,
			final int length2) throws Exception{
		cases++;
		final long observed=pairLength(length1, length2);
		System.out.println(name+"\tobserved="+observed+"\texpected="+expected);
		if(observed!=expected){failures++;}
	}

	/** Calls the private static helper, preserving its underlying failure cause. */
	private long pairLength(final int length1, final int length2) throws Exception{
		try{
			return ((Long)pairLength.invoke(null, length1, length2)).longValue();
		}catch(final InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof Exception){throw (Exception)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw new RuntimeException(cause);
		}
	}

	/** Finds the private production helper without changing production accessibility. */
	private static Method pairLengthMethod() throws NoSuchMethodException{
		final Method method=CountDuplicates.class.getDeclaredMethod(
				"pairLength", int.class, int.class);
		method.setAccessible(true);
		return method;
	}

	/** Production helper under test. */
	private final Method pairLength;
	/** Number of finite fixtures executed. */
	private int cases=0;
	/** Number of fixture results that missed their independent oracle. */
	private int failures=0;

}
