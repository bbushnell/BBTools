package test.ml;


import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

/**
 * Selectable direct-call fixture for RegressionTrainerRowTest.assertClose().
 *
 * @author Sucrose
 */
public final class RegressionTrainerAssertCloseTest{

	/** Runs one named assertClose() case. */
	public static void main(String[] args) throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Expected one case name");}
		run(args[0]);
		System.err.println("RegressionTrainerAssertCloseTest PASS "+args[0]);
	}

	/** Dispatches the requested direct-call comparison case. */
	private static void run(final String name) throws Exception{
		if(name.equals("finite-equal")){
			accept(name, 1, 1);
		}else if(name.equals("within1e-12")){
			accept(name, 1, 1+5e-13);
		}else if(name.equals("beyond1e-12")){
			reject(name, 1, 1+2e-12);
		}else if(name.equals("nan-actual")){
			reject(name, 1, Double.NaN);
		}else if(name.equals("nan-expected")){
			reject(name, Double.NaN, 1);
		}else if(name.equals("both-infinities")){
			reject(name, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
		}else if(name.equals("unequal-infinities")){
			reject(name, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY);
		}else{
			throw new IllegalArgumentException("Unknown case: "+name);
		}
	}

	/** Requires the private helper to accept one comparison. */
	private static void accept(final String name, final double expected,
		final double actual) throws Exception{
		invoke(name, expected, actual);
	}

	/** Requires the private helper to reject one comparison. */
	private static void reject(final String name, final double expected,
		final double actual) throws Exception{
		try{
			invoke(name, expected, actual);
		}catch(AssertionError expectedFailure){
			return;
		}
		throw new AssertionError("Accepted invalid assertClose case: "+name);
	}

	/** Invokes RegressionTrainerRowTest.assertClose() through reflection. */
	private static void invoke(final String name, final double expected,
		final double actual) throws Exception{
		try{
			ASSERT_CLOSE.invoke(null, expected, actual, name);
		}catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof AssertionError){throw (AssertionError)cause;}
			if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw e;
		}
	}

	/** Looks up the private helper once so every case tests the same method. */
	private static Method findAssertClose(){
		try{
			final Method method=RegressionTrainerRowTest.class.getDeclaredMethod(
				"assertClose", double.class, double.class, String.class);
			method.setAccessible(true);
			return method;
		}catch(Exception e){
			throw new AssertionError("Could not reflect RegressionTrainerRowTest.assertClose", e);
		}
	}

	private RegressionTrainerAssertCloseTest(){}

	private static final Method ASSERT_CLOSE=findAssertClose();

}
