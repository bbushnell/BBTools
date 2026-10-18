package test.ml;

import ml.RegressionTrainer;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

/**
 * Direct probes for private relative-error guard helpers in RegressionTrainer fixtures.
 *
 * @author Sucrose
 */
public class RegressionTrainerRelativeGuardTest{

	public static void main(String[] args) throws Exception{
		if(args.length!=1){
			throw new IllegalArgumentException("Expected exactly one case name");
		}
		run(args[0]);
		System.err.println("RegressionTrainerRelativeGuardTest PASS "+args[0]);
	}

	/** Runs one selectable helper case. */
	private static void run(String name) throws Exception{
		switch(name){
			case "simd-finite-within": acceptSimd(name, 1.0, 1.04); break;
			case "simd-finite-beyond": rejectSimd(name, 1.0, 1.06); break;
			case "simd-nan-left": rejectSimd(name, Double.NaN, 1.0); break;
			case "simd-nan-right": rejectSimd(name, 1.0, Double.NaN); break;
			case "simd-both-infinities": rejectSimd(name, Double.POSITIVE_INFINITY,
				Double.POSITIVE_INFINITY); break;
			case "simd-finite-infinity": rejectSimd(name, 1.0,
				Double.POSITIVE_INFINITY); break;
			case "float-finite-within": acceptFloat(name, 1.0, 1.04); break;
			case "float-finite-beyond": rejectFloat(name, 1.0, 1.06); break;
			case "float-nan-left": rejectFloat(name, Double.NaN, 1.0); break;
			case "float-nan-right": rejectFloat(name, 1.0, Double.NaN); break;
			case "float-both-infinities": rejectFloat(name, Double.POSITIVE_INFINITY,
				Double.POSITIVE_INFINITY); break;
			case "float-finite-infinity": rejectFloat(name, 1.0,
				Double.POSITIVE_INFINITY); break;
			case "parallel-finite-within": acceptParallel(name, 1.0, 1.0+5e-10); break;
			case "parallel-finite-beyond": rejectParallel(name, 1.0, 1.0+2e-9); break;
			case "parallel-nan-left": rejectParallel(name, Double.NaN, 1.0); break;
			case "parallel-nan-right": rejectParallel(name, 1.0, Double.NaN); break;
			case "parallel-both-infinities": rejectParallel(name, Double.POSITIVE_INFINITY,
				Double.POSITIVE_INFINITY); break;
			case "parallel-finite-infinity": rejectParallel(name, 1.0,
				Double.POSITIVE_INFINITY); break;
			default: throw new IllegalArgumentException("Unknown case: "+name);
		}
	}

	/** Requires the SIMD consistency guard to accept these values. */
	private static void acceptSimd(String name, double a, double b) throws Exception{
		invokeSimd(a, b);
	}

	/** Requires the SIMD consistency guard to reject these values. */
	private static void rejectSimd(String name, double a, double b) throws Exception{
		reject(name, () -> invokeSimd(a, b));
	}

	/** Requires the float validation guard to accept these values. */
	private static void acceptFloat(String name, double a, double b) throws Exception{
		invokeFloat(a, b);
	}

	/** Requires the float validation guard to reject these values. */
	private static void rejectFloat(String name, double a, double b) throws Exception{
		reject(name, () -> invokeFloat(a, b));
	}

	/** Requires the parallel validation guard to accept these values. */
	private static void acceptParallel(String name, double a, double b) throws Exception{
		invokeParallel(a, b);
	}

	/** Requires the parallel validation guard to reject these values. */
	private static void rejectParallel(String name, double a, double b) throws Exception{
		reject(name, () -> invokeParallel(a, b));
	}

	/** Invokes RegressionTrainerSimdThreadsConsistencyTest's relative helper. */
	private static void invokeSimd(double a, double b) throws Exception{
		final Method m=RegressionTrainerSimdThreadsConsistencyTest.class.getDeclaredMethod(
			"assertRelativeBestValid", String.class, String.class, double.class,
			String.class, double.class);
		m.setAccessible(true);
		m.invoke(null, "SIMD", "scalar", a, "simd-serial", b);
	}

	/** Invokes RegressionTrainerValidateFloatTest's MSE helper. */
	private static void invokeFloat(double a, double b) throws Exception{
		final Method m=RegressionTrainerValidateFloatTest.class.getDeclaredMethod(
			"assertFloatMseClose", double.class, double.class);
		m.setAccessible(true);
		m.invoke(null, a, b);
	}

	/** Invokes RegressionTrainerValidateParallelTest's MSE helper. */
	private static void invokeParallel(double a, double b) throws Exception{
		final Method m=RegressionTrainerValidateParallelTest.class.getDeclaredMethod(
			"assertValidateMseClose", String.class, double.class, double.class);
		m.setAccessible(true);
		m.invoke(null, "direct", a, b);
	}

	/** Requires an AssertionError from one private helper call. */
	private static void reject(String name, ThrowingRunnable runnable) throws Exception{
		try{
			runnable.run();
			throw new AssertionError("Accepted invalid relative guard case: "+name);
		}catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof AssertionError){return;}
			if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
			if(cause instanceof Exception){throw (Exception)cause;}
			throw new RuntimeException(cause);
		}
	}

	/** Runnable that may raise a reflective helper exception. */
	private interface ThrowingRunnable{
		void run() throws Exception;
	}
}
