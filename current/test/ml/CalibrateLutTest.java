package test.ml;

import ml.Calibrate;

import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.util.Arrays;

/**
 * Selectable fixture for Calibrate's isotonic lookup table smoothing.
 * @author Jean
 */
public class CalibrateLutTest{

	public static void main(String[] args) throws Exception{
		if(args.length==1){
			run(args[0]);
			return;
		}
		check(args.length==0, "This fixture accepts only one optional case name");
		run("smoothed_reversal");
		run("already_monotone");
	}

	/** Runs one selectable lookup-table case. */
	private static void run(String name) throws Exception{
		if(name.equals("smoothed_reversal")){
			checkSmoothedReversal();
		}else if(name.equals("already_monotone")){
			checkAlreadyMonotone();
		}else{
			throw new IllegalArgumentException("Unknown fixture mode: "+name);
		}
		System.err.println("CalibrateLutTest PASS "+name);
	}

	/**
	 * Two empirically ordered bins expose that smoothing each PAV block
	 * independently can reverse the finished probabilities.
	 */
	private static void checkSmoothedReversal() throws Exception{
		final Calibrate.Lut lut=build(add(new AccumAccess(16),
			LEFT, 0, 1,
			RIGHT, 1, 10));
		checkKnots("smoothed_reversal", lut, new float[]{LEFT, RIGHT},
			new float[]{2f/13f, 2f/13f});
	}

	/** Already ordered smoothed blocks should remain separate and unchanged. */
	private static void checkAlreadyMonotone() throws Exception{
		final Calibrate.Lut lut=build(add(new AccumAccess(16),
			LEFT, 0, 10,
			RIGHT, 1, 1));
		checkKnots("already_monotone", lut, new float[]{LEFT, RIGHT},
			new float[]{0.5f/11f, 1.5f/2f});
	}

	/** Adds two finite raw-score bins with literal positive and total counts. */
	private static Object add(AccumAccess access, float leftRaw, int leftPos, int leftCount,
			float rightRaw, int rightPos, int rightCount) throws Exception{
		add(access, leftRaw, leftPos, leftCount);
		add(access, rightRaw, rightPos, rightCount);
		return access.finish();
	}

	/** Adds one occupied bin worth of finite raw scores. */
	private static void add(AccumAccess access, float raw, int positives, int count)
			throws Exception{
		check(raw>=0 && raw<1 && raw==raw, "Raw score must be finite in [0,1): "+raw);
		check(positives>=0 && positives<=count, "Invalid positives/count: "+positives+"/"+count);
		for(int i=0; i<positives; i++){access.add(raw, true);}
		for(int i=positives; i<count; i++){access.add(raw, false);}
	}

	/** Runs Calibrate.buildLut on a private Hist snapshot. */
	private static Calibrate.Lut build(Object hist) throws Exception{
		final Method build=Calibrate.class.getDeclaredMethod("buildLut", hist.getClass());
		build.setAccessible(true);
		return (Calibrate.Lut)build.invoke(null, hist);
	}

	/** Verifies literal knot count, raw scores, and probabilities. */
	private static void checkKnots(String name, Calibrate.Lut lut, float[] x, float[] y){
		check(Arrays.equals(x, lut.x) && Arrays.equals(y, lut.y),
			name+" expected x="+Arrays.toString(x)+" y="+Arrays.toString(y)
			+" observed x="+Arrays.toString(lut.x)+" y="+Arrays.toString(lut.y));
		for(int i=1; i<lut.y.length; i++){
			check(lut.y[i-1]<=lut.y[i],
				name+" LUT must be monotone: x="+Arrays.toString(lut.x)
				+" y="+Arrays.toString(lut.y));
		}
	}

	/** Throws a diagnostic assertion failure for a mismatched literal. */
	private static void check(boolean pass, String message){
		if(!pass){throw new AssertionError(message);}
	}

	/** Reflection bridge for Calibrate's private streaming accumulator. */
	private static final class AccumAccess{
		AccumAccess(int bins) throws Exception{
			final Class<?> clazz=Class.forName("ml.Calibrate$Accum");
			final Constructor<?> constructor=clazz.getDeclaredConstructor(int.class);
			constructor.setAccessible(true);
			accum=constructor.newInstance(bins);
			add=clazz.getDeclaredMethod("add", float.class, boolean.class);
			add.setAccessible(true);
			finish=clazz.getDeclaredMethod("finish");
			finish.setAccessible(true);
		}

		void add(float raw, boolean label) throws Exception{
			add.invoke(accum, raw, label);
		}

		Object finish() throws Exception{
			return finish.invoke(accum);
		}

		private final Object accum;
		private final Method add;
		private final Method finish;
	}

	private static final float LEFT=0.03125f;
	private static final float RIGHT=0.09375f;
}
