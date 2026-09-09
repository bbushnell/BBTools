package prot;

/** Focused, dependency-free tests for {@link MagQCExpectedCopyFeatures}. */
public class MagQCExpectedCopyFeaturesTest {

	public static void main(String[] args){
		basicIntegerExample();
		fractionalExpectation();
		zeroDenominator();
		nontrackedIsolation();
		invalidInputsFailClosed();
		System.out.println("MagQCExpectedCopyFeaturesTest: PASS");
	}

	private static void basicIntegerExample(){
		String[] names=MagQCExpectedCopyFeatures.featureNames();
		check(names.length==2 && MagQCExpectedCopyFeatures.FEATURE_OBSERVED_EXPECTED.equals(names[0])
			&& MagQCExpectedCopyFeatures.FEATURE_EXCESS_EXPECTED.equals(names[1]), "feature order");
		check("ratio".equals(MagQCExpectedCopyFeatures.FEATURE_UNIT), "feature units");
		// METHODS Part II.B's illustrative seven-item row: E=7, O=5, X=1.
		MagQCExpectedCopyFeatures.Result r=MagQCExpectedCopyFeatures.compute(
			new int[]{1, 0, 0, 2, 3, 0, 0},
			new double[]{0, 1, 0, 2, 3, 0, 1});
		check(close(r.expected, 7), "integer E");
		check(close(r.observedCapped, 5), "integer O");
		check(close(r.excess, 1), "integer X");
		check(close(r.observedExpected, 5.0/7), "integer O/E");
		check(close(r.excessExpected, 1.0/7), "integer X/E");
	}

	private static void fractionalExpectation(){
		// Fractional population means must not be rounded before partitioning.
		MagQCExpectedCopyFeatures.Result r=MagQCExpectedCopyFeatures.compute(
			new int[]{1, 0}, new double[]{0.5, 0.5});
		check(close(r.expected, 1), "fractional E");
		check(close(r.observedCapped, 0.5), "fractional O");
		check(close(r.excess, 0.5), "fractional X");
		check(close(r.observedExpected, 0.5), "fractional O/E");
		check(close(r.excessExpected, 0.5), "fractional X/E");
	}

	private static void zeroDenominator(){
		MagQCExpectedCopyFeatures.Result r=MagQCExpectedCopyFeatures.compute(
			new int[]{0, 4}, new double[]{0, 0});
		check(r.expected==0, "zero E");
		check(r.observedCapped==0, "zero O");
		check(r.excess==4, "zero X");
		check(r.observedExpected==0 && r.excessExpected==0, "E=0 and X>0 must yield O/E=X/E=0");
	}

	private static void nontrackedIsolation(){
		final double[] expected={1, 2, 4};
		final int[] tracked={0, 1};
		MagQCExpectedCopyFeatures.Mapping mapping=MagQCExpectedCopyFeatures.compileMapping(3, tracked, expected);
		MagQCExpectedCopyFeatures.Result a=mapping.compute(new int[]{1, 0, 0});
		MagQCExpectedCopyFeatures.Result b=mapping.compute(new int[]{1, 0, 99});
		check(close(a.observedExpected, b.observedExpected) && close(a.excessExpected, b.excessExpected),
			"nontracked item with nonzero expectation must be isolated");
		final MagQCExpectedCopyFeatures.Result excessOverOne=MagQCExpectedCopyFeatures.compute(
			new int[]{3, 0}, new double[]{0.5, 0.5});
		check(close(excessOverOne.excessExpected, 2.5) && excessOverOne.excessExpected>1,
			"valid X/E greater than one");
		expectMappingFailure(3, expected, new int[]{0, 0}, "duplicate tracked item");
		expectMappingFailure(3, expected, new int[]{3}, "missing tracked item");
	}

	private static void invalidInputsFailClosed(){
		expectFailure(null, new double[0], "null observed");
		expectFailure(new int[0], null, "null expected");
		expectFailure(new int[]{1}, new double[0], "length mismatch");
		expectFailure(new int[]{-1}, new double[]{1}, "negative observed");
		expectFailure(new int[]{1}, new double[]{-1}, "negative expected");
		expectFailure(new int[]{1}, new double[]{Double.NaN}, "NaN expected");
		expectFailure(new int[]{1}, new double[]{Double.POSITIVE_INFINITY}, "infinite expected");
		expectFailure(new int[]{0, 0}, new double[]{Double.MAX_VALUE, Double.MAX_VALUE}, "expected-sum overflow");
		expectFailure(new int[]{1}, new double[]{Double.MIN_VALUE}, "ratio overflow");
	}

	private static void expectFailure(int[] observed, double[] expected, String name){
		try{
			MagQCExpectedCopyFeatures.compute(observed, expected);
			throw new RuntimeException(name+" did not fail");
		}catch(IllegalArgumentException expectedFailure){
			// intended
		}
	}

	private static void expectMappingFailure(int width, double[] expected, int[] indexes, String name){
		try{
			MagQCExpectedCopyFeatures.compileMapping(width, indexes, expected);
			throw new RuntimeException(name+" did not fail");
		}catch(IllegalArgumentException expectedFailure){
			// intended
		}
	}

	private static boolean close(double a, double b){return Math.abs(a-b)<=1e-12;}
	private static void check(boolean ok, String msg){if(!ok){throw new RuntimeException("FAIL: "+msg);}}
}
