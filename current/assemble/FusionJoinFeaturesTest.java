package assemble;

import java.util.Arrays;

/** Deterministic numeric-contract tests invoked by the diagnostic's test mode. @author Fischl */
final class FusionJoinFeaturesTest {

	/** Tests crossing flank distributions, empty regions, true zeros and extreme counts. */
	static void test(){
		final FusionJoinFeatures encoder=new FusionJoinFeatures();
		final int[][] depths={{0, 1, 100, 100}, {8, 16}, {2, 2, 2, 200}, {0, 4, 8}, {4}};
		final int[] sizes={4, 2, 4, 3, 1};
		final float[] forward=encoder.fill(depths, sizes, 32, 33, 100, 200, 3, 7, 13, 2, 5, 3, 1, 2).clone();
		final int[][] reversed={depths[2], depths[1], depths[0], depths[3], depths[4]};
		final float[] reverse=encoder.fill(reversed, sizes, 32, 33, 200, 100, 7, 3, 13, 2, 5, 3, 1, 2);
		check(Arrays.equals(forward, reverse), "Reciprocal join changes the numeric vector");
		final int[][] unequal={{0, 1, 100}, {8, 16}, {2}, {0, 4, 8}, {4}};
		final float[] shortLeft=encoder.fill(unequal, new int[]{3, 2, 1, 3, 1}, 32, 33,
				34, 32, 3, 7, 9, 2, 5, 3, 1, 2).clone();
		final int[][] unequalReverse={unequal[2], unequal[1], unequal[0], unequal[3], unequal[4]};
		check(Arrays.equals(shortLeft, encoder.fill(unequalReverse, new int[]{1, 2, 3, 3, 1}, 32, 33,
				32, 34, 7, 3, 9, 2, 5, 3, 1, 2)), "Unequal flank sample counts break symmetry");
		unequal[0]=new int[0];
		unequalReverse[2]=unequal[0];
		final float[] oneFlank=encoder.fill(unequal, new int[]{0, 2, 1, 3, 1}, 32, 33,
				10, 32, 3, 7, 6, 1, 5, 3, 1, 2).clone();
		check(Arrays.equals(oneFlank, encoder.fill(unequalReverse, new int[]{1, 2, 0, 3, 1}, 32, 33,
				32, 10, 7, 3, 6, 1, 5, 3, 1, 2)), "One absent flank breaks symmetry");
		check(forward.length==107 && FusionJoinFeatures.names().length==107, "Schema width changed");
		check(value(forward, "flank_min_p50_log")==0, "Crossing flank P50 minimum is incorrect");
		check(close(value(forward, "flank_max_p50_log"), .125), "Flank P50 maximum must use depth2");
		check(close(value(forward, "overlap_p0_log"), .375), "Depth8 must encode as0.375");
		check(close(value(forward, "overlap_p100_log"), .5), "Depth16 must encode as0.5");
		check(close(value(forward, "flank_max_zero_fraction"), .25), "Observed zeros were lost");
		check(close(value(forward, "missing_fraction"), 2.0/13), "Missing words have the wrong denominator");
		check(close(value(forward, "alt_expected_fraction"), .2), "Branch fractions have the wrong denominator");
		final int[][] empty={{0}, {}, {1}, {0}, {}};
		final float[] absent=encoder.fill(empty, new int[]{1, 0, 1, 1, 0}, 32, 20,
				32, 32, 0, 0, 3, 2, 0, 0, 0, 0).clone();
		check(value(absent, "overlap_present")==0 && value(absent, "span_present")==0,
				"Empty regions need explicit absent masks");
		check(value(absent, "enrichment_present")==0 && value(absent, "branch_present")==0,
				"Absent evidence must not manufacture enrichment or branch tests");
		empty[1]=new int[]{0};
		final float[] observed=encoder.fill(empty, new int[]{1, 1, 1, 1, 0}, 32, 32,
				32, 32, 0, 0, 4, 3, 0, 0, 0, 0);
		check(value(observed, "overlap_present")==1 && value(observed, "overlap_zero_fraction")==1,
				"Observed zero depth must differ from an absent region");
		empty[1][0]=Integer.MAX_VALUE;
		final float[] large=encoder.fill(empty, new int[]{1, 1, 1, 1, 0}, 32, 32,
				32, 32, 0, 0, 4, 2, 1, 1, 1, Integer.MAX_VALUE);
		for(float value : large){check(Float.isFinite(value), "Extreme counts produced a nonfinite input");}
		check(value(large, "overlap_p100_log")>3.8, "High depth was silently clipped");
		System.err.println("JOIN_FEATURES_TEST_PASS width=107 symmetry empty zero extremes");
	}

	/** Finds a named column so the fixture checks semantics, not duplicate offset arithmetic. */
	private static float value(final float[] vector, final String name){
		final String[] names=FusionJoinFeatures.names();
		for(int i=0; i<names.length; i++){if(names[i].equals(name)){return vector[i];}}
		throw new AssertionError("Missing feature "+name);
	}

	/** Float tolerances cover only the specified log transform's rounding. */
	private static boolean close(final double a, final double b){return Math.abs(a-b)<1e-6;}

	/** Tests fail regardless of the launcher's assertion setting. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
