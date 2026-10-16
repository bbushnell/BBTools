package test.structures;

import structures.CoverageArray2;

/** Baseline and fixed oracles for wholly negative CoverageArray2 range increments.
 * @author Sucrose
 */
public final class CoverageArray2NegativeRangeCase{

	/** Runs current or fixed metadata and unchanged-value checks. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");

		emptyNegativeRange(expect);
		seededNegativeRange(expect);
		partlyNegativeRange(expect);
		positiveRangeControl(expect);
	}

	/** Records metadata separately from values for an empty wholly negative range. */
	private static void emptyNegativeRange(final String expect){
		final CoverageArray2 ca=new CoverageArray2(1, 5);
		ca.incrementRange(-3, -1);

		checkEquals("empty-negative-get0", 0, ca.get(0));
		checkEquals("empty-negative-sum", 0L, ca.sum());
		checkObservation(expect, "empty_negative_values", "0,0",
				ca.get(0)+","+ca.sum());

		final String expected="baseline".equals(expect) ? "0,-1" : Integer.MAX_VALUE+",-1";
		checkObservation(expect, "empty_negative_metadata", expected,
				ca.minIndex+","+ca.maxIndex);
	}

	/** Records metadata separately from values after a real write. */
	private static void seededNegativeRange(final String expect){
		final CoverageArray2 ca=new CoverageArray2(1, 5);
		ca.incrementRange(2, 2);
		ca.incrementRange(-3, -1);

		checkEquals("seeded-negative-get0", 0, ca.get(0));
		checkEquals("seeded-negative-get2", 1, ca.get(2));
		checkObservation(expect, "seeded_negative_values", "0,1",
				ca.get(0)+","+ca.get(2));

		final String expected="baseline".equals(expect) ? "0,2" : "2,2";
		checkObservation(expect, "seeded_negative_metadata", expected,
				ca.minIndex+","+ca.maxIndex);
	}

	/** Verifies partly negative ranges keep clamping to position zero. */
	private static void partlyNegativeRange(final String expect){
		final CoverageArray2 ca=new CoverageArray2(1, 5);
		ca.incrementRange(-3, 1);

		checkObservation(expect, "partly_negative_values", "1,1,0",
				ca.get(0)+","+ca.get(1)+","+ca.get(2));
		checkObservation(expect, "partly_negative_metadata", "0,1",
				ca.minIndex+","+ca.maxIndex);
	}

	/** Verifies ordinary positive ranges still write values and metadata. */
	private static void positiveRangeControl(final String expect){
		final CoverageArray2 ca=new CoverageArray2(1, 5);
		ca.incrementRange(1, 3, 2);

		checkEquals("positive-range-get0", 0, ca.get(0));
		checkEquals("positive-range-get1", 2, ca.get(1));
		checkEquals("positive-range-get2", 2, ca.get(2));
		checkEquals("positive-range-get3", 2, ca.get(3));
		checkEquals("positive-range-get4", 0, ca.get(4));
		checkObservation(expect, "positive_range_values", "0,2,2,2,0",
				ca.get(0)+","+ca.get(1)+","+ca.get(2)+","+ca.get(3)+","+ca.get(4));

		checkObservation(expect, "positive_range_metadata", "1,3",
				ca.minIndex+","+ca.maxIndex);
	}

	/** Checks an int oracle. */
	private static void checkEquals(final String name, final int expected, final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Checks a long oracle. */
	private static void checkEquals(final String name, final long expected, final long observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Checks and prints one observed value. */
	private static void checkObservation(final String expect, final String name,
			final String expected, final String observed){
		check(expected.equals(observed), name+" expected "+expected+", got "+observed);
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+observed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
