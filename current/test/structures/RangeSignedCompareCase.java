package test.structures;

import structures.Range;

/** Finite oracle for signed Range.compareTo coordinate comparisons.
 * @author Jean
 */
public final class RangeSignedCompareCase{

	/** Runs baseline or fixed expectations for signed endpoint comparisons. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		final boolean fixed="fixed".equals(expect);

		checkObservation(expect, "small_start_low", "-2",
				observeRaw(new Range(1, 9), new Range(3, 4)));
		checkObservation(expect, "small_start_high", "2",
				observeRaw(new Range(3, 4), new Range(1, 9)));
		checkObservation(expect, "small_end_low", "-2",
				observeRaw(new Range(5, 7), new Range(5, 9)));
		checkObservation(expect, "small_end_high", "2",
				observeRaw(new Range(5, 9), new Range(5, 7)));
		checkObservation(expect, "equal_negative", "0",
				observeRaw(new Range(-7, -3), new Range(-7, -3)));
		checkObservation(expect, "equal_extrema", "0",
				observeRaw(new Range(Integer.MIN_VALUE, Integer.MAX_VALUE),
						new Range(Integer.MIN_VALUE, Integer.MAX_VALUE)));

		checkObservation(expect, "start_min_overflow",
				fixed ? "negative" : "positive",
				observeSign(new Range(Integer.MIN_VALUE, Integer.MIN_VALUE+1),
						new Range(1, 2)));
		checkObservation(expect, "start_max_overflow",
				fixed ? "positive" : "negative",
				observeSign(new Range(Integer.MAX_VALUE-1, Integer.MAX_VALUE),
						new Range(-2, -1)));
		checkObservation(expect, "end_min_overflow",
				fixed ? "negative" : "positive",
				observeSign(new Range(Integer.MIN_VALUE, Integer.MIN_VALUE+1),
						new Range(Integer.MIN_VALUE, Integer.MAX_VALUE)));
		checkObservation(expect, "end_max_overflow",
				fixed ? "positive" : "negative",
				observeSign(new Range(Integer.MIN_VALUE, Integer.MAX_VALUE),
						new Range(Integer.MIN_VALUE, Integer.MIN_VALUE+1)));

		checkObservation(expect, "antisymmetry_start", fixed ? "true" : "false",
				antisymmetry(new Range(Integer.MAX_VALUE-1, Integer.MAX_VALUE),
						new Range(-2, -1)));
		checkObservation(expect, "transitivity_start", fixed ? "true" : "false",
				descendingTransitivity(new Range(Integer.MAX_VALUE-1, Integer.MAX_VALUE),
						new Range(0, 1), new Range(-2, -1)));
		checkObservation(expect, "transitivity_end", fixed ? "true" : "false",
				descendingTransitivity(new Range(Integer.MIN_VALUE, Integer.MAX_VALUE),
						new Range(Integer.MIN_VALUE, 0),
						new Range(Integer.MIN_VALUE, Integer.MIN_VALUE+1)));
	}

	/** Calls Range.compareTo and returns the exact integer result. */
	private static String observeRaw(final Range left, final Range right){
		return Integer.toString(left.compareTo(right));
	}

	/** Calls Range.compareTo and returns only the sign of the result. */
	private static String observeSign(final Range left, final Range right){
		final int x=left.compareTo(right);
		return x<0 ? "negative" : x>0 ? "positive" : "zero";
	}

	/** Checks whether comparison signs reverse when argument order is reversed. */
	private static String antisymmetry(final Range left, final Range right){
		final int x=left.compareTo(right);
		final int y=right.compareTo(left);
		return Boolean.toString(Integer.signum(x)==-Integer.signum(y));
	}

	/** Checks high>mid and mid>low implying high>low. */
	private static String descendingTransitivity(final Range high, final Range mid,
			final Range low){
		return Boolean.toString(high.compareTo(mid)>0 && mid.compareTo(low)>0
				&& high.compareTo(low)>0);
	}

	/** Checks and prints a value outcome. */
	private static void checkObservation(final String expect, final String name,
			final String expected, final String observed){
		check(expected.equals(observed), name+" expected "+expected+", got "+observed);
		System.out.println("PASS	"+expect+'	'+name+'	'+observed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
