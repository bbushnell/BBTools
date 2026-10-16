package test.structures;

import structures.Range;

/** Literal scalar cases for Range result-width and intermediate-overflow handling.
 * No case allocates storage proportional to its coordinates.
 * @author Shinobu
 */
public final class RangeResultWidthCase{

	/** Runs the same cases against baseline or repaired Range arithmetic. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final boolean fixed="expect=fixed".equals(args[0]);
		check(fixed || "expect=baseline".equals(args[0]), "Unknown expectation: "+args[0]);
		check(Range.class.desiredAssertionStatus(), "Range assertions must be enabled to check result-width rejection");
		final String mode=fixed ? "fixed" : "baseline";

		length(mode, "ordinary_length", 2, 5, -1, "4");
		length(mode, "negative_length", -5, -2, -1, "4");
		length(mode, "maximum_length", 0, MAX-1, -1, "2147483647");
		length(mode, "clamped_empty", -5, -1, 0, "0");
		length(mode, "clamped_reversed", 10, 20, 5, "-5");
		length(mode, "full_int_length", MIN, MAX, -1, fixed ? ASSERTION : "0");
		length(mode, "one_too_wide_length", 0, MAX, -1, fixed ? ASSERTION : "-2147483648");

		pair(mode, "adjacent_overlap", false, 10, 20, 21, 30, "0");
		pair(mode, "ordinary_overlap", false, 10, 20, 15, 25, "6");
		pair(mode, "reverse_overlap", false, 15, 25, 10, 20, "6");
		pair(mode, "maximum_overlap", false, 0, MAX-1, 0, MAX-1, "2147483647");
		pair(mode, "minimum_gap_result", false, MIN, MIN, 1, 1, "-2147483648");
		pair(mode, "extreme_gap", false, MIN, MIN, MAX, MAX, fixed ? ASSERTION : "2");
		pair(mode, "wide_contained_overlap", false, MIN, MAX, 0, MAX, fixed ? ASSERTION : "-2147483648");
		pair(mode, "wide_reverse_overlap", false, 0, MAX, MIN, MAX, fixed ? ASSERTION : "-2147483648");

		pair(mode, "ordinary_lopsidedness", true, 10, 20, 12, 22, "2");
		pair(mode, "contained_lopsidedness", true, 10, 30, 12, 22, "0");
		pair(mode, "wide_left_intermediate", true, MIN, MAX-1, MAX-1, MAX, fixed ? "1" : "0");
		pair(mode, "wide_right_intermediate", true, MIN, MIN+1, MIN+1, MAX, fixed ? "1" : "0");
		pair(mode, "wide_negative_right", true, MIN, MAX, MIN+1, MIN+1, fixed ? "0" : "1");
		pair(mode, "maximum_lopsidedness", true, MIN, MIN, -1, -1, "2147483647");
		pair(mode, "extreme_lopsidedness", true, MIN, MIN, MAX, MAX, fixed ? ASSERTION : "-1");
	}

	/** Checks length after optional one-sided clamping; negative clamp means none. */
	private static void length(final String mode, final String name, final int a,
			final int b, final int clamp, final String expected){
		final Range range=new Range(a, b);
		if(clamp>=0){range.fixBounds(clamp);}
		String observed;
		try{observed=Integer.toString(range.length());}
		catch(AssertionError e){observed=ASSERTION;}
		result(mode, name, expected, observed);
	}

	/** Checks the selected arithmetic operation, preserving ordered input endpoints. */
	private static void pair(final String mode, final String name, final boolean lopsided,
			final int a1, final int b1, final int a2, final int b2, final String expected){
		check(a1<=b1 && a2<=b2, "Fixture endpoints must be ordered to isolate arithmetic: "+name);
		String observed;
		try{
			observed=Integer.toString(lopsided ? Range.lopsidedness(a1, b1, a2, b2) :
					Range.overlap(a1, b1, a2, b2));
		}catch(AssertionError e){observed=ASSERTION;}
		result(mode, name, expected, observed);
	}

	/** Emits one row only after the independent literal expectation matches. */
	private static void result(final String mode, final String name,
			final String expected, final String observed){
		check(expected.equals(observed), name+" expected "+expected+", observed "+observed);
		System.out.println("PASS\t"+mode+'\t'+name+'\t'+observed);
	}

	/** Keeps the oracle checks active independently of JVM assertion flags. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final int MIN=Integer.MIN_VALUE, MAX=Integer.MAX_VALUE;
	private static final String ASSERTION="AssertionError";
}
