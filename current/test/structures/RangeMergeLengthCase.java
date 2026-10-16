package test.structures;

import structures.Range;

/** Literal oracle for Range.merge assertion overflow on representable merged spans.
 * @author Jean
 */
public final class RangeMergeLengthCase{

	private static final int MAX=Integer.MAX_VALUE;

	/** Runs baseline or fixed expectations for the Range.merge length-sum assertion. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		final boolean fixed="fixed".equals(expect);

		checkObservation(expect, "identical_wide",
				fixed ? "0\t1499999999\t1500000000" : "AssertionError",
				observe(0, 1499999999, 0, 1499999999));
		checkObservation(expect, "overlap_wide_left",
				fixed ? "0\t1999999999\t2000000000" : "AssertionError",
				observe(0, 1499999999, 500000000, 1999999999));
		checkObservation(expect, "overlap_wide_right",
				fixed ? "0\t1999999999\t2000000000" : "AssertionError",
				observe(500000000, 1999999999, 0, 1499999999));
		checkObservation(expect, "max_length_contains",
				fixed ? "0\t"+(MAX-1)+"\t"+MAX : "AssertionError",
				observe(0, MAX-1, 0, 0));

		checkObservation(expect, "ordinary_overlap", "10\t25\t16",
				observe(10, 20, 15, 25));
		checkObservation(expect, "ordinary_adjacent", "10\t25\t16",
				observe(10, 20, 21, 25));
		checkObservation(expect, "ordinary_containment", "10\t30\t21",
				observe(10, 30, 15, 20));
		checkObservation(expect, "ordinary_disjoint", "AssertionError",
				observe(10, 20, 22, 25));
	}

	/** Calls Range.merge and returns merged endpoints and length or the thrown type. */
	private static String observe(final int a1, final int b1, final int a2, final int b2){
		try{
			final Range merged=new Range(a1, b1).merge(new Range(a2, b2));
			return merged.a+"\t"+merged.b+"\t"+merged.length();
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Checks and prints a value-or-exception outcome. */
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
