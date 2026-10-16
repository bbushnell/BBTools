package test.structures;

import structures.Range;

/** Literal endpoint-overflow oracle for static Range.touch and Range.adjacent.
 * @author Sucrose
 */
public final class RangeEndpointCase{

	private static final int MIN=Integer.MIN_VALUE;
	private static final int MAX=Integer.MAX_VALUE;

	/** Runs baseline or fixed expectations for endpoint-overflow cases. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		final boolean fixed="fixed".equals(expect);

		checkObservation(expect, "touch_point_ordinary_contained",
				"true", observeTouchPoint(10, 20, 15));
		checkObservation(expect, "touch_point_ordinary_adjacent",
				"true", observeTouchPoint(10, 20, 21));
		checkObservation(expect, "touch_point_ordinary_separated",
				"false", observeTouchPoint(10, 20, 22));
		checkObservation(expect, "touch_point_max_contained",
				fixed ? "true" : "false", observeTouchPoint(0, MAX, MAX));
		checkObservation(expect, "touch_point_max_neighbor",
				"true", observeTouchPoint(MAX-3, MAX-1, MAX));
		checkObservation(expect, "touch_point_min_contained",
				fixed ? "true" : "false", observeTouchPoint(MIN, 0, MIN));
		checkObservation(expect, "touch_point_min_neighbor",
				"true", observeTouchPoint(MIN+1, MIN+3, MIN));

		checkObservation(expect, "touch_range_ordinary_overlap",
				"true", observeTouchRange(10, 20, 15, 25));
		checkObservation(expect, "touch_range_ordinary_adjacent",
				"true", observeTouchRange(10, 20, 21, 25));
		checkObservation(expect, "touch_range_ordinary_separated",
				"false", observeTouchRange(10, 20, 22, 25));
		checkObservation(expect, "touch_range_max_contained",
				fixed ? "true" : "false", observeTouchRange(0, MAX, MAX, MAX));
		checkObservation(expect, "touch_range_max_neighbor",
				"true", observeTouchRange(MAX-3, MAX-1, MAX, MAX));
		checkObservation(expect, "touch_range_min_contained",
				fixed ? "true" : "false", observeTouchRange(MIN, 0, MIN, MIN));
		checkObservation(expect, "touch_range_min_neighbor",
				"true", observeTouchRange(MIN+1, MIN+3, MIN, MIN));

		checkObservation(expect, "adjacent_point_ordinary_low",
				"true", observeAdjacentPoint(10, 20, 9));
		checkObservation(expect, "adjacent_point_ordinary_contained",
				"false", observeAdjacentPoint(10, 20, 10));
		checkObservation(expect, "adjacent_point_ordinary_separated",
				"false", observeAdjacentPoint(10, 20, 8));
		checkObservation(expect, "adjacent_point_max_wrap",
				fixed ? "false" : "true", observeAdjacentPoint(0, MAX, MIN));
		checkObservation(expect, "adjacent_point_max_neighbor",
				"true", observeAdjacentPoint(MAX-3, MAX-1, MAX));
		checkObservation(expect, "adjacent_point_min_wrap",
				fixed ? "false" : "true", observeAdjacentPoint(MIN, 0, MAX));
		checkObservation(expect, "adjacent_point_min_neighbor",
				"true", observeAdjacentPoint(MIN+1, MIN+3, MIN));

		checkObservation(expect, "adjacent_range_ordinary_high",
				"true", observeAdjacentRange(10, 20, 21, 25));
		checkObservation(expect, "adjacent_range_ordinary_overlap",
				"false", observeAdjacentRange(10, 20, 20, 25));
		checkObservation(expect, "adjacent_range_ordinary_separated",
				"false", observeAdjacentRange(10, 20, 22, 25));
		checkObservation(expect, "adjacent_range_max_wrap",
				fixed ? "false" : "true", observeAdjacentRange(0, MAX, MIN, MIN));
		checkObservation(expect, "adjacent_range_max_neighbor",
				"true", observeAdjacentRange(MAX-3, MAX-1, MAX, MAX));
		checkObservation(expect, "adjacent_range_min_wrap",
				fixed ? "false" : "true", observeAdjacentRange(MIN, 0, MAX, MAX));
		checkObservation(expect, "adjacent_range_min_neighbor",
				"true", observeAdjacentRange(MIN+1, MIN+3, MIN, MIN));
	}

	/** Calls Range.touch(a, b, p) and returns its result or thrown type. */
	private static String observeTouchPoint(final int a, final int b, final int p){
		try{
			return Boolean.toString(Range.touch(a, b, p));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Calls Range.touch(a, b, a2, b2) and returns its result or thrown type. */
	private static String observeTouchRange(final int a, final int b,
			final int a2, final int b2){
		try{
			return Boolean.toString(Range.touch(a, b, a2, b2));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Calls Range.adjacent(a, b, p) and returns its result or thrown type. */
	private static String observeAdjacentPoint(final int a, final int b, final int p){
		try{
			return Boolean.toString(Range.adjacent(a, b, p));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Calls Range.adjacent(a, b, a2, b2) and returns its result or thrown type. */
	private static String observeAdjacentRange(final int a, final int b,
			final int a2, final int b2){
		try{
			return Boolean.toString(Range.adjacent(a, b, a2, b2));
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
