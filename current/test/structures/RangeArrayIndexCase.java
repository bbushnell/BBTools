package test.structures;

import structures.Range;

/** Literal ordinary-array oracle for Range.findIndex.
 * @author Sucrose
 */
public final class RangeArrayIndexCase{

	/** Runs unchanged expectations for finite Range.findIndex array cases. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=unchanged");
		check("expect=unchanged".equals(args[0]), "Unknown expectation: "+args[0]);
		final String expect="unchanged";

		checkObservation(expect, "empty", "-1", observe(5));
		checkObservation(expect, "before_first", "-1", observe(5, ranges()));
		checkObservation(expect, "first_contained", "0", observe(11, ranges()));
		checkObservation(expect, "gap", "-2", observe(15, ranges()));
		checkObservation(expect, "second_contained", "1", observe(21, ranges()));
		checkObservation(expect, "after_last", "-3", observe(30, ranges()));
	}

	/** Calls Range.findIndex and returns its result or thrown type. */
	private static String observe(final int p, final Range... ranges){
		try{
			return Integer.toString(Range.findIndex(p, ranges));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Returns a tiny sorted, disjoint Range array. */
	private static Range[] ranges(){
		return new Range[] {new Range(10, 12), new Range(20, 22)};
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
