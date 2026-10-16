package test.structures;

import java.util.ArrayList;

import structures.Feature;
import structures.Range;

/** Finite oracle for Range.findIndexAfterBinary insertion indexes.
 * @author Jean
 */
public final class RangeAfterIndexCase{

	/** Runs baseline or fixed expectations for findIndexAfterBinary. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		final boolean fixed="fixed".equals(expect);

		checkObservation(expect, "two_before_first_zero", fixed ? "0" : "AssertionError",
				observe(0, 10, 20));
		checkObservation(expect, "two_before_first_nine", fixed ? "0" : "AssertionError",
				observe(9, 10, 20));
		checkObservation(expect, "two_between", "1", observe(15, 10, 20));
		checkObservation(expect, "two_equal_first", "1", observe(10, 10, 20));
		checkObservation(expect, "two_equal_last", "2", observe(20, 10, 20));
		checkObservation(expect, "two_after_last", "2", observe(25, 10, 20));
		checkObservation(expect, "duplicate_equal_first", "2", observe(10, 10, 10, 20));
		checkObservation(expect, "duplicate_equal_second", "3", observe(20, 10, 20, 20, 30));
		checkObservation(expect, "single_before", "0", observe(0, 10));
		checkObservation(expect, "single_equal", "1", observe(10, 10));
		checkObservation(expect, "single_after", "1", observe(20, 10));
		checkObservation(expect, "empty", "-1", observe(0));
	}

	/** Calls Range.findIndexAfterBinary and returns its result or thrown type. */
	private static String observe(final int p, final int... starts){
		try{
			final ArrayList<LiteralFeature> list=new ArrayList<LiteralFeature>(starts.length);
			for(final int start : starts){
				list.add(new LiteralFeature(start));
			}
			return Integer.toString(Range.findIndexAfterBinary(p, list));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Checks and prints a value-or-exception outcome. */
	private static void checkObservation(final String expect, final String name,
			final String expected, final String observed){
		check(expected.equals(observed), name+" expected "+expected+", got "+observed);
		System.out.println("PASS	"+expect+'	'+name+'	'+observed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Literal sorted feature with a start coordinate sufficient for binary search. */
	private static final class LiteralFeature implements Feature{

		/** Constructs a feature whose start and stop are the same point. */
		private LiteralFeature(final int start_){
			start=start_;
		}

		@Override
		public int start(){return start;}

		@Override
		public int stop(){return start;}

		@Override
		public int strand(){return 0;}

		@Override
		public String seqid(){return "literal";}

		@Override
		public String name(){return Integer.toString(start);}

		private final int start;
	}
}
