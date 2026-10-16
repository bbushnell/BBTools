package test.structures;

import structures.Range;

/** Literal equals-contract checks for Range coordinate pairs.
 * @author Jean
 */
public final class RangeEqualsContractCase{

	/** Runs current and proposed Range.equals(Object) expectations. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final boolean fixed="expect=fixed".equals(args[0]);
		check(fixed || "expect=baseline".equals(args[0]), "Unknown expectation: "+args[0]);
		final String expect=(fixed ? "fixed" : "baseline");

		rangeControls(expect);
		rangeSubclassControls(expect);
		rangeHashAndOrderControls(expect);
		rangeBoundaries(expect, fixed);
	}

	/** Checks literal Range equality controls through typed and Object calls. */
	private static void rangeControls(final String expect){
		final Range value=new Range(2, 5);
		checkBoolean(expect, "range-equal-typed", true, value.equals(new Range(2, 5)));
		checkBoolean(expect, "range-equal-object", true, value.equals((Object)new Range(2, 5)));
		checkBoolean(expect, "range-different-start", false, value.equals(new Range(3, 5)));
		checkBoolean(expect, "range-different-end", false, value.equals(new Range(2, 6)));
	}

	/** Checks that Range equality remains compatible with subclasses. */
	private static void rangeSubclassControls(final String expect){
		final Range base=new Range(2, 5);
		final RangeChild child=new RangeChild(2, 5);
		checkBoolean(expect, "range-subclass-argument-typed", true, base.equals(child));
		checkBoolean(expect, "range-subclass-argument-object", true, base.equals((Object)child));
		checkBoolean(expect, "range-subclass-receiver-typed", true, child.equals(base));
		checkBoolean(expect, "range-subclass-receiver-object", true, child.equals((Object)base));
		checkBoolean(expect, "range-subclass-different", false, child.equals(new Range(2, 6)));
	}

	/** Checks Range hash and ordering controls remain coordinate-based. */
	private static void rangeHashAndOrderControls(final String expect){
		checkBoolean(expect, "range-hash-equal", true,
				new Range(2, 5).hashCode()==new Range(2, 5).hashCode());
		checkInteger(expect, "range-compare-equal", 0,
				new Range(2, 5).compareTo(new Range(2, 5)));
		checkInteger(expect, "range-compare-start-low", -2,
				new Range(1, 9).compareTo(new Range(3, 4)));
		checkInteger(expect, "range-compare-start-high", 2,
				new Range(3, 4).compareTo(new Range(1, 9)));
		checkInteger(expect, "range-compare-end-low", -1,
				new Range(2, 5).compareTo(new Range(2, 6)));
		checkInteger(expect, "range-compare-end-high", 1,
				new Range(2, 6).compareTo(new Range(2, 5)));
	}

	/** Checks current Range throws versus proposed Object.equals-compatible false. */
	private static void rangeBoundaries(final String expect, final boolean fixed){
		checkObservation(expect, "range-object-null", fixed ? "false" : "NullPointerException",
				observeRangeObject(null));
		checkObservation(expect, "range-object-unrelated", fixed ? "false" : "ClassCastException",
				observeRangeObject(new Object()));
		checkObservation(expect, "range-typed-null", fixed ? "false" : "NullPointerException",
				observeRangeTypedNull());
	}

	/** Calls Range.equals(Object) and returns a value or thrown type. */
	private static String observeRangeObject(final Object other){
		try{
			return Boolean.toString(new Range(2, 5).equals(other));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Calls Range.equals(Range) with a typed null and returns a value or thrown type. */
	private static String observeRangeTypedNull(){
		try{
			return Boolean.toString(new Range(2, 5).equals((Range)null));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Checks and prints a boolean outcome. */
	private static void checkBoolean(final String expect, final String name,
			final boolean expected, final boolean observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
		pass(expect, name, Boolean.toString(observed));
	}

	/** Checks and prints an integer outcome. */
	private static void checkInteger(final String expect, final String name,
			final int expected, final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
		pass(expect, name, Integer.toString(observed));
	}

	/** Checks and prints a value-or-exception outcome. */
	private static void checkObservation(final String expect, final String name,
			final String expected, final String observed){
		check(expected.equals(observed), name+" expected "+expected+", got "+observed);
		pass(expect, name, observed);
	}

	/** Prints a checked outcome. */
	private static void pass(final String expect, final String name, final String observed){
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+observed);
	}

	/** Fails independently of assertion flags so an oracle cannot silently disappear. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Range subclass used to preserve instanceof-compatible equality. */
	private static final class RangeChild extends Range{

		/** Constructs a child range from inclusive literal coordinates. */
		private RangeChild(final int a_, final int b_){
			super(a_, b_);
		}
	}
}
