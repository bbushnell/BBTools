package test.structures;

import structures.CRange;

/** Finite equality and nonnegative endpoint-boundary oracles for CRange.
 * @author Shinobu
 */
public final class CRangeEdgeCase{

	/** Runs the same literal cases before and after the bounded repairs. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final boolean baseline="expect=baseline".equals(args[0]);
		check(baseline || "expect=fixed".equals(args[0]), "Unknown expectation: "+args[0]);
		final String mode=baseline ? "baseline" : "fixed";
		equalityBoundary(mode, baseline, "null-object", null, false);
		equalityBoundary(mode, baseline, "null-range", null, true);
		equalityBoundary(mode, baseline, "other-type", "not a range", false);
		equalityControls(mode);
		touchBoundary(mode, baseline, "maximum-overlap", Integer.MAX_VALUE-1,
				Integer.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE);
		touchBoundary(mode, baseline, "maximum-left-adjacency", Integer.MAX_VALUE,
				Integer.MAX_VALUE, Integer.MAX_VALUE-1, Integer.MAX_VALUE-1);
		touchBoundary(mode, baseline, "maximum-span", 0, Integer.MAX_VALUE, 0, 0);
		touchControls(mode);
	}

	/** Checks null through both overloads and an unrelated object through Object.equals. */
	private static void equalityBoundary(final String mode, final boolean baseline,
			final String name, final Object other, final boolean typed){
		final CRange range=new CRange(7, 10, 20);
		RuntimeException failure=null;
		boolean equal=false;
		try{
			equal=typed ? range.equals((CRange)other) : range.equals(other);
		}catch(final RuntimeException e){failure=e;}
		if(baseline){
			check(other==null ? failure instanceof NullPointerException :
					failure instanceof ClassCastException, name+" expected baseline type/null failure, got "+failure);
			pass(mode, name, failure.getClass().getSimpleName());
		}else{
			check(failure==null && !equal, name+" expected false without exception, got "+equal+", "+failure);
			pass(mode, name, "false");
		}
	}

	/** Checks equal coordinates, differing fields and payload-independent identity. */
	private static void equalityControls(final String mode){
		final CRange a=new CRange(7, 10, 20, "first");
		final CRange b=new CRange(7, 10, 20, "second");
		check(a.equals(b) && a.equals((Object)b), "Equal coordinates should match through either overload");
		check(a.hashCode()==b.hashCode(), "Equal ranges must have equal hashes");
		check(!a.equals(new CRange(8, 10, 20)), "Contig must participate in equality");
		check(!a.equals(new CRange(7, 11, 20)), "Start must participate in equality");
		check(!a.equals(new CRange(7, 10, 21)), "End must participate in equality");
		pass(mode, "equality-controls", "passed");
	}

	/** All three literal boundary pairs overlap or are directly adjacent. */
	private static void touchBoundary(final String mode, final boolean baseline,
			final String name, final int a1, final int b1, final int a2, final int b2){
		check(a1>=0 && a1<=b1 && a2>=0 && a2<=b2, "Fixture requires valid nonnegative ranges");
		final boolean observed=CRange.touch(a1, b1, a2, b2);
		check(observed==!baseline, name+" expected "+!baseline+", got "+observed);
		pass(mode, name, Boolean.toString(observed));
	}

	/** Checks ordinary overlap, adjacency, separation and maximum-coordinate adjacency. */
	private static void touchControls(final String mode){
		check(CRange.touch(10, 20, 20, 25), "Shared endpoint should touch");
		check(CRange.touch(10, 20, 21, 25), "Consecutive endpoints should touch");
		check(!CRange.touch(10, 20, 22, 25), "A one-base gap must not touch");
		check(!CRange.touch(Integer.MAX_VALUE-2, Integer.MAX_VALUE-2,
				Integer.MAX_VALUE, Integer.MAX_VALUE), "Maximum-coordinate gap must remain separate");
		check(CRange.adjacent(Integer.MAX_VALUE-1, Integer.MAX_VALUE-1,
				Integer.MAX_VALUE, Integer.MAX_VALUE), "Maximum endpoint has a valid left neighbor");
		check(CRange.adjacent(Integer.MAX_VALUE, Integer.MAX_VALUE,
				Integer.MAX_VALUE-1, Integer.MAX_VALUE-1), "Adjacency should hold in the reverse direction");
		pass(mode, "touch-controls", "passed");
	}

	/** Prints a checked outcome. */
	private static void pass(final String mode, final String name, final String value){
		System.out.println("PASS\t"+mode+'\t'+name+'\t'+value);
	}

	/** Keeps fixture checks active regardless of JVM assertion settings. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
