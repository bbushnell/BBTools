package test.structures;

import structures.DoubleList;
import structures.FloatList;

/** Baseline and fixed oracles for DoubleList/FloatList finite mode cases.
 * @author Sucrose
 */
public final class DoubleFloatListCase{

	/** Runs DoubleList and FloatList finite mode cases. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");

		doubleListMode(expect);
		floatListMode(expect);
	}

	/** Verifies DoubleList mode counts the first run once. */
	private static void doubleListMode(final String expect){
		final DoubleList firstRun=new DoubleList(3);
		firstRun.add(1);
		firstRun.add(1);
		firstRun.add(2);
		checkEquals("double-mode-first-run-control", 1, firstRun.mode());

		final DoubleList firstTie=new DoubleList(2);
		firstTie.add(1);
		firstTie.add(2);
		checkEquals("double-mode-first-tie-control", 1, firstTie.mode());

		final DoubleList laterRun=new DoubleList(3);
		laterRun.add(1);
		laterRun.add(2);
		laterRun.add(2);
		final double expected="baseline".equals(expect) ? 1 : 2;
		checkEquals("double-mode-later-run", expected, laterRun.mode());
		pass("double-mode", Double.toString(laterRun.mode()), expect);
	}

	/** Verifies FloatList mode counts the first run once. */
	private static void floatListMode(final String expect){
		final FloatList firstRun=new FloatList(3);
		firstRun.add(1);
		firstRun.add(1);
		firstRun.add(2);
		checkEquals("float-mode-first-run-control", 1, firstRun.mode());

		final FloatList firstTie=new FloatList(2);
		firstTie.add(1);
		firstTie.add(2);
		checkEquals("float-mode-first-tie-control", 1, firstTie.mode());

		final FloatList laterRun=new FloatList(3);
		laterRun.add(1);
		laterRun.add(2);
		laterRun.add(2);
		final float expected="baseline".equals(expect) ? 1 : 2;
		checkEquals("float-mode-later-run", expected, laterRun.mode());
		pass("float-mode", Float.toString(laterRun.mode()), expect);
	}

	/** Checks a double oracle. */
	private static void checkEquals(final String name, final double expected,
			final double observed){
		check(Double.compare(expected, observed)==0,
				name+" expected "+expected+", got "+observed);
	}

	/** Checks a float oracle. */
	private static void checkEquals(final String name, final float expected,
			final float observed){
		check(Float.compare(expected, observed)==0,
				name+" expected "+expected+", got "+observed);
	}

	/** Prints one successful case row. */
	private static void pass(final String name, final String observed, final String expect){
		System.out.println("PASS	"+expect+'	'+name+'	'+observed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
