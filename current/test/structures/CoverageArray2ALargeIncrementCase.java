package test.structures;

import structures.CoverageArray2A;

/** Finite STU-015 oracles for CoverageArray2A large-increment saturation.
 * @author Sucrose
 */
public final class CoverageArray2ALargeIncrementCase{
	
	/** Runs the large-increment overflow checks. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		
		ordinaryLargeIncrement(expect);
		lowerLargeIncrement(expect);
		upperLargeIncrement(expect);
	}
	
	/** Verifies a huge increment from zero saturates without overflowing int. */
	private static void ordinaryLargeIncrement(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 8);
		ca.increment(1, Integer.MAX_VALUE);
		ca.increment(2, Integer.MAX_VALUE);
		checkEquals("ordinary-upper-saturation", 0xFFFF, ca.get(1));
		checkEquals("ordinary-lower-saturation", 0xFFFF, ca.get(2));
		pass("ordinary-large-increment", ca.get(1)+","+ca.get(2), expect);
	}
	
	/** Verifies lower-lane overflow does not disturb the packed upper lane. */
	private static void lowerLargeIncrement(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 8);
		ca.set(1, 3);
		ca.set(2, 1);
		ca.increment(2, Integer.MAX_VALUE);
		
		final int expectedUpper="baseline".equals(expect) ? 32771 : 3;
		final int expectedLower="baseline".equals(expect) ? 0 : 0xFFFF;
		checkEquals("lower-large-preserves-upper", expectedUpper, ca.get(1));
		checkEquals("lower-large-saturates-lower", expectedLower, ca.get(2));
		pass("lower-large-increment", ca.get(1)+","+ca.get(2), expect);
	}
	
	/** Verifies upper-lane overflow does not disturb the packed lower lane. */
	private static void upperLargeIncrement(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 8);
		ca.set(1, 1);
		ca.set(2, 4);
		ca.increment(1, Integer.MAX_VALUE);
		
		final int expectedUpper="baseline".equals(expect) ? 0 : 0xFFFF;
		checkEquals("upper-large-saturates-upper", expectedUpper, ca.get(1));
		checkEquals("upper-large-preserves-lower", 4, ca.get(2));
		pass("upper-large-increment", ca.get(1)+","+ca.get(2), expect);
	}
	
	/** Checks an int oracle. */
	private static void checkEquals(final String name, final int expected,
			final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}
	
	/** Prints one successful case row. */
	private static void pass(final String name, final String observed,
			final String expect){
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+observed);
	}
	
	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	
}
