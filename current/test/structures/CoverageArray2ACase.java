package test.structures;

import structures.CoverageArray2A;

/** Finite STU-010 oracles for packed CoverageArray2A compare-and-set repairs.
 * @author Sucrose
 */
public final class CoverageArray2ACase{
	
	/** Runs the finite packed-lane and range checks. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		
		lowerUpperIncrement(expect);
		zeroIncrement(expect);
		incrementSaturation(expect);
		lanePreservingSet(expect);
		laneClearingSet(expect);
		rangeIncrement(expect);
	}
	
	/** Verifies zero increments are no-ops on lower and upper lanes. */
	private static void zeroIncrement(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 8);
		ca.set(1, 12);
		ca.set(2, 34);
		ca.increment(1, 0);
		ca.increment(2, 0);
		checkEquals("zero-increment-upper", 12, ca.get(1));
		checkEquals("zero-increment-lower", 34, ca.get(2));
		pass("zero-increment", ca.get(1)+","+ca.get(2), expect);
	}
	
	/** Verifies lower and upper increments sharing one packed integer. */
	private static void lowerUpperIncrement(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 8);
		ca.increment(1, 7);
		ca.increment(2, 5);
		checkEquals("upper-increment", 7, ca.get(1));
		checkEquals("lower-increment", 5, ca.get(2));
		
		ca.increment(1);
		ca.increment(2, 2);
		checkEquals("upper-increment-one", 8, ca.get(1));
		checkEquals("lower-increment-two", 7, ca.get(2));
		pass("increment-lanes", ca.get(1)+","+ca.get(2), expect);
	}
	
	/** Verifies increment and set saturation at 16 bits. */
	private static void incrementSaturation(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 8);
		ca.set(1, 0xFFFE);
		ca.increment(1, 5);
		ca.set(2, 70000);
		checkEquals("increment-upper-saturation", 0xFFFF, ca.get(1));
		checkEquals("set-lower-saturation", 0xFFFF, ca.get(2));
		pass("saturation", ca.get(1)+","+ca.get(2), expect);
	}
	
	/** Verifies set() replaces one lane without disturbing the other lane. */
	private static void lanePreservingSet(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 8);
		ca.set(1, 0x1357);
		ca.set(2, 0x2468);
		checkEquals("set-upper", 0x1357, ca.get(1));
		checkEquals("set-lower", 0x2468, ca.get(2));
		
		ca.set(1, 3);
		checkEquals("set-upper-replaced", 3, ca.get(1));
		checkEquals("set-upper-preserves-lower", 0x2468, ca.get(2));
		
		ca.set(2, 4);
		checkEquals("set-lower-replaced", 4, ca.get(2));
		checkEquals("set-lower-preserves-upper", 3, ca.get(1));
		pass("set-preserves-other-lane", ca.get(1)+","+ca.get(2), expect);
	}
	
	/** Verifies set() can clear either lane to zero. */
	private static void laneClearingSet(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 8);
		ca.set(1, 3);
		ca.set(2, 4);
		ca.set(1, 0);
		final int expectedUpper="baseline".equals(expect) ? 3 : 0;
		checkEquals("set-upper-zero", expectedUpper, ca.get(1));
		checkEquals("set-upper-zero-preserves-lower", 4, ca.get(2));
		
		ca.set(2, 0);
		final int expectedLower="baseline".equals(expect) ? 4 : 0;
		checkEquals("set-lower-zero", expectedLower, ca.get(2));
		checkEquals("set-lower-zero-preserves-upper", expectedUpper, ca.get(1));
		pass("set-zero", ca.get(1)+","+ca.get(2), expect);
	}
	
	/** Verifies range increments clamp bounds and update every covered lane. */
	private static void rangeIncrement(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 8);
		ca.incrementRange(-2, 99, 3);
		checkArray("range-clamped", ca, 3, 3, 3, 3, 3, 3, 3, 3);
		
		ca.incrementRange(1, 4, 2);
		checkArray("range-middle", ca, 3, 5, 5, 5, 5, 3, 3, 3);
		pass("range-increment", ca.toString(), expect);
	}
	
	/** Checks exact coverage values. */
	private static void checkArray(final String name, final CoverageArray2A ca,
			final int... expected){
		for(int i=0; i<expected.length; i++){
			checkEquals(name+" element "+i, expected[i], ca.get(i));
		}
	}
	
	/** Checks an int oracle. */
	private static void checkEquals(final String name, final int expected, final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}
	
	/** Prints one successful case row. */
	private static void pass(final String name, final String observed, final String expect){
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+observed);
	}
	
	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
