package test.structures;

import structures.CoverageArray2A;
import structures.IntList;

/** Finite oracle for CoverageArray.incrementRanges amount handling.
 * @author Sucrose
 */
public final class CoverageArrayRangeAmountCase{
	
	/** Runs baseline or fixed expectations for half-open range increments. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		
		positiveAmount3Range(expect);
		overlappingAmount3Ranges(expect);
		amount1Control(expect);
		emptyListControl(expect);
	}
	
	/** Exposes a positive amount above one being ignored. */
	private static void positiveAmount3Range(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 8);
		ca.incrementRanges(ranges(2, 5), 3);
		final int amount="baseline".equals(expect) ? 1 : 3;
		checkArray("positive-amount3-range", ca,
				0, 0, amount, amount, amount, 0, 0, 0);
		pass("positive-amount3-range", ca.toString(), expect);
	}
	
	/** Exposes two overlapping ranges both receiving unit increments. */
	private static void overlappingAmount3Ranges(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 8);
		ca.incrementRanges(ranges(1, 4, 3, 6), 3);
		final int amount="baseline".equals(expect) ? 1 : 3;
		checkArray("overlapping-amount3-ranges", ca,
				0, amount, amount, amount*2, amount, amount, 0, 0);
		pass("overlapping-amount3-ranges", ca.toString(), expect);
	}
	
	/** Preserves the existing amount-one behavior. */
	private static void amount1Control(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 8);
		ca.incrementRanges(ranges(0, 2, 6, 8), 1);
		checkArray("amount1-control", ca, 1, 1, 0, 0, 0, 0, 1, 1);
		pass("amount1-control", ca.toString(), expect);
	}
	
	/** Preserves empty-list no-op behavior. */
	private static void emptyListControl(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 8);
		ca.set(3, 2);
		ca.incrementRanges(new IntList(1), 3);
		checkArray("empty-list-control", ca, 0, 0, 0, 2, 0, 0, 0, 0);
		pass("empty-list-control", ca.toString(), expect);
	}
	
	/** Creates a flat list of half-open start, stop pairs. */
	private static IntList ranges(final int... bounds){
		final IntList ranges=new IntList(bounds.length);
		for(int i=0; i<bounds.length; i++){
			ranges.add(bounds[i]);
		}
		return ranges;
	}
	
	/** Checks exact coverage values. */
	private static void checkArray(final String name, final CoverageArray2A ca,
			final int... expected){
		for(int i=0; i<expected.length; i++){
			checkEquals(name+" element "+i, expected[i], ca.get(i));
		}
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
