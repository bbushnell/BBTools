package test.structures;

import structures.IntLongHashMap;

/** Baseline and fixed oracles for IntLongHashMap thresholded array export.
 * @author Jean
 */
public final class IntLongHashMapThresholdCase{
	
	/** Runs the IntLongHashMap threshold checks. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		
		emptyThresholdZero(expect);
		emptyThresholdOne(expect);
		presentThresholdOne(expect);
	}
	
	/** Verifies that empty cells are not exported as sentinel keys at thresh=0. */
	private static void emptyThresholdZero(final String expect){
		final IntLongHashMap map=new IntLongHashMap(3);
		AssertionError failure=null;
		try{
			map.toArray(0);
		}catch(AssertionError e){failure=e;}
		if("baseline".equals(expect)){
			check(failure!=null && failure.getMessage()!=null
					&& failure.getMessage().startsWith("\nNegative key ("),
					"baseline empty threshold zero should hit the sentinel assertion, got "+failure);
			pass("empty-threshold-zero", "sentinel assertion", expect);
		}else{
			check(failure==null, "fixed empty threshold zero should not throw "+failure);
			checkEquals("empty-threshold-zero-length", 0, map.toArray(0).length);
			pass("empty-threshold-zero", "0", expect);
		}
	}
	
	/** Verifies the empty positive-threshold control. */
	private static void emptyThresholdOne(final String expect){
		final IntLongHashMap map=new IntLongHashMap(3);
		final long[] keys=map.toArray(1);
		checkEquals("empty-threshold-one-length", 0, keys.length);
		pass("empty-threshold-one", Integer.toString(keys.length), expect);
	}
	
	/** Verifies a stored positive value still exports at thresh=1. */
	private static void presentThresholdOne(final String expect){
		final IntLongHashMap map=new IntLongHashMap(3);
		checkEquals("put-present-key", true, map.put(5, 1));
		final long[] keys=map.toArray(1);
		checkEquals("present-threshold-one-length", 1, keys.length);
		checkEquals("present-threshold-one-key", 5, keys[0]);
		pass("present-threshold-one", Long.toString(keys[0]), expect);
	}
	
	/** Checks a long oracle. */
	private static void checkEquals(final String name, final long expected, final long observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}
	
	/** Checks a boolean oracle. */
	private static void checkEquals(final String name, final boolean expected,
			final boolean observed){
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
