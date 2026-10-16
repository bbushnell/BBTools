package test.structures;

import structures.LongHeapMap;

/** Baseline and fixed oracles for LongHeapMap retained-count overflow.
 * @author Jean
 */
public final class LongHeapCountOverflowCase{

	/** Runs finite LongHeapMap count-overflow checks. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");

		ordinaryRetainedIncrements(expect);
		positiveAddControl(expect);
		fullMapTooSmallNewKey(expect);
		incrementOverflowSpareCapacity(expect);
		incrementOverflowAboveFloor(expect);
		incrementOverflowMinimum(expect);
		addOverflowBoundary(expect);
	}

	/** Verifies ordinary retained keys accumulate positive counts. */
	private static void ordinaryRetainedIncrements(final String expect){
		final LongHeapMap map=new LongHeapMap(2);
		checkEquals("ordinary-first", 2, map.increment(17, 2));
		checkEquals("ordinary-second", 5, map.increment(17, 3));
		checkEquals("ordinary-count", 5, map.map.get(17));
		checkEquals("ordinary-size", 1, map.size());
		pass("ordinary-retained-increments", Integer.toString(map.map.get(17)), expect);
	}

	/** Verifies add reports true while a retained key's increment stays positive. */
	private static void positiveAddControl(final String expect){
		final LongHeapMap map=new LongHeapMap(1);
		checkEquals("positive-seed", Integer.MAX_VALUE-2,
				map.increment(23, Integer.MAX_VALUE-2));
		checkEquals("positive-add-retained", true, map.add(23));
		checkEquals("positive-count", Integer.MAX_VALUE-1, map.map.get(23));
		pass("positive-add-control", Integer.toString(map.map.get(23)), expect);
	}

	/** Verifies a full map normally returns zero for a too-small new key. */
	private static void fullMapTooSmallNewKey(final String expect){
		final LongHeapMap map=new LongHeapMap(1);
		checkEquals("full-seed", 7, map.increment(100, 7));
		checkEquals("too-small-rejection", 0, map.increment(99, 1));
		checkEquals("too-small-absent", false, map.contains(99));
		checkEquals("retained-count", 7, map.map.get(100));
		checkEquals("retained-size", 1, map.size());
		pass("full-map-too-small-new-key", Integer.toString(map.map.get(100)), expect);
	}

	/** Exposes retained-count overflow while the heap has spare capacity. */
	private static void incrementOverflowSpareCapacity(final String expect){
		final LongHeapMap map=new LongHeapMap(2);
		checkEquals("increment-spare-overflow-seed", Integer.MAX_VALUE,
				map.increment(11, Integer.MAX_VALUE));
		checkIncrementOverflow("increment-spare-overflow", map, 11, 1, expect);
	}

	/** Exposes retained-count overflow for a full-heap key above the floor. */
	private static void incrementOverflowAboveFloor(final String expect){
		final LongHeapMap map=new LongHeapMap(2);
		checkEquals("increment-above-floor-seed-minimum", 1, map.increment(7, 1));
		checkEquals("increment-above-floor-seed-retained", Integer.MAX_VALUE,
				map.increment(11, Integer.MAX_VALUE));
		checkIncrementOverflow("increment-above-floor-overflow", map, 11, 2, expect);
	}

	/** Exposes retained-count overflow for a full heap's minimum key. */
	private static void incrementOverflowMinimum(final String expect){
		final LongHeapMap map=new LongHeapMap(1);
		checkEquals("increment-minimum-overflow-seed", Integer.MAX_VALUE,
				map.increment(11, Integer.MAX_VALUE));
		checkIncrementOverflow("increment-minimum-overflow", map, 11, 1, expect);
	}

	/** Checks one retained-key increment overflow branch. */
	private static void checkIncrementOverflow(final String name, final LongHeapMap map,
			final long key, final int expectedSize, final String expect){
		AssertionError failure=null;
		int observed=0;
		try{
			observed=map.increment(key, 1);
		}catch(AssertionError e){failure=e;}
		if("baseline".equals(expect)){
			check(failure==null, "baseline "+name+" should silently wrap, got "+failure);
			checkEquals(name+"-return", Integer.MIN_VALUE, observed);
			checkEquals(name+"-count", Integer.MIN_VALUE, map.map.get(key));
			checkEquals(name+"-retained", true, map.contains(key));
			checkEquals(name+"-size", expectedSize, map.size());
			pass(name, Integer.toString(observed), expect);
		}else{
			checkAssertionsEnabled();
			check(failure!=null, "fixed "+name+" should assert on unrepresentable count");
			checkEquals(name+"-count", Integer.MAX_VALUE, map.map.get(key));
			checkEquals(name+"-retained", true, map.contains(key));
			checkEquals(name+"-size", expectedSize, map.size());
			pass(name, "asserted", expect);
		}
	}

	/** Exposes add returning false for an existing key after count overflow. */
	private static void addOverflowBoundary(final String expect){
		final LongHeapMap map=new LongHeapMap(1);
		checkEquals("add-overflow-seed", Integer.MAX_VALUE,
				map.increment(11, Integer.MAX_VALUE));
		AssertionError failure=null;
		boolean retained=false;
		try{
			retained=map.add(11);
		}catch(AssertionError e){failure=e;}
		if("baseline".equals(expect)){
			check(failure==null, "baseline add should silently wrap, got "+failure);
			checkEquals("add-overflow-return", false, retained);
			checkEquals("add-overflow-count", Integer.MIN_VALUE, map.map.get(11));
			checkEquals("add-overflow-retained", true, map.contains(11));
			checkEquals("add-overflow-size", 1, map.size());
			pass("add-overflow-boundary", Boolean.toString(retained), expect);
		}else{
			checkAssertionsEnabled();
			check(failure!=null, "fixed add should assert on unrepresentable count");
			checkEquals("add-overflow-count", Integer.MAX_VALUE, map.map.get(11));
			checkEquals("add-overflow-retained", true, map.contains(11));
			checkEquals("add-overflow-size", 1, map.size());
			pass("add-overflow-boundary", "asserted", expect);
		}
	}

	/** Checks an int oracle. */
	private static void checkEquals(final String name, final int expected, final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Checks a boolean oracle. */
	private static void checkEquals(final String name, final boolean expected,
			final boolean observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Checks that JVM assertions are enabled. */
	private static void checkAssertionsEnabled(){
		check(LongHeapCountOverflowCase.class.desiredAssertionStatus(),
				"Fixed count-overflow checks require -ea");
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
