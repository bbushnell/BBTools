package test.structures;

import structures.IntList;
import structures.LongList;

/** Baseline and fixed oracles for STU-005 through STU-008 IntList/LongList cases.
 * @author Sucrose
 */
public final class IntLongListCase{
	
	/** Runs all IntList and LongList edge cases. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		
		longListRangeAdd(expect);
		intListMode(expect);
		longListMode(expect);
		intListUnique(expect);
		intListMax(expect);
	}
	
	/** Verifies STU-005: range append grows capacity for the entire requested range. */
	private static void longListRangeAdd(final String expect){
		final LongList control=new LongList(4);
		control.add(new long[] {9, 8, 7, 6}, 1, 3);
		checkLongList("range-add-capacity-control", control, 8, 7);
		pass("range-add-capacity-control", control.toString(), expect);
		
		final LongList observed=new LongList(2);
		observed.add(7);
		final Throwable thrown=addRange(observed, new long[] {1, 2, 3, 4, 5}, 0, 5);
		if("baseline".equals(expect)){
			check(thrown instanceof ArrayIndexOutOfBoundsException,
					"baseline large range add should throw, got "+thrown);
			pass("range-add-capacity", thrown.getClass().getSimpleName(), expect);
		}else{
			check(thrown==null, "fixed large range add should not throw "+thrown);
			checkLongList("range-add-capacity", observed, 7, 1, 2, 3, 4, 5);
			pass("range-add-capacity", observed.toString(), expect);
		}
	}
	
	/** Verifies STU-006: IntList mode counts the first run once. */
	private static void intListMode(final String expect){
		checkEquals("int-mode-first-run-control", 1, IntList.mode(new int[] {1, 1, 2}, 3));
		checkEquals("int-mode-first-tie-control", 1, IntList.mode(new int[] {1, 2}, 2));
		
		final int expected="baseline".equals(expect) ? 1 : 2;
		checkEquals("int-mode-later-run", expected, IntList.mode(new int[] {1, 2, 2}, 3));
		checkEquals("int-mode-later-longer-run", expected,
				IntList.mode(new int[] {1, 1, 2, 2, 2}, 5));
		
		final IntList unsorted=new IntList(3);
		unsorted.add(2);
		unsorted.add(1);
		unsorted.add(2);
		checkEquals("int-mode-unsorted-later-run", expected, unsorted.modeUnsorted());
		pass("int-mode", Integer.toString(expected), expect);
	}
	
	/** Verifies STU-006: LongList mode counts the first run once. */
	private static void longListMode(final String expect){
		final LongList firstRun=new LongList(3);
		firstRun.add(1);
		firstRun.add(1);
		firstRun.add(2);
		checkEquals("long-mode-first-run-control", 1, firstRun.mode());
		
		final LongList firstTie=new LongList(2);
		firstTie.add(1);
		firstTie.add(2);
		checkEquals("long-mode-first-tie-control", 1, firstTie.mode());
		
		final LongList laterRun=new LongList(3);
		laterRun.add(1);
		laterRun.add(2);
		laterRun.add(2);
		checkEquals("long-mode-later-run", "baseline".equals(expect) ? 1 : 2,
				laterRun.mode());
		pass("long-mode", Long.toString(laterRun.mode()), expect);
	}
	
	/** Verifies STU-007: duplicate detection and a distinct-list control. */
	private static void intListUnique(final String expect){
		final IntList distinct=new IntList(2);
		distinct.add(5);
		distinct.add(6);
		checkEquals("unique-distinct-control", true, distinct.unique());
		
		final IntList singleton=new IntList(1);
		singleton.add(5);
		checkEquals("unique-singleton-control", true, singleton.unique());
		
		final IntList duplicate=new IntList(2);
		duplicate.add(5);
		duplicate.add(5);
		checkEquals("unique-duplicate", "baseline".equals(expect), duplicate.unique());
		pass("unique", Boolean.toString(duplicate.unique()), expect);
	}
	
	/** Verifies STU-008: max can return Integer.MIN_VALUE. */
	private static void intListMax(final String expect){
		final IntList control=new IntList(3);
		control.add(-9);
		control.add(3);
		control.add(1);
		checkEquals("max-control", 3, control.max());
		
		final IntList min=new IntList(1);
		min.add(Integer.MIN_VALUE);
		final int expected="baseline".equals(expect) ? -Integer.MAX_VALUE : Integer.MIN_VALUE;
		checkEquals("max-min-value", expected, min.max());
		pass("max", Integer.toString(expected), expect);
	}
	
	/** Calls range add and captures any thrown exception. */
	private static Throwable addRange(final LongList list, final long[] values,
			final int from, final int to){
		try{
			list.add(values, from, to);
			return null;
		}catch(Throwable t){return t;}
	}
	
	/** Checks exact list contents. */
	private static void checkLongList(final String name, final LongList observed,
			final long... expected){
		checkEquals(name+" size", expected.length, observed.size);
		for(int i=0; i<expected.length; i++){
			checkEquals(name+" element "+i, expected[i], observed.get(i));
		}
	}
	
	/** Checks an int oracle. */
	private static void checkEquals(final String name, final int expected, final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
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
