package test.structures;

import structures.DoubleList;
import structures.FloatList;

/** Baseline and fixed oracles for DoubleList/FloatList finite alias cases.
 * @author Sucrose
 */
public final class DoubleFloatListAliasCase{

	/** Runs DoubleList and FloatList finite alias cases. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");

		doubleUniqueCounts(expect);
		floatUniqueCounts(expect);
		doubleEmptySelfAppend(expect);
		doubleNonemptySelfAppend(expect);
	}

	/** Verifies DoubleList getUniqueCounts handles distinct and aliased output. */
	private static void doubleUniqueCounts(final String expect){
		final DoubleList distinctInput=doubleCountsInput();
		final DoubleList distinctCounts=new DoubleList(5);
		distinctCounts.add(9);
		distinctCounts.add(7);
		distinctCounts.add(8);
		distinctCounts.add(6);
		distinctInput.getUniqueCounts(distinctCounts);
		checkEquals("double-unique-distinct-input-size", 3, distinctInput.size());
		checkEquals("double-unique-distinct-count-size", 3, distinctCounts.size());
		checkEquals("double-unique-distinct-value-0", 1, distinctInput.get(0));
		checkEquals("double-unique-distinct-value-1", 2, distinctInput.get(1));
		checkEquals("double-unique-distinct-value-2", 3, distinctInput.get(2));
		checkEquals("double-unique-distinct-count-0", 2, distinctCounts.get(0));
		checkEquals("double-unique-distinct-count-1", 1, distinctCounts.get(1));
		checkEquals("double-unique-distinct-count-2", 3, distinctCounts.get(2));

		final DoubleList aliased=doubleCountsInput();
		if("baseline".equals(expect)){
			aliased.getUniqueCounts(aliased);
			checkEquals("double-unique-alias-cleared", 0, aliased.size());
			check("[]".equals(aliased.toStringListView()),
					"Expected cleared aliased DoubleList, got "+aliased);
			pass("double-unique-alias", "cleared", expect);
		}else{
			checkAssertionsEnabled();
			boolean thrown=false;
			try{
				aliased.getUniqueCounts(aliased);
			}catch(final AssertionError expected){
				thrown=true;
			}
			check(thrown, "Expected double getUniqueCounts alias assertion");
			checkDoubleCountsInput("double-unique-alias-preserved", aliased);
			pass("double-unique-alias", "asserted", expect);
		}
		pass("double-unique-distinct", distinctCounts.toStringListView(), expect);
	}

	/** Verifies FloatList getUniqueCounts handles distinct and aliased output. */
	private static void floatUniqueCounts(final String expect){
		final FloatList distinctInput=floatCountsInput();
		final FloatList distinctCounts=new FloatList(5);
		distinctCounts.add(9);
		distinctCounts.add(7);
		distinctCounts.add(8);
		distinctCounts.add(6);
		distinctInput.getUniqueCounts(distinctCounts);
		checkEquals("float-unique-distinct-input-size", 3, distinctInput.size());
		checkEquals("float-unique-distinct-count-size", 3, distinctCounts.size());
		checkEquals("float-unique-distinct-value-0", 1, distinctInput.get(0));
		checkEquals("float-unique-distinct-value-1", 2, distinctInput.get(1));
		checkEquals("float-unique-distinct-value-2", 3, distinctInput.get(2));
		checkEquals("float-unique-distinct-count-0", 2, distinctCounts.get(0));
		checkEquals("float-unique-distinct-count-1", 1, distinctCounts.get(1));
		checkEquals("float-unique-distinct-count-2", 3, distinctCounts.get(2));

		final FloatList aliased=floatCountsInput();
		if("baseline".equals(expect)){
			aliased.getUniqueCounts(aliased);
			checkEquals("float-unique-alias-cleared", 0, aliased.size());
			check("[]".equals(aliased.toStringListView()),
					"Expected cleared aliased FloatList, got "+aliased);
			pass("float-unique-alias", "cleared", expect);
		}else{
			checkAssertionsEnabled();
			boolean thrown=false;
			try{
				aliased.getUniqueCounts(aliased);
			}catch(final AssertionError expected){
				thrown=true;
			}
			check(thrown, "Expected float getUniqueCounts alias assertion");
			checkFloatCountsInput("float-unique-alias-preserved", aliased);
			pass("float-unique-alias", "asserted", expect);
		}
		pass("float-unique-distinct", distinctCounts.toStringListView(), expect);
	}

	/** Verifies empty DoubleList self-append remains a no-op. */
	private static void doubleEmptySelfAppend(final String expect){
		final DoubleList list=new DoubleList(1);
		list.append(list);
		checkEquals("double-empty-self-append-size", 0, list.size());
		pass("double-empty-self-append", Integer.toString(list.size()), expect);
	}

	/** Verifies nonempty DoubleList self-append asserts after repair. */
	private static void doubleNonemptySelfAppend(final String expect){
		if("baseline".equals(expect)){return;}
		checkAssertionsEnabled();

		final DoubleList list=new DoubleList(1);
		list.add(5);
		boolean thrown=false;
		try{
			list.append(list);
		}catch(final AssertionError expected){
			thrown=true;
		}
		check(thrown, "Expected nonempty DoubleList self-append assertion");
		checkEquals("double-nonempty-self-append-preserved-size", 1, list.size());
		checkEquals("double-nonempty-self-append-preserved-value", 5, list.get(0));
		pass("double-nonempty-self-append", "asserted", expect);
	}

	/** Returns a sorted DoubleList with three unique values. */
	private static DoubleList doubleCountsInput(){
		final DoubleList list=new DoubleList(6);
		list.add(1);
		list.add(1);
		list.add(2);
		list.add(3);
		list.add(3);
		list.add(3);
		return list;
	}

	/** Returns a sorted FloatList with three unique values. */
	private static FloatList floatCountsInput(){
		final FloatList list=new FloatList(6);
		list.add(1);
		list.add(1);
		list.add(2);
		list.add(3);
		list.add(3);
		list.add(3);
		return list;
	}

	/** Checks the original DoubleList alias input. */
	private static void checkDoubleCountsInput(final String name, final DoubleList list){
		checkEquals(name+"-size", 6, list.size());
		checkEquals(name+"-0", 1, list.get(0));
		checkEquals(name+"-1", 1, list.get(1));
		checkEquals(name+"-2", 2, list.get(2));
		checkEquals(name+"-3", 3, list.get(3));
		checkEquals(name+"-4", 3, list.get(4));
		checkEquals(name+"-5", 3, list.get(5));
	}

	/** Checks the original FloatList alias input. */
	private static void checkFloatCountsInput(final String name, final FloatList list){
		checkEquals(name+"-size", 6, list.size());
		checkEquals(name+"-0", 1, list.get(0));
		checkEquals(name+"-1", 1, list.get(1));
		checkEquals(name+"-2", 2, list.get(2));
		checkEquals(name+"-3", 3, list.get(3));
		checkEquals(name+"-4", 3, list.get(4));
		checkEquals(name+"-5", 3, list.get(5));
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

	/** Checks an int oracle. */
	private static void checkEquals(final String name, final int expected,
			final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Checks that JVM assertions are enabled. */
	private static void checkAssertionsEnabled(){
		check(DoubleFloatListAliasCase.class.desiredAssertionStatus(),
				"Fixed alias checks require -ea");
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
