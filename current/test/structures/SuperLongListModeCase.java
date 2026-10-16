package test.structures;

import structures.SuperLongList;

/** Finite oracle for SuperLongList overflow mode calculation.
 * @author Sucrose
 */
public final class SuperLongListModeCase{
	
	/** Runs baseline or fixed expectations for overflow-list mode cases. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		final boolean fixed="fixed".equals(expect);
		
		checkCase("single_overflow", fixed ? 5 : 0, 5);
		checkCase("all_unique_overflow", fixed ? 5 : 0, 5, 6, 7);
		checkCase("duplicate_overflow", 5, 5, 5, 6);
		checkCase("histogram_winner", 1, 1, 1, 5);
		checkCase("histogram_overflow_tie", 1, 1, 1, 5, 5);
	}
	
	/** Checks the mode of the supplied values in a histogram with limit 2. */
	private static void checkCase(final String name, final long expected,
			final long... values){
		final SuperLongList list=new SuperLongList(2);
		for(final long value : values){
			list.add(value);
		}
		list.sort();
		final long observed=list.mode();
		check(observed==expected,
				name+" expected mode "+expected+", got "+observed);
		System.out.println("PASS\t"+name+'\t'+observed);
	}
	
	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	
}
