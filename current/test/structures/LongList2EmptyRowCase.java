package test.structures;

import java.util.Arrays;

import structures.LongList2;

/** Finite oracle for LongList2 zero-length entry rows and ordinary insert paths.
 * @author Sucrose
 */
public final class LongList2EmptyRowCase{

	/** Runs baseline or fixed expectations for zero-length LongList2 rows. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		final boolean fixed="fixed".equals(expect);

		checkObservation(expect, "zero_add", fixed ? "1:[5, -1]:-1" :
				"ArrayIndexOutOfBoundsException", observeZeroAdd());
		checkObservation(expect, "zero_set", fixed ? "1:[6, -1]:-1" :
				"ArrayIndexOutOfBoundsException", observeZeroSet());
		checkObservation(expect, "null_row", "1:[5, -1]:-1", observeNullRow());
		checkObservation(expect, "spare_sentinel", "1:[3, 5]:-1", observeSpareSentinel());
		checkObservation(expect, "full_row_growth",
				"1:[3, 5000000000, 7000000000, -1]:-1",
				observeFullRowGrowth());
		checkObservation(expect, "duplicate", "0:[3, -1]:-1", observeDuplicate());
		checkObservation(expect, "ignored_len_getlen", "-1:0:[7, -1]",
				observeIgnoredLenGetLen());
	}

	/** Observes add(new long[0]) followed by insertIntoList(). */
	private static String observeZeroAdd(){
		try{
			final LongList2 list=new LongList2(1);
			list.add(new long[0]);
			final int added=list.insertIntoList(5, 0);
			return observation(added, list, 0);
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Observes set(0, new long[0]) followed by insertIntoList(). */
	private static String observeZeroSet(){
		try{
			final LongList2 list=new LongList2(1);
			list.set(0, new long[0]);
			final int added=list.insertIntoList(6, 0);
			return observation(added, list, 0);
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Observes inserting into a null row appended by loc==size. */
	private static String observeNullRow(){
		final LongList2 list=new LongList2(1);
		final int added=list.insertIntoList(5, 0);
		return observation(added, list, 0);
	}

	/** Observes insertion before an existing INVALID sentinel. */
	private static String observeSpareSentinel(){
		final LongList2 list=new LongList2(1);
		list.set(0, new long[] {3, LongList2.INVALID});
		final int added=list.insertIntoList(5, 0);
		return observation(added, list, 0);
	}

	/** Observes doubling a full row and appending before the new INVALID spare. */
	private static String observeFullRowGrowth(){
		final LongList2 list=new LongList2(1);
		list.set(0, new long[] {3, 5000000000L});
		final int added=list.insertIntoList(7000000000L, 0);
		return observation(added, list, 0);
	}

	/** Observes duplicate detection before an existing INVALID sentinel. */
	private static String observeDuplicate(){
		final LongList2 list=new LongList2(1);
		list.set(0, new long[] {3, LongList2.INVALID});
		final int added=list.insertIntoList(3, 0);
		return observation(added, list, 0);
	}

	/** Observes add(entry, len) ignoring len and getLen() compatibility markers. */
	private static String observeIgnoredLenGetLen(){
		final LongList2 list=new LongList2(1);
		list.add(new long[] {7, LongList2.INVALID}, 99);
		return list.getLen(0)+":"+list.getLen(1)+":"+Arrays.toString(list.get(0));
	}

	/** Returns the inserted count, retained row, and LongList3-compatible length. */
	private static String observation(final int added, final LongList2 list, final int loc){
		return added+":"+Arrays.toString(list.get(loc))+":"+list.getLen(loc);
	}

	/** Checks and prints a value-or-exception outcome. */
	private static void checkObservation(final String expect, final String name,
			final String expected, final String observed){
		check(expected.equals(observed), name+" expected "+expected+", got "+observed);
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+observed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

}
