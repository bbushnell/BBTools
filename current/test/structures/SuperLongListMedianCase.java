package test.structures;

import structures.SuperLongList;

/** Finite oracle for SuperLongList median percentile calculation.
 * @author Sucrose
 */
public final class SuperLongListMedianCase{

	/** Runs baseline or fixed expectations for median cases. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		final boolean fixed="fixed".equals(expect);

		checkCase("empty", 0);
		checkCase("singleton_histogram", fixed ? 3 : 0, 3);
		checkCase("singleton_overflow", fixed ? 7 : 0, 7);
		checkCase("odd_histogram", fixed ? 2 : 1, 1, 2, 3);
		checkCase("odd_overflow", fixed ? 6 : 5, 5, 6, 7);
		checkCase("odd_split", fixed ? 5 : 2, 2, 5, 6);
		checkCase("even_lower_histogram", 1, 1, 2);
		checkCase("even_lower_overflow", 5, 5, 6);
		checkCase("even_lower_split", 2, 2, 5);
		checkCase("duplicate_middle", 5, 5, 5, 6);
		checkHugeHistogramCount();
		checkEmptyLimit0();
	}

	/** Checks the median of the supplied values in a histogram with limit 5. */
	private static void checkCase(final String name, final long expected,
			final long... values){
		final SuperLongList list=new SuperLongList(5);
		for(final long value : values){
			list.add(value);
		}
		list.sort();
		final long observed=list.median();
		check(observed==expected,
				name+" expected median "+expected+", got "+observed);
		System.out.println("PASS\t"+name+'\t'+observed);
	}

	/** Checks overflow-safe median rank at the largest finite histogram count. */
	private static void checkHugeHistogramCount(){
		final SuperLongList list=new SuperLongList(5);
		list.increment(1, Long.MAX_VALUE);
		list.sort();
		final long observed=list.median();
		check(observed==1,
				"huge_histogram_count expected median 1, got "+observed);
		System.out.println("PASS\thuge_histogram_count\t"+observed);
	}

	/** Checks that an empty zero-limit histogram still trips its final assertion. */
	private static void checkEmptyLimit0(){
		final SuperLongList list=new SuperLongList(0);
		try{
			list.median();
		}catch(final AssertionError expected){
			System.out.println("PASS\tempty_limit0\tAssertionError");
			return;
		}
		throw new AssertionError("empty_limit0 expected AssertionError");
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

}
