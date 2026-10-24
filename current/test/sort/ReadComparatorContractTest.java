package test.sort;

import java.nio.charset.StandardCharsets;
import java.util.Comparator;

import sort.ReadComparatorTopological;
import sort.ReadComparatorTopological5Bit;
import stream.Read;

/** Focused comparator contract tests for sort package maintenance.
 * @author Sucrose
 */
public final class ReadComparatorContractTest{

	/** Runs baseline or repaired comparator expectations. */
	public static void main(String[] args){
		check(args.length==1, "Expected expect=baseline|fixed");
		check(args[0].startsWith("expect="), "Expected expect=baseline|fixed, got "+args[0]);
		final String expect=args[0].substring("expect=".length());
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline|fixed, got "+args[0]);

		topologicalMixedMatePresence(expect);
		topological5BitMixedMatePresence(expect);
		topologicalStrictPrefixBeforeMate(expect);
		topological5BitStrictPrefixBeforeMate(expect);
		topologicalCompareMatesFalse(expect);
		topological5BitCompareMatesFalse(expect);
	}

	/** Checks the exact mixed mate-presence transitivity failure in ReadComparatorTopological. */
	private static void topologicalMixedMatePresence(final String expect){
		final Read a=read("AAAAAAAA", "a", 1, "TT");
		final Read b=read("AAAAAAAA", "b", 2, null);
		final Read c=read("AAAAAAAA", "c", 3, "AA");
		checkMixedMatePresence(expect, "topological",
				ReadComparatorTopological.ascending, ReadComparatorTopological.descending,
				a, b, c);
	}

	/** Checks the same mixed mate-presence failure after the 12-symbol topological key. */
	private static void topological5BitMixedMatePresence(final String expect){
		final Read a=read("AAAAAAAAAAAA", "a", 0, "TT");
		final Read b=read("AAAAAAAAAAAA", "b", 0, null);
		final Read c=read("AAAAAAAAAAAA", "c", 0, "AA");
		checkMixedMatePresence(expect, "topological5bit",
				ReadComparatorTopological5Bit.ascending, ReadComparatorTopological5Bit.descending,
				a, b, c);
	}

	/** Checks that strict primary-base prefixes are resolved before mate bases. */
	private static void topologicalStrictPrefixBeforeMate(final String expect){
		final Read a=read("A", "a", 1, "G");
		final Read b=read("AC", "b", 2, "T");
		final Read c=read("AG", "c", 3, "A");
		checkStrictPrefixBeforeMate(expect, "topological-prefix",
				ReadComparatorTopological.ascending, ReadComparatorTopological.descending,
				a, b, c);
	}

	/** Checks the same strict-prefix failure after the 12-symbol topological key. */
	private static void topological5BitStrictPrefixBeforeMate(final String expect){
		final Read a=read("AAAAAAAAAAAA", "a", 0, "G");
		final Read b=read("AAAAAAAAAAAAC", "b", 0, "T");
		final Read c=read("AAAAAAAAAAAAG", "c", 0, "A");
		checkStrictPrefixBeforeMate(expect, "topological5bit-prefix",
				ReadComparatorTopological5Bit.ascending, ReadComparatorTopological5Bit.descending,
				a, b, c);
	}

	/** Checks the public direct comparator can still ignore mates. */
	private static void topologicalCompareMatesFalse(final String expect){
		final Read a=read("AA", "a", 1, "TT");
		final Read b=read("AA", "b", 2, "AA");
		final int aToB=fixed(expect) ? -1 : 1;
		checkCompare(expect, "topological-no-mates-a-b", aToB,
				ReadComparatorTopological.ascending, a, b, false);
		checkCompare(expect, "topological-no-mates-b-a", 0-aToB,
				ReadComparatorTopological.ascending, b, a, false);
		System.out.println("PASS\t"+expect+"\ttopological-no-mates");
	}

	/** Checks the 5-bit public direct comparator can still ignore mates. */
	private static void topological5BitCompareMatesFalse(final String expect){
		final Read a=read("AAAAAAAAAAAA", "a", 0, "TT");
		final Read b=read("AAAAAAAAAAAA", "b", 0, "AA");
		final int aToB=fixed(expect) ? -1 : 1;
		checkCompare(expect, "topological5bit-no-mates-a-b", aToB,
				ReadComparatorTopological5Bit.ascending, a, b, false);
		checkCompare(expect, "topological5bit-no-mates-b-a", 0-aToB,
				ReadComparatorTopological5Bit.ascending, b, a, false);
		System.out.println("PASS\t"+expect+"\ttopological5bit-no-mates");
	}

	/** Checks baseline's {@code A < B < C && A > C} cycle or the repaired total order {@code C < A < B}. */
	private static void checkMixedMatePresence(final String expect, final String label,
			final Comparator<Read> ascending, final Comparator<Read> descending,
			final Read a, final Read b, final Read c){
		if(fixed(expect)){
			checkOrder(expect, label+"-ascending", ascending, c, a, b);
			checkOrder(expect, label+"-descending", descending, b, a, c);
		}else{
			checkCompare(expect, label+"-ascending-a-b", -1, ascending, a, b);
			checkCompare(expect, label+"-ascending-b-c", -1, ascending, b, c);
			checkCompare(expect, label+"-ascending-a-c", 1, ascending, a, c);
			checkCompare(expect, label+"-descending-a-b", 1, descending, a, b);
			checkCompare(expect, label+"-descending-b-c", 1, descending, b, c);
			checkCompare(expect, label+"-descending-a-c", -1, descending, a, c);
			System.out.println("PASS\t"+expect+"\t"+label+"-baseline-cycle");
		}
	}

	/** Checks baseline's {@code A < B < C < A} cycle or the repaired primary order {@code A < B < C}. */
	private static void checkStrictPrefixBeforeMate(final String expect, final String label,
			final Comparator<Read> ascending, final Comparator<Read> descending,
			final Read a, final Read b, final Read c){
		if(fixed(expect)){
			checkOrder(expect, label+"-ascending", ascending, a, b, c);
			checkOrder(expect, label+"-descending", descending, c, b, a);
		}else{
			checkCompare(expect, label+"-ascending-a-b", -1, ascending, a, b);
			checkCompare(expect, label+"-ascending-b-c", -1, ascending, b, c);
			checkCompare(expect, label+"-ascending-c-a", -1, ascending, c, a);
			checkCompare(expect, label+"-descending-a-b", 1, descending, a, b);
			checkCompare(expect, label+"-descending-b-c", 1, descending, b, c);
			checkCompare(expect, label+"-descending-c-a", 1, descending, c, a);
			System.out.println("PASS\t"+expect+"\t"+label+"-baseline-cycle");
		}
	}

	/** Requires {@code first < second < third} and all reverse comparisons. */
	private static void checkOrder(final String expect, final String label,
			final Comparator<Read> comparator, final Read first, final Read second,
			final Read third){
		checkCompare(expect, label+"-first-second", -1, comparator, first, second);
		checkCompare(expect, label+"-second-third", -1, comparator, second, third);
		checkCompare(expect, label+"-first-third", -1, comparator, first, third);
		checkCompare(expect, label+"-second-first", 1, comparator, second, first);
		checkCompare(expect, label+"-third-second", 1, comparator, third, second);
		checkCompare(expect, label+"-third-first", 1, comparator, third, first);
		System.out.println("PASS\t"+expect+"\t"+label+"-order");
	}

	/** Checks one comparison sign. */
	private static void checkCompare(final String expect, final String label,
			final int expected, final Comparator<Read> comparator, final Read a, final Read b){
		final int observed=sign(comparator.compare(a, b));
		check(observed==expected, expect+" "+label+" expected "+expected+" observed "+observed);
	}

	/** Checks one direct three-parameter topological comparison sign. */
	private static void checkCompare(final String expect, final String label,
			final int expected, final ReadComparatorTopological comparator,
			final Read a, final Read b, final boolean compareMates){
		final int observed=sign(comparator.compare(a, b, compareMates));
		check(observed==expected, expect+" "+label+" expected "+expected+" observed "+observed);
	}

	/** Checks one direct three-parameter 5-bit topological comparison sign. */
	private static void checkCompare(final String expect, final String label,
			final int expected, final ReadComparatorTopological5Bit comparator,
			final Read a, final Read b, final boolean compareMates){
		final int observed=sign(comparator.compare(a, b, compareMates));
		check(observed==expected, expect+" "+label+" expected "+expected+" observed "+observed);
	}

	/** Builds a read with optional mate. */
	private static Read read(final String bases, final String id, final long numericID,
			final String mateBases){
		final byte[] baseBytes=bases.getBytes(StandardCharsets.US_ASCII);
		final Read r=new Read(baseBytes, new byte[baseBytes.length], id, numericID, false);
		if(mateBases!=null){
			final byte[] mateBaseBytes=mateBases.getBytes(StandardCharsets.US_ASCII);
			final Read mate=new Read(mateBaseBytes, new byte[mateBaseBytes.length],
					id+"/2", numericID, false);
			r.mate=mate;
			mate.mate=r;
		}
		return r;
	}

	/** Fails loudly if a test expectation is violated. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Returns the sign of a comparator result. */
	private static int sign(final int value){
		return value<0 ? -1 : value>0 ? 1 : 0;
	}

	private static boolean fixed(final String expect){return "fixed".equals(expect);}

}
