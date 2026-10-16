package test.structures;

import structures.CoverageArray2A;

/** Finite oracle for CoverageArray.toGraph final-block emission.
 * @author Sucrose
 */
public final class CoverageArrayToGraphCase{

	/** Runs baseline or fixed expectations for graph block emission. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");

		exactlyOneFullBlock(expect);
		partialFinalBlock(expect);
		multipleBlocks(expect);
		emptyClampedIntervalControl(expect);
	}

	/** Exposes a single full block being accumulated but never emitted. */
	private static void exactlyOneFullBlock(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 4);
		set(ca, 2, 4, 6, 8);
		final double[][] expected="baseline".equals(expect) ? rows() :
				rows(row(2, 5));
		checkRows("exactly-one-full-block", expected, ca.toGraph(4, 0, 3));
		pass("exactly-one-full-block", expected, expect);
	}

	/** Exposes a partial tail being accumulated but never emitted. */
	private static void partialFinalBlock(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 5);
		set(ca, 2, 4, 6, 8, 10);
		final double[][] expected="baseline".equals(expect) ? rows(row(2, 5)) :
				rows(row(2, 5), row(4, 10));
		checkRows("partial-final-block", expected, ca.toGraph(4, 0, 4));
		pass("partial-final-block", expected, expect);
	}

	/** Exposes the last block being dropped after earlier blocks are emitted. */
	private static void multipleBlocks(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 9);
		set(ca, 3, 6, 9, 12, 15, 18, 21, 24, 27);
		final double[][] expected="baseline".equals(expect) ?
				rows(row(1, 6), row(4, 15)) :
				rows(row(1, 6), row(4, 15), row(7, 24));
		checkRows("multiple-blocks", expected, ca.toGraph(3, 0, 8));
		pass("multiple-blocks", expected, expect);
	}

	/** Preserves empty results after clamping to a disjoint interval. */
	private static void emptyClampedIntervalControl(final String expect){
		final CoverageArray2A ca=new CoverageArray2A(0, 6);
		ca.set(4, 4);
		ca.set(5, 6);
		final double[][] expected=rows();
		checkRows("empty-clamped-interval-control", expected, ca.toGraph(2, 6, 9));
		pass("empty-clamped-interval-control", expected, expect);
	}

	/** Sets the leading coverage depths. */
	private static void set(final CoverageArray2A ca, final int... values){
		for(int i=0; i<values.length; i++){
			ca.set(i, values[i]);
		}
	}

	/** Creates expected output rows. */
	private static double[][] rows(final double[]... rows){return rows;}

	/** Creates one expected output row. */
	private static double[] row(final double x, final double y){
		return new double[] {x, y};
	}

	/** Checks every row and column exactly. */
	private static void checkRows(final String name, final double[][] expected,
			final double[][] observed){
		checkEquals(name+" row count", expected.length, observed.length);
		for(int r=0; r<expected.length; r++){
			checkEquals(name+" row "+r+" width", expected[r].length,
					observed[r].length);
			for(int c=0; c<expected[r].length; c++){
				checkEquals(name+" row "+r+" column "+c, expected[r][c],
						observed[r][c]);
			}
		}
	}

	/** Checks an int oracle. */
	private static void checkEquals(final String name, final int expected,
			final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Checks a double oracle. */
	private static void checkEquals(final String name, final double expected,
			final double observed){
		check(Double.compare(expected, observed)==0,
				name+" expected "+expected+", got "+observed);
	}

	/** Prints one successful case row. */
	private static void pass(final String name, final double[][] expected,
			final String expect){
		System.out.println("PASS\t"+expect+'\t'+name+'\t'+toString(expected));
	}

	/** Formats graph rows for literal stdout review. */
	private static String toString(final double[][] rows){
		final StringBuilder sb=new StringBuilder("[");
		for(int r=0; r<rows.length; r++){
			if(r>0){sb.append(", ");}
			sb.append('[').append(rows[r][0]).append(", ").append(rows[r][1]).append(']');
		}
		return sb.append(']').toString();
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

}
