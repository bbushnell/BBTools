package test.structures;

import java.util.ArrayList;

import structures.Feature;
import structures.Range;

/** Finite midpoint-overflow oracle for Range.findIndexAfterBinary.
 * @author Sucrose
 */
public final class RangeMidpointOverflowCase{

	/** Runs baseline or fixed expectations for midpoint-overflow cases. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		final boolean fixed="fixed".equals(expect);

		checkObservation(expect, "empty", "-1", observe(0, 0));
		checkObservation(expect, "single_equal", "1", observe(0, 1));
		checkObservation(expect, "small_first", "1", observe(0, 12));
		checkObservation(expect, "small_between", "6", observe(5, 12));
		checkObservation(expect, "small_last", "12", observe(11, 12));
		checkObservation(expect, "large_low", "11", observe(10, 1500000000));
		checkObservation(expect, "large_overflow_first",
				fixed ? "750000000" : "IndexOutOfBoundsException",
				observe(749999999, 1500000000));
		checkObservation(expect, "large_overflow_near_end",
				fixed ? "1499999999" : "IndexOutOfBoundsException",
				observe(1499999998, 1500000000));
		checkObservation(expect, "large_overflow_after_end",
				fixed ? "1500000000" : "IndexOutOfBoundsException",
				observe(1499999999, 1500000000));
	}

	/** Calls Range.findIndexAfterBinary and returns its result or thrown type. */
	private static String observe(final int p, final int size){
		try{
			final LogicalFeatureList list=new LogicalFeatureList(size);
			return Integer.toString(Range.findIndexAfterBinary(p, list));
		}catch(final Throwable t){return t.getClass().getSimpleName();}
	}

	/** Checks and prints a value-or-exception outcome. */
	private static void checkObservation(final String expect, final String name,
			final String expected, final String observed){
		check(expected.equals(observed), name+" expected "+expected+", got "+observed);
		System.out.println("PASS	"+expect+'	'+name+'	'+observed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Feature list with a huge logical size and O(log n) allocation. */
	private static final class LogicalFeatureList extends ArrayList<LiteralFeature>{

		/** Constructs a list whose feature start equals its index. */
		private LogicalFeatureList(final int size_){
			super(0);
			check(size_>=0, "Negative logical size: "+size_);
			size=size_;
		}

		@Override
		public boolean isEmpty(){return size==0;}

		@Override
		public int size(){return size;}

		@Override
		public LiteralFeature get(final int index){
			if(index<0 || index>=size){
				throw new IndexOutOfBoundsException("Index "+index+" out of bounds [0,"+size+")");
			}
			check(++gets<=128, "Too many logical feature reads: "+gets);
			return new LiteralFeature(index);
		}

		private final int size;
		private int gets;
	}

	/** Literal sorted feature with a start coordinate sufficient for binary search. */
	private static final class LiteralFeature implements Feature{

		/** Constructs a feature whose start and stop are the same point. */
		private LiteralFeature(final int start_){
			start=start_;
		}

		@Override
		public int start(){return start;}

		@Override
		public int stop(){return start;}

		@Override
		public int strand(){return 0;}

		@Override
		public String seqid(){return "literal";}

		@Override
		public String name(){return Integer.toString(start);}

		private final int start;
	}
}
