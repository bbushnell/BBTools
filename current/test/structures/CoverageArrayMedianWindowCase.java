package test.structures;

import java.util.ArrayList;

import dna.Scaffold;
import structures.CoverageArray;
import structures.CoverageArray2;
import structures.CoverageArray2A;
import structures.CoverageArray3;
import structures.CoverageArray3A;
import structures.StandardDeviator;

/** Finite oracle for median() reordering before low-window coverage.
 * @author Sucrose
 */
public final class CoverageArrayMedianWindowCase{

	/** Runs baseline or fixed expectations for median/window interaction. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");

		checkBackend("ca2", CoverageArray2.class, true, expect);
		checkBackend("ca2a", CoverageArray2A.class, false, expect);
		checkBackend("ca3", CoverageArray3.class, true, expect);
		checkBackend("ca3a", CoverageArray3A.class, false, expect);
	}

	/** Checks one coverage backend through the real two-pass statistics path. */
	private static void checkBackend(final String name,
			final Class<? extends CoverageArray> type, final boolean baselineSorts,
			final String expect){
		final ArrayList<Scaffold> scaffolds=new ArrayList<Scaffold>();
		final CoverageArray first=array(type, 0);
		final CoverageArray second=array(type, 1);
		scaffolds.add(scaffold("first-"+name, first));
		scaffolds.add(scaffold("second-"+name, second));

		checkMinMax(name+" before", first, beforeMin(type), beforeMax(type));
		runStandardDeviator(scaffolds);

		final boolean sorted=baselineSorts && "baseline".equals(expect);
		final int[] expected=sorted ? SORTED : ORIGINAL;
		final int windowBases=sorted ? 3 : 0;
		final int min=sorted ? 0 : beforeMin(type);
		final int max=sorted ? 2 : beforeMax(type);

		checkCoverage(name+" first", first, expected);
		checkCoverage(name+" second", second, expected);
		checkEquals(name+" median", 0, first.median());
		checkEquals(name+" window bases", windowBases,
				first.basesUnderAverageCoverage(WINDOW_AVERAGE, WINDOW));
		checkMinMax(name+" after", first, min, max);
		pass(name, first, expect);
	}

	/** Creates one 5-base array with depth at positions 1 and 3. */
	private static CoverageArray array(final Class<? extends CoverageArray> type,
			final int chrom){
		final CoverageArray coverage=CoverageArray.makeArray(chrom, LENGTH, type);
		coverage.set(1, 1);
		coverage.set(3, 1);
		return coverage;
	}

	/** Creates a scaffold carrying the supplied coverage array. */
	private static Scaffold scaffold(final String name,
			final CoverageArray coverage){
		final Scaffold scaffold=new Scaffold(name, LENGTH);
		scaffold.obj0=coverage;
		scaffold.basehits=2;
		return scaffold;
	}

	/** Seeds median and low-window coverage via StandardDeviator.calculateStuff. */
	private static void runStandardDeviator(final ArrayList<Scaffold> scaffolds){
		final StandardDeviator sd=new StandardDeviator(false, 0);
		final boolean error=sd.calculateStuff(1, scaffolds, true, true, false,
				1, HIST_MAX, WINDOW_AVERAGE, WINDOW);
		check(!error, "StandardDeviator reported an error");
	}

	/** Checks exact coverage values. */
	private static void checkCoverage(final String name, final CoverageArray ca,
			final int[] expected){
		for(int i=0; i<expected.length; i++){
			checkEquals(name+" element "+i, expected[i], ca.get(i));
		}
	}

	/** Checks exact touched bounds. */
	private static void checkMinMax(final String name, final CoverageArray ca,
			final int min, final int max){
		checkEquals(name+" minIndex", min, ca.minIndex);
		checkEquals(name+" maxIndex", max, ca.maxIndex);
	}

	/** Returns the expected starting minIndex. */
	private static int beforeMin(final Class<? extends CoverageArray> type){
		return atomic(type) ? 0 : 1;
	}

	/** Returns the expected starting maxIndex. */
	private static int beforeMax(final Class<? extends CoverageArray> type){
		return atomic(type) ? LENGTH-1 : 3;
	}

	/** Returns true for fixed-length atomic coverage arrays. */
	private static boolean atomic(final Class<? extends CoverageArray> type){
		return type==CoverageArray2A.class || type==CoverageArray3A.class;
	}

	/** Checks an exact integer oracle. */
	private static void checkEquals(final String name, final int expected,
			final int observed){
		check(expected==observed, name+" expected "+expected+", got "+observed);
	}

	/** Prints one successful backend row. */
	private static void pass(final String name, final CoverageArray ca,
			final String expect){
		System.out.println("PASS\t"+expect+'\t'+name+'\t'
				+ca.median()+'\t'
				+ca.basesUnderAverageCoverage(WINDOW_AVERAGE, WINDOW)+'\t'
				+ca.minIndex+'\t'+ca.maxIndex+'\t'+ca.toString());
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final int LENGTH=5;
	private static final int HIST_MAX=8;
	private static final int WINDOW=2;
	private static final double WINDOW_AVERAGE=0.5;
	private static final int[] ORIGINAL={0, 1, 0, 1, 0};
	private static final int[] SORTED={1, 1, 0, 0, 0};

}
