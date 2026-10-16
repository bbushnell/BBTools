package test.structures;

import java.util.ArrayList;

import dna.Scaffold;
import structures.CoverageArray2A;
import structures.StandardDeviator;

/** Finite oracle for CoverageArray.devSum parameter-cache reuse.
 * @author Sucrose
 */
public final class CoverageArrayDevSumCacheCase{

	/** Runs baseline or fixed expectations for StandardDeviator cache reuse. */
	public static void main(final String[] args){
		check(args.length==1, "Expected expect=baseline or expect=fixed");
		final String[] split=args[0].split("=", 2);
		check(split.length==2 && "expect".equals(split[0]),
				"Expected expect=baseline or expect=fixed");
		final String expect=split[1];
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");

		minScafReusesAllScaffoldMean(expect);
		sameParameterControl(expect);
		zeroMappedSetterReusesAllScaffoldMean(expect);
	}

	/** Exposes devSum(mean) reuse after minscaf changes the filtered mean. */
	private static void minScafReusesAllScaffoldMean(final String expect){
		final ArrayList<Scaffold> scaffolds=new ArrayList<Scaffold>();
		scaffolds.add(scaffold("long", 4, 2));
		scaffolds.add(scaffold("short", 2, 8));
		seedDevSum(scaffolds);

		final double[] observed=stdev(scaffolds, 4);
		final double expected="baseline".equals(expect) ? 2 : 0;
		checkEquals("filtered mean", 2, observed[0]);
		checkEquals("minscaf devsum reuse", expected, observed[1]);
		pass("minscaf-devsum-reuse", observed, expect);
	}

	/** Preserves reuse when both StandardDeviator passes use the same mean. */
	private static void sameParameterControl(final String expect){
		final ArrayList<Scaffold> scaffolds=new ArrayList<Scaffold>();
		scaffolds.add(scaffold("long", 4, 2));
		seedDevSum(scaffolds);

		final double[] observed=stdev(scaffolds, 4);
		checkEquals("same mean", 2, observed[0]);
		checkEquals("same-parameter stdev", 0, observed[1]);
		pass("same-parameter-control", observed, expect);
	}

	/** Exposes zero-mapped setDevSum() reuse after filtering changes the mean. */
	private static void zeroMappedSetterReusesAllScaffoldMean(final String expect){
		final ArrayList<Scaffold> scaffolds=new ArrayList<Scaffold>();
		scaffolds.add(scaffold("long-zero", 4, 0));
		scaffolds.add(scaffold("short-high", 2, 6));
		seedDevSum(scaffolds);

		final double[] observed=stdev(scaffolds, 4);
		final double expected="baseline".equals(expect) ? 2 : 0;
		checkEquals("zero-mapped filtered mean", 0, observed[0]);
		checkEquals("zero-mapped devsum reuse", expected, observed[1]);
		pass("zero-mapped-setter-reuse", observed, expect);
	}

	/** Creates a scaffold carrying a constant-depth coverage array. */
	private static Scaffold scaffold(final String name, final int length,
			final int depth){
		final Scaffold scaffold=new Scaffold(name, length);
		final CoverageArray2A coverage=new CoverageArray2A(0, length);
		for(int i=0; i<length; i++){
			coverage.set(i, depth);
		}
		scaffold.obj0=coverage;
		scaffold.basehits=(long)length*depth;
		return scaffold;
	}

	/** Seeds aggregate caches through the multi-threaded public chain. */
	private static void seedDevSum(final ArrayList<Scaffold> scaffolds){
		final StandardDeviator sd=new StandardDeviator(false, 0);
		final boolean ok=sd.calculateStuff(1, scaffolds, true, false, false,
				1, -1, -1, -1);
		check(!ok, "StandardDeviator reported an error");
	}

	/** Recomputes filtered standard deviation through the public chain. */
	private static double[] stdev(final ArrayList<Scaffold> scaffolds,
			final int minscaf){
		final StandardDeviator sd=new StandardDeviator(false, 0);
		return sd.standardDeviation(scaffolds, minscaf);
	}

	/** Checks an exact double oracle. */
	private static void checkEquals(final String name, final double expected,
			final double observed){
		check(Double.compare(expected, observed)==0,
				name+" expected "+expected+", got "+observed);
	}

	/** Prints one successful case row. */
	private static void pass(final String name, final double[] observed,
			final String expect){
		System.out.println("PASS\t"+expect+'\t'+name+'\t'
				+observed[0]+'\t'+observed[1]);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

}
