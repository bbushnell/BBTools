package test.ml;

import ml.CellNet;
import ml.Sample;
import ml.SampleValueComparator;
import ml.Seed;

/**
 * Direct comparator-law checks for CellNet, Seed and SampleValueComparator.
 *
 * <p>These cases set fields directly to exercise comparator contracts only.
 * They do not show that NaN predictions can reach live training metrics.
 *
 * @author Sucrose
 */
public final class ComparatorLawTest{

	/** Runs one selectable comparator-law case. */
	public static void main(String[] args){
		if(args.length!=1){throw new IllegalArgumentException("Expected one case name");}
		final String name=args[0];
		final int oldCompareCode=CellNet.compareCode;
		try{
			CellNet.compareCode=CellNet.compareWER;
			run(name);
		}finally{
			CellNet.compareCode=oldCompareCode;
		}
		System.err.println("ComparatorLawTest PASS "+name);
	}

	/** Dispatches one direct comparator case. */
	private static void run(final String name){
		if(name.equals("cellnet-finite")){
			checkCellNetFinite();
		}else if(name.equals("cellnet-selected-modes")){
			checkCellNetSelectedModes();
		}else if(name.equals("cellnet-nan-left")){
			assertLeftWorse("cellnet-nan-left", net(1, Float.NaN, 2, 3, 4),
				net(2, 1, 2, 3, 4));
		}else if(name.equals("cellnet-nan-right")){
			assertLeftBetter("cellnet-nan-right", net(1, 1, 2, 3, 4),
				net(2, Float.NaN, 2, 3, 4));
		}else if(name.equals("cellnet-both-nan-fallback")){
			assertLeftBetter("cellnet-both-nan-fallback", net(1, Float.NaN, 1, 3, 4),
				net(2, Float.NaN, 2, 3, 4));
		}else if(name.equals("cellnet-nan-primary-fallback-tie")){
			assertLeftBetter("cellnet-nan-primary-fallback-tie",
				net(1, Float.NaN, Float.NaN, 1, 4),
				net(2, Float.NaN, Float.NaN, 2, 4));
		}else if(name.equals("cellnet-signed-zero-fallback")){
			assertLeftWorse("cellnet-signed-zero-fallback",
				net(1, -0.0f, 2, 3, 4), net(2, 0.0f, 1, 3, 4));
		}else if(name.equals("cellnet-equality")){
			checkCellNetEquality();
		}else if(name.equals("cellnet-triple")){
			checkCellNetTriple();
		}else if(name.equals("seed-finite")){
			assertAscending("seed-finite", new Seed(1, 1), new Seed(2, 2));
		}else if(name.equals("seed-nan-left")){
			assertAscending("seed-nan-left", new Seed(2, 1), new Seed(1, Float.NaN));
		}else if(name.equals("seed-nan-right")){
			assertAscending("seed-nan-right", new Seed(1, 1), new Seed(2, Float.NaN));
		}else if(name.equals("seed-both-nan-fallback")){
			assertAscending("seed-both-nan-fallback", new Seed(1, Float.NaN),
				new Seed(2, Float.NaN));
		}else if(name.equals("seed-signed-zero-fallback")){
			assertAscending("seed-signed-zero-fallback", new Seed(1, 0.0f),
				new Seed(2, -0.0f));
		}else if(name.equals("seed-equality")){
			checkSeedEquality();
		}else if(name.equals("seed-triple")){
			checkSeedTriple();
		}else if(name.equals("sample-finite")){
			assertSampleAscending("sample-finite", sample(1, 1), sample(2, 2));
		}else if(name.equals("sample-nan-left")){
			assertSampleAscending("sample-nan-left", sample(7, 1), sample(2, Float.NaN));
		}else if(name.equals("sample-nan-right")){
			assertSampleAscending("sample-nan-right", sample(8, 1), sample(3, Float.NaN));
		}else if(name.equals("sample-both-nan-fallback")){
			assertSampleAscending("sample-both-nan-fallback", sample(1, Float.NaN),
				sample(2, Float.NaN));
		}else if(name.equals("sample-signed-zero-fallback")){
			assertSampleAscending("sample-signed-zero-fallback", sample(1, 0.0f),
				sample(2, -0.0f));
		}else if(name.equals("sample-equality")){
			checkSampleEquality();
		}else if(name.equals("sample-triple")){
			checkSampleTriple();
		}else{
			throw new IllegalArgumentException("Unknown case: "+name);
		}
	}

	/** Checks that lower finite primary metrics win in both directions. */
	private static void checkCellNetFinite(){
		assertLeftBetter("cellnet-finite", net(1, 1, 3, 5, 7), net(2, 2, 3, 5, 7));
	}

	/** Checks every selectable primary metric can override opposing fallback metrics. */
	private static void checkCellNetSelectedModes(){
		checkCellNetSelected("wer", CellNet.compareWER,
			net(1, 1, 9, 9, 9), net(2, 2, 1, 1, 1));
		checkCellNetSelected("err", CellNet.compareERR,
			net(1, 9, 1, 9, 9), net(2, 1, 2, 1, 1));
		checkCellNetSelected("fnr", CellNet.compareFNR,
			net(1, 9, 9, 1, 9), net(2, 1, 1, 2, 1));
		checkCellNetSelected("fpr", CellNet.compareFPR,
			net(1, 9, 9, 9, 1), net(2, 1, 1, 1, 2));
	}

	/** Requires the selected primary metric to beat the fallback chain. */
	private static void checkCellNetSelected(final String name, final int compareCode,
			final CellNet a, final CellNet b){
		CellNet.compareCode=compareCode;
		assertLeftBetter("cellnet-selected-"+name, a, b);
	}

	/** Checks NaN reflexivity and equality for distinct matching networks. */
	private static void checkCellNetEquality(){
		final CellNet nan=net(1, Float.NaN, 2, 3, 4);
		final CellNet a=net(2, 1, 2, 3, 4), b=net(2, 1, 2, 3, 4);
		assertEqual("cellnet-nan-self", nan.compareTo(nan));
		assertEqual("cellnet-equal-forward", a.compareTo(b));
		assertEqual("cellnet-equal-reverse", b.compareTo(a));
	}

	/** Checks transitivity over finite best, finite worse, and NaN worst networks. */
	private static void checkCellNetTriple(){
		final CellNet best=net(1, 1, 3, 5, 7);
		final CellNet finite=net(2, 2, 3, 5, 7);
		final CellNet nan=net(3, Float.NaN, 3, 5, 7);
		assertLeftBetter("best-finite", best, finite);
		assertLeftBetter("finite-nan", finite, nan);
		assertLeftBetter("best-nan", best, nan);
	}

	/** Checks NaN reflexivity and equality for distinct matching Seeds. */
	private static void checkSeedEquality(){
		final Seed nan=new Seed(1, Float.NaN);
		final Seed a=new Seed(1, 1), b=new Seed(1, 1);
		assertEqual("seed-nan-self", nan.compareTo(nan));
		assertEqual("seed-equal-forward", a.compareTo(b));
		assertEqual("seed-equal-reverse", b.compareTo(a));
	}

	/** Checks transitivity over finite low, finite high and NaN worst Seeds. */
	private static void checkSeedTriple(){
		final Seed low=new Seed(1, 1);
		final Seed high=new Seed(2, 2);
		final Seed nan=new Seed(3, Float.NaN);
		assertAscending("seed-low-high", low, high);
		assertAscending("seed-high-nan", high, nan);
		assertAscending("seed-low-nan", low, nan);
	}

	/** Checks NaN reflexivity and equality for distinct matching Samples. */
	private static void checkSampleEquality(){
		final Sample nan=sample(1, Float.NaN);
		final Sample a=sample(1, 1), b=sample(1, 1);
		assertEqual("sample-nan-self", SampleValueComparator.COMPARATOR.compare(nan, nan));
		assertEqual("sample-equal-forward", SampleValueComparator.COMPARATOR.compare(a, b));
		assertEqual("sample-equal-reverse", SampleValueComparator.COMPARATOR.compare(b, a));
	}

	/** Checks transitivity for finite Sample values before NaN values. */
	private static void checkSampleTriple(){
		final Sample low=sample(3, 0);
		final Sample high=sample(1, 1);
		final Sample nan=sample(2, Float.NaN);
		assertSampleAscending("sample-low-high", low, high);
		assertSampleAscending("sample-high-nan", high, nan);
		assertSampleAscending("sample-low-nan", low, nan);
	}

	/** Builds a tiny net and installs comparison metrics directly. */
	private static CellNet net(final long seed, final float wer, final float err,
			final float fnr, final float fpr){
		final CellNet net=new CellNet(new int[]{1, 1}, seed, 1f, 0f, 1, null);
		net.weightedErrorRate=wer;
		net.errorRate=err;
		net.fnRate=fnr;
		net.fpRate=fpr;
		return net;
	}

	/** Builds a minimal Sample and installs its primary prediction directly. */
	private static Sample sample(final int id, final float value){
		final Sample sample=new Sample(new float[]{0}, new float[]{0}, 1, id);
		sample.result[0]=value;
		return sample;
	}

	/** Requires {@code a} to beat {@code b}, and {@code b} to lose to {@code a}. */
	private static void assertLeftBetter(final String name, final CellNet a, final CellNet b){
		assertPositive(name, a.compareTo(b));
		assertNegative(name+" reverse", b.compareTo(a));
	}

	/** Requires {@code a} to lose to {@code b}, and {@code b} to beat {@code a}. */
	private static void assertLeftWorse(final String name, final CellNet a, final CellNet b){
		assertNegative(name, a.compareTo(b));
		assertPositive(name+" reverse", b.compareTo(a));
	}

	/** Requires the first Seed to sort before the second in both directions. */
	private static void assertAscending(final String name, final Seed a, final Seed b){
		assertNegative(name, a.compareTo(b));
		assertPositive(name+" reverse", b.compareTo(a));
	}

	/** Requires the first Sample to sort before the second in both directions. */
	private static void assertSampleAscending(final String name, final Sample a, final Sample b){
		assertNegative(name, SampleValueComparator.COMPARATOR.compare(a, b));
		assertPositive(name+" reverse", SampleValueComparator.COMPARATOR.compare(b, a));
	}

	/** Requires a positive comparison result. */
	private static void assertPositive(final String name, final int value){
		if(value<=0){throw new AssertionError(name+" expected >0 but observed "+value);}
	}

	/** Requires a negative comparison result. */
	private static void assertNegative(final String name, final int value){
		if(value>=0){throw new AssertionError(name+" expected <0 but observed "+value);}
	}

	/** Requires equality. */
	private static void assertEqual(final String name, final int value){
		if(value!=0){throw new AssertionError(name+" expected 0 but observed "+value);}
	}

	private ComparatorLawTest(){}

}
