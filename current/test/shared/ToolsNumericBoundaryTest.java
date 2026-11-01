package test.shared;

import java.util.concurrent.atomic.AtomicLongArray;

import shared.Tools;

/** Finite numeric boundary checks for Tools overflow and nearest-index cases.
 * Cases use literal independent expectations and no production callers.
 * @author Jean
 */
public final class ToolsNumericBoundaryTest{

	/** Runs every case with explicit checks under both assertion modes. */
	public static void main(String[] args){
		if(args.length!=0){throw new IllegalArgumentException("Expected no arguments");}
		final ToolsNumericBoundaryTest test=new ToolsNumericBoundaryTest();
		test.check("binarySearch exact keys", ToolsNumericBoundaryTest::binarySearchExactKeys);
		test.check("binarySearch nearest final interval", ToolsNumericBoundaryTest::binarySearchNearestFinalInterval);
		test.check("weightedAverage long constants", ToolsNumericBoundaryTest::weightedAverageLongConstants);
		test.check("long histogram small controls", ToolsNumericBoundaryTest::longHistogramSmallControls);
		test.check("long histogram huge single-bin average", ToolsNumericBoundaryTest::longHistogramHugeSingleBinAverage);
		test.check("long histogram huge single-bin deviation", ToolsNumericBoundaryTest::longHistogramHugeSingleBinDeviation);
		test.check("long histogram huge two-bin average", ToolsNumericBoundaryTest::longHistogramHugeTwoBinAverage);
		test.check("long histogram huge two-bin deviation", ToolsNumericBoundaryTest::longHistogramHugeTwoBinDeviation);
		test.check("long histogram huge percentile", ToolsNumericBoundaryTest::longHistogramHugePercentile);
		test.check("standardDeviationHistogram int constant", ToolsNumericBoundaryTest::standardDeviationHistogramIntConstant);
		test.check("standardDeviationHistogram char constant", ToolsNumericBoundaryTest::standardDeviationHistogramCharConstant);
		test.check("requireNear rejects NaN", ToolsNumericBoundaryTest::requireNearRejectsNaN);
		test.check("AtomicLong histogram huge single-bin deviation",
				ToolsNumericBoundaryTest::standardDeviationHistogramAtomicLongSingleBin);
		test.check("AtomicLong histogram huge two-bin deviation",
				ToolsNumericBoundaryTest::standardDeviationHistogramAtomicLongTwoBin);
		test.check("long unique-kmer histogram huge two-bin deviation",
				ToolsNumericBoundaryTest::standardDeviationHistogramKmerLongTwoBin);
		test.check("AtomicLong unique-kmer histogram huge two-bin deviation",
				ToolsNumericBoundaryTest::standardDeviationHistogramKmerAtomicLongTwoBin);
		test.check("murmur32plus sign mask", ToolsNumericBoundaryTest::murmur32plusSignMask);
		test.check("murmur32plus below upper fold", ToolsNumericBoundaryTest::murmur32plusBelowUpperFold);
		test.check("murmur32plus at upper fold", ToolsNumericBoundaryTest::murmur32plusAtUpperFold);
		test.check("murmur32plus last upper fold", ToolsNumericBoundaryTest::murmur32plusLastUpperFold);
		test.check("murmur32plus masked to upper fold", ToolsNumericBoundaryTest::murmur32plusMaskedToUpperFold);
		test.check("murmur32plus masked to last upper fold",
				ToolsNumericBoundaryTest::murmur32plusMaskedToLastUpperFold);
		System.out.println("ToolsNumericBoundaryTest cases="+test.cases+" failures="+test.failures);
		if(test.failures!=0){throw new AssertionError("Tools numeric boundary failures="+test.failures);}
	}

	/** Exact-key controls for the three-element float search. */
	private static void binarySearchExactKeys(){
		final float[] array={0, 10, 20};
		require(Tools.binarySearch(array, 0)==0, "exact 0 returned wrong index");
		require(Tools.binarySearch(array, 10)==1, "exact 10 returned wrong index");
		require(Tools.binarySearch(array, 20)==2, "exact 20 returned wrong index");
	}

	/** Keys on each side of the final midpoint must choose the nearest element. */
	private static void binarySearchNearestFinalInterval(){
		final float[] array={0, 10, 20};
		require(Tools.binarySearch(array, 14)==1, "key 14 must be nearest to index 1");
		require(Tools.binarySearch(array, 16)==2, "key 16 must be nearest to index 2");
	}

	/** Constant weighted averages must equal that constant, even near Long.MAX_VALUE. */
	private static void weightedAverageLongConstants(){
		require(Tools.weightedAverage(new long[] {4, 4})==4.0, "small constant average changed");
		final double expected=(double)Long.MAX_VALUE;
		final double observed=Tools.weightedAverage(new long[] {Long.MAX_VALUE, Long.MAX_VALUE});
		require(observed==expected, "Long.MAX_VALUE constant average was "+observed);
	}

	/** Small long histogram controls with representable counts and weighted products. */
	private static void longHistogramSmallControls(){
		final long[] singleBin=new long[4];
		singleBin[3]=7;
		require(Tools.averageHistogram(singleBin)==3, "small long single-bin average changed");
		require(Tools.standardDeviationHistogram(singleBin)==0, "small long single-bin deviation was nonzero");

		final long[] doubleBin={0, 5, 5};
		require(Tools.averageHistogram(doubleBin)==1.5, "small long two-bin average changed");
		require(Tools.standardDeviationHistogram(doubleBin)==0.5, "small long two-bin deviation changed");
	}

	/** A huge single populated bin must not overflow the weighted sum. */
	private static void longHistogramHugeSingleBinAverage(){
		require(Tools.averageHistogram(hugeSingleBin())==3, "huge long single-bin average changed");
	}

	/** A huge single populated bin must have zero dispersion. */
	private static void longHistogramHugeSingleBinDeviation(){
		require(Tools.standardDeviationHistogram(hugeSingleBin())==0, "huge long single-bin deviation was nonzero");
	}

	/** Huge adjacent bins must not overflow the count sum for average. */
	private static void longHistogramHugeTwoBinAverage(){
		require(Tools.averageHistogram(hugeDoubleBin())==1.5, "huge long two-bin average changed");
	}

	/** Huge adjacent bins must not overflow the count sum for dispersion. */
	private static void longHistogramHugeTwoBinDeviation(){
		require(Tools.standardDeviationHistogram(hugeDoubleBin())==0.5, "huge long two-bin deviation changed");
	}

	/** Huge adjacent bins must not overflow the count sum for percentile targets. */
	private static void longHistogramHugePercentile(){
		require(Tools.percentileHistogram(new long[] {0, 5, 5}, 0.75)==2,
				"small long 75th percentile changed");
		require(Tools.percentileHistogram(hugeDoubleBin(), 0)==0,
				"huge long zero percentile changed");
		require(Tools.percentileHistogram(hugeDoubleBin(), Double.MIN_VALUE)==0,
				"huge long tiny positive percentile changed");
		require(Tools.percentileHistogram(hugeDoubleBin(), 0.5)==1,
				"huge long 50th percentile changed");
		require(Tools.percentileHistogram(hugeDoubleBin(), Math.nextUp(0.5))==2,
				"huge long next-up 50th percentile changed");
		require(Tools.percentileHistogram(hugeDoubleBin(), 0.75)==2,
				"huge long 75th percentile changed");
		require(Tools.percentileHistogram(hugeDoubleBin(), 1)==2,
				"huge long 100th percentile changed");
		require(Tools.percentileHistogram(hugeDoubleBinWithTrailingZero(), 1)==2,
				"huge long 100th percentile with trailing zero changed");
	}

	/** Histogram with one extreme bin at position 3. */
	private static long[] hugeSingleBin(){
		final long[] histogram=new long[4];
		histogram[3]=Long.MAX_VALUE;
		return histogram;
	}

	/** Histogram with two adjacent extreme bins at positions 1 and 2. */
	private static long[] hugeDoubleBin(){
		return new long[] {0, Long.MAX_VALUE, Long.MAX_VALUE};
	}

	/** Histogram with two adjacent extreme bins followed by an empty bin. */
	private static long[] hugeDoubleBinWithTrailingZero(){
		return new long[] {0, Long.MAX_VALUE, Long.MAX_VALUE, 0};
	}

	/** A one-bin int histogram has zero dispersion regardless of count. */
	private static void standardDeviationHistogramIntConstant(){
		final int[] control=new int[4];
		control[3]=2;
		require(Tools.standardDeviationHistogram(control)==0, "small int constant histogram was nonzero");
		final int[] histogram=new int[3001];
		histogram[3000]=1000000;
		require(Tools.standardDeviationHistogram(histogram)==0, "large int constant histogram was nonzero");
	}

	/** A one-bin char histogram has zero dispersion regardless of count. */
	private static void standardDeviationHistogramCharConstant(){
		final char[] control=new char[6];
		control[5]=3;
		require(Tools.standardDeviationHistogram(control)==0, "small char constant histogram was nonzero");
		final char[] histogram=new char[50001];
		histogram[50000]=(char)65535;
		require(Tools.standardDeviationHistogram(histogram)==0, "large char constant histogram was nonzero");
	}

	/** Verifies that approximate comparisons cannot accept a NaN observation. */
	private static void requireNearRejectsNaN(){
		boolean rejected=false;
		try{
			requireNear(Double.NaN, Math.sqrt(2), 0.000000000001, "NaN should fail");
		}catch(AssertionError e){
			rejected=true;
		}
		require(rejected, "requireNear accepted NaN");
	}

	/** AtomicLong histogram with one extreme bin at position 3. */
	private static void standardDeviationHistogramAtomicLongSingleBin(){
		final AtomicLongArray control=new AtomicLongArray(4);
		control.set(3, 7);
		require(Tools.standardDeviationHistogram(control)==0, "small AtomicLong single-bin deviation was nonzero");
		final AtomicLongArray histogram=new AtomicLongArray(4);
		histogram.set(3, Long.MAX_VALUE);
		require(Tools.standardDeviationHistogram(histogram)==0, "huge AtomicLong single-bin deviation was nonzero");
	}

	/** AtomicLong histogram with two adjacent extreme bins at positions 1 and 2. */
	private static void standardDeviationHistogramAtomicLongTwoBin(){
		final AtomicLongArray control=new AtomicLongArray(new long[] {0, 5, 5});
		require(Tools.standardDeviationHistogram(control)==0.5, "small AtomicLong two-bin deviation changed");
		final AtomicLongArray histogram=new AtomicLongArray(
				new long[] {0, Long.MAX_VALUE, Long.MAX_VALUE});
		require(Tools.standardDeviationHistogram(histogram)==0.5, "huge AtomicLong two-bin deviation changed");
	}

	/** Unique-kmer long histogram with two extreme bins at depths 3 and 6. */
	private static void standardDeviationHistogramKmerLongTwoBin(){
		requireNear(Tools.standardDeviationHistogramKmer(smallKmerDoubleBin()), Math.sqrt(2), 0.000000000001,
				"small long unique-kmer two-bin deviation changed");
		requireNear(Tools.standardDeviationHistogramKmer(hugeKmerDoubleBin()), Math.sqrt(2), 0.000000000001,
				"huge long unique-kmer two-bin deviation changed");
	}

	/** Unique-kmer AtomicLong histogram with two extreme bins at depths 3 and 6. */
	private static void standardDeviationHistogramKmerAtomicLongTwoBin(){
		requireNear(Tools.standardDeviationHistogramKmer(new AtomicLongArray(smallKmerDoubleBin())), Math.sqrt(2),
				0.000000000001, "small AtomicLong unique-kmer two-bin deviation changed");
		requireNear(Tools.standardDeviationHistogramKmer(new AtomicLongArray(hugeKmerDoubleBin())), Math.sqrt(2),
				0.000000000001, "huge AtomicLong unique-kmer two-bin deviation changed");
	}

	/** Input selected by inverting the 32-bit finalizer to raw hash 0x80000000. */
	private static void murmur32plusSignMask(){
		require(Tools.murmur32plus(0x7EC69360)==0,
				"raw hash 0x80000000 must clear the sign bit to zero");
	}

	/** Input selected by inverting the 32-bit finalizer to raw hash 0x7FFFF7FF. */
	private static void murmur32plusBelowUpperFold(){
		require(Tools.murmur32plus(0x5650E5CD)==0x7FFFF7FF,
				"raw hash 0x7FFFF7FF must remain below the upper fold");
	}

	/** Input selected by inverting the 32-bit finalizer to raw hash 0x7FFFF800. */
	private static void murmur32plusAtUpperFold(){
		require(Tools.murmur32plus(0xD4FBAEDA)==0,
				"raw hash 0x7FFFF800 must fold to zero");
	}

	/** Input selected by inverting the 32-bit finalizer to raw hash 0x7FFFFFFF. */
	private static void murmur32plusLastUpperFold(){
		require(Tools.murmur32plus(0x9C0D8E22)==131008,
				"raw hash 0x7FFFFFFF must fold to 131008");
	}

	/** Input selected by inverting the 32-bit finalizer to raw hash 0xFFFFF800. */
	private static void murmur32plusMaskedToUpperFold(){
		require(Tools.murmur32plus(0xC4FD90D1)==0,
				"raw hash 0xFFFFF800 must mask into the upper fold");
	}

	/** Input selected by inverting the 32-bit finalizer to raw hash 0xFFFFFFFF. */
	private static void murmur32plusMaskedToLastUpperFold(){
		require(Tools.murmur32plus(0x331DA083)==131008,
				"raw hash 0xFFFFFFFF must mask into the last upper fold");
	}

	/** Histogram with two low-count equal-copy bins at depths 3 and 6. */
	private static long[] smallKmerDoubleBin(){
		final long[] histogram=new long[7];
		histogram[3]=6;
		histogram[6]=6;
		return histogram;
	}

	/** Histogram with two extreme equal-copy bins at depths 3 and 6. */
	private static long[] hugeKmerDoubleBin(){
		final long[] histogram=new long[7];
		histogram[3]=Long.MAX_VALUE;
		histogram[6]=Long.MAX_VALUE;
		return histogram;
	}

	/** Runs one case and records rather than short-circuiting failures. */
	private void check(String name, Action action){
		cases++;
		try{
			action.run();
			System.out.println(name+" PASS");
		}catch(AssertionError e){
			failures++;
			System.out.println(name+" FAIL: "+e.getMessage());
		}
	}

	/** Checks expectations in both assertion modes. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Checks a finite floating-point expectation within a named absolute tolerance. */
	private static void requireNear(double observed, double expected, double tolerance, String message){
		if(!(Math.abs(observed-expected)<=tolerance)){
			throw new AssertionError(message+"; observed="+observed+", expected="+expected+", tolerance="+tolerance);
		}
	}

	/** One finite Tools numeric contract case. */
	private interface Action{void run();}

	/** Attempted cases. */
	private int cases;
	/** Failed cases. */
	private int failures;
}
