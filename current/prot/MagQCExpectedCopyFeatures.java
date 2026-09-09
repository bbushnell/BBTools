package prot;

/**
 * Computes the two table-derived subnet features defined by METHODS Part II.B.
 *
 * <p>For each tracked item {@code j}, {@code observed[j]} is the whole-bin
 * observed copy count and {@code expected[j]} is the frozen population
 * expectation.  The expectation may be fractional.  The capped observed term
 * is {@code min(observed[j], expected[j])}; excess is
 * {@code max(0, observed[j] - expected[j])}.  The returned ratios use the
 * common denominator {@code E=sum(expected)} and deliberately return zero for
 * {@code E==0} (the project's 0/0 convention).</p>
 *
 * <p>This class has no file or bundle dependency.  A future table consumer can
 * bind its validated, hash-checked rows to these arrays without duplicating the
 * feature arithmetic in the vectorizer or runtime tool.</p>
 */
public final class MagQCExpectedCopyFeatures {

	/** Stable feature order for later subnet/composite bundle metadata. */
	public static final String FEATURE_OBSERVED_EXPECTED="derived_observed_expected";
	public static final String FEATURE_EXCESS_EXPECTED="derived_excess_expected";
	public static final String FEATURE_UNIT="ratio";
	private static final String[] FEATURE_NAMES={FEATURE_OBSERVED_EXPECTED, FEATURE_EXCESS_EXPECTED};

	private MagQCExpectedCopyFeatures(){/* utility class */}

	/** Returns the two feature names in their emitted order. */
	public static String[] featureNames(){return FEATURE_NAMES.clone();}

	/** Immutable result for one subnet and one bin. */
	public static final class Result {
		public final double expected;
		public final double observedCapped;
		public final double excess;
		public final double observedExpected;
		public final double excessExpected;

		private Result(double expected, double observedCapped, double excess){
			if(!finite(expected) || !finite(observedCapped) || !finite(excess)){
				throw new IllegalArgumentException("non-finite derived sum");
			}
			this.expected=expected;
			this.observedCapped=observedCapped;
			this.excess=excess;
			if(expected>0){
				observedExpected=observedCapped/expected;
				excessExpected=excess/expected;
			}else{
				// Explicit METHODS convention: 0/0 is 0, never NaN/Infinity.
				observedExpected=0;
				excessExpected=0;
			}
			if(!finite(observedExpected) || !finite(excessExpected)){
				throw new IllegalArgumentException("non-finite derived ratio");
			}
		}
	}

	/**
	 * A prevalidated table-to-subnet mapping.  Compile one mapping per subnet
	 * (and context row), then reuse it for every bin; the hot path only visits
	 * the ordered tracked items and never allocates or scans untracked items.
	 */
	public static final class Mapping {
		private final int sourceLength;
		private final int[] itemIndexes;
		private final double[] expected;

		private Mapping(int sourceLength, int[] itemIndexes, double[] expected){
			this.sourceLength=sourceLength;
			this.itemIndexes=itemIndexes;
			this.expected=expected;
		}

		/** Computes features from whole-bin observations using this frozen mapping. */
		public Result compute(int[] observedByItem){
			if(observedByItem==null || observedByItem.length!=sourceLength){
				throw new IllegalArgumentException("observed item width mismatch");
			}
			double eSum=0, oSum=0, xSum=0;
			for(int i=0; i<itemIndexes.length; i++){
				final int c=observedByItem[itemIndexes[i]];
				if(c<0){throw new IllegalArgumentException("negative observed count at item "+itemIndexes[i]+": "+c);}
				final double e=expected[i];
				eSum+=e;
				oSum+=Math.min(c, e);
				xSum+=Math.max(0, c-e);
				if(!finite(eSum) || !finite(oSum) || !finite(xSum)){
					throw new IllegalArgumentException("non-finite derived sum");
				}
			}
			return new Result(eSum, oSum, xSum);
		}

		/** Number of tracked items in this mapping. */
		public int size(){return itemIndexes.length;}
	}

	/**
	 * Compiles an exact ordered item mapping against one validated expectation
	 * row.  The full expectation row is checked once here; callers should retain
	 * the returned object for repeated bin calculations.
	 */
	public static Mapping compileMapping(int sourceLength, int[] orderedItemIndexes, double[] expectedByItem){
		if(sourceLength<0 || orderedItemIndexes==null || expectedByItem==null || expectedByItem.length!=sourceLength){
			throw new IllegalArgumentException("invalid item mapping dimensions");
		}
		for(int i=0; i<expectedByItem.length; i++){
			final double e=expectedByItem[i];
			if(e<0 || !finite(e)){throw new IllegalArgumentException("invalid expected count at "+i+": "+e);}
		}
		final int[] indexes=orderedItemIndexes.clone();
		final double[] selected=new double[indexes.length];
		for(int i=0; i<indexes.length; i++){
			final int index=indexes[i];
			if(index<0 || index>=sourceLength){
				throw new IllegalArgumentException("tracked item out of range at "+i+": "+index);
			}
			for(int j=0; j<i; j++){
				if(indexes[j]==index){throw new IllegalArgumentException("duplicate tracked item: "+index);}
			}
			selected[i]=expectedByItem[index];
		}
		return new Mapping(sourceLength, indexes, selected);
	}

	/**
	 * Computes O/E and X/E for one ordered subnet item list.
	 *
	 * @param observed nonnegative whole-bin copy counts, in table item order
	 * @param expected finite nonnegative population means, same order
	 * @return the sums and their two ratios
	 * @throws IllegalArgumentException for null, unequal, negative, or non-finite inputs
	 */
	public static Result compute(int[] observed, double[] expected){
		validateArrays(observed, expected);
		double eSum=0, oSum=0, xSum=0;
		for(int i=0; i<observed.length; i++){
			final int c=observed[i];
			final double e=expected[i];
			eSum+=e;
			oSum+=Math.min(c, e);
			xSum+=Math.max(0, c-e);
			if(!finite(eSum) || !finite(oSum) || !finite(xSum)){
				throw new IllegalArgumentException("non-finite derived sum");
			}
		}
		return new Result(eSum, oSum, xSum);
	}

	/**
	 * Computes the same features for an exact ordered list of tracked item
	 * indexes.  Counts and expectations may contain additional, nontracked
	 * items; those items are ignored rather than silently becoming part of the
	 * subset. Duplicate or out-of-range indexes fail closed. For repeated bins,
	 * use {@link #compileMapping(int, int[], double[])} instead.
	 */
	public static Result computeForItems(int[] observedByItem, double[] expectedByItem, int[] orderedItemIndexes){
		validateArrays(observedByItem, expectedByItem);
		return compileMapping(observedByItem.length, orderedItemIndexes, expectedByItem).compute(observedByItem);
	}

	private static void validateArrays(int[] observed, double[] expected){
		if(observed==null || expected==null){
			throw new IllegalArgumentException("observed and expected are required");
		}
		if(observed.length!=expected.length){
			throw new IllegalArgumentException("observed/expected length mismatch: "+
				observed.length+" != "+expected.length);
		}
		for(int i=0; i<observed.length; i++){
			final int c=observed[i];
			final double e=expected[i];
			if(c<0){throw new IllegalArgumentException("negative observed count at "+i+": "+c);}
			if(e<0 || Double.isNaN(e) || Double.isInfinite(e)){
				throw new IllegalArgumentException("invalid expected count at "+i+": "+e);
			}
		}
	}

	private static boolean finite(double value){return !Double.isNaN(value) && !Double.isInfinite(value);}
}
