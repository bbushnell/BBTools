package prot;

import java.util.Arrays;

/** Exact class-budget reduction of per-query first-rejection cutoffs.
 * The caller measures every deployed candidate and OR-reduces each query first.
 * This class performs no alignment, input sampling or positive-guided selection.
 * @author Keqing
 */
public final class FamilyLengthFit{

	/** Integer ceiling for the selected100-per-million policy, without multiplication overflow. */
	public static long allowed100ppm(final long queries){
		if(queries<=0){throw new IllegalArgumentException("A complete calibration class must contain queries: "+queries);}
		return queries/10000;
	}

	/** Finds the smallest finite amplitude meeting a complete class's explicit count budget.
	 * Sorts cutoffs[0,count) in place; the caller owns this array exclusively during
	 * fitting. No array is retained or copied. An immutable excess fails explicitly.
	 */
	public static Fit fit(final int[] cutoffs, final int count, final int maxAccepted){
		if(cutoffs==null || count<=0 || count>cutoffs.length || maxAccepted<0 || maxAccepted>count){
			throw new IllegalArgumentException("Require a nonempty complete cutoff prefix and count budget in0..count");
		}
		int immutable=0;
		for(int i=0; i<count; i++){
			final int value=cutoffs[i];
			if(value==FamilyLengthAcceptance.IMMUTABLE){immutable++;}
			else if(value<0 || value>MAX_FINITE_BITS){
				throw new IllegalArgumentException("Invalid first-rejection cutoff at query "+i+": "+value);
			}
		}
		if(immutable>maxAccepted){
			throw new IllegalStateException("The anchored length curve is infeasible: "+immutable+
				" immutable accepts exceed class budget "+maxAccepted+" among "+count+" queries");
		}
		Arrays.sort(cutoffs, 0, count);
		final int adjustableBudget=maxAccepted-immutable;
		final int threshold=maxAccepted==count ? 0 : cutoffs[count-1-adjustableBudget];
		assert(threshold>=0) : "Immutable events consume budget before selecting an adjustable quantile";
		final int accepted=countAccepted(cutoffs, count, immutable, threshold);
		final int previous=threshold==0 ? accepted : countAccepted(cutoffs, count, immutable, threshold-1);
		if(accepted>maxAccepted || (threshold>0 && previous<=maxAccepted)){
			throw new IllegalStateException("Selected cutoff must meet the budget and its preceding float must fail");
		}
		return new Fit(count, maxAccepted, immutable, threshold, accepted, previous);
	}

	/** The maximum class requirement is one shared coefficient, not separate deployed policies. */
	public static float sharedAmplitude(final Fit[] classes){
		if(classes==null || classes.length==0){throw new IllegalArgumentException("At least one fitted class is required");}
		int bits=0;
		for(int i=0; i<classes.length; i++){
			if(classes[i]==null){throw new IllegalArgumentException("Missing fitted class "+i);}
			bits=Math.max(bits, classes[i].amplitudeBits);
		}
		return FamilyLengthAcceptance.amplitude(bits);
	}

	/** Immutable events plus cutoffs strictly above the chosen amplitude accept. */
	private static int countAccepted(final int[] sorted, final int count, final int immutable, final int threshold){
		assert(immutable>=0 && immutable<=count && threshold>=0) : "Cutoff order partitions immutable and adjustable query events";
		int low=immutable, high=count;
		while(low<high){
			final int middle=low+((high-low)>>>1);
			if(sorted[middle]<=threshold){low=middle+1;}else{high=middle;}
		}
		return immutable+(count-low);
	}

	/** Small immutable final report; previousAccepted is meaningful only for positive amplitudes. */
	public static final class Fit{
		private Fit(int queries_, int budget_, int immutable_, int bits_, int accepted_, int previous_){
			queries=queries_; budget=budget_; immutable=immutable_; amplitudeBits=bits_;
			accepted=accepted_; previousAccepted=previous_;
		}
		public float amplitude(){return FamilyLengthAcceptance.amplitude(amplitudeBits);}
		public boolean hasPrevious(){return amplitudeBits>0;}
		public final int queries, budget, immutable, amplitudeBits, accepted, previousAccepted;
	}
	private FamilyLengthFit(){}
	private static final int MAX_FINITE_BITS=0x7f7fffff;
}
