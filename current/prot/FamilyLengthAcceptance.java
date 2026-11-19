package prot;

/** Shared overlap-length penalty for the reviewed LONG/SHORT16K5 study.
 * Inputs are already measured on one alignment and include the frozen effective
 * family gates. This class neither loads classifiers nor chooses an error budget.
 * @author Keqing
 */
public final class FamilyLengthAcceptance{

	/** Returns the same alignment-route codes as FamilyTailCascade:2 identity, 3 BLOSUM.
	 * Both spans must reach28. The BLOSUM branch retains its inclusive identity
	 * lower comparison; both upper comparisons remain strict. No gate is clamped.
	 */
	public static int route(final float identity, final float blosum,
			final double identityLower, final double identityUpper, final double blosumUpper,
			final int querySpan, final int referenceSpan, final float amplitude){
		final int length=overlap(querySpan, referenceSpan);
		validateScores(identity, blosum, identityLower, identityUpper, blosumUpper);
		validateAmplitude(amplitude);
		if(length<MIN_OVERLAP){return REJECT;}
		return routeMeasured(identity, blosum, identityLower, identityUpper, blosumUpper, penalty(querySpan, referenceSpan, amplitude));
	}

	/** Penalty added to the identity gates; BLOSUM receives BLOSUM_UNIT times this value. */
	public static double penalty(final int querySpan, final int referenceSpan, final float amplitude){
		final int length=overlap(querySpan, referenceSpan);
		validateAmplitude(amplitude);
		return length<MIN_OVERLAP ? Double.POSITIVE_INFINITY : (double)amplitude*((MAX_OVERLAP-length)/(double)length);
	}

	/** First nonnegative finite-float amplitude rejecting this candidate, encoded
	 * as its positive IEEE bits. Zero means already rejected. IMMUTABLE means an
	 * accepted100-residue overlap that this anchored curve cannot change.
	 * The search uses the deployed predicate, including its double rounding and ties.
	 */
	public static int firstRejectionBits(final float identity, final float blosum,
			final double identityLower, final double identityUpper, final double blosumUpper,
			final int querySpan, final int referenceSpan){
		final int length=overlap(querySpan, referenceSpan);
		validateScores(identity, blosum, identityLower, identityUpper, blosumUpper);
		if(length<MIN_OVERLAP || routeMeasured(identity, blosum, identityLower, identityUpper, blosumUpper, 0)==REJECT){return 0;}
		if(length==MAX_OVERLAP){return IMMUTABLE;}
		final double factor=(MAX_OVERLAP-length)/(double)length;
		if(routeMeasured(identity, blosum, identityLower, identityUpper, blosumUpper, (double)Float.MAX_VALUE*factor)!=REJECT){
			throw new IllegalArgumentException("Finite-float amplitude cannot reject adjustable candidate; identity="+identity+" overlap="+length);
		}
		int low=0, high=MAX_FINITE_BITS;
		while(low<high){
			final int middle=low+((high-low)>>>1);
			final double penalty=(double)Float.intBitsToFloat(middle)*factor;
			if(routeMeasured(identity, blosum, identityLower, identityUpper, blosumUpper, penalty)==REJECT){high=middle;}
			else{low=middle+1;}
		}
		assert(low>0 && routeMeasured(identity, blosum, identityLower, identityUpper, blosumUpper,
			(double)Float.intBitsToFloat(low)*factor)==REJECT &&
			routeMeasured(identity, blosum, identityLower, identityUpper, blosumUpper,
			(double)Float.intBitsToFloat(low-1)*factor)!=REJECT) :
			"Query cutoffs require the first rejecting float, not an arbitrary sufficient threshold";
		return low;
	}

	/** OR-reduces candidate cutoffs for one query without allocating per-pair state.
	 * An unchanged shortcut is represented explicitly by IMMUTABLE as well.
	 */
	public static int merge(final int left, final int right){
		validateCutoff(left);validateCutoff(right);
		return left==IMMUTABLE || right==IMMUTABLE ? IMMUTABLE : Math.max(left, right);
	}

	/** Evaluates a stored query cutoff. Signed zero amplitudes both mean zero. */
	public static boolean accepts(final int cutoffBits, final float amplitude){
		validateCutoff(cutoffBits);validateAmplitude(amplitude);
		final int bits=amplitude==0 ? 0 : Float.floatToRawIntBits(amplitude);
		return cutoffBits==IMMUTABLE || bits<cutoffBits;
	}

	/** Decodes a finite fitted threshold; an immutable event is not an infinity value. */
	public static float amplitude(final int cutoffBits){
		validateCutoff(cutoffBits);
		require(cutoffBits!=IMMUTABLE, "An immutable accept has no finite rejecting amplitude");
		return Float.intBitsToFloat(cutoffBits);
	}

	/** Uses the shorter aligned span; a long span cannot compensate for a short one.
	 * A consensus can exceed the query tail cap, and an indel can make its aligned
	 * span exceed100. The curve's upper bound applies to the minimum, not each span. */
	public static int overlap(final int querySpan, final int referenceSpan){
		final int length=Math.min(querySpan, referenceSpan);
		if(length<=0 || length>MAX_OVERLAP){
			throw new IllegalArgumentException("Measured spans must be positive and their minimum must be1..100: query="+querySpan+" reference="+referenceSpan);
		}
		return length;
	}

	/** Fixed arithmetic order shared by direct scoring and the cutoff search. */
	private static int routeMeasured(final float identity, final float blosum,
			final double lower, final double upper, final double finalBlosum, final double penalty){
		assert(penalty>=0 && Double.isFinite(penalty)) : "Finite nonnegative length penalty preserves monotone acceptance";
		if(identity>upper+penalty){return IDENTITY;}
		return identity>=lower+penalty && blosum>finalBlosum+BLOSUM_UNIT*penalty ? BLOSUM : REJECT;
	}

	/** Native identities and effective gates must remain valid before numerical fitting. */
	private static void validateScores(final float identity, final float blosum,
			final double lower, final double upper, final double finalBlosum){
		if(!Float.isFinite(identity) || identity<0 || identity>100 || !Float.isFinite(blosum)){
			throw new IllegalArgumentException("Length fitting needs measured finite scores, not missing-value sentinels: identity="+identity+" blosum="+blosum);
		}
		if(!Double.isFinite(lower) || lower<0 || lower>100 || !Double.isFinite(upper) || upper<lower || !Double.isFinite(finalBlosum)){
			throw new IllegalArgumentException("Effective family gates must retain finite ordered identity bounds: lower="+lower+" upper="+upper+" blosum="+finalBlosum);
		}
	}

	/** Amplitudes are finite float parameters; no negative or infinity policy is implicit. */
	public static void validateAmplitude(final float value){
		if(!Float.isFinite(value) || value<0){throw new IllegalArgumentException("Length amplitude must be finite and nonnegative: "+value);}
	}

	/** Positive finite float bits preserve numeric order; the sole sentinel is explicit. */
	private static void validateCutoff(final int bits){
		if(bits!=IMMUTABLE && (bits<0 || bits>MAX_FINITE_BITS)){throw new IllegalArgumentException("Invalid encoded length cutoff: "+bits);}
	}
	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}
	private FamilyLengthAcceptance(){}

	public static final int MIN_OVERLAP=28, MAX_OVERLAP=100;
	public static final int REJECT=0, IDENTITY=2, BLOSUM=3;
	public static final int IMMUTABLE=-1;
	public static final double BLOSUM_UNIT=.005;
	private static final int MAX_FINITE_BITS=0x7f7fffff;
}
