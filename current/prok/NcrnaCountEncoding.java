package prok;

import structures.CountFloat;

/** Shared count-to-score mapping for compact endpoint resources and training.
 * Anchor counts are integers. Background values are nonnegative effective counts
 * on the anchor denominator, stored with a declared power-of-two fractional
 * scale. Pseudocount background remains exact metadata, separate from observed
 * counts. Sparse mode removes anchor counts0and1before quantization.
 * All fields are immutable; per-locus LOO quantities are supplied by the caller.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCountEncoding{

	/** Selects a training-only background scale, leaving a factor of two headroom.
	 * Both count codes decode into checked longs; doubles participate only in the
	 * probability/score arithmetic, never as a fallback for an overflowing code.
	 */
	public static NcrnaCountEncoding create(final int bits, final boolean sparse,
			final double alpha, final double largestBackground){
		final CountFloat format=CountFloat.forBits(bits);
		require(Double.isFinite(largestBackground) && largestBackground>=0,
			"Background scaling must use the observed finite nonnegative training-count range");
		final double capacity=format.maximum;
		require(largestBackground<=capacity/2,
			"Observed background counts need a larger count format before a scale can leave training-view headroom");
		final int scale=largestBackground==0 ? 0 : Math.max(0,
			Math.min(52, Math.getExponent(capacity/(2*largestBackground))));
		return new NcrnaCountEncoding(bits, sparse, alpha, scale);
	}

	/** Reconstructs a declared resource mapping without learning from its inputs. */
	public NcrnaCountEncoding(final int bits, final boolean sparse_, final double alpha_, final int scalePower_){
		format=CountFloat.forBits(bits);
		require(Double.isFinite(alpha_) && alpha_>0 && scalePower_>=0 && scalePower_<=52,
			"Compact endpoint metadata requires positive finite smoothing and an explicit bounded binary background scale");
		sparse=sparse_;
		alpha=alpha_;
		scalePower=scalePower_;
	}

	/** Sparse evidence is true-anchor occurrences in this actual training view. */
	public boolean retained(final long anchorCount){
		require(anchorCount>=0, "Anchor support is an unsmoothed observed count, never negative");
		return !sparse || anchorCount>1;
	}

	/** Encodes an effective background count as a nonnegative scaled count. */
	public int backgroundCode(final double count){
		require(Double.isFinite(count) && count>=0, "Observed weighted background must be a finite nonnegative effective count");
		return format.encode(Math.scalb(count, scalePower));
	}

	/** Exact decoding of the declared fixed-point unit; overflow is checked first. */
	public double backgroundCount(final int code){
		return Math.scalb((double)format.decode(code), -scalePower);
	}

	/** Training and producer use the same native count decoder as resource loading.
	 * LOO removes the locus before passing its counts/prior here; held-out and
	 * serving views supply full-training counts. No held-out values choose scale.
	 */
	public float score(final long anchorCount, final double backgroundCount, final double backgroundPrior){
		require(anchorCount>=0 && Double.isFinite(backgroundCount) && backgroundCount>=0,
			"Quantize actual nonnegative remaining-population counts before deriving an endpoint score");
		if(!retained(anchorCount)){return defaultScore(backgroundPrior);}
		return scoreCodes(format.encode((double)anchorCount), backgroundCode(backgroundCount), backgroundPrior);
	}

	/** Serving decodes precisely the two counts used to define a training score. */
	public float scoreCodes(final int anchorCode, final int backgroundCode, final double backgroundPrior){
		final long anchor=format.decode(anchorCode);
		require(!sparse || anchor>1, "A sparse endpoint record must retain more than one true-anchor training occurrence");
		return scoreCounts(anchor, backgroundCount(backgroundCode), backgroundPrior);
	}

	/** Omitted rows have zero observed anchor and background, with exact smoothing. */
	public float defaultScore(final double backgroundPrior){return scoreCounts(0, 0, backgroundPrior);}

	private float scoreCounts(final long anchor, final double background, final double prior){
		require(Double.isFinite(prior) && prior>0,
			"Smoothed wrong-offset background must remain positive even for an omitted endpoint word");
		final double numerator=anchor+alpha, denominator=background+prior;
		final float result=(float)Math.log(numerator/denominator);
		require(Float.isFinite(result), "Decoded positive smoothed counts must yield a finite neural endpoint input");
		return result;
	}

	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}

	public final CountFloat format;
	public final boolean sparse;
	public final double alpha;
	public final int scalePower;
}
