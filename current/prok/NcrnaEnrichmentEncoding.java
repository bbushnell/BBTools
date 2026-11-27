package prok;

import structures.UnsignedFloat;

/** One common conversion of existing neural float-log inputs into compact
 * positive enrichment ratios and back. The per-probe weighted statistic is
 * unchanged; quantization starts from the float the legacy feature API serves.
 * A declared positive default code gives file code0an explicit meaning, so no
 * log of zero or invented missing-key probability enters the network.
 * Immutable; no per-word allocation. Metadata/defaults are training-only.
 * @author Brian, Raiden
 */
public final class NcrnaEnrichmentEncoding{

	public NcrnaEnrichmentEncoding(final int bits, final int defaultCode_){
		format=UnsignedFloat.forBits(bits);
		if(defaultCode_<=0 || !(format.decode(defaultCode_)>0)){
			throw new IllegalArgumentException("An enrichment resource needs a positive canonical decoded default ratio");
		}
		defaultCode=defaultCode_;
		defaultScore=(float)Math.log(format.decode(defaultCode));
	}

	/** Encodes exp(the exact legacy served float), retaining no sign bit. */
	public int rawCode(final float oldLogScore){return rawCode(format, oldLogScore);}

	/** Shared producer/training conversion; zero is never a rounded positive ratio. */
	public static int rawCode(final UnsignedFloat format, final float oldLogScore){
		if(!Float.isFinite(oldLogScore)){throw new IllegalArgumentException("Enrichment encoding requires a finite existing float-log input");}
		double ratio=Math.exp(oldLogScore);
		if(!(ratio>0) || !Double.isFinite(ratio)){
			throw new IllegalArgumentException("Existing log input has no positive finite ratio: "+oldLogScore);
		}
		// A decoded endpoint's float log can round just outside the double ratio
		// range. Accept that exact float representation, never an arbitrary clamp.
		if(ratio<format.minimum && oldLogScore==(float)Math.log(format.minimum)){
			ratio=format.minimum;
		}else if(ratio>format.maximum && oldLogScore==(float)Math.log(format.maximum)){
			ratio=format.maximum;
		}
		return format.encode(ratio);
	}

	/** File zero is shorthand for the explicitly declared default, not ratio zero. */
	public int fileCode(final float oldLogScore){
		final int code=rawCode(oldLogScore);
		return code==defaultCode ? 0 : code;
	}

	/** Serving and saved-vector projection feed the same decoded float-log value. */
	public float scoreCode(final int code){
		if(code==0){return defaultScore;}
		final double ratio=format.decode(code);
		assert(ratio>0) : "A nondefault enrichment code must name a positive ratio before logarithm";
		return (float)Math.log(ratio);
	}

	public float score(final float oldLogScore){return scoreCode(fileCode(oldLogScore));}
	public final UnsignedFloat format;
	public final int defaultCode;
	public final float defaultScore;
}
