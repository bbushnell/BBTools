package prot;

/** Immutable serving-only policy for the capped subnet experiment. Training targets,
 * weights, arithmetic features and feature units are unchanged.
 * @author Brian Bushnell, Yoimiya
 */
public final class MagQCSubnetOutputCaps {

	/** Ordinary predictions have already been blended/decoded, but errors are still fractions. */
	public static MagQCSubnetOutputCaps ordinary(final MagQCPredictionContract prediction){
		if(prediction==null || !prediction.hasCanonicalSubnetEndpointBounds()){
			throw new IllegalArgumentException("Subnet caps require the six-head 0.1–100 percent endpoint decoder");
		}
		return ORDINARY;
	}

	/** Size slots1/6 are unbounded linear size and its residual; all others target [0,1]. */
	public static MagQCSubnetOutputCaps size(){return SIZE;}

	private MagQCSubnetOutputCaps(final boolean size_){size=size_; count=size ? 10 : 6;}

	/** Caps only model slots, retaining any ordinary blend diagnostics in the tail.
	 * Endpoint decoding is d=clamp(10^(3*f-1),0.1,100) percent. Re-encoding that
	 * bounded d gives clamp(f,0,1) exactly. Use the equivalent monotone clamp to
	 * avoid pow overflow and preserve every in-range float bit, rather than adding
	 * a lossy pow/log round trip. The ordinary factory validates those constants.
	 */
	public void apply(final float[] values){
		if(values==null || values.length<count){throw new IllegalArgumentException("Subnet cap scratch is too short");}
		for(int i=0; i<count; i++){
			final float value=values[i];
			if(!Float.isFinite(value)){throw new IllegalArgumentException("Nonfinite subnet prediction at head "+i);}
			if(!size || (i!=1 && i!=6)){values[i]=value<0 ? 0 : value>1 ? 1 : value;}
		}
		assert(!size || count==10) : "Size cap exclusions1/6 follow the ten-head MagQCSizeFeatures target order";
	}

	private final boolean size;
	private final int count;
	private static final MagQCSubnetOutputCaps ORDINARY=new MagQCSubnetOutputCaps(false);
	private static final MagQCSubnetOutputCaps SIZE=new MagQCSubnetOutputCaps(true);
}
