package prok;

/** Shared raw-span inputs and proposal rule for the full-error movement gate.
 * Callers own every array. No sequence, score or feature buffer is modified
 * except the explicit result array; methods allocate no per-proposal objects.
 * @author Brian, Raiden
 */
public final class NcrnaEndpointGateFeatures{

	/*--------------------------------------------------------------*/
	/*----------------       Proposal and Inputs    ----------------*/
	/*--------------------------------------------------------------*/

	/** A reachable tie or unchanged-position winner holds; head21 is only an input. */
	public static int proposal(final float[] scores){
		require(scores!=null && scores.length==22, "Movement proposal requires all22first-stage scores");
		int best=0;
		for(int i=0; i<22; i++){
			require(Float.isFinite(scores[i]), "Movement scores must be finite");
			if(i<21 && scores[i]>scores[best]){best=i;}
		}
		if(best==10){return -1;}
		for(int i=0; i<21; i++){if(i!=best && scores[i]==scores[best]){return -1;}}
		return best;
	}

	/** Length ratios use raw/proposed spans; GC, identity and fuzz stay on the raw span. */
	public static void fill(final float[] scores, final float[] raw, final float spanGC, final float refGC,
			final long proposedLength, final int referenceLength, final float[] result){
		require(scores!=null && scores.length==22 && raw!=null && raw.length==69 && result!=null && result.length==32
			&& referenceLength>0 && spanGC>=0 && spanGC<=1 && refGC>=0 && refGC<=1,
			"Gate layout is22scores,2lengths,4GC,identity,3fuzz; invalid proposed length is retained for audit");
		System.arraycopy(scores, 0, result, 0, 22);
		result[22]=raw[63];
		result[23]=proposedLength/(float)referenceLength;
		result[24]=raw[64];
		result[25]=spanGC;
		result[26]=refGC;
		result[27]=spanGC-refGC;
		result[28]=2f*raw[65]-1f;
		System.arraycopy(raw, 66, result, 29, 3);
		for(final float value:result){require(Float.isFinite(value), "Gate features cannot contain nonfinite values");}
	}

	/** GC/(A+C+G+T), preserving the training reader's recognized-base denominator. */
	public static float gc(final byte[] sequence, final int start, final int stop){
		require(sequence!=null && start>=0 && stop>=start && stop<sequence.length, "GC needs a real inclusive sequence span");
		int recognized=0;
		int gc=0;
		for(int i=start; i<=stop; i++){
			final int base=RrnaResourceIO.base(sequence[i]);
			if(base>=0){recognized++;}
			if(base==1 || base==2){gc++;}
		}
		require(recognized>0, "GC is unavailable for an all-ambiguous span; never invent a zero fraction");
		return gc/(float)recognized;
	}

	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private NcrnaEndpointGateFeatures(){}
	public static final String FEATURE_ORDER="scores0..21,rawLengthRatio,proposedLengthRatio,contigGC,rawSpanGC,referenceGC,spanMinusReferenceGC,rescaledRawIdentity,fuzzMinus,fuzzAt,fuzzPlus";
}
