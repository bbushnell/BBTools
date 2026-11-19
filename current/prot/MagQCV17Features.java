package prot;

/** Observation-only additions for the v1.7 training and serving layouts.
 * @author Brian Bushnell, Yoimiya
 */
final class MagQCV17Features {

	/** Counts tracked protein copies over all recorded calls, including existing ROTHER calls.
	 * Unknown anticodons are already included in total tRNA. This is the reverse of
	 * the v1.4 called-to-tracked operands, not the composite's older mapped/CDS input.
	 * No truth labels, taxonomy reassignment, clipping or inferred missing genes enter it.
	 * The legacy prepared contract does not imply sum(families)==mapped or mapped<=CDS;
	 * therefore this method does not invent either invariant or clamp the result.
	 */
	static double trackedToCalled(final MagQCPreparedBin bin){
		if(bin==null){throw new IllegalArgumentException("Tracked-call input requires observed bin counts");}
		final MagQCVectorMaker.Agg a=bin.agg;
		if(a.cds<0 || a.r16<0 || a.r23<0 || a.r5<0 || a.rother<0 || a.trna<0){
			throw new IllegalArgumentException("Negative called-gene count for "+bin.id);
		}
		final long calls=(long)a.cds+a.r16+a.r23+a.r5+a.rother+a.trna;
		long tracked=0;
		for(int count:bin.families){
			if(count<0){throw new IllegalArgumentException("Negative tracked-family count for "+bin.id);}
			tracked+=count;
		}
		final double value=calls==0 ? 0 : tracked/(double)calls;
		assert(Double.isFinite(value) && value>=0) : "Nonnegative integer observations and the explicit zero fallback must produce a finite ratio";
		return value;
	}

	static final String COLUMN="tracked_to_called";
	static final String DEFINITION="sum_family_copies/(cds+r16+r23+r5+rother+trna);zero_denominator=0;unclipped";
}
