package aligner;

import idaligner.AlignmentStats;
import idaligner.EndClippedAligner;

/** Once-per-consensus CM map composed with an existing caller alignment.
 * Does not align a candidate or consult its exact CM score. A map is built
 * from a separately bound global CM alignment of the caller's consensus.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelConsensusMap {
	public CovarianceModelConsensusMap(CovarianceModel model_, byte[] consensus, byte[] aligned, byte[] rf){
		final CovarianceModelPlacementBands map=new CovarianceModelPlacementBands(model_, consensus, aligned, rf, 0);
		model=model_;consensusLength=consensus.length;low=map.low.clone();high=map.high.clone();
	}

	/** Query=full caller consensus, reference=oriented candidate window.
	 * Window/candidate bounds are zero-based on the same oriented chromosome.
	 * Null means no usable anchors: missing trace or a genuinely clipped query.
	 * Malformed traces throw; they must never become plausible band coordinates. */
	public CovarianceModelPlacementBands compose(AlignmentStats stats, int windowStart,
			int candidateStart, int candidateStop, int radius){
		return compose(stats, windowStart, candidateStart, candidateStop, radius, false);
	}
	/** Experimental D34: unaligned query ends describe intervals, not fabricated
	 * point anchors. Extend toward the retained window edges, then intersect with
	 * the unchanged candidate coordinate domain before adding the usual radius. */
	public CovarianceModelPlacementBands compose(AlignmentStats stats, int windowStart,
			int candidateStart, int candidateStop, int radius, boolean extendClipped){
		require(windowStart>=0 && candidateStart>=0 && candidateStop>=candidateStart
			&& (long)candidateStop-candidateStart<Integer.MAX_VALUE && radius>=0,
			"Anchor composition requires inclusive candidate bounds on one oriented chromosome");
		if(stats==null || stats.matchString==null){return null;}
		require(stats.qLen==consensusLength && stats.rStart>=0 && stats.rStop>=stats.rStart && stats.rStop<stats.rLen,
			"Caller trace must bind the full consensus query and an actual reference interval");
		int queryStart=0, queryEnd=consensusLength;
		if(stats instanceof EndClippedAligner.Result){
			final EndClippedAligner.Result clipped=(EndClippedAligner.Result)stats;
			require(clipped.qStart>=0 && clipped.qStop>=clipped.qStart && clipped.qStop<consensusLength,
				"Clipped query bounds must remain within the mapped consensus");
			queryStart=clipped.qStart;queryEnd=clipped.qStop+1;
			if(!extendClipped && (queryStart!=0 || queryEnd!=consensusLength)){return null;}
		}
		final int[] queryLow=new int[consensusLength+1], queryHigh=new int[consensusLength+1];
		// Every boundary strictly before/after the retained trace can lie
		// anywhere in the corresponding unanchored reference prefix/suffix.
		for(int i=0; i<queryStart; i++){queryHigh[i]=stats.rStart;}
		for(int i=queryEnd+1; i<=consensusLength; i++){
			queryLow[i]=stats.rStop+1;queryHigh[i]=stats.rLen;
		}
		int q=queryStart, r=stats.rStart;queryLow[q]=queryHigh[q]=r;
		// AlignmentStats.setFromMatchString and QuantumAligner: I consumes query
		// only, D consumes reference only; m/S/N consume both. These are long ops.
		for(byte op:stats.matchString){
			if(op=='D'){require(r<=stats.rStop, "Trace deletion exceeds its reference endpoint");queryHigh[q]=++r;}
			else{
				require(op=='m' || op=='S' || op=='N' || op=='I', "Unsupported caller traceback operation: "+(char)op);
				require(q<queryEnd, "Trace consumes query bases beyond the declared retained alignment");
				if(op!='I'){require(r<=stats.rStop, "Trace diagonal exceeds its reference endpoint");r++;}
				q++;queryLow[q]=queryHigh[q]=r;
			}
		}
		require(q==queryEnd && r==stats.rStop+1,
			"Trace must conserve its declared query/reference spans before missing ends become intervals");
		final int length=candidateStop-candidateStart+1;
		final long shift=(long)windowStart-candidateStart;
		final int[] targetLow=new int[low.length], targetHigh=new int[high.length];
		for(int m=0; m<low.length; m++){
			targetLow[m]=clip(queryLow[low[m]]+shift, length);
			targetHigh[m]=clip(queryHigh[high[m]]+shift, length);
		}
		// Refined ends may extend beyond the retained alignment. Keep terminal
		// insertions representable; interior anchors still come only from the trace.
		targetLow[0]=0;targetHigh[model.clen]=length;
		return new CovarianceModelPlacementBands(model, length, targetLow, targetHigh, radius);
	}
	private static int clip(long position, int length){
		assert(length>=0):"Clipping uses the exact nonnegative refined-candidate length";
		return (int)Math.max(0L, Math.min(length, position));
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	boolean belongsTo(CovarianceModel expected){return model==expected;}
	/** Number of CM match positions; positions0..modelLength are boundaries. */
	public int modelLength(){return model.clen;}
	/** Inclusive range of zero-based consensus boundaries corresponding to a CM
	 * boundary. Inserted residues widen it; deleted model columns may share it.
	 * These are frozen alignment coordinates, not learned endpoint means. */
	public int boundaryLow(int boundary){
		require(boundary>=0 && boundary<low.length, "CM boundary must belong to this bound global consensus map");return low[boundary];
	}
	public int boundaryHigh(int boundary){
		require(boundary>=0 && boundary<high.length, "CM boundary must belong to this bound global consensus map");return high[boundary];
	}
	public final int consensusLength;
	private final CovarianceModel model;
	private final int[] low, high;
}
