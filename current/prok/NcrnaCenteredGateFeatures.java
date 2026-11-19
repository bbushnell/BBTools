package prok;

/** Separately versioned centered gate: accept/reject the WHOLE raw-to-proposal
 * move. The original 32-input uncentered layout and its callers stay unchanged.
 * All buffers belong to the invoking worker; no per-proposal allocation.
 * @author Brian, Ganyu
 */
public final class NcrnaCenteredGateFeatures{
	/*--------------------------------------------------------------*/
	/*----------------       Proposal and Inputs    ----------------*/
	/*--------------------------------------------------------------*/

	/** Current best or a reachable tie means no neural adjustment, while a
	 * nonzero median displacement can still be a whole-move gate candidate. */
	public static int adjustment(final float[] scores){
		final int proposed=NcrnaEndpointGateFeatures.proposal(scores);
		return proposed<0 ? 0 : proposed-10;
	}

	/** Both ends' features use the same guarded center. The caller must reject
	 * fixed/invalid-center/missing-feature cases before requesting gate inference.
	 * A single-end proposed span changes only that end of C, not the other end.
	 * Returns the candidate coordinate; equality with the raw end is a hold.
	 */
	public static long fill(final float[] scores, final float[] centeredFeatures,
			final float spanGC, final float referenceGC, final int referenceLength,
			final int rawStart, final int rawStop, final int centerStart, final int centerStop,
			final int end, final float[] scratch32, final float[] result35){
		require(rawStart>=0 && rawStop>=rawStart && centerStart>=0 && centerStop>=centerStart
			&& (end==0 || end==1) && scratch32!=null && scratch32.length==32
			&& result35!=null && result35.length==INPUTS && centeredFeatures!=null && centeredFeatures.length==69
			&& referenceLength>0 && (long)centerStop-centerStart+1<=Integer.MAX_VALUE,
			"Centered whole-move inputs need valid raw/center pairs, an explicit end and separate32/35worker buffers");
		require(Float.floatToRawIntBits(centeredFeatures[63])==Float.floatToRawIntBits((centerStop-centerStart+1)/(float)referenceLength),
			"NcrnaPositionalBoundaryFeatures.fill supplies center length/reference length at63; a raw-span vector cannot teach or serve a centered gate");
		final int neural=adjustment(scores);
		final long proposedLength=(long)centerStop-centerStart+1+(end==0 ? -neural : neural);
		NcrnaEndpointGateFeatures.fill(scores, centeredFeatures, spanGC, referenceGC, proposedLength, referenceLength, scratch32);
		System.arraycopy(scratch32, 0, result35, 0, 32);
		final int raw=end==0 ? rawStart : rawStop;
		final int center=end==0 ? centerStart : centerStop;
		final long applied=(long)center-raw;
		result35[32]=Math.abs(applied);
		result35[33]=applied;
		result35[34]=neural==0 ? 0 : 1;
		for(final float value:result35){require(Float.isFinite(value), "Centered gate inputs must never invent a missing or nonfinite feature");}
		return (long)center+neural;
	}

	/** Teaching target compares against R, not C. Equal distance is negative;
	 * a correction outside the first-stage radius still supplies a real label. */
	public static boolean improves(final int raw, final long proposal, final int truth){
		require(raw>=0 && truth>=0 && proposal>=-10 && proposal<=(long)Integer.MAX_VALUE+10,
			"Whole-move targets use physical int coordinates plus a bounded neural adjustment");
		return Math.abs(proposal-truth)<Math.abs((long)raw-truth);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Fields / Constants     ----------------*/
	/*--------------------------------------------------------------*/

	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private NcrnaCenteredGateFeatures(){}
	public static final int INPUTS=35;
	public static final String SCHEMA="CENTERED_WHOLE_MOVE_35_V1";
	public static final String FEATURE_ORDER="scores0..21,centerLengthRatio,proposedLengthRatio,contigGC,centerSpanGC,referenceGC,spanMinusReferenceGC,rescaledCenterIdentity,fuzzMinus,fuzzAt,fuzzPlus,absoluteAppliedMedianDisplacement,signedAppliedMedianDisplacement,neuralMoved";
}
