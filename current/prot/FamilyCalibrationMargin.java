package prot;

/** Converts an effective-margin increase back to the stored family upper gate.
 * Member lower gates and the global floor stay fixed. No family or policy state
 * is mutated; callers retain the exact margin used to observe the winning score.
 * @author Keqing
 */
public final class FamilyCalibrationMargin{

	/** Raises a stored gate enough to put the effective gate at least increment
	 * above a measured score. The classifier uses a strict score>gate comparison.
	 * Scale1 retains the historical upward-rounded score-plus-increment rule.
	 * Infeasible finite-float output fails rather than becoming an infinity sentinel.
	 */
	public static float raiseUpper(final float lower,final float upper,final float floor,
			final double scale,final float score,final double increment){
		require(Double.isFinite(scale) && scale>0,
			"A positive finite margin is required to invert an effective acceptance gate");
		require(Float.isFinite(score) && Double.isFinite(increment) && increment>0,
			"Gate repair requires a finite measured score and positive finite increment");
		final double current=FamilyTailCascade.evaluationGate(lower,upper,floor,scale);
		final double target=(double)score+increment;
		require(Double.isFinite(target),"Measured score plus repair increment overflowed");
		if(target<=current){return upper;}
		final double base=Math.max(lower,floor);
		final double inverse=scale==1 ? target : base+(target-base)/scale;
		float result=FamilyGateCodec.roundUpper(Math.max((double)upper,inverse));
		// Inverting in double precision and publishing as float can meet a rounding
		// boundary. Check the actual deployed arithmetic, not just the inverse formula.
		if(FamilyTailCascade.evaluationGate(lower,result,floor,scale)<target){result=Math.nextUp(result);}
		require(Float.isFinite(result) && result>=upper &&
			FamilyTailCascade.evaluationGate(lower,result,floor,scale)>=target,
			"Stored float gate cannot represent the requested effective-margin increase");
		return result;
	}

	private static void require(final boolean valid,final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}
	private FamilyCalibrationMargin(){}
}
