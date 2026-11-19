package prok;

/** Shared opt-in confidence rule for scalar and positional endpoint scores.
 * The legacy selector supplies its own historical choice, including ties.
 * Its zero margin returns that exact choice. A positive margin permits it only
 * when it exceeds both the current position and every other available position.
 * The explicit full-error route also compares its non-positional head and
 * requires strict gaps even at zero margin.
 * No scores, coordinates or buffers are modified; all state belongs to callers.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaEndpointConfidence{

	/** Selects an index in an ordered score window, retaining current on doubt.
	 * Every supplied position is available and must have a finite computed score.
	 * Zero margin returns the historical choice; an unchanged historical choice
	 * holds. Missing physical positions use the explicit-mask overload.
	 * Both current and legacy must be available. Equal best positions prevent a
	 * positive-margin move, regardless of the historical selector's tie rule.
	 * @param scores Borrowed per-position network scores, not probabilities
	 * @param current Index of the unchanged caller end
	 * @param legacy Exact index chosen by the existing ungated selector
	 * @param margin Nonnegative finite absolute difference in network-score units
	 * @return Legacy index for zero margin; legacy or current for positive margin
	 */
	public static int select(final float[] scores, final int current, final int legacy, final float margin){
		return select(scores, current, legacy, margin, null);
	}

	/** A null mask means all positions exist. Otherwise false positions are never
	 * read as measurements; no finite signed value is reserved for missingness.
	 * Every available score is validated even when the chosen end will not move.
	 */
	public static int select(final float[] scores, final int current, final int legacy, final float margin,
			final boolean[] available){
		if(scores==null || current<0 || current>=scores.length || legacy<0 || legacy>=scores.length
				|| available!=null && (available.length!=scores.length || !available[current] || !available[legacy])
				|| !Float.isFinite(margin) || margin<0 || !Float.isFinite(scores[current]) || !Float.isFinite(scores[legacy])){
			throw new IllegalArgumentException("Endpoint confidence needs available current/legacy positions and a finite nonnegative score margin");
		}
		float runnerUp=-Float.MAX_VALUE;
		for(int i=0; i<scores.length; i++){
			if(available!=null && !available[i]){continue;}
			final float score=scores[i];
			if(!Float.isFinite(score)){
				throw new IllegalArgumentException("Nonfinite computed endpoint score cannot establish a confidence gap: index="+i);
			}
			if(i!=legacy && score>runnerUp){runnerUp=score;}
		}
		if(margin==0){return legacy;}
		if(legacy==current){return current;}
		assert(Float.isFinite(runnerUp)) : "A proposed move has a distinct finite current position, so runner-up must exist";
		return permits(scores[current], scores[legacy], runnerUp, margin) ? legacy : current;
	}

	/** Selects a reachable offset with an additional non-positional out-of-range head.
	 * The first 2*radius+1 scores are ordered offsets; the last score only vetoes
	 * movement. Even at zero margin, ties with current, another offset or that
	 * head hold the end. All scores must be finite; fixed or unavailable ends
	 * are held by the caller before inference. The borrowed buffer is unchanged.
	 * @param scores Raw scores, including the last non-positional head
	 * @param radius Radius in bases, also the unchanged-end output index
	 * @param margin Finite nonnegative gap in raw network-score units
	 * @return A reachable index, never the out-of-range head
	 */
	public static int selectFull(final float[] scores, final int radius, final float margin){
		return selectFull(scores, radius, margin, margin);
	}

	/** Selects a unique reachable winner with independently controlled head rejection.
	 * Current and runner-up require strict positional gaps. Only the comparison
	 * against the last, non-positional head accepts a finite signed margin.
	 * All input scores must still be finite.
	 * @param scores Borrowed raw offset scores followed by one out-of-range score
	 * @param radius Center index and correction radius in bases
	 * @param positionalMargin Finite nonnegative gap against current and runner-up
	 * @param oorMargin Finite signed gap against the last head
	 * @return A reachable index; current on a positional tie or insufficient gap
	 */
	public static int selectFull(final float[] scores, final int radius, final float positionalMargin, final float oorMargin){
		return selectFull(scores, radius, positionalMargin, oorMargin, true);
	}

	/** An ignored head cannot veto, but its produced score must still be finite.
	 * The signed numeric margin is finite in both modes and unused when ignored.
	 * Positional ties and strict current/runner-up gaps retain the same rules.
	 */
	public static int selectFull(final float[] scores, final int radius, final float positionalMargin,
			final float oorMargin, final boolean useOorHead){
		if(scores==null || radius<1 || (long)radius*2+2!=scores.length
				|| !Float.isFinite(positionalMargin) || positionalMargin<0
				|| !Float.isFinite(oorMargin)){
			throw new IllegalArgumentException("Full-error confidence requires a centered window, a finite nonnegative positional margin and a finite signed head margin");
		}
		final int positions=scores.length-1;
		int best=0;
		for(int i=0; i<scores.length; i++){
			if(!Float.isFinite(scores[i])){
				throw new IllegalArgumentException("Nonfinite full-error score cannot establish a move: index="+i);
			}
			if(i<positions && scores[i]>scores[best]){best=i;}
		}
		if(best==radius){return radius;}
		float runner=-Float.MAX_VALUE;
		for(int i=0; i<positions; i++){
			if(i!=best){runner=Math.max(runner, scores[i]);}
		}
		final double top=scores[best];
		final int selected=top>(double)scores[radius]+positionalMargin && top>(double)runner+positionalMargin
			&& (!useOorHead || top>(double)scores[positions]+oorMargin) ? best : radius;
		assert(selected<positions) : "The last full-error output is a veto, never a coordinate correction";
		return selected;
	}

	/** Same positive-margin rule for a scalar sweep retaining only its top two scores.
	 * Callers preserve their historical zero-margin selector before using this gate.
	 */
	public static boolean permits(final float current, final float best, final float runnerUp, final float margin){
		if(!Float.isFinite(current) || !Float.isFinite(best) || !Float.isFinite(runnerUp)
				|| !Float.isFinite(margin) || margin<0){
			throw new IllegalArgumentException("Endpoint confidence requires finite scores and a nonnegative finite margin");
		}
		return best-current>margin && best-runnerUp>margin;
	}

	private NcrnaEndpointConfidence(){}
}
