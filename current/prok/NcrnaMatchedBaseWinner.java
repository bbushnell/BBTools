package prok;

import idaligner.AlignmentStats;

/** Comparison evidence for the explicit shared-clan experiment, not a new
 * alignment score. Preserve the existing aligners' numerator, including their
 * half contribution for ambiguous bases. Reference span excludes insertion
 * columns, so identity times reference span is not this quantity.
 * @author Ganyu
 */
public final class NcrnaMatchedBaseWinner{
	public static final int IDENTITY=0, MATCHES=1, CIGAR_SCORE=2, MODEL_SCORE=3;

	/** Reproduce Tracer's own denominator, including its substitution clamp.
	 * This normally equals queryLength+deletions; never derive it from rounded id.
	 * Tracer's integer reported M/S fields can truncate two half-valued counts.
	 */
	public static float columns(final int[] positions, final int queryLength){
		if(positions==null || positions.length<7 || positions[1]<positions[0]){
			throw new IllegalArgumentException("The experimental packed observer requires the seven actual Quantum/Scrabble result fields");
		}
		final float matches=numerator(positions[2], queryLength, positions[3]);
		final int insertions=Math.max(0, queryLength+positions[3]-(positions[1]-positions[0]+1));
		final float substitutions=Math.max(0, queryLength-matches-insertions);
		if(positions[4]!=insertions || positions[5]!=(int)matches || positions[6]!=(int)substitutions){
			throw new IllegalArgumentException("Packed MSDI fields must agree with the same Tracer postprocess equation");
		}
		final float exact=matches+substitutions+insertions+positions[3];
		final long reported=reportedColumns(positions);
		return reported==exact ? reported : exact;
	}

	/** Literal integer MSDI fields, retained separately to expose N truncation. */
	public static long reportedColumns(final int[] positions){
		assert(positions!=null && positions.length>=7) : "Reported MSDI requires the original seven-result observation";
		return (long)positions[5]+positions[6]+positions[3]+positions[4];
	}

	public static float columns(final AlignmentStats stats){
		if(stats==null || stats.matches<0 || stats.subs<0 || stats.ins<0 || stats.dels<0 || stats.ns<0){
			throw new IllegalArgumentException("A traced alignment needs its complete nonnegative column counts, including N columns");
		}
		return (long)stats.matches+stats.subs+stats.ins+stats.dels+stats.ns;
	}

	/** B/B2 use the returned verifier identity, promoted to double for arithmetic.
	 * M retains the existing numerator. No threshold or patience changes here. */
	public static double value(final int rule, final float matches, final float identity,
			final float columns, final int modelLength){
		return value(rule, matches, identity, columns, modelLength, 0.54);
	}
	public static double value(final int rule, final float matches, final float identity,
			final float columns, final int modelLength, final double offset){
		assert(Float.isFinite(matches) && matches>=0 && Float.isFinite(identity) && identity>=0 && identity<=1
			&& Float.isFinite(columns) && columns>0 && modelLength>0 && Double.isFinite(offset) && offset>=0 && offset<1)
			: "Shared winner scores require actual finite alignment evidence and the bound model length";
		if(rule==MATCHES){return matches;}
		if(rule==CIGAR_SCORE){return ((double)identity-offset)*columns;}
		if(rule==MODEL_SCORE){return modelScore(identity, modelLength, offset);}
		throw new IllegalArgumentException("Only the explicit M, B and B2 rules consume matched alignment evidence");
	}
	/** Shared by winner selection and the opt-in score-stagnation policy, including
	 * failing alignments; no fabricated match/column evidence is needed. */
	public static double modelScore(final float identity, final int modelLength, final double offset){
		assert(Float.isFinite(identity) && identity>=0 && identity<=1 && modelLength>0
			&& Double.isFinite(offset) && offset>=0 && offset<1)
			: "Model scores require the actual finite identity, full model length and declared offset";
		return ((double)identity-offset)*modelLength;
	}

	public static boolean prefer(final int rule, final float matches, final float identity,
			final float columns, final int modelLength, final float bestMatches, final float bestIdentity,
			final float bestColumns, final int bestModelLength){
		return prefer(rule, matches, identity, columns, modelLength, bestMatches, bestIdentity, bestColumns, bestModelLength, 0.54);
	}
	public static boolean prefer(final int rule, final float matches, final float identity,
			final float columns, final int modelLength, final float bestMatches, final float bestIdentity,
			final float bestColumns, final int bestModelLength, final double offset){
		final double score=value(rule, matches, identity, columns, modelLength, offset);
		final double best=value(rule, bestMatches, bestIdentity, bestColumns, bestModelLength, offset);
		return score>best || score==best && identity>bestIdentity;
	}

	public static String name(final int rule){
		switch(rule){
			case IDENTITY: return "identity";
			case MATCHES: return "matches";
			case CIGAR_SCORE: return "score";
			case MODEL_SCORE: return "modelscore";
			default: throw new IllegalArgumentException("Unknown shared winner rule: "+rule);
		}
	}

	/** Same equation as idaligner.Tracer.postprocess; no additional kernel/trace. */
	public static float numerator(final int rawScore, final int queryLength, final int deletions){
		final long twice=(long)rawScore+queryLength+deletions;
		if(queryLength<=0 || deletions<0 || twice<0 || twice>2L*queryLength){
			throw new IllegalArgumentException("Tracer's match numerator must fit its actual query: score="
				+rawScore+", queryLength="+queryLength+", deletions="+deletions);
		}
		return twice/2f;
	}

	/** Existing traceback/physical-edge identities count an N as half a match. */
	public static float numerator(final AlignmentStats stats){
		if(stats==null || stats.matches<0 || stats.ns<0){
			throw new IllegalArgumentException("A passing traced alignment needs its existing nonnegative match/N counts");
		}
		return stats.matches+0.5f*stats.ns;
	}

	/** Call only after both candidates passed their own threshold and span gates.
	 * An exact count/identity tie retains the earlier candidate's walk position. */
	public static boolean prefer(final float matches, final float identity,
			final float bestMatches, final float bestIdentity){
		assert(Float.isFinite(matches) && matches>=0 && Float.isFinite(bestMatches) && bestMatches>=0
			&& Float.isFinite(identity) && identity>=0 && identity<=1
			&& Float.isFinite(bestIdentity) && bestIdentity>=0 && bestIdentity<=1)
			: "Passing shared candidates must compare actual finite numerators and verifier identities";
		return matches>bestMatches || matches==bestMatches && identity>bestIdentity;
	}

	private NcrnaMatchedBaseWinner(){}
}
