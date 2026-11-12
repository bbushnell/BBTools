package align2;

/** Immutable instance-specific PacBio penalties, in unshifted score units.
 * DEFAULT retains the native constants; experimental instances do not change it.
 * @author Ganyu
 */
public final class PacBioScoreParameters{

	/**
	 * Creates independent penalties with scaled native substitution costs.
	 * Gap arguments accept zero and negative values down to the native deletion
	 * opening penalty (-292). Substitutions use the nearest integer magnitude,
	 * with half ties rounded away from zero; a small scale can round them to zero.
	 * Acceptance here does not establish compatibility with every aligner helper.
	 * @param del4 Fourth deletion-tier penalty, in [-292,0]
	 * @param del5 Fifth deletion-tier penalty, in [-292,0]
	 * @param ins4 Fourth insertion-tier penalty, in [-292,0]
	 * @param ins Insertion-opening penalty, in [-292,0]
	 * @param del Deletion-opening penalty, in [-292,0]
	 * @param subScale Finite multiplier in (0,1] for the four native substitution penalties
	 * @return New immutable parameter set; DEFAULT and native constants are unchanged
	 * @throws IllegalArgumentException If the scale or a resulting penalty is outside its accepted range
	 */
	public static PacBioScoreParameters experimental(int del4, int del5, int ins4, int ins, int del, float subScale){
		//TODO: ALN204 - zero gap penalties are accepted here, but can violate
		//MultiStateAligner9PacBio.maxImperfectScore's assertion. Contract question
		//remains open; inspected prok callers do not use that score-bound helper.
		if(!Float.isFinite(subScale) || subScale<=0 || subScale>1){
			throw new IllegalArgumentException("Experimental substitution scale must be finite in (0,1]");
		}
		return new PacBioScoreParameters(del4, del5, ins4, ins, del,
			scaled(MultiStateAligner9PacBio.POINTS_SUB, subScale), scaled(MultiStateAligner9PacBio.POINTS_SUBR, subScale),
			scaled(MultiStateAligner9PacBio.POINTS_SUB2, subScale), scaled(MultiStateAligner9PacBio.POINTS_SUB3, subScale));
	}
	/** Validates unshifted penalties and caches their packed-score counterparts. */
	private PacBioScoreParameters(int d4, int d5, int i4, int i, int d, int s, int sr, int s2, int s3){
		del4=penalty(d4);del5=penalty(d5);ins4=penalty(i4);ins=penalty(i);del=penalty(d);
		sub=penalty(s);subR=penalty(sr);sub2=penalty(s2);sub3=penalty(s3);
		offDel4=del4<<MultiStateAligner9PacBio.SCOREOFFSET;offDel5=del5<<MultiStateAligner9PacBio.SCOREOFFSET;
		offIns4=ins4<<MultiStateAligner9PacBio.SCOREOFFSET;offIns=ins<<MultiStateAligner9PacBio.SCOREOFFSET;
		offDel=del<<MultiStateAligner9PacBio.SCOREOFFSET;offSub=sub<<MultiStateAligner9PacBio.SCOREOFFSET;
		offSubR=subR<<MultiStateAligner9PacBio.SCOREOFFSET;offSub2=sub2<<MultiStateAligner9PacBio.SCOREOFFSET;
		offSub3=sub3<<MultiStateAligner9PacBio.SCOREOFFSET;
	}
	/** Rounds the scaled magnitude to the nearest integer, with half ties away from zero. */
	private static int scaled(int nativePenalty, float scale){
		assert(nativePenalty<0) : "Scaling a penalty preserves its sign; half ties round away from zero";
		return -Math.round(-nativePenalty*scale);
	}
	/** Rejects rewards and costs below the native deletion-opening penalty. */
	private static int penalty(int value){
		// Initial development domain covers every native scalar penalty and approved
		// D38 candidate; reject large shifts/rewards before allocating native matrices.
		if(value<MultiStateAligner9PacBio.POINTS_DEL || value>0){
			throw new IllegalArgumentException("Experimental PacBio penalties must be in [-292,0]: "+value);
		}
		return value;
	}
	/** Tests all nine unshifted penalties for equality; null never matches. */
	public boolean sameValues(PacBioScoreParameters other){
		return other!=null && del4==other.del4 && del5==other.del5 && ins4==other.ins4 && ins==other.ins && del==other.del
			&& sub==other.sub && subR==other.subR && sub2==other.sub2 && sub3==other.sub3;
	}
	public static final PacBioScoreParameters DEFAULT=new PacBioScoreParameters(
		MultiStateAligner9PacBio.POINTS_DEL4, MultiStateAligner9PacBio.POINTS_DEL5, MultiStateAligner9PacBio.POINTS_INS4,
		MultiStateAligner9PacBio.POINTS_INS, MultiStateAligner9PacBio.POINTS_DEL, MultiStateAligner9PacBio.POINTS_SUB,
		MultiStateAligner9PacBio.POINTS_SUBR, MultiStateAligner9PacBio.POINTS_SUB2, MultiStateAligner9PacBio.POINTS_SUB3);
	/** Unshifted penalties; insertion/deletion extension tiers not listed here remain native. */
	public final int del4, del5, ins4, ins, del, sub, subR, sub2, sub3;
	/** Same penalties shifted by SCOREOFFSET for packed dynamic-programming cells. */
	final int offDel4, offDel5, offIns4, offIns, offDel, offSub, offSubR, offSub2, offSub3;
}
