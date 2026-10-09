package align2;

/** Immutable development-only PacBio penalty choices. Native constants remain
 * public defaults; existing callers share DEFAULT and retain their behavior.
 * @author Ganyu
 */
public final class PacBioScoreParameters {

	public static PacBioScoreParameters experimental(int del4,int del5,int ins4,int ins,int del,float subScale){
		if(!Float.isFinite(subScale) || subScale<=0 || subScale>1){
			throw new IllegalArgumentException("Experimental substitution scale must be finite in (0,1]");
		}
		return new PacBioScoreParameters(del4,del5,ins4,ins,del,
			scaled(MultiStateAligner9PacBio.POINTS_SUB,subScale),scaled(MultiStateAligner9PacBio.POINTS_SUBR,subScale),
			scaled(MultiStateAligner9PacBio.POINTS_SUB2,subScale),scaled(MultiStateAligner9PacBio.POINTS_SUB3,subScale));
	}
	private PacBioScoreParameters(int d4,int d5,int i4,int i,int d,int s,int sr,int s2,int s3){
		del4=penalty(d4);del5=penalty(d5);ins4=penalty(i4);ins=penalty(i);del=penalty(d);
		sub=penalty(s);subR=penalty(sr);sub2=penalty(s2);sub3=penalty(s3);
		offDel4=del4<<MultiStateAligner9PacBio.SCOREOFFSET;offDel5=del5<<MultiStateAligner9PacBio.SCOREOFFSET;
		offIns4=ins4<<MultiStateAligner9PacBio.SCOREOFFSET;offIns=ins<<MultiStateAligner9PacBio.SCOREOFFSET;
		offDel=del<<MultiStateAligner9PacBio.SCOREOFFSET;offSub=sub<<MultiStateAligner9PacBio.SCOREOFFSET;
		offSubR=subR<<MultiStateAligner9PacBio.SCOREOFFSET;offSub2=sub2<<MultiStateAligner9PacBio.SCOREOFFSET;
		offSub3=sub3<<MultiStateAligner9PacBio.SCOREOFFSET;
	}
	private static int scaled(int nativePenalty,float scale){
		assert(nativePenalty<0) : "Scaling a penalty preserves its sign; half ties round away from zero";
		return -Math.round(-nativePenalty*scale);
	}
	private static int penalty(int value){
		// Initial development domain covers every native scalar penalty and approved
		// D38 candidate; reject large shifts/rewards before allocating native matrices.
		if(value<MultiStateAligner9PacBio.POINTS_DEL || value>0){
			throw new IllegalArgumentException("Experimental PacBio penalties must be in [-292,0]: "+value);
		}
		return value;
	}
	public boolean sameValues(PacBioScoreParameters other){
		return other!=null && del4==other.del4 && del5==other.del5 && ins4==other.ins4 && ins==other.ins && del==other.del
			&& sub==other.sub && subR==other.subR && sub2==other.sub2 && sub3==other.sub3;
	}
	public static final PacBioScoreParameters DEFAULT=new PacBioScoreParameters(
		MultiStateAligner9PacBio.POINTS_DEL4,MultiStateAligner9PacBio.POINTS_DEL5,MultiStateAligner9PacBio.POINTS_INS4,
		MultiStateAligner9PacBio.POINTS_INS,MultiStateAligner9PacBio.POINTS_DEL,MultiStateAligner9PacBio.POINTS_SUB,
		MultiStateAligner9PacBio.POINTS_SUBR,MultiStateAligner9PacBio.POINTS_SUB2,MultiStateAligner9PacBio.POINTS_SUB3);
	public final int del4,del5,ins4,ins,del,sub,subR,sub2,sub3;
	final int offDel4,offDel5,offIns4,offIns,offDel,offSub,offSubR,offSub2,offSub3;
}
