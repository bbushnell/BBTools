package prot.cluster;

import java.util.Arrays;

/** Immutable result of one recorded BLOSUM clustering gate. */
public final class TailAlignmentDecision {

	/** Accepted candidates sort by better normalized score, then identity, score and stable ID. */
	public static int compareAccepted(TailAlignmentDecision a,int familyIdA,
			TailAlignmentDecision b,int familyIdB){
		if(!a.accepted() || !b.accepted()){
			throw new IllegalArgumentException("Best-match comparison requires two accepted candidates");
		}
		int c=Double.compare(b.scoreRatio, a.scoreRatio);
		if(c!=0){return c;}
		c=Double.compare(b.identity, a.identity);
		if(c!=0){return c;}
		if(a.score!=b.score){return a.score>b.score ? -1 : 1;}
		return familyIdA<familyIdB ? -1 : familyIdA==familyIdB ? 0 : 1;
	}

	static TailAlignmentDecision reject(Reason reason){
		return new TailAlignmentDecision(reason,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,null);
	}

	static TailAlignmentDecision of(Reason reason,int score,int queryStart,int queryEnd,
			int referenceStart,int referenceEnd,int columns,int paired,int identities,
			int alignedSelfScore,int denominatorSelfScore,double identity,double scoreRatio,
			byte[] path){
		return new TailAlignmentDecision(reason,score,queryStart,queryEnd,referenceStart,
			referenceEnd,span(queryStart,queryEnd),span(referenceStart,referenceEnd),
			columns,paired,identities,alignedSelfScore,denominatorSelfScore,identity,
			scoreRatio,path);
	}

	private TailAlignmentDecision(Reason reason_,int score_,int queryStart_,int queryEnd_,
			int referenceStart_,int referenceEnd_,int querySpan_,int referenceSpan_,
			int columns_,int paired_,int identities_,int alignedSelfScore_,
			int denominatorSelfScore_,double identity_,double scoreRatio_,byte[] path_){
		this(reason_,score_,queryStart_,queryEnd_,referenceStart_,referenceEnd_,querySpan_,referenceSpan_,
			columns_,paired_,identities_,alignedSelfScore_,denominatorSelfScore_,identity_,scoreRatio_,path_,0,false);
	}

	static TailAlignmentDecision of(Reason reason,int score,int queryStart,int queryEnd,
			int referenceStart,int referenceEnd,int columns,int paired,int identities,
			int alignedSelfScore,int denominatorSelfScore,double identity,double scoreRatio,
			byte[] path,double information,boolean informationMeasured){
		return new TailAlignmentDecision(reason,score,queryStart,queryEnd,referenceStart,referenceEnd,
			span(queryStart,queryEnd),span(referenceStart,referenceEnd),columns,paired,identities,
			alignedSelfScore,denominatorSelfScore,identity,scoreRatio,path,information,informationMeasured);
	}

	private TailAlignmentDecision(Reason reason_,int score_,int queryStart_,int queryEnd_,
			int referenceStart_,int referenceEnd_,int querySpan_,int referenceSpan_,
			int columns_,int paired_,int identities_,int alignedSelfScore_,
			int denominatorSelfScore_,double identity_,double scoreRatio_,byte[] path_,
			double information_,boolean informationMeasured_){
		if(reason_==null){throw new IllegalArgumentException("Tail gate reason is required");}
		if(!Double.isFinite(information_) || information_<0 || !informationMeasured_ && information_!=0){
			throw new IllegalArgumentException("Information must be finite, nonnegative and zero when unmeasured");
		}
		information=information_;informationMeasured=informationMeasured_;
		reason=reason_; score=score_; queryStart=queryStart_; queryEnd=queryEnd_;
		referenceStart=referenceStart_; referenceEnd=referenceEnd_; querySpan=querySpan_;
		referenceSpan=referenceSpan_; columns=columns_; paired=paired_;
		identities=identities_; alignedSelfScore=alignedSelfScore_;
		denominatorSelfScore=denominatorSelfScore_; identity=identity_;
		scoreRatio=scoreRatio_;
		path=path_==null ? null : path_.clone();
	}

	private static int span(int start,int end){return start<0 || end<start ? 0 : end-start+1;}

	/** True only when the pair passed every configured gate. */
	public boolean accepted(){return reason==Reason.ACCEPT;}

	/** Defensive path copy for fixture checks and graph construction. */
	public byte[] path(){return path==null ? null : path.clone();}
	byte[] pathUnsafe(){return path;}

	public final Reason reason;
	public final int score, queryStart, queryEnd, referenceStart, referenceEnd;
	public final int querySpan, referenceSpan, columns, paired, identities;
	public final int alignedSelfScore, denominatorSelfScore;
	public final double identity, scoreRatio;
	public final double information;
	public final boolean informationMeasured;
	private final byte[] path;

	public enum Reason {
		ACCEPT,
		NO_ALIGNMENT,
		QUERY_SPAN_SHORT,
		REFERENCE_SPAN_SHORT,
		IDENTITY_LOW,
		SCORE_RATIO_LOW,
		NONPOSITIVE_SCORE,
		PAIRED_SPAN_SHORT
	}

	@Override
	public String toString(){
		return reason+" score="+score+" query="+queryStart+"-"+queryEnd+
			" reference="+referenceStart+"-"+referenceEnd+" columns="+columns+
			" paired="+paired+" identities="+identities+" identity="+identity+
			" alignedSelf="+alignedSelfScore+" denominatorSelf="+denominatorSelfScore+
			" scoreRatio="+scoreRatio;
	}

	@Override
	public int hashCode(){
		int hash=reason.hashCode();
		hash=31*hash+score;
		hash=31*hash+queryStart;
		hash=31*hash+queryEnd;
		hash=31*hash+referenceStart;
		hash=31*hash+referenceEnd;
		hash=31*hash+columns;
		hash=31*hash+paired;
		hash=31*hash+identities;
		hash=31*hash+alignedSelfScore;
		hash=31*hash+denominatorSelfScore;
		long bits=Double.doubleToLongBits(identity);
		hash=31*hash+(int)(bits^(bits>>>32));
		bits=Double.doubleToLongBits(scoreRatio);
		hash=31*hash+(int)(bits^(bits>>>32));
		hash=31*hash+Arrays.hashCode(path);
		if(informationMeasured){
			bits=Double.doubleToLongBits(information);
			hash=31*hash+(int)(bits^(bits>>>32));
			hash=31*hash+1;
		}
		return hash;
	}

	@Override
	public boolean equals(Object other){
		if(this==other){return true;}
		if(!(other instanceof TailAlignmentDecision)){return false;}
		final TailAlignmentDecision t=(TailAlignmentDecision)other;
		return reason==t.reason && score==t.score && queryStart==t.queryStart &&
			queryEnd==t.queryEnd && referenceStart==t.referenceStart &&
			referenceEnd==t.referenceEnd && querySpan==t.querySpan &&
			referenceSpan==t.referenceSpan && columns==t.columns && paired==t.paired &&
			identities==t.identities && alignedSelfScore==t.alignedSelfScore &&
			denominatorSelfScore==t.denominatorSelfScore &&
			Double.doubleToLongBits(identity)==Double.doubleToLongBits(t.identity) &&
			Double.doubleToLongBits(scoreRatio)==Double.doubleToLongBits(t.scoreRatio) &&
			informationMeasured==t.informationMeasured &&
			Double.doubleToLongBits(information)==Double.doubleToLongBits(t.information) &&
			Arrays.equals(path,t.path);
	}
}
