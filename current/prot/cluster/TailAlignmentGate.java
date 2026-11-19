package prot.cluster;

import prot.Blosum62;
import prot.FamilyAlignmentCounts;
import prot.FamilyMatchInformation;
import prot.HbmPositionModel;
import prot.ProteinSubstitutionMatrix;

/** Traceback gate with explicit immutable substitution/gap units.
 * Existing clustering constructors retain BLOSUM62 without global mutation.
 */
public final class TailAlignmentGate {
	/** Explicit right-pinned recipe labels; legacy pin/trim defaults remain separate. */
	public static Mode parseRightMode(String label){
		if("pin_right_trim_left_aligned_query".equals(label)){return Mode.PIN_RIGHT_TRIM_LEFT;}
		if("pin_right_overhang_only_aligned_query".equals(label)){return Mode.PIN_RIGHT_OVERHANG_ONLY;}
		throw new IllegalArgumentException("Unknown explicit right-pinned alignmentmode: "+label);
	}

	public TailAlignmentGate(int maxLength_,int minSpan_,double minIdentity_,double minScoreRatio_){
		this(Mode.PIN_RIGHT_TRIM_LEFT, Denominator.ALIGNED_QUERY, maxLength_, minSpan_,
			minIdentity_, minScoreRatio_);
	}

	public TailAlignmentGate(Mode mode_,Denominator denominator_,int maxLength_,int minSpan_,
			double minIdentity_,double minScoreRatio_){
		this(mode_,denominator_,maxLength_,minSpan_,minIdentity_,minScoreRatio_,null);
	}

	/** Optional information uses the unchanged default BLOSUM path; it never changes the gate. */
	public TailAlignmentGate(Mode mode_,Denominator denominator_,int maxLength_,int minSpan_,
			double minIdentity_,double minScoreRatio_,FamilyMatchInformation information_){
		this(mode_,denominator_,maxLength_,minSpan_,minIdentity_,minScoreRatio_,information_,ProteinSubstitutionMatrix.BLOSUM62);
	}
	/** Explicit matrix binds emissions, gap costs, traceback and self normalization. */
	public TailAlignmentGate(Mode mode_,Denominator denominator_,int maxLength_,int minSpan_,
			double minIdentity_,double minScoreRatio_,FamilyMatchInformation information_,ProteinSubstitutionMatrix matrix_){
		this(mode_,denominator_,maxLength_,minSpan_,minIdentity_,minScoreRatio_,information_,matrix_,0);
	}

	/** Optional paired-residue floor; zero preserves the historical span-only gate. */
	public TailAlignmentGate(Mode mode_,Denominator denominator_,int maxLength_,int minSpan_,
			double minIdentity_,double minScoreRatio_,FamilyMatchInformation information_,ProteinSubstitutionMatrix matrix_,int minPaired_){
		if(minPaired_<0 || minPaired_>maxLength_){throw new IllegalArgumentException("Require 0 <= minPaired <= maxLength");}
		minPaired=minPaired_;
		if(matrix_==null){throw new IllegalArgumentException("Tail gate requires an explicit substitution matrix");}
		matrix=matrix_;
		if(mode_==null || denominator_==null){
			throw new IllegalArgumentException("Tail gate mode and denominator are required");
		}
		if(maxLength_<1 || minSpan_<1 || minSpan_>maxLength_){
			throw new IllegalArgumentException("Require 1 <= minSpan <= maxLength");
		}
		if(!Double.isFinite(minIdentity_) || minIdentity_<0 || minIdentity_>1){
			throw new IllegalArgumentException("minIdentity must be a fraction in [0,1]: "+minIdentity_);
		}
		if(!Double.isFinite(minScoreRatio_)){
			throw new IllegalArgumentException("minScoreRatio must be finite: "+minScoreRatio_);
		}
		mode=mode_; denominator=denominator_; maxLength=maxLength_; minSpan=minSpan_;
		minIdentity=minIdentity_; minScoreRatio=minScoreRatio_;
		information=information_;
	}

	/** Builds the selected-matrix consensus scorer and applies this gate. */
	public TailAlignmentDecision decide(byte[] query,byte[] representative){
		if(query==null || query.length<1){throw new IllegalArgumentException("Tail gate needs a nonempty query");}
		if(representative==null || representative.length<minSpan){
			throw new IllegalArgumentException("Representative is shorter than minSpan");
		}
		final HbmPositionModel model=HbmPositionModel.fromConsensus(representative,matrix);
		if(mode==Mode.PIN_RIGHT_OVERHANG_ONLY){
			return evaluate(query,representative,model.alignTailOverhang(query,maxLength,true));
		}else if(mode==Mode.PIN_RIGHT_TRIM_LEFT || mode==Mode.PIN_RIGHT_FULL_QUERY){
			return evaluate(query,representative,model.alignTail(query,maxLength,true,
				mode==Mode.PIN_RIGHT_TRIM_LEFT));
		}else if(mode==Mode.PIN_LEFT_TRIM_RIGHT){
			return evaluate(query,representative,model.alignHead(query,maxLength,true,true));
		}
		assert(mode==Mode.FULL_QUERY_UNPINNED) : "Unhandled tail alignment mode: "+mode;
		final HbmPositionModel.Result alignment=model.align(query,true);
		if(alignment.score<=0){return TailAlignmentDecision.reject(TailAlignmentDecision.Reason.NO_ALIGNMENT);}
		return evaluateTrace(query,representative,new HbmPositionModel.TailResult(alignment.score,
			0,query.length-1,alignment.start,alignment.end,alignment.path),mode);
	}

	/** Applies this gate with a worker-local scorer already reset to the representative bytes. */
	public TailAlignmentDecision decide(byte[] query,byte[] representative,HbmPositionModel.ConsensusScorer model){
		return decide(query,representative,model,null);
	}

	/** Optional caller-owned development counts; never changes the chosen path. */
	public TailAlignmentDecision decide(byte[] query,byte[] representative,HbmPositionModel.ConsensusScorer model,
			FamilyAlignmentCounts counts){
		if(counts!=null){counts.clear();}
		if(model==null){throw new IllegalArgumentException("Tail gate needs a reusable consensus scorer");}
		if(!matrix.same(model.matrix)){throw new IllegalArgumentException("Consensus scorer and gate use different substitution matrices");}
		if(query==null || query.length<1){throw new IllegalArgumentException("Tail gate needs a nonempty query");}
		if(representative==null || representative.length<minSpan){
			throw new IllegalArgumentException("Representative is shorter than minSpan");
		}
		if(!model.matchesConsensus(representative)){
			throw new IllegalArgumentException("Reusable consensus scorer differs from representative");
		}
		if(mode==Mode.PIN_RIGHT_OVERHANG_ONLY){
			return evaluate(query,representative,model.alignTailOverhang(query,maxLength,true),counts);
		}else if(mode==Mode.PIN_RIGHT_TRIM_LEFT || mode==Mode.PIN_RIGHT_FULL_QUERY){
			return evaluate(query,representative,model.alignTail(query,maxLength,true,
				mode==Mode.PIN_RIGHT_TRIM_LEFT),counts);
		}else if(mode==Mode.PIN_LEFT_TRIM_RIGHT){
			return evaluate(query,representative,model.alignHead(query,maxLength,true,true),counts);
		}
		assert(mode==Mode.FULL_QUERY_UNPINNED) : "Unhandled tail alignment mode: "+mode;
		final HbmPositionModel.Result alignment=model.align(query,true);
		if(alignment.score<=0){return TailAlignmentDecision.reject(TailAlignmentDecision.Reason.NO_ALIGNMENT);}
		return evaluateTrace(query,representative,new HbmPositionModel.TailResult(alignment.score,
			0,query.length-1,alignment.start,alignment.end,alignment.path),mode,counts);
	}

	/** Evaluates a recorded {@link HbmPositionModel#alignTail} trace. */
	public TailAlignmentDecision evaluate(byte[] query,byte[] representative,HbmPositionModel.TailResult alignment){
		return evaluateTrace(query,representative,alignment,mode);
	}
	public TailAlignmentDecision evaluate(byte[] query,byte[] representative,HbmPositionModel.TailResult alignment,
			FamilyAlignmentCounts counts){
		return evaluateTrace(query,representative,alignment,mode,counts);
	}

	private TailAlignmentDecision evaluateTrace(byte[] query,byte[] representative,
			HbmPositionModel.TailResult alignment,Mode mode){
		return evaluateTrace(query,representative,alignment,mode,null);
	}
	private TailAlignmentDecision evaluateTrace(byte[] query,byte[] representative,
			HbmPositionModel.TailResult alignment,Mode mode,FamilyAlignmentCounts counts){
		if(counts!=null){counts.clear();}
		if(query==null || representative==null){
			throw new IllegalArgumentException("Tail gate needs non-null query and representative");
		}
		final boolean overhangOnly=mode==Mode.PIN_RIGHT_OVERHANG_ONLY;
		final boolean requireRightPinned=overhangOnly || mode==Mode.PIN_RIGHT_TRIM_LEFT || mode==Mode.PIN_RIGHT_FULL_QUERY;
		final boolean requireLeftPinned=mode==Mode.PIN_LEFT_TRIM_RIGHT;
		final boolean requireQueryEnd=mode!=Mode.PIN_LEFT_TRIM_RIGHT;
		if(alignment==null){return TailAlignmentDecision.reject(TailAlignmentDecision.Reason.NO_ALIGNMENT);}
		if(alignment.path==null || alignment.path.length<1 || (!overhangOnly && alignment.score<=0) ||
				alignment.queryStart<0 || alignment.queryStart>query.length ||
				alignment.queryEnd<alignment.queryStart-(overhangOnly ? 1 : 0) ||
				alignment.queryEnd>=query.length || alignment.start<0 ||
				alignment.start>representative.length || alignment.end<alignment.start-(overhangOnly ? 1 : 0) ||
				alignment.end>=representative.length ||
				(requireQueryEnd && alignment.queryEnd!=query.length-1) ||
				(requireRightPinned && alignment.end!=representative.length-1) ||
				(overhangOnly && alignment.queryStart!=Math.max(0,query.length-maxLength) &&
					alignment.start!=Math.max(0,representative.length-maxLength)) ||
				(requireLeftPinned && (alignment.queryStart!=0 || alignment.start!=0))){
			throw new IllegalArgumentException("Gate requires a "+(overhangOnly ? "signed" : "positive")+" recorded trace consuming the configured query edge"+
				(requireRightPinned ? " and ending at the reference endpoint" :
					requireLeftPinned ? " and starting at the reference endpoint" : " within the reference"));
		}
		final int querySpan=alignment.queryEnd-alignment.queryStart+1;
		if(querySpan<minSpan){return reject(TailAlignmentDecision.Reason.QUERY_SPAN_SHORT,alignment);}
		final int referenceSpan=alignment.end-alignment.start+1;
		if(referenceSpan<minSpan){return reject(TailAlignmentDecision.Reason.REFERENCE_SPAN_SHORT,alignment);}
		int q=alignment.queryStart, r=alignment.start, identities=0, paired=0, score=0, alignedSelf=0;
		double matchedInformation=0;
		for(byte op : alignment.path){
			if(op=='m'){
				if(q>=query.length || r>=representative.length){
					throw new IllegalArgumentException("Tail trace overran a paired column");
				}
				final byte qa=query[q], ra=representative[r];
				validateResidue(qa, "query", q);
				validateResidue(ra, "representative", r);
				if(counts!=null){counts.paired(qa,ra);}
				if(qa==ra && Blosum62.isStandard(qa)){
					identities++;
					if(information!=null){matchedInformation+=information.weight(qa);}
				}
				paired++;
				score+=HbmPositionModel.SCALE*matrix.score(qa,ra);
				alignedSelf+=HbmPositionModel.SCALE*matrix.score(qa,qa);
				q++; r++;
			}else if(op=='I'){
				if(q>=query.length){throw new IllegalArgumentException("Tail trace overran query insertion");}
				validateResidue(query[q], "query", q);
				alignedSelf+=HbmPositionModel.SCALE*matrix.score(query[q],query[q]);
				score-=HbmPositionModel.SCALE*matrix.linearGap;
				q++;
			}else if(op=='D'){
				if(r>=representative.length){throw new IllegalArgumentException("Tail trace overran representative deletion");}
				validateResidue(representative[r], "representative", r);
				if(counts!=null){counts.deleted(representative[r]);}
				score-=HbmPositionModel.SCALE*matrix.linearGap;
				r++;
			}else{
				throw new IllegalArgumentException("Unknown tail trace operation: "+(char)op);
			}
		}
		if(q!=alignment.queryEnd+1 || r!=alignment.end+1 || score!=alignment.score){
			throw new IllegalArgumentException("Recorded trace does not match its declared alignment");
		}
		final int self=denominator==Denominator.ALIGNED_QUERY ? alignedSelf : selfScore(query);
		if(alignedSelf<=0 || self<=0){throw new IllegalArgumentException("Query self-score must be positive");}
		if(counts!=null){counts.finish();}
		final double identity=identities/(double)alignment.path.length;
		final double scoreRatio=score/(double)self;
		if(paired<minPaired){
			return decision(TailAlignmentDecision.Reason.PAIRED_SPAN_SHORT,alignment,
				identities,paired,alignedSelf,self,identity,scoreRatio,matchedInformation,information!=null);
		}
		if(overhangOnly && score<=0){
			return decision(TailAlignmentDecision.Reason.NONPOSITIVE_SCORE,alignment,
				identities,paired,alignedSelf,self,identity,scoreRatio,matchedInformation,information!=null);
		}
		if(identity<minIdentity){
			return decision(TailAlignmentDecision.Reason.IDENTITY_LOW,alignment,
				identities,paired,alignedSelf,self,identity,scoreRatio,matchedInformation,information!=null);
		}
		if(scoreRatio<minScoreRatio){
			return decision(TailAlignmentDecision.Reason.SCORE_RATIO_LOW,alignment,
				identities,paired,alignedSelf,self,identity,scoreRatio,matchedInformation,information!=null);
		}
		return decision(TailAlignmentDecision.Reason.ACCEPT,alignment,
			identities,paired,alignedSelf,self,identity,scoreRatio,matchedInformation,information!=null);
	}

	private static TailAlignmentDecision reject(TailAlignmentDecision.Reason reason,
			HbmPositionModel.TailResult alignment){
		return TailAlignmentDecision.of(reason,alignment.score,alignment.queryStart,
			alignment.queryEnd,alignment.start,alignment.end,0,0,0,0,0,0,0,null);
	}

	private static TailAlignmentDecision decision(TailAlignmentDecision.Reason reason,
			HbmPositionModel.TailResult alignment,int identities,int paired,
			int alignedSelf,int denominatorSelf,double identity,double ratio,
			double matchedInformation,boolean informationMeasured){
		return TailAlignmentDecision.of(reason,alignment.score,alignment.queryStart,
			alignment.queryEnd,alignment.start,alignment.end,alignment.path.length,
			paired,identities,alignedSelf,denominatorSelf,identity,ratio,alignment.path,
			matchedInformation,informationMeasured);
	}

	private int selfScore(byte[] query){
		int self=0;
		for(int i=0; i<query.length; i++){
			validateResidue(query[i], "query", i);
			self+=HbmPositionModel.SCALE*matrix.score(query[i],query[i]);
		}
		return self;
	}

	private static void validateResidue(byte residue,String side,int index){
		if(residue<0 || residue>Blosum62.X_CODE || residue==20){
			throw new IllegalArgumentException("Invalid "+side+" residue "+residue+" at "+index);
		}
	}

	public enum Mode {PIN_RIGHT_TRIM_LEFT, PIN_RIGHT_FULL_QUERY, PIN_LEFT_TRIM_RIGHT, FULL_QUERY_UNPINNED, PIN_RIGHT_OVERHANG_ONLY}
	public enum Denominator {ALIGNED_QUERY, FULL_QUERY}

	public final Mode mode;
	public final Denominator denominator;
	public final int maxLength, minSpan, minPaired;
	public final double minIdentity, minScoreRatio;
	private final FamilyMatchInformation information;
	public final ProteinSubstitutionMatrix matrix;
}
