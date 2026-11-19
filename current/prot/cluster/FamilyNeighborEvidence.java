package prot.cluster;

import java.util.Arrays;

import ddl.FamilySketchBuilder;
import ddl.MutableFamilyIndex;
import prot.HbmPositionModel;

/** Nearest-family quality evidence against one frozen population, without changing models.
 * Moving inputs are excluded before native DDL shortlist admission. An explicit
 * sibling or other moving input can be scored separately, outside that quota.
 * The caller keeps all model payloads immutable for this helper's lifetime and
 * supplies a separate Workspace to each simultaneously active query worker.
 * @author Brian Bushnell, Fischl
 */
public final class FamilyNeighborEvidence{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds a dense stable-ID-ordered index of existing family-modal sketches.
	 * Pinned geometry and span requirements are retained; identity/ratio cutoffs
	 * are deliberately ignored because even a weak valid neighbor is evidence.
	 */
	public FamilyNeighborEvidence(final Family[] frozenModels, final FamilySketchBuilder builder,
			final TailAlignmentGate pinnedGeometry, final int shortlist_, final int buildThreads){
		if(frozenModels==null || builder==null || pinnedGeometry==null || shortlist_<1 ||
			shortlist_==Integer.MAX_VALUE || buildThreads<1 ||
			pinnedGeometry.mode==TailAlignmentGate.Mode.FULL_QUERY_UNPINNED){
			throw new IllegalArgumentException("Neighbor evidence needs frozen models, pinned geometry and positive bounded limits");
		}
		width=builder.width(); shortlist=shortlist_;
		memberInputs=null;memberRanks=null;
		geometry=new TailAlignmentGate(pinnedGeometry.mode,pinnedGeometry.denominator,
			pinnedGeometry.maxLength,pinnedGeometry.minSpan,0,0,null,pinnedGeometry.matrix);
		models=frozenModels.clone();
		for(Family model:models){validateModel(model);}
		Arrays.sort(models,(a,b)->Integer.compare(a.familyId,b.familyId));
		for(int i=1; i<models.length; i++){
			if(models[i].familyId==models[i-1].familyId){throw new IllegalArgumentException("Duplicate frozen family ID: "+models[i].familyId);}
		}
		view=MutableFamilyIndex.snapshot(width,models.length,new MutableFamilyIndex.SketchSource(){
			@Override public int size(){return models.length;}
			@Override public char[] sketch(final int dense){return models[dense].sketchUnsafe();}
		},buildThreads);
	}

	/** Creates a ranked immutable view without rebuilding the frozen index.
	 * The caller keeps models unchanged; legacy workspaces retain their original owner.
	 */
	public FamilyNeighborEvidence withMemberRanking(final MemberRankInputs inputs){
		if(inputs==null || memberInputs!=null){
			throw new IllegalArgumentException("Ranked evidence requires one original cache and cannot replace an existing binding");
		}
		return new FamilyNeighborEvidence(this,inputs);
	}

	/** Shares only frozen model/index state; target descriptors are captured once in stable-ID order. */
	private FamilyNeighborEvidence(final FamilyNeighborEvidence source,final MemberRankInputs inputs){
		assert(source!=null && inputs!=null) : "Ranked evidence views require an existing immutable index and original cache";
		models=source.models;view=source.view;geometry=source.geometry;width=source.width;shortlist=source.shortlist;
		memberInputs=inputs;memberRanks=new MemberRankStatistics.Snapshot[models.length];
		for(int i=0;i<models.length;i++){
			validateMemberBinding(models[i]);
			memberRanks[i]=models[i].memberRankingSnapshot();
		}
	}

	/** A workspace belongs to this frozen population and can be reused for sequential proposals. */
	public Workspace newWorkspace(){return new Workspace(this);}

	/*--------------------------------------------------------------*/
	/*----------------          Evidence            ----------------*/
	/*--------------------------------------------------------------*/

	/** Scores a query against up to shortlist eligible frozen families plus one explicit neighbor.
	 * Exclusion arguments are stable family IDs, with -1 meaning absent. The query
	 * ID and explicit-neighbor ID are also excluded from the frozen shortlist.
	 * The explicit neighbor is still scored when its ID is one of the moving inputs;
	 * this permits fusion-input and fission-sibling evidence without stale duplicates.
	 * It must have a different ID from the query. No family-size restriction applies.
	 */
	public Result evaluate(final Family query, final int excludeFamilyIdA, final int excludeFamilyIdB,
			final Family explicitNeighbor, final Workspace work){
		validateModel(query);
		if(memberInputs!=null){validateMemberBinding(query);}
		if(excludeFamilyIdA< -1 || excludeFamilyIdB< -1 || work==null || work.owner!=this || work.inUse){
			throw new IllegalArgumentException("Evidence needs valid excluded IDs and exclusive scratch from this frozen population");
		}
		if(explicitNeighbor!=null){
			validateModel(explicitNeighbor);
			if(memberInputs!=null){validateMemberBinding(explicitNeighbor);}
			if(explicitNeighbor.familyId==query.familyId){throw new IllegalArgumentException("An explicit neighbor must not be the query family");}
		}
		work.begin(query,excludeFamilyIdA,excludeFamilyIdB,explicitNeighbor);
		try{
			final int indexed=models.length==0 ? 0 : view.rankIndexed(query.sketchUnsafe(),
				work.candidates.length,work.candidates,null,memberInputs==null ? work : work.ranking,work.indexScratch,work);
			for(int i=0; i<indexed; i++){work.consider(models[work.candidates[i]], false);}
			if(explicitNeighbor!=null){work.consider(explicitNeighbor, true);}
			assert(work.attempted==indexed+(explicitNeighbor==null ? 0 : 1)) :
				"Every shortlisted or explicit neighbor must receive exactly one alignment";
			return new Result(work,indexed,explicitNeighbor==null ? 0 : 1,
				models.length==0 ? 0 : work.indexScratch.postingsVisited());
		}finally{
			work.query=null; work.inUse=false;
		}
	}

	/** An empty eligible population is known zero quality (available=true, ID=-1).
	 * Attempted candidates with no valid alignment are unavailable, also with finite zero.
	 * Explicit fields describe only the supplied directed pair, independently of the
	 * winning neighbor. Missing/invalid pairs have explicitAvailable=false and zero
	 * values. explicitScoreRatio is raw; explicitQuality uses the same clipped-ratio
	 * mean as ordinary neighbor quality. Neither field is a geometric distance.
	 */
	public static final class Result{
		private Result(final Workspace work,final int indexed,final int explicit,final long visited){
			available=work.valid>0 || work.attempted==0; quality=work.bestQuality; nearestFamilyId=work.nearest;
			identity=work.bestIdentity; scoreRatio=work.bestRatio; sharedFraction=work.bestShared;
			indexedCandidates=indexed; explicitCandidates=explicit;
			attemptedAlignments=work.attempted; validAlignments=work.valid; postingsVisited=visited;
			explicitAvailable=work.explicitAvailable;
			explicitIdentity=work.explicitIdentity; explicitScoreRatio=work.explicitScoreRatio;
			explicitSharedFraction=work.explicitSharedFraction; explicitQuality=work.explicitQuality;
			assert(Double.isFinite(quality) && quality>=0 && quality<=1 &&
				validAlignments<=attemptedAlignments && available==(validAlignments>0 || attemptedAlignments==0)) :
				"Local objective evidence must be finite and distinguish missing neighbors from real scores";
			assert(Double.isFinite(explicitIdentity) && Double.isFinite(explicitScoreRatio) &&
				Double.isFinite(explicitSharedFraction) && Double.isFinite(explicitQuality) &&
				(explicitAvailable ? explicitCandidates==1 && validAlignments>0 :
					explicitIdentity==0 && explicitScoreRatio==0 && explicitSharedFraction==0 && explicitQuality==0)) :
				"Directed explicit-pair evidence needs an actual valid alignment; missing values use finite zero sentinels";
		}
		public final boolean available;
		public final double quality,identity,scoreRatio,sharedFraction;
		public final int nearestFamilyId,indexedCandidates,explicitCandidates,attemptedAlignments,validAlignments;
		public final long postingsVisited;
		public final boolean explicitAvailable;
		public final double explicitIdentity,explicitScoreRatio,explicitSharedFraction,explicitQuality;
	}

	/** Worker-owned alignment, candidate and index scratch; callbacks see one fixed query context. */
	public static final class Workspace implements MutableFamilyIndex.Cosine,MutableFamilyIndex.Eligibility{
		private Workspace(final FamilyNeighborEvidence owner_){
			owner=owner_; scorer=new HbmPositionModel.ConsensusScorer(owner.geometry.matrix); candidates=new int[Math.min(owner.shortlist,owner.models.length)];
			queryRank=owner.memberInputs==null ? null : new MemberRankStatistics.View();
			ranking=owner.memberInputs==null ? null : new MutableFamilyIndex.Ranking(){
				@Override public double score(final int dense){return Workspace.this.score(dense);}
				@Override public double score(final int dense,final int hits){
					final MemberRankStatistics.Snapshot target=owner.memberRanks[dense];
					final double ratio=Math.min(queryRank.medianLength,target.medianLength)/Math.max(queryRank.medianLength,target.medianLength);
					final double cosine=Workspace.this.score(dense),similarity=1-Math.abs(queryRank.meanKd-target.meanKd)/9;
					assert(inUse && hits>=0 && ratio>0 && ratio<=1 && cosine>=0 && cosine<=1 && similarity>=0 && similarity<=1) :
						"Frozen member descriptors and native cosine must fit complete-score index bounds";
					return hits*Math.sqrt(Math.sqrt(ratio))+1.5*cosine+owner.memberInputs.hydroWeight*similarity;
				}
				@Override public double lowerBound(final int hits){return 0;}
				@Override public double upperBound(final int hits){return hits+1.5+owner.memberInputs.hydroWeight;}
			};
		}

		/** Sets a complete evidence context before invoking any index callback. */
		private void begin(final Family query_,final int excludeA_,final int excludeB_,final Family explicit){
			assert(!inUse) : "Caller must not reuse neighbor scratch until its previous evaluation returns";
			query=query_; excludeA=excludeA_; excludeB=excludeB_;
			if(owner.memberInputs!=null){query.readMemberRanking(queryRank);}
			explicitId=explicit==null ? -1 : explicit.familyId;
			nearest=-1; attempted=0; valid=0;
			bestQuality=bestIdentity=bestRatio=bestShared=0;
			explicitAvailable=false;
			explicitIdentity=explicitScoreRatio=explicitSharedFraction=explicitQuality=0;
			queryNorm=0;
			for(float value:query.monomersUnsafe()){queryNorm+=value*(double)value;}
			inUse=true;
		}

		/** Excludes moving/self/explicit IDs before positive-hit and zero-hit admission alike. */
		@Override public boolean accept(final int dense){
			assert(inUse && query!=null) : "Eligibility callbacks require the active immutable query context";
			final int id=owner.models[dense].familyId;
			return id!=query.familyId && id!=excludeA && id!=excludeB && id!=explicitId;
		}

		/** Native monomer cosine resolves DDL hit-count ties before stable family-ID order. */
		@Override public double score(final int dense){
			assert(inUse && query!=null) : "Cosine callbacks require the active immutable query context";
			final float[] a=query.monomersUnsafe(),b=owner.models[dense].monomersUnsafe();
			double dot=0,norm=0;
			for(int i=0; i<a.length; i++){dot+=a[i]*(double)b[i]; norm+=b[i]*(double)b[i];}
			return queryNorm==0 || norm==0 ? 0 : clamp(dot/Math.sqrt(queryNorm*norm));
		}

		/** Max mean of identity, clipped BLOSUM ratio and nonzero matches divided by full width. */
		private void consider(final Family target, final boolean explicit){
			assert(inUse && target.familyId!=query.familyId) : "A neighbor alignment must not compare a family to itself";
			assert(!explicit || target.familyId==explicitId) : "Only the separately supplied neighbor may populate explicit-pair evidence";
			attempted++;
			scorer.resetConsensus(target.consensusUnsafe());
			final TailAlignmentDecision decision=owner.geometry.decide(query.consensusUnsafe(),target.consensusUnsafe(),scorer);
			if(!decision.accepted()){return;}
			if(!Double.isFinite(decision.identity) || !Double.isFinite(decision.scoreRatio) ||
				decision.identity<0 || decision.identity>1){
				throw new IllegalStateException("Pinned neighbor alignment produced invalid normalized evidence");
			}
			valid++;
			int shared=0;
			final char[] a=query.sketchUnsafe(),b=target.sketchUnsafe();
			for(int bin=0; bin<owner.width; bin++){if(a[bin]!=0 && a[bin]==b[bin]){shared++;}}
			final double fraction=shared/(double)owner.width,ratio=clamp(decision.scoreRatio);
			final double quality=(decision.identity+ratio+fraction)/3;
			if(explicit){
				explicitAvailable=true; explicitIdentity=decision.identity;
				explicitScoreRatio=decision.scoreRatio;
				explicitSharedFraction=fraction; explicitQuality=quality;
			}
			if(nearest<0 || quality>bestQuality || quality==bestQuality && target.familyId<nearest){
				nearest=target.familyId; bestQuality=quality;
				bestIdentity=decision.identity; bestRatio=ratio; bestShared=fraction;
			}
		}
		private final FamilyNeighborEvidence owner;
		private final HbmPositionModel.ConsensusScorer scorer;
		private final MutableFamilyIndex.Workspace indexScratch=new MutableFamilyIndex.Workspace();
		private final int[] candidates;
		private final MutableFamilyIndex.Ranking ranking;
		private final MemberRankStatistics.View queryRank;
		private Family query;
		private int excludeA,excludeB,explicitId,nearest,attempted,valid;
		private double queryNorm,bestQuality,bestIdentity,bestRatio,bestShared;
		private double explicitIdentity,explicitScoreRatio,explicitSharedFraction,explicitQuality;
		private boolean explicitAvailable;
		private boolean inUse;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Validation          ----------------*/
	/*--------------------------------------------------------------*/

	/** A frozen roster or proposed model must use the same width and valid pinned span. */
	private void validateModel(final Family model){
		if(model==null || model.familyId<0 || model.consensusUnsafe()==null ||
			model.consensusUnsafe().length<geometry.minSpan || model.sketchUnsafe()==null ||
			model.sketchUnsafe().length!=width){
			throw new IllegalArgumentException("Neighbor evidence requires a valid family with matching sketch width and pinned span");
		}
	}
	private static double clamp(final double value){return Math.max(0,Math.min(1,value));}

	/** Ordinals cannot be interpreted through a different population or a separately built cache. */
	private void validateMemberBinding(final Family model){
		if(!model.memberRankingUses(memberInputs)){
			throw new IllegalArgumentException("Neighbor evidence requires the same original member cache for family "+model.familyId);
		}
	}

	private final Family[] models;
	private final int width,shortlist;
	private final TailAlignmentGate geometry;
	private final MutableFamilyIndex.FrozenView view;
	private final MemberRankInputs memberInputs;
	private final MemberRankStatistics.Snapshot[] memberRanks;
}
