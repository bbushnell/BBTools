package prot.cluster;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.ReadWriteLock;

import ddl.FamilySketchBuilder;
import fileIO.ByteStreamWriter;
import map.IntHashMap2;
import prot.HbmPositionModel;
import prot.ProteinSequence;
import structures.ByteBuilder;
import structures.IntList;
import template.Accumulator;
import template.ThreadWaiter;

/** Parallel structural passes over exact original memberships.
 * Each fission wave and fusion pass has an immutable neighbor population.
 * Workers own alignment/refinement scratch and publish proposals only after join.
 * No input member is discarded and no negative sequences are consulted.
 * @author Keqing
 */
public final class LayerCakeClustering implements LayerCakeSchedule.Operations{

	public LayerCakeClustering(final Family[] initial,final ProteinSequence[] sequences_,
			final FamilySketchBuilder sketches_,final TailAlignmentGate memberGate_,final TailAlignmentGate mergeGate_,
			final LocalFamilyObjective objective_,final Options options_,final ByteStreamWriter report_){
		if(initial==null || sequences_==null || sketches_==null || memberGate_==null || mergeGate_==null || objective_==null || options_==null){
			throw new IllegalArgumentException("Structural passes need complete inputs, gates, objective and options");
		}
		if(memberGate_.mode==TailAlignmentGate.Mode.FULL_QUERY_UNPINNED || mergeGate_.mode==TailAlignmentGate.Mode.FULL_QUERY_UNPINNED){
			throw new IllegalArgumentException("Structural evidence and rebuilds require pinned alignment geometry");
		}
		if(options_.localRefineRounds!=0){throw new IllegalArgumentException("Bounded competitive refinement currently requires the nonpublishing local probe engine");}
		if(!memberGate_.matrix.same(mergeGate_.matrix)){
			throw new IllegalArgumentException("Structural member and merge gates require the same substitution matrix");
		}
		sequences=sequences_; sketches=sketches_; memberGate=memberGate_; mergeGate=mergeGate_;
		objective=objective_; options=options_; report=report_;
		final BitSet seen=new BitSet(sequences.length);
		for(int i=0; i<initial.length; i++){
			final Family family=initial[i];
			if(family==null || family.familyId!=i || family.size()<1){
				throw new IllegalArgumentException("Initial structural families must have dense ordered IDs and nonempty memberships");
			}
			for(int j=0; j<family.size(); j++){
				final int ordinal=family.memberIndex(j);
				if(ordinal<0 || ordinal>=sequences.length || sequences[ordinal]==null || seen.get(ordinal)){
					throw new IllegalArgumentException("Initial structural families must have unique valid original input ordinals");
				}
				seen.set(ordinal);
			}
		}
		expected=initial.clone(); nextId=initial.length;
		localOnly=false; constructionOnly=false; firstLocalChildId=initial.length;
		activeIndexByFamilyId=null; qualityIndexByFamilyId=null;
		ensureCache(nextId);
		if(report!=null){report.println("cycle\tphase\twave\ttemperature\tslot_cost\tinput_a\tinput_b\taccepted\treason\tgain_available\tgain\tparent_members\tleft_members\tright_members\tparent_q\tparent_n\tleft_q\tleft_n\tright_q\tright_n");}
	}
	
	public static LayerCakeClustering localProbes(final Family[] active,final Family[] evidence,
			final ProteinSequence[] sequences,final FamilySketchBuilder sketches,final TailAlignmentGate memberGate,
			final TailAlignmentGate mergeGate,final LocalFamilyObjective objective,final Options options,final int firstNewId){
		return localEngine(active,evidence,sequences,sketches,memberGate,mergeGate,objective,options,firstNewId,false,true);
	}

	/** Constructs actual refined candidates without using regression-feature availability as a gate. */
	public static LayerCakeClustering liveProbes(final Family[] active,final Family[] evidence,
			final ProteinSequence[] sequences,final FamilySketchBuilder sketches,final TailAlignmentGate memberGate,
			final TailAlignmentGate mergeGate,final LocalFamilyObjective objective,final Options options,final int firstNewId){
		return localEngine(active,evidence,sequences,sketches,memberGate,mergeGate,objective,options,firstNewId,true,false);
	}

	/** Keeps construction independent of feature availability, but shares the frozen
	 * neighbor index for the subsequent measured-gain regression. Feature measurement
	 * does not publish models or reuse the obsolete fitted acceptance function.
	 */
	public static LayerCakeClustering liveProbesWithFeatures(final Family[] active,final Family[] evidence,
			final ProteinSequence[] sequences,final FamilySketchBuilder sketches,final TailAlignmentGate memberGate,
			final TailAlignmentGate mergeGate,final LocalFamilyObjective objective,final Options options,final int firstNewId){
		return localEngine(active,evidence,sequences,sketches,memberGate,mergeGate,objective,options,firstNewId,true,true);
	}

	private static LayerCakeClustering localEngine(final Family[] active,final Family[] evidence,
			final ProteinSequence[] sequences,final FamilySketchBuilder sketches,final TailAlignmentGate memberGate,
			final TailAlignmentGate mergeGate,final LocalFamilyObjective objective,final Options options,final int firstNewId,
			final boolean constructionOnly,final boolean withFeatures){
		if(active==null || evidence==null || sequences==null || options==null){
			throw new IllegalArgumentException("Local probes need active endpoints, all evidence models, original sequences and options");
		}
		if(constructionOnly && options.localRefineRounds<1){throw new IllegalArgumentException("Live proposals require competitive refinement before scoring");}
		final Family[] sorted=active.clone();
		for(Family family:sorted){
			if(family==null){throw new IllegalArgumentException("Active local endpoints need positive current members and selected stable IDs");}
		}
		Arrays.sort(sorted,(a,b)->Integer.compare(a.familyId,b.familyId));
		final IntHashMap2 evidenceIds=new IntHashMap2(Math.max(16,evidence.length*2));
		int maxId=-1;
		for(int i=0; i<evidence.length; i++){
			final Family family=evidence[i];
			if(family==null || family.familyId<0 || evidenceIds.get(family.familyId)>=0){
				throw new IllegalArgumentException("All selected evidence models need unique stable family IDs");
			}
			evidenceIds.put(family.familyId,i); maxId=Math.max(maxId,family.familyId);
		}
		if(firstNewId<=maxId){
			throw new IllegalArgumentException("Local fission child ID must exceed every immutable selected family ID");
		}
		final BitSet seen=new BitSet(sequences.length);
		final IntHashMap2 activeIds=new IntHashMap2(Math.max(16,sorted.length*2));
		for(int i=0; i<sorted.length; i++){
			final Family family=sorted[i];
			if(family==null || family.size()<1 || activeIds.get(family.familyId)>=0 || evidenceIds.get(family.familyId)<0){
				throw new IllegalArgumentException("Active local endpoints need positive current members and selected stable IDs");
			}
			activeIds.put(family.familyId,i);
			for(int j=0; j<family.size(); j++){
				final int ordinal=family.memberIndex(j);
				if(ordinal<0 || ordinal>=sequences.length || sequences[ordinal]==null || seen.get(ordinal)){
					throw new IllegalArgumentException("Active local endpoints must have unique valid original input ordinals");
				}
				seen.set(ordinal);
			}
		}
		return new LayerCakeClustering(sorted,sequences,sketches,memberGate,mergeGate,objective,options,
			firstNewId,activeIds,withFeatures ?
			new FamilyNeighborEvidence(evidence,sketches,memberGate,options.evidenceNeighbors,options.indexThreads) : null,constructionOnly);
	}
	
	private LayerCakeClustering(final Family[] active,final ProteinSequence[] sequences_,
			final FamilySketchBuilder sketches_,final TailAlignmentGate memberGate_,final TailAlignmentGate mergeGate_,
			final LocalFamilyObjective objective_,final Options options_,final int firstNewId,
			final IntHashMap2 activeIds,final FamilyNeighborEvidence neighbors,final boolean constructionOnly_){
		this(active,sequences_,sketches_,memberGate_,mergeGate_,objective_,options_,firstNewId,
			activeIds,qualitySlots(active,firstNewId),neighbors,constructionOnly_);
	}
	
	private LayerCakeClustering(final Family[] active,final ProteinSequence[] sequences_,
			final FamilySketchBuilder sketches_,final TailAlignmentGate memberGate_,final TailAlignmentGate mergeGate_,
			final LocalFamilyObjective objective_,final Options options_,final int firstNewId,
			final IntHashMap2 activeIds,final IntHashMap2 qualitySlots,
			final FamilyNeighborEvidence neighbors,final boolean constructionOnly_){
		if(sketches_==null || memberGate_==null || mergeGate_==null || objective_==null || options_==null || activeIds==null ||
			qualitySlots==null || neighbors==null && !constructionOnly_){
			throw new IllegalArgumentException("Local probes need complete inputs, gates, objective, options and neighbor evidence");
		}
		if(memberGate_.mode==TailAlignmentGate.Mode.FULL_QUERY_UNPINNED || mergeGate_.mode==TailAlignmentGate.Mode.FULL_QUERY_UNPINNED){
			throw new IllegalArgumentException("Structural evidence and rebuilds require pinned alignment geometry");
		}
		if(!memberGate_.matrix.same(mergeGate_.matrix)){
			throw new IllegalArgumentException("Structural member and merge gates require the same substitution matrix");
		}
		sequences=sequences_; sketches=sketches_; memberGate=memberGate_; mergeGate=mergeGate_;
		objective=objective_; options=options_; report=null;
		expected=active; nextId=firstNewId; probeNeighbors=neighbors;
		localOnly=true;constructionOnly=constructionOnly_;firstLocalChildId=firstNewId;activeIndexByFamilyId=activeIds;
		qualityIndexByFamilyId=qualitySlots; ensureCache(active.length+1);
	}

	@Override public synchronized LayerCakeSchedule.PassResult run(final LayerCakeSchedule.Phase phase,final Family[] frozen,
			final double temperature,final int cycle){
		if(localOnly){throw new IllegalStateException("Local probe engines cannot publish scheduled structural passes");}
		if(probing){throw new IllegalStateException("Independent probes and a publishing schedule require separate engines");}
		if(phase==null || frozen==null || !Arrays.equals(expected,frozen)){
			throw new IllegalArgumentException("A structural pass must consume this engine's last published family population");
		}
		scheduling=true;
		final int slotCost=LayerCakeSchedule.summarize(frozen,options.familyBudget,sequences.length).slotCost;
		final LayerCakeSchedule.PassResult result=phase==LayerCakeSchedule.Phase.FISSION ?
			fission(frozen,temperature,cycle,slotCost) : fusion(frozen,temperature,cycle,slotCost);
		expected=result.families;
		return result;
	}

	/** Shares immutable indexes and gives each worker its own nearest-choice cache;
	 * construction, cohesion caches and all proposal scratch belong to this worker.
	 * The source engine never publishes models while these forks exist.
	 */
	public synchronized LayerCakeClustering forkLocalProbes(){
		return forkLocalProbes(true);
	}

	/** Fission-only live workers can skip the fusion nearest-neighbor selector until a fusion request arrives. */
	public synchronized LayerCakeClustering forkLocalProbes(final boolean withFusionQueries){
		if(!localOnly || scheduling){throw new IllegalStateException("Only a frozen local engine may create independent proposal workers");}
		forked=true;
		if(withFusionQueries && probeFusion==null){
			final int[] counts=new int[expected.length];
			for(int i=0;i<counts.length;i++){counts[i]=expected[i].size();}
			probeFusion=new MutualNearestMerge.PairQueries(expected,counts,sketches,mergeGate,
				options.fusionSizeRatio,options.mergeShortlist,options.indexThreads,memberInputs);
		}
		final LayerCakeClustering worker=new LayerCakeClustering(expected,sequences,sketches,memberGate,mergeGate,
			objective,options,nextId,activeIndexByFamilyId,qualityIndexByFamilyId,probeNeighbors,constructionOnly);
		worker.memberInputs=memberInputs;worker.memberLayer=memberLayer;worker.forked=true;
		if(withFusionQueries){worker.probeFusion=probeFusion.fork();}
		return worker;
	}

	/** Binds a fresh engine before queries, scheduling or forks. Detached proposals
	 * receive descriptors before evidence scoring; borrowed inputs stay unchanged.
	 */
	public synchronized LayerCakeClustering withMemberRanking(final MemberRankInputs inputs,final int layer){
		if(probing || scheduling || forked || memberInputs!=null || probeFusion!=null || featureWorkspace!=null){
			throw new IllegalStateException("Member ranking must bind a fresh engine before queries, scheduling or worker forks");
		}
		if(inputs==null || !inputs.belongsTo(sequences) || layer<0){
			throw new IllegalArgumentException("Engine ranking needs the identical original population cache and a nonnegative layer");
		}
		final MemberRankStatistics.View view=new MemberRankStatistics.View();
		for(Family family:expected){
			if(!family.memberRankingUses(inputs)){
				throw new IllegalArgumentException("Every active proposal parent must use the same original ordinal cache");
			}
			family.readMemberRanking(view);
			if(view.layer>layer){throw new IllegalArgumentException("Proposal measurement layer cannot precede its parent descriptors");}
		}
		// Local factories may already own an index; the ranked view shares that immutable index.
		final FamilyNeighborEvidence ranked=probeNeighbors==null ? null : probeNeighbors.withMemberRanking(inputs);
		memberInputs=inputs;memberLayer=layer;probeNeighbors=ranked;
		return this;
	}

	/** Measures actual rejected models without publishing their cohesion under old IDs.
	 * Feasible controls reuse the evidence already measured during construction.
	 */
	public synchronized ProbeFeatures measureLocalFeatures(final Probe probe){
		return measureLocalFeatures(probe,null);
	}

	/** Optional numeric capture preserves raw evidence without retaining rejected models.
	 * The live caller invokes this before any residual-pool scoring or publication.
	 */
	public synchronized ProbeFeatures measureLocalFeatures(final Probe probe, final WideProbeFeatures capture){
		if(!localOnly || probe==null){throw new IllegalArgumentException("Feature measurement requires this frozen local engine and a probe");}
		if(capture!=null && !constructionOnly){throw new IllegalArgumentException("Wide capture currently requires live construction mode");}
		if(probeNeighbors==null){throw new IllegalStateException("Feature measurement requires an explicitly enabled frozen neighbor context");}
		if(probe.features!=null){return probe.features;}
		if(constructionOnly ? probe.candidateFirst==null : probe.audit==null || !probe.audit.available()){return null;}
		final Family a=expected[activeIndexByFamilyId.get(probe.inputA)];
		final Family b=probe.inputB<0 ? null : expected[activeIndexByFamilyId.get(probe.inputB)];
		final boolean fission=probe.phase==LayerCakeSchedule.Phase.FISSION;
		final Decision d=new Decision(probe.inputA,probe.inputB);
		if(constructionOnly){
			d.parent=fission ? a : probe.candidateFirst;
			d.left=fission ? probe.candidateFirst : a;
			d.right=fission ? probe.candidateSecond : b;
		}else{
			final Audit candidate=probe.audit;
			d.parent=fission ? a : candidate.union;
			d.left=fission ? candidate.left : a;
			d.right=fission ? candidate.right : b;
		}
		if(featureScorer==null){
			featureScorer=new HbmPositionModel.ConsensusScorer(memberGate.matrix);
			featurePlacement=new TailAlignmentGate(memberGate.mode,memberGate.denominator,memberGate.maxLength,1,0,0,null,memberGate.matrix);
			featureWorkspace=probeNeighbors.newWorkspace();
		}
		final HbmPositionModel.ConsensusScorer scorer=featureScorer;
		final TailAlignmentGate placement=featurePlacement;
		final ConsensusQualityReport.Snapshot p=measureQuality(d.parent,placement,scorer);
		final ConsensusQualityReport.Snapshot l=measureQuality(d.left,placement,scorer);
		final ConsensusQualityReport.Snapshot r=measureQuality(d.right,placement,scorer);
		if(capture!=null){capture.quality(0,p);capture.quality(1,l);capture.quality(2,r);}
		if(p==null || l==null || r==null){return null;}
		d.parentQ=normalized(p);d.leftQ=normalized(l);d.rightQ=normalized(r);
		final int excludeA=a.familyId,excludeB=fission ? d.right.familyId : b.familyId;
		final FamilyNeighborEvidence.Workspace scratch=featureWorkspace;
		final FamilyNeighborEvidence.Result pn=probeNeighbors.evaluate(d.parent,excludeA,excludeB,null,scratch);
		final FamilyNeighborEvidence.Result ln=probeNeighbors.evaluate(d.left,excludeA,excludeB,d.right,scratch);
		final FamilyNeighborEvidence.Result rn=probeNeighbors.evaluate(d.right,excludeA,excludeB,d.left,scratch);
		if(capture!=null){capture.neighbor(0,pn);capture.neighbor(1,ln);capture.neighbor(2,rn);}
		if(!pn.available || !ln.available || !rn.available){return null;}
		d.parentN=pn.quality;d.leftN=ln.quality;d.rightN=rn.quality;d.leftPair=ln;d.rightPair=rn;
		// Live construction uses temperature zero because it never decides by the
		// old objective. The fitted predictor is the unscaled saved size-pressure
		// term (temperature one), not an accidentally zeroed feature column.
		return new ProbeFeatures(probe.phase,d,probe.slotCost,constructionOnly ? 1 : probe.temperature,objective.parameters);
	}

	/** Unthresholded pinned cohesion of every original member; null means a missing placement. */
	private ConsensusQualityReport.Snapshot measureQuality(final Family family,final TailAlignmentGate placement,
			final HbmPositionModel.ConsensusScorer scorer){
		assert(family!=null && family.size()>0) : "Counterfactual cohesion must measure the complete nonempty proposed membership";
		final byte[] consensus=family.consensusUnsafe();scorer.resetConsensus(consensus);
		double identity=0,ratio=0;
		for(int i=0;i<family.size();i++){
			final TailAlignmentDecision decision=placement.decide(sequences[family.memberIndex(i)].enc,consensus,scorer);
			if(!decision.accepted()){return null;}
			identity+=decision.identity;ratio+=decision.scoreRatio;
		}
		return ConsensusQualityReport.available(family.familyId,family.size(),identity,ratio);
	}

	/** Constructs exactly one independent proposal against the constructor's immutable roster.
	 * No model is published, no child recurses and objective weights/margin never censor a probe.
	 * Calls are serialized; immutable neighbor and fusion-selection contexts are reused. Returned
	 * original models are borrowed read-only; new child/union payloads belong to the caller.
	 * Each independent fission uses the same fresh right ID, because overlays must never be chained.
	 */
	public synchronized Probe probe(final LayerCakeSchedule.Phase phase,final int leftId,final int rightId,final double temperature){
		return probeByIndex(phase,leftId,rightId,temperature,false);
	}
	
	public synchronized Probe probeByFamilyId(final LayerCakeSchedule.Phase phase,final int leftFamilyId,final int rightFamilyId,
			final double temperature){
		return probeByFamilyId(phase,leftFamilyId,rightFamilyId,temperature,false);
	}

	public synchronized Probe probeByFamilyId(final LayerCakeSchedule.Phase phase,final int leftFamilyId,final int rightFamilyId,
			final double temperature,final boolean summarizeRebuild){
		return probeByFamilyId(phase,leftFamilyId,rightFamilyId,temperature,summarizeRebuild,false);
	}

	/** Optional audit candidates remain private to this independent request, never to a publishing pass. */
	public synchronized Probe probeByFamilyId(final LayerCakeSchedule.Phase phase,final int leftFamilyId,final int rightFamilyId,
			final double temperature,final boolean summarizeRebuild,final boolean auditRejected){
		if(!localOnly){throw new IllegalStateException("Stable-ID local probes require a local probe engine");}
		if(options.localRefineRounds>0 && !auditRejected){throw new IllegalArgumentException("Competitive diagnostic refinement requires actual rejected-proposal capture");}
		if(phase==null || phase==LayerCakeSchedule.Phase.FISSION && rightFamilyId!=-1){
			throw new IllegalArgumentException("A stable-ID probe needs a phase and no fission partner");
		}
		final int left=activeIndexByFamilyId.get(leftFamilyId);
		final int right=phase==LayerCakeSchedule.Phase.FISSION ? -1 : activeIndexByFamilyId.get(rightFamilyId);
		if(left<0 || phase==LayerCakeSchedule.Phase.FUSION && right<0){
			throw new IllegalArgumentException("Local probe endpoints must be active stable family IDs");
		}
		return probeByIndex(phase,left,right,temperature,true,summarizeRebuild,auditRejected);
	}

	/** The pass coordinator reserves a different right-child ID for every fission.
	 * Returned candidates are detached and may have strict member vetoes; only measured
	 * gain and the committing pass determine whether they are published.
	 */
	public synchronized Probe constructByFamilyId(final LayerCakeSchedule.Phase phase,final int leftFamilyId,
			final int rightFamilyId,final int reservedChildId){
		if(!localOnly || !constructionOnly || phase==null ||
			(phase==LayerCakeSchedule.Phase.FISSION ? reservedChildId<firstLocalChildId || reservedChildId>=Integer.MAX_VALUE-8 : reservedChildId!=-1)){
			throw new IllegalArgumentException("Live construction needs a dedicated engine and a bounded freshly reserved split ID");
		}
		final int saved=nextId;
		if(phase==LayerCakeSchedule.Phase.FISSION){nextId=reservedChildId;}
		try{
			final Probe result=probeByFamilyId(phase,leftFamilyId,rightFamilyId,0,true,true);
			if(memberInputs!=null){
				verifyCandidateRanking(result.candidateFirst);verifyCandidateRanking(result.candidateSecond);
			}
			return result;
		}
		finally{nextId=saved;}
	}

	/** Runs once after refinement and before evidence, including usable audit candidates.
	 * These candidates conserve construction members; live publication later measures winners.
	 */
	private void bindCandidate(final Family candidate,final int layer){
		if(memberInputs==null || candidate==null){return;}
		if(candidate.size()<1 || candidate.hasMemberRanking()){
			throw new IllegalStateException("A completed proposal must own fresh nonempty membership before descriptor binding");
		}
		candidate.enableMemberRanking(memberInputs,layer,null);
	}

	/** Live candidate handoff verifies early binding rather than attaching a second time. */
	private void verifyCandidateRanking(final Family candidate){
		if(candidate==null){return;}
		if(!candidate.memberRankingUses(memberInputs)){
			throw new IllegalStateException("A completed live candidate must already use the original member cache");
		}
		final MemberRankStatistics.Snapshot snapshot=candidate.memberRankingSnapshot();
		if(snapshot.members!=candidate.size() || snapshot.layer!=memberLayer){
			throw new IllegalStateException("Live candidate descriptors must cover its final construction membership and measurement layer");
		}
	}

	/** Each publishing wave/pass measures fresh models later than its frozen parents. */
	private void advanceMemberLayer(){
		if(memberInputs==null){return;}
		if(memberLayer==Integer.MAX_VALUE){throw new IllegalStateException("Structural member-ranking layer cannot advance beyond Integer.MAX_VALUE");}
		memberLayer=Math.addExact(memberLayer,1);
	}

	/** Both modes retain the same index geometry; only the optional shortlist pivot differs. */
	private FamilyNeighborEvidence neighborEvidence(final Family[] models){
		final FamilyNeighborEvidence evidence=new FamilyNeighborEvidence(models,sketches,memberGate,options.evidenceNeighbors,options.indexThreads);
		return memberInputs==null ? evidence : evidence.withMemberRanking(memberInputs);
	}
	
	private Probe probeByIndex(final LayerCakeSchedule.Phase phase,final int leftId,final int rightId,
			final double temperature,final boolean reportStableIds){
		return probeByIndex(phase,leftId,rightId,temperature,reportStableIds,false);
	}

	private Probe probeByIndex(final LayerCakeSchedule.Phase phase,final int leftId,final int rightId,
			final double temperature,final boolean reportStableIds,final boolean summarizeRebuild){
		return probeByIndex(phase,leftId,rightId,temperature,reportStableIds,summarizeRebuild,false);
	}

	private Probe probeByIndex(final LayerCakeSchedule.Phase phase,final int leftId,final int rightId,
			final double temperature,final boolean reportStableIds,final boolean summarizeRebuild,final boolean auditRejected){
		if(auditRejected && !localOnly){throw new IllegalStateException("Rejected-model audits require a nonpublishing local engine");}
		if(scheduling){throw new IllegalStateException("Independent probes cannot follow a publishing structural pass");}
		if(phase==null || !Double.isFinite(temperature) || temperature<0 || leftId<0 || leftId>=expected.length ||
			(phase==LayerCakeSchedule.Phase.FISSION ? rightId!=-1 : rightId<0 || rightId>=expected.length || rightId==leftId)){
			throw new IllegalArgumentException("A probe needs a valid initial parent or pair, finite temperature and no fission partner");
		}
		probing=true;
		if(probeSlotCost<0){probeSlotCost=LayerCakeSchedule.summarize(expected,options.familyBudget,sequences.length).slotCost;}
		final int slotCost=probeSlotCost;
		final int a=phase==LayerCakeSchedule.Phase.FISSION ? leftId : Math.min(leftId,rightId);
		final int b=phase==LayerCakeSchedule.Phase.FISSION ? -1 : Math.max(leftId,rightId);
		final int reportA=reportStableIds ? expected[a].familyId : a, reportB=b<0 ? -1 : (reportStableIds ? expected[b].familyId : b);
		MutualNearestMerge.Pair pair=null;
		if(phase==LayerCakeSchedule.Phase.FISSION){
			if(!fissionEligible(expected[a])){return new Probe(phase,reportA,reportB,slotCost,temperature,"below_fission_size_floor");}
			if(nextId>=Integer.MAX_VALUE-8){throw new IllegalStateException("No fresh fission ID remains for the independent overlay");}
		}else{
			if(Math.max(expected[a].size(),expected[b].size())/(double)Math.min(expected[a].size(),expected[b].size())>options.fusionSizeRatio){
				return new Probe(phase,reportA,reportB,slotCost,temperature,"fusion_size_ratio");
			}
			if(probeFusion==null){
				final int[] counts=new int[expected.length];
				for(int i=0; i<counts.length; i++){counts[i]=expected[i].size();}
				probeFusion=new MutualNearestMerge.PairQueries(expected,counts,sketches,mergeGate,
					options.fusionSizeRatio,options.mergeShortlist,options.indexThreads,memberInputs);
			}
			final MutualNearestMerge.PairChoice choice=probeFusion.query(a,b);
			selectionAlignments=Math.addExact(selectionAlignments,choice.attemptedAlignments);
			if(!choice.reciprocal){
				return new Probe(phase,reportA,reportB,slotCost,temperature,"fusion_not_reciprocal");
			}
			if(!choice.gatesPassed){
				return new Probe(phase,reportA,reportB,slotCost,temperature,"fusion_similarity_gate");
			}
			pair=choice.pair;
			if(pair==null || pair.leftIndex!=a || pair.rightIndex!=b){throw new IllegalStateException("A reciprocal probe must retain the requested ordered endpoints");}
		}
		if(probeNeighbors==null && !constructionOnly){
			probeNeighbors=neighborEvidence(expected);
		}
		final Work work=new Work(expected,phase==LayerCakeSchedule.Phase.FISSION ? new int[]{a} : null,
			pair==null ? null : new MutualNearestMerge.Pair[]{pair},temperature,0,1,nextId,slotCost,true,probeNeighbors,summarizeRebuild,auditRejected);
		work.execute();
		return new Probe(phase,reportA,reportB,slotCost,temperature,work.decisions[0],objective.parameters);
	}

	/** Accepted children alone enter the next wave; each wave rebuilds its full current neighbor context. */
	private LayerCakeSchedule.PassResult fission(final Family[] initial,final double temperature,final int cycle,final int slotCost){
		Family[] current=initial;
		IntList eligible=new IntList();
		for(int i=0; i<current.length; i++){if(fissionEligible(current[i])){eligible.add(i);}}
		long accepted=0;
		for(int wave=1; eligible.size>0; wave++){
			final long end=(long)nextId+eligible.size;
			if(end>Integer.MAX_VALUE-8){throw new IllegalStateException("Stable fission family-ID space exhausted");}
			final int firstRightId=nextId; nextId=(int)end; ensureCache(nextId);
			advanceMemberLayer();
			final Work work=new Work(current,eligible.toArray(),null,temperature,cycle,wave,firstRightId,slotCost);
			work.execute();
			int kept=0;
			for(Decision decision:work.decisions){if(decision.accepted){kept++;}}
			final Family[] next=Arrays.copyOf(current,Math.addExact(current.length,kept));
			final IntList children=new IntList(kept*2);
			int tail=current.length;
			for(int j=0; j<work.decisions.length; j++){
				final Decision decision=work.decisions[j]; writeDecision(work,decision);
				if(!decision.accepted){continue;}
				final int at=work.jobs[j];
				next[at]=decision.left; next[tail]=decision.right;
				publishQuality(decision.left,decision.leftQ); publishQuality(decision.right,decision.rightQ);
				if(fissionEligible(decision.left)){children.add(at);}
				if(fissionEligible(decision.right)){children.add(tail);}
				tail++;
			}
			assert(tail==next.length) : "Every accepted binary fission adds exactly one new family";
			accepted+=kept; acceptedFissions+=kept;
			children.sort(); current=next; eligible=children;
			if(wave==Integer.MAX_VALUE && eligible.size>0){throw new IllegalStateException("Fission wave counter exhausted");}
		}
		return new LayerCakeSchedule.PassResult(current,accepted);
	}

	/** One selection and one publication only; rejected pairs never fall back to a second-best partner. */
	private LayerCakeSchedule.PassResult fusion(final Family[] frozen,final double temperature,final int cycle,final int slotCost){
		if(frozen.length<2){return new LayerCakeSchedule.PassResult(frozen,0);}
		final int[] counts=new int[frozen.length];
		for(int i=0; i<counts.length; i++){counts[i]=frozen[i].size();}
		final MutualNearestMerge.Result selected=MutualNearestMerge.selectWithinSizeRatio(frozen,counts,sketches,
			mergeGate,options.fusionSizeRatio,options.mergeShortlist,options.threads,options.indexThreads,memberInputs);
		selectionAlignments=Math.addExact(selectionAlignments,selected.attemptedAlignments);
		if(selected.pairs.length==0){return new LayerCakeSchedule.PassResult(frozen,0);}
		advanceMemberLayer();
		final Work work=new Work(frozen,null,selected.pairs,temperature,cycle,1,0,slotCost);
		work.execute();
		final boolean[] removed=new boolean[frozen.length],claimed=new boolean[frozen.length];
		final Family[] replacements=new Family[frozen.length];
		int accepted=0;
		for(int i=0; i<selected.pairs.length; i++){
			final MutualNearestMerge.Pair pair=selected.pairs[i];
			if(claimed[pair.leftIndex] || claimed[pair.rightIndex]){
				throw new IllegalStateException("A frozen fusion pass cannot reuse either endpoint of any proposed pair");
			}
			claimed[pair.leftIndex]=claimed[pair.rightIndex]=true;
			final Decision decision=work.decisions[i]; writeDecision(work,decision);
			if(!decision.accepted){continue;}
			replacements[pair.leftIndex]=decision.parent; removed[pair.rightIndex]=true;
			publishQuality(decision.parent,decision.parentQ); accepted++;
		}
		final Family[] next=new Family[frozen.length-accepted];
		for(int i=0,j=0; i<frozen.length; i++){
			if(!removed[i]){next[j++]=replacements[i]==null ? frozen[i] : replacements[i];}
		}
		acceptedFusions+=accepted;
		return new LayerCakeSchedule.PassResult(next,accepted);
	}

	private boolean fissionEligible(final Family family){
		return family.size()>=options.fissionMinSize && family.size()>=2L*options.minChildSize;
	}

	private static IntHashMap2 qualitySlots(final Family[] active,final int firstNewId){
		final IntHashMap2 map=new IntHashMap2(Math.max(16,active.length*2+2));
		for(int i=0; i<active.length; i++){map.put(active[i].familyId,i);}
		map.put(firstNewId,active.length);
		return map;
	}

	/** Cache growth occurs on the coordinator only, between worker joins. */
	private void ensureCache(final int needed){
		if(needed<=quality.length){return;}
		final int capacity=(int)Math.min(Integer.MAX_VALUE-8,Math.max(needed,Math.max(16L,quality.length*2L)));
		quality=Arrays.copyOf(quality,capacity); qualityState=Arrays.copyOf(qualityState,capacity);
	}
	/** A worker may memoize its uniquely owned old input; accepted new models are published after join. */
	private void publishQuality(final Family family,final double q){
		final int slot=qualitySlot(family);
		assert(slot<quality.length && Double.isFinite(q) && q>=0 && q<=1) : "Accepted model quality needs a reserved cache slot and normalized evidence";
		quality[slot]=q; qualityState[slot]=VALID;
	}
	private int qualitySlot(final Family family){
		if(!localOnly){return family.familyId;}
		final int slot=qualityIndexByFamilyId.get(family.familyId);
		if(slot<0){throw new IllegalStateException("No local quality slot reserved for family "+family.familyId);}
		return slot;
	}

	/** Preserves explicit missingness; rejected proposals do not acquire invented zero-valued evidence. */
	private void writeDecision(final Work work,final Decision decision){
		if(report==null){return;}
		final ByteBuilder row=reportRow.clear();
		row.append(work.cycle).tab().append(work.pairs==null ? "fission" : "fusion").tab().append(work.wave)
			.tab().appendSlow(work.temperature).tab().append(work.slotCost).tab().append(decision.inputA).tab().append(decision.inputB)
			.tab().append(decision.accepted ? 1 : 0).tab().append(decision.reason).tab().append(decision.scored ? 1 : 0).tab();
		if(decision.scored){
			row.appendSlow(decision.gain).tab().append(decision.parent.size()).tab().append(decision.left.size()).tab().append(decision.right.size())
				.tab().appendSlow(decision.parentQ).tab().appendSlow(decision.parentN)
				.tab().appendSlow(decision.leftQ).tab().appendSlow(decision.leftN)
				.tab().appendSlow(decision.rightQ).tab().appendSlow(decision.rightN);
		}else{row.append("NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA");}
		report.println(row);
	}

	/** One frozen population and an independent result slot per proposal. */
	private final class Work implements Accumulator<Worker>{
		Work(final Family[] families_,final int[] jobs_,final MutualNearestMerge.Pair[] pairs_,final double temperature_,
				final int cycle_,final int wave_,final int firstRightId_,final int slotCost_){
			this(families_,jobs_,pairs_,temperature_,cycle_,wave_,firstRightId_,slotCost_,false,null,false,false);
		}
		Work(final Family[] families_,final int[] jobs_,final MutualNearestMerge.Pair[] pairs_,final double temperature_,
				final int cycle_,final int wave_,final int firstRightId_,final int slotCost_,final boolean calibration_,
				final FamilyNeighborEvidence existingNeighbors,final boolean summarizeRebuild_,final boolean auditRejected_){
			families=families_; jobs=jobs_; pairs=pairs_; temperature=temperature_; cycle=cycle_; wave=wave_; firstRightId=firstRightId_;slotCost=slotCost_;
			calibration=calibration_; summarizeRebuild=summarizeRebuild_||auditRejected_; auditRejected=auditRejected_;
			decisions=new Decision[pairs==null ? jobs.length : pairs.length];
			measurementLayer=memberLayer;
			neighbors=constructionOnly ? null : existingNeighbors==null ?
				neighborEvidence(families) : existingNeighbors;
		}
		void execute(){
			final ArrayList<Worker> workers=new ArrayList<Worker>();
			for(int i=0; i<Math.min(options.threads,decisions.length); i++){workers.add(new Worker(this));}
			try{
				if(localOnly && workers.size()==1){workers.get(0).run();accumulate(workers.get(0));}
				else if(!ThreadWaiter.startAndWait(workers,this)){fail(new IllegalStateException("Structural proposal workers did not finish"));}
			}catch(Throwable cause){
				fail(cause);
				final ArrayList<Worker> started=new ArrayList<Worker>();
				for(Worker worker:workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
				ThreadWaiter.waitForThreadsToFinish(started);
			}
			if(failure!=null || done!=decisions.length){throw new IllegalStateException("Structural pass failed before publishing models",failure);}
			for(Decision decision:decisions){if(decision==null){throw new IllegalStateException("Structural pass omitted a proposal result");}}
			if(pairs==null){fissionProposals+=done;}else{fusionProposals+=done;}
		}
		int claim(){return (int)Math.min(next.getAndAdd(BATCH),decisions.length);}
		synchronized void fail(final Throwable cause){if(failure==null){failure=cause;}}
		@Override public void accumulate(final Worker worker){
			if(!worker.success){fail(new IllegalStateException("Structural worker ended before its completion marker"));}
			done=Math.addExact(done,worker.done);
			qualityAlignments=Math.addExact(qualityAlignments,worker.qualityAlignments);
			neighborAlignments=Math.addExact(neighborAlignments,worker.neighborAlignments);
			rebuildAlignments=Math.addExact(rebuildAlignments,worker.rebuildAlignments);
			memberPostings=Math.addExact(memberPostings,worker.memberPostings);
			neighborPostings=Math.addExact(neighborPostings,worker.neighborPostings);
			fissionRetainedEdges=Math.addExact(fissionRetainedEdges,worker.retainedEdges);
			fissionUnionEdges=Math.addExact(fissionUnionEdges,worker.unionEdges);
		}
		@Override public boolean success(){return failure==null && done==decisions.length;}
		@Override public ReadWriteLock rwlock(){return null;}
		final Family[] families;
		final int[] jobs;
		final MutualNearestMerge.Pair[] pairs;
		final Decision[] decisions;
		final FamilyNeighborEvidence neighbors;
		final double temperature;
		final int cycle,wave,firstRightId,slotCost,measurementLayer;
		final boolean calibration,summarizeRebuild,auditRejected;
		final AtomicLong next=new AtomicLong();
		volatile Throwable failure;
		int done;
	}

	/** Caches only its disjoint input families; accepted replacements enter the cache after every worker joins. */
	private final class Worker extends Thread{
		Worker(final Work work_){work=work_;}
		@Override public void run(){
			try(ConsensusRefiner owner=new ConsensusRefiner(1,options.consensusTrimLeft,options.trimDepth,null,true);
					ConsensusRefiner refiner=owner.fork();
					ClusteringMetrics metrics=new ClusteringMetrics(null,null,true,1,Long.MAX_VALUE)){
				final FamilySketchBuilder builder=sketches.fork();
				final FamilyFission splitter=new FamilyFission(options.minChildSize,options.memberNeighbors,options.memberMinHits,1);
				final FamilyNeighborEvidence.Workspace neighborWork=work.neighbors==null ? null : work.neighbors.newWorkspace();
				for(int start=work.claim(); start<work.decisions.length && work.failure==null; start=work.claim()){
					final int end=(int)Math.min((long)start+BATCH,work.decisions.length);
					for(int i=start; i<end && work.failure==null; i++){
						work.decisions[i]=work.pairs==null ? split(i,builder,splitter,refiner,metrics,neighborWork) :
							merge(i,builder,refiner,metrics,neighborWork);
						done++;
					}
				}
				success=work.failure==null;
			}catch(Throwable cause){work.fail(cause);}
		}

		private Decision split(final int index,final FamilySketchBuilder builder,final FamilyFission splitter,
				final ConsensusRefiner refiner,final ClusteringMetrics metrics,final FamilyNeighborEvidence.Workspace neighborWork){
			final Family parent=work.families[work.jobs[index]];
			final Decision result=new Decision(parent.familyId,-1);
			if(!constructionOnly && !cacheQuality(parent)){result.reason="parent_quality_unavailable";return result;}
			metrics.begin("consensus",work.cycle,parent.size());
			final int rightId=work.firstRightId+index;
			final long before=metrics.constructionAlignments();
			final FamilyFission.Result proposal=splitter.split(parent,sequences,builder,memberGate,refiner,metrics,rightId,
				work.summarizeRebuild,work.auditRejected,options.localRefineRounds);
			result.refinement=proposal.refinement;
			result.leftRebuild=proposal.leftUpdate==null ? null : proposal.leftUpdate.rebuildSummary();
			result.rightRebuild=proposal.rightUpdate==null ? null : proposal.rightUpdate.rebuildSummary();
			if(proposal.audit!=null){
				final FamilyFission.Audit audit=proposal.audit;
				if(audit.available()){
					Probe.verifyMembers(parent,audit.left,audit.right);
					bindCandidate(audit.left,work.measurementLayer);bindCandidate(audit.right,work.measurementLayer);
				}
				result.audit=new Audit(true,audit.left,audit.right,null,audit.leftRebuild,audit.rightRebuild,null,
					audit.unavailableReason,audit.extraRebuilds,audit.extraAccepted,audit.extraRejected,audit.extraAlignments);
			}
			memberPostings+=proposal.postingsVisited;retainedEdges+=proposal.retainedEdges;unionEdges+=proposal.unionEdges;
			metrics.end(parent.size(),parent.size(),0,1,1);
			rebuildAlignments+=metrics.constructionAlignments()-before;
			if(proposal.refinement!=null){rebuildAlignments+=proposal.refinement.rebuildAlignments;}
			if(!proposal.accepted()){result.reason=proposal.status.name();return result;}
			if(proposal.left.familyId!=parent.familyId || proposal.right.familyId!=rightId ||
				proposal.left.size()<options.minChildSize || proposal.right.size()<options.minChildSize ||
				(long)proposal.left.size()+proposal.right.size()!=parent.size()){
				throw new IllegalStateException("A fission proposal must conserve its parent and use the reserved child identities");
			}
			// Accepted fissions have no audit payload; each detached child is attached exactly once.
			assert(proposal.audit==null) : "FamilyFission.Result keeps accepted children and rejected audit candidates mutually exclusive";
			bindCandidate(proposal.left,work.measurementLayer);bindCandidate(proposal.right,work.measurementLayer);
			if(constructionOnly){
				result.parent=parent;result.left=proposal.left;result.right=proposal.right;
				result.constructed=true;result.reason="live_proposal";return result;
			}
			if(!hasQuality(proposal.left,proposal.leftUpdate) || !hasQuality(proposal.right,proposal.rightUpdate)){
				result.reason="child_quality_unavailable";return result;
			}
			result.parent=parent;result.left=proposal.left;result.right=proposal.right;
			result.parentQ=quality[qualitySlot(parent)];result.leftQ=normalized(proposal.leftUpdate.quality());result.rightQ=normalized(proposal.rightUpdate.quality());
			score(result,parent.familyId,rightId,neighborWork);
			return result;
		}

		private Decision merge(final int index,final FamilySketchBuilder builder,final ConsensusRefiner refiner,
				final ClusteringMetrics metrics,final FamilyNeighborEvidence.Workspace neighborWork){
			final MutualNearestMerge.Pair pair=work.pairs[index];
			final Family a=work.families[pair.leftIndex],b=work.families[pair.rightIndex];
			final Decision result=new Decision(a.familyId,b.familyId);
			if(!constructionOnly && (!cacheQuality(a) || !cacheQuality(b))){result.reason="input_quality_unavailable";return result;}
			final IntList members=new IntList(Math.addExact(a.size(),b.size()));
			for(int i=0; i<a.size(); i++){members.add(a.memberIndex(i));}
			for(int i=0; i<b.size(); i++){members.add(b.memberIndex(i));}
			members.sort();
			Family merged=new Family(a.familyId,a.id,members.array[0],a.consensusUnsafe());
			for(int i=1; i<members.size; i++){
				if(members.array[i]<=members.array[i-1]){throw new IllegalStateException("Fusion inputs must have disjoint original memberships");}
				merged.addMember(members.array[i]);
			}
			metrics.begin("consensus",work.cycle,merged.size());
			ConsensusRefiner.Update update=refiner.updateAfterConcurrentGreedy(merged,sequences,memberGate,metrics,
				true,work.summarizeRebuild,work.auditRejected);
			result.unionRebuild=update.rebuildSummary();
			metrics.end(merged.size(),merged.size(),0,1,1);
			rebuildAlignments+=update.alignments();
			if(options.localRefineRounds>0){
				if(update.accepted()+update.unchanged()<1){
					final byte[] actual=update.rejectedProposal();
					if(actual==null){result.reason="initial_union_candidate_unavailable";return result;}
					merged.setConsensus(actual);
				}
				final CompetitiveRefinement.Result refined=CompetitiveRefinement.refine(merged,null,members.toArray(),
					sequences,builder,memberGate,refiner,options.minChildSize,options.localRefineRounds);
				result.refinement=refined.summary;
				rebuildAlignments+=refined.summary.rebuildAlignments;
				if(!refined.available()){result.reason=refined.unavailableReason;return result;}
				merged=refined.first;update=refined.firstUpdate;result.unionRebuild=update.rebuildSummary();
			}
			if(update.accepted()+update.unchanged()<1){
				result.reason="union_rebuild_rejected";
				if(work.auditRejected){
					final byte[] proposed=update.rejectedProposal();
					Family candidate=null;
					if(proposed!=null){
						candidate=new Family(merged.familyId,merged.id,merged.memberIndex(0),merged.representative());
						for(int i=1;i<merged.size();i++){candidate.addMember(merged.memberIndex(i));}
						candidate.setConsensus(proposed);
						candidate.refreshStats(sequences,builder); candidate.markSketchRebuilt();
						Probe.verifyMembers(candidate,a,b);
						bindCandidate(candidate,work.measurementLayer);
					}
					result.audit=new Audit(false,null,null,candidate,null,null,result.unionRebuild,
						candidate==null ? "union_proposal_unavailable" : null,0,0,0,0);
				}
				return result;
			}
			if(constructionOnly){
				merged.refreshStats(sequences,builder);merged.markSketchRebuilt();
				bindCandidate(merged,work.measurementLayer);
				result.parent=merged;result.left=a;result.right=b;
				result.constructed=true;result.reason="live_proposal";return result;
			}
			if(!hasQuality(merged,update)){result.reason="union_quality_unavailable";return result;}
			merged.refreshStats(sequences,builder); merged.markSketchRebuilt();
			bindCandidate(merged,work.measurementLayer);
			result.parent=merged;result.left=a;result.right=b;
			result.parentQ=normalized(update.quality());result.leftQ=quality[qualitySlot(a)];result.rightQ=quality[qualitySlot(b)];
			score(result,a.familyId,b.familyId,neighborWork);
			return result;
		}

		/** All three models use one external population; each child also sees the other child explicitly. */
		private void score(final Decision result,final int excludeA,final int excludeB,final FamilyNeighborEvidence.Workspace scratch){
			final FamilyNeighborEvidence.Result p=work.neighbors.evaluate(result.parent,excludeA,excludeB,null,scratch);
			final FamilyNeighborEvidence.Result a=work.neighbors.evaluate(result.left,excludeA,excludeB,result.right,scratch);
			final FamilyNeighborEvidence.Result b=work.neighbors.evaluate(result.right,excludeA,excludeB,result.left,scratch);
			neighborAlignments+=(long)p.attemptedAlignments+a.attemptedAlignments+b.attemptedAlignments;
			neighborPostings+=p.postingsVisited+a.postingsVisited+b.postingsVisited;
			if(!p.available || !a.available || !b.available){result.reason="neighbor_quality_unavailable";return;}
			result.parentN=p.quality;result.leftN=a.quality;result.rightN=b.quality;
			result.leftPair=a;result.rightPair=b;
			result.scored=true;
			if(work.calibration){result.reason="calibration_proposal";return;}
			final LocalFamilyObjective.Evidence parent=new LocalFamilyObjective.Evidence(result.parent.size(),result.parentQ,result.parentN);
			final LocalFamilyObjective.Evidence left=new LocalFamilyObjective.Evidence(result.left.size(),result.leftQ,result.leftN);
			final LocalFamilyObjective.Evidence right=new LocalFamilyObjective.Evidence(result.right.size(),result.rightQ,result.rightN);
			result.gain=work.pairs==null ? objective.splitGain(parent,left,right,work.temperature,work.slotCost) :
				objective.fusionGain(left,right,parent,work.temperature,work.slotCost);
			result.accepted=objective.accepts(result.gain);
			result.reason=result.accepted ? "local_gain_accepted" : "local_gain_below_margin";
		}

		/** Each old family belongs to one proposal per wave/pass, so lazy cache writes have unique owners. */
		private boolean cacheQuality(final Family family){
			final int slot=qualitySlot(family);
			if(qualityState[slot]!=UNKNOWN){return qualityState[slot]==VALID;}
			final byte[] consensus=family.consensusUnsafe(); scorer.resetConsensus(consensus);
			double identity=0,ratio=0;
			for(int i=0; i<family.size(); i++){
				qualityAlignments++;
				final TailAlignmentDecision decision=placement.decide(sequences[family.memberIndex(i)].enc,consensus,scorer);
				if(!decision.accepted()){qualityState[slot]=MISSING;return false;}
				identity+=decision.identity;ratio+=decision.scoreRatio;
			}
			publishQuality(family,normalized(ConsensusQualityReport.available(family.familyId,family.size(),identity,ratio)));
			return true;
		}
		final Work work;
		final HbmPositionModel.ConsensusScorer scorer=new HbmPositionModel.ConsensusScorer(memberGate.matrix);
		final TailAlignmentGate placement=new TailAlignmentGate(memberGate.mode,memberGate.denominator,memberGate.maxLength,1,0,0,null,memberGate.matrix);
		int done;
		long qualityAlignments,neighborAlignments,rebuildAlignments;
		long memberPostings,neighborPostings,retainedEdges,unionEdges;
		boolean success;
	}

	private static boolean hasQuality(final Family family,final ConsensusRefiner.Update update){
		final ConsensusQualityReport.Snapshot q=update.quality();
		if(q==null || q.familyId()!=family.familyId || q.members()!=family.size()){
			throw new IllegalStateException("Rebuild quality must bind the proposed family identity and complete member count");
		}
		return q.available();
	}
	private static double normalized(final ConsensusQualityReport.Snapshot q){
		if(q==null || !q.available()){throw new IllegalArgumentException("Missing cohesion cannot become a local objective value");}
		return .5*(q.identityMean()+Math.max(0,Math.min(1,q.scoreRatioMean())));
	}

	private static final class Decision{
		Decision(final int a,final int b){inputA=a;inputB=b;}
		final int inputA,inputB;
		Family parent,left,right;
		FamilyNeighborEvidence.Result leftPair,rightPair;
		Audit audit;
		CompetitiveRefinement.Summary refinement;
		ConsensusRefiner.RebuildSummary leftRebuild,rightRebuild,unionRebuild;
		double parentQ,leftQ,rightQ,parentN,leftN,rightN,gain;
		boolean scored,accepted,constructed;
		String reason="not_evaluated";
	}

	/** Request-local counterfactual models, separate from legal production replacements.
	 * Consumers must score and discard these models before retaining a batch result.
	 */
	public static final class Audit{
		private Audit(final boolean fission_,final Family left_,final Family right_,final Family union_,
				final ConsensusRefiner.RebuildSummary leftRebuild_,final ConsensusRefiner.RebuildSummary rightRebuild_,
				final ConsensusRefiner.RebuildSummary unionRebuild_,final String unavailableReason_,
				final long rebuilds_,final long accepted_,final long rejected_,final long alignments_){
			fission=fission_;left=left_;right=right_;union=union_;
			leftRebuild=leftRebuild_;rightRebuild=rightRebuild_;unionRebuild=unionRebuild_;
			unavailableReason=unavailableReason_;rebuilds=rebuilds_;accepted=accepted_;rejected=rejected_;alignments=alignments_;
			if((fission ? union!=null : left!=null || right!=null) ||
				rebuilds<0 || accepted<0 || rejected<0 || alignments<0 || accepted>rebuilds || rejected>rebuilds-accepted ||
				available()!=(unavailableReason==null)){
				throw new IllegalStateException("Audit models and work counts must identify a complete detached replacement or explicit unavailability");
			}
		}
		public boolean available(){return fission ? left!=null && right!=null : union!=null;}
		public final boolean fission;
		public final Family left,right,union;
		public final ConsensusRefiner.RebuildSummary leftRebuild,rightRebuild,unionRebuild;
		public final String unavailableReason;
		public final long rebuilds,accepted,rejected,alignments;
	}

	/** One raw legal edit, or an explicit rejection with no production replacement/evidence payload.
	 * Feasibility means the native proposal/rebuild/evidence gates passed, never objective acceptance.
	 * Optional audit models are request-local and never make a rejected probe feasible.
	 */
	public static final class Probe{
		private Probe(final LayerCakeSchedule.Phase phase_,final int a,final int b,final int cost,final double t,final String reason_){
			phase=phase_;inputA=a;inputB=b;slotCost=cost;temperature=t;reason=reason_;parent=left=right=null;features=null;
			leftRebuild=rightRebuild=unionRebuild=null;
			audit=null;
			refinement=null;
			candidateFirst=candidateSecond=null;
		}
		private Probe(final LayerCakeSchedule.Phase phase_,final int a,final int b,final int cost,final double t,
				final Decision decision,final LocalFamilyObjective.Parameters parameters){
			phase=phase_;inputA=a;inputB=b;slotCost=cost;temperature=t;reason=decision.reason;
			parent=decision.scored ? decision.parent : null;
			left=decision.scored ? decision.left : null;right=decision.scored ? decision.right : null;
			leftRebuild=decision.leftRebuild; rightRebuild=decision.rightRebuild; unionRebuild=decision.unionRebuild;
			audit=decision.audit;
			refinement=decision.refinement;
			candidateFirst=decision.constructed ? (phase==LayerCakeSchedule.Phase.FISSION ? decision.left : decision.parent) :
				audit!=null && audit.available() ? (audit.fission ? audit.left : audit.union) : null;
			candidateSecond=candidateFirst==null ? null : decision.constructed ?
				(phase==LayerCakeSchedule.Phase.FISSION ? decision.right : null) : audit.right;
			features=decision.scored ? new ProbeFeatures(phase,decision,cost,t,parameters) : null;
			if(features!=null && audit!=null){throw new IllegalStateException("Feasible controls reuse ordinary models and cannot retain rejected audit candidates");}
			if(features!=null){verifyMembers(parent,left,right);}
			if(decision.constructed){verifyMembers(decision.parent,decision.left,decision.right);}
		}
		public boolean feasible(){return features!=null;}
		public Family[] replacements(){
			if(!feasible()){throw new IllegalStateException("A rejected probe has no legal replacement population");}
			return phase==LayerCakeSchedule.Phase.FISSION ? new Family[]{left,right} : new Family[]{parent};
		}
		public int[] removals(){return inputB<0 ? new int[]{inputA} : new int[]{inputA,inputB};}
		/** Count conservation alone cannot prove that the same original tails survived the edit. */
		private static void verifyMembers(final Family parent,final Family left,final Family right){
			final int[] before=parent.memberIndexes(),after=new int[Math.addExact(left.size(),right.size())];
			if(before.length!=after.length){throw new IllegalStateException("A probe must conserve the affected original-tail count");}
			for(int i=0; i<left.size(); i++){after[i]=left.memberIndex(i);}
			for(int i=0; i<right.size(); i++){after[left.size()+i]=right.memberIndex(i);}
			Arrays.sort(before);Arrays.sort(after);
			if(!Arrays.equals(before,after)){throw new IllegalStateException("A probe must conserve the exact affected original-tail ordinals");}
			for(int i=1; i<before.length; i++){
				if(before[i]<=before[i-1]){throw new IllegalStateException("A probe cannot duplicate any affected original-tail ordinal");}
			}
		}
		public final LayerCakeSchedule.Phase phase;
		public final int inputA,inputB,slotCost;
		public final double temperature;
		public final String reason;
		public final Family parent,left,right;
		public final ProbeFeatures features;
		public final Audit audit;
		public final CompetitiveRefinement.Summary refinement;
		/** Actual live candidates, independent of old objective/strict-veto eligibility. */
		public final Family candidateFirst,candidateSecond;
		public final ConsensusRefiner.RebuildSummary leftRebuild,rightRebuild,unionRebuild;
	}

	/** Unweighted predictors, oriented as gains for the requested operation. No provisional LQ/LS/LP or margin is used. */
	public static final class ProbeFeatures{
		private ProbeFeatures(final LayerCakeSchedule.Phase phase,final Decision d,final int cost,final double t,
				final LocalFamilyObjective.Parameters parameters){
			parent=new LocalFamilyObjective.Evidence(d.parent.size(),d.parentQ,d.parentN);
			left=new LocalFamilyObjective.Evidence(d.left.size(),d.leftQ,d.leftN);
			right=new LocalFamilyObjective.Evidence(d.right.size(),d.rightQ,d.rightN);
			leftPair=d.leftPair;rightPair=d.rightPair;
			final double sign=phase==LayerCakeSchedule.Phase.FISSION ? 1 : -1;
			cohesionDelta=sign*(left.members*left.cohesion+right.members*right.cohesion-parent.members*parent.cohesion);
			separationDelta=sign*(left.members*(left.cohesion-left.neighbor)+right.members*(right.cohesion-right.neighbor)-
				parent.members*(parent.cohesion-parent.neighbor));
			balance=1-Math.abs((long)left.members-right.members)/(double)parent.members;
			sizePressure=sign*t*parent.members*parameters.sizeForm.pressure(parent.members,parameters.largeBoundary)*balance;
			slotDelta=-sign*cost;
			if(!Double.isFinite(cohesionDelta) || !Double.isFinite(separationDelta) || !Double.isFinite(sizePressure)){
				throw new IllegalArgumentException("Raw probe predictors overflowed; use a finite bounded temperature and size form");
			}
		}
		public final LocalFamilyObjective.Evidence parent,left,right;
		public final FamilyNeighborEvidence.Result leftPair,rightPair;
		public final double cohesionDelta,separationDelta,sizePressure,slotDelta,balance;
	}

	/** Explicit controls for the opt-in implementation; ordinary clustering never constructs this engine. */
	public static final class Options{
		public Options(final int threads_,final int indexThreads_,final int fissionMinSize_,final int minChildSize_,
				final int memberNeighbors_,final int memberMinHits_,final int evidenceNeighbors_,final int mergeShortlist_,
				final double fusionSizeRatio_,final boolean consensusTrimLeft_,final float trimDepth_,final int familyBudget_){
			this(threads_,indexThreads_,fissionMinSize_,minChildSize_,memberNeighbors_,memberMinHits_,evidenceNeighbors_,
				mergeShortlist_,fusionSizeRatio_,consensusTrimLeft_,trimDepth_,familyBudget_,0);
		}
		public Options(final int threads_,final int indexThreads_,final int fissionMinSize_,final int minChildSize_,
				final int memberNeighbors_,final int memberMinHits_,final int evidenceNeighbors_,final int mergeShortlist_,
				final double fusionSizeRatio_,final boolean consensusTrimLeft_,final float trimDepth_,final int familyBudget_,
				final int localRefineRounds_){
			this(threads_,indexThreads_,fissionMinSize_,minChildSize_,memberNeighbors_,memberMinHits_,evidenceNeighbors_,
				mergeShortlist_,fusionSizeRatio_,consensusTrimLeft_,trimDepth_,familyBudget_,localRefineRounds_,
				DEFAULT_PRICED_FRACTION);
		}
		public Options(final int threads_,final int indexThreads_,final int fissionMinSize_,final int minChildSize_,
				final int memberNeighbors_,final int memberMinHits_,final int evidenceNeighbors_,final int mergeShortlist_,
				final double fusionSizeRatio_,final boolean consensusTrimLeft_,final float trimDepth_,final int familyBudget_,
				final int localRefineRounds_,final double pricedFraction_){
			if(threads_<1 || indexThreads_<1 || indexThreads_>threads_ || fissionMinSize_<1 || minChildSize_<1 || familyBudget_<1 ||
				memberNeighbors_<1 || memberMinHits_<1 || evidenceNeighbors_<1 || mergeShortlist_<1 || mergeShortlist_==Integer.MAX_VALUE ||
				!Double.isFinite(fusionSizeRatio_) || fusionSizeRatio_<1 || !Float.isFinite(trimDepth_) || trimDepth_<0 || trimDepth_>1 ||
				localRefineRounds_<0 || localRefineRounds_>2 ||
				!Double.isFinite(pricedFraction_) || pricedFraction_<=0 || pricedFraction_>1){
				throw new IllegalArgumentException("Invalid structural worker, proposal, shortlist, size-ratio, trimming or priced-rank controls");
			}
			threads=threads_;indexThreads=indexThreads_;fissionMinSize=fissionMinSize_;minChildSize=minChildSize_;
			memberNeighbors=memberNeighbors_;memberMinHits=memberMinHits_;evidenceNeighbors=evidenceNeighbors_;
			mergeShortlist=mergeShortlist_;fusionSizeRatio=fusionSizeRatio_;consensusTrimLeft=consensusTrimLeft_;trimDepth=trimDepth_;familyBudget=familyBudget_;
			localRefineRounds=localRefineRounds_;pricedFraction=pricedFraction_;
		}
		public final int threads,indexThreads,fissionMinSize,minChildSize,memberNeighbors,memberMinHits,evidenceNeighbors,mergeShortlist;
		public final int familyBudget;
		public final int localRefineRounds;
		public final double fusionSizeRatio,pricedFraction;
		public final boolean consensusTrimLeft;
		public final float trimDepth;
		public static final double DEFAULT_PRICED_FRACTION=0.99;
	}

	private final ProteinSequence[] sequences;
	private final FamilySketchBuilder sketches;
	private final TailAlignmentGate memberGate,mergeGate;
	private final LocalFamilyObjective objective;
	private final Options options;
	private final ByteStreamWriter report;
	private final ByteBuilder reportRow=new ByteBuilder();
	private Family[] expected;
	private int nextId;
	private final boolean localOnly;
	private final boolean constructionOnly;
	private final int firstLocalChildId;
	/** Hashed stable-ID lookups, populated once and shared read-only by local forks. */
	private final IntHashMap2 activeIndexByFamilyId,qualityIndexByFamilyId;
	private boolean probing,scheduling,forked;
	/** Set before forks, then borrowed read-only by every request-local worker. */
	private MemberRankInputs memberInputs;
	private int memberLayer;
	private int probeSlotCost=-1;
	private FamilyNeighborEvidence probeNeighbors;
	/** Lazily allocated per fork; feature probes reuse scorer/index scratch across requests. */
	private HbmPositionModel.ConsensusScorer featureScorer;
	private TailAlignmentGate featurePlacement;
	private FamilyNeighborEvidence.Workspace featureWorkspace;
	private MutualNearestMerge.PairQueries probeFusion;
	private double[] quality=new double[0];
	private byte[] qualityState=new byte[0];
	public long fissionProposals,fusionProposals,acceptedFissions,acceptedFusions,qualityAlignments,neighborAlignments,rebuildAlignments,selectionAlignments;
	public long memberPostings,neighborPostings,fissionRetainedEdges,fissionUnionEdges;
	private static final byte UNKNOWN=0,VALID=1,MISSING=2;
	// A family proposal is already a coarse unit; one-at-a-time claims balance widely varying sizes.
	private static final int BATCH=1;
}
