package prot.cluster;

import java.util.Arrays;
import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import prot.ProteinSequence;
import prot.HbmPositionModel;
import structures.ByteBuilder;

/** Transactional consensus proposals: construction never filters or removes members.
 * Every current member must pass against both the current and proposed consensus.
 * A failing proposal leaves the published consensus intact and records its failures.
 * @author Keqing
 */
public final class ConsensusRefiner implements GreedyClusterer.ConsensusUpdater, AutoCloseable {

	public ConsensusRefiner(int rounds, boolean trimLeft, float trimDepthFraction, String events, boolean overwrite){
		if(rounds<1 || !Float.isFinite(trimDepthFraction) || trimDepthFraction<0 || trimDepthFraction>1){
			throw new IllegalArgumentException("Consensus rounds must be positive and trim fraction in [0,1]");
		}
		this.rounds=rounds; this.trimLeft=trimLeft; this.trimDepthFraction=trimDepthFraction;
		writer=events==null ? null : new ByteStreamWriter(events, overwrite, false, true);
		if(writer!=null){
			writer.start(); writer.println("proposal\tfamily_id\tmember_index\tstatus\treason\told_length\tnew_length");
		}
	}

	/** A worker owns its scorer and event buffer; it never writes the shared output stream. */
	private ConsensusRefiner(ConsensusRefiner parent){
		this(parent, parent.rounds);
	}

	/** A worker with a stage-local round count still reports through its parent driver. */
	private ConsensusRefiner(ConsensusRefiner parent,int rounds_){
		if(rounds_<1){throw new IllegalArgumentException("Consensus rounds must be positive");}
		rounds=rounds_; trimLeft=parent.trimLeft; trimDepthFraction=parent.trimDepthFraction;
		writer=null; worker=true; recordEvents=parent.writer!=null;
	}

	public ConsensusRefiner fork(){
		if(closed || worker){throw new IllegalStateException("Fork an open driver refiner only");}
		return new ConsensusRefiner(this);
	}

	/** Creates a worker with a stage-local round count and the same event-capture contract. */
	public ConsensusRefiner fork(int rounds){
		if(closed || worker){throw new IllegalStateException("Fork an open driver refiner only");}
		return new ConsensusRefiner(this, rounds);
	}

	/** Creates an independent writerless driver with identical refinement settings. */
	public ConsensusRefiner freshDriver(){
		return freshDriver(rounds);
	}

	/** Creates an independent writerless driver with identical trimming and a new round count. */
	public ConsensusRefiner freshDriver(int rounds){
		if(closed || worker){throw new IllegalStateException("Copy an open driver refiner only");}
		return new ConsensusRefiner(rounds, trimLeft, trimDepthFraction, null, true);
	}

	public int rounds(){return rounds;}

	/** Mutates only the worker-owned family, returning counters/events for ordered publication. */
	Update updateFamily(Family family,ProteinSequence[] sequences,TailAlignmentGate gate,ClusteringMetrics metrics){
		return updateFamily(family,sequences,gate,metrics,false);
	}

	Update updateFamily(Family family,ProteinSequence[] sequences,TailAlignmentGate gate,ClusteringMetrics metrics,
			boolean capture){
		return updateFamily(family,family.membershipView(),sequences,gate,metrics,false,false,capture);
	}

	/** Constructs one online proposal from exactly one published membership prefix.
	 * The caller owns this unpublished family copy and later attempts its version CAS.
	 * Later joins share membership but neither enter this proposal nor veto it.
	 */
	public Update updateOnline(Family family,Family.MembershipView view,ProteinSequence[] sequences,
			TailAlignmentGate gate,ClusteringMetrics metrics){
		if(gate.mode==TailAlignmentGate.Mode.FULL_QUERY_UNPINNED){
			throw new IllegalArgumentException("Concurrent greedy construction requires right-pinned placement");
		}
		return updateFamily(family,view,sequences,gate,metrics,true,true,false);
	}

	/** The first rebuild after asynchronous admission can still contain late joins.
	 * It uses inclusive placement for construction until global reassignment restores
	 * the strict old-model membership invariant; all new proposals remain strictly gated.
	 */
	public Update updateAfterConcurrentGreedy(Family family,ProteinSequence[] sequences,
			TailAlignmentGate gate,ClusteringMetrics metrics){
		return updateAfterConcurrentGreedy(family,sequences,gate,metrics,false);
	}

	public Update updateAfterConcurrentGreedy(Family family,ProteinSequence[] sequences,
			TailAlignmentGate gate,ClusteringMetrics metrics,boolean capture){
		return updateAfterConcurrentGreedy(family,sequences,gate,metrics,capture,false);
	}

	public Update updateAfterConcurrentGreedy(Family family,ProteinSequence[] sequences,
			TailAlignmentGate gate,ClusteringMetrics metrics,boolean capture,boolean summarizeRebuild){
		return updateAfterConcurrentGreedy(family,sequences,gate,metrics,capture,summarizeRebuild,false);
	}

	/** Optionally retains a structurally valid proposal vetoed by strict member checks.
	 * The published family and ordinary quality/counters retain their normal semantics.
	 */
	public Update updateAfterConcurrentGreedy(Family family,ProteinSequence[] sequences,
			TailAlignmentGate gate,ClusteringMetrics metrics,boolean capture,boolean summarizeRebuild,
			boolean captureRejectedProposal){
		if(gate.mode==TailAlignmentGate.Mode.FULL_QUERY_UNPINNED){
			throw new IllegalArgumentException("Concurrent greedy construction requires right-pinned placement");
		}
		return updateFamily(family,family.membershipView(),sequences,gate,metrics,false,true,capture,summarizeRebuild,captureRejectedProposal);
	}

	private Update updateFamily(Family family,Family.MembershipView view,ProteinSequence[] sequences,
			TailAlignmentGate gate,ClusteringMetrics metrics,boolean online,boolean inclusivePlacement,boolean capture){
		return updateFamily(family,view,sequences,gate,metrics,online,inclusivePlacement,capture,false,false);
	}

	private Update updateFamily(Family family,Family.MembershipView view,ProteinSequence[] sequences,
			TailAlignmentGate gate,ClusteringMetrics metrics,boolean online,boolean inclusivePlacement,
			boolean capture,boolean summarizeRebuild,boolean captureRejectedProposal){
		if(!worker){throw new IllegalStateException("Detached updates require a private worker refiner");}
		if(summarizeRebuild && rounds!=1){throw new IllegalStateException("Rebuild summaries are single-proposal diagnostics");}
		proposals=accepted=rejected=unchanged=rejectedMemberChecks=incompatibleSnapshots=0;
		events=recordEvents ? new ArrayList<Event>() : null;
		summary=null; summarize=summarizeRebuild;
		rejectedProposal=null; retainRejectedProposal=captureRejectedProposal;
		try{
			final long before=metrics.constructionAlignments();
			final ConsensusQualityReport.Snapshot quality=update(family,view,sequences,gate,metrics,
				online,inclusivePlacement,capture);
			final RebuildSummary rebuildSummary=summary==null ? null : summary.toSummary();
			return new Update(proposals,accepted,rejected,unchanged,rejectedMemberChecks,
				metrics.constructionAlignments()-before,incompatibleSnapshots,events,quality,rebuildSummary,
				rejectedProposal);
		}finally{
			events=null; summary=null; summarize=false;
			rejectedProposal=null; retainRejectedProposal=false;
		}
	}

	/** Replays local proposal IDs after the preceding families' proposals, preserving serial bytes. */
	void merge(Update update){
		if(closed || worker || update==null){throw new IllegalStateException("Only the driver merges completed family updates");}
		final long offset=proposals;
		if(update.events!=null){
			for(Event event : update.events){writeEvent(Math.addExact(offset,event.proposal),event);}
		}
		proposals=Math.addExact(proposals,update.proposals);
		accepted=Math.addExact(accepted,update.accepted);
		rejected=Math.addExact(rejected,update.rejected);
		unchanged=Math.addExact(unchanged,update.unchanged);
		rejectedMemberChecks=Math.addExact(rejectedMemberChecks,update.rejectedMemberChecks);
		incompatibleSnapshots=Math.addExact(incompatibleSnapshots,update.incompatibleSnapshots);
	}

	@Override
	public void update(Family family, ProteinSequence[] sequences, TailAlignmentGate gate,
			ClusteringMetrics metrics, boolean online){
		if(family==null){throw new IllegalStateException("Consensus refinement needs a family");}
		update(family,family.membershipView(),sequences,gate,metrics,online,false,false);
	}

	public ConsensusQualityReport.Snapshot updateObserved(Family family, ProteinSequence[] sequences,
			TailAlignmentGate gate, ClusteringMetrics metrics, boolean online){
		if(family==null){throw new IllegalStateException("Consensus refinement needs a family");}
		return update(family,family.membershipView(),sequences,gate,metrics,online,false,true);
	}

	private ConsensusQualityReport.Snapshot update(Family family,Family.MembershipView view,
			ProteinSequence[] sequences,TailAlignmentGate gate,ClusteringMetrics metrics,
			boolean online,boolean inclusivePlacement,boolean capture){
		if(closed || family==null || view==null || !view.viewOf(family) || view.size()==0){
			throw new IllegalStateException("Consensus refinement needs an open accumulator and assigned members");
		}
		final int size=view.size();
		// Construction admits positive traces, not new family members. Only the
		// proposed NEW model uses the configured gate, for every member in this S.
		final TailAlignmentGate placement=inclusivePlacement ? new TailAlignmentGate(gate.mode,
			gate.denominator,gate.maxLength,1,0,0,null,gate.matrix) : gate;
		final int limit=online ? 1 : rounds;
		ConsensusQualityReport.Snapshot quality=null;
		for(int round=0; round<limit; round++){
			final byte[] old=family.consensusUnsafe();
			if(summarize){summary=new RebuildSummaryBuilder(family.familyId,size,old.length);}
			final FamilyConsensus accumulator=new FamilyConsensus(old, gate.mode, trimLeft, trimDepthFraction, gate.matrix);
			double oldIdentitySum=0, oldScoreRatioSum=0;
			for(int i=0; i<size; i++){
				final int input=view.memberIndex(i);
				metrics.alignment(true);
				final TailAlignmentDecision decision=decide(placement,sequences[input].enc,old);
				if(!accumulator.add(sequences[input].enc, decision)){
					if(inclusivePlacement){
						assert(decision.reason==TailAlignmentDecision.Reason.NO_ALIGNMENT ||
							placement.mode==TailAlignmentGate.Mode.PIN_RIGHT_OVERHANG_ONLY &&
							(decision.reason==TailAlignmentDecision.Reason.NONPOSITIVE_SCORE && decision.score<=0 ||
							 decision.score<0 &&
							 (decision.reason==TailAlignmentDecision.Reason.QUERY_SPAN_SHORT && decision.querySpan==0 ||
							  decision.reason==TailAlignmentDecision.Reason.REFERENCE_SPAN_SHORT && decision.referenceSpan==0))) :
							"Inclusive placement must admit every positive paired trace; rejected reason="+decision.reason;
						proposals++; rejected++; incompatibleSnapshots++;
						rebuildReason("snapshot_no_positive_placement",1);
						event(family.familyId,input,"rejected","snapshot_no_positive_placement",old.length,old.length);
						return capture ? ConsensusQualityReport.unavailable(family.familyId,size,
							"snapshot_no_positive_placement") : null;
					}
					throw new IllegalStateException("Published consensus no longer accepts assigned member "+input+
						" of "+family.id+": "+decision+"; no member was silently excluded");
				}
				if(capture){
					oldIdentitySum+=decision.identity;
					oldScoreRatioSum+=decision.scoreRatio;
				}
			}
			if(accumulator.acceptedMembers()!=size){
				throw new IllegalStateException("Every assigned member must contribute to consensus construction: "+family.id);
			}
			if(capture){quality=ConsensusQualityReport.available(family.familyId,size,oldIdentitySum,oldScoreRatioSum);}
			proposals++;
			final byte[] proposed;
			try{proposed=accumulator.consensus();}
			catch(FamilyConsensus.EmptyConsensusException empty){
				rejected++;
				proposedLength(0);
				rebuildReason("consensus_empty_after_trim",1);
				event(family.familyId, -1, "rejected", "consensus_empty_after_trim", old.length, 0);
				return quality;
			}
			proposedLength(proposed.length);
			if(Arrays.equals(old, proposed) && !inclusivePlacement){
				unchanged++;
				rebuildReason("same_sequence",1);
				event(family.familyId, -1, "unchanged", "same_sequence", old.length, proposed.length);
				return quality;
			}
			long failed=0;
			double newIdentitySum=0, newScoreRatioSum=0;
			if(proposed.length<gate.minSpan){
				failed=size;
				rebuildReason("consensus_span_short",failed);
				event(family.familyId, -1, "rejected", "consensus_span_short", old.length, proposed.length);
			}else{
				for(int i=0; i<size; i++){
					final int input=view.memberIndex(i);
					metrics.alignment(true);
					final TailAlignmentDecision decision=decide(gate,sequences[input].enc,proposed);
					if(!decision.accepted()){
						failed++;
						rebuildReason(decision.reason.name(),1);
						event(family.familyId, input, "rejected", decision.reason.name(), old.length, proposed.length);
					}else if(capture){
						newIdentitySum+=decision.identity;
						newScoreRatioSum+=decision.scoreRatio;
					}
				}
			}
			if(failed>0){
				rejected++; rejectedMemberChecks+=failed;
				if(retainRejectedProposal && proposed.length>0 && proposed.length>=gate.minSpan){
					rejectedProposal=proposed.clone();
				}
				return quality;
			}
			if(Arrays.equals(old,proposed)){
				unchanged++;
				rebuildReason("same_sequence",1);
				event(family.familyId,-1,"unchanged","same_sequence",old.length,proposed.length);
				return quality;
			}
			family.setConsensus(proposed);
			if(capture){quality=ConsensusQualityReport.available(family.familyId,size,newIdentitySum,newScoreRatioSum);}
			accepted++;
			// Online concurrent callers count publication only after their version CAS wins.
			if(!(online && inclusivePlacement)){metrics.consensusUpdated();}
			rebuildReason("all_members_pass",1);
			event(family.familyId, -1, "accepted", "all_members_pass", old.length, proposed.length);
		}
		return quality;
	}

	private void event(int family, int member, String status, String reason, int oldLength, int newLength){
		assert(proposals>0) : "Consensus event IDs identify the proposal just constructed";
		if(worker){
			if(events!=null){events.add(new Event(proposals,family,member,status,reason,oldLength,newLength));}
		}else if(writer!=null){
			writer.println(new ByteBuilder().append(proposals).tab().append(family).tab().append(member).tab()
				.append(status).tab().append(reason).tab().append(oldLength).tab().append(newLength));
		}
	}

	private TailAlignmentDecision decide(TailAlignmentGate gate,byte[] query,byte[] reference){
		if(!worker){return gate.decide(query,reference);}
		// The gate arrives with the update; a worker never borrows another worker's scorer.
		if(scorer==null || !scorer.matrix.same(gate.matrix)){
			scorer=new HbmPositionModel.ConsensusScorer(gate.matrix);
		}
		scorer.resetConsensus(reference);
		return gate.decide(query,reference,scorer);
	}

	private void writeEvent(long proposal,Event event){
		if(writer!=null){
			writer.println(new ByteBuilder().append(proposal).tab().append(event.family).tab().append(event.member).tab()
				.append(event.status).tab().append(event.reason).tab().append(event.oldLength).tab().append(event.newLength));
		}
	}

	private void proposedLength(int length){
		if(summary!=null){summary.proposedLength=length;}
	}

	private void rebuildReason(String reason,long count){
		if(summary!=null){summary.add(reason,count);}
	}

	/** One completed family; no mutable worker scratch escapes with this result. */
	public static final class Update{
		Update(long proposals,long accepted,long rejected,long unchanged,long rejectedMemberChecks,
				long alignments,long incompatibleSnapshots,ArrayList<Event> events,
				ConsensusQualityReport.Snapshot quality,RebuildSummary rebuildSummary,byte[] rejectedProposal){
			this.proposals=proposals; this.accepted=accepted; this.rejected=rejected; this.unchanged=unchanged;
			this.rejectedMemberChecks=rejectedMemberChecks; this.alignments=alignments; this.events=events;
			this.incompatibleSnapshots=incompatibleSnapshots; this.quality=quality;
			this.rebuildSummary=rebuildSummary;
			this.rejectedProposal=rejectedProposal;
		}
		public long proposals(){return proposals;}
		public long accepted(){return accepted;}
		public long rejected(){return rejected;}
		public long unchanged(){return unchanged;}
		public long rejectedMemberChecks(){return rejectedMemberChecks;}
		public long alignments(){return alignments;}
		public long incompatibleSnapshots(){return incompatibleSnapshots;}
		public ConsensusQualityReport.Snapshot quality(){return quality;}
		public RebuildSummary rebuildSummary(){return rebuildSummary;}
		/** Actual strict-veto proposal, copied for the caller; null if disabled or unavailable. */
		public byte[] rejectedProposal(){return rejectedProposal==null ? null : rejectedProposal.clone();}
		private final long proposals,accepted,rejected,unchanged,rejectedMemberChecks,alignments,incompatibleSnapshots;
		private final ArrayList<Event> events;
		private final ConsensusQualityReport.Snapshot quality;
		private final RebuildSummary rebuildSummary;
		private final byte[] rejectedProposal;
	}

	public static final class RebuildSummary{
		private RebuildSummary(int familyId_,int members_,int oldLength_,int proposedLength_,
				long proposals_,long accepted_,long unchanged_,long rejected_,long rejectedMemberChecks_,
				long incompatibleSnapshots_,String[] reasons_,long[] reasonCounts_){
			familyId=familyId_; members=members_; oldLength=oldLength_; proposedLength=proposedLength_;
			proposals=proposals_; accepted=accepted_; unchanged=unchanged_; rejected=rejected_;
			rejectedMemberChecks=rejectedMemberChecks_; incompatibleSnapshots=incompatibleSnapshots_;
			reasons=reasons_; reasonCounts=reasonCounts_;
		}
		public ByteBuilder appendReasonCounts(ByteBuilder bb){
			if(reasons.length<1){return bb.append("none");}
			for(int i=0; i<reasons.length; i++){
				if(i>0){bb.append(';');}
				bb.append(reasons[i]).append(':').append(reasonCounts[i]);
			}
			return bb;
		}
		public final int familyId,members,oldLength,proposedLength;
		public final long proposals,accepted,unchanged,rejected,rejectedMemberChecks,incompatibleSnapshots;
		private final String[] reasons;
		private final long[] reasonCounts;
	}

	private final class RebuildSummaryBuilder{
		RebuildSummaryBuilder(int familyId_,int members_,int oldLength_){
			familyId=familyId_; members=members_; oldLength=oldLength_; proposedLength=-1;
		}
		void add(String reason,long count){
			if(count<1){throw new IllegalArgumentException("Rebuild reason counts must be positive");}
			for(int i=0; i<reasonCount; i++){
				if(reasons[i].equals(reason)){counts[i]=Math.addExact(counts[i],count); return;}
			}
			if(reasonCount>=reasons.length){
				reasons=Arrays.copyOf(reasons,reasons.length*2);
				counts=Arrays.copyOf(counts,counts.length*2);
			}
			reasons[reasonCount]=reason; counts[reasonCount]=count; reasonCount++;
		}
		RebuildSummary toSummary(){
			return new RebuildSummary(familyId,members,oldLength,proposedLength,proposals,
				accepted,unchanged,rejected,rejectedMemberChecks,incompatibleSnapshots,
				Arrays.copyOf(reasons,reasonCount),Arrays.copyOf(counts,reasonCount));
		}
		String[] reasons=new String[4];
		long[] counts=new long[4];
		final int familyId,members,oldLength;
		int proposedLength,reasonCount;
	}

	private static final class Event{
		Event(long proposal,int family,int member,String status,String reason,int oldLength,int newLength){
			this.proposal=proposal; this.family=family; this.member=member; this.status=status;
			this.reason=reason; this.oldLength=oldLength; this.newLength=newLength;
		}
		final long proposal;
		final int family,member,oldLength,newLength;
		final String status,reason;
	}

	@Override
	public void close(){
		if(closed){return;}
		closed=true;
		if(writer!=null && writer.poisonAndWait()){throw new IllegalStateException("Consensus proposal diagnostics failed to write");}
	}

	public long proposals(){return proposals;}
	public long accepted(){return accepted;}
	public long rejected(){return rejected;}
	public long unchanged(){return unchanged;}
	public long rejectedMemberChecks(){return rejectedMemberChecks;}
	public long incompatibleSnapshots(){return incompatibleSnapshots;}
	private final int rounds;
	private final boolean trimLeft;
	private final float trimDepthFraction;
	private final ByteStreamWriter writer;
	private boolean worker,recordEvents;
	private HbmPositionModel.ConsensusScorer scorer;
	private RebuildSummaryBuilder summary;
	private ArrayList<Event> events;
	private boolean summarize,retainRejectedProposal;
	private byte[] rejectedProposal;
	private long proposals, accepted, rejected, unchanged, rejectedMemberChecks,incompatibleSnapshots;
	private boolean closed;
}
