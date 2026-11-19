package prot.cluster;

import java.util.Arrays;

import ddl.FamilySketchBuilder;
import prot.HbmPositionModel;
import prot.ProteinSequence;

/** Bounded competitive reassignment and rebuilding of detached local candidates.
 * Every original member stays in exactly one child. Relaxed placement ties and
 * failures preserve its prior child; strict rebuild vetoes retain only actual
 * proposed bytes. No baseline family, member list or cached quality is mutated.
 * @author Brian Bushnell, Fischl
 */
public final class CompetitiveRefinement{

	/** Performs one mandatory reassignment/rebuild, then at most one repeat if members moved.
	 * first/second are actual initial candidate consensuses, not stale seeds. A null
	 * second requests fusion refinement. Original members must be exactly the unique
	 * union of the supplied memberships. The caller retains that original denominator
	 * for final pool evaluation even when strict gates veto some members.
	 */
	public static Result refine(final Family first, final Family second, final int[] originalMembers,
			final ProteinSequence[] sequences, final FamilySketchBuilder sketches,
			final TailAlignmentGate strictGate, final ConsensusRefiner worker,
			final int minChildSize, final int maxRounds){
		require(first!=null && originalMembers!=null && originalMembers.length>0 && sequences!=null &&
			sketches!=null && strictGate!=null && worker!=null, "Refinement needs candidates and a complete frozen member store");
		require(maxRounds>=1 && maxRounds<=2 && minChildSize>0 && worker.rounds()==1,
			"Competitive refinement is bounded to one or two rounds with one rebuild per nonempty child");
		require(strictGate.mode!=TailAlignmentGate.Mode.FULL_QUERY_UNPINNED &&
			(second==null || first.familyId!=second.familyId), "Refinement requires pinned geometry and distinct child IDs");
		final int[] members=originalMembers.clone(); Arrays.sort(members);
		for(int i=0; i<members.length; i++){
			require(members[i]>=0 && members[i]<sequences.length && sequences[members[i]]!=null &&
				sequences[members[i]].enc.length>0 && (i==0 || members[i]!=members[i-1]),
				"Original refinement members must be unique valid input ordinals");
		}
		byte[] prior=new byte[members.length], next=new byte[members.length]; Arrays.fill(prior, (byte)-1);
		mark(first, (byte)0, members, prior);
		if(second!=null){mark(second, (byte)1, members, prior);}
		for(byte side:prior){require(side>=0, "Candidate memberships must account for every original member");}
		final Work counts=new Work();
		if(first.size()==0 || second!=null && second.size()==0){return unavailable("empty_child", null, null, counts);}
		if(first.size()<minChildSize || second!=null && second.size()<minChildSize){
			return unavailable("below_child_floor", null, null, counts);
		}
		if(!usable(first, strictGate) || second!=null && !usable(second, strictGate)){
			return unavailable("initial_consensus_unavailable", null, null, counts);
		}
		final TailAlignmentGate placement=new TailAlignmentGate(strictGate.mode, strictGate.denominator,
			strictGate.maxLength, 1, 0, 0, null, strictGate.matrix);
		final HbmPositionModel.ConsensusScorer aScorer=new HbmPositionModel.ConsensusScorer(strictGate.matrix);
		final HbmPositionModel.ConsensusScorer bScorer=new HbmPositionModel.ConsensusScorer(strictGate.matrix);
		Family a=first, b=second;
		ConsensusRefiner.Update aUpdate=null, bUpdate=null;
		for(int round=0; round<maxRounds; round++){
			counts.rounds++;
			aScorer.resetConsensus(a.consensusUnsafe());
			if(b!=null){bScorer.resetConsensus(b.consensusUnsafe());}
			int aSize=0, bSize=0, moved=0;
			for(int i=0; i<members.length; i++){
				final byte[] query=sequences[members[i]].enc;
				final TailAlignmentDecision left=placement.decide(query, a.consensusUnsafe(), aScorer);
				counts.assignmentAlignments++;
				final TailAlignmentDecision right=b==null ? null : placement.decide(query, b.consensusUnsafe(), bScorer);
				if(b!=null){counts.assignmentAlignments++;}
				byte side=prior[i];
				if(left.accepted() && right!=null && right.accepted()){
					// Equal placeholder IDs disable only the final ID tie-break. Exact
					// native metric ties preserve membership rather than draining a child.
					final int compare=TailAlignmentDecision.compareAccepted(left, 0, right, 0);
					if(compare!=0){side=(byte)(compare<0 ? 0 : 1);}
				}else if(left.accepted()){side=0;}
				else if(right!=null && right.accepted()){side=1;}
				next[i]=side; moved+=(side==prior[i] ? 0 : 1);
				if(side==0){aSize++;}else{bSize++;}
			}
			counts.moved+=moved;
			assert(aSize+bSize==members.length) : "Relaxed placement cannot discard an original construction member";
			if(aSize==0 || b!=null && bSize==0){return unavailable("empty_child", aUpdate, bUpdate, counts);}
			if(aSize<minChildSize || b!=null && bSize<minChildSize){
				return unavailable("below_child_floor", aUpdate, bUpdate, counts);
			}
			final Family newA=assigned(a, (byte)0, members, next, sequences);
			final Family newB=b==null ? null : assigned(b, (byte)1, members, next, sequences);
			try(ClusteringMetrics metrics=new ClusteringMetrics(null, null, true, 1, Long.MAX_VALUE)){
				metrics.begin("consensus", round, members.length);
				aUpdate=worker.updateAfterConcurrentGreedy(newA, sequences, strictGate, metrics, true, true, true);
				bUpdate=newB==null ? null : worker.updateAfterConcurrentGreedy(newB, sequences, strictGate, metrics, true, true, true);
				counts.rebuildAlignments+=metrics.constructionAlignments();
				metrics.end(members.length, members.length, 0, newB==null ? 1 : 2, newB==null ? 1 : 2);
			}
			if(!adopt(newA, aUpdate, strictGate)){
				return unavailable("first_proposal_unavailable", aUpdate, bUpdate, counts);
			}
			if(newB!=null && !adopt(newB, bUpdate, strictGate)){
				return unavailable("second_proposal_unavailable", aUpdate, bUpdate, counts);
			}
			a=newA; b=newB;
			final byte[] swap=prior; prior=next; next=swap;
			if(moved==0){break;}
		}
		final FamilySketchBuilder local=sketches.fork();
		a.refreshStats(sequences, local); a.markSketchRebuilt();
		if(b!=null){b.refreshStats(sequences, local); b.markSketchRebuilt();}
		aScorer.resetConsensus(a.consensusUnsafe());
		if(b!=null){bScorer.resetConsensus(b.consensusUnsafe());}
		for(int ordinal:members){
			final TailAlignmentDecision left=strictGate.decide(sequences[ordinal].enc, a.consensusUnsafe(), aScorer);
			counts.assignmentAlignments++;
			final TailAlignmentDecision right=b==null ? null : strictGate.decide(sequences[ordinal].enc, b.consensusUnsafe(), bScorer);
			if(b!=null){counts.assignmentAlignments++;}
			if(!left.accepted() && (right==null || !right.accepted())){counts.forcedRetention++;}
		}
		return new Result(a, b, aUpdate, bUpdate, null, counts);
	}

	/** Checks membership provenance before any model construction or alignment. */
	private static void mark(final Family model, final byte side, final int[] members, final byte[] assignment){
		for(int i=0; i<model.size(); i++){
			final int at=Arrays.binarySearch(members, model.memberIndex(i));
			require(at>=0 && assignment[at]<0, "Candidate memberships must be disjoint and drawn only from the original member list");
			assignment[at]=side;
		}
	}

	/** Fresh membership storage prevents writes through Family's sharing copy constructor. */
	private static Family assigned(final Family previous, final byte side, final int[] members,
			final byte[] assignment, final ProteinSequence[] sequences){
		Family result=null;
		for(int i=0; i<members.length; i++){
			if(assignment[i]!=side){continue;}
			if(result==null){result=new Family(previous.familyId, previous.id, members[i], sequences[members[i]].enc);}
			else{result.addMember(members[i]);}
		}
		assert(result!=null) : "Empty assignment sides must be rejected before rebuilding";
		result.setConsensus(previous.consensusUnsafe());
		return result;
	}

	/** A strict veto may retain an actual valid proposed sequence, never the previous seed. */
	private static boolean adopt(final Family candidate, final ConsensusRefiner.Update update, final TailAlignmentGate gate){
		assert(update!=null) : "Every candidate must have a completed native rebuild outcome";
		if(update.accepted()+update.unchanged()==0){
			final byte[] proposed=update.rejectedProposal();
			if(proposed==null || proposed.length<gate.minSpan){return false;}
			candidate.setConsensus(proposed);
		}
		return usable(candidate, gate);
	}
	private static boolean usable(final Family model, final TailAlignmentGate gate){
		return model.consensusUnsafe()!=null && model.consensusUnsafe().length>=gate.minSpan;
	}
	private static Result unavailable(final String reason, final ConsensusRefiner.Update first,
			final ConsensusRefiner.Update second, final Work counts){
		return new Result(null, null, first, second, reason, counts);
	}

	/** Detached final models and native final-side rebuild outcomes; availability is not strict acceptance. */
	public static final class Result{
		private Result(final Family first_, final Family second_, final ConsensusRefiner.Update firstUpdate_,
				final ConsensusRefiner.Update secondUpdate_, final String reason, final Work counts){
			assert((reason==null)==(first_!=null)) : "Unavailable refinement cannot expose a stale model as a usable replacement";
			first=first_; second=second_; firstUpdate=firstUpdate_; secondUpdate=secondUpdate_;
			unavailableReason=reason; summary=new Summary(counts);
		}
		public boolean available(){return unavailableReason==null;}
		/** A valid rejected proposal stays available for audit but does not silently pass strict rebuilding. */
		public boolean strictAccepted(){
			return available() && firstUpdate.accepted()+firstUpdate.unchanged()>0 &&
				(second==null || secondUpdate.accepted()+secondUpdate.unchanged()>0);
		}
		public final Family first, second;
		public final ConsensusRefiner.Update firstUpdate, secondUpdate;
		public final String unavailableReason;
		public final Summary summary;
	}

	/** Scalar work only: rounds attempted, cumulative moves, and final strict-both-vetoed members.
	 * forcedRetention is meaningful only for an available Result. Assignment work
	 * includes the final strict check, which never moves or drops any member.
	 */
	public static final class Summary{
		private Summary(final Work counts){
			rounds=counts.rounds; moved=counts.moved; forcedRetention=counts.forcedRetention;
			assignmentAlignments=counts.assignmentAlignments; rebuildAlignments=counts.rebuildAlignments;
			assert(rounds>=0 && rounds<=2 && moved>=0 && forcedRetention>=0 &&
				assignmentAlignments>=0 && rebuildAlignments>=0) : "Bounded refinement reports finite nonnegative work counts";
		}
		public final int rounds;
		public final long moved, forcedRetention, assignmentAlignments, rebuildAlignments;
	}
	private static final class Work{
		int rounds;
		long moved, forcedRetention, assignmentAlignments, rebuildAlignments;
	}
	private static void require(final boolean valid, final String message){if(!valid){throw new IllegalArgumentException(message);}}
	private CompetitiveRefinement(){}
}
