package prot;

import java.util.Arrays;
import java.util.HashSet;

/** Highest-score eligible stop-tail family on a caller-supplied shortlist.
 * Construction gates are shared with HbmTailProfileRefiner; no full-protein
 * coverage/core threshold or calibrated detection threshold is substituted.
 * The library is detached and immutable. Each worker owns its Scratch.
 * @author Keqing
 */
public final class HbmTailCompetitiveGate {

	public HbmTailCompetitiveGate(String[] ids, byte[][] references, HbmPositionModel[] profiles){
		this(ids, references, profiles, HbmTailProfileRefiner.TailPolicy.LEGACY_CONSTRUCTION);
	}

	public HbmTailCompetitiveGate(String[] ids, byte[][] references, HbmPositionModel[] profiles,
			HbmTailProfileRefiner.TailPolicy policy_){
		require(ids!=null && references!=null && profiles!=null && ids.length>0
			&& ids.length==references.length && ids.length==profiles.length, "Tail library dimensions must agree");
		require(policy_!=null, "Tail competition needs a named acceptance policy");
		familyIds=ids.clone(); refs=new byte[ids.length][]; models=profiles.clone();
		final HashSet<String> seen=new HashSet<String>();
		for(int i=0; i<ids.length; i++){
			require(ids[i]!=null && !ids[i].isEmpty() && seen.add(ids[i]), "Tail family IDs must be nonempty and unique");
			for(int j=0; j<ids[i].length(); j++){
				require(ids[i].charAt(j)>32 && ids[i].charAt(j)<127, "Tail family IDs must be printable ASCII tokens for deterministic ties");
			}
			validateSequence(references[i]);
			require(profiles[i]!=null, "Missing frozen tail profile: "+ids[i]);
			refs[i]=references[i].clone(); models[i].requireConsensus(refs[i]);
		}
		policy=policy_;
	}

	/** Scores every candidate, then tests traces in score/lexical-family order.
	 * Returns -1 for an orphan, including an empty shortlist. A high-scoring
	 * ineligible family cannot hide a lower-scoring eligible family.
	 * All counters and prior winner evidence are reset for each query.
	 */
	public int select(byte[] query, int[] shortlist, int count, Scratch work){
		validateSequence(query);
		require(shortlist!=null && work!=null && count>=0 && count<=shortlist.length
			&& count<=work.order.length && work.seen.length==models.length, "Tail shortlist/scratch dimensions must agree");
		work.reset();
		int positive=0;
		for(int i=0; i<count; i++){
			final int rank=shortlist[i];
			require(rank>=0 && rank<models.length && work.seen[rank]!=work.generation, "Duplicate or out-of-range tail candidate: "+rank);
			work.seen[rank]=work.generation;
			final HbmPositionModel.TailResult result=models[rank].alignTail(query, HbmTailProfileRefiner.CAP, false);
			work.scored++;
			if(result==null){work.nonpositive++; continue;}
			int at=positive;
			while(at>0 && HbmCompetitiveGate.better(result.score, familyIds[rank], work.scores[at-1], familyIds[work.order[at-1]])){
				work.order[at]=work.order[at-1]; work.scores[at]=work.scores[at-1]; at--;
			}
			work.order[at]=rank; work.scores[at]=result.score; positive++;
		}
		assert(positive+work.nonpositive==count) : "Every supplied family must have one recorded score outcome";
		for(int i=0; i<positive; i++){
			final int rank=work.order[i];
			final HbmPositionModel.TailResult result=models[rank].alignTail(query, HbmTailProfileRefiner.CAP, true);
			work.recorded++;
			require(result!=null && result.score==work.scores[i], "Recording a frozen tail trace changed its optimum score");
			if(HbmTailProfileRefiner.rejection(query, refs[rank], result, policy)==null){
				work.winner=rank; work.score64=result.score; work.alignment=result;
				return rank;
			}
			work.rejected++;
		}
		return -1;
	}

	private static void validateSequence(byte[] sequence){
		require(sequence!=null && sequence.length>=HbmTailProfileRefiner.MIN_LENGTH
			&& sequence.length<=HbmTailProfileRefiner.CAP, "Tail construction sequences must contain100..200 residues");
		for(byte residue : sequence){require(residue>=0 && residue<20, "Tail construction excludes stop/ambiguous residues");}
	}
	private static void require(boolean condition, String message){if(!condition){throw new IllegalArgumentException(message);}}

	/** Caller-owned reusable state; alignment is nonnull only for an accepted winner. */
	public static final class Scratch {
		public Scratch(int families, int capacity){
			require(families>0 && capacity>0 && capacity<=families, "Tail scratch capacity must fit its family roster");
			seen=new int[families]; order=new int[capacity]; scores=new int[capacity];
		}
		private void reset(){
			if(generation==Integer.MAX_VALUE){Arrays.fill(seen, 0); generation=0;}
			generation++; winner=-1; score64=Integer.MIN_VALUE; alignment=null;
			scored=nonpositive=recorded=rejected=0;
		}
		public int winner=-1, score64=Integer.MIN_VALUE;
		public int scored, nonpositive, recorded, rejected;
		public HbmPositionModel.TailResult alignment;
		private final int[] seen, order, scores;
		private int generation;
	}
	private final String[] familyIds;
	private final byte[][] refs;
	private final HbmPositionModel[] models;
	private final HbmTailProfileRefiner.TailPolicy policy;
}
