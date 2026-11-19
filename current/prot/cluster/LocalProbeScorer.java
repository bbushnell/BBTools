package prot.cluster;

import prot.HbmPositionModel;

/** Worker-local competition between the one or two replacement models of a local probe.
 * The original tail remains the alignment query. Shared-register hits may come
 * from the reverse tail index, but filtering and alignment retain their native
 * query/reference direction. Model and query bytes are borrowed read-only.
 * No baseline model, membership or index is changed by this scorer.
 * @author Brian Bushnell, Keqing
 */
public final class LocalProbeScorer{

	public LocalProbeScorer(final TailAlignmentGate gate_,final CandidateFilter filter_,
			final int width_,final int minimumLength_){
		require(gate_!=null && filter_!=null,"Local scoring requires the native alignment gate and candidate filter");
		require(minimumLength_>=gate_.minSpan,"Local query length floor cannot be shorter than the alignment span floor");
		filter_.validateWidth(width_);
		gate=gate_; filter=filter_; width=width_; minimumLength=minimumLength_;
		scorer=new HbmPositionModel.ConsensusScorer(gate.matrix);
	}

	/** Borrows the fixed replacements of one probe until endProbe. Each consensus
	 * matrix is initialized only when that model first passes the candidate filter.
	 * Models must remain immutable during the probe; a later probe may reuse changed
	 * Family objects. Callers without a probe binding retain per-candidate resetting.
	 */
	public void beginProbe(final Family first,final Family second){
		require(!probeBound,"End the previous local probe before binding new replacement models");
		checkModel(first);
		if(second!=null){checkModel(second);}
		require(second==null || first.familyId!=second.familyId,
			"Two probe replacements need distinct IDs for native tie breaking");
		firstReplacement=first; secondReplacement=second;
		firstReady=secondReady=false; probeBound=true;
		query=null; best=null;
	}

	/** Releases borrowed models and invalidates both matrices, also after a failure. */
	public void endProbe(){
		firstReplacement=secondReplacement=null;
		firstReady=secondReady=probeBound=false;
		query=null; best=null;
	}

	/** Begins one independent tail query, clearing the preceding winner and counters.
	 * False means the tail fails the length floor and must not be passed to consider.
	 * The caller keeps the encoded bytes immutable until the next prepare call.
	 */
	public boolean prepare(final byte[] encoded){
		query=null; best=null; winner=-1; eligible=false;
		considered=aligned=hitFiltered=hydroFiltered=0;
		require(encoded!=null,"Local scoring needs an encoded tail");
		if(encoded.length>=minimumLength){
			queryFilter.reset(encoded);
			eligible=true;
		}
		query=encoded;
		return eligible;
	}

	/** Reference occupancy determines the native reassignment hit floor.
	 * A zero floor requires the reverse index to enumerate zero-hit tails too.
	 */
	public int requiredHits(final Family replacement){
		checkModel(replacement);
		return filter.requiredHits(false,filter.passHitsFlex ? replacement.occupiedBins() : 0);
	}

	/** Considers a replacement once using its exact shared nonzero-register count.
	 * The two-model limit enforces the local-probe cost shape; baseline neighbors
	 * never enter this competition. Best-match ties use the existing native rule.
	 */
	public void consider(final Family replacement,final int exactHits){
		require(query!=null && eligible,"Prepare an eligible original tail before considering a replacement");
		checkModel(replacement);
		require(!probeBound || replacement==firstReplacement || replacement==secondReplacement,
			"A bound local probe may consider only its borrowed replacement models");
		require(exactHits>=0 && exactHits<=width,"Shared-register count must lie within the frozen sketch width");
		require(considered<2,"A local binary probe may score only its one or two replacement models");
		considered++;
		final int rejected=filter.rejection(false,exactHits,
			filter.passHitsFlex ? replacement.occupiedBins() : 0,queryFilter,replacement);
		if(rejected==CandidateFilter.LOW_HITS){hitFiltered++; return;}
		if(rejected==CandidateFilter.HYDRO){hydroFiltered++; return;}
		assert(rejected==CandidateFilter.PASS) : "Every candidate-screen reason must be handled before alignment";
		final byte[] reference=replacement.consensusUnsafe();
		aligned++;
		final HbmPositionModel.ConsensusScorer model=scorer(replacement,reference);
		final TailAlignmentDecision decision=gate.decide(query,reference,model);
		if(decision.accepted() && (best==null ||
			TailAlignmentDecision.compareAccepted(decision,replacement.familyId,best,winner)<0)){
			best=decision; winner=replacement.familyId;
		}
	}

	/** Stable family ID of the best accepted replacement, or -1 when none passed. */
	public int winner(){checkPrepared(); return winner;}
	public TailAlignmentDecision decision(){checkPrepared(); return best;}
	public int candidates(){checkPrepared(); return considered;}
	public int alignments(){checkPrepared(); return aligned;}
	public int lowHitRejections(){checkPrepared(); return hitFiltered;}
	public int hydroRejections(){checkPrepared(); return hydroFiltered;}

	/** Selects a worker-owned matrix after filtering. Bound models initialize once;
	 * unbound callers reset for every pair. The gate verifies the retained bytes.
	 */
	private HbmPositionModel.ConsensusScorer scorer(final Family replacement, final byte[] reference){
		if(!probeBound){scorer.resetConsensus(reference); return scorer;}
		if(replacement==firstReplacement){
			if(!firstReady){scorer.resetConsensus(reference); firstReady=true;}
			return scorer;
		}
		assert(replacement==secondReplacement) : "consider validates the model before selecting a probe matrix";
		if(secondScorer==null){secondScorer=new HbmPositionModel.ConsensusScorer(gate.matrix);}
		if(!secondReady){secondScorer.resetConsensus(reference); secondReady=true;}
		return secondScorer;
	}

	private void checkModel(final Family model){
		require(model!=null && model.familyId>=0 && model.sketchUnsafe()!=null &&
			model.sketchUnsafe().length==width,"Replacement model must match the frozen pool sketch geometry");
		require(model.consensusUnsafe()!=null && model.consensusUnsafe().length>=gate.minSpan,
			"Replacement consensus must satisfy the native reference span precondition");
	}
	private void checkPrepared(){require(query!=null,"Local scoring has no prepared query");}
	private static void require(final boolean valid,final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}

	private final TailAlignmentGate gate;
	private final CandidateFilter filter;
	private final int width,minimumLength;
	private final HbmPositionModel.ConsensusScorer scorer;
	private HbmPositionModel.ConsensusScorer secondScorer;
	private Family firstReplacement,secondReplacement;
	private final CandidateFilter.Context queryFilter=new CandidateFilter.Context();
	private byte[] query;
	private TailAlignmentDecision best;
	private int winner=-1,considered,aligned,hitFiltered,hydroFiltered;
	private boolean eligible;
	private boolean probeBound,firstReady,secondReady;
}
