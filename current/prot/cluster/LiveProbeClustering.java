package prot.cluster;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.ReadWriteLock;

import ddl.FamilySketchBuilder;
import ddl.TailPoolIndex;
import fileIO.ByteStreamWriter;
import prot.ProteinSequence;
import structures.ByteBuilder;
import structures.IntList;
import structures.IntHashMap;
import template.Accumulator;
import template.ThreadWaiter;

/** One frozen live-gain pass: construct and measure in parallel, then publish a
 * deterministic overlap-adjusted selection. Returned state can feed the next
 * phase directly. Single-phase passes use measured gains; combined passes can
 * use measured or predicted ranking with exact diagnostic coverage accounting.
 * @author Brian Bushnell, Keqing
 */
public final class LiveProbeClustering{

	/** Compatibility entrypoint for a zero additive margin, using fraction one.
	 * Nonzero additive margins are rejected because they break split/merge symmetry.
	 * A strict inequality
	 * is used. All original input ordinals and evaluated model bytes are preserved
	 * by the state publisher; no mutable candidate enters the next pass early.
	 */
	public static Result run(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final LayerCakeSchedule.Phase phase,final double margin){
		return run(state,sequences,builder,memberGate,mergeGate,filter,minimumLength,objective,options,phase,margin,Integer.MAX_VALUE);
	}

	/** Optional inclusive bound on BOTH current fusion endpoint sizes. Neighbor search
	 * still uses the full active roster and unchanged ratio/gates; only completed
	 * mutual pairs are filtered. Zero excludes all fusion pairs. Fission is unchanged.
	 */
	public static Result run(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final LayerCakeSchedule.Phase phase,final double margin,
			final int maximumFusionMembers){
		if(margin!=0){throw new IllegalArgumentException("Additive gain margins are unsupported; use a symmetric slot fraction");}
		return runWithFraction(state,sequences,builder,memberGate,mergeGate,filter,minimumLength,objective,
			options,phase,1,maximumFusionMembers);
	}

	/** Reads B from this phase's input population and applies only the signed fraction*B threshold. */
	public static Result runWithFraction(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final LayerCakeSchedule.Phase phase,final double slotFraction,
			final int maximumFusionMembers){
		if(maximumFusionMembers<0){throw new IllegalArgumentException("Maximum fusion endpoint membership must be nonnegative");}
		if(state==null || sequences==null || builder==null || memberGate==null || mergeGate==null || filter==null ||
			objective==null || options==null || phase==null || !Double.isFinite(slotFraction) || slotFraction<0 ||
			minimumLength<memberGate.minSpan || sequences.length!=state.sequenceCount() || options.localRefineRounds<1){
			throw new IllegalArgumentException("Live passes need a matching frozen state, bounded refinement, native gates and finite nonnegative slot fraction");
		}
		final ClusterSketchCodec.Header h=state.sourceHeader();
		if(!builder.matchesGeometry(h.width,h.geometry,h.sketchSeed,h.positionMode,h.partitions,h.overlap,h.blacklistSha80)){
			throw new IllegalArgumentException("Live model and query sketch geometry/blacklist must remain bound");
		}
		final long start=System.nanoTime();
		final Family[] active=state.active().clone();
		Arrays.sort(active,(a,b)->Integer.compare(a.familyId,b.familyId));
		final double slot=LayerCakeSchedule.summarize(active,options.familyBudget,sequences.length).slotCost;
		final double threshold=(phase==LayerCakeSchedule.Phase.FISSION ? slot : -slot)*slotFraction;
		if(!Double.isFinite(threshold)){throw new IllegalArgumentException("Live signed slot margin overflows");}
		checkRankingPopulation(state,sequences);
		final Requests requests=generate(active,builder,mergeGate,options,phase,state.firstNewId(),maximumFusionMembers,
			null,state.memberRankingInputs());
		return runRequests(state,sequences,builder,memberGate,mergeGate,filter,minimumLength,objective,
			options,phase,requests,slot,threshold,null,start);
	}

	/** Measures one frozen mixed batch, saves every raw row, fits once and commits once.
	 * ALL selected models receive fission rows, including explicitly unavailable
	 * empty/small models. The caller must set fissionMinSize=1; minimum child size still
	 * governs construction. Fusion neighbors are selected from the full active set.
	 * Initial B is shared: split threshold=+B, fusion threshold=-B. Exact measured
	 * marginal gain alone determines acceptance. Top-budget coverage is diagnostic.
	 */
	public static CombinedResult runCombined(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final int maximumFusionMembers,
			final String rawReport,final boolean overwrite){
		return runCombined(state,sequences,builder,memberGate,mergeGate,filter,minimumLength,objective,
			options,maximumFusionMembers,rawReport,overwrite,false,0,0,true,null);
	}

	/** Selects measured or predicted ranking with explicit signed phase thresholds. A diagnostic
	 * fit cannot exclude proposals in measured mode, including those lacking features.
	 * Numerical fit insufficiency is unavailable only in measured mode; predicted mode
	 * still requires a fit. Invalid data and unexpected errors always fail loudly.
	 * Signed measured gain is compared only with the supplied phase threshold.
	 */
	public static CombinedResult runCombined(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final int maximumFusionMembers,
			final String rawReport,final boolean overwrite,final boolean rankByPrediction,
			final double fissionThreshold,final double fusionThreshold){
		return runCombined(state,sequences,builder,memberGate,mergeGate,filter,minimumLength,objective,
			options,maximumFusionMembers,rawReport,overwrite,rankByPrediction,fissionThreshold,fusionThreshold,false,null);
	}

	/** Optional wide sidecar does not alter selection. Predictors are captured before
	 * pool scoring, and every requested row is written before fitting can fail.
	 */
	public static CombinedResult runCombined(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final int maximumFusionMembers,
			final String rawReport,final boolean overwrite,final boolean rankByPrediction,
			final double fissionThreshold,final double fusionThreshold,final WideProbeFeatures.Output wide){
		return runCombined(state,sequences,builder,memberGate,mergeGate,filter,minimumLength,objective,
			options,maximumFusionMembers,rawReport,overwrite,rankByPrediction,fissionThreshold,fusionThreshold,false,wide);
	}

	/** Shared measurement/publication path; the default overload derives symmetric unit-fraction thresholds. */
	private static CombinedResult runCombined(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final int maximumFusionMembers,
			final String rawReport,final boolean overwrite,final boolean rankByPrediction,
			final double fissionThreshold,final double fusionThreshold,final boolean deriveSlotThresholds,final WideProbeFeatures.Output wide){
		if(state==null || sequences==null || builder==null || memberGate==null || mergeGate==null || filter==null ||
			objective==null || options==null || options.fissionMinSize!=1 || options.localRefineRounds<1 ||
			minimumLength<memberGate.minSpan || sequences.length!=state.sequenceCount() || maximumFusionMembers<0 || rawReport==null ||
			!Double.isFinite(fissionThreshold) || !Double.isFinite(fusionThreshold)){
			throw new IllegalArgumentException("Combined probes need complete frozen inputs, fissionMinSize=1, refinement and a raw report path");
		}
		final ClusterSketchCodec.Header header=state.sourceHeader();
		if(!builder.matchesGeometry(header.width,header.geometry,header.sketchSeed,header.positionMode,header.partitions,header.overlap,header.blacklistSha80)){
			throw new IllegalArgumentException("Combined model and query sketch geometry must remain bound");
		}
		final long start=System.nanoTime();
		final Family[] active=state.active().clone(); Arrays.sort(active,(a,b)->Integer.compare(a.familyId,b.familyId));
		final Family[] requested=state.models().clone(); Arrays.sort(requested,(a,b)->Integer.compare(a.familyId,b.familyId));
		final IntHashMap requestedIndex=new IntHashMap(mapCapacity(requested.length));
		for(int i=0;i<requested.length;i++){requestedIndex.put(requested[i].familyId,i);}
		final LayerCakeSchedule.Summary before=LayerCakeSchedule.summarize(active,options.familyBudget,sequences.length);
		final double splitLimit=deriveSlotThresholds ? before.slotCost : fissionThreshold;
		final double mergeLimit=deriveSlotThresholds ? -before.slotCost : fusionThreshold;
		checkRankingPopulation(state,sequences);
		final Requests fusions=generate(active,builder,mergeGate,options,LayerCakeSchedule.Phase.FUSION,state.firstNewId(),maximumFusionMembers,
			null,state.memberRankingInputs());
		final int total=Math.addExact(requested.length,fusions.size());
		final long next=(long)state.firstNewId()+active.length;
		if(next>=Integer.MAX_VALUE-8){throw new IllegalStateException("Combined split reservations exhaust stable family IDs");}
		final int[] left=new int[total],right=new int[total],children=new int[total];Arrays.fill(children,-1);
		for(int i=0,child=state.firstNewId();i<requested.length;i++){
			left[i]=i;right[i]=-1;if(requested[i].size()>0){children[i]=child++;}
		}
		for(int i=0;i<fusions.size();i++){
			left[requested.length+i]=requestedIndex.get(active[fusions.left[i]].familyId);
			right[requested.length+i]=requestedIndex.get(active[fusions.right[i]].familyId);
		}
		final Requests requests=new Requests(left,right,state.firstNewId(),(int)next,children);
		final Outcome[] outcomes=new Outcome[total];
		if(active.length>0){
			final TailPoolIndex pool=LocalProbePool.build(state.residualOrdinals(),sequences,builder,options.threads);
			final LocalProbeEvaluator evaluator=new LocalProbeEvaluator(sequences,pool,builder,memberGate,filter,minimumLength);
			final LayerCakeClustering engine=LayerCakeClustering.liveProbesWithFeatures(active,state.models(),sequences,builder,
				memberGate,mergeGate,objective,options,state.firstNewId());
			bindRanking(engine,state);
			new Work(requested,requests,null,engine,evaluator,sequences.length,outcomes,options.threads,true,wide,before.slotCost).run();
		}else{
			for(int i=0;i<requested.length;i++){
				outcomes[i]=new Outcome(null,new Row(i,requested[i].familyId,-1,-1,"zero_member_parent",false,null,LayerCakeSchedule.Phase.FISSION));
				if(wide!=null){outcomes[i].row.wide=new WideProbeFeatures(requested[i],null,null,sequences.length,requested.length,before.slotCost,wide.layer,true);}
			}
		}
		final Row[] rows=new Row[total];
		int fitted=0,measured=0;
		for(int i=0;i<total;i++){
			rows[i]=outcomes[i].row;
			rows[i].threshold=rows[i].phase==LayerCakeSchedule.Phase.FISSION ? splitLimit : mergeLimit;
			if(rows[i].available){measured++;}
			if(rows[i].available && rows[i].features!=null){fitted++;}
		}
		// This write completes before any fit or selection. Even a required fit's
		// rank failure preserves measurements without publishing partial state.
		writeCombinedRows(rows,rawReport,overwrite,before.slotCost,false);
		if(wide!=null){writeWideRows(rows,wide,overwrite);}
		final double[] cohesion=new double[fitted],separation=new double[fitted],pressure=new double[fitted];
		final long[] gains=new long[fitted];
		for(int i=0,at=0;i<total;i++){
			final Row row=rows[i];
			if(!row.available || row.features==null){continue;}
			cohesion[at]=row.features.cohesionDelta;separation[at]=row.features.separationDelta;pressure[at]=row.features.sizePressure;
			gains[at++]=row.gain;
		}
		DdlClustProbeRegression.LiveFit fit=null;
		String fitUnavailableReason=null;
		try{fit=DdlClustProbeRegression.fitLive(cohesion,separation,pressure,gains);}
		catch(DdlClustProbeRegression.FitUnavailableException unavailable){
			if(rankByPrediction){throw unavailable;}
			fitUnavailableReason=unavailable.reason();
		}
		final int selectable=rankByPrediction ? fitted : measured;
		final LiveProbeState.Proposal[] proposals=new LiveProbeState.Proposal[selectable];
		final LiveProbeSelection.Operation[] operations=new LiveProbeSelection.Operation[selectable];
		final double[] predictions=rankByPrediction ? new double[selectable] : null,thresholds=new double[selectable];
		for(int i=0,at=0;i<total;i++){
			final Row row=rows[i];
			if(!row.available){continue;}
			if(fit!=null && row.features!=null){
				row.prediction=fit.predict(row.features.cohesionDelta,row.features.separationDelta,row.features.sizePressure);
				row.predictionAvailable=true;
			}
			if(rankByPrediction && !row.predictionAvailable){continue;}
			proposals[at]=outcomes[i].proposal;operations[at]=proposals[at].selectionOperation();thresholds[at]=row.threshold;
			if(rankByPrediction){predictions[at]=row.prediction;}
			at++;
		}
		final BudgetGate budget=new BudgetGate(state,proposals,options.familyBudget,rows);
		if(budget.currentTop!=before.topMembers){throw new IllegalStateException("Initial budget histogram disagrees with the independent frozen census");}
		final LiveProbeSelection.Result selection=rankByPrediction ? LiveProbeSelection.selectPredicted(operations,predictions,thresholds,budget) :
			LiveProbeSelection.selectMeasured(operations,thresholds,budget);
		final LiveProbeState committed=state.commit(proposals,selection,requests.nextUnusedId);
		long gain=0;
		for(int i=0;i<selection.size();i++){
			final Row row=rows[selection.requestId(i)];
			if(row.accepted){throw new IllegalStateException("Combined selection duplicated a request");}
			row.accepted=true;row.commitOrder=i;row.marginal=selection.marginalGain(i);
			if(row.predictionAvailable){
				// Exact unique rescue claims discount both diagnostics and measured gain.
				row.predictedMarginal=row.prediction-Math.subtractExact(row.gain,row.marginal);
			}
			gain=Math.addExact(gain,row.marginal);
		}
		final long topAfter=LayerCakeSchedule.summarize(committed.active(),options.familyBudget,sequences.length).topMembers;
		if(gain!=(long)committed.assigned()-state.assigned() || gain!=committed.lastMarginalGain() ||
			topAfter!=budget.currentTop){
			throw new IllegalStateException("Published coverage must match exact signed measured gain and the diagnostic top-budget ledger");
		}
		return new CombinedResult(committed,rows,fit,fitUnavailableReason,fitted,rankByPrediction,splitLimit,mergeLimit,
			before.slotCost,before.topMembers,topAfter,gain,System.nanoTime()-start);
	}

	/** Exhaustive eligible parents or reciprocal size-compatible pairs; no response-dependent sampling. */
	static Requests generateFissions(final LiveProbeState state,final int[] parentFamilyIds,
			final LayerCakeClustering.Options options){
		if(state==null || parentFamilyIds==null || options==null){
			throw new IllegalArgumentException("Explicit fission waves need a state, parent IDs and options");
		}
		final Family[] active=state.active().clone(); Arrays.sort(active,(a,b)->Integer.compare(a.familyId,b.familyId));
		final IntHashMap index=new IntHashMap(mapCapacity(active.length));
		for(int i=0;i<active.length;i++){index.put(active[i].familyId,i);}
		final long floor=Math.max((long)options.fissionMinSize,2L*options.minChildSize);
		final IntList left=new IntList();
		for(int id:parentFamilyIds){
			final int at=index.get(id);
			if(at>=0 && active[at].size()>=floor){left.add(at);}
		}
		final int[] right=new int[left.size];Arrays.fill(right,-1);
		final int next=Math.addExact(state.firstNewId(),left.size);
		if(next>=Integer.MAX_VALUE-8){throw new IllegalStateException("Live split reservations exhaust stable family IDs");}
		return new Requests(left.toArray(),right,state.firstNewId(),next);
	}

	static Requests generateFusions(final LiveProbeState state,final FamilySketchBuilder builder,
			final TailAlignmentGate mergeGate,final LayerCakeClustering.Options options){
		return generateFusions(state,builder,mergeGate,options,null);
	}

	static Requests generateFusions(final LiveProbeState state,final FamilySketchBuilder builder,
			final TailAlignmentGate mergeGate,final LayerCakeClustering.Options options,final int[] eligibleFamilyIds){
		if(state==null || builder==null || mergeGate==null || options==null){
			throw new IllegalArgumentException("Fusion generation needs a state, sketches, gates and options");
		}
		final Family[] active=state.active().clone(); Arrays.sort(active,(a,b)->Integer.compare(a.familyId,b.familyId));
		return generate(active,builder,mergeGate,options,LayerCakeSchedule.Phase.FUSION,state.firstNewId(),Integer.MAX_VALUE,
			eligibleFamilyIds,state.memberRankingInputs());
	}

	static Result runRequests(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final LayerCakeSchedule.Phase phase,final Requests requests,
			final double slot,final double threshold,final WideProbeFeatures.Output wide,final long start){
		return runRequests(state,sequences,builder,memberGate,mergeGate,filter,minimumLength,objective,options,
			phase,requests,slot,threshold,0,wide,start,null);
	}

	/** Size pressure affects selection only; measurements and published coverage remain raw counts. */
	static Result runRequests(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final LayerCakeSchedule.Phase phase,final Requests requests,
			final double slot,final double threshold,final double sizeLambda,final WideProbeFeatures.Output wide,final long start){
		return runRequests(state,sequences,builder,memberGate,mergeGate,filter,minimumLength,objective,options,
			phase,requests,slot,threshold,sizeLambda,wide,start,null);
	}

	static Result runRequests(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final LayerCakeSchedule.Phase phase,final Requests requests,
			final double slot,final double threshold,final double sizeLambda,final WideProbeFeatures.Output wide,
			final long start,final LocalProbePool.ResidualCache residualCache){
		if(!Double.isFinite(sizeLambda) || sizeLambda<0){throw new IllegalArgumentException("sizelambda must be finite and nonnegative");}
		if(requests==null || phase==null || !Double.isFinite(slot) || slot<0 || !Double.isFinite(threshold)){
			throw new IllegalArgumentException("Explicit live passes need requests, phase and a finite signed threshold");
		}
		if(state==null || sequences==null || builder==null || memberGate==null || mergeGate==null || filter==null ||
			objective==null || options==null || minimumLength<memberGate.minSpan || sequences.length!=state.sequenceCount() ||
			options.localRefineRounds<1){
			throw new IllegalArgumentException("Live passes need a matching frozen state, bounded refinement and native gates");
		}
		final ClusterSketchCodec.Header h=state.sourceHeader();
		if(!builder.matchesGeometry(h.width,h.geometry,h.sketchSeed,h.positionMode,h.partitions,h.overlap,h.blacklistSha80)){
			throw new IllegalArgumentException("Live model and query sketch geometry/blacklist must remain bound");
		}
		if(residualCache!=null && phase!=LayerCakeSchedule.Phase.FISSION){
			throw new IllegalArgumentException("Residual caches are scoped to recursive fission waves");
		}
		if(residualCache!=null && !residualCache.matchesContext(sequences,builder,options.threads)){
			throw new IllegalArgumentException("Residual cache belongs to a different frozen input store or sketch recipe");
		}
		checkRankingPopulation(state,sequences);
		final Family[] active=state.active().clone();
		Arrays.sort(active,(a,b)->Integer.compare(a.familyId,b.familyId));
		final Outcome[] outcomes=new Outcome[requests.size()];
		if(requests.size()>0){
			final TailPoolIndex pool;
			try(PhaseCostLog.Span interval=PhaseCostLog.begin("residual_index")){
				pool=residualCache==null ? LocalProbePool.build(state.residualOrdinals(),sequences,builder,options.threads) :
					residualCache.build(state);
			}
			final LocalProbeEvaluator evaluator=new LocalProbeEvaluator(sequences,pool,builder,memberGate,filter,minimumLength);
			final Work work;
			try(PhaseCostLog.Span interval=PhaseCostLog.begin("proposal_setup")){
				final LayerCakeClustering engine=wide==null ?
					LayerCakeClustering.liveProbes(active,state.models(),sequences,builder,memberGate,mergeGate,objective,options,state.firstNewId()) :
					LayerCakeClustering.liveProbesWithFeatures(active,state.models(),sequences,builder,memberGate,mergeGate,objective,options,state.firstNewId());
				bindRanking(engine,state);
				work=new Work(active,requests,phase,engine,evaluator,sequences.length,outcomes,options.threads,wide!=null,wide,slot);
			}
			try(PhaseCostLog.Span interval=PhaseCostLog.begin("proposal_evaluation")){
				work.run();
			}
		}
		try(PhaseCostLog.Span interval=PhaseCostLog.begin("selection_publication")){
			final LiveProbeState.Proposal[] proposals;
			final LiveProbeSelection.Operation[] operations;
			final double[] offsets;
			final Row[] rows;
			int count=0;
			try(PhaseCostLog.Span inputs=PhaseCostLog.begin("selection_inputs")){
				for(Outcome outcome:outcomes){if(outcome.proposal!=null){count++;}}
				proposals=new LiveProbeState.Proposal[count];
				operations=new LiveProbeSelection.Operation[count];
				offsets=sizeLambda==0 ? null : new double[count];
				for(int i=0,at=0;i<outcomes.length;i++){
					if(outcomes[i].proposal!=null){
						proposals[at]=outcomes[i].proposal;operations[at]=proposals[at].selectionOperation();
						if(offsets!=null){
							final Row row=outcomes[i].row;
							row.sizePressure=sizePressure(row.affected,slot,sizeLambda)*(phase==LayerCakeSchedule.Phase.FISSION ? 1 : -1);
							offsets[at]=row.sizePressure;
						}
						at++;
					}
				}
				rows=new Row[outcomes.length];
				for(int i=0;i<rows.length;i++){rows[i]=outcomes[i].row;}
			}
			final LiveProbeSelection.Result selection;
			try(PhaseCostLog.Span select=PhaseCostLog.begin("selection")){
				final BudgetGate budget=new BudgetGate(state,proposals,options.familyBudget,rows);
				selection=offsets==null ? LiveProbeSelection.select(operations,threshold,budget) :
					LiveProbeSelection.selectAdjusted(operations,offsets,threshold,budget);
			}
			final LiveProbeState committed;
			try(PhaseCostLog.Span publish=PhaseCostLog.begin("publication")){
				committed=state.commit(proposals,selection,requests.nextUnusedId);
			}
			long sum=0;
			try(PhaseCostLog.Span annotate=PhaseCostLog.begin("accepted_rows")){
				for(int i=0;i<selection.size();i++){
					final int id=selection.requestId(i);
					if(id<0 || id>=rows.length || rows[id].accepted){throw new IllegalStateException("Live selection must name unique measured request rows");}
					rows[id].accepted=true;rows[id].marginal=selection.marginalGain(i);rows[id].commitOrder=i;
					sum=Math.addExact(sum,selection.marginalGain(i));
				}
			}
			if(sum!=(long)committed.assigned()-state.assigned() || sum!=committed.lastMarginalGain()){
				throw new IllegalStateException("Published covered-tail delta must equal the actual overlap-adjusted accepted gains");
			}
			return new Result(committed,rows,count,selection.size(),slot,threshold,sizeLambda,phase,sum,System.nanoTime()-start);
		}
	}

	/** Quadratic pressure in tail-count units; zero B explicitly disables pressure.
	 * A split parent and its reverse union use the same conserved original member count.
	 */
	public static double sizePressure(final int affected,final double slot,final double lambda){
		if(affected<0 || !Double.isFinite(slot) || slot<0 || !Double.isFinite(lambda) || lambda<0){
			throw new IllegalArgumentException("Size pressure requires nonnegative original membership, finite B and lambda");
		}
		if(slot==0 || lambda==0){return 0;}
		final double ratio=affected/slot,pressure=lambda*slot*ratio*ratio;
		if(!Double.isFinite(pressure)){throw new IllegalArgumentException("Quadratic size pressure overflowed");}
		return pressure;
	}

	/** Cache identity is checked even for passes with no proposals, before publication. */
	private static void checkRankingPopulation(final LiveProbeState state,final ProteinSequence[] sequences){
		final MemberRankInputs inputs=state.memberRankingInputs();
		if(inputs!=null && !inputs.belongsTo(sequences)){
			throw new IllegalArgumentException("A ranked live pass must retain its identical original protein population");
		}
	}

	/** Both combined and ordinary workers inherit the state envelope before forking. */
	private static void bindRanking(final LayerCakeClustering engine,final LiveProbeState state){
		if(state.memberRankingInputs()!=null){engine.withMemberRanking(state.memberRankingInputs(),state.memberRankingLayer());}
	}

	private static Requests generate(final Family[] active,final FamilySketchBuilder builder,
			final TailAlignmentGate mergeGate,final LayerCakeClustering.Options options,
			final LayerCakeSchedule.Phase phase,final int firstNewId,final int maximumFusionMembers,
			final int[] eligibleFusionFamilies,final MemberRankInputs inputs){
		assert(active!=null && firstNewId>=0) : "Live requests require a frozen ordered population and an unused ID boundary";
		final IntList left=new IntList(),right=new IntList();
		if(phase==LayerCakeSchedule.Phase.FISSION){
			final long floor=Math.max((long)options.fissionMinSize,2L*options.minChildSize);
			for(int i=0;i<active.length;i++){if(active[i].size()>=floor){left.add(i);right.add(-1);}}
		}else if(active.length>1){
			final int[] counts=new int[active.length];
			for(int i=0;i<counts.length;i++){counts[i]=active[i].size();}
			final IntHashMap eligible=eligibleFusionFamilies==null ? null :
				new IntHashMap(Math.max(4,eligibleFusionFamilies.length*2));
			if(eligible!=null){for(int id:eligibleFusionFamilies){eligible.put(id,1);}}
			final MutualNearestMerge.Result pairs=MutualNearestMerge.selectWithinSizeRatio(active,counts,builder,
				mergeGate,options.fusionSizeRatio,options.mergeShortlist,options.threads,options.indexThreads,inputs);
			for(MutualNearestMerge.Pair pair:pairs.pairs){
				if(counts[pair.leftIndex]<=maximumFusionMembers && counts[pair.rightIndex]<=maximumFusionMembers){
					final int a=active[pair.leftIndex].familyId,b=active[pair.rightIndex].familyId;
					if(eligible==null || eligible.get(a)>=0 && eligible.get(b)>=0){
						left.add(pair.leftIndex);right.add(pair.rightIndex);
					}
				}
			}
		}
		final int next=Math.addExact(firstNewId,phase==LayerCakeSchedule.Phase.FISSION ? left.size : 0);
		if(next>=Integer.MAX_VALUE-8){throw new IllegalStateException("Live split reservations exhaust stable family IDs");}
		return new Requests(left.toArray(),right.toArray(),firstNewId,next);
	}

	static final class Requests{
		Requests(int[] left_,int[] right_,int first_,int next_){
			this(left_,right_,first_,next_,null);
		}
		Requests(int[] left_,int[] right_,int first_,int next_,int[] children_){
			assert(left_.length==right_.length) : "Each live request has exactly one parent or pair";
			left=left_;right=right_;firstId=first_;nextUnusedId=next_;
			children=children_;
		}
		int size(){return left.length;}
		final int[] left,right,children;
		final int firstId,nextUnusedId;
	}

	private static final class Work implements Accumulator<Worker>{
		Work(Family[] active_,Requests requests_,LayerCakeSchedule.Phase phase_,LayerCakeClustering engine,
				LocalProbeEvaluator evaluator,int sequenceCount_,Outcome[] output_,int threads){
			this(active_,requests_,phase_,engine,evaluator,sequenceCount_,output_,threads,false);
		}
		Work(Family[] active_,Requests requests_,LayerCakeSchedule.Phase phase_,LayerCakeClustering engine,
				LocalProbeEvaluator evaluator,int sequenceCount_,Outcome[] output_,int threads,boolean features_){
			this(active_,requests_,phase_,engine,evaluator,sequenceCount_,output_,threads,features_,null,0);
		}
		Work(Family[] active_,Requests requests_,LayerCakeSchedule.Phase phase_,LayerCakeClustering engine,
				LocalProbeEvaluator evaluator,int sequenceCount_,Outcome[] output_,int threads,boolean features_,
				WideProbeFeatures.Output wide_,double slot_){
			active=active_;requests=requests_;phase=phase_;sequenceCount=sequenceCount_;output=output_;
			features=features_;wide=wide_;slot=slot_;
			for(int i=0;i<Math.min(threads,requests.size());i++){
				workers.add(new Worker(this,engine.forkLocalProbes(phase!=LayerCakeSchedule.Phase.FISSION),evaluator.newWorker()));
			}
		}
		void run(){
			try{if(!ThreadWaiter.startAndWait(workers,this)){fail(new IllegalStateException("Live workers did not finish"));}}
			catch(Throwable problem){
				fail(problem);final ArrayList<Worker> started=new ArrayList<Worker>();
				for(Worker worker:workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
				ThreadWaiter.waitForThreadsToFinish(started);
			}
			if(failure!=null || completed!=output.length){throw new IllegalStateException("Live pass failed before any state publication",failure);}
			for(Outcome outcome:output){if(outcome==null){throw new IllegalStateException("A live request was omitted");}}
		}
		int claim(){return (int)Math.min(next.getAndIncrement(),output.length);}
		synchronized void fail(Throwable problem){if(failure==null){failure=problem;}}
		@Override public void accumulate(Worker worker){
			if(!worker.success){fail(new IllegalStateException("Live worker exited without completion"));}
			completed=Math.addExact(completed,worker.completed);
		}
		@Override public boolean success(){return failure==null && completed==output.length;}
		@Override public ReadWriteLock rwlock(){return null;}
		final Family[] active;
		final Requests requests;
		final LayerCakeSchedule.Phase phase;
		final int sequenceCount;
		final boolean features;
		final WideProbeFeatures.Output wide;
		final double slot;
		final Outcome[] output;
		final AtomicLong next=new AtomicLong();
		final ArrayList<Worker> workers=new ArrayList<Worker>();
		volatile Throwable failure;
		int completed;
	}

	private static final class Worker extends Thread{
		Worker(Work owner_,LayerCakeClustering engine_,LocalProbeEvaluator.Worker evaluator_){owner=owner_;engine=engine_;evaluator=evaluator_;}
		@Override public void run(){
			try{
				for(int id=owner.claim();id<owner.output.length && owner.failure==null;id=owner.claim()){
					final int ai=owner.requests.left[id],bi=owner.requests.right[id];
					final Family a=owner.active[ai],b=bi<0 ? null : owner.active[bi];
					final LayerCakeSchedule.Phase phase=owner.phase==null ?
						(b==null ? LayerCakeSchedule.Phase.FISSION : LayerCakeSchedule.Phase.FUSION) : owner.phase;
					final int child=owner.requests.children!=null ? owner.requests.children[id] :
						phase==LayerCakeSchedule.Phase.FISSION ? owner.requests.firstId+id : -1;
					if(a.size()==0){
						assert(owner.features && phase==LayerCakeSchedule.Phase.FISSION && child==-1) : "Only combined zero-model reporting bypasses construction";
						final Row empty=new Row(id,a.familyId,-1,-1,"zero_member_parent",false,null,phase);
						if(owner.wide!=null){empty.wide=new WideProbeFeatures(a,null,null,owner.sequenceCount,owner.active.length,owner.slot,owner.wide.layer,true);}
						owner.output[id]=new Outcome(null,empty);completed++;continue;
					}
					final LayerCakeClustering.Probe probe=engine.constructByFamilyId(phase,a.familyId,b==null ? -1 : b.familyId,child);
					final Row row=new Row(id,a.familyId,b==null ? -1 : b.familyId,child,probe.reason,
						probe.audit!=null && probe.audit.available(),probe.refinement,phase);
					// Capture before evaluator access: neither pool outcomes nor later
					// acceptance/overlap decisions can leak into these predictors.
					if(owner.wide!=null){
						final long featureStart=System.nanoTime();
						final boolean split=phase==LayerCakeSchedule.Phase.FISSION;
						row.wide=new WideProbeFeatures(split ? a : probe.candidateFirst,
							split ? probe.candidateFirst : a,split ? probe.candidateSecond : b,
							owner.sequenceCount,owner.active.length,owner.slot,owner.wide.layer,split);
						row.features=engine.measureLocalFeatures(probe,row.wide);
						row.wide.deltas(row.features);
						row.featureNanos=System.nanoTime()-featureStart;
					}
					LiveProbeState.Proposal proposal=null;
					if(probe.candidateFirst!=null){
						final int[] affected=new int[Math.addExact(a.size(),b==null ? 0 : b.size())];
						for(int i=0;i<a.size();i++){affected[i]=a.memberIndex(i);}
						if(b!=null){for(int i=0;i<b.size();i++){affected[a.size()+i]=b.memberIndex(i);}}
						final LiveProbeAssignments collector=new LiveProbeAssignments(owner.sequenceCount,probe.candidateFirst.familyId,
							probe.candidateSecond==null ? -1 : probe.candidateSecond.familyId);
						final LocalProbeEvaluator.Result scored=evaluator.evaluate(probe.candidateFirst,probe.candidateSecond,affected,collector);
						final LiveProbeAssignments.Snapshot assignments=collector.snapshot();
						if(scored.affected!=affected.length || assignments.rescuedCount()!=scored.poolAccepted ||
							assignments.retainedCount()!=scored.affectedAccepted){
							throw new IllegalStateException("Retained winner evidence must match the exact measured local counts");
						}
						proposal=new LiveProbeState.Proposal(id,b==null ? new int[]{a.familyId} : new int[]{a.familyId,b.familyId},
							probe.candidateFirst,probe.candidateSecond,assignments,scored.netGain);
						row.available=true;row.gain=scored.netGain;row.rescued=scored.poolAccepted;
						row.retained=scored.affectedAccepted;row.affected=affected.length;
						if(owner.features && owner.wide==null){row.features=engine.measureLocalFeatures(probe);}
					}
					owner.output[id]=new Outcome(proposal,row);completed++;
				}
				success=owner.failure==null;
			}catch(Throwable problem){owner.fail(problem);}
		}
		final Work owner;
		final LayerCakeClustering engine;
		final LocalProbeEvaluator.Worker evaluator;
		int completed;
		boolean success;
	}

	private static final class Outcome{
		Outcome(LiveProbeState.Proposal p,Row r){proposal=p;row=r;}
		final LiveProbeState.Proposal proposal;
		final Row row;
	}
	private static final class Row{
		Row(int id_,int a_,int b_,int child_,String reason_,boolean veto_,CompetitiveRefinement.Summary refinement_,LayerCakeSchedule.Phase phase_){
			id=id_;a=a_;b=b_;child=child_;reason=reason_;veto=veto_;refinement=refinement_;
			phase=phase_;
		}
		final int id,a,b,child;
		final String reason;
		final boolean veto;
		final CompetitiveRefinement.Summary refinement;
		final LayerCakeSchedule.Phase phase;
		LayerCakeClustering.ProbeFeatures features;
		WideProbeFeatures wide;
		long featureNanos;
		boolean available,accepted;
		boolean predictionAvailable,topBudgetDecreased,emptyChildVeto;
		double threshold,prediction,predictedMarginal,sizePressure;
		int commitOrder=-1,affected;
		long gain,marginal,rescued,retained;
	}

	/** Exact final membership gate, separate from model payload and prediction.
	 * The callback preserves top-budget diagnostics, rejects fissions that lose a
	 * final child after rescue overlap, and throws on invalid conservation evidence.
	 */
	public static final class BudgetGate implements LiveProbeSelection.AcceptanceGate{
		public BudgetGate(final LiveProbeState state,final LiveProbeState.Proposal[] proposals,final int budget){
			this(state,proposals,budget,null);
		}
		private BudgetGate(final LiveProbeState state,final LiveProbeState.Proposal[] proposals,final int budget,final Row[] rows_){
			if(state==null || proposals==null || budget<1){throw new IllegalArgumentException("A budget gate needs a frozen state, proposals and positive budget");}
			source=state;models=state.models();input=proposals.clone();rows=rows_;limit=budget;used=new boolean[models.length];
			proposalRows=new int[input.length][];affectedSizes=new long[input.length];
			int maximum=0;
			for(int i=0;i<models.length;i++){maximum=Math.max(maximum,models[i].size());}
			for(int i=0;i<input.length;i++){
				final LiveProbeState.Proposal proposal=input[i];
				if(proposal==null){throw new IllegalArgumentException("A budget proposal cannot be missing");}
				final int[] parents=proposal.parentIds(),parentRows=new int[parents.length];
				long affected=0;
				for(int p=0;p<parents.length;p++){
					final int at=row(parents[p]);
					if(at<0 || models[at].size()<1){throw new IllegalArgumentException("Budget parents must belong to the frozen active roster");}
					parentRows[p]=at;affected+=models[at].size();
				}
				proposalRows[i]=parentRows;affectedSizes[i]=affected;
				winnerCounts(proposal,null,counts);
				maximum=Math.max(maximum,Math.max(counts[0],counts[1]));
			}
			sizes=new SizeLedger(maximum);
			for(Family model:models){sizes.add(model.size(),1);}
			currentTop=sizes.top(limit);
		}
		@Override
		public boolean tryAccept(final int operationIndex,final long measuredMarginalGain){
			final LiveProbeState.Proposal proposal=input[operationIndex];
			final int[] parentRows=proposalRows[operationIndex];
			final long affected=affectedSizes[operationIndex];
			for(int at:parentRows){
				if(used[at]){throw new IllegalStateException("Budget gate received an already consumed parent");}
			}
			winnerCounts(proposal,claimed,counts);
			if((long)counts[0]+counts[1]-affected!=measuredMarginalGain){
				throw new IllegalStateException("Budget child winners must match the selector's exact marginal coverage");
			}
			if(proposal.second!=null && (counts[0]==0 || counts[1]==0)){
				if(rows!=null){rows[proposal.requestId].emptyChildVeto=true;}
				emptyChildRejected++;return false;
			}
			for(int at:parentRows){sizes.add(models[at].size(),-1);}
			sizes.add(counts[0],1);if(proposal.second!=null){sizes.add(counts[1],1);}
			final long after=sizes.top(limit);
			if(after<currentTop){
				if(rows!=null){rows[proposal.requestId].topBudgetDecreased=true;}
				decreased++;
			}
			for(int at:parentRows){used[at]=true;}
			for(int i=0;i<proposal.assignments.rescuedCount();i++){claimed.set(proposal.assignments.rescuedOrdinal(i));}
			currentTop=after;return true;
		}
		public long topMembers(){return currentTop;}
		public int rejected(){return emptyChildRejected;}
		public int decreased(){return decreased;}
		private int row(final int parent){return source.rowOf(parent);}
		/** Counts actual winning retained rows plus unclaimed rescue rows by child, never construction sizes. */
		private static void winnerCounts(final LiveProbeState.Proposal proposal,final BitSet claimed,final int[] counts){
			Arrays.fill(counts,0);
			final LiveProbeAssignments.Snapshot rows=proposal.assignments;
			for(int i=0;i<rows.retainedCount();i++){addWinner(proposal,rows.retainedWinner(i),counts);}
			for(int i=0;i<rows.rescuedCount();i++){
				if(claimed==null || !claimed.get(rows.rescuedOrdinal(i))){addWinner(proposal,rows.rescuedWinner(i),counts);}
			}
		}
		private static void addWinner(final LiveProbeState.Proposal proposal,final int winner,final int[] counts){
			final int side=winner==proposal.first.familyId ? 0 : proposal.second!=null && winner==proposal.second.familyId ? 1 : -1;
			if(side<0){throw new IllegalArgumentException("Budget evidence names an unknown child winner");}
			counts[side]=Math.addExact(counts[side],1);
		}
		private final Family[] models;
		private final LiveProbeState source;
		private final LiveProbeState.Proposal[] input;
		private final int[][] proposalRows;
		private final long[] affectedSizes;
		private final BitSet claimed=new BitSet();
		private final boolean[] used;
		private final int[] counts=new int[2];
		private final SizeLedger sizes;
		private final int limit;
		private final Row[] rows;
		private long currentTop;
		private int decreased,emptyChildRejected;
	}

	private static int mapCapacity(final int entries){
		return Math.max(4, entries>=(Integer.MAX_VALUE-8)/2 ? Integer.MAX_VALUE-8 : entries*2+1);
	}

	/** Primitive Fenwick histogram of positive family sizes and their weighted sum. */
	private static final class SizeLedger{
		SizeLedger(final int maximum){
			if(maximum<0 || maximum>=Integer.MAX_VALUE-8){throw new IllegalArgumentException("Family-size histogram exceeds array capacity");}
			counts=new int[maximum+1];sums=new long[maximum+1];
		}
		void add(final int size,final int change){
			if(size==0){return;}
			if(size<0 || size>=counts.length || (change!=1 && change!=-1)){throw new IllegalArgumentException("Invalid family-size histogram update");}
			totalCount=Math.addExact(totalCount,change);totalSum=Math.addExact(totalSum,(long)size*change);
			for(int at=size;at<counts.length;){
				counts[at]=Math.addExact(counts[at],change);sums[at]=Math.addExact(sums[at],(long)size*change);
				assert(counts[at]>=0 && sums[at]>=0) : "Only sizes currently present in the frozen or accepted roster may be removed";
				final long next=(long)at+(at&-at);if(next>=counts.length){break;}at=(int)next;
			}
		}
		long top(final int budget){
			if(totalCount<=budget){return totalSum;}
			int remove=totalCount-budget,at=0;
			long bottom=0;
			for(int step=Integer.highestOneBit(counts.length-1);step>0;step>>>=1){
				final long next=(long)at+step;
				if(next<counts.length && counts[(int)next]<=remove){
					at=(int)next;remove-=counts[at];bottom+=sums[at];
				}
			}
			assert(remove==0 || at+1<counts.length) : "The bottom-count search must stop inside the nonempty size range";
			return totalSum-bottom-(long)remove*(at+1);
		}
		final int[] counts;
		final long[] sums;
		int totalCount;
		long totalSum;
	}

	/** Writes all requests in stable order. Missing measured labels are NA, never zero.
	 * Columns through feature_ns are audit/label metadata; only named f_ columns are predictors.
	 */
	private static void writeWideRows(final Row[] rows,final WideProbeFeatures.Output output,final boolean overwrite){
		final ByteStreamWriter writer=new ByteStreamWriter(output.path,overwrite,false,true);writer.start();
		try{
			writer.println(WIDE_HEADER);
			final ByteBuilder line=new ByteBuilder();
			for(Row row:rows){
				if(row.wide==null){throw new IllegalStateException("Every request must retain an explicit numeric feature snapshot");}
				line.clear().append(WideProbeFeatures.SCHEMA).tab().append(output.context).tab().append(row.id).tab()
					.append(row.phase.name().toLowerCase(Locale.ROOT)).tab().append(row.a).tab().append(row.b).tab()
					.append(row.child).tab().append(row.reason).tab().append(row.available ? 1 : 0).tab();
				if(row.available){line.append(row.gain);}else{line.append("NA");}
				line.tab().append(row.featureNanos);row.wide.append(line);writer.println(line);
			}
		}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Wide feature output failed: "+output.path);}}
	}

	/** Writes raw evidence before fitting; final reports append prediction/selection columns to the same raw fields. */
	private static void writeCombinedRows(final Row[] rows,final String path,final boolean overwrite,final double slot,final boolean decisions){
		if(path==null){throw new IllegalArgumentException("Combined rows require an output path");}
		final ByteStreamWriter writer=new ByteStreamWriter(path,overwrite,false,true);writer.start();
		try{
			writer.println(decisions ? COMBINED_DECISION_HEADER : COMBINED_RAW_HEADER);
			final ByteBuilder line=new ByteBuilder();
			for(Row row:rows){
				line.clear().append(row.id).tab().append(row.phase.name().toLowerCase(Locale.ROOT)).tab().append(row.a).tab().append(row.b).tab()
					.append(row.child).tab();
				LiveProbeLayer.appendSlotCost(line,slot).tab().appendSlow(row.threshold).tab().append(row.available ? 1 : 0).tab()
					.append(row.reason).tab().append(row.features!=null ? 1 : 0).tab();
				if(row.features!=null){line.appendSlow(row.features.cohesionDelta).tab().appendSlow(row.features.separationDelta).tab().appendSlow(row.features.sizePressure);}
				else{line.append("NA\tNA\tNA");}
				if(row.available){line.tab().append(row.gain).tab().append(row.rescued).tab().append(row.retained).tab().append(row.affected);}
				else{line.append("\tNA\tNA\tNA\tNA");}
				line.tab();if(row.available && row.refinement!=null){line.append(row.refinement.forcedRetention);}else{line.append("NA");}
				if(decisions){
					line.tab();if(row.predictionAvailable){line.appendSlow(row.prediction);}else{line.append("NA");}
					line.tab().append(row.accepted ? 1 : 0).tab();
					if(row.accepted){
						line.append(row.commitOrder).tab();
						if(row.predictionAvailable){line.appendSlow(row.predictedMarginal);}else{line.append("NA");}
						line.tab().append(row.marginal);
					}else{line.append("NA\tNA\tNA");}
					line.tab().append(row.topBudgetDecreased ? 1 : 0).tab().append(row.emptyChildVeto ? 1 : 0);
				}
				writer.println(line);
			}
		}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Combined report failed: "+path);}}
	}

	/** Completed one-fit/one-commit result; raw evidence remains on disk even when fitting cannot complete. */
	public static final class CombinedResult{
		private CombinedResult(final LiveProbeState state_,final Row[] rows_,final DdlClustProbeRegression.LiveFit fit_,
				final String fitUnavailableReason_,final int featureObservations_,final boolean rankByPrediction_,
				final double fissionThreshold_,final double fusionThreshold_,
				final int slot_,final long topBefore_,final long topAfter_,final long gain_,final long elapsed_){
			state=state_;rows=rows_;fit=fit_;slot=slot_;topBefore=topBefore_;topAfter=topAfter_;gain=gain_;elapsed=elapsed_;
			fitUnavailableReason=fitUnavailableReason_;featureObservations=featureObservations_;rankByPrediction=rankByPrediction_;
			fissionThreshold=fissionThreshold_;fusionThreshold=fusionThreshold_;
			for(Row row:rows){
				final int phase=row.phase==LayerCakeSchedule.Phase.FISSION ? 0 : 1;
				candidates[phase]++;if(row.available){measured[phase]++;}if(row.predictionAvailable){fitted[phase]++;}
				if(row.accepted){accepted[phase]++;}
				if(row.topBudgetDecreased){budgetDecreased++;}
				if(row.emptyChildVeto){emptyChildRejected++;}
			}
			assert((fit==null ? fitted()==0 && fitUnavailableReason!=null && !rankByPrediction :
				fitted()==fit.observations() && fitted()==featureObservations && fitUnavailableReason==null) &&
				topBefore>=0 && topAfter>=0) : "Combined fit and publication must conserve their explicit available-row and signed coverage ledgers";
		}
		public LiveProbeState state(){return state;}
		public boolean fitAvailable(){return fit!=null;}
		public boolean rankedByPrediction(){return rankByPrediction;}
		public double fissionThreshold(){return fissionThreshold;}
		public double fusionThreshold(){return fusionThreshold;}
		/** Null only when a fit is available; missing diagnostics never use NaN. */
		public String fitUnavailableReason(){return fitUnavailableReason;}
		public DdlClustProbeRegression.LiveFit fit(){
			if(fit==null){throw new IllegalStateException("Diagnostic fit unavailable: "+fitUnavailableReason);}
			return fit;
		}
		public int candidates(){return candidates[0]+candidates[1];}
		public int fissionCandidates(){return candidates[0];}
		public int fusionCandidates(){return candidates[1];}
		public int measured(){return measured[0]+measured[1];}
		public int fissionMeasured(){return measured[0];}
		public int fusionMeasured(){return measured[1];}
		public int fitted(){return fitted[0]+fitted[1];}
		public int fissionFitted(){return fitted[0];}
		public int fusionFitted(){return fitted[1];}
		public int accepted(){return accepted[0]+accepted[1];}
		public int fissionAccepted(){return accepted[0];}
		public int fusionAccepted(){return accepted[1];}
		/** Compatibility accessor: top-budget changes no longer veto selected operations. */
		public int topBudgetRejected(){return 0;}
		public int topBudgetDecreased(){return budgetDecreased;}
		public int emptyChildRejected(){return emptyChildRejected;}
		public int slotCost(){return slot;}
		public long topMembersBefore(){return topBefore;}
		public long topMembersAfter(){return topAfter;}
		public long netGain(){return gain;}
		public long elapsedNanos(){return elapsed;}
		public void writeReport(final String path,final boolean overwrite){writeCombinedRows(rows,path,overwrite,slot,true);}
		public void writeFit(final String path,final boolean overwrite){
			if(path==null){throw new IllegalArgumentException("Combined fit needs an output path");}
			final ByteStreamWriter writer=new ByteStreamWriter(path,overwrite,false,true);writer.start();
			try{
				writer.println("metric\tvalue");
				if(fit==null){
					writer.println("status\tunavailable");
					writer.println("reason\t"+fitUnavailableReason);
					writer.println("observations\t"+featureObservations);
					return;
				}
				final String[] names={"observations","intercept","cohesion_coefficient","separation_coefficient","size_coefficient","sse","rmse","r_squared"};
				final double[] values={fit.observations(),fit.intercept(),fit.cohesionCoefficient(),fit.separationCoefficient(),fit.sizeCoefficient(),fit.sse(),fit.rmse(),fit.rSquared()};
				final ByteBuilder line=new ByteBuilder();
				for(int i=0;i<names.length;i++){line.clear().append(names[i]).tab().appendSlow(values[i]);writer.println(line);}
			}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Combined fit output failed: "+path);}}
		}
		private final LiveProbeState state;
		private final Row[] rows;
		private final DdlClustProbeRegression.LiveFit fit;
		private final String fitUnavailableReason;
		private final int featureObservations;
		private final boolean rankByPrediction;
		private final double fissionThreshold,fusionThreshold;
		private final int slot;
		private final long topBefore,topAfter,gain,elapsed;
		private final int[] candidates=new int[2],measured=new int[2],fitted=new int[2],accepted=new int[2];
		private int budgetDecreased,emptyChildRejected;
	}
	public static final String COMBINED_RAW_HEADER="request_id\tphase\tparent_a\tparent_b\tchild_id\tslot_cost\tthreshold\tavailable\treason\tfeature_available"+
		"\tweighted_cohesion_delta\tweighted_separation_delta\tsize_pressure\tmeasured_gain\trescued\tretained\taffected\tforced_retention";
	public static final String COMBINED_DECISION_HEADER=COMBINED_RAW_HEADER+"\tprediction\taccepted\tcommit_order\tpredicted_marginal_gain\tmeasured_marginal_gain\ttop_budget_decreased\tempty_child_veto";

	/** Published state plus scalar per-request decisions; rejected model payloads are released. */
	public static final class Result{
		private Result(LiveProbeState state_,Row[] rows_,int available_,int accepted_,double slot_,double threshold_,double sizeLambda_,
				LayerCakeSchedule.Phase phase_,long gain_,long elapsed_){
			state=state_;rows=rows_;available=available_;accepted=accepted_;slot=slot_;threshold=threshold_;
			sizeLambda=sizeLambda_;
			phase=phase_;gain=gain_;elapsed=elapsed_;
			assert(accepted<=available && available<=rows.length && Double.isFinite(slot) && slot>=0 && elapsed>=0) :
				"Live result counts must describe the completed pass with a finite nonnegative effective price";
		}
		public LiveProbeState state(){return state;}
		public int candidates(){return rows.length;}
		public int available(){return available;}
		public int accepted(){return accepted;}
		public int emptyChildRejected(){
			int count=0;for(Row row:rows){if(row.emptyChildVeto){count++;}}
			return count;
		}
		public double slotCost(){return slot;}
		public double sizeLambda(){return sizeLambda;}
		public double threshold(){return threshold;}
		public long netGain(){return gain;}
		public long elapsedNanos(){return elapsed;}
		public void writeReport(final String path,final boolean overwrite){
			if(path==null){throw new IllegalArgumentException("Live decisions require an output path");}
			final ByteStreamWriter writer=new ByteStreamWriter(path,overwrite,false,true);writer.start();
			try{
				final ByteBuilder summary=new ByteBuilder();
				summary.append("#live_probe_pass\tphase=").append(phase.name().toLowerCase(Locale.ROOT))
					.append("\tcandidates=").append(rows.length).append("\tavailable=").append(available)
					.append("\taccepted=").append(accepted).append("\tslot_cost=");
				LiveProbeLayer.appendSlotCost(summary,slot)
					.append("\tthreshold=").appendSlow(threshold).append("\tcovered_before=").append(state.assigned()-gain)
					.append("\tcovered_after=").append(state.assigned()).append("\tresidual=").append(state.residualOrdinals().length)
					.append("\tfiltered=").append(state.filtered()).append("\tmodels=").append(state.models().length)
					.append("\tactive_models=").append(state.active().length).append("\tnext_family_id=").append(state.firstNewId())
					.append("\tnet_gain=").append(gain).append("\tempty_child_rejected=").append(emptyChildRejected());
				writer.println(summary);
				writer.println("request_id\tphase\tparent_a\tparent_b\tchild_id\tslot_cost\tthreshold\tavailable\treason\tstrict_member_veto\tmeasured_gain\taccepted\tcommit_order\tmarginal_gain\trescued\tretained\taffected\tforced_retention"+(sizeLambda>0 ? SIZE_DECISION_COLUMNS : "")+"\tempty_child_veto");
				final ByteBuilder line=new ByteBuilder();
				for(Row row:rows){
					line.clear().append(row.id).tab().append(phase.name().toLowerCase(Locale.ROOT)).tab().append(row.a).tab()
						.append(row.b).tab().append(row.child).tab();
					LiveProbeLayer.appendSlotCost(line,slot).tab().appendSlow(threshold).tab()
						.append(row.available ? 1 : 0).tab().append(row.reason).tab()
						.append(row.veto ? 1 : 0).tab();
					if(row.available){line.append(row.gain);}else{line.append("NA");}
					line.tab().append(row.accepted ? 1 : 0).tab();
					if(row.accepted){line.append(row.commitOrder).tab().append(row.marginal);}else{line.append("NA\tNA");}
					if(row.available){line.tab().append(row.rescued).tab().append(row.retained).tab().append(row.affected);}
					else{line.append("\tNA\tNA\tNA");}
					line.tab();if(row.refinement==null){line.append("NA");}else{line.append(row.refinement.forcedRetention);}
					appendPressure(line,row);
					line.tab().append(row.emptyChildVeto ? 1 : 0);
					writer.println(line);
				}
			}finally{if(writer.poisonAndWait()){throw new IllegalStateException("Live decision output failed: "+path);}}
		}
		void appendAcceptedFissionChildren(final IntList children){
			if(children==null){throw new IllegalArgumentException("Accepted-child collection needs an output list");}
			if(phase!=LayerCakeSchedule.Phase.FISSION){throw new IllegalStateException("Only fission passes have two child-frontier endpoints");}
			for(Row row:rows){
				if(row.accepted){
					children.add(row.a);
					children.add(row.child);
				}
			}
		}
		void writeLayerRows(final ByteStreamWriter writer,final int wave,final int firstRequestId){
			if(writer==null || wave<1 || firstRequestId<0){
				throw new IllegalArgumentException("Layer rows need a writer, positive wave and nonnegative request base");
			}
			final ByteBuilder line=new ByteBuilder();
			for(Row row:rows){
				final int layerRequest=Math.addExact(firstRequestId,row.id);
				line.clear().append(layerRequest).tab().append(wave).tab().append(row.id).tab()
					.append(row.phase.name().toLowerCase(Locale.ROOT)).tab().append(row.a).tab()
					.append(row.b).tab().append(row.child).tab();
				LiveProbeLayer.appendSlotCost(line,slot).tab()
					.appendSlow(threshold).tab().append(row.available ? 1 : 0).tab()
					.append(row.reason).tab().append(row.veto ? 1 : 0).tab();
				if(row.available){line.append(row.gain);}else{line.append("NA");}
				line.tab().append(row.accepted ? 1 : 0).tab();
				if(row.accepted){line.append(row.commitOrder).tab().append(row.marginal);}
				else{line.append("NA\tNA");}
				if(row.available){line.tab().append(row.rescued).tab().append(row.retained).tab().append(row.affected);}
				else{line.append("\tNA\tNA\tNA");}
				line.tab();if(row.refinement==null){line.append("NA");}else{line.append(row.refinement.forcedRetention);}
				appendPressure(line,row);
				line.tab().append(row.emptyChildVeto ? 1 : 0);
				writer.println(line);
			}
		}
		/** Unavailable measurements are not invented zeros; uncommitted rows have no marginal score. */
		private void appendPressure(final ByteBuilder line,final Row row){
			if(sizeLambda==0){return;}
			if(!row.available){line.append("\tNA\tNA\tNA");return;}
			line.tab().appendSlow(row.sizePressure).tab().appendSlow(row.gain+row.sizePressure).tab();
			if(row.accepted){line.appendSlow(row.marginal+row.sizePressure);}else{line.append("NA");}
		}
		void writeWideRows(final ByteStreamWriter writer,final WideProbeFeatures.Output output,
				final int wave,final int firstRequestId){
			if(writer==null || output==null || wave<1 || firstRequestId<0){
				throw new IllegalArgumentException("Wide layer rows need a writer, output and stable request base");
			}
			final ByteBuilder line=new ByteBuilder();
			for(Row row:rows){
				if(row.wide==null){throw new IllegalStateException("Every layer request must retain an explicit numeric feature snapshot");}
				final int layerRequest=Math.addExact(firstRequestId,row.id);
				line.clear().append(WideProbeFeatures.SCHEMA).tab().append(output.context).tab().append(layerRequest).tab()
					.append(row.phase.name().toLowerCase(Locale.ROOT)).tab().append(row.a).tab()
					.append(row.b).tab().append(row.child).tab().append(row.reason).tab()
					.append(row.available ? 1 : 0).tab();
				if(row.available){line.append(row.gain);}else{line.append("NA");}
				line.tab().append(row.featureNanos);row.wide.append(line);writer.println(line);
			}
		}
		private final LiveProbeState state;
		private final Row[] rows;
		private final int available,accepted;
		private final double slot;
		private final double threshold,sizeLambda;
		private final LayerCakeSchedule.Phase phase;
		private final long gain,elapsed;
	}
	static final String LAYER_DECISION_HEADER="layer_request_id\twave\trequest_id\tphase\tparent_a\tparent_b\tchild_id\tslot_cost\tthreshold\tavailable\treason\tstrict_member_veto\tmeasured_gain\taccepted\tcommit_order\tmarginal_gain\trescued\tretained\taffected\tforced_retention";
	static final String SIZE_DECISION_COLUMNS="\tsize_pressure\tselection_score\tmarginal_score";
	static final String WIDE_HEADER;
	static{
		final ByteBuilder line=new ByteBuilder("schema\tcontext\trequest_id\tphase\tparent_a\tparent_b\tchild_id\tconstruction_reason\tlabel_available\tmeasured_gain\tfeature_ns");
		for(String name:WideProbeFeatures.names()){line.tab().append("f_").append(name);}
		WIDE_HEADER=line.toString();
	}
	private LiveProbeClustering(){}
}
