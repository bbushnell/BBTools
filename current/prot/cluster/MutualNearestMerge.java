package prot.cluster;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.ReadWriteLock;

import ddl.FamilySketchBuilder;
import ddl.MutableFamilyIndex;
import prot.HbmPositionModel;
import template.Accumulator;
import template.ThreadWaiter;

/** Selects disjoint reciprocal neighbors from one frozen family population.
 * All families are indexed, including large families which may block a small
 * family's merge. Only small families need queries. Selection is by pinned
 * identity, then score ratio, then stable ID within the indexed shortlist.
 * The opt-in size-ratio mode queries every family, restricting targets before
 * shortlisting to nonself families within the configured size ratio.
 * Merge thresholds apply after nearest selection, in both directions, without
 * fallback. This class never merges, rebuilds or changes memberships.
 * @author Brian Bushnell, Fischl
 */
public final class MutualNearestMerge{

	/*--------------------------------------------------------------*/
	/*----------------             API              ----------------*/
	/*--------------------------------------------------------------*/

	private MutualNearestMerge(){}

	/** Selects one pass; the caller must keep model payloads unchanged until return.
	 * originalMembers is aligned with models and counts original tails, never the
	 * number of consensus records. shortlist excludes self. Returned indices refer
	 * to the caller's model order; pair order follows stable family IDs.
	 */
	public static Result select(final Family[] models,final int[] originalMembers,
			final FamilySketchBuilder builder,final TailAlignmentGate mergeGate,
			final int maxSize,final int shortlist,final int threads){
		return select(models,originalMembers,builder,mergeGate,maxSize,shortlist,threads,1);
	}

	/** Separately bounds immutable index construction and nearest-neighbor query workers. */
	public static Result select(final Family[] models,final int[] originalMembers,
			final FamilySketchBuilder builder,final TailAlignmentGate mergeGate,
			final int maxSize,final int shortlist,final int threads,final int frozenIndexThreads){
		final Prepared work=new Prepared(models,originalMembers,builder,mergeGate,maxSize,shortlist,threads,frozenIndexThreads,0,null);
		work.run();
		return work.result();
	}

	/** Optional family-query pivot; counts/eligibility, alignment and reciprocity are unchanged. */
	public static Result selectMemberRanked(final Family[] models,final int[] originalMembers,
			final FamilySketchBuilder builder,final TailAlignmentGate mergeGate,
			final int maxSize,final int shortlist,final int threads,final int frozenIndexThreads,
			final MemberRankInputs inputs){
		if(inputs==null){throw new IllegalArgumentException("Member-ranked family queries require their original ordinal cache");}
		final Prepared work=new Prepared(models,originalMembers,builder,mergeGate,maxSize,shortlist,threads,frozenIndexThreads,0,inputs);
		work.run();
		return work.result();
	}

	/** Selects one frozen fusion proposal pass among all size-compatible families.
	 * Both counts are positive original-member counts; maxSizeRatio is inclusive.
	 * Incompatible families cannot consume shortlist capacity or block a compatible
	 * neighbor. Returned pairs still require both directional alignment gates;
	 * the caller must separately accept each union under the shared local objective.
	 */
	public static Result selectWithinSizeRatio(final Family[] models,final int[] originalMembers,
			final FamilySketchBuilder builder,final TailAlignmentGate mergeGate,
			final double maxSizeRatio,final int shortlist,final int threads,final int frozenIndexThreads){
		return selectWithinSizeRatio(models,originalMembers,builder,mergeGate,maxSizeRatio,shortlist,threads,frozenIndexThreads,null);
	}

	/** Non-null inputs select the retained-member pivot; null retains the legacy control. */
	public static Result selectWithinSizeRatio(final Family[] models,final int[] originalMembers,
			final FamilySketchBuilder builder,final TailAlignmentGate mergeGate,
			final double maxSizeRatio,final int shortlist,final int threads,final int frozenIndexThreads,
			final MemberRankInputs inputs){
		if(!Double.isFinite(maxSizeRatio) || maxSizeRatio<1){
			throw new IllegalArgumentException("Fusion size ratio must be finite and at least one");
		}
		final Prepared work=new Prepared(models,originalMembers,builder,mergeGate,
			Integer.MAX_VALUE,shortlist,threads,frozenIndexThreads,maxSizeRatio,inputs);
		work.run();
		return work.result();
	}

	/** Lazy pair validation against one full frozen candidate population.
	 * Each requested input's nearest choice is computed once and cached, including
	 * absence or gate failure. Models must remain unchanged for this object's lifetime.
	 * Construction builds only the index; query never starts a whole-population run.
	 */
	public static final class PairQueries{
		/** Uses the complete size-compatible candidate pool and existing nearest-choice path. */
		public PairQueries(final Family[] models,final int[] originalMembers,
				final FamilySketchBuilder builder,final TailAlignmentGate mergeGate,
				final double maxSizeRatio,final int shortlist,final int indexThreads){
			this(models,originalMembers,builder,mergeGate,maxSizeRatio,shortlist,indexThreads,null);
		}

		/** A ranked session freezes member descriptors once and shares them across worker forks. */
		public PairQueries(final Family[] models,final int[] originalMembers,
				final FamilySketchBuilder builder,final TailAlignmentGate mergeGate,
				final double maxSizeRatio,final int shortlist,final int indexThreads,final MemberRankInputs inputs){
			if(!Double.isFinite(maxSizeRatio) || maxSizeRatio<1){
				throw new IllegalArgumentException("Fusion size ratio must be finite and at least one");
			}
			prepared=new Prepared(models,originalMembers,builder,mergeGate,
				Integer.MAX_VALUE,shortlist,1,indexThreads,maxSizeRatio,inputs);
			inputToDense=new int[prepared.models.length];
			for(int dense=0; dense<inputToDense.length; dense++){
				inputToDense[prepared.denseToInput[dense]]=dense;
			}
			queried=new boolean[inputToDense.length];
			worker=new Worker(prepared);
		}

		/** Private query outcomes and scratch over the same immutable index and roster. */
		public PairQueries fork(){return new PairQueries(this);}
		private PairQueries(final PairQueries source){
			prepared=new Prepared(source.prepared);
			inputToDense=source.inputToDense;
			queried=new boolean[inputToDense.length];
			worker=new Worker(prepared);
		}

		/** Computes at most two new nearest choices, returning only this request's work.
		 * Inputs are positions in the constructor's model array, not stable family IDs.
		 * Reversed and overlapping requests reuse cached endpoints. Calls are serialized
		 * because the scorer, rank scratch and outcome cache have one owner.
		 */
		public synchronized PairChoice query(final int firstInput,final int secondInput){
			if(firstInput<0 || firstInput>=queried.length || secondInput<0 ||
				secondInput>=queried.length || firstInput==secondInput){
				throw new IllegalArgumentException("A pair query needs two distinct valid input positions");
			}
			if(failure!=null){throw new IllegalStateException("A failed pair-query session cannot reuse partial outcomes",failure);}
			final long attempted=worker.attempted,valid=worker.valid;
			int newQueries=0;
			try{
				if(!queried[firstInput]){compute(firstInput); newQueries++;}
				if(!queried[secondInput]){compute(secondInput); newQueries++;}
			}catch(RuntimeException | Error problem){
				failure=problem;
				throw problem;
			}
			final boolean reciprocal=prepared.nearest[firstInput]==secondInput && prepared.nearest[secondInput]==firstInput;
			final boolean gatesPassed=reciprocal && prepared.gatePass[firstInput] && prepared.gatePass[secondInput];
			Pair pair=null;
			if(gatesPassed){
				final boolean firstLeft=prepared.models[firstInput].familyId<prepared.models[secondInput].familyId;
				pair=new Pair(firstLeft ? firstInput : secondInput,firstLeft ? secondInput : firstInput);
			}
			return new PairChoice(pair,reciprocal,gatesPassed,worker.attempted-attempted,worker.valid-valid,newQueries);
		}

		/** Marks an endpoint cached only after the same selector used by full passes succeeds. */
		private void compute(final int input){
			assert(!queried[input]) : "A cached endpoint must never repeat nearest-neighbor alignment work";
			worker.processQuery(inputToDense[input]);
			queried[input]=true;
		}
		private final Prepared prepared;
		private final Worker worker;
		private final int[] inputToDense;
		private final boolean[] queried;
		private Throwable failure;
	}

	/** Immutable requested-pair decision and incremental work counters.
	 * gatesPassed is true only for a reciprocal pair passing both directional gates;
	 * pair is null otherwise. Cached gate rejections or absence return zero new work.
	 */
	public static final class PairChoice{
		private PairChoice(final Pair pair_,final boolean reciprocal_,final boolean gatesPassed_,
				final long attempted_,final long valid_,final int newQueries_){
			assert(newQueries_>=0 && newQueries_<=2 && attempted_>=0 && valid_>=0 && valid_<=attempted_ &&
				(!gatesPassed_ || reciprocal_) && (pair_!=null)==gatesPassed_) :
				"One cached pair request can compute at most two endpoints and publish only a reciprocal gate-passing pair";
			pair=pair_; reciprocal=reciprocal_; gatesPassed=gatesPassed_;
			attemptedAlignments=attempted_; spanValidAlignments=valid_; newQueries=newQueries_;
		}
		public final Pair pair;
		public final boolean reciprocal,gatesPassed;
		public final long attemptedAlignments,spanValidAlignments;
		public final int newQueries;
	}

	/** One accepted pair; left has the lower stable family ID, not necessarily the lower input index. */
	public static final class Pair{
		private Pair(final int left,final int right){leftIndex=left; rightIndex=right;}
		public final int leftIndex,rightIndex;
	}

	/** Completed choices and counters; absence is -1, never a nonfinite numeric sentinel. */
	public static final class Result{
		private Result(final Prepared work,final Pair[] pairs_,final int reciprocal,final int rejected){
			pairs=pairs_; eligibleQueries=work.eligibleDense.length;
			attemptedAlignments=work.attempted; spanValidAlignments=work.valid;
			reciprocalPairs=reciprocal; gateRejectedPairs=rejected;
			nearest=work.nearest; identities=work.identities; ratios=work.ratios; gatePass=work.gatePass;
			familyIds=new int[work.models.length];
			for(int i=0; i<familyIds.length; i++){familyIds[i]=work.models[i].familyId;}
		}
		/** Returns -1 for a legacy unqueried large family or a query without an eligible valid neighbor. */
		public int nearestIndexFor(final int input){return nearest[input];}
		public int nearestFamilyIdFor(final int input){return nearest[input]<0 ? -1 : familyIds[nearest[input]];}
		/** Zero for no neighbor; call nearestIndexFor to distinguish a real zero-identity alignment. */
		public double nearestIdentityFor(final int input){return identities[input];}
		public double nearestRatioFor(final int input){return ratios[input];}
		public boolean nearestPassesGate(final int input){return gatePass[input];}
		public int inputs(){return nearest.length;}
		public final Pair[] pairs;
		public final int eligibleQueries,reciprocalPairs,gateRejectedPairs;
		public final long attemptedAlignments,spanValidAlignments;
		private final int[] nearest,familyIds;
		private final double[] identities,ratios;
		private final boolean[] gatePass;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Frozen work           ----------------*/
	/*--------------------------------------------------------------*/

	/** Dense IDs follow stable-ID order, avoiding arrays proportional to a sparse maximum ID. */
	private static final class Prepared implements Accumulator<Worker>,MutableFamilyIndex.SketchSource{
		/** Only result arrays are mutable; never rebuild or copy the frozen posting index per worker. */
		Prepared(final Prepared source){
			assert(source!=null) : "A query fork requires a completely constructed frozen index";
			models=source.models;originalMembers=source.originalMembers;builder=source.builder;gate=source.gate;
			memberRanks=source.memberRanks;
			hydroWeight=source.hydroWeight;
			maxSize=source.maxSize;shortlist=source.shortlist;maxSizeRatio=source.maxSizeRatio;ratioRestricted=source.ratioRestricted;
			denseToInput=source.denseToInput;eligibleDense=source.eligibleDense;view=source.view;workerCount=1;
			nearest=new int[models.length];Arrays.fill(nearest,-1);
			identities=new double[models.length];ratios=new double[models.length];gatePass=new boolean[models.length];
		}
		Prepared(final Family[] models_,final int[] counts,final FamilySketchBuilder builder_,
				final TailAlignmentGate gate_,final int maxSize_,final int shortlist_,final int threads,
				final int frozenIndexThreads,final double maxSizeRatio_,final MemberRankInputs inputs){
			if(models_==null || counts==null || models_.length!=counts.length || builder_==null || gate_==null ||
				maxSize_<1 || shortlist_<1 || shortlist_==Integer.MAX_VALUE || threads<1 || frozenIndexThreads<1 ||
				gate_.mode==TailAlignmentGate.Mode.FULL_QUERY_UNPINNED){
				throw new IllegalArgumentException("Mutual merging needs matching models/counts, a pinned gate, positive limits and workers");
			}
			models=models_.clone(); originalMembers=counts.clone(); builder=builder_; gate=gate_;
			memberRanks=inputs==null ? null : new MemberRankStatistics.Snapshot[models.length];
			hydroWeight=inputs==null ? MemberRankInputs.DEFAULT_HYDRO_WEIGHT : inputs.hydroWeight;
			maxSize=maxSize_; shortlist=shortlist_; maxSizeRatio=maxSizeRatio_; ratioRestricted=maxSizeRatio>0;
			final long[] ordered=new long[models.length];
			int eligible=0;
			for(int i=0; i<models.length; i++){
				final Family model=models[i];
				if(model==null || model.familyId<0 || originalMembers[i]<1 || model.consensusUnsafe()==null ||
					model.consensusUnsafe().length<gate.minSpan || model.sketchUnsafe()==null ||
					model.sketchUnsafe().length!=builder.width()){
					throw new IllegalArgumentException("Invalid frozen model, original-tail count, span or sketch width at index "+i);
				}
				if(memberRanks!=null){
					if(!model.memberRankingUses(inputs)){
						throw new IllegalArgumentException("Family query descriptors must use the same original ordinal cache at index "+i);
					}
					memberRanks[i]=model.memberRankingSnapshot();
				}
				ordered[i]=((long)model.familyId<<32)|(i&0xffffffffL);
				if(ratioRestricted || originalMembers[i]<maxSize){eligible++;}
			}
			Arrays.sort(ordered);
			denseToInput=new int[models.length]; eligibleDense=new int[eligible];
			int kept=0,previous=-1;
			for(int dense=0; dense<ordered.length; dense++){
				final int input=(int)ordered[dense],id=models[input].familyId;
				if(id==previous){throw new IllegalArgumentException("Duplicate stable family ID: "+id);}
				previous=id; denseToInput[dense]=input;
				if(ratioRestricted || originalMembers[input]<maxSize){eligibleDense[kept++]=dense;}
			}
			assert(kept==eligible) : "Every eligible family must get exactly one frozen nearest-neighbor query";
			// Models are already immutable for this pass: build CSR once instead
			// of inserting a mutable delta and repeatedly compacting it. With
			// no small queries, the empty snapshot avoids all target-row work.
			try(PhaseCostLog.Span phase=PhaseCostLog.begin("family_index")){
				view=frozenIndexThreads==1 ? MutableFamilyIndex.snapshot(builder.width(),eligible==0 ? 0 : models.length,this) :
					MutableFamilyIndex.snapshot(builder.width(),eligible==0 ? 0 : models.length,this,frozenIndexThreads);
			}
			workerCount=Math.min(threads,eligible);
			nearest=new int[models.length]; Arrays.fill(nearest,-1);
			identities=new double[models.length]; ratios=new double[models.length]; gatePass=new boolean[models.length];
		}

		/** Dense snapshot source borrows rows under select's frozen-payload contract. */
		@Override public int size(){return models.length;}
		@Override public char[] sketch(final int dense){return models[denseToInput[dense]].sketchUnsafe();}

		/** Joins all started workers before returning or propagating a failure. */
		void run(){
			try(PhaseCostLog.Span phase=PhaseCostLog.begin("nearest_neighbor_comparison")){
				runWorkers();
			}
		}

		/** Indexed shortlists and their alignments, not an exhaustive all-pairs DP pass. */
		private void runWorkers(){
			final ArrayList<Worker> workers=new ArrayList<Worker>(workerCount);
			for(int i=0; i<workerCount; i++){workers.add(new Worker(this));}
			try{
				if(!ThreadWaiter.startAndWait(workers,this)){fail(new IllegalStateException("Mutual-neighbor workers did not complete"));}
			}catch(Throwable cause){
				fail(cause);
				final ArrayList<Worker> started=new ArrayList<Worker>();
				for(Worker worker:workers){if(worker.getState()!=Thread.State.NEW){started.add(worker);}}
				ThreadWaiter.waitForThreadsToFinish(started);
			}
			if(failure!=null){throw new IllegalStateException("Mutual-neighbor pass failed before pair publication",failure);}
			if(done!=eligibleDense.length){throw new IllegalStateException("Mutual-neighbor pass omitted eligible queries");}
		}

		/** Each reciprocal edge contributes exactly once; no connected-component closure is taken. */
		Result result(){
			assert(done==eligibleDense.length && failure==null) : "Pairs require all choices from the same completed frozen pass";
			final ArrayList<Pair> pairs=new ArrayList<Pair>();
			int reciprocal=0,rejected=0;
			for(int dense:eligibleDense){
				final int a=denseToInput[dense],b=nearest[a];
				if(b<0 || !ratioRestricted && originalMembers[b]>=maxSize ||
					models[a].familyId>=models[b].familyId || nearest[b]!=a){continue;}
				assert(!ratioRestricted || withinRatio(originalMembers[a],originalMembers[b],maxSizeRatio)) :
					"Filtered reciprocal choices must both belong to the inclusive size-compatible search population";
				reciprocal++;
				if(gatePass[a] && gatePass[b]){pairs.add(new Pair(a,b));}
				else{rejected++;}
			}
			assert(reciprocal==pairs.size()+rejected) : "Every size-eligible reciprocal pair must pass or fail the bidirectional gate";
			return new Result(this,pairs.toArray(new Pair[pairs.size()]),reciprocal,rejected);
		}

		/** Batch claims keep synchronization outside the per-alignment hot path. */
		int claim(){return (int)Math.min(next.getAndAdd(BATCH),eligibleDense.length);}
		synchronized void fail(final Throwable cause){if(failure==null){failure=cause;}}
		@Override public void accumulate(final Worker worker){
			if(!worker.success){fail(new IllegalStateException("Mutual-neighbor worker ended before completing its batch"));}
			done=Math.addExact(done,worker.done);
			attempted=Math.addExact(attempted,worker.attempted); valid=Math.addExact(valid,worker.valid);
		}
		@Override public boolean success(){return failure==null && done==eligibleDense.length;}
		@Override public ReadWriteLock rwlock(){return null;}// Ordered joins establish visibility.
		final Family[] models;
		final MemberRankStatistics.Snapshot[] memberRanks;
		final int[] originalMembers,denseToInput,eligibleDense,nearest;
		final double[] identities,ratios;
		final boolean[] gatePass;
		final FamilySketchBuilder builder;
		final TailAlignmentGate gate;
		final MutableFamilyIndex.FrozenView view;
		final int maxSize,shortlist,workerCount;
		final double maxSizeRatio,hydroWeight;
		final boolean ratioRestricted;
		final AtomicLong next=new AtomicLong();
		volatile Throwable failure;
		int done;
		long attempted,valid;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Workers            ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns its sketch builder, scorer and index scratch; writes only its claimed result slots. */
	private static final class Worker extends Thread implements MutableFamilyIndex.Cosine,MutableFamilyIndex.Eligibility{
		Worker(final Prepared parent_){
			parent=parent_; scorer=new HbmPositionModel.ConsensusScorer(parent.gate.matrix);
			builder=parent.builder.fork(); querySketch=new char[builder.width()];
			ranking=parent.memberRanks==null ? this : new MutableFamilyIndex.Ranking(){
				@Override public double score(final int dense){return Worker.this.score(dense);}
				@Override public double score(final int dense,final int hits){
					final MemberRankStatistics.Snapshot query=parent.memberRanks[parent.denseToInput[queryDense]];
					final MemberRankStatistics.Snapshot target=parent.memberRanks[parent.denseToInput[dense]];
					final double ratio=Math.min(query.medianLength,target.medianLength)/Math.max(query.medianLength,target.medianLength);
					final double cosine=Worker.this.score(dense),similarity=1-Math.abs(query.meanKd-target.meanKd)/9;
					assert(hits>=0 && ratio>0 && ratio<=1 && cosine>=0 && cosine<=1 && similarity>=0 && similarity<=1) :
						"Frozen member medians/KD and native cosine must obey the complete-score index bounds";
					return hits*Math.sqrt(Math.sqrt(ratio))+1.5*cosine+parent.hydroWeight*similarity;
				}
				@Override public double lowerBound(final int hits){return 0;}
				@Override public double upperBound(final int hits){return hits+1.5+parent.hydroWeight;}
			};
			candidateIds=new int[(int)Math.min((long)parent.shortlist+(parent.ratioRestricted ? 0 : 1),Math.max(1,parent.models.length))];
			selectionGate=new TailAlignmentGate(parent.gate.mode,parent.gate.denominator,
				parent.gate.maxLength,parent.gate.minSpan,0,0,null,parent.gate.matrix,parent.gate.minPaired);
		}
		@Override public void run(){
			try{
				for(int start=parent.claim(); start<parent.eligibleDense.length && parent.failure==null; start=parent.claim()){
					final int stop=(int)Math.min((long)start+BATCH,parent.eligibleDense.length);
					for(int i=start; i<stop && parent.failure==null; i++){processQuery(parent.eligibleDense[i]); done++;}
				}
				success=parent.failure==null;
			}catch(Throwable cause){parent.fail(cause);}
		}

		/** Nearest is selected without merge thresholds so a failing closer model still blocks fallback. */
		void processQuery(final int dense){
			final int input=parent.denseToInput[dense];
			queryDense=dense;
			final byte[] query=parent.models[input].consensusUnsafe();
			builder.sketchEncoded(query,querySketch);
			Arrays.fill(monomers,0); norm=0;
			for(byte aa:query){if(aa>=0 && aa<20){monomers[aa]++;}}
			for(double count:monomers){norm+=count*count;}
			final int candidates=parent.ratioRestricted ?
				parent.view.rankIndexed(querySketch,candidateIds.length,candidateIds,null,ranking,workspace,this) :
				parent.view.rankIndexed(querySketch,candidateIds.length,candidateIds,null,ranking,workspace);
			int kept=0,winner=-1;
			double identity=0,ratio=0;
			for(int i=0; i<candidates && kept<parent.shortlist; i++){
				final int otherDense=candidateIds[i];
				if(otherDense==dense){continue;}
				kept++; attempted++;
				final int other=parent.denseToInput[otherDense];
				final byte[] reference=parent.models[other].consensusUnsafe();
				scorer.resetConsensus(reference);
				final TailAlignmentDecision decision=selectionGate.decide(query,reference,scorer);
				if(!decision.accepted()){continue;}
				valid++;
				if(winner<0 || decision.identity>identity || decision.identity==identity &&
					(decision.scoreRatio>ratio || decision.scoreRatio==ratio &&
					parent.models[other].familyId<parent.models[winner].familyId)){
					winner=other; identity=decision.identity; ratio=decision.scoreRatio;
				}
			}
			assert(kept<=parent.shortlist && winner!=input) : "Self exclusion must precede the nonself shortlist quota and winner selection";
			parent.nearest[input]=winner; parent.identities[input]=identity; parent.ratios[input]=ratio;
			parent.gatePass[input]=winner>=0 && identity>=parent.gate.minIdentity && ratio>=parent.gate.minScoreRatio;
		}

		/** Query-local restriction is evaluated before either positive-hit or zero-hit heap admission. */
		@Override public boolean accept(final int dense){
			assert(queryDense>=0 && queryDense<parent.models.length) : "The worker must select its frozen query before candidate filtering";
			return dense!=queryDense && withinRatio(parent.originalMembers[parent.denseToInput[queryDense]],
				parent.originalMembers[parent.denseToInput[dense]],parent.maxSizeRatio);
		}

		/** Same native monomer cosine used by the ordinary indexed neighbor report. */
		@Override public double score(final int dense){
			final float[] target=parent.models[parent.denseToInput[dense]].monomersUnsafe();
			double dot=0,targetNorm=0;
			for(int i=0; i<20; i++){dot+=monomers[i]*target[i]; targetNorm+=target[i]*(double)target[i];}
			return norm==0 || targetNorm==0 ? 0 : Math.min(1,Math.max(0,dot/Math.sqrt(norm*targetNorm)));
		}
		final Prepared parent;
		final MutableFamilyIndex.Cosine ranking;
		final FamilySketchBuilder builder;
		final TailAlignmentGate selectionGate;
		final HbmPositionModel.ConsensusScorer scorer;
		final MutableFamilyIndex.Workspace workspace=new MutableFamilyIndex.Workspace();
		final char[] querySketch;
		final int[] candidateIds;
		final double[] monomers=new double[20];
		double norm;
		int done,queryDense=-1;
		long attempted,valid;
		boolean success;
	}

	/** Division avoids overflowing a large configured ratio times an original-member count. */
	private static boolean withinRatio(final int a,final int b,final double ratio){
		assert(a>0 && b>0 && Double.isFinite(ratio) && ratio>=1) : "Validated positive counts and a finite ratio define fusion eligibility";
		return Math.max(a,b)/(double)Math.min(a,b)<=ratio;
	}
	private static final int BATCH=32;
}
