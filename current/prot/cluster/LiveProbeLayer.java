package prot.cluster;

import java.math.BigDecimal;
import java.util.Arrays;
import java.util.Locale;

import ddl.FamilySketchBuilder;
import fileIO.ByteStreamWriter;
import prot.ProteinSequence;
import structures.ByteBuilder;
import structures.IntList;

/** One native layer: recursive accepted-child fission waves, then one fusion pass.
 * B is fixed from the layer-start population. Reports are written pass-by-pass
 * as soon as each wave publishes, while only the current state is retained.
 * @author Brian Bushnell, Yelan
 */
public final class LiveProbeLayer{

	/*--------------------------------------------------------------*/
	/*----------------             API              ----------------*/
	/*--------------------------------------------------------------*/

	/** Compatibility entry point for the fraction-priced layer API. */
	public static Result run(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final double slotFraction,final Outputs outputs){
		return run(state,sequences,builder,memberGate,mergeGate,filter,minimumLength,objective,options,
			slotFraction,PriceMode.FRACTION,-1,DEFAULT_EDGE_FRACTION,outputs);
	}

	/** Uses one layer-start price for every wave; floor mode receives the below-edge prior price. */
	public static Result run(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final double slotFraction,final PriceMode priceMode,
			final int setupSlotCost,final double edgeFraction,final Outputs outputs){
		return run(state,sequences,builder,memberGate,mergeGate,filter,minimumLength,objective,options,
			slotFraction,priceMode,setupSlotCost,edgeFraction,BelowMode.CONSTANT,outputs);
	}

	/** Below-edge scaling changes effective price only; source family counts remain integers. */
	public static Result run(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final double slotFraction,final PriceMode priceMode,
			final int setupSlotCost,final double edgeFraction,final BelowMode belowMode,final Outputs outputs){
		return run(state,sequences,builder,memberGate,mergeGate,filter,minimumLength,objective,options,
			slotFraction,priceMode,setupSlotCost,edgeFraction,belowMode,0,outputs);
	}

	/** Applies a fixed size objective to measured proposals without changing their raw gains. */
	public static Result run(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final double slotFraction,final PriceMode priceMode,
			final int setupSlotCost,final double edgeFraction,final BelowMode belowMode,final double sizeLambda,final Outputs outputs){
		return run(state,sequences,builder,memberGate,mergeGate,filter,minimumLength,objective,options,
			slotFraction,priceMode,setupSlotCost,edgeFraction,belowMode,sizeLambda,1,1,outputs);
	}

	/** Rank-capped layer mode starts from the largest fission parents and only
	 * admits smallest-family fusion endpoints. Rank ties are resolved by ID.
	 */
	public static Result run(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final double slotFraction,final PriceMode priceMode,
			final int setupSlotCost,final double edgeFraction,final BelowMode belowMode,final double sizeLambda,
			final double fissionTopFraction,final double fusionBottomFraction,final Outputs outputs){
		return run(state,sequences,builder,memberGate,mergeGate,filter,minimumLength,objective,options,
			slotFraction,priceMode,setupSlotCost,edgeFraction,belowMode,sizeLambda,fissionTopFraction,
			fusionBottomFraction,outputs,null);
	}

	/** An explicit fission frontier disables accepted-child recursion because accepted
	 * splits consume fusion-product eligibility until a later fusion restores it.
	 */
	public static Result run(final LiveProbeState state,final ProteinSequence[] sequences,
			final FamilySketchBuilder builder,final TailAlignmentGate memberGate,final TailAlignmentGate mergeGate,
			final CandidateFilter filter,final int minimumLength,final LocalFamilyObjective objective,
			final LayerCakeClustering.Options options,final double slotFraction,final PriceMode priceMode,
			final int setupSlotCost,final double edgeFraction,final BelowMode belowMode,final double sizeLambda,
			final double fissionTopFraction,final double fusionBottomFraction,final Outputs outputs,
			final int[] fissionParentIds){
		if(!Double.isFinite(sizeLambda) || sizeLambda<0){throw new IllegalArgumentException("sizelambda must be finite and nonnegative");}
		validateEligibilityFraction(fissionTopFraction,"fissiontopfraction");
		validateEligibilityFraction(fusionBottomFraction,"fusionbottomfraction");
		if(fissionParentIds!=null && fissionTopFraction!=1){
			throw new IllegalArgumentException("Explicit fission parents require fissiontopfraction=1");
		}
		if(state==null || sequences==null || builder==null || memberGate==null || mergeGate==null || filter==null ||
			objective==null || options==null || !Double.isFinite(slotFraction) || slotFraction<0 ||
			minimumLength<memberGate.minSpan || sequences.length!=state.sequenceCount() || options.localRefineRounds<1){
			throw new IllegalArgumentException("Layer live probing needs a matching state, finite fraction, gates and refinement");
		}
		if(priceMode==null || belowMode==null || !Double.isFinite(edgeFraction) || edgeFraction<=0 || setupSlotCost< -1 ||
			priceMode==PriceMode.FLOOR && setupSlotCost<0){
			throw new IllegalArgumentException("Layer pricing needs a finite positive edge and an explicit nonnegative setup price in floor mode");
		}
		final Outputs out=outputs==null ? new Outputs(null,null,null,null,true) : outputs;
		final long start=System.nanoTime();
		final Family[] initial=state.active().clone();
		if(initial.length<1){throw new IllegalArgumentException("A live layer needs at least one active family for priced B");}
		Arrays.sort(initial,(a,b)->Integer.compare(a.familyId,b.familyId));
		final int budgetSlot=LayerCakeSchedule.summarize(initial,options.familyBudget,sequences.length).slotCost;
		final int pricedRank=priceMode==PriceMode.FLOOR ? options.familyBudget :
			LayerCakeSchedule.pricedRank(initial.length,options.familyBudget,options.pricedFraction);
		final double slot=priceMode==PriceMode.FLOOR ?
			conditionalSlotCost(initial.length,options.familyBudget,budgetSlot,setupSlotCost,edgeFraction,belowMode) :
			pricedRank==options.familyBudget ? budgetSlot : LayerCakeSchedule.summarize(initial,pricedRank,sequences.length).slotCost;
		final double fissionThreshold=slot*slotFraction, fusionThreshold=-fissionThreshold;
		if(!Double.isFinite(fissionThreshold) || !Double.isFinite(fusionThreshold)){
			throw new IllegalArgumentException("Live layer signed thresholds overflowed");
		}
		LiveProbeState current=state;
		LiveProbeState afterFission=null;
		final IntList frontier=new IntList(initial.length);
		if(fissionParentIds==null){addLargestFamilies(initial,fissionTopFraction,frontier);}
		else{for(int id:fissionParentIds){if(id<0){throw new IllegalArgumentException("Explicit fission parents must be nonnegative");}frontier.add(id);}}
		final boolean recursiveFission=fissionParentIds==null && fissionTopFraction==1;
		int waves=0, fissionCandidates=0, fissionAvailable=0, acceptedFissions=0, acceptedFusions=0, fusionCandidates=0, fusionAvailable=0;
		long fissionGain=0, fusionGain=0, fissionElapsed=0, fusionElapsed=0;
		int nextRequest=0;
		final Writers writers=out.open(sizeLambda>0);
		try{
			try(final LocalProbePool.ResidualCache fissionPool=recursiveFission ?
						new LocalProbePool.ResidualCache(sequences,builder,options.threads) : null;
					PhaseCostLog.Span phase=PhaseCostLog.begin("fission_waves")){
				for(;frontier.size>0;){
					frontier.sort();
					final LiveProbeClustering.Requests requests=LiveProbeClustering.generateFissions(current,frontier.toArray(),options);
					if(requests.size()==0){break;}
					final LiveProbeClustering.Result result=LiveProbeClustering.runRequests(current,sequences,builder,
						memberGate,mergeGate,filter,minimumLength,objective,options,LayerCakeSchedule.Phase.FISSION,
						requests,slot,fissionThreshold,sizeLambda,out.wide,System.nanoTime(),fissionPool);
					waves++;
					fissionCandidates=Math.addExact(fissionCandidates,result.candidates());
					fissionAvailable=Math.addExact(fissionAvailable,result.available());
					acceptedFissions=Math.addExact(acceptedFissions,result.accepted());
					fissionGain=Math.addExact(fissionGain,result.netGain());fissionElapsed=Math.addExact(fissionElapsed,result.elapsedNanos());
					if(writers.fission!=null){result.writeLayerRows(writers.fission,waves,nextRequest);}
					if(writers.wide!=null){result.writeWideRows(writers.wide,out.wide,waves,nextRequest);}
					nextRequest=Math.addExact(nextRequest,result.candidates());
					if(writers.wave!=null){writePass(writers.wave,waves,LayerCakeSchedule.Phase.FISSION,result);}
					frontier.clear();
					if(recursiveFission){result.appendAcceptedFissionChildren(frontier);}
					current=result.state();
				}
			}
			afterFission=current;
			try(PhaseCostLog.Span phase=PhaseCostLog.begin("fusion")){
				final int[] fusionEligible=rankedFamilyIds(current.active(),fusionBottomFraction,false);
				final LiveProbeClustering.Requests requests;
				try(PhaseCostLog.Span candidates=PhaseCostLog.begin("fusion_candidates")){
					requests=LiveProbeClustering.generateFusions(current,builder,mergeGate,options,fusionEligible);
				}
				final LiveProbeClustering.Result fusion=LiveProbeClustering.runRequests(current,sequences,builder,
					memberGate,mergeGate,filter,minimumLength,objective,options,LayerCakeSchedule.Phase.FUSION,
					requests,slot,fusionThreshold,sizeLambda,out.wide,System.nanoTime());
				fusionCandidates=fusion.candidates();fusionAvailable=fusion.available();acceptedFusions=fusion.accepted();
				fusionGain=fusion.netGain();fusionElapsed=fusion.elapsedNanos();
				if(writers.fusion!=null){fusion.writeLayerRows(writers.fusion,1,nextRequest);}
				if(writers.wide!=null){fusion.writeWideRows(writers.wide,out.wide,1,nextRequest);}
				if(writers.wave!=null){writePass(writers.wave,1,LayerCakeSchedule.Phase.FUSION,fusion);}
				current=fusion.state();
			}
		}finally{writers.close();}
		return new Result(current,afterFission,initial.length,pricedRank,options.pricedFraction,slot,slotFraction,
			priceMode,setupSlotCost,budgetSlot,edgeFraction,belowMode,sizeLambda,
			waves,fissionCandidates,fusionCandidates,fissionAvailable,fusionAvailable,acceptedFissions,
			acceptedFusions,fissionGain,fusionGain,fissionElapsed,fusionElapsed,System.nanoTime()-start);
	}

	/** Above the exact decimal population edge use the budget-rank size; equality and below use the prior B.
	 * This is a conditional fallback, not max(prior,current), and either valid price may be zero.
	 */
	public static int conditionalSlotCost(final int positiveFamilies,final int familyBudget,
			final int budgetSlotCost,final int setupSlotCost,final double edgeFraction){
		return (int)conditionalSlotCost(positiveFamilies,familyBudget,budgetSlotCost,setupSlotCost,edgeFraction,BelowMode.CONSTANT);
	}

	/** The ratio is not clamped within the edge band; these are the declared linear/sqrt equations. */
	public static double conditionalSlotCost(final int positiveFamilies,final int familyBudget,
			final int budgetSlotCost,final int setupSlotCost,final double edgeFraction,final BelowMode belowMode){
		if(positiveFamilies<1 || familyBudget<1 || budgetSlotCost<0 || setupSlotCost<0 ||
			belowMode==null || !Double.isFinite(edgeFraction) || edgeFraction<=0 ||
			(positiveFamilies<familyBudget ? budgetSlotCost!=0 : budgetSlotCost==0)){
			throw new IllegalArgumentException("Conditional pricing requires positive counts, a finite positive edge and valid nonnegative rank/setup prices");
		}
		if(aboveEdge(positiveFamilies,familyBudget,edgeFraction)){return budgetSlotCost;}
		final double ratio=positiveFamilies/(double)familyBudget;
		final double price=belowMode==BelowMode.SQRT ? setupSlotCost*Math.sqrt(ratio) :
			(belowMode==BelowMode.LINEAR || belowMode==BelowMode.PRIORSCALED ? setupSlotCost*ratio : setupSlotCost);
		if(!Double.isFinite(price) || price<0){throw new IllegalArgumentException("Effective below-edge price must remain finite and nonnegative");}
		return price;
	}

	/** Exact decimal edge predicate shared by cake dispatch and native layer pricing. */
	public static boolean aboveEdge(final int positiveFamilies,final int familyBudget,final double edgeFraction){
		if(positiveFamilies<1 || familyBudget<1 || !Double.isFinite(edgeFraction) || edgeFraction<=0){
			throw new IllegalArgumentException("A population edge comparison needs positive counts and a finite positive edge");
		}
		final BigDecimal edge=BigDecimal.valueOf(edgeFraction).multiply(BigDecimal.valueOf(familyBudget));
		return BigDecimal.valueOf(positiveFamilies).compareTo(edge)>0;
	}

	/** Preserves legacy integer report text while retaining every fractional price digit. */
	public static ByteBuilder appendSlotCost(final ByteBuilder line,final double price){
		if(line==null || !Double.isFinite(price) || price<0){throw new IllegalArgumentException("A reported slot price must be finite and nonnegative");}
		if(price==Math.rint(price) && price<Long.MAX_VALUE){line.append((long)price);}
		else{line.appendSlow(price);}
		return line;
	}

	private static void validateEligibilityFraction(final double fraction,final String label){
		if(!Double.isFinite(fraction) || fraction<0 || fraction>1){
			throw new IllegalArgumentException(label+" must be finite and in [0,1]");
		}
	}

	private static void addLargestFamilies(final Family[] initial,final double fraction,final IntList frontier){
		final int[] ids=rankedFamilyIds(initial,fraction,true);
		for(int id:ids){frontier.add(id);}
	}

	private static int[] rankedFamilyIds(final Family[] families,final double fraction,final boolean largest){
		final int count=(int)Math.ceil(families.length*fraction);
		final Family[] ranked=families.clone();
		Arrays.sort(ranked,(a,b)->{
			final int size=largest ? Integer.compare(b.size(),a.size()) : Integer.compare(a.size(),b.size());
			return size!=0 ? size : Integer.compare(a.familyId,b.familyId);
		});
		final int[] ids=new int[count];
		for(int i=0;i<count;i++){ids[i]=ranked[i].familyId;}
		return ids;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Reports            ----------------*/
	/*--------------------------------------------------------------*/

	private static void writePass(final ByteStreamWriter writer,final int wave,final LayerCakeSchedule.Phase phase,
			final LiveProbeClustering.Result result){
		final ByteBuilder line=new ByteBuilder();
		line.append(phase.name().toLowerCase(Locale.ROOT)).tab().append(wave).tab().append(result.candidates()).tab()
			.append(result.available()).tab().append(result.accepted()).tab();
		appendSlotCost(line,result.slotCost()).tab()
			.appendSlow(result.threshold()).tab().append(result.netGain()).tab().append(result.elapsedNanos()).tab()
			.append(result.state().assigned()).tab().append(result.state().residualOrdinals().length)
			.tab().append(result.state().models().length).tab().append(result.state().active().length)
			.tab().append(result.state().firstNewId());
		writer.println(line);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Immutable records      ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Outputs{
		public Outputs(final String fissionReport_,final String fusionReport_,final String waveReport_,
				final WideProbeFeatures.Output wide_,final boolean overwrite_){
			fissionReport=fissionReport_;fusionReport=fusionReport_;waveReport=waveReport_;wide=wide_;overwrite=overwrite_;
			final String widePath=wide==null ? null : wide.path;
			if(same(fissionReport,fusionReport) || same(fissionReport,waveReport) || same(fusionReport,waveReport) ||
				same(fissionReport,widePath) || same(fusionReport,widePath) || same(waveReport,widePath)){
				throw new IllegalArgumentException("Layer reports need distinct output paths");
			}
		}
		Writers open(final boolean pressure){
			ByteStreamWriter fission=null,fusion=null,waveWriter=null,wideWriter=null;
			try{
				fission=open(fissionReport);fusion=open(fusionReport);waveWriter=open(waveReport);
				wideWriter=wide==null ? null : open(wide.path);
				final String decisionHeader=LiveProbeClustering.LAYER_DECISION_HEADER+
					(pressure ? LiveProbeClustering.SIZE_DECISION_COLUMNS : "")+"\tempty_child_veto";
				if(fission!=null){fission.println(decisionHeader);}
				if(fusion!=null){fusion.println(decisionHeader);}
				if(waveWriter!=null){waveWriter.println(WAVE_HEADER);}
				if(wideWriter!=null){wideWriter.println(LiveProbeClustering.WIDE_HEADER);}
				final Writers writers=new Writers(fission,fusion,waveWriter,wideWriter);
				fission=fusion=waveWriter=wideWriter=null;
				return writers;
			}catch(Throwable failure){
				closeOpen(fission,failure);closeOpen(fusion,failure);closeOpen(waveWriter,failure);closeOpen(wideWriter,failure);
				if(failure instanceof RuntimeException){throw (RuntimeException)failure;}
				if(failure instanceof Error){throw (Error)failure;}
				throw new RuntimeException(failure);
			}
		}
		private ByteStreamWriter open(final String path){
			if(path==null){return null;}
			final ByteStreamWriter writer=new ByteStreamWriter(path,overwrite,false,true);
			writer.start();return writer;
		}
		private static boolean same(final String a,final String b){return a!=null && b!=null && a.equals(b);}
		private static void closeOpen(final ByteStreamWriter writer,final Throwable failure){
			if(writer!=null){
				try{
					if(writer.poisonAndWait()){
						failure.addSuppressed(new IllegalStateException("Live layer report output failed while recovering from open failure"));
					}
				}catch(Throwable cleanup){failure.addSuppressed(cleanup);}
			}
		}
		public final String fissionReport,fusionReport,waveReport;
		public final WideProbeFeatures.Output wide;
		public final boolean overwrite;
	}

	public static final class Result{
		private Result(final LiveProbeState state_,final LiveProbeState afterFission_,final int positiveCount_,
				final int pricedRank_,final double pricedFraction_,final double slot_,final double fraction_,
				final PriceMode priceMode_,final int setupSlot_,final int budgetSlot_,final double edgeFraction_,final BelowMode belowMode_,final double sizeLambda_,
				final int waves_,final int fissionCandidates_,final int fusionCandidates_,
				final int fissionAvailable_,final int fusionAvailable_,final int acceptedFissions_,
				final int acceptedFusions_,final long fissionGain_,final long fusionGain_,final long fissionElapsed_,
				final long fusionElapsed_,final long elapsed_){
			state=state_;afterFission=afterFission_;positiveCount=positiveCount_;pricedRank=pricedRank_;
			pricedFraction=pricedFraction_;slot=slot_;fraction=fraction_;waves=waves_;
			priceMode=priceMode_;setupSlot=setupSlot_;budgetSlot=budgetSlot_;edgeFraction=edgeFraction_;belowMode=belowMode_;
			sizeLambda=sizeLambda_;
			fissionCandidates=fissionCandidates_;fusionCandidates=fusionCandidates_;
			fissionAvailable=fissionAvailable_;fusionAvailable=fusionAvailable_;
			acceptedFissions=acceptedFissions_;acceptedFusions=acceptedFusions_;
			fissionGain=fissionGain_;fusionGain=fusionGain_;
			fissionElapsed=fissionElapsed_;fusionElapsed=fusionElapsed_;elapsed=elapsed_;
			assert(state!=null && afterFission!=null && positiveCount>0 && pricedRank>0 && Double.isFinite(slot) && slot>=0 && budgetSlot>=0 && waves>=0) :
				"Layer result needs both publication states and nonnegative counters; an unfilled setup may have zero price";
		}
		public LiveProbeState state(){return state;}
		public LiveProbeState afterFissionState(){return afterFission;}
		public int initialPositiveCount(){return positiveCount;}
		public int initialPricedRank(){return pricedRank;}
		public double initialPricedFraction(){return pricedFraction;}
		public double initialSlotCost(){return slot;}
		public PriceMode priceMode(){return priceMode;}
		public int setupSlotCost(){return setupSlot;}
		public int initialBudgetSlotCost(){return budgetSlot;}
		public double edgeFraction(){return edgeFraction;}
		public BelowMode belowMode(){return belowMode;}
		public double sizeLambda(){return sizeLambda;}
		public double slotFraction(){return fraction;}
		public double fissionThreshold(){return slot*fraction;}
		public double fusionThreshold(){return -fissionThreshold();}
		public int fissionWaves(){return waves;}
		public int fissionCandidates(){return fissionCandidates;}
		public int fusionCandidates(){return fusionCandidates;}
		public int fissionAvailable(){return fissionAvailable;}
		public int fusionAvailable(){return fusionAvailable;}
		public int acceptedFissions(){return acceptedFissions;}
		public int acceptedFusions(){return acceptedFusions;}
		public long fissionNetGain(){return fissionGain;}
		public long fusionNetGain(){return fusionGain;}
		public long fissionElapsedNanos(){return fissionElapsed;}
		public long fusionElapsedNanos(){return fusionElapsed;}
		public long elapsedNanos(){return elapsed;}
		private final LiveProbeState state,afterFission;
		private final int positiveCount,pricedRank,waves,fissionCandidates,fusionCandidates,fissionAvailable,fusionAvailable;
		private final int acceptedFissions,acceptedFusions;
		private final double pricedFraction,fraction,slot;
		private final PriceMode priceMode;
		private final BelowMode belowMode;
		private final int setupSlot,budgetSlot;
		private final double edgeFraction,sizeLambda;
		private final long fissionGain,fusionGain,fissionElapsed,fusionElapsed,elapsed;
	}

	private static final class Writers{
		Writers(final ByteStreamWriter fission_,final ByteStreamWriter fusion_,final ByteStreamWriter wave_,
				final ByteStreamWriter wide_){
			fission=fission_;fusion=fusion_;wave=wave_;wide=wide_;
		}
		void close(){
			boolean error=false;
			if(fission!=null){error|=fission.poisonAndWait();}
			if(fusion!=null){error|=fusion.poisonAndWait();}
			if(wave!=null){error|=wave.poisonAndWait();}
			if(wide!=null){error|=wide.poisonAndWait();}
			if(error){throw new IllegalStateException("Live layer report output failed");}
		}
		final ByteStreamWriter fission,fusion,wave,wide;
	}

	private static final String WAVE_HEADER="phase\twave\tcandidates\tavailable\taccepted\tslot_cost\tthreshold\tnet_gain\telapsed_ns\tcovered_after\tresidual\tmodels\tactive_models\tnext_family_id";

	public enum PriceMode{FLOOR,FRACTION}
	public enum BelowMode{CONSTANT,LINEAR,SQRT,PRIOR,PRIORSCALED}
	public static final double DEFAULT_EDGE_FRACTION=1.01;
	public static final double DEFAULT_SIZE_LAMBDA=0.25;
	private LiveProbeLayer(){}
}
