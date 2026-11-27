package prok;

import java.util.Arrays;

import ml.CellNet;

/** Pooled-member positional inference after ordinary alignment trimming.
 * Each caller worker owns its network copies and feature scratch. Both ends
 * are predicted from the SAME post-trim span, optionally centered first by
 * supported training medians. An invalid center retains the raw pair without
 * inference; invalid NN geometry retains the applicable raw or centered base.
 * No detection score, selected model or genomic bases are changed.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaPositionalBoundaryRefiner implements NcrnaBoundaryRefiner{
	/** Control endpoints are state, not nonfinite numeric thresholds. */
	public enum GatePolicy{THRESHOLD, NEVER, ALWAYS}

	/** Immutable routing installed on this exact loaded family before workers start. */
	public static final class Resources implements NcrnaBoundaryRefiner.Factory{
		public Resources(final NcrnaFamily family_, final RrnaPositionalKmerTable.Table five,
				final RrnaPositionalKmerTable.Table three, final CellNet fiveNet, final CellNet threeNet, final float margin_){
			this(family_, new RrnaPositionalKmerTable.Table[]{five}, new RrnaPositionalKmerTable.Table[]{three}, fiveNet, threeNet, margin_);
		}

		/** Binds either one ordinary table/end or three ordered K9 probes/end. */
		public Resources(final NcrnaFamily family_, final RrnaPositionalKmerTable.Table[] five,
				final RrnaPositionalKmerTable.Table[] three, final CellNet fiveNet, final CellNet threeNet, final float margin_){
			this(family_, five, three, null, null, fiveNet, threeNet, margin_);
		}

		/** Explicit full-error route; both ends borrow three signed-background probes. */
		public Resources(final NcrnaFamily family_, final RrnaPositionalKmerTable.Table[] five,
				final RrnaPositionalKmerTable.Table[] three, final NcrnaWeightedEndpointTable[] fiveBackground,
				final NcrnaWeightedEndpointTable[] threeBackground, final CellNet fiveNet, final CellNet threeNet, final float margin_){
			this(family_, five, three, fiveBackground, threeBackground, fiveNet, threeNet, margin_, margin_, true);
		}

		/** Full-error head margins are finite and signed; a separate flag can ignore the head.
		 * Existing constructors retain one common positional/head confidence margin.
		 */
		public Resources(final NcrnaFamily family_, final RrnaPositionalKmerTable.Table[] five,
				final RrnaPositionalKmerTable.Table[] three, final NcrnaWeightedEndpointTable[] fiveBackground,
				final NcrnaWeightedEndpointTable[] threeBackground, final CellNet fiveNet, final CellNet threeNet,
				final float margin_, final float oorMargin_, final boolean useOorHead_){
			this(family_, five, three, fiveBackground, threeBackground, fiveNet, threeNet,
				margin_, oorMargin_, useOorHead_, null, null, 0, 0);
		}

		/** Optional second-stage gates replace margin selection on full-error proposals.
		 * Each gate consumes32raw inputs and emits one score; equality at its finite
		 * cutoff accepts. The last first-stage head remains an input, never a veto.
		 */
		public Resources(final NcrnaFamily family_, final RrnaPositionalKmerTable.Table[] five,
				final RrnaPositionalKmerTable.Table[] three, final NcrnaWeightedEndpointTable[] fiveBackground,
				final NcrnaWeightedEndpointTable[] threeBackground, final CellNet fiveNet, final CellNet threeNet,
				final float margin_, final float oorMargin_, final boolean useOorHead_, final CellNet fiveGate, final CellNet threeGate,
				final float fiveCutoff, final float threeCutoff){
			this(family_, five, three, fiveBackground, threeBackground, fiveNet, threeNet,
				margin_, oorMargin_, useOorHead_, fiveGate, threeGate, fiveCutoff, threeCutoff, null);
		}

		/** Optional model offsets precede both full-error feature vectors. Existing
		 * gates cannot consume centered scores without a separately trained resource.
		 */
		public Resources(final NcrnaFamily family_, final RrnaPositionalKmerTable.Table[] five,
				final RrnaPositionalKmerTable.Table[] three, final NcrnaWeightedEndpointTable[] fiveBackground,
				final NcrnaWeightedEndpointTable[] threeBackground, final CellNet fiveNet, final CellNet threeNet,
				final float margin_, final float oorMargin_, final boolean useOorHead_, final CellNet fiveGate, final CellNet threeGate,
				final float fiveCutoff, final float threeCutoff, final NcrnaEndpointModelOffsets modelOffsets_){
			this(family_, five, three, fiveBackground, threeBackground, fiveNet, threeNet, margin_, oorMargin_, useOorHead_,
				fiveGate, threeGate, fiveCutoff, threeCutoff, modelOffsets_, GatePolicy.THRESHOLD, GatePolicy.THRESHOLD);
		}

		/** Explicit whole-move controls retain finite numeric cutoffs in every mode. */
		public Resources(final NcrnaFamily family_, final RrnaPositionalKmerTable.Table[] five,
				final RrnaPositionalKmerTable.Table[] three, final NcrnaWeightedEndpointTable[] fiveBackground,
				final NcrnaWeightedEndpointTable[] threeBackground, final CellNet fiveNet, final CellNet threeNet,
				final float margin_, final float oorMargin_, final boolean useOorHead_, final CellNet fiveGate, final CellNet threeGate,
				final float fiveCutoff, final float threeCutoff, final NcrnaEndpointModelOffsets modelOffsets_,
				final GatePolicy fivePolicy, final GatePolicy threePolicy){
			if(family_==null || family_.library==null || family_.modelNames==null || family_.library.length==0
					|| family_.modelNames.length!=family_.library.length || !Float.isFinite(margin_) || margin_<0){
				throw new IllegalArgumentException("Positional resources bind a loaded family and finite nonnegative confidence margin");
			}
			// Validates the shared table geometry, including finalization and end order.
			if(five==null || three==null || five.length!=three.length || (five.length!=1 && five.length!=3)){
				throw new IllegalArgumentException("Positional resources require one or three ordered tables per end");
			}
			final boolean triple=five.length==3;
			final boolean full=fiveBackground!=null;
			final boolean gated=fiveGate!=null;
			centeredGate=gated && modelOffsets_!=null;
			if(fivePolicy==null || threePolicy==null
					|| !centeredGate && (fivePolicy!=GatePolicy.THRESHOLD || threePolicy!=GatePolicy.THRESHOLD)){
				throw new IllegalArgumentException("Explicit never/always policies belong only to separately trained centered35 gates");
			}
			if(gated!=(threeGate!=null) || (gated ? !full || margin_!=0 || useOorHead_ || oorMargin_!=0
					|| !Float.isFinite(fiveCutoff) || !Float.isFinite(threeCutoff)
					: fiveCutoff!=0 || threeCutoff!=0)){
				throw new IllegalArgumentException("Movement gates require both nets, declared per-end cutoffs and full-error zero-margin/ignored-head proposals");
			}
			if(!Float.isFinite(oorMargin_) || !full && (oorMargin_!=margin_ || !useOorHead_)){
				throw new IllegalArgumentException("Only full-error resources may select a separate finite signed head margin or ignore the head");
			}
			if(full!=(threeBackground!=null) || full && !triple){
				throw new IllegalArgumentException("Full-error inference requires three signed-background probes at both ends");
			}
			final NcrnaPositionalBoundaryFeatures geometry=full
				? new NcrnaPositionalBoundaryFeatures(five, three, fiveBackground, threeBackground, false)
				: triple ? new NcrnaPositionalBoundaryFeatures(five, three, false)
				: new NcrnaPositionalBoundaryFeatures(five[0], three[0]);
			if(gated && geometry.radius!=10){
				throw new IllegalArgumentException("Movement gates bind the radius-ten69/22first stage and their declared32/35-input layout");
			}
			if(modelOffsets_!=null && (!full || geometry.radius!=10 || !modelOffsets_.belongsTo(family_))){
				throw new IllegalArgumentException("Model offsets bind this exact radius-ten full-error family");
			}
			final NcrnaFamilyRegistry.Entry entry=NcrnaFamilyRegistry.get(family_.name);
			final String member=entry.group==null ? entry.id : entry.member;
			if(!member.equals(five[0].model) || !member.equals(three[0].model)){
				throw new IllegalArgumentException("Positional tables must bind the registered caller member: family="+entry.id+", member="+member);
			}
			RrnaEndpointCallerFeatures.Resources.validateNet(fiveNet, geometry.values[0].length, five[0].sites+(full ? 1 : 0));
			RrnaEndpointCallerFeatures.Resources.validateNet(threeNet, geometry.values[1].length, three[0].sites+(full ? 1 : 0));
			family=family_;
			tables=new RrnaPositionalKmerTable.Table[]{five[0], three[0]};
			probes=triple ? new RrnaPositionalKmerTable.Table[][]{five.clone(), three.clone()} : null;
			backgrounds=full ? new NcrnaWeightedEndpointTable[][]{fiveBackground.clone(), threeBackground.clone()} : null;
			plateaus=triple ? new float[family.library.length] : null;
			if(triple){
				if(family.models==null || family.models.length!=family.library.length){
					throw new IllegalArgumentException("Optional HBM features require one actual paired model per loaded consensus");
				}
				for(int model=0; model<plateaus.length; model++){
					if(family.models[model]==null || !Arrays.equals(family.models[model].original, family.library[model])){
						throw new IllegalArgumentException("Optional positional HBM and consensus disagree at model "+model);
					}
					plateaus[model]=TrnaBoundaryFeatures.modelMaxCoverage(family.models[model]);
					if(plateaus[model]<=0){throw new IllegalArgumentException("Optional positional HBM must have positive training coverage");}
				}
			}
			nets=new CellNet[]{fiveNet, threeNet};
			gateNets=gated ? new CellNet[]{fiveGate, threeGate} : null;
			gateCutoffs=gated ? new float[]{fiveCutoff, threeCutoff} : null;
			gatePolicies=gated ? new GatePolicy[]{fivePolicy, threePolicy} : null;
			modelGC=gated ? new float[family.library.length] : null;
			if(gated){
				RrnaEndpointCallerFeatures.Resources.validateNet(fiveGate, centeredGate ? 35 : 32, 1);
				RrnaEndpointCallerFeatures.Resources.validateNet(threeGate, centeredGate ? 35 : 32, 1);
				for(int model=0; model<modelGC.length; model++){
					modelGC[model]=NcrnaEndpointGateFeatures.gc(family.library[model], 0, family.library[model].length-1);
				}
			}
			margin=margin_;
			oorMargin=oorMargin_;
			useOorHead=useOorHead_;
			modelOffsets=modelOffsets_;
		}
		@Override
		public NcrnaBoundaryRefiner create(final NcrnaFamily family_){
			if(family_!=family){throw new IllegalArgumentException("Positional factory cannot be reused against another loaded family/library");}
			return new NcrnaPositionalBoundaryRefiner(this);
		}
		/** Immutable common radius, for the trace's positional output schema. */
		public int radius(){return tables[0].radius;}
		/** Actual primary table word length, for explicitly bound replay resources. */
		public int wordLength(){return tables[0].k;}
		/** Declared confidence margin; an ungated trace requires zero. */
		public float margin(){return margin;}
		/** Effective finite signed head margin; unused when the head is ignored. */
		public float oorMargin(){return oorMargin;}
		/** True when the extra first-stage output participates in margin selection. */
		public boolean usesOorHead(){return useOorHead;}
		/** Full-error output includes a non-positional head unsupported by trace v1. */
		public boolean fullErrors(){return backgrounds!=null;}
		/** True only when both movement gates replace the legacy margin selector. */
		public boolean gated(){return gateNets!=null;}
		/** Centered gates accept/reject the entire raw-to-reported movement. */
		public boolean centeredGate(){return centeredGate;}
		/** Immutable optional model axis, also available to fixed-vector parity checks. */
		public NcrnaEndpointModelOffsets modelOffsets(){return modelOffsets;}
		private final NcrnaFamily family;
		private final RrnaPositionalKmerTable.Table[] tables;
		private final RrnaPositionalKmerTable.Table[][] probes;
		private final NcrnaWeightedEndpointTable[][] backgrounds;
		private final float[] plateaus;
		private final CellNet[] nets;
		/** Templates and immutable per-end thresholds; both are absent without gates. Workers copy only the nets. */
		private final CellNet[] gateNets;
		private final float[] gateCutoffs, modelGC;
		private final GatePolicy[] gatePolicies;
		private final float margin;
		private final float oorMargin;
		private final boolean useOorHead;
		private final NcrnaEndpointModelOffsets modelOffsets;
		private final boolean centeredGate;
	}

	private NcrnaPositionalBoundaryRefiner(final Resources resources_){
		resources=resources_;
		features=resources.probes==null ? new NcrnaPositionalBoundaryFeatures(resources.tables[0], resources.tables[1])
			: resources.fullErrors() ? new NcrnaPositionalBoundaryFeatures(resources.probes[0], resources.probes[1],
				resources.backgrounds[0], resources.backgrounds[1], false)
			: new NcrnaPositionalBoundaryFeatures(resources.probes[0], resources.probes[1], false);
		nets=new CellNet[]{resources.nets[0].copy(false), resources.nets[1].copy(false)};
		gateNets=resources.gateNets==null ? null : new CellNet[]{resources.gateNets[0].copy(false), resources.gateNets[1].copy(false)};
		gateValues=gateNets==null ? null : new float[resources.centeredGate ? 35 : 32];
		gateScratch=resources.centeredGate ? new float[32] : null;
		inferenceScores=gateNets==null && resources.modelOffsets==null ? null : new float[22];
		center=resources.modelOffsets==null ? null : new int[2];
	}
	@Override
	public boolean applies(final int model){return model>=0 && model<resources.family.library.length;}
	@Override
	public boolean afterAlignmentTrim(){return true;}
	@Override
	public void setBoundaryTrace(final NcrnaBoundaryTrace.Family trace_){
		if(trace_!=null && resources.fullErrors()){
			throw new IllegalArgumentException("Boundary trace v1 cannot represent the full-error veto head or strict zero-margin policy");
		}
		if(trace_!=null && resources.margin!=0){throw new IllegalArgumentException("Trace the declared ungated positional inference");}
		trace=trace_;
	}

	/** Optional worker-confined diagnostic. Ordinary calls leave it null and allocate
	 * no diagnostic copies. The observer owns its copied raw features/scores and
	 * cannot change the network's buffers; callback exceptions propagate to caller.
	 */
	public void setScoreObserver(final ScoreObserver observer_){observer=observer_;}

	/** Receives actual inference inputs, every output and the selected offset index. */
	public interface ScoreObserver{
		void accept(int end, float[] input, float[] scores, int selected);
	}

	/** Receives copied raw gate inputs and its score/decision, before joint geometry holds. */
	public interface GateObserver{
		void accept(int end, float[] input, float score, int proposed, boolean accepted);
	}
	/** Optional worker-confined parity observer; ordinary calls leave it null. */
	public void setGateObserver(final GateObserver observer_){gateObserver=observer_;}
	/** Separate whole-coordinate observer keeps the legacy class-index contract intact. */
	public interface WholeMoveObserver{
		void accept(int end, float[] input, float score, long proposal, boolean accepted);
	}
	public void setWholeMoveObserver(final WholeMoveObserver observer_){wholeMoveObserver=observer_;}

	/** Coordinates are sense-oriented even on a reverse-strand genomic call. */
	@Override
	public boolean refine(final Orf orf, final byte[] bases, final int model, final float gc, final float identity,
			final boolean quantumIdentity, final int fixedMask, final int minLen, final int maxLen){
		if(orf==null || !applies(model) || !resources.family.name.equals(orf.ncrnaFamily)
				|| !resources.family.modelNames[model].equals(orf.trnaModel)){
			throw new IllegalArgumentException("Positional refinement requires the actual selected model of its bound family");
		}
		final int start=orf.start, stop=orf.stop;
		if(resources.modelOffsets!=null){
			if(!resources.modelOffsets.center(model, start, stop, fixedMask, bases.length, minLen, maxLen, center)){return false;}
			orf.start=center[0];
			orf.stop=center[1];
		}
		final int baseStart=orf.start, baseStop=orf.stop;
		gateSpanGC=-1;
		if(trace!=null){
			traceStart=start;
			traceStop=stop;
		}
		final int usable=resources.probes==null
			? features.fill(bases, baseStart, baseStop, resources.family.library[model], gc, fixedMask)
			: features.fill(bases, baseStart, baseStop, resources.family.library[model], resources.family.models[model], resources.plateaus[model], gc, fixedMask);
		if(resources.centeredGate){
			// Missing features, rejected decisions and invalid joint geometry all hold
			// the ordinary path anchor, never a partially applied median center.
			orf.start=start;
			orf.stop=stop;
			final long chosenStart=(usable&1)==0 ? start : predictWhole(0, bases, start, stop, baseStart, baseStop, model);
			final long chosenStop=(usable&2)==0 ? stop : predictWhole(1, bases, start, stop, baseStart, baseStop, model);
			if(chosenStart>=0 && chosenStop<bases.length && chosenStart<=Integer.MAX_VALUE && chosenStop>=0){
				RrnaEndpointCallerFeatures.apply(orf, (int)chosenStart, (int)chosenStop, bases.length, minLen, maxLen);
			}
			return orf.start!=start || orf.stop!=stop;
		}
		final int chosenStart=(usable&1)==0 ? baseStart : predict(0, baseStart, bases, baseStart, baseStop, model);
		final int chosenStop=(usable&2)==0 ? baseStop : predict(1, baseStop, bases, baseStart, baseStop, model);
		RrnaEndpointCallerFeatures.apply(orf, chosenStart, chosenStop, bases.length, minLen, maxLen);
		return orf.start!=start || orf.stop!=stop;
	}

	/** Both ends borrow the same center; no truth or selected neighbour enters inference. */
	private long predictWhole(final int end, final byte[] bases, final int rawStart, final int rawStop,
			final int centerStart, final int centerStop, final int model){
		assert(resources.centeredGate && gateScratch.length==32 && gateValues.length==35)
			: "Only separately trained centered35 resources may gate a whole raw-to-proposal move";
		final CellNet net=nets[end];
		net.applyInput(features.values[end]);
		net.feedForward();
		for(int i=0; i<22; i++){inferenceScores[i]=net.getOutput(i);}
		if(gateSpanGC<0){
			gateSpanGC=NcrnaEndpointGateFeatures.gc(bases, centerStart, centerStop);
		}
		final long proposal=NcrnaCenteredGateFeatures.fill(inferenceScores, features.values[end], gateSpanGC,
			resources.modelGC[model], resources.family.library[model].length,
			rawStart, rawStop, centerStart, centerStop, end, gateScratch, gateValues);
		final int raw=end==0 ? rawStart : rawStop;
		if(observer!=null){
			observer.accept(end, features.values[end].clone(), inferenceScores.clone(),
				NcrnaCenteredGateFeatures.adjustment(inferenceScores)+features.radius);
		}
		if(proposal==raw){return raw;}
		final CellNet gate=gateNets[end];
		gate.applyInput(gateValues);
		gate.feedForward();
		final float score=gate.getOutput(0);
		if(!Float.isFinite(score)){throw new IllegalArgumentException("Whole-move gates must emit finite scores even at never/always control endpoints");}
		final GatePolicy policy=resources.gatePolicies[end];
		final boolean accepted=policy==GatePolicy.ALWAYS || policy==GatePolicy.THRESHOLD && score>=resources.gateCutoffs[end];
		if(wholeMoveObserver!=null){wholeMoveObserver.accept(end, gateValues.clone(), score, proposal, accepted);}
		return accepted ? proposal : raw;
	}

	/** One forward evaluation yields every position; confidence never rescans the network. */
	private int predict(final int end, final int raw, final byte[] bases, final int start, final int stop, final int model){
		final CellNet net=nets[end];
		final float[] observedInput=observer==null ? null : features.values[end].clone();
		net.applyInput(features.values[end]);
		net.feedForward();
		final float[] scores=inferenceScores==null ? net.getOutput() : inferenceScores;
		if(inferenceScores!=null){for(int i=0; i<22; i++){scores[i]=net.getOutput(i);}}
		if(trace!=null){
			trace.scores(end, NcrnaBoundaryTrace.INDEPENDENT, traceStart, traceStop, scores.clone());
		}
		final int selected;
		if(gateNets==null){
			selected=resources.fullErrors() ? NcrnaEndpointConfidence.selectFull(scores, features.radius,
				resources.margin, resources.oorMargin, resources.useOorHead)
				: NcrnaEndpointConfidence.select(scores, features.radius, RrnaResourceIO.argmax(scores), resources.margin);
		}else{
			final int proposed=NcrnaEndpointGateFeatures.proposal(scores);
			if(proposed<0){selected=features.radius;}
			else{
				if(gateSpanGC<0){
					gateSpanGC=NcrnaEndpointGateFeatures.gc(bases, start, stop);
				}
				final long length=(long)stop-start+1+(end==0 ? 10-proposed : proposed-10);
				NcrnaEndpointGateFeatures.fill(scores, features.values[end], gateSpanGC,
					resources.modelGC[model], length, resources.family.library[model].length, gateValues);
				final float[] observedGate=gateObserver==null ? null : gateValues.clone();
				final CellNet gate=gateNets[end];
				gate.applyInput(gateValues);
				gate.feedForward();
				final float score=gate.getOutput(0);
				if(!Float.isFinite(score)){throw new IllegalArgumentException("Movement gate must emit a finite score");}
				final boolean accepted=score>=resources.gateCutoffs[end];
				selected=accepted ? proposed : features.radius;
				if(gateObserver!=null){gateObserver.accept(end, observedGate, score, proposed, accepted);}
			}
		}
		if(observer!=null){observer.accept(end, observedInput, scores.clone(), selected);}
		return RrnaEndpointVector.endpointFromClass(raw, selected, features.radius);
	}

	private final Resources resources;
	private final NcrnaPositionalBoundaryFeatures features;
	private final CellNet[] nets;
	private final CellNet[] gateNets;
	/** Reused gate inputs and full-error scores; centered inference also reuses scores. */
	private final float[] gateValues, gateScratch, inferenceScores;
	/** Present only for centered resources; no per-call coordinate allocation. */
	private final int[] center;
	/** GC is in [0,1]; -1 is a finite, noncolliding unmeasured cache value. */
	private float gateSpanGC=-1;
	private GateObserver gateObserver;
	private WholeMoveObserver wholeMoveObserver;
	private NcrnaBoundaryTrace.Family trace;
	private ScoreObserver observer;
	private int traceStart, traceStop;
}
