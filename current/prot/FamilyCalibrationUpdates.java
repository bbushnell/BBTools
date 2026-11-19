package prot;

/** Worker-local acceptance maxima against one immutable classifier parent.
 * Merge only after workers finish; publish a new roster at the round barrier.
 * The borrowed parent array and its records must remain immutable.
 * @author Brian Bushnell, Keqing
 */
public final class FamilyCalibrationUpdates{

	public FamilyCalibrationUpdates(FamilyReference[] parent,int reject,int accept,float idFloor,float blFloor,FamilyAlignmentMetric metric_){
		require(parent!=null && reject>=0 && reject<accept && accept<=1025,
			"Calibration needs a parent roster and 0<=DDL reject<accept<=1025");
		require(Float.isFinite(idFloor) && idFloor>=30 && idFloor<=100 && Float.isFinite(blFloor) && blFloor>=0,
			"Calibration floors must match finite deployed identity/BLOSUM gates");
		require(metric_!=null,"Calibration updates must retain the scoring metric");
		families=parent;ddlReject=reject;ddlAccept=accept;identityFloor=idFloor;blosumFloor=blFloor;
		metric=metric_;
		identity=new double[parent.length];blosum=new double[parent.length];
		for(int i=0;i<parent.length;i++){
			require(parent[i]!=null,"A calibration parent cannot have missing roster entries");
			identity[i]=parent[i].identityUpper;blosum[i]=parent[i].blosumFinal;
		}
	}

	/** Historical callers calibrate the legacy glocal full-query metric. */
	public FamilyCalibrationUpdates(FamilyReference[] parent,int reject,int accept,float idFloor,float blFloor){
		this(parent,reject,accept,idFloor,blFloor,FamilyAlignmentMetric.LEGACY);
	}

	/** Stages 1/2/3 are the existing DDL/identity/BLOSUM terminal accepts.
	 * Scores are the original Java floats. Apply their historical double increments
	 * before any upward float rounding, exactly as FamilyPanelEscalation does.
	 */
	public void observe(int ordinal,int stage,int ddl,float identityScore,float blosumScore){
		require(metric.same(FamilyAlignmentMetric.LEGACY),
			"Raw float calibration observations are only valid for legacy alignment scoring");
		observeState(ordinal,stage,ddl,identityScore,blosumScore,FamilyAlignmentMetric.State.MEASURED);
	}

	/** Observation-aware path for metric-bound packages; failed alignments cannot enter maxima. */
	public void observe(int ordinal,int stage,int ddl,FamilyAlignmentMetric.Observation observation){
		require(observation!=null,"Metric-bound calibration needs an explicit alignment observation");
		observeState(ordinal,stage,ddl,observation.measured() ? (float)observation.identityPct() : 0f,
			observation.measured() ? (float)observation.blosum() : 0f,observation.state);
	}

	/** Measured callback path for callers that receive scalar scores plus explicit state from a classifier. */
	public void observe(int ordinal,int stage,int ddl,float identityScore,float blosumScore,
			FamilyAlignmentMetric.State alignmentState){
		observeState(ordinal,stage,ddl,identityScore,blosumScore,alignmentState);
	}

	private void observeState(int ordinal,int stage,int ddl,float identityScore,float blosumScore,
			FamilyAlignmentMetric.State alignmentState){
		observeAtMargin(ordinal,stage,ddl,identityScore,blosumScore,alignmentState,1);
	}

	/** Scores come from this parent's explicitly scaled evaluation view. Proposals
	 * are converted back into the parent's stored upper-gate coordinate system,
	 * so the existing immutable publication path retains all other policy fields. */
	public void observeAtMargin(int ordinal,int stage,int ddl,float identityScore,float blosumScore,
			FamilyAlignmentMetric.State alignmentState,double marginScale){
		require(Double.isFinite(marginScale) && marginScale>0,
			"Calibration observations require their positive finite evaluation margin");
		require(ordinal>=0 && ordinal<families.length && stage>=1 && stage<=3,
			"Only accepted winners from this parent may propose gate updates");
		require(alignmentState!=null,"Aligned calibration observations need explicit state");
		final FamilyReference family=families[ordinal];
		require(family.members>0 && ddl>=0 && ddl<=1024,"Calibration winner must be live with a valid DDL count");
		if(stage==1){
			require(ddl>=ddlAccept,"Shortcut winner must reach the parent DDL acceptance threshold");
			return;
		}
		require(family.alignmentCalibration==FamilyReference.AlignmentCalibration.ENABLED && alignmentState==FamilyAlignmentMetric.State.MEASURED,
			"Aligned calibration winners require measured scores for an alignment-enabled family");
		require(ddl>=ddlReject && ddl<ddlAccept && Float.isFinite(identityScore) && identityScore>=0 &&
			identityScore<=100 && Float.isFinite(blosumScore),"Invalid aligned calibration winner scores");
		final double idGate=marginScale==1 ? Math.max(identityFloor,family.identityUpper) :
			FamilyTailCascade.evaluationGate(family.identityLower,family.identityUpper,identityFloor,marginScale);
		final double blGate=marginScale==1 ? Math.max(blosumFloor,family.blosumFinal) :
			FamilyTailCascade.evaluationGate(family.blosumLower,family.blosumFinal,blosumFloor,marginScale);
		if(stage==2){
			require(identityScore>idGate,"Identity winner must exceed the immutable parent gate");
			identity[ordinal]=Math.max(identity[ordinal],marginScale==1 ? (double)identityScore+1.0 :
				FamilyCalibrationMargin.raiseUpper(family.identityLower,family.identityUpper,identityFloor,marginScale,identityScore,1.0));
		}else{
			require(identityScore<=idGate && identityScore>=family.identityLower && blosumScore>blGate,
				"BLOSUM winner must follow the immutable parent identity/BLOSUM acceptance path");
			blosum[ordinal]=Math.max(blosum[ordinal],marginScale==1 ? (double)blosumScore+.005 :
				FamilyCalibrationMargin.raiseUpper(family.blosumLower,family.blosumFinal,blosumFloor,marginScale,blosumScore,.005));
		}
	}

	/** Creates zero-observation worker state bound to the same immutable parent. */
	public FamilyCalibrationUpdates newAccumulator(){
		return new FamilyCalibrationUpdates(families,ddlReject,ddlAccept,identityFloor,blosumFloor,metric);
	}

	/** Package-internal identity check binds maxima to both roster and deployed gates. */
	boolean boundTo(FamilyTailCascade.SearchConfig config){
		return config!=null && config.margins==null && families==config.families && ddlReject==config.ddlReject && ddlAccept==config.ddlAccept &&
			identityFloor==config.identityFloor && blosumFloor==config.blosumFloor && metric.same(config.alignmentMetric);
	}

	/** Maximum reduction is independent of worker and panel order. Parents may not be mixed. */
	public void add(FamilyCalibrationUpdates other){
		require(other!=null && other.families==families && other.ddlReject==ddlReject && other.ddlAccept==ddlAccept &&
			other.identityFloor==identityFloor && other.blosumFloor==blosumFloor && other.metric.same(metric),
			"Merge calibration maxima only from the same immutable parent and deployed policy");
		for(int i=0;i<families.length;i++){
			identity[i]=Math.max(identity[i],other.identity[i]);blosum[i]=Math.max(blosum[i],other.blosum[i]);
		}
	}

	public double identityAt(int ordinal){return identity[ordinal];}
	public double blosumAt(int ordinal){return blosum[ordinal];}
	public boolean changed(int ordinal){
		return identity[ordinal]>families[ordinal].identityUpper || blosum[ordinal]>families[ordinal].blosumFinal;
	}

	/** Upward rounding happens once when publishing. Unchanged immutable records are shared. */
	public FamilyReference[] updated(){
		final FamilyReference[] result=families.clone();
		for(int i=0;i<result.length;i++){
			if(!changed(i)){continue;}
			final FamilyReference f=families[i];
			final float id=FamilyGateCodec.roundUpper(identity[i]),bl=FamilyGateCodec.roundUpper(blosum[i]);
			assert(f.members>0 && id>=f.identityUpper && bl>=f.blosumFinal) :
				"Winner maxima may only raise live-family acceptance gates; lower gates and payloads remain fixed";
			result[i]=f.explicitAlignmentCalibration ?
				new FamilyReference(f.id,f.rank,f.members,f.residues,f.validWords,f.consensus,f.ddl,f.spectrum,
					f.cosineLower,f.identityLower,f.blosumLower,f.cosineAccept,id,bl,
					f.consensusRaw,f.scalarMask,f.scalarBounds,f.alignmentCalibration) :
				new FamilyReference(f.id,f.rank,f.members,f.residues,f.validWords,f.consensus,f.ddl,f.spectrum,
					f.cosineLower,f.identityLower,f.blosumLower,f.cosineAccept,id,bl,
					f.consensusRaw,f.scalarMask,f.scalarBounds);
		}
		return result;
	}

	private static void require(boolean valid,String why){if(!valid){throw new IllegalArgumentException(why);}}
	private final FamilyReference[] families;
	private final double[] identity,blosum;
	private final FamilyAlignmentMetric metric;
	private final int ddlReject,ddlAccept;
	private final float identityFloor,blosumFloor;
}
