package prot.cluster;

import java.util.Arrays;

import structures.ByteBuilder;

/** Numeric evidence captured before residual-pool scoring, with explicit missingness.
 * Parent/left/right describe the same three models in either direction: a split
 * replaces parent with left/right; a fusion replaces left/right with parent.
 * No model, sequence, membership array or measured outcome is retained here.
 * @author Brian Bushnell, Fischl
 */
public final class WideProbeFeatures{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** A context names a committed run manifest containing source, inputs, seeds and
	 * configuration. It is audit metadata, never a numeric predictor.
	 */
	public static final class Output{
		public Output(final String path_, final String context_, final int layer_){
			if(path_==null || path_.isEmpty() || context_==null || context_.isEmpty() || layer_<0){
				throw new IllegalArgumentException("Wide capture requires an output, manifest context token and nonnegative layer");
			}
			for(int i=0; i<context_.length(); i++){
				if(context_.charAt(i)<=32 || context_.charAt(i)>=127){
					throw new IllegalArgumentException("Feature context must be a printable ASCII token referencing run provenance");
				}
			}
			path=path_;context=context_;layer=layer_;
		}
		public final String path,context;
		public final int layer;
	}

	/** Snapshot model-local predictors while the request still owns its candidates. */
	public WideProbeFeatures(final Family parent, final Family left, final Family right,
			final int population, final int modelCount, final double slot, final int layer, final boolean fission){
		if(population<0 || modelCount<0 || !Double.isFinite(slot) || slot<0 || layer<0){
			throw new IllegalArgumentException("Feature census and scheduled layer must be nonnegative");
		}
		values[0]=population;values[1]=Math.log1p(population);values[2]=modelCount;
		values[3]=slot;values[4]=layer;values[5]=fission ? 1 : 0;
		captureModel(0,parent);captureModel(1,left);captureModel(2,right);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Extraction          ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses current construction memberships, never the codec's source reference count. */
	private void captureModel(final int model, final Family family){
		if(family==null){return;}
		final int at=modelOffset(model);
		final byte[] sequence=family.consensusUnsafe();
		values[at]=1;values[at+1]=family.size();values[at+2]=Math.log1p(family.size());values[at+3]=sequence.length;
		final int[] counts=new int[20];
		final boolean[] dimers=new boolean[400];
		int valid=0,pairs=0,distinct=0;
		for(int i=0; i<sequence.length; i++){
			final int residue=sequence[i];
			if(residue<0 || residue>=20){continue;}
			counts[residue]++;valid++;
			if(i>0 && sequence[i-1]>=0 && sequence[i-1]<20){
				final int pair=sequence[i-1]*20+residue;
				pairs++;if(!dimers[pair]){dimers[pair]=true;distinct++;}
			}
		}
		double entropy=0;
		if(valid>0){
			for(int count:counts){
				if(count>0){final double p=count/(double)valid;entropy-=p*Math.log(p)/LOG2;}
			}
		}
		values[at+4]=entropy;
		values[at+5]=pairs==0 ? 0 : distinct/(double)pairs;
		values[at+6]=(sequence.length-valid)/(double)sequence.length;
		final char[] sketch=family.sketchUnsafe();
		if(sketch!=null){
			if(sketch.length==0){throw new IllegalArgumentException("Present model sketches need positive width");}
			values[at+7]=1;values[at+8]=family.occupiedBins();values[at+9]=family.occupiedBins()/(double)sketch.length;
		}
		final double[] scalars=family.scalarsUnsafe();
		for(int i=0; i<5; i++){
			final boolean available=(family.scalarMask()&(1<<i))!=0;
			values[at+10+2*i]=available ? 1 : 0;
			values[at+11+2*i]=available ? scalars[i] : 0;
		}
		assert(valid<=sequence.length && distinct<=pairs) : "Composition features count only actual standard residues and adjacent pairs";
	}

	/** Stores each independent cohesion snapshot, even if another model is unavailable. */
	public void quality(final int model, final ConsensusQualityReport.Snapshot quality){
		final int at=modelOffset(model)+20;
		if(quality==null || !quality.available()){return;}
		values[at]=1;values[at+1]=quality.identityMean();values[at+2]=quality.scoreRatioMean();
		values[at+3]=.5*(quality.identityMean()+Math.max(0,Math.min(1,quality.scoreRatioMean())));
	}

	/** Keeps nearest quality and the explicit directed sibling/input pair separately. */
	public void neighbor(final int model, final FamilyNeighborEvidence.Result neighbor){
		final int at=modelOffset(model)+24;
		if(neighbor==null){return;}
		if(neighbor.available){
			values[at]=1;values[at+1]=neighbor.quality;values[at+2]=neighbor.identity;
			values[at+3]=neighbor.scoreRatio;values[at+4]=neighbor.sharedFraction;
		}
		if(model>0 && neighbor.explicitAvailable){
			final int pair=PAIR_START+(model-1)*4;
			values[pair]=1;values[pair+1]=neighbor.explicitIdentity;
			values[pair+2]=neighbor.explicitScoreRatio;values[pair+3]=neighbor.explicitSharedFraction;
		}
	}

	/** Original three-feature fit remains optional; masked raw blocks survive its absence. */
	public void deltas(final LayerCakeClustering.ProbeFeatures features){
		if(features==null){return;}
		final int at=DELTA_START;
		values[at]=1;values[at+1]=features.cohesionDelta;values[at+2]=features.separationDelta;
		values[at+3]=features.sizePressure;values[at+4]=features.balance;
		values[at+5]=features.cohesionDelta/features.parent.members;
		values[at+6]=features.separationDelta/features.parent.members;
		values[at+7]=Math.max(features.left.members,features.right.members)/(double)Math.min(features.left.members,features.right.members);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Stable TSV schema      ----------------*/
	/*--------------------------------------------------------------*/

	/** Appends only predictors. Outcomes, IDs, statuses and timings belong outside this range. */
	public void append(final ByteBuilder line){
		if(line==null){throw new IllegalArgumentException("Feature output needs a byte builder");}
		for(double value:values){
			if(!Double.isFinite(value)){throw new IllegalStateException("Nonfinite pre-probe feature cannot enter training data");}
			line.tab().appendSlow(value);
		}
	}

	/** Defensive numeric view for model consumers and independent fixtures. */
	public double[] values(){return values.clone();}
	public static String[] names(){return NAMES.clone();}

	/** Header ordering is generated from the same fixed block sizes as extraction. */
	private static String[] makeNames(){
		final String[] names=new String[DIMENSIONS];
		final String[] global={"input_tails","log_input_tails","model_count","slot_cost","layer","is_fission"};
		System.arraycopy(global,0,names,0,global.length);
		final String[] model={"present","members","log_members","length","entropy_bits","distinct_dimer_fraction","unknown_fraction",
			"sketch_available","occupied_bins","occupied_fraction","hydrophobicity_available","hydrophobicity",
			"charge_available","charge","aromaticity_available","aromaticity","small_fraction_available","small_fraction",
			"pi_available","pi","quality_available","identity_mean","score_ratio_mean","cohesion",
			"neighbor_available","neighbor_quality","neighbor_identity","neighbor_score_ratio","neighbor_shared_fraction"};
		final String[] prefixes={"parent_","left_","right_"};
		for(int side=0; side<3; side++){
			for(int i=0; i<MODEL_WIDTH; i++){names[modelOffset(side)+i]=prefixes[side]+model[i];}
		}
		final String[] pair={"available","identity","score_ratio","shared_fraction"};
		for(int side=0; side<2; side++){
			for(int i=0; i<4; i++){names[PAIR_START+side*4+i]=(side==0 ? "left_to_right_" : "right_to_left_")+pair[i];}
		}
		final String[] delta={"deltas_available","weighted_cohesion_delta","weighted_separation_delta","size_pressure","balance",
			"per_member_cohesion_delta","per_member_separation_delta","pair_size_ratio"};
		System.arraycopy(delta,0,names,DELTA_START,delta.length);
		assert(!Arrays.asList(names).contains(null)) : "Every stored numeric predictor needs a stable versioned column name";
		return names;
	}

	private static int modelOffset(final int model){
		if(model<0 || model>=3){throw new IllegalArgumentException("Feature model must be parent, left or right");}
		return GLOBAL_WIDTH+model*MODEL_WIDTH;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	public static final String SCHEMA="wide_probe_v1";
	private static final int GLOBAL_WIDTH=6,MODEL_WIDTH=29,PAIR_START=GLOBAL_WIDTH+3*MODEL_WIDTH,DELTA_START=PAIR_START+8;
	public static final int DIMENSIONS=DELTA_START+8;
	private static final double LOG2=Math.log(2);
	private static final String[] NAMES=makeNames();
	private final double[] values=new double[DIMENSIONS];
}
