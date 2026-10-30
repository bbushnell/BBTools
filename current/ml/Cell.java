package ml;

import java.util.Arrays;
import shared.Random;

import shared.Tools;
import shared.KillSwitch;
import simd.Vector;

/**
 * One neuron row in a {@link CellNet}: a shared activation slot, one activation
 * function, a bias, and either dense incoming weights or sparse incoming
 * weights plus {@link #inputs} indices.
 *
 * @author Brian Bushnell
 * @date 2025
 */
public class Cell extends Source{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs a dense row wrapper around the owning layer's shared value buffers. */
	public Cell(int id_, int activationType, int lpos_, int layer_, int maxLayer_,
			int prevLayerStart_, int nextLayerStart_, int wid,
			float[] values_, float[] eOverNetArray_){
		this(id_, activationType, lpos_, layer_, maxLayer_, prevLayerStart_, nextLayerStart_,
			wid, values_, eOverNetArray_, true);
	}

	/**
	 * Creates a row bound to the owning layer's shared activation and gradient arrays.
	 * @param id_ global cell id
	 * @param activationType activation function id
	 * @param lpos_ position within this layer's shared arrays
	 * @param layer_ layer index for this row
	 * @param maxLayer_ output layer index
	 * @param prevLayerStart_ global id of the previous layer's first cell
	 * @param nextLayerStart_ global id of the next layer's first cell
	 * @param wid width of the owning layer
	 * @param values_ borrowed activation array for the owning layer
	 * @param eOverNetArray_ borrowed gradient array for the owning layer
	 * @param dense_ immutable dense versus sparse incoming-edge representation
	 */
	public Cell(int id_, int activationType, int lpos_, int layer_, int maxLayer_,
			int prevLayerStart_, int nextLayerStart_, int wid,
			float[] values_, float[] eOverNetArray_, boolean dense_){
		dense=dense_;
		id=id_;
//		type=activationType;
		function=Function.getFunction(activationType);
		lpos=lpos_;
		layer=layer_;
		maxLayer=maxLayer_;
		values=values_;
		eOverNetArray=eOverNetArray_;
		assert(values.length==wid);
		prevLayerStart=prevLayerStart_;
		nextLayerStart=nextLayerStart_;

		//Initialize later
//		inputs=in;
//		outputs=out;
//		weights=(inputs==null || inputs.length<1 ? null : new float[inputs.length]);
//		deltas=(inputs==null || inputs.length<1 ? null : new float[inputs.length]);

		assert(lpos>=0 & lpos<wid);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Non-Mutators         ----------------*/
	/*--------------------------------------------------------------*/

	/** Dot-products all previous-layer values with dense incoming weights and stores activation. */
	public void summateDense(float[] valuesIn){
		sum=bias;
		assert(valuesIn.length==weights.length) : valuesIn.length+", "+weights.length;
		sum+=Vector.fma(weights, valuesIn);
		final float v=(float)activation(sum);
		setValue(v);
	}

	/** Dot-products sorted-unique sparse inputs; full-width rows are identity maps and use dense FMA. */
	public void summateSparse(float[] valuesIn, int edgeBlockSize){
		sum=bias;
		sum+=(weights.length==valuesIn.length ? Vector.fma(weights, valuesIn) :
			Vector.fma(weights, valuesIn, inputs, edgeBlockSize, Vector.SIMD_FMA_SPARSE));
		final float v=(float)activation(sum);
		setValue(v);
	}

	/** Returns half squared error for the current output versus one target. */
	public float calcError(float ideal){
		float e=ideal-value();
		return 0.5f*e*e;
	}

	/** Assertion-backed structural sanity check for shared buffers and row arrays. */
	synchronized public boolean check(){
//		assert(false);

		if(!dense){
			assert(layer==maxLayer || outputs!=null) :
				layer+", "+maxLayer+", "+lpos+", "+id();
		}

		if(value!=values[lpos] && (value!=-1 && values[lpos]!=0)){
			assert(false) : id+", "+layer+", "+lpos+": "+value+", "+values[lpos]+", "+Arrays.toString(values);
			return false;
		}
		if(eOverNet!=eOverNetArray[lpos]){
			assert(false) : id+", "+layer+", "+lpos+": "+eOverNet+", "+eOverNetArray[lpos]+", "+Arrays.toString(values);
			return false;
		}
		if(layer>0){
			assert(weights!=null);
			assert(dense==(inputs==null));
//			assert(weights2.length==weights.length);
			assert(deltas==null || deltas.length==weights.length);
//			for(int i=0; i<weights2.length; i++) {
//				Edge e=outputs.get(i);
//				if(e.weight2!=weights2[i]) {
//					assert(false) : id+", "+layer+", "+lpos+": "+e.weight2+", "+weights2[i];
//					return false;
//				}
//			}
		}else{
//			assert(inputs==null || outputs.isEmpty())  : id+", "+layer+", "+lpos+": "+weights2+": "+outputs+"\n"+toString()+"\n";
			assert(inputs==null || outputs==null)  : id+", "+layer+", "+lpos+": "+deltas+": "+outputs+"\n"+toString()+"\n";
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Mutators           ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Applies accumulated bias/edge deltas; input rows are unchanged, and deltas
	 * are cleared while retaining old weights for dense zero slots and candidate
	 * magnitudes at or below {@link Float#MIN_NORMAL}.
	 */
	void applyUpdates(float invSamples, float alpha){//only called once per epoch
		if(layer<1){return;}
		adjustBias(invSamples, alpha);
		for(int i=0; i<weights.length; i++){
			final float w=weights[i];

			final float d3=(float)(deltas[i]*invSamples*alpha);
			float w4=w+d3;
			{
				final float absW4=Math.abs(w4);
				if(absW4>Math.abs(w) && absW4>edgeAmplitudeIncreaseThresh){
					//Very sensitive; strong values break convergence.
					//0.98 seems OK though.  Goal is to keep edges low.
					//.94 is too extreme.  .96 seems OK.
					w4=w+edgeAmplitudeIncreaseMult*d3;//Slow down increase in edge magnitude
				}
			}
			deltas[i]=0;

			if((w!=0 || !dense) && Math.abs(w4)>Float.MIN_NORMAL){weights[i]=w4;}
		}
	}

	/** Accumulates one loss value into this cell's epoch error. */
	public void addError(float e){
		assert(error>=0);
		error+=e;
	}

	/** Clears the accumulated epoch error. */
	public void clearError(){
		error=0;
	}

	/** Clears accumulated bias, error and, for non-input layers, weight deltas. */
	public void clearTemp(){
//		bias2=0;
		biasDelta=0;
		error=0;
		if(layer<1){return;}
//		Arrays.fill(weights2, 0);
		assert(deltas!=null) : id()+", "+layer;
		Arrays.fill(deltas, 0);
	}

	/** Sets the bias; ignoreAssertion is retained for disabled legacy difference checks. */
	public void setBias(float b2, boolean ignoreAssertion){
		assert(layer>0);
//		float dif=Tools.absdif(bias, b2);
//		assert(ignoreAssertion || dif<0.1f || dif<0.2f*bias) : dif+", "+bias+", "+b2; //just checking for bugs.
//		if(!(ignoreAssertion || dif<0.1f || dif<0.2f*bias)) {System.out.print("*");}
//		if(!ignoreAssertion) {System.err.println(bias);}
		bias=b2;
	}

	/** Applies and clears the accumulated bias delta using epoch and learning-rate scaling. */
	private void adjustBias(float invSamples, float alpha){
		assert(layer>0);
		assert(alpha>0) : alpha;
//		System.err.println(String.format("b"+id()+": %.5f -> %.5f * %.2f = %.5f",
//				bias, bias2, invSamples, bias2*invSamples));
//		float b=((float)bias2)*invSamples;
		float bFromDelta=(float)(bias+biasDelta*invSamples*alpha*biasAlphaMult);
//		assert(Tools.absdif(b, bFromDelta)<0.00001) : "b="+b+", bFD="+bFromDelta+", bias="+bias+", bias2="+bias2+", bD="+biasDelta+", invS="+invSamples+", a="+alpha+", bAM="+biasAlphaMult
//			+", layer="+layer+" cid="+id;
		setBias(bFromDelta, false);
//		bias2=0;
		biasDelta=0;
	}

	/** Adds another training copy's accumulated epoch error into this cell. */
	public void addError(Cell c2){
		error+=c2.error;
	}

	@Override
	public void setValue(float v){
		assert(layer==0 || v==(float)activation(sum)) : v+", "+activation(sum)+", "+sum;
//		assert(value==values[lpos]);
		values[lpos]=value=v;
	}

	/*--------------------------------------------------------------*/

	/** Applies this row's activation function to a pre-activation sum. */
	public final double activation(double x){
		return function.activate(x);
	}

	/** Returns this activation's derivative from precomputed x and f(x). */
	public final double derivativeXFX(double x, double fx){
		double d=function.derivativeXFX(x, fx);
		assert(!Double.isNaN(d)) : x+", "+fx+", "+d;
		return d;
	}

	/*--------------------------------------------------------------*/

	//TODO: Vectorize?  Final layer is usually small though
	/** Accumulates final-layer gradients for rows whose source index equals weight position. */
	public void updateEdgesFinalLayerDense(float target, float[] valuesIn, float weightMult){
		//assert(check());
		final float v=value();
		eTotalOverOut=calcETotalOverOut(v, target, weightMult);
		outOverNet=(float)derivativeXFX(sum, v);
//		final double eTotalOverOut_X_outOverNet=eTotalOverOut*outOverNet;
		eOverNet=eOverNetArray[lpos]=eTotalOverOut*outOverNet;

		for(int i=0; i<weights.length; i++){
//			assert(e.source.check());
			final float netOverWeight=valuesIn[i];//e.source.value();
//			assert(valuesIn[i]==e.source.value) : valuesIn[i]+", "+e.source.value+", "+valuesIn.length+", "+inputs.size();
			final float eTotalOverWeight=eOverNet*netOverWeight;

//			final double w=weights[i];
//			final double incr=w-alpha*eTotalOverWeight;
//			weights2[i]+=incr;

			deltas[i]-=eTotalOverWeight;
			//TODO: Alpha is not really needed here;
			//the result could be accumulated and multiplied by alpha at the end
		}
		{
			assert(layer>0);
//			bias2+=bias-alpha*eOverNet*biasAlphaMult;
			biasDelta-=eOverNet;
		}
		//assert(check());
	}

	/*--------------------------------------------------------------*/

	/** Classifies output connectivity once, during CellNet.makeWeightMatrices(). */
	void prepareOutputGradient(final int inputWidth){
		sparseOutputGradient=false;
		if(inputs==null){return;}
		assert(inputs.length==weights.length) : KillSwitch.assertDie(
			"Sparse output-gradient indexing requires one source per weight; cell="+id);
		for(int i=0; i<inputs.length; i++){
			final int source=inputs[i];
			assert(source>=0 && source<inputWidth && (i==0 || inputs[i-1]<source)) : KillSwitch.assertDie(
				"Sparse output-gradient indexing requires ordered unique in-range sources; positional indexing "+
				"previously trained the wrong edges; cell="+id+", edge="+i+", source="+source);
			sparseOutputGradient|=(source!=i);
		}
	}

	/** Accumulates final-layer gradients for sparse rows with non-positional input indices. */
	void updateEdgesFinalLayerSparse(final float target, final float[] valuesIn, final float weightMult){
		assert(inputs!=null && inputs.length==weights.length) : KillSwitch.assertDie(
			"Cached sparse output-gradient mode must match the topology prepared with the training matrices");
		final float v=value();
		eTotalOverOut=calcETotalOverOut(v, target, weightMult);
		outOverNet=(float)derivativeXFX(sum, v);
		eOverNet=eOverNetArray[lpos]=eTotalOverOut*outOverNet;
		for(int i=0; i<weights.length; i++){
			final int source=inputs[i];
			assert(source>=0 && source<valuesIn.length) : KillSwitch.assertDie(
				"Sparse output-gradient fix indexes source activations, not weight positions; source="+source);
			final float netOverWeight=valuesIn[source];
			final float eTotalOverWeight=eOverNet*netOverWeight;
			deltas[i]-=eTotalOverWeight;
		}
		biasDelta-=eOverNet;
	}

	/** Returns the connectivity classification prepared with the training matrices. */
	boolean needsIndexedOutputGradient(){return sparseOutputGradient;}

	/** Backpropagates through a hidden dense row, using precomputed SPECIAL_FMA sums when enabled. */
	public void updateEdgesHiddenLayerDense(float[] valuesIn, float[] eOverNetNext, float[] weightsOut){
		//assert(check());
//		eTotalOverOut=0;
		final float v=value();
		assert(v==values[lpos]) : v+", "+values[lpos]+//Also fires on NaN, but that shouldn't happen...
			"\n"+Arrays.toString(values)+
			"\n"+Arrays.toString(weights)+"\n";
		outOverNet=(float)derivativeXFX(sum, v);

		assert(dense);
		if(!CellNet.SPECIAL_FMA){
			eTotalOverOut=Vector.fma(weightsOut, eOverNetNext);
		}

//		if(CellNet.SIMD && eOverNetNext.length>=16) {
//			eTotalOverOut=shared.Vector.fma(weightsOut, eOverNetNext);
//		}else {
//			for(int dest=0; dest<eOverNetNext.length; dest++){
//				final float netOverOut=weightsOut[dest];
//				final float eOverNetDest=eOverNetNext[dest];
//				final float eOverOut=eOverNetDest*netOverOut;
//
//				eTotalOverOut+=eOverOut;
//			}
//		}

		eOverNetArray[lpos]=eOverNet=eTotalOverOut*outOverNet;

//		for(int source=0; source<valuesPrev.length; source++) {
//			final float netOverWeight=valuesPrev[source];
//			final float eTotalOverWeight=eOverNet*netOverWeight;
//
////			final double incr=weights[source]-alpha*eTotalOverWeight;
////			weights2[source]+=incr;
//
//			deltas[source]-=eTotalOverWeight;
//		}
		Vector.addProduct(deltas, valuesIn, -eOverNet);

//		if(layer>0){
		{
			assert(layer>0) : layer;
//			bias2+=bias-alpha*eOverNet*biasAlphaMult;//Bias adjusts slower than edges
			biasDelta-=eOverNet;
		}
		//assert(check());
	}

	/** Backpropagates sorted-unique sparse rows; full-width rows are identity maps using dense overloads. */
	public void updateEdgesHiddenLayerSparse(float[] valuesIn, float[] eOverNetNext,
			float[] weightsOut, int edgeBlockSize){
		//assert(check());
//		eTotalOverOut=0;
		final float v=value();
		assert(v==values[lpos]) : v+", "+values[lpos]+//Also fires on NaN, but that shouldn't happen...
			"\n"+Arrays.toString(values)+
			"\n"+Arrays.toString(weights)+"\n";
		outOverNet=(float)derivativeXFX(sum, v);

		assert(!dense);
		eTotalOverOut=(weightsOut.length==eOverNetNext.length ? Vector.fma(weightsOut, eOverNetNext) :
			Vector.fma(weightsOut, eOverNetNext, outputs, 1, false));

//		if(CellNet.SIMD && eOverNetNext.length>=16) {
//			eTotalOverOut=shared.Vector.fma(weightsOut, eOverNetNext);
//		}else {
//			for(int dest=0; dest<eOverNetNext.length; dest++){
//				final float netOverOut=weightsOut[dest];
//				final float eOverNetDest=eOverNetNext[dest];
//				final float eOverOut=eOverNetDest*netOverOut;
//
//				eTotalOverOut+=eOverOut;
//			}
//		}

		eOverNetArray[lpos]=eOverNet=eTotalOverOut*outOverNet;

//		for(int source=0; source<valuesPrev.length; source++) {
//			final float netOverWeight=valuesPrev[source];
//			final float eTotalOverWeight=eOverNet*netOverWeight;
//
////			final double incr=weights[source]-alpha*eTotalOverWeight;
////			weights2[source]+=incr;
//
//			deltas[source]-=eTotalOverWeight;
//		}
		if(deltas.length==valuesIn.length){Vector.addProduct(deltas, valuesIn, -eOverNet);
		}else{Vector.addProduct(deltas, valuesIn, inputs, -eOverNet, edgeBlockSize);}

//		if(layer>0){
		{
			assert(layer>0) : layer;
//			bias2+=bias-alpha*eOverNet*biasAlphaMult;//Bias adjusts slower than edges
			biasDelta-=eOverNet;
		}
		//assert(check());
	}

	/*--------------------------------------------------------------*/

	/** Adds one thread-local training copy's accumulated deltas into this row. */
	public void accumulate(Cell c2){
		//assert(check());
		//assert(c2.check());
		error+=c2.error;
//		bias2+=c2.bias2;
		biasDelta+=c2.biasDelta;
//		if(layer>0) {
//			final double[] c2w2=c2.weights2;
//			for(int i=0; i<weights2.length; i++) {
//				weights2[i]+=c2w2[i];
//			}
//		}
		if(layer>0){
//			final float[] c2d=c2.deltas;
//			for(int i=0; i<deltas.length; i++) {
//				deltas[i]+=c2d[i];
//			}
			Vector.add(deltas, c2.deltas);
		}
	}

	/*--------------------------------------------------------------*/

	/** Adds bounded random perturbations to the bias and nonzero weights. */
	public void anneal(float strength, Random randy){
		if(layer<1){return;}
		if(annealBias){
			final float abs=Math.abs(bias);
			float xb=strength*(randy.nextFloat()-randy.nextFloat())*biasAnnealMult;
			if(abs<0.2f){xb=xb*(Tools.max(abs*5, 0.2f));} //Weaker anneal for weaker bias //TODO: Make lower limit even weaker
			setBias(bias+xb, true);
		}
		for(int i=0; i<weights.length; i++){
//			final Edge e=inputs.get(i);
			final float w=weights[i];
			if(w!=0){
				//			assert(w==e.weight());
				final float abs=Math.abs(w);
				float xe=strength*(randy.nextFloat()-randy.nextFloat());
				if(abs<lowWeightAnnealCutoff){xe=xe*(Tools.max(abs*lowWeightAnnealMult, lowWeightAnnealCutoff));} //Weaker anneal for weaker weight
				if(Math.abs(w+xe)>abs){xe*=edgeAmplitudeIncreaseMult;} //Weaker anneal when it increases absolute magnitude of weight
				final float w2=w+xe;
				//			e.setWeight(w2);
				if(Math.abs(w2)>Float.MIN_NORMAL){weights[i]=w2;}
				//			assert(w2==e.weight());
			}
		}
	}

	/** Copies row state into this existing row, preserving id, position and shared layer arrays. */
	public void setFrom(Cell c, boolean copyDelta){
		assert(dense==c.dense) : KillSwitch.assertDie(
			"Cell.setFrom reuses weight/index arrays and requires the same representation: destination dense="+
			dense+", source dense="+c.dense);
		sparseOutputGradient=c.sparseOutputGradient;
		eTotalOverOut=c.eTotalOverOut;
		outOverNet=c.outOverNet;
		bias=c.bias;
//		type=c.type;
		function=c.function;
		value=c.value;
		sum=c.sum;
		error=c.error;
//		if(layer==0){return;}
		if(weights==null){
			weights=(c.weights==null ? null : Arrays.copyOf(c.weights, c.weights.length));
			assert(inputs==null);
			inputs=(c.inputs==null ? null : Arrays.copyOf(c.inputs, c.inputs.length));
		}else{
			Vector.copy(weights, c.weights);
			if(!dense){Vector.copy(inputs, c.inputs);}
		}

		if(outputs==null){
			outputs=(c.outputs==null ? null : Arrays.copyOf(c.outputs, c.outputs.length));
		}else{
			Vector.copy(outputs, c.outputs);
		}

		biasDelta=0;
		if(copyDelta){
			biasDelta=c.biasDelta;
			deltas=(c.deltas==null ? null : Arrays.copyOf(c.deltas, c.deltas.length));
		}else if(deltas==null){
			deltas=(weights==null ? null : new float[weights.length]);//(c.deltas==null ? null : new float[c.deltas.length]);
		}else{
			Arrays.fill(deltas, 0);
		}
		assert(Tools.equals(inputs, c.inputs));
		assert(Tools.equals(outputs, c.outputs));
		assert(Tools.equals(weights, c.weights));
//		assert(false) : outputs+", "+c.outputs;

		//assert(check());
		//assert(c.check());
//		assert(c.weights.length==weights.length) : "\n"+inputs+"\n"+c.inputs;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Overrides           ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns a diagnostic row summary with bias, accumulated error and live edge weights. */
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
//		sb.append("C"+id()+": v="+String.format("%.4f, b=%.4f, b2=%.4f, e=%.5f", value(), bias, bias2, error));
		sb.append("C"+id()+": v="+String.format("%.4f, b=%.4f, b2=%.4f, e=%.5f", value(), bias, biasDelta, error));
		if(!terminal()){
			sb.append(", Edges: {");
			if(dense){
				int prevBase=id-lpos-weights.length;
				for(int i=0; i<weights.length; i++){
					//TODO: Could use source array or source value array here
//					return String.format("C"+(prevBase+i)+"->C"+id+",%.3f,w?=%.4f,w2=%.4f,d=%.4f; ", -1, weights[i], weights2[i], deltas[i]);
					sb.append(String.format("C"+(prevBase+i)+"->C"+id+",%.3f,w?=%.4f,d=%.4f; ", -1, weights[i], deltas[i]));
				}
			}else{
				for(int i=0; i<weights.length; i++){
					final int prev=inputs[i]+prevLayerStart;
					sb.append(String.format("C"+(prev)+"->C"+id+",%.3f,w?=%.4f,d=%.4f; ", -1, weights[i], deltas[i]));
				}
			}
			sb.setLength(sb.length()-2);
			sb.append("}");
		}
		return sb.toString();
	}

	/** Returns true for input-layer cells, which do not have incoming edges. */
	@Override
	public boolean terminal(){
//		return inputs==null || inputs.isEmpty();
		return layer==0;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the global cell identifier. */
	public int id(){return id;}
	/** Returns the current bias. */
	public float bias(){return bias;}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Applies boundary/class weighting to one signed raw output error, unless weightMult is zero. */
	public static float toWeightedError(double rawError, float v, float target, float weightMult){
		if(weightMult==0){return (float)rawError;}
//		assert(rawError>=0) : rawError;
		float mult=toErrorMult(v, target, weightMult);
		double incr=toErrorIncr(rawError, v, target);
		assert(incr==0 || incr>=0 == rawError>=0) : incr+", "+rawError;
		float ret=(float)((rawError+incr)*mult);
		assert(ret>=0 || rawError<0) : "raw="+rawError+", v="+v+", t="+target+
			", wm="+weightMult+", m="+mult+", incr="+incr+", ret="+ret;
		return ret;
	}

	/** Returns the signed boundary increment paired with the sign of rawError. */
	public static double toErrorIncr(double rawError, float v, float target){
		double incr=0;
		if(target>cutoffForTraining){
			if(v<=cutoffForTraining+spread){incr=fnErrorIncr;}
		}else if(target<cutoffForTraining){
			if(v>=cutoffForTraining-spread){incr=fpErrorIncr;}
		}
		double ret=incr*(rawError>=0 ? 1 : -1);
		return ret;
	}

//	public static float toErrorMult(float v, float target) {
//		if(true) {return toErrorMult(v, target, 1);}
//		final float mult;
////		assert(cutoff-spread==lowThresh);
////		assert(cutoff+spread==highThresh);
////		if(v>target) {
////			mult=positiveErrorMult;
////			if(v>lowThresh && target<lowThresh) {
////				mult*=falsePositiveErrorMult;
////			}
////		}else{
////			mult=negativeErrorMult;
////			if(v<highThresh && target>=highThresh) {
////				mult*=falseNegativeErrorMult;
////			}
////		}
//		if(v>target) {
//			if(v>cutoff-spread && target<cutoff+spread) {
//				mult=falsePositiveErrorMult;
//			}else{
//				mult=positiveErrorMult;
//			}
//		}else{
//			if(v<cutoff+spread && target>=cutoff-spread) {
//				mult=falseNegativeErrorMult;
//			}else{
//				mult=negativeErrorMult;
//			}
//		}
//		return mult;
//	}

	//TODO: multFraction comes from TrainerThread.weightMult() and is basically always 1
	//It should probably be eliminated
	//But test first; it's a way of preventing fpem from being too strong early via "minweightepoch=500" or whatever
	/** Returns zero for exact matches or the selected error multiplier times multFraction. */
	public static float toErrorMult(float v, float target, float multFraction){
		if(v==target){return 0;}
		final float mult;
		final boolean positiveError=v>target;
		final boolean positiveGoal=target>cutoffForTraining;
		final boolean negativeGoal=target<cutoffForTraining;
		final boolean excess=(positiveError == positiveGoal);
		boolean offsides=false;
//		final boolean offsides=(positiveGoal && v<)
		if(positiveError){
			if(positiveGoal){
				assert(excess);
				mult=excessPositiveErrorMult;
			}else if(v>cutoffForTraining-spread){//offsides; false positive
				mult=falsePositiveErrorMult;
				offsides=true;
			}else{
				mult=positiveErrorMult;
			}
		}else{//Negative error
			if(negativeGoal){
				assert(excess);
				mult=excessNegativeErrorMult;
			}else if(v<cutoffForTraining+spread){//offsides; false negative
				mult=falseNegativeErrorMult;
				offsides=true;
			}else{
				mult=negativeErrorMult;
			}
		}
		//multFraction=0.5, for example, returns halfway between 1.0 and mult
//		float ret=((mult-1)*multFraction)+1;

		float ret=mult*multFraction;//Changed to this because the old formula was giving a negative multiplier
		assert(ret>0) : "v="+v+", target="+target+"\n"+
				"mult="+mult+", ret="+ret+", mf="+multFraction+
				", pe="+positiveError+", pg="+positiveGoal+", ng="+negativeGoal+
				", excess="+excess+", offsides="+offsides;
		return ret;
	}

	/** Returns the configured weighted output error for {@code v-target}. */
	public static float calcETotalOverOut(float v, float target, float weightMult){
		float eTotalOverOut=v-target; //Aka out-target
		final float ret=toWeightedError(eTotalOverOut, v, target, weightMult);
		return ret;
	}

	/** Sets the annealing cutoff and reciprocal multiplier for small weights. */
	public static void setLowWeightAnnealCutoff(float c){
		assert(c>=0 && c<=1);
		lowWeightAnnealCutoff=c;
		lowWeightAnnealMult=Tools.max(1f, 1f/Tools.max(lowWeightAnnealCutoff, 0.0001f));
	}

	/*--------------------------------------------------------------*/

	/** Returns this row's activation-function name. */
	public final String typeString(){
		return function.name();
	}

	/*--------------------------------------------------------------*/

	//Fake legacy method
	/** Legacy double-array entry point retained only to fail accidental calls. */
	public void updateEdgesHiddenLayerDense(float alpha, float[] valuesIn, double[] eOverNetNext, float[] weightsOut){
		throw new RuntimeException("Wrong method, for legacy double[] version in class CellNetDouble");
	}

	//Fake legacy method
	/** Legacy double-array constructor retained only to fail accidental calls. */
	public Cell(int size, int type, int i, int layerNum, int prevWidth, int width, int nextWidth, float[] lvals,
			double[] eons){
		throw new RuntimeException("Legacy constructor for CellNetDouble");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Configured backprop signal with respect to this row's output activation. */
	public float eTotalOverOut;
	/** Activation derivative at the current pre-activation sum and output. */
	public float outOverNet;
	/** Configured backprop signal with respect to this row's pre-activation sum. */
	public float eOverNet;
	/** Weighted input sum before activation. */
	public double sum=0;

	/** Bias added to the weighted input sum. */
	public float bias;
//	private double bias2;
	/** Accumulated bias delta for the next epoch update. */
	private double biasDelta;
	/** Accumulated reporting loss for this row. */
	public double error;
	/** Immutable representation inherited from the owning network at construction. */
	public final boolean dense;
	/** Global cell identifier. */
	private final int id;
	/** Position within the current layer. */
	final int lpos;//position within layer
	/** Layer containing this cell. */
	final int layer;
	/** Output-layer index. */
	final int maxLayer;
	/** Index for loading serialized network weights sequentially. */
	int nextWeight=0; //For use when loading networks

//	int type;//TODO: Change to Function.
	/** Activation function implementation. */
	public Function function;

	/*--------------------------------------------------------------*/

	/** Global id of the first cell in the previous layer. */
	final int prevLayerStart;
	/** Global id of the first cell in the next layer. */
	final int nextLayerStart;

	//Lpos (layer position) of inputs
	/** Sparse incoming source positions, or null for dense rows. */
	public int[] inputs;
	/** Prepared with the training matrices; false preserves the original positional output loop. */
	private boolean sparseOutputGradient=false;
	/** Sparse outgoing destination positions used by hidden-row backpropagation. */
	public int[] outputs;
	/** Incoming edge weights. */
	public float[] weights;
	/** Accumulated incoming-edge deltas for the next epoch update. */
	public float[] deltas;

	/*--------------------------------------------------------------*/

	/** Shared activation values for all cells in this layer. */
	private final float[] values;
	/** Shared eOverNet values for all cells in this layer. */
	public final float[] eOverNetArray;

	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/

	/** Highest activation type eligible for the uniform hidden-layer fallback. */
	public static int MAX_TYPE=Function.TANH;
	/** Default activation type for input and hidden rows. */
	public static int defaultActivationType=Function.SIG;
	/** Default output-row activation type. */
	public static int finalLayerType=Function.RSLOG;
	/**
	 * Probability of uniform hidden-layer activation changes when weighted
	 * Function.TYPE_RATES_CUM selection is disabled.
	 */
	public static float randomTypeRate=0.0f;

	/** Learning-rate multiplier for bias updates. */
	static float biasAlphaMult=1f;
	/** Annealing strength multiplier for bias perturbations. */
	static float biasAnnealMult=0.5f;
	/** True to perturb bias values during annealing. */
	static boolean annealBias=true;

	/** Weight magnitude below which annealing perturbations are damped. */
	private static float lowWeightAnnealCutoff=0.2f;
	/** Reciprocal-based annealing multiplier for low-magnitude weights. */
	private static float lowWeightAnnealMult=1f/lowWeightAnnealCutoff;

	/** Classification threshold for binary-style training weights. */
	public static float cutoffForTraining=0.5f;
	/** True after a custom training cutoff has been set. */
	public static boolean setCutoffForTraining=false;
	/** True to infer a classification midpoint from observed labels. */
	public static boolean useMidpoint=false;

	/** Multiplier for ordinary positive prediction errors. */
	static float positiveErrorMult=1.0f;
	/** Multiplier for near-threshold false positives. */
	static float falsePositiveErrorMult=10.5f;
	//0.2 is best for binary classification; otherwise 1.0 is probably better
	//Should be paired with adjusting the Sample pivot function.
	/** Multiplier for excess positive predictions beyond the goal. */
	static float excessPositiveErrorMult=0.2f;

	/** Multiplier for ordinary negative prediction errors. */
	static float negativeErrorMult=1.0f;
	/** Multiplier for near-threshold false negatives. */
	static float falseNegativeErrorMult=10.5f;
	/** Multiplier for excess negative predictions below the goal. */
	static float excessNegativeErrorMult=0.2f;

	//Optimal BBMerge settings at the time
//	static float positiveErrorMult=1.65f;
//	static float falsePositiveErrorMult=12.5f;
//	static float excessPositiveErrorMult=1f;
//
//	static float negativeErrorMult=0.825f;
//	static float falseNegativeErrorMult=2.7f;
//	static float excessNegativeErrorMult=1f;

	/** Raw-error increment near the false-negative boundary. */
	static float fnErrorIncr=0.01f;
	/** Raw-error increment near the false-positive boundary. */
	static float fpErrorIncr=0.00f;
	/** Half-width around the cutoff for boundary-specific weighting. */
	static float spread=0.050f;
	/** Multiplier for edge changes that would grow a large magnitude. */
	static float edgeAmplitudeIncreaseMult=0.98f;
	/** Magnitude above which edge increases are damped. */
	static float edgeAmplitudeIncreaseThresh=0.1f;

}
