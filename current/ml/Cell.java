package ml;

import java.util.Arrays;
import shared.Random;

import shared.Tools;
import simd.Vector;

/**
 * Represents a neural network cell with advanced computational and training capabilities.
 * Handles neural network cell operations including activation, error computation, edge updates,
 * and weight management for both dense and sparse network architectures.
 * Supports multiple activation functions through dynamic function selection.
 *
 * @author Brian Bushnell
 * @date 2025
 */
public class Cell extends Source {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs a neural network cell with specified configuration parameters.
	 * Initializes cell identification, activation function, layer positioning, and shared arrays.
	 *
	 * @param id_ Unique identifier for this cell
	 * @param activationType Type of activation function to use
	 * @param lpos_ Position within the layer
	 * @param layer_ Current layer number
	 * @param maxLayer_ Maximum layer number in the network
	 * @param prevLayerStart_ Starting index of the previous layer
	 * @param nextLayerStart_ Starting index of the next layer
	 * @param wid Width of the current layer
	 * @param values_ Shared array storing activation values for all cells in layer
	 * @param eOverNetArray_ Shared array storing error over net values for all cells in layer
	 */
	public Cell(int id_, int activationType, int lpos_, int layer_, int maxLayer_, 
			int prevLayerStart_, int nextLayerStart_, int wid,
			float[] values_, float[] eOverNetArray_) {
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
	
	/**
	 * Computes weighted sum and activation for dense network connections.
	 * Uses vectorized operations for efficient computation of dot product between weights and inputs.
	 * @param valuesIn Input values from the previous layer
	 */
	public void summateDense(float[] valuesIn) {
		sum=bias;
		assert(valuesIn.length==weights.length) : valuesIn.length+", "+weights.length;
		sum+=Vector.fma(weights, valuesIn);
		final float v=(float)activation(sum);
		setValue(v);
	}
	
	/**
	 * Computes weighted sum and activation for sparse network connections.
	 * Uses sparse matrix multiplication with specified edge block size for efficiency.
	 * @param valuesIn Input values from the previous layer
	 * @param edgeBlockSize Block size for sparse operations
	 */
	public void summateSparse(float[] valuesIn, int edgeBlockSize) {
		sum=bias;
		sum+=Vector.fma(weights, valuesIn, inputs, edgeBlockSize, Vector.SIMD_FMA_SPARSE);
		final float v=(float)activation(sum);
		setValue(v);
	}
	
	/**
	 * Calculates squared error between current cell value and ideal target value.
	 * @param ideal Target value for this cell
	 * @return Squared error (0.5 * (ideal - actual)^2)
	 */
	public float calcError(float ideal) {
		float e=ideal-value();
		return 0.5f*e*e;
	}
	
	/**
	 * Performs comprehensive validation of cell state and connectivity.
	 * Verifies array consistency, weight initialization, and network topology constraints.
	 * @return true if cell passes all validation checks, false otherwise
	 */
	synchronized public boolean check() {
//		assert(false);
		
		if(!CellNet.DENSE) {
			assert(layer==maxLayer || outputs!=null) :
				layer+", "+maxLayer+", "+lpos+", "+id();
		}
		
		if(value!=values[lpos] && (value!=-1 && values[lpos]!=0)) {
			assert(false) : id+", "+layer+", "+lpos+": "+value+", "+values[lpos]+", "+Arrays.toString(values);
			return false;
		}
		if(eOverNet!=eOverNetArray[lpos]) {
			assert(false) : id+", "+layer+", "+lpos+": "+eOverNet+", "+eOverNetArray[lpos]+", "+Arrays.toString(values);
			return false;
		}
		if(layer>0) {
			assert(weights!=null);
			assert(CellNet.DENSE==(inputs==null));
//			assert(weights2.length==weights.length);
			assert(deltas==null || deltas.length==weights.length);
//			for(int i=0; i<weights2.length; i++) {
//				Edge e=outputs.get(i);
//				if(e.weight2!=weights2[i]) {
//					assert(false) : id+", "+layer+", "+lpos+": "+e.weight2+", "+weights2[i];
//					return false;
//				}
//			}
		}else {
//			assert(inputs==null || outputs.isEmpty())  : id+", "+layer+", "+lpos+": "+weights2+": "+outputs+"\n"+toString()+"\n";
			assert(inputs==null || outputs==null)  : id+", "+layer+", "+lpos+": "+deltas+": "+outputs+"\n"+toString()+"\n";
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Mutators           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Applies accumulated weight and bias updates after an epoch of training.
	 * Includes amplitude control to prevent excessive weight growth that can break convergence.
	 * @param invSamples Inverse of the number of samples processed
	 * @param alpha Learning rate multiplier
	 */
	void applyUpdates(float invSamples, float alpha) {//only called once per epoch
		if(layer<1) {return;}
		adjustBias(invSamples, alpha);
		for(int i=0; i<weights.length; i++) {
			final float w=weights[i];
			
			final float d3=(float)(deltas[i]*invSamples*alpha);
			float w4=w+d3;
			{
				final float absW4=Math.abs(w4);
				if(absW4>Math.abs(w) && absW4>edgeAmplitudeIncreaseThresh) {
					//Very sensitive; strong values break convergence.
					//0.98 seems OK though.  Goal is to keep edges low.
					//.94 is too extreme.  .96 seems OK.
					w4=w+edgeAmplitudeIncreaseMult*d3;//Slow down increase in edge magnitude
				}
			}
			deltas[i]=0;

			if((w!=0 || !CellNet.DENSE) && Math.abs(w4)>Float.MIN_NORMAL) {weights[i]=w4;}
		}
	}
	
	/** Accumulates error value to the cell's total error.
	 * @param e Error value to add */
	public void addError(float e) {
		assert(error>=0);
		error+=e;
	}
	
	/** Resets the cell's accumulated error to zero */
	public void clearError() {
		error=0;
	}
	
	/**
	 * Clears temporary training variables including bias delta and weight deltas
	 */
	public void clearTemp() {
//		bias2=0;
		biasDelta=0;
		error=0;
		if(layer<1) {return;}
//		Arrays.fill(weights2, 0);
		assert(deltas!=null) : id()+", "+layer;
		Arrays.fill(deltas, 0);
	}
	
	/**
	 * Sets the bias value for this cell.
	 * @param b2 New bias value
	 * @param ignoreAssertion Whether to skip validation assertions
	 */
	public void setBias(float b2, boolean ignoreAssertion) {
		assert(layer>0);
//		float dif=Tools.absdif(bias, b2);
//		assert(ignoreAssertion || dif<0.1f || dif<0.2f*bias) : dif+", "+bias+", "+b2; //just checking for bugs.
//		if(!(ignoreAssertion || dif<0.1f || dif<0.2f*bias)) {System.out.print("*");}
//		if(!ignoreAssertion) {System.err.println(bias);}
		bias=b2;
	}
	
	/**
	 * Updates bias using accumulated bias delta with learning rate scaling.
	 * @param invSamples Inverse of the number of samples processed
	 * @param alpha Learning rate
	 */
	private void adjustBias(float invSamples, float alpha) {
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
	
	/** Adds error from another cell to this cell's error accumulator.
	 * @param c2 Source cell to add error from */
	public void addError(Cell c2) {
		error+=c2.error;
	}
	
	@Override
	public void setValue(float v) {
		assert(layer==0 || v==(float)activation(sum)) : v+", "+activation(sum)+", "+sum;
//		assert(value==values[lpos]);
		values[lpos]=value=v;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Applies the cell's activation function to the input value.
	 * @param x Input value to activate
	 * @return Activated output value
	 */
	public final double activation(double x) {
		return function.activate(x);
	}
	
	/**
	 * Calculates the derivative of the activation function for backpropagation.
	 * @param x Input value
	 * @param fx Function output value
	 * @return Derivative value for gradient calculation
	 */
	public final double derivativeXFX(double x, double fx) {
		double d=function.derivativeXFX(x, fx);
		assert(!Double.isNaN(d)) : x+", "+fx+", "+d;
		return d;
	}
	
	/*--------------------------------------------------------------*/
	
	//TODO: Vectorize?  Final layer is usually small though
	/**
	 * Updates weights and gradients for output layer cells in dense networks.
	 * Computes error derivatives and accumulates weight updates based on target values.
	 *
	 * @param target Target output value for this cell
	 * @param valuesIn Input values from previous layer
	 * @param weightMult Weight multiplier for error calculation
	 */
	void updateEdgesFinalLayerDense(float target, float[] valuesIn, float weightMult) {
		//assert(check());
		final float v=value();
		eTotalOverOut=calcETotalOverOut(v, target, weightMult);
		outOverNet=(float)derivativeXFX(sum, v);
//		final double eTotalOverOut_X_outOverNet=eTotalOverOut*outOverNet;
		eOverNet=eOverNetArray[lpos]=eTotalOverOut*outOverNet;
		
		for(int i=0; i<weights.length; i++) {
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
	
	/**
	 * Updates weights and gradients for hidden layer cells in dense networks.
	 * Propagates error backwards from next layer using chain rule.
	 *
	 * @param valuesIn Input values from previous layer
	 * @param eOverNetNext Error over net values from next layer
	 * @param weightsOut Weights connecting to next layer
	 */
	public void updateEdgesHiddenLayerDense(float[] valuesIn, float[] eOverNetNext, float[] weightsOut) {
		//assert(check());
//		eTotalOverOut=0;
		final float v=value();
		assert(v==values[lpos]) : v+", "+values[lpos]+//Also fires on NaN, but that shouldn't happen...
			"\n"+Arrays.toString(values)+
			"\n"+Arrays.toString(weights)+"\n";
		outOverNet=(float)derivativeXFX(sum, v);

		assert(CellNet.DENSE);
		if(!CellNet.SPECIAL_FMA) {
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
	
	/**
	 * Updates weights and gradients for hidden layer cells in sparse networks.
	 * Uses sparse matrix operations for efficient computation with limited connections.
	 *
	 * @param valuesIn Input values from previous layer
	 * @param eOverNetNext Error over net values from next layer
	 * @param weightsOut Weights connecting to next layer
	 * @param edgeBlockSize Block size for sparse operations
	 */
	public void updateEdgesHiddenLayerSparse(float[] valuesIn, float[] eOverNetNext,
			float[] weightsOut, int edgeBlockSize) {
		//assert(check());
//		eTotalOverOut=0;
		final float v=value();
		assert(v==values[lpos]) : v+", "+values[lpos]+//Also fires on NaN, but that shouldn't happen...
			"\n"+Arrays.toString(values)+
			"\n"+Arrays.toString(weights)+"\n";
		outOverNet=(float)derivativeXFX(sum, v);

		assert(!CellNet.DENSE);
		eTotalOverOut=Vector.fma(weightsOut, eOverNetNext, outputs, 1, false);
		
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
		Vector.addProduct(deltas, valuesIn, inputs, -eOverNet, edgeBlockSize);
		
//		if(layer>0){
		{
			assert(layer>0) : layer;
//			bias2+=bias-alpha*eOverNet*biasAlphaMult;//Bias adjusts slower than edges
			biasDelta-=eOverNet;
		}
		//assert(check());
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Accumulates training state from another cell including error, bias delta, and weight deltas.
	 * Used for parallel training where multiple threads process the same network structure.
	 * @param c2 Source cell to accumulate from
	 */
	public void accumulate(Cell c2) {
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
		if(layer>0) {
//			final float[] c2d=c2.deltas;
//			for(int i=0; i<deltas.length; i++) {
//				deltas[i]+=c2d[i];
//			}
			Vector.add(deltas, c2.deltas);
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Applies simulated annealing by adding random noise to weights and bias.
	 * Uses scaled random perturbations to help escape local minima during training.
	 * @param strength Strength of annealing perturbations
	 * @param randy Random number generator for reproducible annealing
	 */
	public void anneal(float strength, Random randy) {
		if(layer<1){return;}
		if(annealBias) {
			final float abs=Math.abs(bias);
			float xb=strength*(randy.nextFloat()-randy.nextFloat())*biasAnnealMult;
			if(abs<0.2f) {xb=xb*(Tools.max(abs*5, 0.2f));} //Weaker anneal for weaker bias //TODO: Make lower limit even weaker
			setBias(bias+xb, true);
		}
		for(int i=0; i<weights.length; i++) {
//			final Edge e=inputs.get(i);
			final float w=weights[i];
			if(w!=0) {
				//			assert(w==e.weight());
				final float abs=Math.abs(w);
				float xe=strength*(randy.nextFloat()-randy.nextFloat());
				if(abs<lowWeightAnnealCutoff) {xe=xe*(Tools.max(abs*lowWeightAnnealMult, lowWeightAnnealCutoff));} //Weaker anneal for weaker weight
				if(Math.abs(w+xe)>abs) {xe*=edgeAmplitudeIncreaseMult;} //Weaker anneal when it increases absolute magnitude of weight
				final float w2=w+xe;
				//			e.setWeight(w2);
				if(Math.abs(w2)>Float.MIN_NORMAL) {weights[i]=w2;}
				//			assert(w2==e.weight());
			}
		}
	}
	
	/**
	 * Copies complete cell state from another cell including weights, biases, and connections.
	 * Used for network cloning and parallel training setup.
	 * @param c Source cell to copy from
	 * @param copyDelta Whether to copy delta values for gradient accumulation
	 */
	public void setFrom(Cell c, boolean copyDelta) {
		eTotalOverOut=c.eTotalOverOut;
		outOverNet=c.outOverNet;
		bias=c.bias;
//		type=c.type;
		function=c.function;
		value=c.value;
		sum=c.sum;
		error=c.error;
//		if(layer==0){return;}
		if(weights==null) {
			weights=(c.weights==null ? null : Arrays.copyOf(c.weights, c.weights.length));
			assert(inputs==null);
			inputs=(c.inputs==null ? null : Arrays.copyOf(c.inputs, c.inputs.length));
		}else{
			Vector.copy(weights, c.weights);
			if(!CellNet.DENSE) {Vector.copy(inputs, c.inputs);}
		}
		
		if(outputs==null) {
			outputs=(c.outputs==null ? null : Arrays.copyOf(c.outputs, c.outputs.length));
		}else{
			Vector.copy(outputs, c.outputs);
		}

		biasDelta=0;
		if(copyDelta) {
			biasDelta=c.biasDelta;
			deltas=(c.deltas==null ? null : Arrays.copyOf(c.deltas, c.deltas.length));
		}else if(deltas==null){
			deltas=(weights==null ? null : new float[weights.length]);//(c.deltas==null ? null : new float[c.deltas.length]);
		}else {
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
	
	@Override
	public String toString() {
		StringBuilder sb=new StringBuilder();
//		sb.append("C"+id()+": v="+String.format("%.4f, b=%.4f, b2=%.4f, e=%.5f", value(), bias, bias2, error));
		sb.append("C"+id()+": v="+String.format("%.4f, b=%.4f, b2=%.4f, e=%.5f", value(), bias, biasDelta, error));
		if(!terminal()) {
			sb.append(", Edges: {");
			if(CellNet.DENSE) {
				int prevBase=id-lpos-weights.length;
				for(int i=0; i<weights.length; i++) {
					//TODO: Could use source array or source value array here
//					return String.format("C"+(prevBase+i)+"->C"+id+",%.3f,w?=%.4f,w2=%.4f,d=%.4f; ", -1, weights[i], weights2[i], deltas[i]);
					sb.append(String.format("C"+(prevBase+i)+"->C"+id+",%.3f,w?=%.4f,d=%.4f; ", -1, weights[i], deltas[i]));
				}
			}else {
				for(int i=0; i<weights.length; i++) {
					final int prev=inputs[i]+prevLayerStart;
					sb.append(String.format("C"+(prev)+"->C"+id+",%.3f,w?=%.4f,d=%.4f; ", -1, weights[i], deltas[i]));
				}
			}
			sb.setLength(sb.length()-2);
			sb.append("}");
		}
		return sb.toString();
	}

	@Override
	public boolean terminal() {
//		return inputs==null || inputs.isEmpty();
		return layer==0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the unique identifier for this cell */
	public int id() {return id;}
	/** Returns the current bias value for this cell */
	public float bias() {return bias;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Converts raw error to weighted error using configurable class weighting.
	 * Applies different error multipliers based on prediction accuracy and target class.
	 *
	 * @param rawError Raw error value before weighting
	 * @param v Current predicted value
	 * @param target Target value
	 * @param weightMult Weight multiplier for error scaling
	 * @return Weighted error value
	 */
	public static float toWeightedError(double rawError, float v, float target, float weightMult) {
		if(weightMult==0) {return (float)rawError;}
//		assert(rawError>=0) : rawError;
		float mult=toErrorMult(v, target, weightMult);
		double incr=toErrorIncr(rawError, v, target);
		assert(incr==0 || incr>=0 == rawError>=0) : incr+", "+rawError;
		float ret=(float)((rawError+incr)*mult);
		assert(ret>=0 || rawError<0) : "raw="+rawError+", v="+v+", t="+target+
			", wm="+weightMult+", m="+mult+", incr="+incr+", ret="+ret;
		return ret;
	}
	
	/**
	 * Calculates additional error increment for predictions near classification boundaries.
	 * Adds penalty for false positives and false negatives near the decision threshold.
	 *
	 * @param rawError Raw error value
	 * @param v Current predicted value
	 * @param target Target value
	 * @return Error increment to add
	 */
	public static double toErrorIncr(double rawError, float v, float target) {
		double incr=0;
		if(target>cutoffForTraining) {
			if(v<=cutoffForTraining+spread) {incr=fnErrorIncr;}
		}else if(target<cutoffForTraining) {
			if(v>=cutoffForTraining-spread) {incr=fpErrorIncr;}
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
	/**
	 * Calculates error multiplier based on prediction type and target class.
	 * Applies different multipliers for correct predictions, false positives, and false negatives.
	 *
	 * @param v Current predicted value
	 * @param target Target value
	 * @param multFraction Fraction for scaling the multiplier effect
	 * @return Error multiplier to apply
	 */
	public static float toErrorMult(float v, float target, float multFraction) {
		if(v==target) {return 0;}
		final float mult;
		final boolean positiveError=v>target;
		final boolean positiveGoal=target>cutoffForTraining;
		final boolean negativeGoal=target<cutoffForTraining;
		final boolean excess=(positiveError == positiveGoal);
		boolean offsides=false;
//		final boolean offsides=(positiveGoal && v<)
		if(positiveError) {
			if(positiveGoal) {
				assert(excess);
				mult=excessPositiveErrorMult;
			}else if(v>cutoffForTraining-spread){//offsides; false positive
				mult=falsePositiveErrorMult;
				offsides=true;
			}else {
				mult=positiveErrorMult;
			}
		}else{//Negative error
			if(negativeGoal) {
				assert(excess);
				mult=excessNegativeErrorMult;
			}else if(v<cutoffForTraining+spread){//offsides; false negative
				mult=falseNegativeErrorMult;
				offsides=true;
			}else {
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
	
	/**
	 * Calculates total error over output for gradient computation.
	 * Combines raw error difference with weighted error adjustments.
	 *
	 * @param v Current output value
	 * @param target Target output value
	 * @param weightMult Weight multiplier for error weighting
	 * @return Total error over output value
	 */
	public static float calcETotalOverOut(float v, float target, float weightMult) {
		float eTotalOverOut=v-target; //Aka out-target
		final float ret=toWeightedError(eTotalOverOut, v, target, weightMult);
		return ret;
	}
	
	/**
	 * Sets the threshold below which weights receive reduced annealing strength.
	 * Prevents small weights from being perturbed too strongly during annealing.
	 * @param c Cutoff threshold for low weight annealing (0.0 to 1.0)
	 */
	public static void setLowWeightAnnealCutoff(float c) {
		assert(c>=0 && c<=1);
		lowWeightAnnealCutoff=c;
		lowWeightAnnealMult=Tools.max(1f, 1f/Tools.max(lowWeightAnnealCutoff, 0.0001f));
	}
	
	/*--------------------------------------------------------------*/
	
	/** Returns the name of the activation function used by this cell.
	 * @return String name of the activation function */
	public final String typeString() {
		return function.name();
	}
	
	/*--------------------------------------------------------------*/
	
	//Fake legacy method
	/**
	 * Legacy method placeholder that throws RuntimeException to indicate deprecated double[] version
	 */
	public void updateEdgesHiddenLayerDense(float alpha, float[] valuesIn, double[] eOverNetNext, float[] weightsOut) {
		throw new RuntimeException("Wrong method, for legacy double[] version in class CellNetDouble");
	}

	//Fake legacy method
	/**
	 * Legacy constructor placeholder that throws RuntimeException to indicate CellNetDouble compatibility
	 */
	public Cell(int size, int type, int i, int layerNum, int prevWidth, int width, int nextWidth, float[] lvals,
			double[] eons) {
		throw new RuntimeException("Legacy constructor for CellNetDouble");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Total error propagated from output layer for gradient computation */
	public float eTotalOverOut;
	/** Derivative of activation function for backpropagation */
	public float outOverNet;
	/** Error over net value for this cell in backpropagation */
	public float eOverNet;
	/** Weighted sum of inputs before activation function is applied */
	public double sum=0;
	
	/** Bias value added to weighted sum before activation */
	public float bias;
//	private double bias2;
	/** Accumulated bias gradient for weight updates */
	private double biasDelta;
	/** Accumulated error for this cell during training */
	public double error;
	/** Unique identifier for this cell */
	private final int id;
	/** Position index within the current layer */
	final int lpos;//position within layer
	/** Layer number containing this cell */
	final int layer;
	/** Maximum layer number in the neural network */
	final int maxLayer;
	/** Index for loading network weights sequentially */
	int nextWeight=0; //For use when loading networks
	
//	int type;//TODO: Change to Function.
	/** Activation function implementation for this cell */
	Function function;
	
	/*--------------------------------------------------------------*/

	/** Starting index of the previous layer in global cell array */
	final int prevLayerStart;
	/** Starting index of the next layer in global cell array */
	final int nextLayerStart;
	
	//Lpos (layer position) of inputs
	/** Array of input connection indices for sparse networks */
	public int[] inputs;
	/** Array of output connection indices for sparse networks */
	public int[] outputs;
	/** Weights for incoming connections from previous layer */
	public float[] weights;
	/** Accumulated weight gradients for gradient descent updates */
	float[] deltas;
	
	/*--------------------------------------------------------------*/
	
	/** Shared array containing activation values for all cells in this layer */
	private final float[] values; 
	/** Shared array containing error over net values for all cells in this layer */
	public final float[] eOverNetArray;
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Maximum activation function type value */
	public static int MAX_TYPE=Function.TANH;
	/** Default activation function type for new cells */
	public static int defaultActivationType=Function.SIG;
	/** Activation function type used for output layer cells */
	public static int finalLayerType=Function.RSLOG;
	/** Rate of random activation function selection during network generation */
	public static float randomTypeRate=0.0f;

	/** Learning rate multiplier for bias updates */
	static float biasAlphaMult=1f;
	/** Annealing strength multiplier for bias perturbations */
	static float biasAnnealMult=0.5f;
	/** Whether to apply annealing to bias values */
	static boolean annealBias=true;
	
	/** Threshold below which weights receive reduced annealing strength */
	private static float lowWeightAnnealCutoff=0.2f;
	/** Multiplier for annealing strength on low weights */
	private static float lowWeightAnnealMult=1f/lowWeightAnnealCutoff;
	
	/** Classification threshold for binary training (typically 0.5) */
	static float cutoffForTraining=0.5f;
	/** Whether custom training cutoff has been set */
	static boolean setCutoffForTraining=false;
	/** Whether to use midpoint for classification decisions */
	static boolean useMidpoint=false;

	/** Error multiplier for positive prediction errors */
	static float positiveErrorMult=1.0f;
	/** Error multiplier for false positive predictions */
	static float falsePositiveErrorMult=10.5f;
	//0.2 is best for binary classification; otherwise 1.0 is probably better
	//Should be paired with adjusting the Sample pivot function.
	/** Error multiplier for excessive positive predictions beyond target */
	static float excessPositiveErrorMult=0.2f;

	/** Error multiplier for negative prediction errors */
	static float negativeErrorMult=1.0f;
	/** Error multiplier for false negative predictions */
	static float falseNegativeErrorMult=10.5f;
	/** Error multiplier for excessive negative predictions beyond target */
	static float excessNegativeErrorMult=0.2f;

	//Optimal BBMerge settings at the time
//	static float positiveErrorMult=1.65f;
//	static float falsePositiveErrorMult=12.5f;
//	static float excessPositiveErrorMult=1f;
//
//	static float negativeErrorMult=0.825f;
//	static float falseNegativeErrorMult=2.7f;
//	static float excessNegativeErrorMult=1f;
	
	/** Additional error increment for false negative predictions */
	static float fnErrorIncr=0.01f;
	/** Additional error increment for false positive predictions */
	static float fpErrorIncr=0.00f;
	/** Spread around cutoff threshold for classification boundaries */
	static float spread=0.050f;
	/**
	 * Multiplier to reduce weight amplitude increases that can break convergence
	 */
	static float edgeAmplitudeIncreaseMult=0.98f;
	/** Threshold above which weight amplitude increases are dampened */
	static float edgeAmplitudeIncreaseThresh=0.1f;
	
}
