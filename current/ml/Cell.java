package ml;

import java.util.Arrays;
import java.util.Random;

import shared.Tools;
import shared.Vector;

public class Cell extends Source {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	public void summateDense(float[] valuesIn) {
		sum=bias;
		assert(valuesIn.length==weights.length) : valuesIn.length+", "+weights.length;
		sum+=Vector.fma(weights, valuesIn);
		final float v=(float)activation(sum);
		setValue(v);
	}
	
	public void summateSparse(float[] valuesIn, int edgeBlockSize) {
		sum=bias;
		sum+=Vector.fma(weights, valuesIn, inputs, edgeBlockSize, Vector.SIMD_FMA_SPARSE);
		final float v=(float)activation(sum);
		setValue(v);
	}
	
	public float calcError(float ideal) {
		float e=ideal-value();
		return 0.5f*e*e;
	}
	
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
	
	public void addError(float e) {
		assert(error>=0);
		error+=e;
	}
	
	public void clearError() {
		error=0;
	}
	
	public void clearTemp() {
//		bias2=0;
		biasDelta=0;
		error=0;
		if(layer<1) {return;}
//		Arrays.fill(weights2, 0);
		assert(deltas!=null) : id()+", "+layer;
		Arrays.fill(deltas, 0);
	}
	
	public void setBias(float b2, boolean ignoreAssertion) {
		assert(layer>0);
//		float dif=Tools.absdif(bias, b2);
//		assert(ignoreAssertion || dif<0.1f || dif<0.2f*bias) : dif+", "+bias+", "+b2; //just checking for bugs.
//		if(!(ignoreAssertion || dif<0.1f || dif<0.2f*bias)) {System.out.print("*");}
//		if(!ignoreAssertion) {System.err.println(bias);}
		bias=b2;
	}
	
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
	
	public final double activation(double x) {
		return function.activate(x);
	}
	
	public final double derivativeXFX(double x, double fx) {
		double d=function.derivativeXFX(x, fx);
		assert(!Double.isNaN(d)) : x+", "+fx+", "+d;
		return d;
	}
	
	/*--------------------------------------------------------------*/
	
	//TODO: Vectorize?  Final layer is usually small though
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

	public int id() {return id;}
	public float bias() {return bias;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	public static float calcETotalOverOut(float v, float target, float weightMult) {
		float eTotalOverOut=v-target; //Aka out-target
		final float ret=toWeightedError(eTotalOverOut, v, target, weightMult);
		return ret;
	}
	
	public static void setLowWeightAnnealCutoff(float c) {
		assert(c>=0 && c<=1);
		lowWeightAnnealCutoff=c;
		lowWeightAnnealMult=Tools.max(1f, 1f/Tools.max(lowWeightAnnealCutoff, 0.0001f));
	}
	
	/*--------------------------------------------------------------*/
	
	public final String typeString() {
		return function.name();
	}
	
	/*--------------------------------------------------------------*/
	
	//Fake legacy method
	public void updateEdgesHiddenLayerDense(float alpha, float[] valuesIn, double[] eOverNetNext, float[] weightsOut) {
		throw new RuntimeException("Wrong method, for legacy double[] version in class CellNetDouble");
	}

	//Fake legacy method
	public Cell(int size, int type, int i, int layerNum, int prevWidth, int width, int nextWidth, float[] lvals,
			double[] eons) {
		throw new RuntimeException("Legacy constructor for CellNetDouble");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public float eTotalOverOut;
	public float outOverNet;
	public float eOverNet;
	public double sum=0;
	
	public float bias;
//	private double bias2;
	private double biasDelta;
	public double error;
	private final int id;
	final int lpos;//position within layer
	final int layer;
	final int maxLayer;
	int nextWeight=0; //For use when loading networks
	
//	int type;//TODO: Change to Function.
	Function function;
	
	/*--------------------------------------------------------------*/

	final int prevLayerStart;
	final int nextLayerStart;
	
	//Lpos (layer position) of inputs
	public int[] inputs;
	public int[] outputs;
	public float[] weights;
	float[] deltas;
	
	/*--------------------------------------------------------------*/
	
	private final float[] values; 
	public final float[] eOverNetArray;
	
	/*--------------------------------------------------------------*/
	/*----------------           Statics            ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int MAX_TYPE=Function.TANH;
	public static int defaultActivationType=Function.SIG;
	public static int finalLayerType=Function.RSLOG;
	public static float randomTypeRate=0.0f;

	static float biasAlphaMult=1f;
	static float biasAnnealMult=0.5f;
	static boolean annealBias=true;
	
	private static float lowWeightAnnealCutoff=0.2f;
	private static float lowWeightAnnealMult=1f/lowWeightAnnealCutoff;
	
	static float cutoffForTraining=0.5f;
	static boolean setCutoffForTraining=false;
	static boolean useMidpoint=false;

	static float positiveErrorMult=1.0f;
	static float falsePositiveErrorMult=10.5f;
	//0.2 is best for binary classification; otherwise 1.0 is probably better
	//Should be paired with adjusting the Sample pivot function.
	static float excessPositiveErrorMult=0.2f;

	static float negativeErrorMult=1.0f;
	static float falseNegativeErrorMult=10.5f;
	static float excessNegativeErrorMult=0.2f;

	//Optimal BBMerge settings at the time
//	static float positiveErrorMult=1.65f;
//	static float falsePositiveErrorMult=12.5f;
//	static float excessPositiveErrorMult=1f;
//
//	static float negativeErrorMult=0.825f;
//	static float falseNegativeErrorMult=2.7f;
//	static float excessNegativeErrorMult=1f;
	
	static float fnErrorIncr=0.01f;
	static float fpErrorIncr=0.00f;
	static float spread=0.050f;
	static float edgeAmplitudeIncreaseMult=0.98f;
	static float edgeAmplitudeIncreaseThresh=0.1f;
	
}
