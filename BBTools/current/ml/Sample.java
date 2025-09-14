package ml;

import java.util.Arrays;

import structures.ByteBuilder;

public class Sample implements Comparable<Sample> {
	
	public Sample(float[] in_, float[] out_, float weight_, int id_) {
		in=in_;
		goal=out_;
		result=new float[goal.length];
		weight=weight_;
		id=id_;
		positive=(goal[0]>=0.5f);
		assert(weight>0) : weight;
	}
	
	@Override
	public int compareTo(Sample o) {
		final float a=pivot, b=o.pivot;
		return a>b ? -1 : b>a ? 1 : id-o.id;
	}
	
	public boolean checkPivot() {
		return pivot==calcPivot();
	}
	
	synchronized void setPivot() {
		pivot=calcPivot();
	}
	
	synchronized float calcPivot() {
		final float v=result[0];
		final boolean positiveError=v>goal[0];
		final boolean excess=(positiveError == positive);
		final float mult=(excess ? excessPivotMult*0.5f : 0.5f);
		return (errorMagnitude+weightedErrorMagnitude)*mult-epoch*EPOCH_MULT;
//		return (errorMagnitude+weightedErrorMagnitude)*0.5f-epoch*EPOCH_MULT;
	}
	
	public String toString() {
//		String s="S%d\t%s\t%s\tep=%d\tg=%4f\tr=%4f\tem=%6f\tev=%.6f\tpv=%.6f";
		String s="S%d\t%s\t%s\tep=%d\tg=%4f\tr=%4f\tem=%6f\tpv=%.6f";
		
		
		boolean gol=(goal[0]>=0.5f);
		boolean pred=(result[0]>=0.5f);
		String type=(gol && pred) ? "TP" : (!gol && !pred) ? "TN" : (!gol && pred) ? "FP" : (gol && !pred) ? "FN" : "??";
		String sign=(positive ? "+" : "-");

//		s=String.format(s, id, sign, type, epoch, goal[0], result[0], errorMagnitude, errorValue, calcPivot());
		s=String.format(s, id, sign, type, epoch, goal[0], result[0], errorMagnitude, calcPivot());
		return s+"\t"+Arrays.toString(in);
	}
	
	public ByteBuilder toBytes() {
		return toBytes(new ByteBuilder());
	}
	
	public ByteBuilder toBytes(ByteBuilder bb) {
		for(float f : in) {bb.append(f, 6).tab();}
		for(float f : goal) {bb.append(f, 6).tab();}
		bb.trimLast(1);
		bb.nl();
		return bb;
	}
	
//	synchronized boolean positive() {
//		return goal[0]>=0.5f;
//	}
	
	public void calcError(float weightMult){
		double error=0;
		for(int i=0; i<result.length; i++){
			float r=result[i];
			float g=goal[i];
			float e=calcError(g, r);
			assert(e>=0);
			error+=e;
		}
		errorMagnitude=(float)error;
		assert(error>=0);
		weightedErrorMagnitude=Cell.toWeightedError(error, result[0], goal[0], weightMult);
		assert(weightedErrorMagnitude>=0);
	}
	
	public synchronized int epoch() {return epoch;}
	public synchronized int lastTID() {return lastTID;}
	public synchronized void setEpoch(long x) {
		epoch=(int)x;
	}
	
	public synchronized void setLastTID(int x) {
		lastTID=x;
	}
	
	public static final float calcError(float goal, float pred) {
		float e=goal-pred;
		return 0.5f*e*e;
	}

	final boolean positive;
	float errorMagnitude=1;
	float weightedErrorMagnitude=1;
//	float errorValue=1;//Unused, commented for efficiency
	private int epoch=-1;
	private int lastTID=-1;
	float pivot=0;
	
	final float[] in;
	final float[] goal;
	final float[] result;//Can't be volatile
	final float weight;//TODO
	final int id;

	//0.2f is good for binary classifiers
	public static float excessPivotMult=0.2f;
	public static final float EPOCH_MULT=1/256f;
}
