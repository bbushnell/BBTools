package hiseq;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Tools;
import stream.Read;

/**
 * Tracks base and quality composition per cycle.
 * 
 * @author Brian Bushnell
 * @date August 16, 2018
 *
 */
public class CycleTracker {
	
	public void add(Read r){
		if(r==null || r.length()<1){return;}
		addBases(r.bases);
		addQuality(r.quality);
	}
	
	public void add(CycleTracker ct){
//		assert(false) : length+", "+ct.length;
		final long[][] matrix=ct.acgtnq;
		if(matrix[0]==null){return;}
		if(acgtnq[0]==null || ct.length>length){
			resize(ct.length);
		}
		for(int i=0; i<matrix.length; i++){
			for(int j=0; j<matrix[i].length; j++){
				acgtnq[i][j]+=matrix[i][j];
			}
		}
	}
	
	//Do this before using results
	public void process(){
//		System.err.println("Processing); length="+length);
		if(length<1){return;}
		long[] cycleSum=new long[length];
		for(int i=0; i<5; i++){
			long[] array=acgtnq[i];
			for(int j=0; j<array.length; j++){
				cycleSum[j]+=array[j];
			}
		}
		maxes=new float[6];
		averages=new float[6];
		for(int i=0; i<6; i++){
			cycleAverages[i]=new float[length];
			long[] array=acgtnq[i];
			for(int j=0; j<length; j++){
				float f=array[j]/(float)cycleSum[j];
				cycleAverages[i][j]=f;
				maxes[i]=Tools.max(f, maxes[i]);
			}
		}
		
		long[] sum=new long[6];
		long sumsum=0;
		for(int i=0; i<5; i++){
			long x=shared.Vector.sum(acgtnq[i]);
			sum[i]=x;
			sumsum+=x;
		}
		sum[5]=shared.Vector.sum(acgtnq[5]);
		for(int i=0; i<6; i++){
			averages[i]=sum[i]/(float)sumsum;
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private void addBases(byte[] bases){
		if(length<bases.length){resize(bases.length);}
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumberACGTN[b];
			assert(x>=0) : "This program does not support degenerate base codes: "+new String(bases);
			acgtnq[x][i]++;
		}
	}
	
	private void addQuality(byte[] quals){
		if(quals==null){return;}
		if(length<quals.length){resize(quals.length);}
		final long[] sum=acgtnq[5];
		for(int i=0; i<quals.length; i++){
			byte q=quals[i];
			sum[i]+=q;
		}
	}
	
	private void resize(int newLen){
		assert(newLen>length);
		length=newLen;
		if(acgtnq[0]==null){
			for(int i=0; i<acgtnq.length; i++){
				acgtnq[i]=new long[newLen];
			}
		}else{
			for(int i=0; i<acgtnq.length; i++){
				acgtnq[i]=Arrays.copyOf(acgtnq[i], newLen);
			}
		}
	}
	
	public float max(byte base) {
		final int x=AminoAcid.baseToNumberACGTN[base];
		return(maxes[x]);
	}
	
	public float max(char base) {
		final int x=AminoAcid.baseToNumberACGTN[base];
		return(maxes[x]);
	}
	
	public float avg(byte base) {
		final int x=AminoAcid.baseToNumberACGTN[base];
		return(averages[x]);
	}
	
	public float avg(char base) {
		final int x=AminoAcid.baseToNumberACGTN[base];
		return(averages[x]);
	}
	
	public void clear() {
		for(long[] array : acgtnq) {Arrays.fill(array, 0);}
		for(float[] array : cycleAverages) {Arrays.fill(array, 0);}
		Arrays.fill(maxes, 0);
		Arrays.fill(averages, 0);
	}
	
	/*--------------------------------------------------------------*/

	public long[][] acgtnq=new long[6][];
	public float[][] cycleAverages=new float[6][];
	public float[] maxes;
	public float[] averages;
	public int length=0;
	
}
