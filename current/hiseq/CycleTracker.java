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
	
	/**
	 * Adds a read to the cycle tracker, updating base and quality statistics.
	 * Skips null or empty reads.
	 * @param r The read to add for cycle tracking
	 */
	public void add(Read r){
		if(r==null || r.length()<1){return;}
		addBases(r.bases);
		addQuality(r.quality);
	}
	
	/**
	 * Merges statistics from another CycleTracker into this instance.
	 * Resizes internal arrays if the other tracker has longer sequences.
	 * @param ct The CycleTracker to merge into this instance
	 */
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
	/**
	 * Processes accumulated data to compute cycle averages and statistics.
	 * Must be called before accessing max() or avg() results.
	 * Calculates per-cycle base frequencies and overall base distributions.
	 */
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
	
	/**
	 * Adds base composition data from a sequence to the tracker.
	 * Resizes arrays if necessary to accommodate sequence length.
	 * @param bases The base sequence to add
	 */
	private void addBases(byte[] bases){
		if(length<bases.length){resize(bases.length);}
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumberACGTN[b];
			assert(x>=0) : "This program does not support degenerate base codes: "+new String(bases);
			acgtnq[x][i]++;
		}
	}
	
	/**
	 * Adds quality score data to the tracker.
	 * Skips processing if quality array is null.
	 * @param quals The quality scores to add
	 */
	private void addQuality(byte[] quals){
		if(quals==null){return;}
		if(length<quals.length){resize(quals.length);}
		final long[] sum=acgtnq[5];
		for(int i=0; i<quals.length; i++){
			byte q=quals[i];
			sum[i]+=q;
		}
	}
	
	/**
	 * Resizes internal arrays to accommodate longer sequences.
	 * Creates new arrays or extends existing ones using Arrays.copyOf().
	 * @param newLen The new array length (must be greater than current length)
	 */
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
	
	/**
	 * Returns the maximum frequency observed for the specified base across all cycles.
	 * @param base The base to query (A, C, G, T, or N)
	 * @return Maximum frequency for the base across all cycles
	 */
	public float max(byte base) {
		final int x=AminoAcid.baseToNumberACGTN[base];
		return(maxes[x]);
	}
	
	/**
	 * Returns the maximum frequency observed for the specified base across all cycles.
	 * @param base The base to query (A, C, G, T, or N)
	 * @return Maximum frequency for the base across all cycles
	 */
	public float max(char base) {
		final int x=AminoAcid.baseToNumberACGTN[base];
		return(maxes[x]);
	}
	
	/**
	 * Returns the average frequency for the specified base across all sequences.
	 * @param base The base to query (A, C, G, T, or N)
	 * @return Average frequency for the base across all sequences
	 */
	public float avg(byte base) {
		final int x=AminoAcid.baseToNumberACGTN[base];
		return(averages[x]);
	}
	
	/**
	 * Returns the average frequency for the specified base across all sequences.
	 * @param base The base to query (A, C, G, T, or N)
	 * @return Average frequency for the base across all sequences
	 */
	public float avg(char base) {
		final int x=AminoAcid.baseToNumberACGTN[base];
		return(averages[x]);
	}
	
	/** Resets all internal arrays and statistics to zero.
	 * Clears base counts, cycle averages, maxes, and averages. */
	public void clear() {
		for(long[] array : acgtnq) {Arrays.fill(array, 0);}
		for(float[] array : cycleAverages) {Arrays.fill(array, 0);}
		Arrays.fill(maxes, 0);
		Arrays.fill(averages, 0);
	}
	
	/*--------------------------------------------------------------*/

	/** 2D array storing base and quality counts: [ACGTN+quality][cycle] */
	public long[][] acgtnq=new long[6][];
	/** Per-cycle averages for each base type and quality */
	public float[][] cycleAverages=new float[6][];
	/** Maximum frequencies observed for each base type */
	public float[] maxes;
	/** Overall average frequencies for each base type */
	public float[] averages;
	/** Current maximum sequence length being tracked */
	public int length=0;
	
}
