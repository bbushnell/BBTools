package hiseq;

import java.util.Arrays;

import dna.AminoAcid;
import shared.Tools;
import stream.Read;

/**
 * Tracks base and quality composition across sequencing cycles for DNA reads.
 * Accumulates statistics about nucleotide bases (A, C, G, T, N) and their quality
 * scores across different read positions. Rows contain A/C/G/T/N counts and
 * numeric quality sums; quality means divide by observed bases. Instances and
 * exposed arrays require exclusive ownership while adding, processing or clearing.
 *
 * @author Brian Bushnell
 * @date August 16, 2018
 */
public class CycleTracker{

	/** Adds a nonempty read's bases and optional qualities; null/empty reads are ignored. */
	public void add(Read r){
		if(r==null || r.length()<1){return;}
		addBases(r.bases);
		addQuality(r.quality);
	}

	/** Adds another quiescent tracker's raw counters without copying its derived arrays. */
	public void add(CycleTracker ct){
//		assert(false) : length+", "+ct.length;
		final long[][] matrix=ct.acgtnq;
		if(matrix[0]==null){return;}//ct empty (never added to); nothing to merge
		//resize's assert(newLen>length) holds: matrix[0]!=null here => ct.length>=1, and acgtnq[0]==null
		//iff our length==0, so the null case has ct.length>=1>0; the other case is guarded by ct.length>length.
		if(acgtnq[0]==null || ct.length>length){
			resize(ct.length);
		}
		for(int i=0; i<matrix.length; i++){
			for(int j=0; j<matrix[i].length; j++){
				acgtnq[i][j]+=matrix[i][j];
			}
		}
	}

	/** Rebuilds derived fractions/quality means; retained cycles with no bases yield zero. */
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
				//HIL-101: clear() retains storage; unused cycles after reset have no observations.
				float f=cycleSum[j]==0 ? 0 : array[j]/(float)cycleSum[j];
				cycleAverages[i][j]=f;
				maxes[i]=Tools.max(f, maxes[i]);
			}
		}

		long[] sum=new long[6];
		long sumsum=0;
		for(int i=0; i<5; i++){
			long x=simd.Vector.sum(acgtnq[i]);
			sum[i]=x;
			sumsum+=x;
		}
		sum[5]=simd.Vector.sum(acgtnq[5]);
		for(int i=0; i<6; i++){
			averages[i]=sumsum==0 ? 0 : sum[i]/(float)sumsum;
		}
	}

	/*--------------------------------------------------------------*/

	/** Counts supported A/C/G/T/N codes at their zero-based cycles, growing storage as needed. */
	private void addBases(byte[] bases){
		if(length<bases.length){resize(bases.length);}
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int x=AminoAcid.baseToNumberACGTN[b];
			assert(x>=0) : "This program does not support degenerate base codes: "+new String(bases);
			acgtnq[x][i]++;
		}
	}

	/** Adds numeric qualities to the quality-sum row; absent qualities contribute no sum. */
	private void addQuality(byte[] quals){
		if(quals==null){return;}
		if(length<quals.length){resize(quals.length);}
		final long[] sum=acgtnq[5];
		for(int i=0; i<quals.length; i++){
			byte q=quals[i];
			sum[i]+=q;
		}
	}

	/** Extends every raw row together; newLen must exceed the retained extent. */
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

	/** Maximum cycle fraction for a supported base, after process() has populated results. */
	public float max(byte base){
		final int x=AminoAcid.baseToNumberACGTN[base];
		return(maxes[x]);
	}

	/** Character overload of max(byte); requires populated derived results. */
	public float max(char base){
		final int x=AminoAcid.baseToNumberACGTN[base];
		return(maxes[x]);
	}

	/** Fraction of all observed bases matching base, after process() has populated results. */
	public float avg(byte base){
		final int x=AminoAcid.baseToNumberACGTN[base];
		return(averages[x]);
	}

	/** Character overload of avg(byte); requires populated derived results. */
	public float avg(char base){
		final int x=AminoAcid.baseToNumberACGTN[base];
		return(averages[x]);
	}

	/** Clears raw and existing derived values while retaining allocated cycle storage. */
	public void clear(){
		//[hiseq/CycleTracker#001] FIXED (G11 2026-07-01, Brian-authorized): acgtnq[i] is null until resize()
		//(via an add), and cycleAverages[i]/maxes/averages are null until process() -- so the old unguarded
		//Arrays.fill(...) NPE'd when clear() ran on a fresh/added-but-unprocessed tracker. Now null-guarded so
		//clear() is valid at any lifecycle stage (a null array is already "empty" -> filling it is a no-op).
		//Reachability at time of fix: LATENT -- the only caller is MicroTile.clear(), which itself has NO
		//callers tree-wide (dead path). Fixed anyway per the existence rule (provably NPE-prone logic).
		for(long[] array : acgtnq){if(array!=null){Arrays.fill(array, 0);}}
		for(float[] array : cycleAverages){if(array!=null){Arrays.fill(array, 0);}}
		if(maxes!=null){Arrays.fill(maxes, 0);}
		if(averages!=null){Arrays.fill(averages, 0);}
	}

	/*--------------------------------------------------------------*/

	/** A/C/G/T/N counts and quality sums; rows are allocated by the first nonempty addition. */
	public long[][] acgtnq=new long[6][];
	/** Per-cycle fractions/mean qualities, rebuilt by process(). */
	public float[][] cycleAverages=new float[6][];
	/** Maximum per-cycle value for each row, allocated by process(). */
	public float[] maxes;
	/** Whole-read-collection value for each row, allocated by process(). */
	public float[] averages;
	/** Retained raw-row extent; clear() zeroes observations without shrinking this storage. */
	public int length=0;

}
