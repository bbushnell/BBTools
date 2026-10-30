package simd;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import dna.AminoAcid;
import ml.Cell;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;

/** Java-8-compatible facade for optional vector kernels and scalar fallbacks.
 * Capability probes avoid directly linking Vector API types into this class.
 * Most methods select a fixed-width backend using {@link Shared#SIMD} and a
 * length threshold; capability fields describe available hardware/API support,
 * not whether a particular call will use that backend. The 128-bit and 64-bit
 * capability flags do not select alternate implementations of 256-bit kernels.
 * Arrays are borrowed unless a method explicitly returns a new copy. Callers
 * must coordinate mutations and changes to process-wide execution flags.
 * @author Brian Bushnell
 * @date Sep 12, 2013
 *
 */
public final class Vector{

	/**
	 * Performance benchmark comparing scalar and SIMD implementations.
	 * Tests sum and max operations on integer arrays using both execution paths.
	 * @param args Command-line arguments (unused)
	 */
	public static void main(String[] args){
		int[] array=new int[19999];
		for(int i=0; i<array.length; i++){
			array[i]=(byte)(i&127);
		}
		Timer t=new Timer();

		int loops=80000;

		long sum=0, sum2=0, max=0;
		for(int outer=0; outer<4; outer++){
			t.start();
			sum=0; sum2=0;
			Shared.SIMD=false;
			for(int i=0; i<loops; i++){
				sum=sum(array);
				sum2+=sum;
			}
			System.err.println(sum2);
			t.stop("Scalar: ");

			t.start();
			sum=0; sum2=0;
			Shared.SIMD=true;
			for(int i=0; i<loops; i++){
				sum=sum(array);
				sum2+=sum;
			}
			System.err.println(sum2);
			t.stop("Vector: ");
		}

		for(int outer=0; outer<4; outer++){
			t.start();
			max=0; sum=0; sum2=0;
			Shared.SIMD=false;
			for(int i=0; i<loops; i++){
				max=max(array);
				sum2+=max;
			}
			System.err.println(max+", "+sum2);
			t.stop("Scalar: ");

			t.start();
			max=0; sum=0; sum2=0;
			Shared.SIMD=true;
			for(int i=0; i<loops; i++){
				max=max(array);
				sum2+=max;
			}
			System.err.println(max+", "+sum2);
			t.stop("Vector: ");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Preserves the historical public no-argument constructor; operations are static. */
	public Vector(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns "c+=a[i]*b[i]" where a and b are equal-length arrays.
	 * @param a A vector to multiply.
	 * @param b A vector to multiply.
	 * @return Sum of products of vector elements.
	 */
	public static final float fma(final float[] a, final float[] b){
		assert(a.length==b.length);
		if(Shared.SIMD && Shared.SIMD_FMA && a.length>=MINLEN32){return SIMD.fma(a, b);}
		float c=0;
		for(int i=0; i<a.length; i++){c+=a[i]*b[i];}
		return c;
	}

	/**
	 * Computes one query dot product against every contiguous row of a matrix.
	 * The matrix may have padding between rows; only the first {@code q.length}
	 * values of each row participate.  The SIMD reduction is not bit-identical
	 * to the scalar sequential summation, and results are not bitwise
	 * reproducible across JIT states; callers comparing results should use a
	 * numerical tolerance.
	 * @param q Query vector.
	 * @param mat Row-major matrix with row {@code r} starting at {@code r*stride}.
	 * @param stride Matrix row stride, at least {@code q.length}.
	 * @param nRows Number of rows to score.
	 * @param out Destination for one score per row; the first {@code nRows} entries are written.
	 */
	public static final void dotRows(final float[] q, final float[] mat, final int stride,
			final int nRows, final float[] out){
		assert(q!=null && mat!=null && out!=null);
		assert(stride>=q.length && stride>0 && nRows>=0 && out.length>=nRows);
		assert((long)nRows*stride<=mat.length || nRows==0);
		if(Shared.SIMD && q.length>=MINLEN32){
			SIMD.dotRows(q, mat, stride, nRows, out);
			return;
		}
		for(int row=0; row<nRows; row++){
			final int base=row*stride;
			float sum=0;
			for(int d=0; d<q.length; d++){sum+=q[d]*mat[base+d];}
			out[row]=sum;
		}
	}

	/**
	 * Selects the top {@code n} score indices, sorted by score descending and
	 * then index ascending. NaN scores sort after all non-NaN scores.
	 * @param scores Scores indexed by family row.
	 * @param n Number of indices to write.
	 * @param outIdx Destination with room for at least {@code n} indices.
	 */
	public static final void topN(final float[] scores, final int n, final int[] outIdx){
		if(scores==null || outIdx==null || n<0 || n>scores.length || n>outIdx.length){
			throw new IllegalArgumentException("Invalid topN arguments: scores="+
				(scores==null ? -1 : scores.length)+" n="+n+" outIdx="+(outIdx==null ? -1 : outIdx.length));
		}
		if(n==0){return;}
		int size=0;
		for(int idx=0; idx<scores.length; idx++){
			if(size==n){
				final int worst=outIdx[size-1];
				if(!betterTop(scores[idx], idx, scores[worst], worst)){continue;}
			}
			int pos=(size<n ? size : size-1);
			while(pos>0 && betterTop(scores[idx], idx, scores[outIdx[pos-1]], outIdx[pos-1])){pos--;}
			if(size<n){size++;}
			for(int j=size-1; j>pos; j--){outIdx[j]=outIdx[j-1];}
			outIdx[pos]=idx;
		}
	}

	/** Tests whether A should displace B during the ascending-index top-N scan.
	 * Two NaNs never displace one another, preserving their encounter order.
	 * @param scoreA Candidate score
	 * @param indexA Candidate index
	 * @param scoreB Existing score
	 * @param indexB Existing index
	 * @return True for a better non-NaN score, or an earlier equal-score index */
	private static boolean betterTop(final float scoreA, final int indexA, final float scoreB, final int indexB){
		final boolean nanA=Float.isNaN(scoreA), nanB=Float.isNaN(scoreB);
		return nanB ? !nanA : (!nanA && (scoreA>scoreB || (scoreA==scoreB && indexA<indexB)));
	}

	/**
	 * Computes the indexed dot product a[i]*b[bSet[i]] without modifying the inputs.
	 * Equal lengths alone do not imply an identity mapping. Callers that already
	 * know the mapping is identity may use the dense overload.
	 * @param a Multipliers, with the same length as bSet
	 * @param b Indexed source values
	 * @param bSet Valid indices into b
	 * @param blockSize Positive contiguous index-block size; use one for arbitrary
	 * mappings. Multiples of eight permit vector loads from each block's first index.
	 * @param allowSimd Whether this call may use sparse SIMD; false preserves scalar order
	 * @return Sum of products; SIMD fusion/reduction order can change rounding
	 */
	public static final float fma(final float[] a, final float[] b, final int[] bSet,
		final int blockSize, final boolean allowSimd){
		assert(a.length==bSet.length);
		//FIXED [simd/Vector#002]: equal array lengths did not prove bSet was identity.
		if(Shared.SIMD && Shared.SIMD_FMA && allowSimd && a.length>=MINLEN32
				&& a.length<b.length && blockSize>0 && ((blockSize&7)==0)){
			return SIMD.fmaSparse(a, b, bSet);
		}
		float c=0;
		for(int i=0; i<a.length; i++){c+=a[i]*b[bSet[i]];}
		return c;
	}

	/**
	 * Executes forward propagation for a neural network layer.
	 * Computes weighted sums, applies bias, and calculates activation values.
	 * @param layer Array of neural network cells to update
	 * @param valuesIn Input values from previous layer
	 */
	public static final void feedForward(Cell[] layer, float[] valuesIn){
		//		assert(layer.length==valuesIn.length);
		if(Shared.SIMD && Shared.SIMD_FEED_FORWARD && valuesIn.length>=MINLEN32){
			SIMD.feedForward(layer, valuesIn);
			return;
		}

		for(int cnum=0; cnum<layer.length; cnum++){
			Cell c=layer[cnum];
			float[] weights=c.weights;
			float sum=c.bias;
			assert(valuesIn.length==weights.length) : valuesIn.length+", "+weights.length;
			sum+=Vector.fma(valuesIn, weights);
			c.sum=sum;
			final float v=(float)c.activation(sum);
			c.setValue(v);
		}
	}

	/**
	 * Dense forward propagation.
	 * @param layer Array of neural network cells to update
	 * @param valuesIn Input values from previous layer
	 */
	public static final void feedForwardDense(Cell[] layer, float[] valuesIn){
		feedForwardDense(layer, valuesIn, true);//Back-compat: callers with no per-job control fall through to the global gate
	}

	/** Computes dense cell activations with a per-job gate for the batched kernel.
	 * Disabling simdFF skips that kernel; individual fallback dot products still
	 * honor the independent Shared.SIMD_FMA flag and may use SIMD.
	 * @param layer Cells whose sums and output values are updated
	 * @param valuesIn Inputs matching every cell's dense weight-array length
	 * @param simdFF Whether this call may use the batched feed-forward kernel */
	public static final void feedForwardDense(Cell[] layer, float[] valuesIn, boolean simdFF){
		//		assert(layer.length==valuesIn.length);
		if(simdFF && Shared.SIMD && Shared.SIMD_FEED_FORWARD && valuesIn.length>=MINLEN32){//Batched-kernel gate; fallback dot products separately honor SIMD_FMA.
			SIMD.feedForward(layer, valuesIn);
			return;
		}

		for(int cnum=0; cnum<layer.length; cnum++){
			Cell c=layer[cnum];
			float[] weights=c.weights;
			float sum=c.bias;
			assert(valuesIn.length==weights.length) : valuesIn.length+", "+weights.length;
			sum+=Vector.fma(valuesIn, weights);
			c.sum=sum;
			final float v=(float)c.activation(sum);
			c.setValue(v);
		}
	}

	/**
	 * Backpropagation fused multiply-add for neural network training.
	 * Computes error gradients using matrix multiplication between weights and errors.
	 *
	 * @param layer Layer cells to update with error gradients
	 * @param eOverNetNext Error derivatives from next layer
	 * @param weightsOutLnum Weight matrix connecting this layer to next layer
	 */
	public static void backPropFma(Cell[] layer, float[] eOverNetNext, float[][] weightsOutLnum){
		if(Shared.SIMD && Shared.SIMD_BACKPROP && eOverNetNext.length>=MINLEN32){
			SIMD.backPropFma(layer, eOverNetNext, weightsOutLnum);
			return;
		}
		for(int i=0; i<layer.length; i++){
			Cell cell=layer[i];
			cell.eTotalOverOut=Vector.fma(weightsOutLnum[i], eOverNetNext);
		}
	}

	/**
	 * Performs "a+=incr" where a and incr are equal-length arrays.
	 * @param a A vector to increment.
	 * @param incr Increment amount.
	 */
	public static final void add(final float[] a, final float[] incr){
		assert(a.length==incr.length);
		if(Shared.SIMD && a.length>=MINLEN32){SIMD.add(a, incr); return;}
		for(int i=0; i<a.length; i++){a[i]+=incr[i];}
	}

	/**
	 * Computes sum of absolute differences between two float vectors.
	 * @param a First vector
	 * @param b Second vector
	 * @return Sum of |a[i] - b[i]| for all elements
	 */
	public static final float absDifFloat(float[] a, float[] b){
		if(Shared.SIMD && a.length>=MINLEN32){return SIMD.absDifFloat(a, b);}
		assert(a.length==b.length);
		float sum=0;
		for(int i=0; i<a.length; i++){
			sum+=Math.abs(a[i]-b[i]);
		}
		return (float)sum;
	}

	/**
	 * SIMD-accelerated DDL bucket comparison.
	 * Compares two DynamicDemiLog maxArray char arrays, counting lower/equal/higher
	 * buckets while skipping empty-empty pairs.
	 * @param a First DDL maxArray
	 * @param b Second DDL maxArray
	 * @return int[]{lower, equal, higher, bothEmpty}
	 */
	public static int[] compareDDL(char[] a, char[] b){
		if(Shared.SIMD){return SIMDLogLog.compareDetailed(a, b);}
		return compareDDLScalar(a, b);
	}

	/**
	 * Hybrid SIMD-accelerated DDL comparison: pre-converted short[] query vs char[] ref.
	 * Pre-convert the query once via toShortArray(), then compare against many char[] refs.
	 * @param queryShort Pre-converted query maxArray (from toShortArray)
	 * @param refArray Reference DDL maxArray (native char[])
	 * @return int[]{lower, equal, higher, bothEmpty}
	 */
	public static int[] compareDDL(short[] queryShort, char[] refArray){
		if(Shared.SIMD){return SIMDLogLog.compareDetailedShort(queryShort, refArray, queryShort.length);}
		return compareDDLScalarHybrid(queryShort, refArray);
	}

	/**
	 * Compares pre-converted DDL arrays as unsigned 16-bit bucket values.
	 * Both inputs retain char bit patterns in short storage. Pre-converted arrays
	 * can be reused across comparisons; each call allocates its four-count result.
	 * @param queryShort Pre-converted query maxArray
	 * @param refShort Reference storage with at least queryShort.length elements
	 * @return int[]{lower, equal, higher, bothEmpty}
	 */
	public static int[] compareDDL(short[] queryShort, short[] refShort){
		if(Shared.SIMD){return SIMDLogLog.compareDetailedShort(queryShort, refShort, queryShort.length);}
		return compareDDLScalarShort(queryShort, refShort);
	}

	/** Copies DDL bucket bit patterns into new short storage for repeated comparisons.
	 * @param arr Non-null source bucket array
	 * @return New array of equal length with the same 16-bit patterns */
	public static short[] toShortArray(final char[] arr){
		final short[] out=new short[arr.length];
		for(int i=0; i<arr.length; i++){
			out[i]=(short)arr[i];
		}
		return out;
	}

	/** Compares unsigned short patterns against char buckets, separating empty pairs.
	 * @param a Query buckets
	 * @param b Reference storage with at least a.length elements
	 * @return New counts in lower, equal, higher, both-empty order */
	private static int[] compareDDLScalarHybrid(final short[] a, final char[] b){
		int lower=0, equal=0, higher=0, bothEmpty=0;
		for(int i=0; i<a.length; i++){
			final int ai=a[i]&0xFFFF, bi=b[i];
			if(ai==0 && bi==0){bothEmpty++; continue;}
			int dif=ai-bi;
			int nbit=(dif>>>31);
			int hbit=((-dif)>>>31);
			int ebit=1-nbit-hbit;
			lower+=nbit;
			higher+=hbit;
			equal+=ebit;
		}
		return new int[]{lower, equal, higher, bothEmpty};
	}

	/** Compares unsigned short patterns, excluding empty pairs from equality counts.
	 * @param a Query buckets
	 * @param b Reference storage with at least a.length elements
	 * @return New counts in lower, equal, higher, both-empty order */
	private static int[] compareDDLScalarShort(final short[] a, final short[] b){
		int lower=0, equal=0, higher=0, bothEmpty=0;
		for(int i=0; i<a.length; i++){
			final int ai=a[i]&0xFFFF, bi=b[i]&0xFFFF;
			if(ai==0 && bi==0){bothEmpty++; continue;}
			int dif=ai-bi;
			int nbit=(dif>>>31);
			int hbit=((-dif)>>>31);
			int ebit=1-nbit-hbit;
			lower+=nbit;
			higher+=hbit;
			equal+=ebit;
		}
		return new int[]{lower, equal, higher, bothEmpty};
	}

	/** Compares char buckets, excluding empty pairs from equality counts.
	 * @param a Query buckets
	 * @param b Reference storage with at least a.length elements
	 * @return New counts in lower, equal, higher, both-empty order */
	private static int[] compareDDLScalar(final char[] a, final char[] b){
		int lower=0, equal=0, higher=0, bothEmpty=0;
		for(int i=0; i<a.length; i++){
			final int ai=a[i], bi=b[i];
			if(ai==0 && bi==0){bothEmpty++; continue;}
			int dif=ai-bi;
			int nbit=(dif>>>31);
			int hbit=((-dif)>>>31);
			int ebit=1-nbit-hbit;
			lower+=nbit;
			higher+=hbit;
			equal+=ebit;
		}
		return new int[]{lower, equal, higher, bothEmpty};
	}

	/**
	 * GC-content compensation for k-mer frequency arrays.
	 * Normalizes values based on GC content to reduce compositional bias.
	 *
	 * @param a Array of k-mer counts
	 * @param k K-mer length
	 * @param gcmap GC content mapping for each k-mer
	 * @return GC-compensated frequency array
	 */
	public static float[] compensate(long[] a, int k, int[] gcmap){
		final float[] aSum=new float[k+1];
		final float inv=1f/(k+1);

		for(int i=0; i<a.length; i++){
			int gc=gcmap[i];
			aSum[gc]+=a[i];
		}

		for(int i=0; i<aSum.length; i++){
			aSum[i]=inv/Math.max(aSum[i], 1);
		}

		float[] comp=new float[a.length];
		for(int i=0; i<a.length; i++){
			int gc=gcmap[i];
			comp[i]=a[i]*aSum[gc];
		}
		//This just needs to add to approximately 1.
		//    	assert(Tools.sum(comp)==1) : "k="+k+", "+Tools.sum(comp)+"\n"+Arrays.toString(aSum)+"\n"+Arrays.toString(gcmap)+"\n"+Arrays.toString(a)+"\n";
		return comp;
	}

	/**
	 * GC-content compensation for k-mer frequency arrays (integer version).
	 * Normalizes values based on GC content to reduce compositional bias.
	 *
	 * @param a Array of k-mer counts
	 * @param k K-mer length
	 * @param gcmap GC content mapping for each k-mer
	 * @return GC-compensated frequency array
	 */
	public static float[] compensate(int[] a, int k, int[] gcmap){
		final float[] aSum=new float[k+1];
		final float inv=1f/(k+1);

		for(int i=0; i<a.length; i++){
			int gc=gcmap[i];
			aSum[gc]+=a[i];
		}

		for(int i=0; i<aSum.length; i++){
			aSum[i]=inv/Math.max(aSum[i], 1);
		}

		float[] comp=new float[a.length];
		for(int i=0; i<a.length; i++){
			int gc=gcmap[i];
			comp[i]=a[i]*aSum[gc];
		}
		//This just needs to add to approximately 1.
		//    	assert(Tools.sum(comp)==1) : "k="+k+", "+Tools.sum(comp)+"\n"+Arrays.toString(aSum)+"\n"+Arrays.toString(gcmap)+"\n"+Arrays.toString(a)+"\n";
		return comp;
	}



	/**
	 * Computes GC-compensated absolute difference between two k-mer frequency arrays.
	 * Applies GC normalization before calculating sum of absolute differences.
	 *
	 * @param a First k-mer frequency array
	 * @param b Second k-mer frequency array
	 * @param k K-mer length
	 * @param gcmap GC content mapping for each k-mer
	 * @return GC-compensated absolute difference (0-1 range)
	 */
	public static float absDifComp(long[] a, long[] b, int k, int[] gcmap){
		//    	if(Shared.SIMD && a.length>Vector.MINLEN32) {return SIMD.absDifComp(a, b, k, gcmap);}
		float[] af=compensate(a, k, gcmap);
		float[] bf=compensate(b, k, gcmap);
		float ret=Vector.absDifFloat(af, bf);
		return Tools.mid(0, 1, (Float.isFinite(ret) && ret>0 ? ret : 0));
	}

	/**
	 * Computes cosine difference (1 - cosine similarity) between two vectors.
	 * @param a First vector
	 * @param b Second vector
	 * @return Cosine difference value (0 = identical, higher = more different)
	 */
	public static float cosineDifference(int[] a, int[] b){
		float inva=1f/Math.max(1, sum(a));
		float invb=1f/Math.max(1, sum(b));
		float ret=1-cosineSimilarity(a, b, inva, invb);
		return (Float.isFinite(ret) && ret>0 ? ret : 0);
	}

	/**
	 * Computes cosine similarity between two normalized vectors.
	 * Uses pre-computed inverse normalization factors for efficiency.
	 *
	 * @param a First vector
	 * @param b Second vector
	 * @param inva Inverse normalization factor for vector a
	 * @param invb Inverse normalization factor for vector b
	 * @return Cosine similarity (1 = identical, 0 = orthogonal)
	 */
	public static float cosineSimilarity(int[] a, int[] b, float inva, float invb){
		assert(a.length==b.length);
		if(Shared.SIMD && a.length>=MINLEN32){return SIMD.cosineSimilarity(a, b, inva, invb);}
		float dotProduct=0f;
		float normVec1=0f;
		float normVec2=0f;

		for(int i=0; i<a.length; i++){
			float ai=a[i]*inva, bi=b[i]*invb;
			dotProduct+=ai*bi;
			normVec1+=ai*ai;
			normVec2+=bi*bi;
		}

		normVec1=Math.max(normVec1, 1e-15f);
		normVec2=Math.max(normVec2, 1e-15f);
		return (float)(dotProduct/(Math.sqrt(normVec1)*Math.sqrt(normVec2)));
	}

	/**
	 * Computes cosine similarity between float arrays
	 * @param a First float array
	 * @param b Second float array
	 * @return Cosine similarity (1.0 = identical, 0.0 = orthogonal)
	 */
	public static float cosineSimilarity(float[] a, float[] b){
		assert(a.length==b.length);
		if(Shared.SIMD && a.length>=MINLEN32){return SIMD.cosineSimilarity(a, b);}
		float dotProduct=0f;
		float normVec1=0f;
		float normVec2=0f;
		for(int i=0; i<a.length; i++){
			float ai=a[i];
			float bi=b[i];
			dotProduct+=ai*bi;
			normVec1+=ai*ai;
			normVec2+=bi*bi;
		}

		if(normVec1==0 || normVec2==0){return 0;} //Avoid NaN
		return (float)(dotProduct/(Math.sqrt(normVec1)*Math.sqrt(normVec2)));
	}

	/**
	 * Performs "a[i]+=b[i]*mult" where a and b are equal-length arrays.
	 * @param a A vector to increment.
	 * @param b Increment amount.
	 * @param mult Increment multiplier.
	 */
	public static final void addProduct(final float[] a, final float[] b, final float mult){
		assert(a.length==b.length);
		if(Shared.SIMD && a.length>=MINLEN32){SIMD.addProduct(a, b, mult); return;}
		for(int i=0; i<a.length; i++){a[i]+=b[i]*mult;}
	}

	/**
	 * Adds b[bSet[i]]*mult into each destination entry.
	 * Equal source and destination lengths do not imply identity indices. Callers
	 * that know the mapping is identity may use the dense overload.
	 * @param a Destination, with the same length as bSet
	 * @param b Source values; must not alias destination writes that affect later reads
	 * @param bSet Valid indices into b
	 * @param mult Increment multiplier
	 * @param blockSize Positive contiguous index-block size; use one for arbitrary mappings
	 */
	public static final void addProduct(final float[] a, final float[] b, final int[] bSet,
			final float mult, final int blockSize){
		assert(a.length==bSet.length);
		//FIXED [simd/Vector#002]: retain indexed semantics even for full-width rows.
		if(Shared.SIMD && a.length>=MINLEN32 && SIMD_MULT_SPARSE && a.length<b.length
				&& blockSize>0 && ((blockSize&7)==0)){
			SIMD.addProductSparse(a, b, bSet, mult);
			return;
		}
		for(int i=0; i<a.length; i++){a[i]+=b[bSet[i]]*mult;}
	}

	/**
	 * Copies the shorter array's length from source to destination.
	 * Any remaining destination suffix is unchanged. Both arrays must be non-null;
	 * passing the same array for both arguments is allowed.
	 * @param dest Destination array, modified in place
	 * @param source Source array, otherwise unchanged
	 */
	public static void copy(final float[] dest, final float[] source){
		//		assert(a.length==b.length);
		if(SIMDCOPY && Shared.SIMD && dest.length>=MINLEN32){SIMD.copy(dest, source); return;}
		for(int i=0, max=Tools.min(dest.length, source.length); i<max; i++){dest[i]=source[i];}
	}

	/**
	 * Copies the shorter array's length from source to destination using scalar code.
	 * Any remaining destination suffix is unchanged. Both arrays must be non-null;
	 * passing the same array for both arguments is allowed.
	 * @param dest Destination array, modified in place
	 * @param source Source array, otherwise unchanged
	 */
	public static void copy(final int[] dest, final int[] source){
		//		assert(a.length==b.length);
		//		if(SIMDCOPY && Shared.SIMD && dest.length>=MINLEN32) {SIMD.copy(dest, source); return;}//TODO
		for(int i=0, max=Tools.min(dest.length, source.length); i<max; i++){dest[i]=source[i];}
	}

	/** Counts equal bytes in two valid, equally sized inclusive ranges.
	 * Equality is literal; ambiguous bases are not treated as wildcards.
	 * @param s1 First sequence
	 * @param s2 Second sequence
	 * @param a1 First index in s1, inclusive
	 * @param b1 Last index in s1, inclusive
	 * @param a2 First index in s2, inclusive
	 * @param b2 Last index in s2, inclusive
	 * @return Number of equal byte pairs; neither sequence is modified */
	public static final int countMatches(final byte[] s1, final byte[] s2,
			final int a1, final int b1, final int a2, final int b2){
		assert(b1-a1==b2-a2) : a1+"-"+b1+", "+a2+"-"+b2+", len="+s1.length+", "+(b1-a1)+"!="+(b2-a2);
		if(Shared.SIMD && b1-a1+1>=MINLEN8){return SIMDByte256.countMatches(s1, s2, a1, b1, a2, b2);}
		int matches=0;
		for(int i=a1, j=a2; j<=b2; i++, j++){
			final byte x=s1[i], y=s2[j];
			final int m=((x==y) ? 1 : 0);//Does not take into account capitalization or undefined bases
			matches+=m;
		}
		assert(matches>=0 && matches<=b1-a1+1);
		return matches;
	}

	/**
	 * Counts unequal bytes in two valid, equally sized inclusive ranges.
	 * Equality is literal; ambiguous bases are not treated as wildcards.
	 *
	 * @param s1 First sequence
	 * @param s2 Second sequence
	 * @param a1 Start position in first sequence
	 * @param b1 End position in first sequence, inclusive
	 * @param a2 Start position in second sequence
	 * @param b2 End position in second sequence, inclusive
	 * @return Number of mismatching bases
	 */
	public static final int countMismatches(final byte[] s1, final byte[] s2,
			final int a1, final int b1, final int a2, final int b2){
		return (b1-a1+1)-countMatches(s1, s2, a1, b1, a2, b2);
	}

	/**
	 * Finds first occurrence of a symbol in array range.
	 *
	 * @param a Array to search
	 * @param symbol Symbol to find
	 * @param from Start position (inclusive)
	 * @param to End position (exclusive)
	 * @return Index of first occurrence, or 'to' if not found
	 */
	public static final int find(final byte[] a, final byte symbol, final int from, final int to){
		//		if(Shared.SIMD && to-from>=MINLEN8) {return SIMDByte.find(a, symbol, from, to);}
		int len=from;
		while(len<to && a[len]!=symbol){len++;}
		return len;
	}


	/**
	 * Sums float elements, returning a double without modifying the input.
	 * Scalar code accumulates in double. SIMD retains float-vector accumulation
	 * for ordinary results and recomputes a nonfinite result in double. Finite
	 * inputs cannot overflow the return value, but SIMD can round differently.
	 * @param array Array to sum
	 * @return Sum of all elements, or 0 if array is null
	 */
	public static double sum(final float[] array){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN32){return SIMD.sum(array, 0, array.length-1);}
		double x=0;
		for(float y : array){x+=y;}
		return x;
	}

	/**
	 * Computes sum of all elements in byte array.
	 * @param array Array to sum
	 * @return Sum of all elements, or 0 if array is null
	 */
	public static long sum(final byte[] array){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN8){return SIMDByte256.sum(array, 0, array.length-1);}
		long x=0;
		for(byte y : array){x+=y;}
		return x;
	}

	/**
	 * Computes sum of all elements in char array.
	 * @param array Array to sum
	 * @return Sum of all elements, or 0 if array is null
	 */
	public static long sum(final char[] array){
		if(array==null){return 0;}
		//		if(Shared.SIMD && array.length>=SMINLEN) {return SIMD.sum(array, 0, array.length-1);}
		long x=0;
		for(char y : array){x+=y;}
		return x;
	}

	/**
	 * Computes sum of all elements in short array.
	 * @param array Array to sum
	 * @return Sum of all elements, or 0 if array is null
	 */
	public static long sum(final short[] array){
		if(array==null){return 0;}
		//		if(Shared.SIMD && array.length>=SMINLEN) {return SIMD.sum(array, 0, array.length-1);}
		long x=0;
		for(short y : array){x+=y;}
		return x;
	}

	/**
	 * Computes sum of all elements in int array.
	 * @param array Array to sum
	 * @return Sum of all elements, or 0 if array is null
	 */
	public static long sum(final int[] array){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN32){return SIMD.sum(array, 0, array.length-1);}
		long x=0;
		for(int y : array){x+=y;}
		return x;
	}

	/**
	 * Sums double elements without modifying the array.
	 * SIMD reduction order can change rounding or overflow behavior.
	 * @param array Array to sum
	 * @return Sum of all elements, or 0 if array is null
	 */
	public static double sum(final double[] array){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN64){return SIMD.sum(array, 0, array.length-1);}
		double x=0;
		for(double y : array){x+=y;}
		return x;
	}

	/**
	 * Sums all long elements with ordinary two's-complement overflow.
	 * @param array Array to sum
	 * @return Sum of all elements, or 0 if array is null
	 */
	public static long sum(final long[] array){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN64){return SIMD.sum(array, 0, array.length-1);}
		long x=0;
		for(long y : array){x+=y;}
		return x;
	}

	/**
	 * Computes sum of elements in int array within specified range.
	 *
	 * @param array Array to sum
	 * @param from Start index (inclusive)
	 * @param to End index (inclusive)
	 * @return Sum of elements in range, or 0 if array is null
	 */
	public static long sum(final int[] array, final int from, final int to){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN32){return SIMD.sum(array, from, to);}//FIXED [simd/Vector#001]: SIMD path passed (0, array.length-1), ignoring from/to, so it silently summed the WHOLE array.  SIMD.sum(int[],from,to) honors arbitrary ranges.
		long x=0;
		for(int i=from; i<=to; i++){x+=array[i];}
		return x;
	}

	/**
	 * Sums an inclusive long range with ordinary two's-complement overflow.
	 *
	 * @param array Array to sum
	 * @param from Start index (inclusive)
	 * @param to End index (inclusive)
	 * @return Sum of elements in range, or 0 if array is null
	 */
	public static long sum(final long[] array, final int from, final int to){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN64){return SIMD.sum(array, from, to);}//FIXED [simd/Vector#001]: SIMD path passed (0, array.length-1), ignoring from/to, so it silently summed the WHOLE array.  SIMD.sum(long[],from,to) honors arbitrary ranges.
		long x=0;
		for(int i=from; i<=to; i++){x+=array[i];}
		return x;
	}

	/**
	 * Finds the maximum signed integer without modifying the array.
	 * @param array Array to search; may be null, but must otherwise be nonempty
	 * @return Maximum value, or 0 if array is null
	 */
	public static final int max(final int[] array){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN32){return SIMD.max(array, 0, array.length-1);}
		int max=array[0];
		for(int i=1; i<array.length; i++){
			final int x=array[i];
			max=(x>max ? x : max);
		}
		return max;
	}

	/**
	 * Finds a float maximum using first-value scalar comparison semantics.
	 * A first-element NaN is retained; later NaNs are ignored. Equal values retain
	 * the earlier value, including the sign of a zero maximum. Input is unchanged.
	 * @param array Array to search; may be null, but must otherwise be nonempty
	 * @return Maximum value, or 0 if array is null
	 */
	public static final float max(final float[] array){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN32){return SIMD.max(array, 0, array.length-1);}
		float max=array[0];
		for(int i=1; i<array.length; i++){
			final float x=array[i];
			max=(x>max ? x : max);
		}
		return max;
	}

	/**
	 * Finds the maximum signed long without modifying the array.
	 * @param array Array to search; may be null, but must otherwise be nonempty
	 * @return Maximum value, or 0 if array is null
	 */
	public static final long max(final long[] array){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN32){return SIMD.max(array, 0, array.length-1);}
		long max=array[0];
		for(int i=1; i<array.length; i++){
			final long x=array[i];
			max=(x>max ? x : max);
		}
		return max;
	}

	/**
	 * Count mismatches or Ns.
	 * @param query Query sequence bases
	 * @param ref Reference sequence bases
	 * @param maxSubs Maximum allowed substitutions (for early termination)
	 * @param rStart Starting position in reference (0-based)
	 * @return Number of substitutions found
	 */
	public static int align(byte[] query, byte[] ref, final int maxSubs, final int rStart){
		if(query==null){return 0;}
		if(Shared.SIMD && query.length>=8){return SIMDByte256.alignSIMD(query, ref, maxSubs, rStart);}
		int subs=0;
		for(int i=0, j=rStart; i<query.length && subs<=maxSubs; i++, j++){
			subs+=(query[i]!=ref[j] || query[i]=='N') ? 1 : 0;
		}
		return subs;
	}

	/** Conservatively filters fused-reference candidates without changing either array.
	 * Reference uppercase N bases count separately as padding. Both paths permit
	 * maxNs+maxSubs padding bases so potentially valid paid clipping survives.
	 * The caller must still combine actual clipping cost and substitutions before
	 * accepting an alignment; passing this prefilter is only a necessary condition.
	 * @param query Query bytes; null returns zero
	 * @param ref Reference containing the complete query-sized interval at rStart
	 * @param maxSubs Nonnegative substitution/paid-clipping budget
	 * @param maxNs Nonnegative free-padding allowance
	 * @param rStart Nonnegative reference offset for query base zero
	 * @return Substitutions, or maxSubs+1 when the selected path exceeds a budget */
	public static int alignFused(byte[] query, byte[] ref, final int maxSubs,
		final int maxNs, final int rStart){
		if(query==null){return 0;}
		if(Shared.SIMD && query.length>=8){
			return SIMDByte256.alignFusedSIMD(query, ref, maxSubs, maxNs, rStart);
		}

		// Scalar Fallback
		// Match the SIMD necessary-condition budget; widening prevents allowance overflow.
		final long limitNs=(long)maxNs+maxSubs;
		int subs=0;
		int ns=0;
		for(int i=0, j=rStart; i<query.length; i++, j++){
			byte q=query[i];
			byte r=ref[j];
			if(r=='N'){ns++;}
			else if(q!=r){subs++;}
			if(subs>maxSubs || ns>limitNs){return maxSubs+1;}
		}
		return subs;
	}

	/** Appends every symbol position in the half-open range.
	 * Null arrays leave positions unchanged. Range bounds are unchecked.
	 * @param array Byte array to search
	 * @param from First index to inspect, inclusive
	 * @param to Stop index, exclusive
	 * @param symbol Byte to find
	 * @param positions Existing destination for matching indexes
	 * @return Final size of positions, including pre-existing entries */
	public static final int findSymbols(final byte[] array,
		final int from, final int to, final byte symbol, final IntList positions){
		if(array==null){return positions.size();}
		if(Shared.SIMD && array.length>=MINLEN8){
			return SIMDByte256.findSymbols(array, from, to, symbol, positions);
		}
		for(int i=from; i<to; i++){
			if(array[i]==symbol){positions.add(i);}
		}
		return positions.size();
	}

	/** Counts occurrences of a byte in the half-open range.
	 * Null arrays return zero. Range bounds are unchecked.
	 * @param array Byte array to search
	 * @param from First index to inspect, inclusive
	 * @param to Stop index, exclusive
	 * @param symbol Byte to count
	 * @return Number of matching positions in the requested range */
	public static final int countSymbols(final byte[] array,
		final int from, final int to, final byte symbol){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN8){
			return SIMDByte256.countSymbols(array, from, to, symbol);
		}
		int positions=0;
		for(int i=from; i<to; i++){
			if(array[i]==symbol){positions++;}
		}
		return positions;
	}

	/** Finds the last matching symbol before limit.
	 * Null buffers return -1. The search range is [0, limit); limit is unchecked.
	 * @param buffer Buffer to scan
	 * @param limit Stop index, exclusive
	 * @param symbol Byte to find
	 * @return Last matching index before limit, or -1 if absent */
	public static int findLastSymbol(final byte[] buffer, final int limit, final byte symbol){
		if(buffer==null){return -1;}
		if(Shared.SIMD && limit>=MINLEN8){
			return SIMDByte256.findLastSymbol(buffer, limit, symbol);
		}
		for(int i=limit-1; i>=0; i--){
			if(buffer[i]==symbol){
				return i;
			}
		}
		return -1;
	}

	/** Appends FASTA record-boundary newline positions.
	 * A boundary is a newline immediately followed by {@code >}. Null arrays leave
	 * positions unchanged. Range bounds are unchecked; the first tested header byte
	 * is {@code from+1} so its preceding newline can be inspected.
	 * @param array Byte array to search
	 * @param from Lower bound for the preceding newline
	 * @param to Stop index for the header marker, exclusive
	 * @param positions Existing destination for boundary newline indexes
	 * @return Final size of positions, including pre-existing entries */
	public static final int findFastaHeaders(final byte[] array,
		final int from, final int to, final IntList positions){
		if(array==null){return positions.size();}
		if(Shared.SIMD && array.length>=MINLEN8){
			return SIMDByte256.findFastaHeaders(array, from, to, positions);
		}
		for(int i=from+1; i<to; i++){
			if(array[i]==carrot && array[i-1]==slashn){
				positions.add(i-1); // Store position of \n
			}
		}
		return positions.size();
	}

	/**
	 * Convert a FASTA record (header + sequence) to a Read.
	 * @param record Byte array containing ">header\nsequence" (may span multiple lines)
	 * @param readID Numeric ID for the read
	 * @param pairnum Pair number (0 or 1)
	 * @param bb Reusable sequence scratch, cleared and populated for this record
	 * @param newlines Reusable newline offsets, cleared when the SIMD scan is used
	 * @return Read object
	 */
	public static Read fastaRecordToRead(byte[] record, long readID, int pairnum,
		ByteBuilder bb, IntList newlines){
		if(record==null || record.length==0){return null;}

		// Find end of header (first newline)
		int headerEnd=0;
		while(headerEnd<record.length && record[headerEnd]!='\n'){headerEnd++;}

		// Extract header (skip leading >)
		byte[] header=Arrays.copyOfRange(record, 1, headerEnd);

		// Extract and concatenate sequence lines
		headerEnd++;
		final int remaining=record.length-headerEnd;
		bb.clear();
		bb.ensureExtra(remaining);
		if(Shared.SIMD && remaining>=MINLEN8){
			newlines.clear();
			SIMDByte256.findSymbols(record, headerEnd, record.length, slashn, newlines);

			// Copy segments between newlines
			int start=headerEnd;
			for(int i=0; i<newlines.size(); i++){
				int nlPos=newlines.get(i);
				int len=nlPos-start;
				if(len>0){
					System.arraycopy(record, start, bb.array, bb.length, len);
					bb.length+=len;
				}
				start=nlPos+1; // Skip the newline
			}

			// Copy final segment after last newline
			int len=record.length-start;
			if(len>0){
				System.arraycopy(record, start, bb.array, bb.length, len);
				bb.length+=len;
			}
		}else{
			for(int i=headerEnd; i<record.length; i++){
				if(record[i]!='\n'){
					bb.array[bb.length++]=record[i];
				}
			}
		}
		//		assert(false) : bb.length+", "+headerEnd+", "+record.length;
		Read r=new Read(bb.toBytes(), null, new String(header, StandardCharsets.US_ASCII), readID);
		r.setPairnum(pairnum);
		return r;
	}

	/** Adds delta to every byte in place; null arrays are no-ops.
	 * @param array Mutable byte array
	 * @param delta Value to add with byte overflow semantics */
	public static void add(byte[] array, byte delta){
		if(array==null){return;}
		if(Shared.SIMD && array.length>=MINLEN8){
			SIMDByte256.add(array, delta);
			return;
		}
		for(int i=0; i<array.length; i++){array[i]+=delta;}
	}

	/** Adds delta, stores each byte as {@code (byte)Math.max(cap, sum)},
	 * and reports the smallest int sum with an initial cap of 127.
	 * Stored bytes can wrap if the capped int sum is outside the byte range.
	 * Null arrays return zero; empty arrays return 127.
	 * @param array Mutable byte array
	 * @param delta Value to add before capping
	 * @param cap Minimum int value before the narrowed store
	 * @return {@code (byte)Math.min(127, min(array[i]+delta))} */
	public static byte addAndCapMin(final byte[] array, final byte delta, final int cap){
		if(array==null){return 0;}
		if(Shared.SIMD && array.length>=MINLEN8){
			return SIMDByte256.addAndCapMin(array, delta, cap);
		}
		int min=127;
		for(int i=0; i<array.length; i++){
			int b=array[i]+delta;
			min=Math.min(min, b);
			array[i]=(byte)Math.max(cap, b);
		}
		return (byte)min;
	}

	/** Applies a quality offset in place for supported uppercase A/C/G/T/N bases.
	 * A/C/G/T qualities are offset and floored at two; uppercase N receives zero.
	 * Vector lanes wrap byte addition before flooring, so A/C/G/T inputs must keep
	 * {@code quals[i]+delta} within signed-byte range for scalar/SIMD equivalence.
	 * Mixed-case, IUPAC and non-N undefined bases retain the legacy divergence:
	 * vector blocks zero only uppercase N while scalar tails zero all undefined bases.
	 * Base bytes must be nonnegative table indexes. Null quality arrays are no-ops.
	 * @param quals Mutable quality array
	 * @param bases Base array with at least quals.length entries
	 * @param delta Offset to add to each quality */
	public static void applyQualOffset(final byte[] quals, final byte[] bases, final int delta){
		if(quals==null){return;}
		if(Shared.SIMD && quals.length>=MINLEN8){
			SIMDByte256.applyQualOffset(quals, bases, delta);
			return;
		}
		for(int i=0; i<quals.length; i++){
			byte b=bases[i];
			int q=quals[i]+delta;
			q=(AminoAcid.baseToNumber[b]<0 ? 0 : Math.max(2, q));
			quals[i]=(byte)q;
		}
	}

	/** Restricts supported uppercase A/C/G/T/N qualities to the called-base range in place.
	 * A/C/G/T qualities are capped to 2..Read.MAX_CALLED_QUALITY; uppercase N receives zero.
	 * Mixed-case, IUPAC and non-N undefined bases retain the legacy divergence:
	 * vector blocks zero only uppercase N while scalar tails zero all undefined bases.
	 * Base bytes must be nonnegative table indexes. Null quality arrays are no-ops.
	 * @param quals Mutable quality array
	 * @param bases Base array with at least quals.length entries */
	public static void capQuality(final byte[] quals, final byte[] bases){
		if(quals==null){return;}
		final byte max=Read.MAX_CALLED_QUALITY;
		if(Shared.SIMD && quals.length>=MINLEN8){
			SIMDByte256.capQuality(quals, bases, max);
			return;
		}
		for(int i=0; i<quals.length; i++){
			byte b=bases[i];
			int q=quals[i];
			q=(AminoAcid.baseToNumber[b]<0 ? 0 : Math.min(Math.max(2, q), max));
			quals[i]=(byte)q;
		}
	}

	/** Converts U/u to T/t in place; null arrays are no-ops.
	 * @param bases Mutable nonnegative table-index base array */
	public static void uToT(byte[] bases){
		if(bases==null){return;}
		if(Shared.SIMD && bases.length>=MINLEN8){
			SIMDByte256.uToT(bases);
			return;
		}
		for(int i=0; i<bases.length; i++){
			bases[i]=AminoAcid.uToT[bases[i]];
		}
	}

	/** Converts '.', '-', X, x, and n to N in place.
	 * Bytes must be nonnegative table indexes. Null arrays are accepted. The current
	 * implementation always reports success.
	 * @param array Mutable byte array
	 * @return true */
	public static boolean dotDashXToN(byte[] array){
		if(array==null){return true;}
		if(Shared.SIMD && array.length>=MINLEN8){
			return SIMDByte256.dotDashXToN(array);
		}
		for(int i=0; i<array.length; i++){
			array[i]=AminoAcid.dotDashXToNocall[array[i]];
		}
		return true;
	}

	/** Checks whether the array contains an amino-acid byte that is not an IUPAC base.
	 * Bytes must be nonnegative table indexes. Null arrays return true.
	 * @param array Byte array to inspect
	 * @return true if any byte is amino-only */
	public static boolean isProtein(byte[] array){
		if(array==null){return true;}
		if(Shared.SIMD && array.length>=MINLEN8){
			return SIMDByte256.isProtein(array);
		}
		boolean protein=false;
		for(int i=0; i<array.length; i++){
			byte b=array[i];
			boolean nuc=AminoAcid.baseToNumberExtended[b]>=0;
			boolean amino=AminoAcid.acidToNumberExtended[b]>=0;
			protein|=(amino && !nuc);
		}
		return protein;
	}

	/** Converts lowercase ASCII letters to N in place.
	 * Bytes must be nonnegative table indexes. Null arrays are accepted. The current
	 * implementation always reports success.
	 * @param array Mutable byte array
	 * @return true */
	public static boolean lowerCaseToN(byte[] array){
		if(array==null){return true;}
		if(Shared.SIMD && array.length>=MINLEN8){
			return SIMDByte256.lowerCaseToN(array);
		}
		for(int i=0; i<array.length; i++){
			array[i]=AminoAcid.lowerCaseToNocall[array[i]];
		}
		return true;
	}

	/** Uppercases lowercase ASCII letters in place and leaves other bytes unchanged.
	 * Bytes must be nonnegative table indexes. Null arrays are accepted. The current
	 * implementation always reports success.
	 * @param array Mutable byte array
	 * @return true */
	public static boolean toUpperCase(byte[] array){
		if(array==null){return true;}
		if(Shared.SIMD && array.length>=MINLEN8){
			return SIMDByte256.toUpperCase(array);
		}
		boolean success=true;
		for(int i=0; i<array.length; i++){
			array[i]=AminoAcid.toUpperCase[array[i]];
		}
		return success;
	}

	/** Checks whether every byte is an ASCII letter.
	 * Null arrays return true.
	 * @param array Byte array to inspect
	 * @return true if every inspected byte is A-Z or a-z */
	public static boolean allLetters(byte[] array){
		if(array==null){return true;}
		if(Shared.SIMD && array.length>=MINLEN8){
			return SIMDByte256.allLetters(array);
		}
		final byte A='A', Z='Z';
		final byte mask=~32;
		boolean success=true;
		for(int i=0; i<array.length; i++){//Could do lookup instead
			int b=(array[i]&mask);
			success&=(b>=A && b<=Z);
		}
		return success;
	}

	/** Masks IUPAC ambiguity letters (B/D/H/K/M/N/R/S/V/W/X/Y in either case) to N.
	 * Called bases A/C/G/T/U retain their case; other letters and nonletters are
	 * unchanged. Uppercasing, U-to-T conversion and junk repair are separate policies.
	 * Null arrays are accepted without mutation.
	 * @param array Mutable nonnegative table-index base bytes
	 * @return true on normal completion */
	public static boolean iupacToN(byte[] array){
		if(array==null){return true;}
		if(Shared.SIMD && array.length>=MINLEN8){
			return SIMDByte256.iupacToN(array);
		}
		for(int i=0; i<array.length; i++){
			array[i]=AminoAcid.iupacToNocall[array[i]];
		}
		return true;
	}

	/** Checks whether every byte is an extended IUPAC nucleotide base.
	 * Bytes must be nonnegative table indexes. Null arrays return true.
	 * @param array Byte array to inspect
	 * @return true if every byte maps through AminoAcid.baseToNumberExtended */
	public static boolean isNucleotide(byte[] array){
		if(array==null){return true;}
		if(Shared.SIMD && array.length>=MINLEN8){
			return SIMDByte256.isNucleotide(array);
		}
		boolean success=true;
		for(int i=0; i<array.length; i++){
			success&=(AminoAcid.baseToNumberExtended[array[i]]>=0);;
		}
		return success;
	}

	/** Adds delta to each quality byte and appends the results.
	 * Null quality arrays leave bb unchanged.
	 * @param quals Quality bytes to append
	 * @param bb Destination builder; length advances by quals.length
	 * @param delta Value to add to each appended byte */
	public static void addAndAppend(byte[] quals, ByteBuilder bb, int delta){
		if(quals==null){return;}
		if(Shared.SIMD && quals.length>=MINLEN8){
			SIMDByte256.addAndAppend(quals, bb, delta);
			return;
		}
		final int qlen=quals.length;
		bb.ensureExtra(qlen);
		final byte[] array=bb.array;
		for(int i=0, j=bb.length; i<qlen; i++, j++){
			array[j]=(byte)(quals[i]+delta);
		}
		bb.length+=qlen;
	}

	/** Adds delta while appending quality bytes in reverse order.
	 * Null quality arrays leave bb unchanged.
	 * @param quals Quality bytes to append from end to start
	 * @param bb Destination builder; length advances by quals.length
	 * @param delta Value to add to each appended byte */
	public static void addAndAppendReversed(byte[] quals, ByteBuilder bb, int delta){
		if(quals==null){return;}
		if(Shared.SIMD && quals.length>=MINLEN8){
			SIMDByte256.addAndAppendReversed(quals, bb, delta);
			return;
		}
		final int qlen=quals.length;
		bb.ensureExtra(qlen);
		final byte[] array=bb.array;
		for(int i=qlen-1, j=bb.length; i>=0; i--, j++){
			array[j]=(byte)(quals[i]+delta);
		}
		bb.length+=qlen;
	}

	/** Appends synthetic qualities selected from base definition status.
	 * Fully defined bases append qFake; all others append qUndef. Base bytes must be
	 * nonnegative table indexes. Null base arrays leave bb unchanged.
	 * @param bases Base bytes to classify
	 * @param bb Destination builder; length advances by bases.length
	 * @param qFake Quality byte for fully defined bases
	 * @param qUndef Quality byte for ambiguous or undefined bases */
	public static void appendFake(byte[] bases, ByteBuilder bb, int qFake, int qUndef){
		if(bases==null){return;}
		if(Shared.SIMD && bases.length>=MINLEN8){
			SIMDByte256.appendFake(bases, bb, qFake, qUndef);
			return;
		}
		final int blen=bases.length;
		bb.ensureExtra(blen);
		final byte[] array=bb.array;
		for(int i=0, j=bb.length; i<blen; i++, j++){
			array[j]=(byte)(AminoAcid.isFullyDefined(bases[i]) ? qFake : qUndef);
		}
		bb.length+=blen;
	}

	/** Appends the low byte of each String UTF-16 code unit to a ByteBuilder.
	 * @param bb Destination builder
	 * @param s Non-null String to append
	 * @return bb */
	public static ByteBuilder append(ByteBuilder bb, String s){
		if(varHandles){
			return VarHandler.appendString(bb, s);
		}else{
			final int len=s.length();
			bb.expand(len);
			for(int i=0, j=bb.length; i<len; i++, j++){
				bb.array[j]=(byte)s.charAt(i);
			}
			bb.length+=len;
			return bb;
		}
	}

	/** Reverses the whole array in place; null arrays are no-ops.
	 * @param array Mutable byte array */
	public static void reverseInPlace(final byte[] array){
		if(array==null){return;}
		if(Shared.SIMD && array.length>=MINLEN8){
			SIMDByte256.reverseInPlace(array, array.length);
			return;
		}
		final int max=array.length/2, last=array.length-1;
		for(int i=0; i<max; i++){
			byte temp=array[i];
			array[i]=array[last-i];
			array[last-i]=temp;
		}
	}

	/** Reverses a prefix in place; null arrays are no-ops.
	 * Elements at or beyond length are left unchanged. Length is unchecked.
	 * @param array Mutable byte array
	 * @param length Prefix length */
	public static void reverseInPlace(final byte[] array, final int length){
		if(array==null){return;}
		if(Shared.SIMD && array.length>=MINLEN8){
			SIMDByte256.reverseInPlace(array, length);
			return;
		}
		final int max=length/2, last=length-1;
		for(int i=0; i<max; i++){
			byte temp=array[i];
			array[i]=array[last-i];
			array[last-i]=temp;
		}
	}

	/** Reverse-complements an entire nucleotide array in place.
	 * @param bases Mutable base array */
	public static void reverseComplementInPlace(final byte[] bases){
		reverseComplementInPlace(bases, bases.length);
	}

	/** Reverse-complements a prefix using the extended IUPAC complement mapping.
	 * Elements at or beyond length remain unchanged; null input is a no-op.
	 * @param bases Mutable table-indexable base bytes
	 * @param length Prefix length from zero through bases.length */
	public static void reverseComplementInPlace(final byte[] bases, final int length){
		if(true){reverseComplementInPlaceFast(bases, length); return;}//Superseded by reverseComplementInPlaceFast, which now self-detects non-ACGTN per block and falls back to scalar - yielding an identical proper complement without this method's separate up-front isACGTN scan.  Old body kept below for trivial revert.
		if(bases==null){return;}
		if(Shared.SIMD && bases.length>=MINLEN8 && SIMDByte256.isACGTN(bases, length)){
			SIMDByte256.reverseComplementInPlace(bases, length);
			return;
		}
		final int last=length-1;
		final int max=length/2;
		for(int i=0; i<max; i++){
			byte a=bases[i];
			byte b=bases[last-i];
			bases[i]=AminoAcid.baseToComplementExtended[b];
			bases[last-i]=AminoAcid.baseToComplementExtended[a];
		}
		if((length&1)==1){//Odd length; process middle
			bases[max]=AminoAcid.baseToComplementExtended[bases[max]];
		}
	}

	/** Reverse-complements an entire nucleotide array with extended IUPAC mappings.
	 * @param bases Mutable table-indexable base bytes */
	public static void reverseComplementInPlaceFast(final byte[] bases){
		reverseComplementInPlaceFast(bases, bases.length);
	}

	/** Reverse-complements a prefix using the extended IUPAC complement mapping.
	 * Elements at or beyond length remain unchanged; null input is a no-op.
	 * @param bases Mutable table-indexable base bytes
	 * @param length Prefix length from zero through bases.length */
	public static void reverseComplementInPlaceFast(final byte[] bases, final int length){
		if(bases==null){return;}
		if(Shared.SIMD && bases.length>=MINLEN8){
			SIMDByte256.reverseComplementInPlace(bases, length);
			return;
		}
		final int last=length-1;
		final int max=length/2;
		for(int i=0; i<max; i++){
			byte a=bases[i];
			byte b=bases[last-i];
			bases[i]=AminoAcid.baseToComplementExtended[b];
			bases[last-i]=AminoAcid.baseToComplementExtended[a];
		}
		if((length&1)==1){//Odd length; process middle
			bases[max]=AminoAcid.baseToComplementExtended[bases[max]];
		}
	}

	/** Probes the fixed backend's integer lane count.
	 * @return Lane count, or one if loading or probing the backend fails */
	private static synchronized int intLanes(){
		try{return SIMD.intLanes();}
		catch(Throwable e){return 1;}
	}

	/** Counts error-free sampled k-mers independently in each fixed SIMD lane.
	 * The caller supplies a scalar implementation when the backend is unavailable.
	 * @param errorBuffer Zero/one error indicators, interleaved by position then lane;
	 * must contain queryLen times INT_LANES entries
	 * @param results At least INT_LANES output entries, overwritten with lane counts
	 * @param k K-mer length from one through32
	 * @param queryLen Nonnegative number of query positions
	 * @param step Positive power-of-two spacing between sampled window starts
	 * @param wildcardMask Bits requiring error-free positions; bit zero is the oldest base
	 * @throws UnsupportedOperationException If SIMD is disabled or unavailable */
	public static void countErrorFreeKmersBatch(int[] errorBuffer, int[] results, int k, int queryLen, int step, int wildcardMask){
		if(Shared.SIMD && INT_LANES>1){
			SIMD.countErrorFreeKmersBatch(errorBuffer, results, k, queryLen, step, wildcardMask);
		}else{
			// Scalar fallback logic not implemented here; calling class should check INT_LANES
			throw new UnsupportedOperationException("Scalar fallback for countErrorFreeKmersBatch should be handled by caller");
		}
	}

	/** Finds a key in a coherent linear-probing table.
	 * Probes from initial to the end, then wraps to zero at most once. Stored keys
	 * must precede the first empty cell in their probe chains; SIMD checks blocks.
	 * @param keys Nonempty table that remains stable during the search
	 * @param key Key to find; the empty marker is always reported absent
	 * @param initial Initial cell in the range zero through keys.length-1
	 * @param invalid Empty-cell marker
	 * @return Matching cell index, or -1 when absent */
	public static int findKey(final int[] keys, final int key, final int initial, final int invalid){
		assert(keys!=null);
		if(key==invalid){return -1;}
		if(Shared.SIMD && keys.length>=MINLEN32){
			return SIMD.findKey(keys, key, initial, invalid);
		}
		final int limit=keys.length;
		int cell=initial;
		for(; cell<limit && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		if(cell<limit){return keys[cell]==key ? cell : -1;}
		cell=0;
		for(; cell<initial && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		return keys[cell]==key ? cell : -1;
	}

	/** Finds an existing key or the first empty cell in wrapped probe order.
	 * The table must contain the key or an empty cell. Exhaustion returns -1 in
	 * the SIMD backend; the scalar path first fails an assertion when enabled.
	 * @param keys Nonempty table that remains stable during the search
	 * @param key Key distinct from the empty marker
	 * @param initial Initial cell in the range zero through keys.length-1
	 * @param invalid Empty-cell marker
	 * @return Matching or empty cell index, or -1 on untrapped exhaustion */
	public static int findKeyOrInvalid(final int[] keys, final int key, final int initial, final int invalid){
		assert(keys!=null);
		assert(key!=invalid) : "Collision - this should have been intercepted.";
		if(Shared.SIMD && keys.length>=MINLEN32){
			return SIMD.findKeyOrInvalid(keys, key, initial, invalid);
		}
		final int limit=keys.length;
		int cell=initial;
		for(; cell<limit && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		if(cell<limit){return cell;}
		cell=0;
		for(; cell<initial && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		assert(cell<initial) : "Overflow at size "+limit+": key="+key+", initial="+initial+", invalid="+invalid;
		return cell<initial ? cell : -1;
	}

	/** Searches a wrapped linear-probing table without SIMD dispatch.
	 * Unlike findKey, this helper does not reject key==invalid before searching.
	 * @param keys Nonempty stable table
	 * @param key Key to find; callers normally exclude the empty marker
	 * @param initial Initial cell in the range zero through keys.length-1
	 * @param invalid Empty-cell marker that terminates an unsuccessful search
	 * @return Matching cell index, or -1 when absent */
	public static int findKeyScalar(final int[] keys, final int key, final int initial, final int invalid){
		final int limit=keys.length;
		int cell=initial;
		for(; cell<limit && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		if(cell<limit){return keys[cell]==key ? cell : -1;}
		cell=0;
		for(; cell<initial && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		return keys[cell]==key ? cell : -1;
	}

	/** Searches for a matching or empty cell without SIMD dispatch.
	 * Exhausting a full table fails an assertion, or returns -1 without assertions.
	 * @param keys Nonempty stable table containing the key or an empty cell
	 * @param key Key to find
	 * @param initial Initial cell in the range zero through keys.length-1
	 * @param invalid Empty-cell marker
	 * @return Matching or empty cell index */
	public static int findKeyOrInvalidScalar(final int[] keys, final int key, final int initial, final int invalid){
		final int limit=keys.length;
		int cell=initial;
		for(; cell<limit && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		if(cell<limit){return cell;}
		cell=0;
		for(; cell<initial && keys[cell]!=key && keys[cell]!=invalid; cell++){}
		assert(cell<initial) : "Overflow at size "+limit+": key="+key+", initial="+initial+", invalid="+invalid;
		return cell<initial ? cell : -1;
	}

	/** Replaces every exact occurrence of one integer value in place.
	 * This does not rehash a table; callers must preserve their own probe invariants.
	 * @param keys Mutable array
	 * @param oldKey Value to replace
	 * @param newKey Replacement value */
	public static void changeAll(final int[] keys, final int oldKey, final int newKey){
		assert(keys!=null);
		if(Shared.SIMD && keys.length>=MINLEN32){
			SIMD.changeAll(keys, oldKey, newKey);
			return;
		}
		for(int i=0; i<keys.length; i++){
			if(keys[i]==oldKey){keys[i]=newKey;}
		}
	}

	/** Checks whether the optional Vector API class is visible to this loader.
	 * @return False when ByteVector cannot be found; other linkage errors propagate */
	private static synchronized boolean vectorLoaded(){
		try{Class.forName("jdk.incubator.vector.ByteVector");}
		catch(ClassNotFoundException e){return false;}
		return true;
	}

	/** Asks the optional backend for its preferred vector width.
	 * @return Reported width in bits, or zero if backend initialization/probing fails */
	private static synchronized int maxSimdWidth(){
		try{return SIMD.maxVectorLength();}
		catch(Throwable e){return 0;}
	}

	/** Probes the optional String backend without consulting partially initialized Shared fields.
	 * @return True when the backend links and its private String lookup succeeds */
	private static boolean varHandlesAvailable(){
		try{return VarHandler.AVAILABLE;}
		catch(LinkageError e){return false;}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Static Fields        ----------------*/
	/*--------------------------------------------------------------*/

	/** Minimum array length for 8-bit SIMD operations */
	public static final int MINLEN8=8;//Due to dual SIMD
	/** Minimum array length for 16-bit SIMD operations */
	public static final int MINLEN16=16;
	/** Minimum array length for 32-bit SIMD operations (optimized for 16) */
	public static final int MINLEN32=16;//16 or 32 are optimal; 0, 24, and 48 are worse.
	/** Minimum array length for 64-bit SIMD operations */
	public static final int MINLEN64=8;
	/**
	 * Whether to use SIMD for copy operations (disabled due to power usage concerns)
	 */
	public static boolean SIMDCOPY=false;//Does not seem to affect speed, but could increase power usage.
	/**
	 * Whether to use SIMD for sparse multiplication (grants speedup, currently broken at ebs=1)
	 */
	public static boolean SIMD_MULT_SPARSE=true;//Grants a speedup, and same results (but currently broken at ebs=1)
	/**
	 * Whether to use SIMD for sparse FMA operations (grants speedup, slightly different results)
	 */
	public static boolean SIMD_FMA_SPARSE=true;//Grants a speedup, slightly different results
	/** Whether the optional Vector API probe can load in this JVM. */
	private static final boolean vectorLoaded=vectorLoaded();
	/** Maximum vector width reported by the probe, or zero when unavailable. */
	public static final int maxSimdWidth=vectorLoaded ? maxSimdWidth() : 0;
	/** Whether the probed width supports the fixed 256-bit kernels. */
	public static final boolean simd256=maxSimdWidth>=256;
	/** Whether the probed width is at least 128 bits; not a backend selector. */
	public static final boolean simd128=maxSimdWidth>=128;
	/** Whether the probed width is at least 64 bits; not a backend selector. */
	public static final boolean simd64=maxSimdWidth>=64;
	/** Fixed-backend integer lanes, or one when unavailable; independent of the runtime SIMD switch. */
	public static final int INT_LANES=simd256 ? intLanes() : 1;
	/** Whether the optional backend can access String storage; false when unavailable or unlinked. */
	public static final boolean varHandles=varHandlesAvailable();
	/** Carriage-return delimiter used by record scanning. */
	private static final byte slashr='\r';
	/** Line-feed delimiter used by record scanning. */
	private static final byte slashn='\n';
	/** FASTA record marker. */
	private static final byte carrot='>';
	/** FASTQ quality-header marker. */
	private static final byte plus='+';
	/** FASTQ record marker. */
	private static final byte at='@';

}
