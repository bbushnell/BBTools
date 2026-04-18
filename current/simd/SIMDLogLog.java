package simd;

import jdk.incubator.vector.ShortVector;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * AVX256-vectorized DynamicDemiLog bucket comparison.
 * Compares two char[2048] maxArrays, counting lower/equal/higher buckets
 * while skipping empty-empty pairs.  Branchless, no tail loop for
 * power-of-2 array sizes (2048/16 = 128 vector ops).
 *
 * ~16x faster than scalar compareDetailed from SIMD alone;
 * ~50x faster than MinHashSketch comparison (branchless + no sort).
 *
 * @author Brian Bushnell, Ady
 * @date April 2026
 */
public class SIMDLogLog {

	/*--------------------------------------------------------------*/
	/*----------------         Species Setup        ----------------*/
	/*--------------------------------------------------------------*/

	/** 256-bit short vectors: 16 lanes of 16-bit values. */
	@SuppressWarnings("restriction")
	private static final VectorSpecies<Short> SSPECIES=ShortVector.SPECIES_256;
	/** Number of short lanes per vector (16 for 256-bit). */
	private static final int SWIDTH=SSPECIES.length();

	/** 256-bit int vectors: 8 lanes of 32-bit values. */
	@SuppressWarnings("restriction")
	private static final VectorSpecies<Integer> ISPECIES=IntVector.SPECIES_256;
	/** Number of int lanes per vector (8 for 256-bit). */
	private static final int IWIDTH=ISPECIES.length();

	/** Bias for unsigned short comparison: XOR flips sign bit,
	 *  converting unsigned ordering to signed ordering. */
	private static final ShortVector SBIAS=ShortVector.broadcast(SSPECIES, Short.MIN_VALUE);
	/** Zero vector for empty-bucket detection. */
	private static final ShortVector SZERO=ShortVector.zero(SSPECIES);

	/** Zero vector for int empty-bucket detection. */
	private static final IntVector IZERO=IntVector.zero(ISPECIES);

	/*--------------------------------------------------------------*/
	/*----------------      16-bit Comparison       ----------------*/
	/*--------------------------------------------------------------*/

	/** Vectorized three-way comparison of DDL maxArray char arrays.
	 *  Equivalent to DynamicDemiLog.compareDetailed() but ~16x faster.
	 *  Excludes empty-empty pairs (both values zero).
	 *  @param arrayA First DDL maxArray
	 *  @param arrayB Second DDL maxArray
	 *  @return int[]{lower, equal, higher, bothEmpty} */
	public static int[] compareDetailed(char[] arrayA, char[] arrayB){
		assert(arrayA.length==arrayB.length);
		final int len=arrayA.length;
		final short[] sa=charToShort(arrayA);
		final short[] sb=charToShort(arrayB);
		return compareDetailedShort(sa, sb, len);
	}

	/** Core SIMD comparison on pre-converted short arrays.
	 *  Use toShortArray() to convert once, then compare many times.
	 *  @param sa First array (char bits reinterpreted as short)
	 *  @param sb Second array (char bits reinterpreted as short)
	 *  @param len Number of elements to compare
	 *  @return int[]{lower, equal, higher, bothEmpty} */
	public static int[] compareDetailedShort(short[] sa, short[] sb, int len){
		int lower=0, equal=0, higher=0, bothEmpty=0;
		final int limit=SSPECIES.loopBound(len);

		for(int i=0; i<limit; i+=SWIDTH){
			final ShortVector va=ShortVector.fromArray(SSPECIES, sa, i);
			final ShortVector vb=ShortVector.fromArray(SSPECIES, sb, i);

			//Detect empty-empty pairs (both zero)
			final VectorMask<Short> aZero=va.compare(VectorOperators.EQ, SZERO);
			final VectorMask<Short> bZero=vb.compare(VectorOperators.EQ, SZERO);
			final VectorMask<Short> bothZero=aZero.and(bZero);
			final VectorMask<Short> active=bothZero.not();

			//Unsigned comparison via XOR bias: flipping sign bit converts
			//unsigned char ordering (0..65535) to signed short ordering
			final ShortVector vaU=va.lanewise(VectorOperators.XOR, SBIAS);
			final ShortVector vbU=vb.lanewise(VectorOperators.XOR, SBIAS);

			final VectorMask<Short> ltMask=vaU.compare(VectorOperators.LT, vbU).and(active);
			final VectorMask<Short> eqMask=va.compare(VectorOperators.EQ, vb).and(active);
			final VectorMask<Short> gtMask=vaU.compare(VectorOperators.GT, vbU).and(active);

			lower+=ltMask.trueCount();
			equal+=eqMask.trueCount();
			higher+=gtMask.trueCount();
			bothEmpty+=bothZero.trueCount();
		}

		//Scalar tail for non-power-of-2 lengths (never fires for 2048)
		for(int i=limit; i<len; i++){
			final int a=sa[i]&0xFFFF, b=sb[i]&0xFFFF;
			if(a==0 && b==0){bothEmpty++; continue;}
			if(a<b){lower++;}
			else if(a>b){higher++;}
			else{equal++;}
		}
		return new int[]{lower, equal, higher, bothEmpty};
	}

	/*--------------------------------------------------------------*/
	/*----------------      32-bit Comparison       ----------------*/
	/*--------------------------------------------------------------*/

	/** Vectorized three-way comparison for future 32-bit DDL arrays.
	 *  32-bit values are non-negative (6-bit NLZ + 24-bit mantissa + 2-bit history),
	 *  so signed int comparison is correct without bias.
	 *  @param arrayA First 32-bit DDL array
	 *  @param arrayB Second 32-bit DDL array
	 *  @return int[]{lower, equal, higher, bothEmpty} */
	public static int[] compareDetailed(int[] arrayA, int[] arrayB){
		assert(arrayA.length==arrayB.length);
		final int len=arrayA.length;
		int lower=0, equal=0, higher=0, bothEmpty=0;
		final int limit=ISPECIES.loopBound(len);

		for(int i=0; i<limit; i+=IWIDTH){
			final IntVector va=IntVector.fromArray(ISPECIES, arrayA, i);
			final IntVector vb=IntVector.fromArray(ISPECIES, arrayB, i);

			final VectorMask<Integer> aZero=va.compare(VectorOperators.EQ, IZERO);
			final VectorMask<Integer> bZero=vb.compare(VectorOperators.EQ, IZERO);
			final VectorMask<Integer> bothZero=aZero.and(bZero);
			final VectorMask<Integer> active=bothZero.not();

			//No bias needed: 32-bit DDL values are non-negative, so
			//signed comparison matches unsigned ordering
			final VectorMask<Integer> ltMask=va.compare(VectorOperators.LT, vb).and(active);
			final VectorMask<Integer> eqMask=va.compare(VectorOperators.EQ, vb).and(active);
			final VectorMask<Integer> gtMask=va.compare(VectorOperators.GT, vb).and(active);

			lower+=ltMask.trueCount();
			equal+=eqMask.trueCount();
			higher+=gtMask.trueCount();
			bothEmpty+=bothZero.trueCount();
		}

		for(int i=limit; i<len; i++){
			final int a=arrayA[i], b=arrayB[i];
			if(a==0 && b==0){bothEmpty++; continue;}
			if(a<b){lower++;}
			else if(a>b){higher++;}
			else{equal++;}
		}
		return new int[]{lower, equal, higher, bothEmpty};
	}

	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Converts char[] to short[] for SIMD loading.
	 *  Same bit pattern; only the signedness interpretation changes.
	 *  Convert once, compare many times to avoid repeated copies.
	 *  @param arr char array (e.g., DDL maxArray)
	 *  @return short array with identical bit patterns */
	public static short[] toShortArray(char[] arr){
		return charToShort(arr);
	}

	/** Internal conversion from char[] to short[]. */
	private static short[] charToShort(char[] arr){
		final short[] out=new short[arr.length];
		for(int i=0; i<arr.length; i++){
			out[i]=(short)arr[i];
		}
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Scalar Fallback       ----------------*/
	/*--------------------------------------------------------------*/

	/** Scalar compareDetailed, identical to DynamicDemiLog.compareDetailed().
	 *  For correctness verification and non-SIMD fallback.
	 *  @param arrayA First DDL maxArray
	 *  @param arrayB Second DDL maxArray
	 *  @return int[]{lower, equal, higher, bothEmpty} */
	public static int[] compareDetailedScalar(char[] arrayA, char[] arrayB){
		int lower=0, equal=0, higher=0, bothEmpty=0;
		for(int i=0; i<arrayA.length; i++){
			final int a=arrayA[i], b=arrayB[i];
			if(a==0 && b==0){bothEmpty++; continue;}
			int dif=a-b;
			int nbit=(dif>>>31);
			int hbit=((-dif)>>>31);
			int ebit=1-nbit-hbit;
			lower+=nbit;
			higher+=hbit;
			equal+=ebit;
			//No bothEmpty increment here since we continue above
		}
		return new int[]{lower, equal, higher, bothEmpty};
	}
}
