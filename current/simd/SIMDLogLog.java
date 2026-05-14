package simd;

import java.util.Arrays;
import java.util.Random;

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
	 *  Zero-allocation: loads char[] directly into ShortVector via fromCharArray.
	 *  @param arrayA First DDL maxArray
	 *  @param arrayB Second DDL maxArray
	 *  @return int[]{lower, equal, higher, bothEmpty} */
	public static int[] compareDetailed(char[] arrayA, char[] arrayB){
		assert(arrayA.length==arrayB.length);
		final int len=arrayA.length;
		int lower=0, equal=0, higher=0, bothEmpty=0;
		final int limit=SSPECIES.loopBound(len);

		for(int i=0; i<limit; i+=SWIDTH){
			final ShortVector va=ShortVector.fromCharArray(SSPECIES, arrayA, i);
			final ShortVector vb=ShortVector.fromCharArray(SSPECIES, arrayB, i);

			final VectorMask<Short> aZero=va.compare(VectorOperators.EQ, SZERO);
			final VectorMask<Short> bZero=vb.compare(VectorOperators.EQ, SZERO);
			final VectorMask<Short> bothZero=aZero.and(bZero);
			final VectorMask<Short> active=bothZero.not();

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

		for(int i=limit; i<len; i++){
			final int a=arrayA[i], b=arrayB[i];
			if(a==0 && b==0){bothEmpty++; continue;}
			if(a<b){lower++;}
			else if(a>b){higher++;}
			else{equal++;}
		}
		return new int[]{lower, equal, higher, bothEmpty};
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

	/** Hybrid SIMD comparison: pre-converted short[] query vs char[] reference.
	 *  Loads query via fromArray (fast direct load) and ref via fromCharArray.
	 *  Pre-convert the query once, then compare against many char[] refs.
	 *  @param sa Query array (pre-converted via toShortArray)
	 *  @param arrayB Reference DDL maxArray (native char[])
	 *  @param len Number of elements to compare
	 *  @return int[]{lower, equal, higher, bothEmpty} */
	public static int[] compareDetailedShort(short[] sa, char[] arrayB, int len){
		int lower=0, equal=0, higher=0, bothEmpty=0;
		final int limit=SSPECIES.loopBound(len);

		for(int i=0; i<limit; i+=SWIDTH){
			final ShortVector va=ShortVector.fromArray(SSPECIES, sa, i);
			final ShortVector vb=ShortVector.fromCharArray(SSPECIES, arrayB, i);

			final VectorMask<Short> aZero=va.compare(VectorOperators.EQ, SZERO);
			final VectorMask<Short> bZero=vb.compare(VectorOperators.EQ, SZERO);
			final VectorMask<Short> bothZero=aZero.and(bZero);
			final VectorMask<Short> active=bothZero.not();

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

		for(int i=limit; i<len; i++){
			final int a=sa[i]&0xFFFF, b=arrayB[i];
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
		}
		return new int[]{lower, equal, higher, bothEmpty};
	}

	/*--------------------------------------------------------------*/
	/*----------------        Verification          ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		final int trials=1000;
		final int len=2048;
		final Random rand=new Random(42);
		int passed=0, failed=0;

		for(int t=0; t<trials; t++){
			char[] a=new char[len];
			char[] b=new char[len];
			for(int i=0; i<len; i++){
				int r=rand.nextInt(100);
				if(r<30){a[i]=0;}
				else if(r<50){a[i]=(char)rand.nextInt(256);}
				else{a[i]=(char)rand.nextInt(65536);}

				r=rand.nextInt(100);
				if(r<30){b[i]=0;}
				else if(r<50){b[i]=(char)rand.nextInt(256);}
				else{b[i]=(char)rand.nextInt(65536);}
			}

			int[] simd=compareDetailed(a, b);
			int[] scalar=compareDetailedScalar(a, b);
			short[] sa=charToShort(a), sb=charToShort(b);
			int[] fromShort=compareDetailedShort(sa, sb, len);
			int[] hybrid=compareDetailedShort(sa, b, len);

			boolean ok=Arrays.equals(simd, scalar) && Arrays.equals(simd, fromShort)
				&& Arrays.equals(simd, hybrid);
			if(ok){
				passed++;
			}else{
				failed++;
				if(failed<=3){
					System.err.println("FAIL trial "+t+
						": simd="+Arrays.toString(simd)+
						" scalar="+Arrays.toString(scalar)+
						" short="+Arrays.toString(fromShort));
				}
			}

			int sum=simd[0]+simd[1]+simd[2]+simd[3];
			if(sum!=len){
				System.err.println("FAIL trial "+t+": sum="+sum+" != len="+len);
				failed++;
			}
		}

		System.err.println("SIMDLogLog fromCharArray verification: "+
			passed+"/"+trials+" passed, "+failed+" failed.");
		if(failed==0){System.err.println("PASS");}
	}
}
