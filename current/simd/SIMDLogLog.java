package simd;

import java.util.Arrays;
import java.util.Random;

import jdk.incubator.vector.ShortVector;
import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/** 256-bit vectorized DynamicDemiLog bucket comparison.
 * Compares equal-length bucket arrays, counting lower, equal and higher
 * nonempty buckets while also counting skipped empty-empty pairs. Array lengths
 * divisible by {@link #SWIDTH} have no scalar tail in the 16-bit paths; that is
 * multiples of 16 for the fixed 256-bit short-vector species.
 * @author Brian Bushnell, Ady
 * @date April 2026
 */
public class SIMDLogLog{

	/**
	 * Runs randomized parity checks against the scalar implementation.
	 * @param args unused
	 */
	public static void main(final String[] args){
		final int trials=1000;
		final int len=2048;
		final Random rand=new Random(42);
		int passed=0, failed=0;

		for(int t=0; t<trials; t++){
			final char[] a=new char[len];
			final char[] b=new char[len];
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

			final int[] simd=compareDetailed(a, b);
			final int[] scalar=compareDetailedScalar(a, b);
			final short[] sa=charToShort(a), sb=charToShort(b);
			final int[] fromShort=compareDetailedShort(sa, sb, len);
			final int[] hybrid=compareDetailedShort(sa, b, len);

			final boolean ok=Arrays.equals(simd, scalar) && Arrays.equals(simd, fromShort)
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

			final int sum=simd[0]+simd[1]+simd[2]+simd[3];
			if(sum!=len){
				System.err.println("FAIL trial "+t+": sum="+sum+" != len="+len);
				failed++;
			}
		}

		System.err.println("SIMDLogLog fromCharArray verification: "+
			passed+"/"+trials+" passed, "+failed+" failed.");
		if(failed==0){System.err.println("PASS");}
		else{throw new AssertionError("SIMDLogLog failures="+failed);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	/** Allows construction for historical public no-arg compatibility. */
	public SIMDLogLog(){}

	/*--------------------------------------------------------------*/
	/*----------------      16-bit Comparison       ----------------*/
	/*--------------------------------------------------------------*/

	/** Vectorized three-way comparison of DDL maxArray char arrays.
	 * Excludes empty-empty pairs from the lower/equal/higher counts. Returns a
	 * new {@code int[4]} and loads char arrays directly without a full-array
	 * conversion.
	 * @param arrayA First DDL maxArray
	 * @param arrayB Second DDL maxArray; must have the same length as arrayA
	 * @return int[]{lower, equal, higher, bothEmpty} */
	public static int[] compareDetailed(final char[] arrayA, final char[] arrayB){
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

			//Comprehension: lt/eq/gt partition each ACTIVE lane exactly once (so lower+equal+higher+bothEmpty==lanes, checked by
			//main()'s sum==len). va==vb => vaU==vbU => lt&gt both false, eq true; va<vb(unsigned) => vaU<vbU => only lt.
			//The XOR bias is required for lt/gt (unsigned char order) but NOT for eq (equality is bias-invariant),
			//which is why eqMask uses raw va/vb.
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
	 * Use toShortArray() to convert once, then compare many times.
	 * @param sa First array, with char bits reinterpreted as short
	 * @param sb Second array, with char bits reinterpreted as short
	 * @param len Number of leading elements to compare; must be between zero
	 * and both array lengths, inclusive
	 * @return int[]{lower, equal, higher, bothEmpty} */
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

		//Scalar tail for lengths not divisible by SWIDTH.
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
	 * Loads query via fromArray and ref via fromCharArray. Pre-convert the query
	 * once, then compare against many char[] refs.
	 * @param sa Query array pre-converted via toShortArray
	 * @param arrayB Reference DDL maxArray in native char[] form
	 * @param len Number of leading elements to compare; must be between zero
	 * and both array lengths, inclusive
	 * @return int[]{lower, equal, higher, bothEmpty} */
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

	/** Vectorized three-way signed comparison for 32-bit DDL arrays.
	 * @param arrayA First 32-bit DDL array
	 * @param arrayB Second 32-bit DDL array; must have the same length as arrayA
	 * @return int[]{lower, equal, higher, bothEmpty} */
	public static int[] compareDetailed(final int[] arrayA, final int[] arrayB){
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

			//Int lanes use signed ordering.
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
	 * Same bit pattern; only the signedness interpretation changes.
	 * Convert once, compare many times to avoid repeated copies.
	 * @param arr char array, such as a DDL maxArray
	 * @return short array with identical bit patterns */
	public static short[] toShortArray(final char[] arr){
		return charToShort(arr);
	}

	/** Internal conversion from char[] to short[].
	 * @param arr char array to convert
	 * @return short array with identical bit patterns
	 */
	private static short[] charToShort(final char[] arr){
		final short[] out=new short[arr.length];
		for(int i=0; i<arr.length; i++){
			out[i]=(short)arr[i];
		}
		return out;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Scalar Fallback       ----------------*/
	/*--------------------------------------------------------------*/

	/** Scalar oracle for verifying SIMD result parity.
	 * This class still runs Vector API static initializers, so this method is
	 * not a module-absent fallback.
	 * @param arrayA First DDL maxArray
	 * @param arrayB Second DDL maxArray; must have the same length as arrayA
	 * @return int[]{lower, equal, higher, bothEmpty} */
	public static int[] compareDetailedScalar(final char[] arrayA, final char[] arrayB){
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
	/*----------------          Constants           ----------------*/
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

}
