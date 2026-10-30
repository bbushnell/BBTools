package test.simd;

import java.util.Arrays;

import cardinality.DynamicDemiLog;
import shared.Shared;
import simd.Vector;

/** Correctness test for DDL public facade comparisons over realistic and edge inputs.
 * Verifies SIMD results match scalar baseline across random DynamicDemiLog arrays,
 * empty arrays, identical arrays, populated/empty pairs, and unsigned char values.
 * @author Brian Bushnell, Ady, Sucrose
 * @date April 2026
 */
public final class SIMDLogLogTest{

	/** Runs finite DDL correctness checks.
	 * @param args Unused
	 */
	public static void main(final String[] args){
		if(args.length>0){throw new IllegalArgumentException("No arguments expected");}
		if(!Vector.simd256){throw new IllegalStateException("DDL SIMD test requires 256-bit Vector support");}
		final boolean oldSimd=Shared.SIMD;
		int pass=0, fail=0;
		try{
			//Test 1: Random DDLs from actual hashing.
			{
				final DynamicDemiLog a=DynamicDemiLog.create(2048, 31, 12345L, 0f, true);
				final DynamicDemiLog b=DynamicDemiLog.create(2048, 31, 12345L, 0f, true);
				final byte[] seqA=makeRandomBases(500000, 42);
				final byte[] seqB=makeRandomBases(500000, 99);
				a.hash(seqA, null);
				b.hash(seqB, null);
				if(checkAll("Random DDLs", a.maxArray(), b.maxArray())){pass++;}else{fail++;}
			}

			//Test 2: Identical DDLs.
			{
				final DynamicDemiLog a=DynamicDemiLog.create(2048, 31, 12345L, 0f, true);
				final byte[] seq=makeRandomBases(200000, 7);
				a.hash(seq, null);
				if(checkAll("Identical DDLs", a.maxArray(), a.maxArray())){pass++;}else{fail++;}
			}

			//Test 3: All-zero arrays.
			{
				final char[] z=new char[2048];
				if(checkAll("All-zero", z, z)){pass++;}else{fail++;}
			}

			//Test 4: One populated, one empty.
			{
				final DynamicDemiLog a=DynamicDemiLog.create(2048, 31, 12345L, 0f, true);
				final byte[] seq=makeRandomBases(200000, 3);
				a.hash(seq, null);
				final char[] z=new char[2048];
				if(checkAll("One empty", a.maxArray(), z)){pass++;}else{fail++;}
			}

			//Test 5: High values (> 32767) validate unsigned comparison.
			{
				final char[] hi=new char[2048];
				final char[] lo=new char[2048];
				for(int i=0; i<2048; i++){
					hi[i]=(char)(40000+(i%100));
					lo[i]=(char)(100+(i%100));
				}
				if(checkAll("High values (>32767)", hi, lo)){pass++;}else{fail++;}
			}

			//Test 6: Mixed high/low values.
			{
				final char[] a=new char[2048];
				final char[] b=new char[2048];
				final java.util.Random rng=new java.util.Random(777);
				for(int i=0; i<2048; i++){
					a[i]=(char)rng.nextInt(65536);
					b[i]=(char)rng.nextInt(65536);
				}
				if(checkAll("Mixed high/low", a, b)){pass++;}else{fail++;}
			}
		}finally{
			Shared.SIMD=oldSimd;
		}

		System.err.println("\n"+pass+" passed, "+fail+" failed.");
		if(fail>0){System.exit(1);}
	}

	/** Checks all public DDL facade overloads for one input pair.
	 * @param name Case name
	 * @param a First DDL maxArray
	 * @param b Second DDL maxArray
	 * @return True if the three SIMD facade overloads match scalar output
	 */
	private static boolean checkAll(final String name, final char[] a, final char[] b){
		Shared.SIMD=false;
		final int[] scalar=Vector.compareDDL(a, b);
		final short[] sa=Vector.toShortArray(a);
		final short[] sb=Vector.toShortArray(b);

		Shared.SIMD=true;
		boolean ok=true;
		ok&=check(name+" char-char", scalar, Vector.compareDDL(a, b));
		ok&=check(name+" short-char", scalar, Vector.compareDDL(sa, b));
		ok&=check(name+" short-short", scalar, Vector.compareDDL(sa, sb));
		return ok;
	}

	/** Checks one expected/observed comparison result.
	 * @param name Case name
	 * @param scalar Expected scalar counts
	 * @param simd Observed SIMD counts
	 * @return True if the four counts are identical
	 */
	private static boolean check(final String name, final int[] scalar, final int[] simd){
		final boolean ok=Arrays.equals(scalar, simd);
		System.err.println((ok ? "PASS" : "FAIL")+": "+name);
		if(!ok){
			System.err.println("  Scalar: "+Arrays.toString(scalar));
			System.err.println("  SIMD:   "+Arrays.toString(simd));
		}else{
			System.err.println("  L="+scalar[0]+" E="+scalar[1]
					+" H="+scalar[2]+" Empty="+scalar[3]);
		}
		return ok;
	}

	/** Makes a deterministic random A/C/G/T sequence.
	 * @param len Number of bases
	 * @param seed Random seed
	 * @return Newly allocated base array
	 */
	private static byte[] makeRandomBases(final int len, final long seed){
		final byte[] bases=new byte[len];
		final java.util.Random rng=new java.util.Random(seed);
		final byte[] acgt={'A','C','G','T'};
		for(int i=0; i<len; i++){
			bases[i]=acgt[rng.nextInt(4)];
		}
		return bases;
	}

	/** Prevents construction of this standalone test. */
	private SIMDLogLogTest(){}
}
