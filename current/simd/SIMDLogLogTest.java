package simd;

import java.util.Arrays;
import cardinality.DynamicDemiLog;

/**
 * Correctness test for SIMDLogLog.compareDetailed().
 * Verifies SIMD results match scalar baseline across:
 * - Random DDLs with realistic value distributions
 * - Edge cases (all-zero, identical, one-empty)
 * - High char values (> 32767) that require unsigned comparison
 *
 * @author Brian Bushnell, Ady
 * @date April 2026
 */
public class SIMDLogLogTest {

	public static void main(String[] args){
		int pass=0, fail=0;

		//Test 1: Random DDLs from actual hashing
		{
			DynamicDemiLog a=DynamicDemiLog.create(2048, 31, 12345L, 0f, true);
			DynamicDemiLog b=DynamicDemiLog.create(2048, 31, 12345L, 0f, true);
			//Hash some fake sequences
			byte[] seqA=makeRandomBases(500000, 42);
			byte[] seqB=makeRandomBases(500000, 99);
			a.hash(seqA, null);
			b.hash(seqB, null);

			int[] scalar=SIMDLogLog.compareDetailedScalar(a.maxArray(), b.maxArray());
			int[] simd=SIMDLogLog.compareDetailed(a.maxArray(), b.maxArray());
			if(check("Random DDLs", scalar, simd)){pass++;}else{fail++;}
		}

		//Test 2: Identical DDLs
		{
			DynamicDemiLog a=DynamicDemiLog.create(2048, 31, 12345L, 0f, true);
			byte[] seq=makeRandomBases(200000, 7);
			a.hash(seq, null);

			int[] scalar=SIMDLogLog.compareDetailedScalar(a.maxArray(), a.maxArray());
			int[] simd=SIMDLogLog.compareDetailed(a.maxArray(), a.maxArray());
			if(check("Identical DDLs", scalar, simd)){pass++;}else{fail++;}
		}

		//Test 3: All-zero arrays
		{
			char[] z=new char[2048];
			int[] scalar=SIMDLogLog.compareDetailedScalar(z, z);
			int[] simd=SIMDLogLog.compareDetailed(z, z);
			if(check("All-zero", scalar, simd)){pass++;}else{fail++;}
		}

		//Test 4: One populated, one empty
		{
			DynamicDemiLog a=DynamicDemiLog.create(2048, 31, 12345L, 0f, true);
			byte[] seq=makeRandomBases(200000, 3);
			a.hash(seq, null);
			char[] z=new char[2048];

			int[] scalar=SIMDLogLog.compareDetailedScalar(a.maxArray(), z);
			int[] simd=SIMDLogLog.compareDetailed(a.maxArray(), z);
			if(check("One empty", scalar, simd)){pass++;}else{fail++;}
		}

		//Test 5: High values (> 32767) — validates unsigned comparison
		{
			char[] hi=new char[2048];
			char[] lo=new char[2048];
			for(int i=0; i<2048; i++){
				hi[i]=(char)(40000+(i%100)); //40000-40099, all > 32767
				lo[i]=(char)(100+(i%100));   //100-199, all < 32767
			}
			int[] scalar=SIMDLogLog.compareDetailedScalar(hi, lo);
			int[] simd=SIMDLogLog.compareDetailed(hi, lo);
			if(check("High values (>32767)", scalar, simd)){pass++;}else{fail++;}
		}

		//Test 6: Mixed high/low values
		{
			char[] a=new char[2048];
			char[] b=new char[2048];
			java.util.Random rng=new java.util.Random(777);
			for(int i=0; i<2048; i++){
				a[i]=(char)rng.nextInt(65536);
				b[i]=(char)rng.nextInt(65536);
			}
			int[] scalar=SIMDLogLog.compareDetailedScalar(a, b);
			int[] simd=SIMDLogLog.compareDetailed(a, b);
			if(check("Mixed high/low", scalar, simd)){pass++;}else{fail++;}
		}

		//Test 7: Performance benchmark
		{
			DynamicDemiLog a=DynamicDemiLog.create(2048, 31, 12345L, 0f, true);
			DynamicDemiLog b=DynamicDemiLog.create(2048, 31, 12345L, 0f, true);
			byte[] seqA=makeRandomBases(2000000, 42);
			byte[] seqB=makeRandomBases(2000000, 99);
			a.hash(seqA, null);
			b.hash(seqB, null);

			short[] sa=SIMDLogLog.toShortArray(a.maxArray());
			short[] sb=SIMDLogLog.toShortArray(b.maxArray());

			//Warmup
			for(int i=0; i<10000; i++){
				SIMDLogLog.compareDetailedScalar(a.maxArray(), b.maxArray());
				SIMDLogLog.compareDetailedShort(sa, sb, sa.length);
			}

			int iters=1000000;
			long t0=System.nanoTime();
			for(int i=0; i<iters; i++){
				SIMDLogLog.compareDetailedScalar(a.maxArray(), b.maxArray());
			}
			long scalarNs=System.nanoTime()-t0;

			t0=System.nanoTime();
			for(int i=0; i<iters; i++){
				SIMDLogLog.compareDetailedShort(sa, sb, sa.length);
			}
			long simdNs=System.nanoTime()-t0;

			System.err.println("\nPerformance ("+iters+" iterations):");
			System.err.println("  Scalar: "+String.format("%.1f", scalarNs/1e6)+" ms"
				+" ("+String.format("%.0f", (double)scalarNs/iters)+" ns/op)");
			System.err.println("  SIMD:   "+String.format("%.1f", simdNs/1e6)+" ms"
				+" ("+String.format("%.0f", (double)simdNs/iters)+" ns/op)");
			System.err.println("  Speedup: "+String.format("%.1fx", (double)scalarNs/simdNs));
		}

		System.err.println("\n"+pass+" passed, "+fail+" failed.");
		if(fail>0){System.exit(1);}
	}

	static boolean check(String name, int[] scalar, int[] simd){
		boolean ok=Arrays.equals(scalar, simd);
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

	static byte[] makeRandomBases(int len, long seed){
		byte[] bases=new byte[len];
		java.util.Random rng=new java.util.Random(seed);
		byte[] acgt={'A','C','G','T'};
		for(int i=0; i<len; i++){
			bases[i]=acgt[rng.nextInt(4)];
		}
		return bases;
	}
}
