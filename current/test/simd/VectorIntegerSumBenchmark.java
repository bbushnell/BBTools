package test.simd;

import shared.Shared;
import simd.Vector;

/** Measures repeated integer sums with a retained, independently checked result.
 * Mutates one array element per call so repeated results cannot be hoisted.
 * Fixed-size small values keep both the historical and repaired sums in range.
 * @author Shinobu */
public final class VectorIntegerSumBenchmark{

	/** Runs warmup and a fixed amount of timed single-threaded work.
	 * @param args loops=positive integer and simd=t|f
	 */
	public static void main(final String[] args){
		int loops=65536;
		boolean simd=false;
		for(String arg : args){
			if(arg.startsWith("loops=")){loops=Integer.parseInt(arg.substring(6));}
			else if(arg.equals("simd=t")){simd=true;}
			else if(arg.equals("simd=f")){simd=false;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(loops<1){throw new IllegalArgumentException("loops must be positive: "+loops);}
		if(simd && !Vector.simd256){
			throw new IllegalStateException("SIMD benchmark requires 256-bit Vector support");
		}
		final boolean oldSimd=Shared.SIMD;
		try{
			Shared.SIMD=simd;
			final int[] values=new int[LENGTH];
			for(int i=0; i<values.length; i++){values[i]=i&1023;}
			final long warmup=repeat(values, 8192);
			check(warmup, 8192);
			final long start=System.nanoTime();
			final long checksum=repeat(values, loops);
			final long elapsed=System.nanoTime()-start;
			check(checksum, loops);
			System.out.println("INTEGER_SUM_BENCH simd="+simd+" length="+LENGTH+" loops="+loops+
					" elements="+((long)LENGTH*loops)+" wall_ns="+elapsed+
					" checksum="+checksum+" warmup_checksum="+warmup);
		}finally{
			Shared.SIMD=oldSimd;
		}
	}

	/** Sums after changing one element on each iteration.
	 * @param values Mutable input array
	 * @param loops Number of calls
	 * @return Accumulated results
	 */
	private static long repeat(final int[] values, final int loops){
		long checksum=0;
		for(int i=0; i<loops; i++){
			values[0]=i&1023;
			checksum+=Vector.sum(values);
		}
		return checksum;
	}

	/** Checks the arithmetic-series oracle without invoking the sum implementation.
	 * @param actual Observed accumulated sums
	 * @param loops Number of calls
	 */
	private static void check(final long actual, final int loops){
		final long base=(long)(LENGTH/1024)*1023*1024/2;
		final long remainder=loops%1024;
		final long expected=base*loops+(long)(loops/1024)*1023*1024/2+
				remainder*(remainder-1)/2;
		if(actual!=expected){
			throw new AssertionError("Retained sum checksum expected="+expected+" actual="+actual);
		}
	}

	/** Prevents construction of this standalone benchmark. */
	private VectorIntegerSumBenchmark(){}

	/** Array length, an exact multiple of the 1024-value input period. */
	private static final int LENGTH=65536;
}
