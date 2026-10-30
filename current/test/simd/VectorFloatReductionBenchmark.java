package test.simd;

import shared.Shared;
import simd.Vector;

/** Measures repeated floating reductions with retained, exact per-call results.
 * Mutates preallocated integer-valued float arrays before each call so sums and
 * maxima both vary while every expected result remains exactly representable.
 * @author Sucrose */
public final class VectorFloatReductionBenchmark{

	/** Runs warmup and a fixed amount of timed single-threaded work.
	 * @param args op=sum|max, loops=positive integer, length=positive integer, and simd=t|f
	 */
	public static void main(final String[] args){
		String op=OP_SUM;
		int loops=65536, length=65536;
		boolean simd=false;
		for(final String arg : args){
			if(arg.startsWith("op=")){op=arg.substring(3);}
			else if(arg.startsWith("loops=")){loops=Integer.parseInt(arg.substring(6));}
			else if(arg.startsWith("length=")){length=Integer.parseInt(arg.substring(7));}
			else if(arg.equals("simd=t")){simd=true;}
			else if(arg.equals("simd=f")){simd=false;}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		if(!validOp(op)){throw new IllegalArgumentException("Unknown op: "+op);}
		if(loops<1){throw new IllegalArgumentException("loops must be positive: "+loops);}
		if(length<1){throw new IllegalArgumentException("length must be positive: "+length);}
		if(simd && !Vector.simd256){
			throw new IllegalStateException("SIMD benchmark requires 256-bit Vector support");
		}
		final Workload workload=OP_SUM.equals(op) ? sumWorkload(length) : maxWorkload(length);
		final boolean oldSimd=Shared.SIMD;
		try{
			Shared.SIMD=simd;
			final long warmup=repeat(op, workload, WARMUP_LOOPS);
			final long start=System.nanoTime();
			final long checksum=repeat(op, workload, loops);
			final long elapsed=System.nanoTime()-start;
			System.out.println("FLOAT_REDUCTION_BENCH op="+op+" simd="+simd+
					" length="+length+" loops="+loops+" buffers="+BANKS+
					" calls="+loops+" elements="+((long)length*loops)+
					" wall_ns="+elapsed+" checksum="+checksum+
					" warmup_checksum="+warmup);
		}finally{
			Shared.SIMD=oldSimd;
		}
	}

	/** Repeats the selected reduction over a mutable bank.
	 * @param op Operation name
	 * @param workload Mutable benchmark state
	 * @param loops Number of calls
	 * @return Retained checksum of exact reduction results and touched values
	 */
	private static long repeat(final String op, final Workload workload, final int loops){
		long checksum=0;
		for(int i=0; i<loops; i++){
			final int bank=i&(BANKS-1);
			if(OP_SUM.equals(op)){
				final int changed=mutateSum(workload, bank);
				final double actual=Vector.sum(workload.arrays[bank]);
				if(Double.doubleToLongBits(actual)!=Double.doubleToLongBits(workload.sums[bank])){
					throw new AssertionError("sum bank="+bank+" expected="+
							workload.sums[bank]+" actual="+actual);
				}
				checksum=mix(mix(checksum, Double.doubleToLongBits(actual)), changed);
			}else{
				final int expected=mutateMax(workload, bank);
				final float actual=Vector.max(workload.arrays[bank]);
				if(Float.floatToIntBits(actual)!=Float.floatToIntBits((float)expected)){
					throw new AssertionError("max bank="+bank+" expected="+
							expected+" actual="+actual);
				}
				checksum=mix(mix(checksum, Float.floatToIntBits(actual)), expected);
			}
		}
		return checksum;
	}

	/** Mutates one float in a summation bank and updates its exact sum.
	 * @param workload Mutable benchmark state
	 * @param bank Bank index
	 * @return New value written at the cursor
	 */
	private static int mutateSum(final Workload workload, final int bank){
		final float[] values=workload.arrays[bank];
		final int pos=workload.cursors[bank];
		int nextPos=pos+1;
		if(nextPos>=values.length){nextPos=0;}
		workload.cursors[bank]=nextPos;
		final int old=(int)values[pos];
		final int value=(old>=MAX_SUM_VALUE ? MIN_SUM_VALUE : old+1);
		values[pos]=value;
		workload.sums[bank]+=value-old;
		return value;
	}

	/** Moves one positive maximum to the next cursor position.
	 * @param workload Mutable benchmark state
	 * @param bank Bank index
	 * @return New exact maximum
	 */
	private static int mutateMax(final Workload workload, final int bank){
		final float[] values=workload.arrays[bank];
		final int pos=workload.cursors[bank];
		values[pos]=MAX_BACKGROUND;
		int nextPos=pos+1;
		if(nextPos>=values.length){nextPos=0;}
		workload.cursors[bank]=nextPos;
		final int value=(workload.maxima[bank]>=MAX_SUM_VALUE ? 1 : workload.maxima[bank]+1);
		workload.maxima[bank]=value;
		values[nextPos]=value;
		return value;
	}

	/** Builds the exact mutable summation workload.
	 * @param length Array length
	 * @return Mutable benchmark state
	 */
	private static Workload sumWorkload(final int length){
		final Workload workload=new Workload(length);
		for(int bank=0; bank<BANKS; bank++){
			final float[] values=workload.arrays[bank];
			double sum=0;
			for(int i=0; i<values.length; i++){
				final int value=((i+bank*7)&31)-16;
				values[i]=value;
				sum+=value;
			}
			workload.sums[bank]=sum;
			workload.cursors[bank]=bank%length;
		}
		return workload;
	}

	/** Builds the mutable maximum workload with one positive maximum per bank.
	 * @param length Array length
	 * @return Mutable benchmark state
	 */
	private static Workload maxWorkload(final int length){
		final Workload workload=new Workload(length);
		for(int bank=0; bank<BANKS; bank++){
			final float[] values=workload.arrays[bank];
			for(int i=0; i<values.length; i++){values[i]=MAX_BACKGROUND;}
			final int pos=bank%length;
			final int value=1+(bank&31);
			values[pos]=value;
			workload.cursors[bank]=pos;
			workload.maxima[bank]=value;
		}
		return workload;
	}

	/** Checks whether an operation name is supported.
	 * @param op Operation name
	 * @return True if supported
	 */
	private static boolean validOp(final String op){
		return OP_SUM.equals(op) || OP_MAX.equals(op);
	}

	/** Mixes a retained long value into the checksum.
	 * @param checksum Previous checksum
	 * @param value Retained value
	 * @return Updated checksum
	 */
	private static long mix(long checksum, final long value){
		checksum^=value+0x9E3779B97F4A7C15L;
		checksum*=0xBF58476D1CE4E5B9L;
		return Long.rotateLeft(checksum, 11);
	}

	/** Prevents construction of this standalone benchmark. */
	private VectorFloatReductionBenchmark(){}

	/** Mutable arrays and expected reductions. */
	private static final class Workload{
		/** Preallocates the input bank and per-bank cursors.
		 * @param length Array length
		 */
		private Workload(final int length){
			arrays=new float[BANKS][length];
			cursors=new int[BANKS];
			sums=new double[BANKS];
			maxima=new int[BANKS];
		}

		/** Mutable reduction inputs. */
		final float[][] arrays;
		/** Per-bank mutation cursors, reset before they reach length. */
		final int[] cursors;
		/** Exact expected sums. */
		final double[] sums;
		/** Exact expected maxima. */
		final int[] maxima;
	}

	/** Number of mutable input buffers. */
	private static final int BANKS=128;
	/** Warmup loops run before timed work. */
	private static final int WARMUP_LOOPS=8192;
	/** Minimum summation input. */
	private static final int MIN_SUM_VALUE=-127;
	/** Maximum ordinary input and varying maximum. */
	private static final int MAX_SUM_VALUE=127;
	/** Background below every generated maximum. */
	private static final int MAX_BACKGROUND=-1000;
	/** Operation selecting Vector.sum(float[]). */
	private static final String OP_SUM="sum";
	/** Operation selecting Vector.max(float[]). */
	private static final String OP_MAX="max";
}
