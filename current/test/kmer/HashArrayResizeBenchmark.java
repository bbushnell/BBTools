package test.kmer;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.Locale;

import kmer.HashArray2D;

/**
 * Matched single-thread timing of multi-value table construction and resizing.
 * Consecutive keys avoid the old resize-to-victim failure on both revisions.
 * rounds=0 calibrates at least 15 seconds of accumulated current-thread CPU time
 * during construction/insertion; replay the printed rounds on both revisions.
 * The seconds column reports kernel wall time, excluding input-array construction
 * and the full oracle. CPU accounting must be supported and enabled in this JVM.
 * The runner must also retain full process wall, user CPU and system CPU with
 * /usr/bin/time, including those costs and JVM startup.
 * @author Sucrose
 */
public final class HashArrayResizeBenchmark {

	/** Calibrates or replays fixed work, verifying every table before discarding it. */
	public static void main(String[] args){
		int requested=0;
		for(String arg : args){
			if(!arg.startsWith("rounds=")){throw new IllegalArgumentException("Expected rounds=: "+arg);}
			requested=Integer.parseInt(arg.substring(7));
		}
		if(requested<0){throw new IllegalArgumentException("Round count must be nonnegative: "+requested);}
		final ThreadMXBean cpu=cpuClock();
		final int[][] inputs=makeInputs();
		long kernelNanos=0, kernelCpuNanos=0, checksum=0;
		int rounds=0;
		do{
			final long cpuStart=cpu.getCurrentThreadCpuTime();
			final long start=System.nanoTime();
			final HashArray2D table=new HashArray2D(SCHEDULE, Long.MAX_VALUE);
			for(int key=0; key<KEYS; key++){table.set(key, inputs[key], VALUES);}
			kernelNanos+=System.nanoTime()-start;
			final long cpuEnd=cpu.getCurrentThreadCpuTime();
			if(cpuStart<0 || cpuEnd<cpuStart){
				throw new IllegalStateException("Current-thread CPU accounting became unavailable or decreased: "+cpuStart+", "+cpuEnd);
			}
			kernelCpuNanos+=cpuEnd-cpuStart;
			checksum+=verify(table);
			rounds++;
		}while(requested==0 ? kernelCpuNanos<TARGET_CPU_NANOS : rounds<requested);
		blackhole=checksum;
		System.err.printf(Locale.ROOT, "kernelCpuSeconds=%.6f%n", kernelCpuNanos*1e-9);
		System.out.printf(Locale.ROOT,
			"mode\trounds\tseconds\tchecksum%nresize\t%d\t%.6f\t%d%n",
			rounds, kernelNanos*1e-9, blackhole);
	}

	/** Requires a usable Java8 CPU clock; wall time cannot substitute for calibration work. */
	private static ThreadMXBean cpuClock(){
		final ThreadMXBean cpu=ManagementFactory.getThreadMXBean();
		if(!cpu.isCurrentThreadCpuTimeSupported()){
			throw new IllegalStateException("Resize calibration requires current-thread CPU-time support");
		}
		if(!cpu.isThreadCpuTimeEnabled()){
			try{cpu.setThreadCpuTimeEnabled(true);}
			catch(UnsupportedOperationException | SecurityException e){
				throw new IllegalStateException("Could not enable current-thread CPU timing for resize calibration", e);
			}
		}
		if(!cpu.isThreadCpuTimeEnabled()){
			throw new IllegalStateException("Current-thread CPU timing remains disabled for resize calibration");
		}
		return cpu;
	}

	/** Builds distinct eight-value lists once; tables borrow them without changing their contents. */
	private static int[][] makeInputs(){
		final int[][] inputs=new int[KEYS][];
		for(int key=0; key<KEYS; key++){
			final int[] values=new int[VALUES+1];
			for(int j=0; j<VALUES; j++){values[j]=key*VALUES+j+1;}
			values[VALUES]=-1;
			inputs[key]=values;
		}
		return inputs;
	}

	/**
	 * Checks all keys, all literal values, terminators, final capacity and both
	 * storage sizes. Expected values are computed from key/index, not input arrays.
	 * Checks remain enabled under both -ea and -da.
	 */
	private static long verify(final HashArray2D table){
		if(table.size()!=KEYS || table.victims().size()!=0){
			throw new AssertionError("Expected 32768 primary keys and no victims: "+table.size()+", "+table.victims().size());
		}
		if(table.arrayLength()!=67469){
			throw new AssertionError("Expected final modulus 67409 plus 60 probe cells: "+table.arrayLength());
		}
		final int[] singleton=new int[1];
		long sum=0;
		for(int key=0; key<KEYS; key++){
			final int[] actual=table.getValues(key, singleton);
			if(actual==null || actual.length!=VALUES+1){throw new AssertionError("Wrong value-list length at key "+key);}
			for(int j=0; j<VALUES; j++){
				final int expected=key*VALUES+j+1;
				if(actual[j]!=expected){
					throw new AssertionError("Changed value at key "+key+", index "+j+": "+actual[j]+" versus "+expected);
				}
				sum+=actual[j];
			}
			if(actual[VALUES]!=-1){throw new AssertionError("Lost terminator at key "+key);}
		}
		if(table.contains(KEYS)){throw new AssertionError("Unexpected key beyond the inserted range");}
		return sum;
	}

	/** Repeated tables use the identical schedule on both source revisions. */
	private static final int[] SCHEDULE={127, 257, 521, 1049, 2099, 4201, 8419, 16843, 33703, 67409};
	private static final int KEYS=32768;
	private static final int VALUES=8;
	private static final long TARGET_CPU_NANOS=15000000000L;
	/** Retains consumed table contents beyond the measured and verified work. */
	private static volatile long blackhole;
}
