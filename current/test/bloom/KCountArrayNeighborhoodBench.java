package test.bloom;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.lang.reflect.Field;
import java.util.Locale;
import java.util.Random;

import bloom.KCountArray;
import bloom.KCountArray7MTA;
import shared.Shared;

/**
 * Fixed-work, single-thread timing of canonical left/right neighborhood sums.
 * Run once per fresh JVM with operations= and warmup= positive multiples of 512.
 * One operation is one readLeft or readRight call, not a pair or a cell lookup.
 * Each cycle visits 256 seeded keys, left then right for each key. Setup and
 * independent component-sum checks are outside both timed loops. Low counts
 * permit collisions without overflow; this is not an overflow correctness test.
 * The runner chooses a common workload and verifies at least ten main-thread
 * CPU seconds per measured sample; a short pilot is only calibration.
 *
 * @author Brian Bushnell, Sucrose
 */
public final class KCountArrayNeighborhoodBench{

	/** Parses explicit call counts, prepares one table and emits one measured sample. */
	public static void main(final String[] args) throws Exception{
		long operations=0, warmup=0;
		for(final String arg : args){
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected name=value: "+arg);}
			final String name=arg.substring(0, equals);
			final long count=Long.parseLong(arg.substring(equals+1));
			if(count<=0 || count%CYCLE_CALLS!=0){
				throw new IllegalArgumentException("Call count must be a positive multiple of "+CYCLE_CALLS+": "+arg);
			}
			if(name.equalsIgnoreCase("operations") && operations==0){operations=count;}
			else if(name.equalsIgnoreCase("warmup") && warmup==0){warmup=count;}
			else{throw new IllegalArgumentException("Unknown or duplicate argument: "+arg);}
		}
		if(operations==0 || warmup==0){
			throw new IllegalArgumentException("Required: operations=<calls> warmup=<calls>; multiples of 512.");
		}
		final long counter=readField(KCountArray7MTA.class, "counter", null);
		final long seedMask=readField(KCountArray7MTA.class, "SEEDMASK", null);
		if(counter!=0 || seedMask!=0 || KCountArray.verbose){
			throw new IllegalStateException("Requires fresh JVM, default mask seeds and verbose=false; counter="+
					counter+" seedMask="+seedMask+" verbose="+KCountArray.verbose);
		}
		final KCountArray7MTA table=new KCountArray7MTA(TARGET_CELLS, 32, 1, null, 0);
		final long partitions=readField(KCountArray.class, "numArrays", table);
		final long[] keys=populate(table);
		final long[] oracle=verify(table, keys, partitions);
		final long expectedWarmup=Math.multiplyExact(oracle[0], warmup/CYCLE_CALLS);
		final long expectedSample=Math.multiplyExact(oracle[0], operations/CYCLE_CALLS);
		System.out.println("CONFIG k="+K+" bits="+table.cellBits+" hashes=1 canonical=true queries="+QUERY_COUNT+
				" query_rng=java.util.Random query_seed="+QUERY_SEED+" mask_counter="+counter+" seed_mask="+seedMask+
				" target_cells="+TARGET_CELLS+" actual_cells="+table.cells+" partitions="+partitions+
				" cells_per_partition="+(table.cells/partitions)+" max_value="+table.maxValue+
				" setup_increment=1+(query_index&3) inserted_mass="+INSERTED_MASS+
				" worker_threads=1 shared_threads="+Shared.threads()+" available_processors="+
				Runtime.getRuntime().availableProcessors()+" assertions="+KCountArray7MTA.class.desiredAssertionStatus()+
				" operation_unit=neighbor_call cycle_calls="+CYCLE_CALLS+" cycle_checksum="+oracle[0]+
				" workload_fingerprint="+Long.toHexString(oracle[1])+" operations="+operations+" warmup="+warmup+
				" expected_sample_checksum="+expectedSample+" expected_warmup_checksum="+expectedWarmup);
		checkOutput();
		final ThreadMXBean bean=ManagementFactory.getThreadMXBean();
		if(!bean.isCurrentThreadCpuTimeSupported()){
			throw new IllegalStateException("Main-thread CPU time is required; wall time alone is insufficient.");
		}
		final boolean wasEnabled=bean.isThreadCpuTimeEnabled();
		try{
			if(!wasEnabled){bean.setThreadCpuTimeEnabled(true);}
			time(table, keys, warmup, expectedWarmup, "warmup", bean);
			time(table, keys, operations, expectedSample, "sample", bean);
		}finally{
			if(!wasEnabled){bean.setThreadCpuTimeEnabled(false);}
		}
	}

	/** Prevents instantiation of the standalone driver. */
	private KCountArrayNeighborhoodBench(){}

	/** Seeds all eight canonical neighbors per query with small public increments. */
	private static long[] populate(final KCountArray7MTA table){
		final long[] keys=new long[QUERY_COUNT];
		final Random random=new Random(QUERY_SEED);
		long inserted=0;
		for(int i=0; i<keys.length; i++){
			final long key=random.nextLong()&KEY_MASK;
			keys[i]=key;
			final int amount=1+(i&3);
			for(int base=0; base<4; base++){
				final long left=(key>>>2)|((long)base<<(2*(K-1)));
				final long right=((key<<2)&KEY_MASK)|base;
				table.increment(KCountArray.makeCanonical2(left, K), amount);
				table.increment(KCountArray.makeCanonical2(right, K), amount);
				inserted+=2*amount;
			}
		}
		assert(inserted==INSERTED_MASS) : "Low-count verification bound depends on total setup mass; observed="+inserted;
		return keys;
	}

	/** Records component counts and returns the cycle sum and a noncryptographic workload fingerprint. */
	private static long[] verify(final KCountArray7MTA table, final long[] keys, final long partitions){
		final int[] left=new int[4], right=new int[4];
		long cycle=0, fingerprint=mix(mix(QUERY_SEED, table.cells), partitions);
		for(int i=0; i<keys.length; i++){
			final long key=keys[i];
			table.readAllLeft(key, K, true, left);
			table.readAllRight(key, K, true, right);
			final long leftSum=sum(left), rightSum=sum(right);
			final int actualLeft=table.readLeft(key, K, true), actualRight=table.readRight(key, K, true);
			if(actualLeft!=leftSum || actualRight!=rightSum){
				throw new AssertionError("Component sums disagree at query "+i+": left="+actualLeft+
						" expected="+leftSum+" right="+actualRight+" expected="+rightSum);
			}
			cycle+=leftSum+rightSum;
			fingerprint=mix(fingerprint, key);
			for(int base=0; base<4; base++){
				fingerprint=mix(mix(fingerprint, left[base]), right[base]);
			}
			System.out.println("QUERY index="+i+" key="+key+" left="+left[0]+","+left[1]+","+left[2]+","+left[3]+
					" right="+right[0]+","+right[1]+","+right[2]+","+right[3]+
					" expected_left="+leftSum+" expected_right="+rightSum);
		}
		return new long[]{cycle, fingerprint};
	}

	/** Widens component addition and checks the known total-insertion bound before any timing. */
	private static long sum(final int[] counts){
		long result=0;
		for(final int count : counts){
			if(count<1 || count>INSERTED_MASS){
				throw new AssertionError("Every queried neighbor was seeded; count outside [1,"+INSERTED_MASS+"]: "+count);
			}
			result+=count;
		}
		return result;
	}

	/** Times one complete loop; setup, output and checksum validation remain outside clocks. */
	private static void time(final KCountArray7MTA table, final long[] keys, final long operations,
			final long expected, final String phase, final ThreadMXBean bean){
		final long cpuStart=bean.getCurrentThreadCpuTime();
		if(cpuStart<0){throw new IllegalStateException("Main-thread CPU time is unavailable.");}
		final long wallStart=System.nanoTime();
		final long checksum=run(table, keys, operations);
		final long wallNanos=System.nanoTime()-wallStart;
		final long cpuEnd=bean.getCurrentThreadCpuTime();
		final long cpuNanos=cpuEnd-cpuStart;
		if(cpuEnd<0 || cpuNanos<=0 || wallNanos<=0 || checksum!=expected){
			throw new AssertionError("Invalid "+phase+": checksum="+checksum+" expected="+expected+
					" wall_ns="+wallNanos+" thread_cpu_ns="+cpuNanos);
		}
		System.out.printf(Locale.ROOT,
				"RESULT phase=%s operations=%d checksum=%d wall_ns=%d thread_cpu_ns=%d wall_seconds=%.9f thread_cpu_seconds=%.9f%n",
				phase, operations, checksum, wallNanos, cpuNanos, wallNanos*1e-9, cpuNanos*1e-9);
		checkOutput();
	}

	/** Executes exactly the requested calls, reusing keys and retaining every returned sum. */
	private static long run(final KCountArray7MTA table, final long[] keys, final long operations){
		assert(keys.length==QUERY_COUNT) : "Masked cycle indexing requires exactly "+QUERY_COUNT+" keys; got "+keys.length;
		long checksum=0;
		for(long op=0; op<operations; op++){
			final long key=keys[(int)(op>>>1)&(QUERY_COUNT-1)];
			checksum+=((op&1)==0 ? table.readLeft(key, K, true) : table.readRight(key, K, true));
		}
		return checksum;
	}

	/** Reads seed/geometry state without changing the table's private fields. */
	private static long readField(final Class<?> owner, final String name, final Object instance) throws Exception{
		final Field field=owner.getDeclaredField(name);
		field.setAccessible(true);
		return ((Number)field.get(instance)).longValue();
	}

	/** Combines workload words with intentional modulo-2^64 wraparound; not a cryptographic hash. */
	private static long mix(final long hash, final long word){return (hash^word)*0x100000001b3L;}

	/** Makes a failed or incomplete stdout capture a visible failure. */
	private static void checkOutput(){
		if(System.out.checkError()){throw new IllegalStateException("Could not write complete benchmark output.");}
	}

	/** Fixed nucleotide length, query count and number of calls per full left/right cycle. */
	private static final int K=31, QUERY_COUNT=256, CYCLE_CALLS=2*QUERY_COUNT;
	/** Requested cell count; actual prime-partitioned geometry is reported at runtime. */
	private static final long TARGET_CELLS=1L<<20;
	/** Repeatable java.util.Random seed for the raw query keys. */
	private static final long QUERY_SEED=20261106;
	/** Mask retaining exactly 31 two-bit nucleotide symbols. */
	private static final long KEY_MASK=(1L<<(2*K))-1;
	/** Eight increments per query, with amounts 1,2,3,4 repeated 64 times. */
	private static final int INSERTED_MASS=8*(QUERY_COUNT/4)*10;
}
