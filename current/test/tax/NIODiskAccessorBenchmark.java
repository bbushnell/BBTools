package test.tax;

import java.nio.file.Files;
import java.nio.file.Path;

import parse.Parse;
import tax.DiskAccessionTable;
import tax.NIODiskAccessor;

/** Repeated NIO disk-accession lookup workload with a fixed mixed query batch.
 * Table construction and expected-result construction are outside timing.
 * @author Jean */
public final class NIODiskAccessorBenchmark{

	/** Measures a tiny positioned-read table with one or four reader threads.
	 * @param args threads=1|4, batch=200, rounds=100000, warmup=auto
	 * @throws Exception If table setup, accessor construction or cleanup fails */
	public static void main(String[] args) throws Exception{
		int threads=1, batch=200, rounds=100000, warmup=-1;
		for(String arg : args){
			String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Expected key=value: "+arg);}
			String key=split[0].toLowerCase(), value=split[1];
			if(key.equals("threads") || key.equals("t")){threads=Parse.parseIntKMG(value);}
			else if(key.equals("batch")){batch=Parse.parseIntKMG(value);}
			else if(key.equals("rounds")){rounds=Parse.parseIntKMG(value);}
			else if(key.equals("warmup")){
				warmup=value.equalsIgnoreCase("auto") ? -1 : Parse.parseIntKMG(value);
			}else{throw new IllegalArgumentException("Unknown parameter "+arg);}
		}
		if(threads!=1 && threads!=4){throw new IllegalArgumentException("threads must be 1 or 4: "+threads);}
		if(batch<1 || rounds<1){throw new IllegalArgumentException("batch and rounds must be positive");}
		if(warmup<0){warmup=warmupRounds(batch);}
		if(warmup<1){throw new IllegalArgumentException("warmup must be positive");}
		if(threads==4 && batch<MIN_FOUR_THREAD_BATCH){
			throw new IllegalArgumentException("batch must be at least "+
					MIN_FOUR_THREAD_BATCH+" for four NIO workers");
		}

		final Path root=Files.createTempDirectory("nio-disk-benchmark-");
		final Path table=root.resolve("accession.bin");
		try{
			buildFourKeyTable(table);
			final long[] keys=keys(batch);
			final int[] expected=expected(keys);
			final long expectedSum=sum(expected);
			NIODiskAccessor accessor=null;
			try{
				accessor=new NIODiskAccessor(table.toString(), threads);
				check(accessor.getBatch(keys), expected);
				run(accessor, keys, warmup, expectedSum);
				final long start=System.nanoTime();
				final long checksum=run(accessor, keys, rounds, expectedSum);
				final long elapsed=System.nanoTime()-start;
				System.out.println("checksum="+checksum+" expected="+
						(expectedSum*rounds)+" batch="+batch+" rounds="+rounds+
						" threads="+threads+" warmup="+warmup+" compute_ns="+elapsed);
			}finally{
				if(accessor!=null){accessor.close();}
			}
		}finally{
			Files.deleteIfExists(table);
			Files.deleteIfExists(root);
		}
	}

	/** Prevents construction of this standalone benchmark. */
	private NIODiskAccessorBenchmark(){}

	/** Builds the four-key disk table used by every benchmark run.
	 * @param path Output table path
	 * @throws Exception If table construction fails */
	private static void buildFourKeyTable(final Path path) throws Exception{
		try(DiskAccessionTable.Builder builder=DiskAccessionTable.build(path.toString(), 4)){
			builder.add(1, 11);
			builder.add(2, 22);
			builder.add(3, 33);
			builder.add(4, 44);
		}
	}

	/** Creates a fixed present/missing query batch.
	 * @param batch Number of queries per batch
	 * @return Reused query keys */
	private static long[] keys(final int batch){
		final long[] keys=new long[batch];
		for(int i=0; i<keys.length; i++){keys[i]=PATTERN[i%PATTERN.length];}
		return keys;
	}

	/** Creates exact expected taxIDs for a query batch.
	 * @param keys Query keys
	 * @return Expected taxID result for every query */
	private static int[] expected(final long[] keys){
		final int[] expected=new int[keys.length];
		for(int i=0; i<keys.length; i++){
			final long key=keys[i];
			expected[i]=(int)(key>=1 && key<=4 ? key*11 : -1);
		}
		return expected;
	}

	/** Runs repeated batch lookups and validates each batch checksum.
	 * @param accessor Accessor under measurement
	 * @param keys Reused query batch
	 * @param rounds Number of repeated batches
	 * @param expectedSum Expected checksum per batch
	 * @return Sum over all result values
	 */
	private static long run(final NIODiskAccessor accessor, final long[] keys,
			final int rounds, final long expectedSum){
		long sum=0;
		for(int round=0; round<rounds; round++){
			final int[] result=accessor.getBatch(keys);
			final long batchSum=sum(result);
			if(batchSum!=expectedSum){
				throw new AssertionError("Batch checksum mismatch: "+batchSum+
						" != "+expectedSum);
			}
			sum+=batchSum;
		}
		return sum;
	}

	/** Sums a batch of lookup results.
	 * @param result Lookup results
	 * @return Sum of all result values */
	private static long sum(final int[] result){
		long sum=0;
		for(int value : result){sum+=value;}
		return sum;
	}

	/** Checks the pre-timing exact-result control.
	 * @param result Observed lookup results
	 * @param expected Expected lookup results */
	private static void check(final int[] result, final int[] expected){
		if(result.length!=expected.length){throw new AssertionError("Result length mismatch");}
		for(int i=0; i<result.length; i++){
			if(result[i]!=expected[i]){
				throw new AssertionError("Result mismatch at "+i+": "+
						result[i]+" != "+expected[i]);
			}
		}
	}

	/** Selects enough untimed batches to compile the hot lookup path.
	 * @param batch Number of queries per batch
	 * @return Untimed warmup round count */
	private static int warmupRounds(final int batch){
		return Math.max(MIN_WARMUP_ROUNDS,
				(MIN_WARMUP_LOOKUPS+batch-1)/batch);
	}

	/** Fixed mix of present and absent query keys. */
	private static final long[] PATTERN={1, 2, 99, 3, 100, 4, 101, 102};
	/** Smallest batch that current NIODiskAccessor settings split over four workers. */
	private static final int MIN_FOUR_THREAD_BATCH=200;
	/** Minimum untimed batches before measurement. */
	private static final int MIN_WARMUP_ROUNDS=32;
	/** Minimum untimed key lookups before measurement. */
	private static final int MIN_WARMUP_LOOKUPS=50000;
}
