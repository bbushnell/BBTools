package test.jgi;

import jgi.GradeMergedReads;

/** Direct microbenchmark for {@link GradeMergedReads#parseInsert(String)}. */
public final class GradeMergedReadsInsertBenchmark{

	public static void main(String[] args){
		final long iterations=parseIterations(args);
		final String[] headers=new String[HEADER_COUNT];
		final int[] expected=new int[HEADER_COUNT];
		fillHeaders(headers, expected);
		final long expectedPerIteration=sum(expected);
		final long expectedWarmup=expectedChecksum(expectedPerIteration, WARMUP_ITERATIONS);
		final long warmupChecksum=run(headers, WARMUP_ITERATIONS);
		if(warmupChecksum!=expectedWarmup){
			throw new AssertionError("Warmup checksum mismatch: expected "+expectedWarmup+", observed "+warmupChecksum);
		}
		final long expectedTimed=expectedChecksum(expectedPerIteration, iterations);

		final long start=System.nanoTime();
		final long checksum=run(headers, iterations);
		final long elapsed=System.nanoTime()-start;
		if(checksum!=expectedTimed){
			throw new AssertionError("Timed checksum mismatch: expected "+expectedTimed+", observed "+checksum);
		}
		System.out.println("GradeMergedReadsInsertBenchmark headers="+HEADER_COUNT+
				" iterations="+iterations+
				" calls="+expectedCalls(iterations)+
				" checksum="+checksum+
				" elapsedNanos="+elapsed);
	}

	private static long parseIterations(final String[] args){
		long iterations=DEFAULT_ITERATIONS;
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			if(arg.startsWith("iterations=")){
				iterations=Long.parseLong(arg.substring("iterations=".length()));
			}else{
				throw new IllegalArgumentException("Unknown argument: "+arg);
			}
		}
		if(iterations<1){
			throw new IllegalArgumentException("iterations must be positive: "+iterations);
		}
		return iterations;
	}

	private static void fillHeaders(final String[] headers, final int[] expected){
		for(int i=0; i<headers.length; i++){
			final int insert=100+(int)((i*37L+13)%900L);
			expected[i]=insert;
			headers[i]="f_"+(i&7)+"_c_"+(1000+i)+"_s_"+(i&1)+"_p_"+(100000+i*31)+
					"_i_"+insert+"_r_"+(insert+500);
		}
	}

	private static long sum(final int[] expected){
		long sum=0;
		for(int i=0; i<expected.length; i++){
			sum+=expected[i];
		}
		return sum;
	}

	private static long run(final String[] headers, final long iterations){
		long sum=0;
		for(long i=0; i<iterations; i++){
			for(int j=0; j<headers.length; j++){
				sum+=GradeMergedReads.parseInsert(headers[j]);
			}
		}
		return sum;
	}

	private static long expectedChecksum(final long expectedPerIteration, final long iterations){
		if(iterations>Long.MAX_VALUE/expectedPerIteration){
			throw new IllegalArgumentException("iterations would overflow checksum: "+iterations);
		}
		return expectedPerIteration*iterations;
	}

	private static long expectedCalls(final long iterations){
		// 1024 headers keeps this multiplication exact for any accepted checksum.
		return iterations*HEADER_COUNT;
	}

	private static final int HEADER_COUNT=1024;
	private static final long WARMUP_ITERATIONS=1024;
	private static final long DEFAULT_ITERATIONS=1024;
}
