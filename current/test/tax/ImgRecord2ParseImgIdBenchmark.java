package test.tax;

import parse.Parse;
import tax.ImgRecord2;

/** Measures valid IMG header ID parsing over a deterministic prebuilt corpus.
 * Corpus construction and known-ID result checks precede warmup and measurement.
 * @author Sucrose */
public final class ImgRecord2ParseImgIdBenchmark{

	/** Measures repeated parsing of a small varied valid-IMG header corpus.
	 * @param args Optional iterations=N (default 67,108,864) and doAssertions=t/f (default t) */
	public static void main(String[] args){
		int iterations=64*1024*1024;
		boolean doAssertions=true;
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Expected key=value: "+arg);}
			if(split[0].equalsIgnoreCase("iterations")){iterations=Parse.parseIntKMG(split[1]);}
			else if(split[0].equalsIgnoreCase("doAssertions")){
				doAssertions=Parse.parseBoolean(split[1]);
			}else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		if(iterations<1){throw new IllegalArgumentException("iterations must be positive");}

		final String[] headers=new String[HEADERS];
		final long[] ids=new long[HEADERS];
		for(int i=0; i<headers.length; i++){
			ids[i]=expectedId(i);
			headers[i]=makeHeader(i, ids[i]);
			final long actual=ImgRecord2.parseImgId(headers[i], doAssertions);
			if(actual!=ids[i]){
				throw new AssertionError("IMG header oracle mismatch: "+headers[i]+
						" -> "+actual+", expected "+ids[i]);
			}
		}
		final long expected=expected(ids, iterations);
		if(run(headers, WARMUP, doAssertions)!=(expected(ids, WARMUP))){
			throw new AssertionError("Warmup checksum mismatch");
		}
		final long start=System.nanoTime();
		final long sum=run(headers, iterations, doAssertions);
		final long elapsed=System.nanoTime()-start;
		if(sum!=expected){throw new AssertionError("Timed checksum mismatch");}
		System.out.println("checksum="+sum+" expected="+expected+
				" iterations="+iterations+" headers="+headers.length+
				" doAssertions="+doAssertions+" compute_ns="+elapsed);
	}

	/** Prevents construction of this standalone benchmark. */
	private ImgRecord2ParseImgIdBenchmark(){}

	/** Creates one valid IMG header with varied delimiters and token positions.
	 * @param index Header slot
	 * @param id Expected IMG ID
	 * @return Header containing exactly one valid IMG token */
	private static String makeHeader(final int index, final long id){
		switch(index&3){
			case 0: return "read|img|"+id+"|tid|7";
			case 1: return "read~img~"+id+"~tid~7";
			case 2: return "read|prefix|img|"+id+" description";
			default: return "read~prefix~img~0000"+id+" description";
		}
	}

	/** Returns mostly realistic 10-12 digit IMG IDs, with small controls retained.
	 * @param index Header slot
	 * @return Expected IMG ID */
	private static long expectedId(final int index){
		if(index==0){return 1;}
		if(index==1){return 123456789;}
		return 1000000000L+((long)index*98245165307L)%900000000000L;
	}

	/** Parses prebuilt valid IMG headers.
	 * @param headers Power-of-two header corpus
	 * @param iterations Total parser calls
	 * @param doAssertions Parser assertion mode to benchmark
	 * @return Sum of all parsed IMG IDs with Java long wraparound */
	private static long run(final String[] headers, final int iterations,
			final boolean doAssertions){
		assert(headers.length>0 && (headers.length&(headers.length-1))==0) :
				"Benchmark indexing requires a nonempty power-of-two header array.";
		long sum=0;
		final int mask=headers.length-1;
		for(int i=0; i<iterations; i++){
			sum+=ImgRecord2.parseImgId(headers[i&mask], doAssertions);
		}
		return sum;
	}

	/** Computes the checksum expected from {@link #run(String[], int, boolean)}.
	 * @param ids Expected IMG IDs in corpus order
	 * @param iterations Total parser calls
	 * @return Expected deterministic checksum with Java long wraparound */
	private static long expected(final long[] ids, final int iterations){
		long cycle=0, tail=0;
		for(long id : ids){cycle+=id;}
		for(int i=0; i<(iterations&(ids.length-1)); i++){tail+=ids[i];}
		return cycle*(iterations/ids.length)+tail;
	}

	/** Number of valid headers in the prebuilt corpus. */
	private static final int HEADERS=16;

	/** Warmup parser calls run before timing starts. */
	private static final int WARMUP=1<<20;
}
