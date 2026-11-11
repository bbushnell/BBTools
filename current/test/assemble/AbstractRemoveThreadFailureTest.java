package test.assemble;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import assemble.AbstractRemoveThread;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/**
 * Injected backing-shard failure contracts for AbstractRemoveThread.
 * Uses a tiny real table set and corrupts one shard to make the private
 * RemoveThread2 fail through the production worker path.
 * @author Brian Bushnell, Jean
 */
public final class AbstractRemoveThreadFailureTest{

	/** Checks that count-pruning survives and worker exceptions are not lost. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, AbstractRemoveThreadFailureTest.class, false);
		try{
			String input=null;
			Path root=null;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				if(pair[0].equalsIgnoreCase("in")){input=pair[1];}
				else if(pair[0].equalsIgnoreCase("outdir")){root=Paths.get(pair[1]).toAbsolutePath();}
				else{throw new IllegalArgumentException("Unknown test option: "+arg);}
			}
			if(input==null || input.length()==0 || root==null){
				throw new IllegalArgumentException("Specify in=existing.fa and new outdir=DIR.");
			}
			Files.createDirectory(root);

			int failed=0;
			final StringBuilder checks=new StringBuilder("case\tpassed\tdetail\n");
			failed+=control(checks, input);
			failed+=workerException(checks, input);

			Files.write(root.resolve("checks.tsv"), checks.toString().getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
			System.out.println("ABSTRACT_REMOVE_THREAD_FAILURE_TEST checks=2 failed="+failed);
			if(failed>0){throw new AssertionError("AbstractRemoveThread worker-failure contract failures; inspect "+root);}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** One ordinary shard should be pruned to completion. */
	private static int control(final StringBuilder checks, final String input) throws Exception{
		final KmerTableSetU tables=tables(input);
		final Kmer kmer=insertOneKmer(tables, 0);
		final int before=tables.getCount(kmer);
		final long removed=AbstractRemoveThread.process(1, 2, Integer.MAX_VALUE, tables, false);
		final int after=tables.getCount(kmer);
		final int liveThreads=liveRemoveThreads();
		return check(checks, "control", before==1 && removed==1 && after<1 && liveThreads==0,
				"removed="+removed+", before="+before+", after="+after+", liveThreads="+liveThreads);
	}

	/**
	 * A corrupt shard should make process() fail, but only after all workers join
	 * and the independent live shards have finished pruning.
	 */
	private static int workerException(final StringBuilder checks, final String input) throws Exception{
		final KmerTableSetU tables=tables(input);
		final Kmer first=insertOneKmer(tables, 0);
		final Kmer second=insertOneKmer(tables, 2);
		final int before=countPositive(tables, first, second);
		nullShard(tables, 1);

		try{
			final long removed=AbstractRemoveThread.process(3, 2, Integer.MAX_VALUE, tables, false);
			final int after=countPositive(tables, first, second);
			final int liveThreads=liveRemoveThreads();
			return check(checks, "worker-null-shard", false,
					"returned removed="+removed+", before="+before+", after="+after+", liveThreads="+liveThreads);
		}catch(Throwable t){
			final Throwable cause=t.getCause()==null ? t : t.getCause();
			if(!(cause instanceof NullPointerException) || !isWorkerNullShardFailure(cause)){
				throw new AssertionError("Unexpected exception before expected worker failure", t);
			}
			final int after=countPositive(tables, first, second);
			final int liveThreads=liveRemoveThreads();
			return check(checks, "worker-null-shard", before==2 && after==0 && liveThreads==0,
					"threw="+cause.getClass().getSimpleName()+", before="+before+", after="+after+", liveThreads="+liveThreads);
		}
	}

	/** Creates a small explicit unlimited table set with exactly three shards. */
	private static KmerTableSetU tables(final String input){
		final String[] options={"k=63", "packed=t", "rcomp=f", "ways=3", "initialsize=16", "prealloc=f",
				"prefilter=f", "minprob=0", "minprobmain=f", "interleaved=f", "t=1", "in="+input};
		return new AllocatedKmerTableSetU(options);
	}

	/** Inserts one k-mer into the requested way. */
	private static Kmer insertOneKmer(final KmerTableSetU tables, final int way){
		for(int seed=0; seed<10000; seed++){
			final Kmer kmer=kmerForSeed(seed);
			if(tables.kmerToWay(kmer)==way){
				if(tables.getTable(kmer).incrementAndReturnNumCreated(kmer)!=1){
					throw new AssertionError("Expected a fresh k-mer for seed "+seed);
				}
				return kmer;
			}
		}
		throw new AssertionError("Could not find a k-mer for way "+way);
	}

	/** Builds a deterministic 63-base key. */
	private static Kmer kmerForSeed(final int seed){
		final Kmer kmer=new Kmer(63, false);
		long x=0x9E3779B97F4A7C15L ^ seed;
		for(int i=0; i<63; i++){
			x^=(x<<13);
			x^=(x>>>7);
			x^=(x<<17);
			kmer.addRight(BASES[(int)(x&3)]);
		}
		return kmer;
	}

	/** Counts how many of the tracked k-mers are still present. */
	private static int countPositive(final KmerTableSetU tables, final Kmer first, final Kmer second){
		int positive=0;
		if(tables.getCount(first)>0){positive++;}
		if(tables.getCount(second)>0){positive++;}
		return positive;
	}

	/** Replaces one backing shard with null so the private worker throws there. */
	private static void nullShard(final KmerTableSetU tables, final int way) throws Exception{
		final Field field=KmerTableSetU.class.getDeclaredField("tables");
		field.setAccessible(true);
		final Object[] shards=(Object[])field.get(tables);
		shards[way]=null;
	}

	/** Returns true only for the expected private unlimited-table removal worker. */
	private static boolean isWorkerNullShardFailure(final Throwable cause){
		for(StackTraceElement element : cause.getStackTrace()){
			if("assemble.AbstractRemoveThread$RemoveThread2".equals(element.getClassName()) &&
					"processNextTable".equals(element.getMethodName())){
				return true;
			}
		}
		return false;
	}

	/** Counts removal workers that survived process() in this fresh JVM. */
	private static int liveRemoveThreads(){
		int live=0;
		for(Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread instanceof AbstractRemoveThread && thread.isAlive()){live++;}
		}
		return live;
	}

	/** Exposes the protected allocation step without loading reads. */
	private static final class AllocatedKmerTableSetU extends KmerTableSetU{
		AllocatedKmerTableSetU(final String[] args){
			super(args, 0);
			allocateTables();
		}
	}

	/** Keeps every scientific expectation active without depending on JVM assertions. */
	private static int check(final StringBuilder checks, final String name,
			final boolean passed, final String detail){
		checks.append(name).append('\t').append(passed ? 1 : 0).append('\t').append(detail).append('\n');
		return passed ? 0 : 1;
	}

	private static final char[] BASES={'A', 'C', 'G', 'T'};

}
