package test.assemble;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import assemble.Tadpole;
import assemble.Tadpole1;
import assemble.Tadpole2;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import ukmer.Kmer;

/** Finite insert-distance contracts on two disjoint five-edge paths.
 * Loads real count tables, verifies each word and all four successors directly,
 * then compares walkers with coordinate-derived distances. Does not change them.
 * @author Brian Bushnell, Shinobu
 */
public final class TadpoleInsertDistanceTest{

	/** Runs one representation in a fresh directory; failed contracts remain in checks.tsv. */
	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, TadpoleInsertDistanceTest.class, false);
		try{
			int k=31;
			long iterations=0, warmup=1000000;
			Path root=null;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String value=pair.length>1 ? pair[1] : null;
				if(pair[0].equalsIgnoreCase("k")){k=Integer.parseInt(value);}
				else if(pair[0].equalsIgnoreCase("iterations")){iterations=Long.parseLong(value);}
				else if(pair[0].equalsIgnoreCase("warmup")){warmup=Long.parseLong(value);}
				else if(pair[0].equalsIgnoreCase("outdir")){root=Paths.get(value).toAbsolutePath();}
				else{throw new IllegalArgumentException("Unknown test option: "+arg);}
			}
			if(root==null || (k!=31 && k!=63)){throw new IllegalArgumentException("Specify new outdir=DIR and k=31 or63.");}
			if(iterations<0 || iterations%10!=0 || warmup<0 || warmup%10!=0){
				throw new IllegalArgumentException("Iteration counts must be nonnegative multiples of ten.");
			}
			Files.createDirectory(root);
			final String[] sequences=k==31 ? SHORT_PATHS : LONG_PATHS;
			final byte[][] paths={sequences[0].getBytes(StandardCharsets.US_ASCII), sequences[1].getBytes(StandardCharsets.US_ASCII)};
			final StringBuilder fasta=new StringBuilder();
			for(int path=0; path<paths.length; path++){
				if(paths[path].length!=k+5){throw new IllegalStateException("Fixture must have exactly five edges.");}
				for(int copy=0; copy<4; copy++){fasta.append('>').append(path).append('_').append(copy).append('\n').append(sequences[path]).append('\n');}
			}
			final Path input=root.resolve("paths.fa");
			Files.write(input, fasta.toString().getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
			final String[] options={"in="+input, "k="+k, "packed=t", "t=1", "buildthreads=1",
					"ways=1", "initialsize=128", "prealloc=f", "prefilter=0", "minprob=0", "minprobmain=f",
					"mcs=1", "mce=1", "shave=f", "rinse=f", "interleaved=f", "pigz=f", "unpigz=f"};
			Files.write(root.resolve("table.args"), java.util.Arrays.asList(options), StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
			final Tadpole assembler=Tadpole.makeTadpole(options, true);
			try{
				assembler.loadKmers(new Timer());
				for(byte[] path : paths){verifyPath(assembler, path, k);}
				final StringBuilder checks=new StringBuilder("case\texpected\tactual\tpassed\n");
				int failed=0;
				failed+=check(checks, "reachable-with-tail", assembler, word(paths[0], 0, k), word(paths[0], 1, k), 100, 1);
				failed+=check(checks, "terminal-target", assembler, word(paths[0], 0, k), word(paths[0], 5, k), 100, 5);
				failed+=check(checks, "disconnected-counted-target", assembler, word(paths[0], 0, k), word(paths[1], 0, k), 100, -1);
				failed+=check(checks, "identity-with-tail", assembler, word(paths[0], 0, k), word(paths[0], 0, k), 100, 0);
				failed+=check(checks, "terminal-identity", assembler, word(paths[0], 5, k), word(paths[0], 5, k), 100, 0);
				failed+=check(checks, "insufficient-limit", assembler, word(paths[0], 0, k), word(paths[0], 5, k), 2, -1);
				failed+=check(checks, "exclusive-limit", assembler, word(paths[0], 0, k), word(paths[0], 1, k), 1, -1);
				failed+=check(checks, "zero-limit", assembler, word(paths[0], 0, k), word(paths[0], 0, k), 0, -1);
				Files.write(root.resolve("checks.tsv"), checks.toString().getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
				if(iterations>0){benchmark(assembler, paths, k, warmup, iterations, root);}
				System.out.println("INSERT_DISTANCE k="+k+" checks=8 failed="+failed);
				if(failed>0){throw new AssertionError("Insert-distance contract failures; inspect "+root);}
			}finally{assembler.clearData();}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Measures a repeatable mixed query loop after warmup; functional failures are still reported afterward. */
	private static void benchmark(final Tadpole assembler, final byte[][] paths, final int k,
			final long warmup, final long iterations, final Path root) throws Exception{
		final Kmer[] from=new Kmer[5], reached=new Kmer[5];
		for(int i=0; i<5; i++){from[i]=word(paths[0], i, k); reached[i]=word(paths[0], i+1, k);}
		final Kmer disconnected=word(paths[1], 0, k), scratch=new Kmer(k);
		final int[] counts=new int[4];
		final long warmupSum=measureMixed(assembler, from, reached, disconnected, scratch, counts, warmup);
		final long start=System.nanoTime();
		final long sum=measureMixed(assembler, from, reached, disconnected, scratch, counts, iterations);
		final long nanos=System.nanoTime()-start;
		final String result="k\titerations\twarmup_iterations\tnanos\tchecksum\twarmup_checksum\n"+
				k+"\t"+iterations+"\t"+warmup+"\t"+nanos+"\t"+sum+"\t"+warmupSum+"\n";
		Files.write(root.resolve("benchmark.tsv"), result.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
		System.out.println("INSERT_BENCHMARK k="+k+" iterations="+iterations+" nanos="+nanos+" checksum="+sum);
	}

	/** Alternates reachable and disconnected targets at five starts, with no per-query allocation.
	 * Each correct pair contributes one plus minus one, hence a zero checksum.
	 * The long walker receives reusable scratch; templates and target words stay unchanged. */
	private static long measureMixed(final Tadpole assembler, final Kmer[] from, final Kmer[] reached,
			final Kmer disconnected, final Kmer scratch, final int[] counts, final long iterations){
		long sum=0;
		for(long i=0; i<iterations; i++){
			final int index=(int)((i>>>1)%5);
			final Kmer source=from[index], target=(i&1)==0 ? reached[index] : disconnected;
			if(assembler instanceof Tadpole1){
				sum+=((Tadpole1)assembler).measureInsert(source.array1()[0], source.array2()[0],
						target.array1()[0], target.array2()[0], 100, counts);
			}else{
				scratch.setFrom(source);
				sum+=((Tadpole2)assembler).measureInsert(scratch, target, 100, counts);
			}
		}
		return sum;
	}

	/** Verifies a unique depth-four chain and its terminal absence without invoking a walker. */
	private static void verifyPath(final Tadpole assembler, final byte[] path, final int k){
		for(int start=0; start<=5; start++){
			if(count(assembler, word(path, start, k))!=4){throw new IllegalStateException("Fixture words must be unique at depth four.");}
			final byte[] next=new byte[k];
			System.arraycopy(path, start+1, next, 0, k-1);
			for(byte base : ALPHABET){
				next[k-1]=base;
				final int expected=start<5 && base==path[start+k] ? 4 : 0;
				if(Math.max(0, count(assembler, word(next, 0, k)))!=expected){
					throw new IllegalStateException("Fixture has an unexpected successor at "+start+" base="+(char)base);
				}
			}
		}
	}

	/** Queries the actual table directly; absent-key sentinels remain visible to the caller. */
	private static int count(final Tadpole assembler, final Kmer key){
		return assembler instanceof Tadpole1 ? ((Tadpole1)assembler).getCount(key.array1()[0], key.array2()[0]) :
				((Tadpole2)assembler).tables().getCount(key);
	}

	/** Builds a fresh oriented query, since the long walker mutates its first query. */
	private static Kmer word(final byte[] bases, final int start, final int k){
		final Kmer key=new Kmer(k);
		if(key.kbig!=k){throw new IllegalStateException("Fixture requires exact packed K.");}
		for(int i=start; i<start+k; i++){key.addRight(bases[i]);}
		assert(key.len()==k) : "Each fixture word contains exactly K defined bases.";
		return key;
	}

	/** Records a desired result independently of assertions; no exception is counted as correct output. */
	private static int check(final StringBuilder checks, final String name, final Tadpole assembler,
			final Kmer from, final Kmer to, final int limit, final int expected){
		final int[] counts=new int[4];
		final int actual=assembler instanceof Tadpole1 ? ((Tadpole1)assembler).measureInsert(
				from.array1()[0], from.array2()[0], to.array1()[0], to.array2()[0], limit, counts) :
				((Tadpole2)assembler).measureInsert(from, to, limit, counts);
		checks.append(name).append('\t').append(expected).append('\t').append(actual).append('\t').append(actual==expected).append('\n');
		return actual==expected ? 0 : 1;
	}

	/** Literal paths proposed by Jean; their lengths and all successor counts are checked before use. */
	private static final String[] SHORT_PATHS={"GAACGGTATCCGGTAGCGAATCCGAAAATGCGACAC", "TGCTCTTTATTACCCAATTGGGAGAACAAATATTCG"};
	private static final String[] LONG_PATHS={"TAACCCTCTCGATCCCGATTGGGGGGCTACGTGACATATCTTCTACTGTGGAGCAGCTGTCCGCTACC",
			"GATTGCCTGCGTCATAAGGACTAAGCGTTCTGACTTTGGCCTGCATGTGAGCTGAGTCCGAGCTACAC"};
	private static final byte[] ALPHABET={'A', 'C', 'G', 'T'};
}
