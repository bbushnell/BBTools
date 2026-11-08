package test.bloom;

import java.io.InputStream;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Properties;
import java.util.concurrent.TimeUnit;

import bloom.KCountArray;
import bloom.ReadCounter;
import kmer.KmerTableSet;
import shared.Shared;
import stream.FASTQ;
import stream.Read;
import ukmer.Kmer;

/** Checks recovery of long-window quality eligibility after a low-quality prefix.
 * Each child counts one real FASTQ record using the public ReadCounter API.
 * Literal expected counts are checked against independently recomputed double
 * window products, never the production rolling float/inverse recurrence.
 * Includes a precision-only boundary case whose fixed oracle is zero windows.
 * Parent invocation selects the assertion mode inherited by all fresh children.
 * @author Brian Bushnell, Sucrose */
public final class ReadCounterProbabilityTest{

	/** Runs all cases even if an earlier case fails; child selectors are internal. */
	public static void main(final String[] args) throws Exception{
		if(args.length==2 && args[0].startsWith("child=") && args[1].startsWith("outdir=")){
			child(args[0].substring(6), Paths.get(args[1].substring(7)));
			return;
		}
		if(args.length!=1 || !args[0].startsWith("outdir=")){
			throw new IllegalArgumentException("Usage: outdir=<new directory>");
		}
		final Path root=Paths.get(args[0].substring(7)).toAbsolutePath();
		check(!Files.exists(root), "Fresh output prevents stale evidence from satisfying a new capture.");
		Files.createDirectories(root);
		int failures=0;
		for(String name : CASES){failures+=runChild(root, name);}
		if(failures>0){throw new AssertionError("ReadCounter probability contracts failed "+failures+" of "+CASES.length+" cases.");}
		System.out.println("READ_COUNTER_PROBABILITY_TEST_PASS cases="+CASES.length);
	}

	/** Bounds only the child this fixture owns; a killed child never counts as cleanup success. */
	private static int runChild(final Path root, final String name) throws Exception{
		final Path dir=Files.createDirectory(root.resolve(name));
		final ArrayList<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.add(ReadCounterProbabilityTest.class.desiredAssertionStatus() ? "-ea" : "-da");
		command.add("-Xmx256m");
		command.add("-XX:+ExitOnOutOfMemoryError");
		try{
			Class.forName("jdk.incubator.vector.ByteVector", false, ReadCounterProbabilityTest.class.getClassLoader());
			command.add("--add-modules=jdk.incubator.vector");
		}catch(ClassNotFoundException | LinkageError unavailable){
			//Older or unresolved-module parents use the scalar path; do not request an absent module.
		}
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(ReadCounterProbabilityTest.class.getName());
		command.add("child="+name);
		command.add("outdir="+dir);
		Files.write(dir.resolve("command.txt"), command, StandardCharsets.UTF_8);
		final Process process=new ProcessBuilder(command).redirectInput(Paths.get("/dev/null").toFile())
				.redirectOutput(dir.resolve("stdout.txt").toFile()).redirectError(dir.resolve("stderr.txt").toFile()).start();
		boolean completed=false;
		try{completed=process.waitFor(15, TimeUnit.SECONDS);}
		finally{
			if(process.isAlive()){
				process.destroy();
				if(!process.waitFor(2, TimeUnit.SECONDS)){
					process.destroyForcibly();
					check(process.waitFor(5, TimeUnit.SECONDS), "Fixture-owned child must terminate after forced destruction.");
				}
			}
		}
		final Properties result=new Properties();
		if(Files.exists(dir.resolve("result.properties"))){
			try(InputStream input=Files.newInputStream(dir.resolve("result.properties"))){result.load(input);}
		}
		final int exit=process.exitValue();
		final boolean pass=completed && exit==0 && "true".equals(result.getProperty("pass")) &&
				name.equals(result.getProperty("case")) &&
				Boolean.toString(ReadCounterProbabilityTest.class.desiredAssertionStatus()).equals(result.getProperty("assertions"));
		final String summary="case="+name+" completed="+completed+" exit="+exit+" result="+result+" pass="+pass;
		Files.write(dir.resolve("outcome.txt"), Arrays.asList(summary), StandardCharsets.UTF_8);
		System.out.println(summary);
		return pass ? 0 : 1;
	}

	/** Saves the literal input and observed counters/table state even when the expected count fails. */
	private static void child(final String name, final Path dir) throws Exception{
		check(name.equals("recovery") || name.equals("high-control") || name.equals("threshold-boundary"), "Unknown child case: "+name);
		final boolean recovery=name.equals("recovery");
		final boolean boundary=name.equals("threshold-boundary");
		final int expected=boundary ? 0 : recovery ? 1 : 129;
		final float cutoff=boundary ? 1f-8388f/16777216f : 0.5f;
		final byte[] bases=new byte[boundary ? K : LENGTH], qualities=new byte[bases.length];
		Arrays.fill(bases, (byte)'A');
		Arrays.fill(qualities, (byte)40);
		if(recovery){Arrays.fill(qualities, 0, K, (byte)2);}
		if(boundary){Arrays.fill(qualities, 0, 123, (byte)93);}
		final byte[] encoded=qualities.clone();
		for(int i=0; i<encoded.length; i++){encoded[i]+=33;}
		final byte[] input=("@"+name+"\n"+new String(bases, StandardCharsets.US_ASCII)+"\n+\n"+
				new String(encoded, StandardCharsets.US_ASCII)+"\n").getBytes(StandardCharsets.US_ASCII);
		final Path path=dir.resolve("input.fq");
		final Properties result=new Properties();
		KCountArray table=null;
		long rawKey=0;
		Throwable failure=null;
		try{
			configure();
			if(boundary){Read.setMaxCalledQuality(93); ReadCounter.minProb=cutoff;}
			Files.write(path, input);
			result.setProperty("cutoff", Float.toString(cutoff));
			final double low=KmerTableSet.PROB_CORRECT[2], high=KmerTableSet.PROB_CORRECT[40];
			result.setProperty("probCorrectQ2", Double.toString(low));
			result.setProperty("probCorrectQ40", Double.toString(high));
			check(low>0 && low<0.5 && high>0.9998 && high<=1,
					"Q2 must reject any mixed window; 128 Q40 bases must stay safely above probability 0.5.");
			check(Read.capQuality((byte)2)==2 && Read.capQuality((byte)40)==40,
					"Normal defined-base validation must preserve the fixture Q2/Q40 values.");
			if(boundary){
				check(Read.capQuality((byte)93)==93 && KmerTableSet.PROB_CORRECT[93]==1f &&
						KmerTableSet.PROB_CORRECT[40]==1f-1678f/16777216f,
						"Boundary proof requires actual Q93 factor1 and Q40 factor1-1678/2^24, without quality capping.");
				result.setProperty("probCorrectQ93", Float.toString(KmerTableSet.PROB_CORRECT[93]));
			}
			final int oracle=eligibleWindows(qualities, cutoff, dir);
			result.setProperty("oracleCount", Integer.toString(oracle));
			check(oracle==expected, "Independent per-window products must match the literal count "+expected+", got "+oracle);
			final Kmer key=new Kmer(K);
			check(key.kbig==K && key.k==32 && !key.hollow && Kmer.getMult(K)==4,
					"Packed k128 must represent exactly four full 32-base words, without rounding or hollow hashing.");
			for(int i=0; i<4; i++){check(key.perWordK(i)==32, "Every k128 word must span 32 bases.");}
			for(int i=0; i<K; i++){key.addRight((byte)'A');}
			rawKey=key.xor();//Address only: this hash does not determine quality eligibility.
			result.setProperty("kbig", Integer.toString(key.kbig));
			result.setProperty("wordCount", Integer.toString(Kmer.getMult(K)));
			result.setProperty("rawKey", Long.toString(rawKey));
			table=KCountArray.makeNew(4096, 8, 1);
			check(table.cellBits==8 && table.maxValue==255, "The largest expected depth129 must not saturate this table.");
			final ReadCounter counter=new ReadCounter(K, true, false, false, false, false);
			final KCountArray returned=counter.count(path.toString(), null, table);
			check(returned==table, "The public count overload must return the borrowed destination.");
		}catch(Throwable caught){failure=caught;}
		try{
			if(table!=null){
				table.shutdown();//Own destination only; never close the counter's input on its behalf.
				final int depth=table.read(rawKey);
				final long occupied=table.cellsUsed(1);
				result.setProperty("actualDepth", Integer.toString(depth));
				result.setProperty("occupiedCells", Long.toString(occupied));
				result.setProperty("allocatedCells", Long.toString(table.cells));
				final long[] histogram=table.transformToFrequency();
				final ArrayList<String> rows=new ArrayList<String>();
				for(int i=0; i<histogram.length; i++){
					if(histogram[i]>0){rows.add(i+"\t"+histogram[i]);}
				}
				Files.write(dir.resolve("histogram.tsv"), rows, StandardCharsets.UTF_8);
				check(depth==expected && occupied==(expected==0 ? 0 : 1) && (expected==0 || histogram[expected]==1) &&
						ReadCounter.keysCounted==expected && ReadCounter.increments==expected && ReadCounter.readsProcessed==1,
						"Expected "+expected+" eligible windows, "+(expected==0 ? 0 : 1)+" occupied keys and one input read; depth="+depth+
						", occupied="+occupied+", keys="+ReadCounter.keysCounted+", increments="+ReadCounter.increments+
						", reads="+ReadCounter.readsProcessed);
			}
		}catch(Throwable observed){if(failure==null){failure=observed;}else{failure.addSuppressed(observed);}}
		try{check(Arrays.equals(input, Files.readAllBytes(path)), "Counter must not change the saved literal FASTQ input.");}
		catch(Throwable changed){if(failure==null){failure=changed;}else{failure.addSuppressed(changed);}}
		result.setProperty("keysCounted", Long.toString(ReadCounter.keysCounted));
		result.setProperty("increments", Long.toString(ReadCounter.increments));
		result.setProperty("readsProcessed", Long.toString(ReadCounter.readsProcessed));
		result.setProperty("case", name);
		result.setProperty("expectedCount", Integer.toString(expected));
		result.setProperty("assertions", Boolean.toString(ReadCounterProbabilityTest.class.desiredAssertionStatus()));
		result.setProperty("pass", Boolean.toString(failure==null));
		result.setProperty("failure", failure==null ? "none" : failure.toString());
		try(OutputStream output=Files.newOutputStream(dir.resolve("result.properties"))){result.store(output, "Observed long-window counts; literal oracle unchanged across baseline/fixed");}
		if(failure!=null){throw new AssertionError("Probability contract failed: "+name, failure);}
		//Natural return leaves a leaked non-daemon input worker visible to the parent deadline.
	}

	/** Pins process-global settings inside an isolated child; no other invocation shares them. */
	private static void configure(){
		Locale.setDefault(Locale.ROOT);
		Shared.setThreads(1);
		FASTQ.ASCII_OFFSET=33;
		FASTQ.DETECT_QUALITY=FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		Read.VALIDATE_IN_CONSTRUCTOR=Read.CHANGE_QUALITY=true;
		Read.setMinCalledQuality(2);
		Read.setMaxCalledQuality(50);
		Kmer.PACKED=Kmer.FULL_MIX=true;
		Kmer.MASK_CORE=false;
		ReadCounter.maxShortKmerLength=31;
		ReadCounter.MAX_COUNT_THREADS=1;
		ReadCounter.maxReads=-1;
		ReadCounter.minQuality=0;
		ReadCounter.minProb=0.5f;
		ReadCounter.KEEP_DUPLICATE_KMERS=ReadCounter.CANONICAL=true;
		ReadCounter.KMERS_PER_READ=ReadCounter.IDMASK=0;
		ReadCounter.BUFFERED=ReadCounter.SKETCH_MODE=ReadCounter.STORE_HASHED=ReadCounter.PREJOIN=false;
		check(ReadCounter.keysCounted==0 && ReadCounter.increments==0 && ReadCounter.readsProcessed==0,
				"Fresh children must begin with zero counting totals; do not mask state contamination by resetting them.");
	}

	/** Recomputes each complete window independently in double precision and saves all products. */
	private static int eligibleWindows(final byte[] qualities, final float cutoff, final Path dir) throws Exception{
		assert((qualities.length==LENGTH || qualities.length==K) && K==128) : "Literal oracle uses 129 recovery/control windows or one boundary window.";
		final ArrayList<String> rows=new ArrayList<String>();
		rows.add("start\tprobability\teligible");
		int count=0;
		for(int start=0; start+K<=qualities.length; start++){
			double probability=1;
			for(int pos=start; pos<start+K; pos++){probability*=(double)KmerTableSet.PROB_CORRECT[qualities[pos]];}
			check(probability>0 && Double.isFinite(probability), "Independent window products must not underflow or become nonfinite.");
			final boolean eligible=probability>=cutoff;
			if(eligible){count++;}
			rows.add(start+"\t"+probability+"\t"+eligible);
		}
		Files.write(dir.resolve("oracle.tsv"), rows, StandardCharsets.UTF_8);
		return count;
	}

	/** Keeps the literal oracle active with either JVM assertion setting. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final int K=128, LENGTH=256;
	private static final String[] CASES={"recovery", "high-control", "threshold-boundary"};
}
