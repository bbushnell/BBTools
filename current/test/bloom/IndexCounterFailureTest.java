package test.bloom;

import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Properties;
import java.util.concurrent.TimeUnit;

import bloom.IndexCounter;
import bloom.KCountArray;
import bloom.KCountArray3;
import bloom.KmerCountAbstract;
import dna.ChromosomeArray;
import dna.Data;
import shared.Shared;

/** Checks caller-visible rejection of real IndexCounter worker failures.
 * A single worker counts trimmed NACG, then an invalid declared chromosome count
 * causes a bounds failure before reference loading. Direct and factory entry points
 * retain separate control/fault cases. Partial tables are observed, never rolled back.
 * Each case runs in a fresh bounded child, inheriting the parent's assertion mode.
 * @author Brian Bushnell, Sucrose */
public final class IndexCounterFailureTest{

	/** Runs every child even when an earlier baseline call incorrectly returns success. */
	public static void main(final String[] args) throws Exception{
		if(args.length==2 && args[0].startsWith("child=") && args[1].startsWith("outdir=")){
			child(args[0].substring(6), Paths.get(args[1].substring(7)));
			return;
		}
		if(args.length!=1 || !args[0].startsWith("outdir=")){
			throw new IllegalArgumentException("Usage: outdir=<new directory>");
		}
		final Path root=Paths.get(args[0].substring(7)).toAbsolutePath();
		check(!Files.exists(root), "Fresh output prevents stale evidence from satisfying a new run.");
		Files.createDirectories(root);
		int failures=0;
		for(String name : CASES){failures+=runChild(root, name);}
		if(failures>0){throw new AssertionError("IndexCounter failure contracts failed "+failures+" of "+CASES.length+" cases.");}
		System.out.println("INDEX_COUNTER_FAILURE_TEST_PASS cases="+CASES.length);
	}

	/** Terminates only an overdue fixture-owned process; forced exit never passes. */
	private static int runChild(final Path root, final String name) throws Exception{
		final Path dir=Files.createDirectory(root.resolve(name));
		final ArrayList<String> command=new ArrayList<String>();
		command.add(Paths.get(System.getProperty("java.home"), "bin", "java").toString());
		command.add(IndexCounterFailureTest.class.desiredAssertionStatus() ? "-ea" : "-da");
		command.add("-Xmx256m");
		command.add("-XX:+ExitOnOutOfMemoryError");
		try{
			Class.forName("jdk.incubator.vector.ByteVector", false, IndexCounterFailureTest.class.getClassLoader());
			command.add("--add-modules=jdk.incubator.vector");
		}catch(ClassNotFoundException | LinkageError unavailable){
			//Do not request an optional module unavailable to the parent.
		}
		command.add("-cp");
		command.add(System.getProperty("java.class.path"));
		command.add(IndexCounterFailureTest.class.getName());
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
				Boolean.toString(IndexCounterFailureTest.class.desiredAssertionStatus()).equals(result.getProperty("assertions"));
		final String summary="case="+name+" completed="+completed+" exit="+exit+" result="+result+" pass="+pass;
		Files.write(dir.resolve("outcome.txt"), Arrays.asList(summary), StandardCharsets.UTF_8);
		System.out.println(summary);
		return pass ? 0 : 1;
	}

	/** Observes failure and partial work before checking the fixed rejection contract. */
	private static void child(final String name, final Path dir) throws Exception{
		check(Arrays.asList(CASES).contains(name), "Unknown child case: "+name);
		final boolean factory=name.startsWith("factory");
		final boolean fault=name.endsWith("fault");
		final Snapshot saved=new Snapshot();
		final Properties result=new Properties();
		final ArrayList<String> failures=new ArrayList<String>();
		final ChromosomeArray chromosome=new ChromosomeArray(1, Shared.PLUS);
		chromosome.set(0, "NACG");
		chromosome.resize(chromosome.maxIndex+1);
		check(chromosome.minIndex==0 && chromosome.maxIndex==3 && chromosome.array.length==4 &&
				Arrays.equals(chromosome.array, BASES),
				"NACG must occupy exactly four indexed bases, including Data.getChromosome's leading N.");
		Files.write(dir.resolve("input.txt"), BASES);
		Files.write(dir.resolve("oracle.txt"), Arrays.asList("forwardKey=6", "observedDepth=1", "occupiedCells=1",
				"callerMustReject="+fault, "tableMode="+(factory ? "hashed-presence" : "exact-direct"),
				"failedCounters=diagnostic-only", "rollback=not-promised"), StandardCharsets.UTF_8);
		final ObservedFactoryCounter observer=factory ? new ObservedFactoryCounter() : null;
		KCountArray table=null, returned=null;
		Throwable callerFailure=null, observationFailure=null;
		try{
			Shared.setThreads(1);
			Data.numChroms=fault ? 2 : 1;
			Data.chromosomePlusMatrix=new ChromosomeArray[]{null, chromosome};
			KmerCountAbstract.CANONICAL=true;
			KmerCountAbstract.maxShortKmerLength=31;
			KmerCountAbstract.readsProcessed=KmerCountAbstract.keysCounted=KmerCountAbstract.increments=0;
			check(Shared.threads()==1 && Data.chromosomePlusMatrix.length==2,
					"One real worker must finish chromosome1 before chromosome2 fails bounds, without loading a null slot.");
			if(!factory){table=new KCountArray3(2048, 4);}
			try{
				returned=factory ? observer.makeKcaFromIndex(4096, 1, 1) : new IndexCounter(3, false).countFromIndex(table);
			}catch(Throwable caught){callerFailure=caught;}
			if(factory){table=observer.observed;}
			result.setProperty("actualReadsProcessed", Long.toString(KmerCountAbstract.readsProcessed));
			result.setProperty("actualKeysCounted", Long.toString(KmerCountAbstract.keysCounted));
			result.setProperty("actualIncrements", Long.toString(KmerCountAbstract.increments));
			if(table==null){
				failures.add("No actual table retained for partial-state observation.");
			}else{
				observeTable(table, factory, dir, result, failures);
			}
			if(fault){
				if(callerFailure==null || returned!=null){failures.add("Worker failure must reach the caller, not return a successful partial table.");}
				final Class<? extends Throwable> expected=IndexCounterFailureTest.class.desiredAssertionStatus() ?
						AssertionError.class : ArrayIndexOutOfBoundsException.class;
				final boolean cause=hasCause(callerFailure, expected);
				result.setProperty("expectedCause", expected.getName());
				result.setProperty("causePreserved", Boolean.toString(cause));
				if(!cause){failures.add("Caller failure chain must retain the real chromosome bounds/assertion cause.");}
			}else if(callerFailure!=null || returned==null || returned!=table){
				failures.add("The ordinary one-chromosome control must return its actual populated table.");
			}
		}catch(Throwable caught){observationFailure=caught; failures.add("Setup/observation failed: "+caught);}
		finally{saved.restore();}
		final boolean restored=saved.matches();
		final boolean preserved=Arrays.equals(chromosome.array, BASES) && Arrays.equals(Files.readAllBytes(dir.resolve("input.txt")), BASES);
		if(!restored){failures.add("All saved globals, including the buffer count after restoring threads, must be restored.");}
		if(!preserved){failures.add("Installed chromosome and saved input bytes must remain unchanged.");}
		result.setProperty("case", name);
		result.setProperty("assertions", Boolean.toString(IndexCounterFailureTest.class.desiredAssertionStatus()));
		result.setProperty("configured", "k=3 rcomp=false canonical=true threads=1 matrixLength=2 chromosome1=NACG numChroms="+(fault ? 2 : 1));
		result.setProperty("globalsRestored", Boolean.toString(restored));
		result.setProperty("inputPreserved", Boolean.toString(preserved));
		result.setProperty("callerRejected", Boolean.toString(callerFailure!=null));
		result.setProperty("tableReturned", Boolean.toString(returned!=null));
		result.setProperty("callerFailure", callerFailure==null ? "none" : callerFailure.toString());
		result.setProperty("pass", Boolean.toString(failures.isEmpty()));
		writeFailure(dir.resolve("caller-exception.txt"), callerFailure);
		writeFailure(dir.resolve("observation-exception.txt"), observationFailure);
		Files.write(dir.resolve("failures.txt"), failures, StandardCharsets.UTF_8);
		try(OutputStream output=Files.newOutputStream(dir.resolve("result.properties"))){
			result.store(output, "Real worker/caller failure observations; partial state is not rollback");
		}
		if(!failures.isEmpty()){throw new AssertionError("IndexCounter failure contract: "+name+": "+failures);}
		//Natural return exposes leaked non-daemon workers; do not force successful process exit.
	}

	/** Saves physical cells and checks the one completed NACG window without assuming absent hashed keys. */
	private static void observeTable(final KCountArray table, final boolean factory, final Path dir,
			final Properties result, final ArrayList<String> failures) throws Exception{
		Files.write(dir.resolve("table-cells.txt"), table.toContentsString().toString().getBytes(StandardCharsets.UTF_8));
		final int depth=table.read(6);
		final long occupied=table.cellsUsed(1);
		final long[] histogram=table.transformToFrequency();
		final ArrayList<String> rows=new ArrayList<String>();
		rows.add("depth\tcells");
		for(int i=0; i<histogram.length; i++){if(histogram[i]!=0){rows.add(i+"\t"+histogram[i]);}}
		Files.write(dir.resolve("histogram.tsv"), rows, StandardCharsets.UTF_8);
		result.setProperty("depth6", Integer.toString(depth));
		result.setProperty("occupiedCells", Long.toString(occupied));
		result.setProperty("allocatedCells", Long.toString(table.cells));
		result.setProperty("cellBits", Integer.toString(table.cellBits));
		result.setProperty("histogramBins", Integer.toString(histogram.length));
		if(depth!=1 || occupied!=1 || histogram.length<2 || histogram[1]!=1){
			failures.add("Chromosome1 must leave exactly one occupied depth1 cell for forward ACG key6.");
		}
		if(!factory){
			check(table.cells==2048 && table.cellBits==4, "The exact direct oracle covers all2048 four-bit cells.");
			for(int key=0; key<2048; key++){
				if(table.read(key)!=(key==6 ? 1 : 0)){failures.add("Unexpected exact direct count at key "+key);}
			}
		}else if(table.cellBits!=1){failures.add("Factory control must retain its presence-bitmap contract.");}
	}

	/** Retains a factory-owned table but delegates all worker behavior unchanged to the real counter. */
	private static final class ObservedFactoryCounter extends IndexCounter{
		ObservedFactoryCounter(){super(3, false);}
		@Override
		public KCountArray countFromIndex(final KCountArray counts) throws Exception{
			observed=counts;
			return super.countFromIndex(counts);
		}
		private KCountArray observed;
	}

	/** Saves the exact shared state touched by this fixture, without resetting unrelated index metadata. */
	private static final class Snapshot{
		void restore(){
			Data.numChroms=numChroms;
			Data.chromosomePlusMatrix=chromosomes;
			KmerCountAbstract.CANONICAL=canonical;
			KmerCountAbstract.maxShortKmerLength=maxShort;
			KmerCountAbstract.readsProcessed=reads;
			KmerCountAbstract.keysCounted=keys;
			KmerCountAbstract.increments=increments;
			Shared.setThreads(threads);
			Shared.setBuffers(buffers);
		}
		boolean matches(){
			return Data.numChroms==numChroms && Data.chromosomePlusMatrix==chromosomes &&
					KmerCountAbstract.CANONICAL==canonical && KmerCountAbstract.maxShortKmerLength==maxShort &&
					KmerCountAbstract.readsProcessed==reads && KmerCountAbstract.keysCounted==keys &&
					KmerCountAbstract.increments==increments && Shared.threads()==threads &&
					Shared.READ_BUFFER_NUM_BUFFERS()==buffers;
		}
		private final int threads=Shared.threads(), buffers=Shared.READ_BUFFER_NUM_BUFFERS();
		private final int numChroms=Data.numChroms, maxShort=KmerCountAbstract.maxShortKmerLength;
		private final ChromosomeArray[] chromosomes=Data.chromosomePlusMatrix;
		private final boolean canonical=KmerCountAbstract.CANONICAL;
		private final long reads=KmerCountAbstract.readsProcessed, keys=KmerCountAbstract.keysCounted,
				increments=KmerCountAbstract.increments;
	}

	/** Requires the real lookup failure in the caller chain, not an unrelated assertion or worker stderr alone. */
	private static boolean hasCause(Throwable failure, final Class<? extends Throwable> type){
		for(; failure!=null; failure=failure.getCause()){
			if(type.isInstance(failure)){
				for(StackTraceElement frame : failure.getStackTrace()){
					if(frame.getClassName().equals("dna.Data") && frame.getMethodName().equals("getChromosome")){return true;}
				}
			}
		}
		return false;
	}

	/** Retains a complete caller/observation chain, including suppressed failures, or an explicit absence. */
	private static void writeFailure(final Path path, final Throwable failure) throws Exception{
		try(PrintStream output=new PrintStream(Files.newOutputStream(path), true, "UTF-8")){
			if(failure==null){
				output.println("none");
			}else{
				failure.printStackTrace(output);
			}
		}
	}

	/** Keeps fixture setup checks active in both assertion modes. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final byte[] BASES="NACG".getBytes(StandardCharsets.US_ASCII);
	private static final String[] CASES={"direct-control", "direct-fault", "factory-control", "factory-fault"};
}
