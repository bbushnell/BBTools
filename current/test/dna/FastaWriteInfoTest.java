package test.dna;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import dna.ChromosomeArray;
import dna.Data;
import dna.FastaToChromArrays2;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Shared;

/** Finite metadata-completion checks for the legacy writeInfo entry point.
 * Requires a fresh JVM and work directory. Retains a tiny serialized chromosome
 * and outputs; process exit confines static Data/converter state. info-failure
 * tests retained raw-file errors. producer-failure loads a serialized chromosome
 * whose nonempty bounds exceed its backing array, then observes writer ownership
 * before rescuing only its newly created info.txt writer. No reflection or public
 * API changes are needed. A distinct watchdog exit97 is always unexpected failure.
 * @author Sucrose
 */
public final class FastaWriteInfoTest{

	/** Runs case=valid|info-failure|producer-failure with {@code work=<fresh-directory>}. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path work=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("work=")){work=Paths.get(arg.substring(5));}
			else{throw new IllegalArgumentException(arg);}
		}
		require(work!=null && ("valid".equals(mode) || "info-failure".equals(mode) || "producer-failure".equals(mode)),
				"Require case=valid|info-failure|producer-failure and work=");
		Files.createDirectory(work);
		Shared.setThreads(1);
		Data.setPath(work.toString());
		Data.CHROMGZ=false;
		FastaToChromArrays2.MID_PADDING=7;
		final Path genome=Paths.get(Data.ROOT_GENOME, "17");
		Files.createDirectories(genome);
		final Path input=Paths.get(Data.chromFname(1, 17));
		final ChromosomeArray original=new ChromosomeArray(1, Shared.PLUS);
		original.set(0, "NACGTN");
		original.resize(6);
		ReadWrite.write(original, input.toString(), false);
		final ChromosomeArray decoded=ReadWrite.read(ChromosomeArray.class, input.toString(), false);
		require(validChromosome(decoded), "Serialized fixture must independently round-trip its six literal bases");
		final Path info=genome.resolve("info.txt"), summary=genome.resolve("summary.txt");
		require(!Files.exists(info) && !Files.exists(summary), "Fresh fixture must not contain metadata");
		if("valid".equals(mode)){
			FastaToChromArrays2.writeInfo(17, 1, "tiny", "fixture", false, true);
			requireLoaded();
			equalMetadata(info, Arrays.asList("#Chromosome sizes", "#Version\t5",
					"#chrom\tscaffolds\tcontigs\tlength\tdefined\tundefined\tstartPad\tstopPad",
					"1\t1\t1\t6\t4\t2\t1\t1"));
			equalMetadata(summary, Arrays.asList("#Summary", "#Version\t5", "chroms\t1", "bases\t6",
					"defined\t4", "undefined\t2", "contigs\t1", "scaffolds\t1", "interpad\t7",
					"name\ttiny", "source\tfixture", "scafprefixes\ttrue"));
		}else if("info-failure".equals(mode)){infoFailure(info, summary, work);}
		else{producerFailure(input, info, summary, work);}
		System.out.println("PASS "+mode);
	}

	/** Creates a readable serialized object with invalid storage bounds, outside the tested catch.
	 * The leading N satisfies Data.getChromosome's assertion in both assertion modes. */
	private static void producerFailure(Path input, Path info, Path summary, Path work) throws Exception{
		final ChromosomeArray malformed=new ChromosomeArray(1, Shared.PLUS);
		malformed.array=new byte[]{'N'};
		malformed.minIndex=0;
		malformed.maxIndex=5;
		ReadWrite.write(malformed, input.toString(), false);
		try(ObjectInputStream in=new ObjectInputStream(Files.newInputStream(input))){
			final Object decoded=in.readObject();
			require(decoded instanceof ChromosomeArray && malformedChromosome((ChromosomeArray)decoded),
					"Independent deserialization must preserve the deliberately malformed chromosome");
		}
		final byte[] savedInput=Files.readAllBytes(input);
		final Set<Thread> before=Thread.getAllStackTraces().keySet();
		startWatchdog(work);
		RuntimeException failure=null;
		int liveWriters=-1;
		boolean loaded=false;
		try{
			try{FastaToChromArrays2.writeInfo(17, 1, "tiny", "fixture", false, true);}
			catch(RuntimeException problem){failure=problem;}
			// Inspect the cache directly; do not mask a skipped load by calling its getter.
			loaded=Data.GENOME_BUILD==17 && Data.chromosomePlusMatrix!=null &&
					Data.chromosomePlusMatrix.length==2 && malformedChromosome(Data.chromosomePlusMatrix[1]);
			liveWriters=ownedInfoWriters(before, info.toString()).size();
			System.out.println("malformedLoaded="+loaded+" exception="+(failure==null ? "none" : failure.getClass().getSimpleName())+
					" liveInfoWriters="+liveWriters+" summaryExists="+Files.exists(summary));
		}finally{
			rescueInfoWriters(before, info.toString());
		}
		require(loaded && failure instanceof ArrayIndexOutOfBoundsException &&
				hasFrame(failure, "dna.ChromosomeArray", "get") &&
				hasFrame(failure, "dna.FastaToChromArrays2", "countInfo") &&
				hasFrame(failure, "dna.FastaToChromArrays2", "writeInfo"),
				"Expected the original counting exception after malformed chromosome loading: "+failure);
		require(Arrays.equals(savedInput, Files.readAllBytes(input)), "Metadata generation changed its serialized input");
		require(!Files.exists(summary), "Producer failure published a new completion summary");
		equalMetadata(info, Arrays.asList("#Chromosome sizes", "#Version\t5",
				"#chrom\tscaffolds\tcontigs\tlength\tdefined\tundefined\tstartPad\tstopPad"));
		System.out.println("producerExceptionVerified=true inputPreserved=true summaryAbsent=true");
		if(liveWriters!=0){
			System.out.println("FAIL_PRODUCER_CLEANUP liveInfoWriters="+liveWriters);
			throw new AssertionError("writeInfo must join its started writer before propagating a producer exception");
		}
	}

	/** Exact malformed object shape which is readable but fails at the second counted base. */
	private static boolean malformedChromosome(ChromosomeArray chromosome){
		return chromosome!=null && chromosome.chromosome==1 && chromosome.strand==Shared.PLUS &&
				chromosome.minIndex==0 && chromosome.maxIndex==5 && Arrays.equals(chromosome.array, new byte[]{'N'});
	}

	/** Checks direct stack frames rather than accepting a wrapper with the intended failure as its cause. */
	private static boolean hasFrame(Throwable failure, String type, String method){
		for(StackTraceElement frame : failure.getStackTrace()){
			if(frame.getClassName().equals(type) && frame.getMethodName().equals(method)){return true;}
		}
		return false;
	}

	/** Selects only live info writers created after the snapshot at this fixture's exact destination. */
	private static ArrayList<TextStreamWriter> ownedInfoWriters(Set<Thread> before, String info){
		final ArrayList<TextStreamWriter> owned=new ArrayList<TextStreamWriter>();
		for(Thread thread : Thread.getAllStackTraces().keySet()){
			if(!before.contains(thread) && thread instanceof TextStreamWriter && thread.isAlive()){
				final TextStreamWriter writer=(TextStreamWriter)thread;
				if(info.equals(writer.fname)){owned.add(writer);}
			}
		}
		return owned;
	}

	/** Rescues only this call's writer after the producer unwinds, preserving pre-rescue observations.
	 * Tiny ordinary buffers fit the queue; watchdog97 covers unexpected poison/join nontermination. */
	private static void rescueInfoWriters(Set<Thread> before, String info) throws InterruptedException{
		final ArrayList<TextStreamWriter> owned=ownedInfoWriters(before, info);
		for(TextStreamWriter writer : owned){writer.poison();}
		for(TextStreamWriter writer : owned){writer.join(5000);}
		for(TextStreamWriter writer : owned){
			require(!writer.isAlive(), "Fixture rescue could not join its info writer");
			require(!writer.errorState, "Fixture rescue encountered unexpected info output failure");
		}
		require(ownedInfoWriters(before, info).isEmpty(), "Fixture rescue left an owned info writer alive");
		System.out.println("rescuedInfoWriters="+owned.size());
	}

	/** Remains armed until this fresh JVM exits; watchdog failure is never an expected baseline result. */
	private static void startWatchdog(final Path work){
		final Thread watchdog=new Thread(new Runnable(){
			@Override
			public void run(){
				try{Thread.sleep(30000);}
				catch(InterruptedException interrupted){Thread.currentThread().interrupt();}
				try{
					Files.write(work.resolve("watchdog-failure.txt"), "WATCHDOG_FAILURE exit=97\n".getBytes(StandardCharsets.US_ASCII));
				}catch(IOException failure){
					System.err.println("WATCHDOG_FAILURE: could not save marker: "+failure);
				}finally{
					System.err.println("WATCHDOG_FAILURE exit=97: writeInfo test did not terminate normally");
					Runtime.getRuntime().halt(97);
				}
			}
		}, "FastaWriteInfoTest-watchdog");
		watchdog.setDaemon(true);
		watchdog.start();
	}

	/** Verifies the device independently, then requires evidence from the actual metadata writer. */
	private static void infoFailure(Path info, Path summary, Path work) throws Exception{
		final Path full=Paths.get("/dev/full");
		require(Files.exists(full), "This controlled output failure requires Linux /dev/full");
		Files.createSymbolicLink(info, full);
		require(Files.isSameFile(info, full), "info.txt must resolve to the failing device");
		boolean rejected=false;
		try(OutputStream probe=Files.newOutputStream(info)){
			try{probe.write(1);}
			catch(IOException expected){rejected=true;}
		}
		require(rejected, "Independent raw write must fail before invoking the converter");

		final PrintStream previous=System.err;
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		RuntimeException failure=null;
		try(PrintStream capture=new PrintStream(bytes, true, "UTF-8")){
			System.setErr(capture);
			try{FastaToChromArrays2.writeInfo(17, 1, "tiny", "fixture", false, true);}
			catch(RuntimeException problem){failure=problem;}
			finally{System.setErr(previous);}
		}
		final String diagnostics=new String(bytes.toByteArray(), StandardCharsets.UTF_8);
		previous.print(diagnostics);
		Files.write(work.resolve("writer.stderr"), bytes.toByteArray());
		requireLoaded();
		final boolean writerFailed=diagnostics.contains("java.io.IOException") &&
				diagnostics.contains("fileIO.TextStreamWriter.run");
		final boolean linked=Files.isSymbolicLink(info) && Files.isSameFile(info, full);
		final boolean summaryExists=Files.exists(summary);
		System.out.println("deviceRejects="+rejected+" writerFailureObserved="+writerFailed+
				" chromosomeLoaded=true targetPreserved="+linked+" exception="+
				(failure==null ? "none" : failure.getClass().getSimpleName())+" summaryExists="+summaryExists);
		require(writerFailed && linked, "Actual info writer must reach the verified failing target");
		require(failure!=null, "writeInfo must report its metadata write failure");
		require(!summaryExists, "Failed info output must not publish a new completion summary");
	}

	/** Checks the cache directly so the fixture cannot itself load a chromosome skipped by writeInfo. */
	private static void requireLoaded(){
		require(Data.GENOME_BUILD==17 && Data.chromosomePlusMatrix!=null &&
				Data.chromosomePlusMatrix.length==2 && validChromosome(Data.chromosomePlusMatrix[1]),
				"writeInfo must have loaded the serialized chromosome before the output verdict");
	}

	private static boolean validChromosome(ChromosomeArray chromosome){
		return chromosome!=null && chromosome.chromosome==1 && chromosome.minIndex==0 && chromosome.maxIndex==5 &&
				Arrays.equals(chromosome.array, "NACGTN".getBytes(StandardCharsets.US_ASCII));
	}

	/** Removes only the one variable generated-date line, comparing every other line to literal expected data. */
	private static void equalMetadata(Path path, List<String> expected) throws Exception{
		final ArrayList<String> actual=new ArrayList<String>();
		int dates=0;
		for(String line : Files.readAllLines(path, StandardCharsets.US_ASCII)){
			if(line.startsWith("#Generated on\t")){dates++;}
			else{actual.add(line);}
		}
		System.out.println(path.getFileName()+" generatedDates="+dates+" actual="+actual);
		require(dates==1 && expected.equals(actual), "Unexpected metadata in "+path.getFileName());
	}

	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
