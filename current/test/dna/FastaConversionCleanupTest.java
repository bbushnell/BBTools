package test.dna;

import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.lang.management.ManagementFactory;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import dna.ChromosomeArray;
import dna.FastaToChromArrays2;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.Shared;

/** Finite producer-failure probe for reference conversion resource ownership.
 * Run each case in a fresh JVM with a fresh work directory on Linux and
 * {@code -XX:+UnlockExperimentalVMOptions -XX:+UseEpsilonGC}. Non-collecting GC
 * prevents collection from closing an abandoned input before descriptor observation.
 * The test changes global converter/I/O settings; JVM exit confines those changes.
 *
 * Collision uses three records and disables merging so the input has not reached
 * empty EOF when the first existing chromosome file rejects overwrite. Observations
 * precede rescue. Only new TextStreamWriter instances at the two exact metadata
 * paths are poisoned and joined after the calling producer has unwound. An orphaned
 * input has no public close handle; process exit releases it. No reflection, unrelated
 * worker intervention, or production cleanup hook is used.
 *
 * A daemon watchdog remains armed through JVM exit. Unexpected nontermination writes
 * WATCHDOG_FAILURE and halts only this test JVM with exit 97. That is never an expected
 * baseline failure. A valid collision with leaked resources reports
 * FAIL_RESOURCE_CLEANUP only after successful writer rescue, then throws AssertionError.
 * @author Sucrose
 */
public final class FastaConversionCleanupTest{

	/** Runs {@code case=collision|write|metadata work=<fresh-directory>} once per JVM. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path work=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("work=")){work=Paths.get(arg.substring(5));}
			else{throw new IllegalArgumentException(arg);}
		}
		require(("collision".equals(mode) || "write".equals(mode) || "metadata".equals(mode)) && work!=null,
				"Require case=collision|write|metadata and a fresh work= directory");
		requireEpsilon();
		Files.createDirectory(work);
		work=work.toRealPath();
		final Path input=work.resolve("reference.fa"), genome=work.resolve("ref/genome/17");
		Files.write(input, FASTA.getBytes(StandardCharsets.US_ASCII));
		Files.createDirectories(genome);
		checkDescriptorObservation(input);
		Shared.setThreads(1);
		ReadWrite.RAWMODE=false;
		ReadWrite.USE_PIGZ=ReadWrite.USE_GZIP=false;
		ReadWrite.USE_UNPIGZ=ReadWrite.USE_GUNZIP=ReadWrite.USE_UNBGZIP=false;
		ReadWrite.ALLOW_NATIVE_BGZF=false;
		require(ReadWrite.countActiveThreads()==0, "Fresh JVM must not contain earlier object writers");
		final boolean collision=mode.equals("collision"), writeChroms=!mode.equals("metadata");
		final Path first=genome.resolve("chr1.chrom");
		if(collision){Files.write(first, SENTINEL);}
		final String info=genome.resolve("info.txt").toString();
		final String scaffolds=genome.resolve("scaffolds.txt.gz").toString();
		final Set<Thread> before=Thread.getAllStackTraces().keySet();
		startWatchdog(work);
		final String[] convert={"in="+input, "build=17", "path="+work, "name=tiny",
				"nodisk=f", "writechroms="+writeChroms, "retain=t", "genscaffoldinfo=t",
				"startpad=0", "midpad=0", "stoppad=0", "minscaf=1", "mergescaffolds=f",
				"overwrite=f", "gz=f", "waitforwriting=t", "scafprefixes=f"};
		ArrayList<ChromosomeArray> result=null;
		RuntimeException failure=null;
		int inputDescriptors=-1, liveWriters=-1;
		try{
			try{result=FastaToChromArrays2.main2(convert);}
			catch(RuntimeException problem){failure=problem;}
			// These are the tested observations, before any fixture cleanup.
			inputDescriptors=countInputDescriptors(input);
			liveWriters=ownedWriters(before, info, scaffolds).size();
			System.out.println("case="+mode+" exception="+(failure==null ? "none" : failure.getClass().getSimpleName())+
					" inputDescriptors="+inputDescriptors+" liveMetadataWriters="+liveWriters+
					" activeObjectWriters="+ReadWrite.countActiveThreads()+
					" summaryExists="+Files.exists(genome.resolve("summary.txt")));
		}finally{
			rescueWriters(before, info, scaffolds);
		}
		require(ReadWrite.countActiveThreads()==0, "Conversion left counted chromosome writes unfinished");
		if(collision){
			require(failure!=null && ("Tried to overwrite existing file "+first+", but overwrite=false.").equals(failure.getMessage()),
					"Expected the first chromosome overwrite collision, observed "+failure);
			require(Arrays.equals(SENTINEL, Files.readAllBytes(first)), "Collision changed the existing chromosome sentinel");
			require(!Files.exists(genome.resolve("summary.txt")), "Failed conversion published a completion summary");
			require(!Files.exists(genome.resolve("chr2.chrom")) && !Files.exists(genome.resolve("chr3.chrom")),
					"First-chromosome collision unexpectedly wrote later chromosomes");
			System.out.println("collisionVerified=true sentinelPreserved=true summaryAbsent=true");
		}else{
			if(failure!=null){throw failure;}
			checkSuccess(result, genome, input, writeChroms);
		}
		if(inputDescriptors!=0 || liveWriters!=0){
			System.out.println("FAIL_RESOURCE_CLEANUP inputDescriptors="+inputDescriptors+" liveMetadataWriters="+liveWriters);
			throw new AssertionError("Converter must close its input and join metadata writers before returning or throwing");
		}
		System.out.println("PASS "+mode);
	}

	/** Requires the runner's non-collecting-GC option rather than silently accepting a GC-masked pass. */
	private static void requireEpsilon(){
		boolean enabled=false;
		for(String arg : ManagementFactory.getRuntimeMXBean().getInputArguments()){
			if(arg.equals("-XX:+UseEpsilonGC")){enabled=true;}
			else if(arg.equals("-XX:-UseEpsilonGC")){enabled=false;}
		}
		require(enabled, "Descriptor observations require an Epsilon-GC JVM; source remains Java 8 compatible");
	}

	/** Independently checks that this process can observe an owned open input and its closure. */
	private static void checkDescriptorObservation(Path input) throws IOException{
		require(countInputDescriptors(input)==0, "Input unexpectedly open before calibration");
		try(InputStream in=Files.newInputStream(input)){
			require(in.read()=='>', "Literal FASTA must begin with a header");
			require(countInputDescriptors(input)==1, "Descriptor scanner failed its known-open control");
		}
		require(countInputDescriptors(input)==0, "Descriptor scanner failed its known-closed control");
	}

	/** Counts exact canonical input targets; only disappearing proc entries may be ignored. */
	private static int countInputDescriptors(Path input) throws IOException{
		int count=0;
		try(DirectoryStream<Path> entries=Files.newDirectoryStream(Paths.get("/proc/self/fd"))){
			for(Path entry : entries){
				try{
					if(Files.readSymbolicLink(entry).toString().equals(input.toString())){count++;}
				}catch(NoSuchFileException disappeared){
					// Another JVM descriptor can close during enumeration.
				}
			}
		}
		return count;
	}

	/** Selects only newly created live metadata writers at the exact two fixture destinations. */
	private static ArrayList<TextStreamWriter> ownedWriters(Set<Thread> before, String info, String scaffolds){
		final ArrayList<TextStreamWriter> owned=new ArrayList<TextStreamWriter>();
		for(Thread thread : Thread.getAllStackTraces().keySet()){
			if(!before.contains(thread) && thread instanceof TextStreamWriter && thread.isAlive()){
				final TextStreamWriter writer=(TextStreamWriter)thread;
				if(info.equals(writer.fname) || scaffolds.equals(writer.fname)){owned.add(writer);}
			}
		}
		return owned;
	}

	/** Ordinary tiny buffers fit the queue; poison then bounded joins rescue only this call's writers.
	 * The watchdog covers an unexpected blocking poison or failed join, not the test oracle. */
	private static void rescueWriters(Set<Thread> before, String info, String scaffolds) throws InterruptedException{
		final ArrayList<TextStreamWriter> owned=ownedWriters(before, info, scaffolds);
		for(TextStreamWriter writer : owned){writer.poison();}
		for(TextStreamWriter writer : owned){writer.join(5000);}
		for(TextStreamWriter writer : owned){
			require(!writer.isAlive(), "Fixture rescue could not join its metadata writer: "+writer.fname);
			require(!writer.errorState, "Fixture rescue encountered unexpected metadata output failure: "+writer.fname);
		}
		require(ownedWriters(before, info, scaffolds).isEmpty(), "Fixture rescue left an owned metadata writer alive");
		System.out.println("rescuedMetadataWriters="+owned.size());
	}

	/** Checks retained bases, independent serialized reads, and complete literal metadata. */
	private static void checkSuccess(ArrayList<ChromosomeArray> result, Path genome, Path input, boolean writeChroms) throws Exception{
		require(result!=null && result.size()==3, "Expected three retained chromosomes with merging disabled");
		for(int i=0; i<3; i++){
			checkChromosome(result.get(i), i);
			final Path file=genome.resolve("chr"+(i+1)+".chrom");
			if(writeChroms){
				try(ObjectInputStream in=new ObjectInputStream(Files.newInputStream(file))){
					final Object decoded=in.readObject();
					require(decoded instanceof ChromosomeArray, "Serialized output is not a chromosome");
					checkChromosome((ChromosomeArray)decoded, i);
				}
			}else{require(!Files.exists(file), "writechroms=false wrote a chromosome file");}
		}
		final String info=withoutDate(new String(Files.readAllBytes(genome.resolve("info.txt")), StandardCharsets.US_ASCII));
		require(info.equals("#Chromosome sizes\n#Version\t5\n#chrom\tscaffolds\tcontigs\tlength\tdefined\tundefined\tstartPad\tstopPad\n"+
				"1\t1\t1\t4\t4\t0\t0\t0\n2\t1\t1\t4\t4\t0\t0\t0\n3\t1\t1\t4\t4\t0\t0\t0\n"), "Chromosome metadata differs: "+info);
		final StringBuilder scaffoldText=new StringBuilder();
		try(InputStream in=new GZIPInputStream(Files.newInputStream(genome.resolve("scaffolds.txt.gz")))){
			for(int b; (b=in.read())!=-1;){scaffoldText.append((char)b);}
		}
		final String scaffolds=withoutDate(scaffoldText.toString());
		require(scaffolds.equals("#Scaffold names\n#Version\t5\n#chrom\tid\tstart\tlength\tname\n"+
				"1\t1\t0\t4\tone\n2\t2\t0\t4\ttwo\n3\t3\t0\t4\tthree\n"), "Scaffold metadata differs: "+scaffolds);
		final String summary=withoutDate(new String(Files.readAllBytes(genome.resolve("summary.txt")), StandardCharsets.US_ASCII));
		require(summary.equals("#Summary\n#Version\t5\nchroms\t3\nbases\t12\ndefined\t12\nundefined\t0\ncontigs\t3\nscaffolds\t3\n"+
				"interpad\t0\nname\ttiny\nsource\t"+input+"\nbytes\t"+Files.size(input)+
				"\nlast modified\t"+input.toFile().lastModified()+"\n"), "Summary metadata differs: "+summary);
		System.out.println("validatedChromosomes=3 validatedBases=12 serializedChromosomes="+(writeChroms ? 3 : 0)+" metadataExact=true");
	}

	/** Checks identifiers, strand, inclusive bounds, capacity and all four independently specified bases. */
	private static void checkChromosome(ChromosomeArray ca, int index){
		require(ca.chromosome==index+1 && ca.strand==Shared.PLUS && ca.minIndex==0 && ca.maxIndex==3 &&
				Arrays.equals(ca.array, BASES[index].getBytes(StandardCharsets.US_ASCII)), "Chromosome differs at index "+index);
	}

	/** Removes exactly one nondeterministic generated-date line; every other byte remains in the oracle. */
	private static String withoutDate(String text){
		final StringBuilder kept=new StringBuilder();
		int dates=0;
		for(String line : text.split("\n", -1)){
			if(line.startsWith("#Generated on\t")){dates++;}
			else{kept.append(line).append('\n');}
		}
		require(dates==1 && text.endsWith("\n"), "Expected one generated-date line and terminal newline");
		kept.setLength(kept.length()-1);//Drop only the split's final empty-field newline.
		return kept.toString();
	}

	/** Bounds any unexpected hang in this fresh test process; never reports expected baseline failure. */
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
					System.err.println("WATCHDOG_FAILURE exit=97: conversion test did not terminate normally");
					Runtime.getRuntime().halt(97);
				}
			}
		}, "FastaConversionCleanupTest-watchdog");
		watchdog.setDaemon(true);
		watchdog.start();
	}

	/** Keeps all oracles enabled under both assertions settings. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** The third header lets the second scaffold be prepared without consuming EOF. */
	private static final String FASTA=">one\nACGT\n>two\nTGCA\n>three\nGATT\n";
	/** Literal successful chromosome contents, independent of the converter's metadata. */
	private static final String[] BASES={"ACGT", "TGCA", "GATT"};
	/** Existing chromosome bytes must survive an overwrite=false rejection unchanged. */
	private static final byte[] SENTINEL="keep-existing-chromosome\n".getBytes(StandardCharsets.US_ASCII);
}
