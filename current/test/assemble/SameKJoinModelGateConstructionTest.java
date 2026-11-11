package test.assemble;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.io.StringWriter;
import java.lang.management.GarbageCollectorMXBean;
import java.lang.management.ManagementFactory;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.List;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;

/** Observes real model-gate constructor ownership in one fresh Linux JVM per case.
 * Uses tiny literal inputs, BF1 and exact-path descriptor access flags. A collection
 * during failure observation invalidates setup because GC can hide an unclosed input.
 * All failed constructions must leave zero descriptors, on every source revision.
 * Output failure is an ordinary blocked-parent open, not injected startup failure.
 * Generator cases use saved feature/model artifacts and one real directory collision,
 * gated by an externally created FIFO. Unreturned BF1 readers remain observable until
 * this isolated JVM exits; the fixture never closes a descriptor behind its owner's back.
 * Reflection preserves production visibility. This is not an assembly or speed test.
 * @author Brian Bushnell, Sucrose
 */
public final class SameKJoinModelGateConstructionTest{

	/** Accepts constructor cases or generator-healthy|generator-output-collision.
	 * Generator cases additionally require features= and expected=; collision needs fifo=.
	 * outdir= must be new, and the FIFO must already exist outside it. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, SameKJoinModelGateConstructionTest.class, false);
		try{
			String name=null;
			Path root=null;
			Path features=null, expected=null, fifo=null;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(pair.length==2 && !pair[1].isEmpty(), "Expected nonempty flag=value: "+arg);
				if(pair[0].equalsIgnoreCase("case")){name=pair[1];}
				else if(pair[0].equalsIgnoreCase("outdir")){root=Paths.get(pair[1]).toAbsolutePath().normalize();}
				else if(pair[0].equalsIgnoreCase("features")){features=Paths.get(pair[1]).toAbsolutePath().normalize();}
				else if(pair[0].equalsIgnoreCase("expected")){expected=Paths.get(pair[1]).toAbsolutePath().normalize();}
				else if(pair[0].equalsIgnoreCase("fifo")){fifo=Paths.get(pair[1]).toAbsolutePath().normalize();}
				else{throw new IllegalArgumentException("Unknown fixture option: "+arg);}
			}
			final boolean generator="generator-healthy".equals(name) || "generator-output-collision".equals(name);
			check(root!=null && (generator || "healthy".equals(name) || "metadata-refusal".equals(name) ||
					"output-open-failure".equals(name)), "Required a supported case and outdir=NEW.");
			if(generator){
				check(features!=null && Files.isRegularFile(features) && expected!=null && Files.isDirectory(expected),
						"Generator controls require saved regular features and the seven expected outputs.");
				check("generator-output-collision".equals(name)==(fifo!=null), "Only the collision case requires fifo=.");
			}else{check(features==null && expected==null && fifo==null, "Generator options do not belong to constructor cases.");}
			Files.createDirectory(root);
			final boolean[] original=backendFlags();
			final int threads=Shared.threads();
			final boolean passed;
			try{
				ByteFile.FORCE_MODE_BF1=true;
				ByteFile.FORCE_MODE_BF2=ByteFile.FORCE_MODE_BF3=ByteFile.FORCE_MODE_BF4=false;
				Shared.setThreads(1);
				final SameKJoinModelGateConstructionTest fixture=new SameKJoinModelGateConstructionTest();
				passed=generator ? fixture.runGenerator(name, root, features, expected, fifo) : fixture.run(name, root);
			}finally{
				ByteFile.FORCE_MODE_BF1=original[0]; ByteFile.FORCE_MODE_BF2=original[1];
				ByteFile.FORCE_MODE_BF3=original[2]; ByteFile.FORCE_MODE_BF4=original[3];
				Shared.setThreads(threads);
			}
			check(Arrays.equals(original, backendFlags()) && Shared.threads()==threads,
					"Fixture must restore every forced-reader flag and the thread setting.");
			final StringBuilder settings=new StringBuilder("setting\tbefore\tduring\tafter\n");
			for(int i=0; i<4; i++){
				settings.append("BF").append(i+1).append('\t').append(original[i]).append('\t')
						.append(i==0).append('\t').append(backendFlags()[i]).append('\n');
			}
			settings.append("threads\t").append(threads).append("\t1\t").append(Shared.threads()).append('\n');
			write(root.resolve("settings.tsv"), settings.toString());
			if(!passed){throw new AssertionError(generator ? "Model-gate fixture generator left its input open." :
					"Model gate left resources open after failed construction.");}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves private entry points and collector counters before the observed interval. */
	private SameKJoinModelGateConstructionTest() throws Exception{
		check(Files.isDirectory(FD_ROOT), "Exact descriptor observations require Linux /proc/self/fd.");
		final Class<?> gate=Class.forName("assemble.SameKJoinModelGate");
		constructor=gate.getDeclaredConstructor(String.class, String.class, String.class, int.class, float.class);
		constructor.setAccessible(true);
		finish=gate.getDeclaredMethod("finish"); finish.setAccessible(true);
		close=gate.getDeclaredMethod("close"); close.setAccessible(true);
		input=gate.getDeclaredField("input"); input.setAccessible(true);
		output=gate.getDeclaredField("output"); output.setAccessible(true);
		collectors=ManagementFactory.getGarbageCollectorMXBeans();
		check(!collectors.isEmpty(), "GC counters are needed to reject collection-masked leaks.");
		gcCounts();
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs unchanged generator code, preserving the accepted seven-file control.
	 * A FIFO controller owns only its coordination handle; it closes that before FD
	 * observation. No target reader handle returns on failure, so JVM exit is its rescue. */
	private boolean runGenerator(final String name, final Path root, final Path features,
			final Path expected, final Path fifo) throws Exception{
		check(SameKJoinModelGateConstructionTest.class.desiredAssertionStatus() && ByteStreamWriter.class.desiredAssertionStatus(),
				"The ordinary existing-output refusal must run with assertions enabled.");
		final boolean collision=fifo!=null;
		if(collision){
			check((((Number)Files.getAttribute(fifo, "unix:mode")).intValue()&0170000)==0010000,
					"fifo= must name a real Linux FIFO, not a regular file.");
			check(countFds(fifo, false)==0 && countFds(fifo, true)==0, "FIFO is already owned by this JVM.");
			// Linux O_RDWR opens a FIFO without waiting for another endpoint; load this path before observation.
			try(RandomAccessFile warm=new RandomAccessFile(fifo.toFile(), "rw")){check(warm.getFD().valid(), "FIFO warmup open failed.");}
		}
		final byte[] original=Files.readAllBytes(features);
		final Path warmup=root.resolve("generator-warmup"), generated=root.resolve("generated");
		SameKJoinModelGateTest.main(new String[]{"generate", features.toString(), warmup.toString()});
		checkGenerated(warmup, expected);
		final FifoBlocker blocker=collision ? new FifoBlocker(generated, fifo) : null;
		final Path observed=collision ? fifo : features;
		check(countFds(observed, false)==0 && countFds(observed, true)==0, "Observed input was not closed after warmup.");
		final long[] before=gcCounts();
		Throwable failure=null;
		try{
			if(blocker!=null){blocker.start();}
			SameKJoinModelGateTest.main(new String[]{"generate", observed.toString(), generated.toString()});
		}catch(RuntimeException|Error e){failure=e;}
		finally{if(blocker!=null){blocker.finish();}}
		final int readers=countFds(observed, false), writers=countFds(observed, true);
		final long[] after=gcCounts();
		final StringBuilder gc=new StringBuilder("collector\tbefore\tafter\n");
		for(int i=0; i<before.length; i++){
			gc.append(collectors.get(i).getName()).append('\t').append(before[i]).append('\t').append(after[i]).append('\n');
		}
		write(root.resolve("gc.tsv"), gc.toString());
		write(root.resolve("observation.tsv"), "reader_fds_after\twriter_fds_after\n"+readers+"\t"+writers+"\n");
		if(failure!=null){
			final StringWriter stack=new StringWriter();
			final PrintWriter writer=new PrintWriter(stack);
			failure.printStackTrace(writer); writer.flush();
			write(root.resolve("construction-error.txt"), stack.toString());
		}
		check(Arrays.equals(before, after), "GC during generator observation could hide its unclosed reader.");
		check(Arrays.equals(original, Files.readAllBytes(features)), "Accepted feature input was changed.");
		if(collision){
			final Path blocked=generated.resolve("badid.tsv");
			check(blocker.ready && blocker.closed && !blocker.isAlive(), "FIFO coordination did not complete and close.");
			check(failure!=null && failure.getClass()==AssertionError.class && failure.getCause()==null &&
					failure.getSuppressed().length==0 && failure.getMessage()!=null &&
					failure.getMessage().contains(blocked.toString()) && failure.getMessage().contains("overwrite=false"),
					"Expected the real existing-output AssertionError, observed: "+failure);
			boolean constructorFrame=false, generatorFrame=false;
			for(StackTraceElement frame : failure.getStackTrace()){
				constructorFrame|=frame.getClassName().equals("fileIO.ByteStreamWriter") && frame.getMethodName().equals("<init>");
				generatorFrame|=frame.getClassName().equals("test.assemble.SameKJoinModelGateTest") && frame.getMethodName().equals("generate");
			}
			check(constructorFrame && generatorFrame, "Refusal did not come from the intended generator writer acquisition.");
			check(Files.isDirectory(blocked) && countFds(blocked, true)==0, "Blocker changed or acquired an output descriptor.");
			try(DirectoryStream<Path> entries=Files.newDirectoryStream(blocked)){
				check(!entries.iterator().hasNext(), "Collision directory must stay empty.");
			}
			checkGeneratorFile(generated, expected, "constant.bbnet");
			checkGeneratorFile(generated, expected, "badscore.bbnet");
			check(entryCount(generated)==3, "Collision must precede all variant writes.");
			write(root.resolve("fifo.tsv"), "controller_joined\thelper_closed\tblocker_directory\tunreturned_reader_policy\ntrue\ttrue\ttrue\tisolated-jvm-exit\n");
		}else{
			check(failure==null, "Healthy generator failed: "+failure);
			checkGenerated(generated, expected);
		}
		check(writers==0, "FIFO coordination retained a writable descriptor.");
		final boolean passed=readers==0;
		final String report="case\texpected_reader_fds\tobserved_reader_fds\texpected_writer_fds\tobserved_writer_fds\tpassed\n"+
				name+"\t0\t"+readers+"\t0\t"+writers+"\t"+passed+"\n";
		write(root.resolve("checks.tsv"), report); System.out.print(report);
		System.out.println("MODEL_GATE_GENERATOR case="+name+" assertions=true setup=true passed="+passed);
		return passed;
	}

	/** Requires exactly the seven accepted outputs and no open writer for any of them. */
	private static void checkGenerated(final Path actual, final Path expected) throws Exception{
		check(entryCount(actual)==GENERATED.length && entryCount(expected)==GENERATED.length,
				"Generator output inventory must contain exactly seven files.");
		for(String file : GENERATED){checkGeneratorFile(actual, expected, file);}
	}

	/** Compares complete saved bytes, not regenerated expectations or counts. */
	private static void checkGeneratorFile(final Path actual, final Path expected, final String file) throws Exception{
		check(Files.isRegularFile(actual.resolve(file)) && Files.isRegularFile(expected.resolve(file)) &&
				Arrays.equals(Files.readAllBytes(actual.resolve(file)), Files.readAllBytes(expected.resolve(file))),
				"Generator output differs from accepted bytes: "+file);
		check(countFds(actual.resolve(file), true)==0, "Generator writer remains open: "+file);
	}

	/** Counts directory entries without accepting an extra report or unfinished output. */
	private static int entryCount(final Path directory) throws Exception{
		int count=0;
		try(DirectoryStream<Path> entries=Files.newDirectoryStream(directory)){for(Path ignored : entries){count++;}}
		return count;
	}

	/** Coordinates a real filesystem collision, without replacing a reader or throwing into the target. */
	private static final class FifoBlocker extends Thread{
		/** Owns only the FIFO writer handshake, never the target's reader or output files. */
		FifoBlocker(final Path generated_, final Path fifo_){
			super("model-gate-fixture-fifo");
			assert(generated_!=null && fifo_!=null) : "FIFO coordination requires the target's new directory and a prepared pipe.";
			generated=generated_; fifo=fifo_; setDaemon(true);
		}
		/** Installs the ordinary collision before allowing the target's FIFO input open to return. */
		@Override
		public void run(){
			try{
				while(!finished && !Files.isRegularFile(generated.resolve("badscore.bbnet"))){Thread.sleep(10);}
				if(finished){return;}
				Files.createDirectory(generated.resolve("badid.tsv"));
				handle=new RandomAccessFile(fifo.toFile(), "rw");
				ready=true;
			}catch(Exception|Error failure){KillSwitch.exceptionKill(failure);}// A failed controller must not leave main blocked in FIFO open.
		}
		/** Joins the controller and closes only its own FIFO handle before observing the target reader. */
		void finish() throws Exception{
			finished=true;
			boolean interrupted=false;
			while(isAlive()){
				try{join();}catch(InterruptedException e){interrupted=true;}
			}
			try{if(handle!=null){handle.close();} closed=true;}
			finally{if(interrupted){Thread.currentThread().interrupt();}}
		}
		/** Target output directory and externally prepared FIFO; neither is deleted here. */
		final Path generated, fifo;
		/** Main signals early failure so the controller cannot wait forever for a model file. */
		volatile boolean finished;
		/** Own read/write coordination handle; finish closes it before descriptor observation. */
		RandomAccessFile handle;
		/** Published by thread completion, then recorded after finish has joined the controller. */
		boolean ready, closed;
	}

	/** Warms the same real constructor/scanner, then records one finite case. */
	private boolean run(final String name, final Path root) throws Exception{
		final Path warmup=root.resolve("warmup");
		Files.createDirectory(warmup);
		healthy(warmup);
		final int readers;
		if(name.equals("healthy")){healthy(root); readers=0;}
		else{readers=failedConstruction(root, name.equals("metadata-refusal"));}
		final boolean passed=readers==0;
		final String report="case\texpected_reader_fds\tobserved_reader_fds\texpected_writer_fds\tobserved_writer_fds\tpassed\n"+
				name+"\t0\t"+readers+"\t0\t0\t"+passed+"\n";
		write(root.resolve("checks.tsv"), report);
		System.out.print(report);
		System.out.println("MODEL_GATE_CONSTRUCTION case="+name+" assertions="+
				SameKJoinModelGateConstructionTest.class.desiredAssertionStatus()+" setup=true passed="+passed);
		return passed;
	}

	/** Confirms live ownership, then empty-population finish, terminal close and exact header bytes. */
	private void healthy(final Path root) throws Exception{
		final Path features=root.resolve("features.tsv"), model=root.resolve("model.bbnet"), scores=root.resolve("scores.tsv");
		write(features, features(false)); write(model, model());
		check(countFds(features, false)==0 && countFds(scores, true)==0, "Fresh control paths already have descriptors.");
		final Object gate=create(features, model, scores);
		ByteStreamWriter writer=null;
		try{
			writer=(ByteStreamWriter)output.get(gate);
			check(input.get(gate) instanceof ByteFile1, "The bounded constructor witness requires BF1 without reader workers.");
			check(countFds(features, false)==1 && countFds(scores, true)==1 && writer.isAlive(),
					"Healthy construction must expose one feature reader and one started score writer to the scanner.");
			invoke(finish, gate);
		}finally{check(!((Boolean)invoke(close, gate)), "Healthy gate close reported an I/O error.");}
		check(countFds(features, false)==0 && countFds(scores, true)==0 && !writer.isAlive(),
				"Healthy finish/close must release both paths and join its writer.");
		checkBytes(features, features(false)); checkBytes(model, model()); checkBytes(scores, SCORE_HEADER);
		write(root.resolve("observation.tsv"), "reader_fds_after\twriter_fds_after\n0\t0\n");
	}

	/** Observes constructor refusal before any GC; cannot rescue an object that never returned. */
	private int failedConstruction(final Path root, final boolean metadata) throws Exception{
		final Path features=root.resolve("features.tsv"), model=root.resolve("model.bbnet"), blocker=root.resolve("blocker");
		final Path scores=metadata ? root.resolve("scores.tsv") : blocker.resolve("scores.tsv");
		write(features, features(metadata)); write(model, model());
		if(!metadata){write(blocker, BLOCKER);}
		check(countFds(features, false)==0 && countFds(scores, true)==0, "Failure observation needs unopened paths.");
		final long[] before=gcCounts();
		RuntimeException failure=null;
		Object unexpected=null;
		try{unexpected=create(features, model, scores);}
		catch(RuntimeException e){failure=e;}
		finally{if(unexpected!=null){check(!((Boolean)invoke(close, unexpected)), "Unexpected gate could not be closed.");}}
		final int readers=countFds(features, false), writers=countFds(scores, true);
		final long[] after=gcCounts();
		final StringBuilder gc=new StringBuilder("collector\tbefore\tafter\n");
		for(int i=0; i<before.length; i++){
			gc.append(collectors.get(i).getName()).append('\t').append(before[i]).append('\t').append(after[i]).append('\n');
		}
		write(root.resolve("gc.tsv"), gc.toString());
		write(root.resolve("observation.tsv"), "reader_fds_after\twriter_fds_after\n"+readers+"\t"+writers+"\n");
		if(failure!=null){
			final StringWriter trace=new StringWriter();
			final PrintWriter writer=new PrintWriter(trace);
			failure.printStackTrace(writer); writer.flush();
			write(root.resolve("construction-error.txt"), trace.toString());
		}
		check(Arrays.equals(before, after), "GC during failure observation could conceal the missing input close.");
		check(unexpected==null && failure!=null, "Expected construction to refuse the prepared failure input.");
		if(metadata){
			check(failure.getClass()==IllegalArgumentException.class && failure.getCause()==null &&
					"Wrong same-K feature schema.".equals(failure.getMessage()), "Expected the real metadata rejection: "+failure);
		}else{
			check(failure.getClass()==RuntimeException.class && failure.getCause() instanceof FileNotFoundException &&
					failure.getCause().getMessage()!=null && failure.getCause().getMessage().contains(scores.toString()),
					"Expected blocked-parent output creation with its original FileNotFoundException cause: "+failure);
			checkBytes(blocker, BLOCKER);
		}
		check(failure.getSuppressed().length==0, "Ordinary plain-file cleanup must not introduce suppressed errors.");
		check(writers==0 && !Files.exists(scores), "Rejected construction must not open or create score output.");
		checkBytes(features, features(metadata)); checkBytes(model, model());
		return readers;
	}

	/** Calls the real constructor and exposes its exact unchecked failure, without wrapping it. */
	private Object create(final Path features, final Path model, final Path scores) throws Exception{
		try{return constructor.newInstance(features.toString(), model.toString(), scores.toString(), 32, 0.5f);}
		catch(InvocationTargetException e){throwOriginal(e); return null;}
	}

	/** Invokes real finish/close; checked reflection setup errors remain fixture failures. */
	private static Object invoke(final Method method, final Object target) throws Exception{
		try{return method.invoke(target);}
		catch(InvocationTargetException e){throwOriginal(e); return null;}
	}

	/** Preserves the exact original RuntimeException/Error and its cause chain. */
	private static void throwOriginal(final InvocationTargetException wrapped) throws Exception{
		final Throwable original=wrapped.getCause();
		if(original instanceof RuntimeException){throw (RuntimeException)original;}
		if(original instanceof Error){throw (Error)original;}
		throw wrapped;
	}

	/** Snapshots every collector; unknown counters cannot exclude hidden collection. */
	private long[] gcCounts(){
		final long[] counts=new long[collectors.size()];
		for(int i=0; i<counts.length; i++){
			counts[i]=collectors.get(i).getCollectionCount();
			check(counts[i]>=0, "Unavailable GC count: "+collectors.get(i).getName());
		}
		return counts;
	}

	/** Counts only matching canonical paths with the requested access; absent output paths are valid. */
	private static int countFds(final Path file, final boolean writable) throws Exception{
		final Path real=file.toFile().getCanonicalFile().toPath();
		int count=0;
		try(DirectoryStream<Path> descriptors=Files.newDirectoryStream(FD_ROOT)){
			for(Path fd : descriptors){
				final Path target;
				try{target=Files.readSymbolicLink(fd);}
				catch(NoSuchFileException disappeared){continue;}
				if(target.equals(real)){
					int flags=-1;
					for(String line : Files.readAllLines(Paths.get("/proc/self/fdinfo", fd.getFileName().toString()), StandardCharsets.US_ASCII)){
						if(line.startsWith("flags:")){
							check(flags==-1, "Duplicate descriptor access flags: "+fd);
							flags=Integer.parseInt(line.substring(6).trim(), 8);
						}
					}
					check(flags>=0 && (flags&3)!=3, "Missing or invalid descriptor access mode: "+fd);
					final int access=flags&3;
					if(access==2 || access==(writable ? 1 : 0)){count++;}
				}
			}
		}
		return count;
	}

	/** Literal empty-population schema; malformed mode changes only its schema name. */
	private static String features(final boolean malformed){
		final StringBuilder text=new StringBuilder(malformed ? "#schema=wrong" : "#schema=samek27-v1");
		text.append(" unique=0 observations=0.0 squares=0.0\nkind\tevent");
		for(int i=0; i<27; i++){text.append("\tf").append(i);}
		return text.append("\n#controls population=0 sampled=0 seed=1 requested=0\n").toString();
	}

	/** One linear output with constant0.8 bias and27zero weights; no training required. */
	private static String model(){
		final StringBuilder text=new StringBuilder("##bbnet\n#version 1\n#concise\n#dense\n#density 1\n#blocksize 1\n#seed 1\n#layers 2\n#dims 27 1\n#edges 27\n#coding decimal\n##layer 1\nC28 LINEAR 0.8");
		for(int i=0; i<27; i++){text.append(" 0");}
		return text.append('\n').toString();
	}

	/** Captures all process-wide reader-force settings in a fixed order. */
	private static boolean[] backendFlags(){
		return new boolean[]{ByteFile.FORCE_MODE_BF1, ByteFile.FORCE_MODE_BF2, ByteFile.FORCE_MODE_BF3, ByteFile.FORCE_MODE_BF4};
	}

	/** Requires complete tiny fixture byte equality, never just an expected line count. */
	private static void checkBytes(final Path path, final String expected) throws Exception{
		check(Arrays.equals(Files.readAllBytes(path), expected.getBytes(StandardCharsets.US_ASCII)), "Fixture bytes changed: "+path);
	}

	/** Creates a small fresh record without replacing earlier evidence. */
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Keeps setup and ownership oracles enabled regardless of JVM assertion switches. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Exact constructor and lifecycle methods; no production visibility changes. */
	private final Constructor<?> constructor;
	private final Method finish, close;
	/** Existing owned-stream fields used only to verify the healthy control. */
	private final Field input, output;
	/** Stable collector roster for the before/after observation. */
	private final List<GarbageCollectorMXBean> collectors;
	/** Linux descriptor namespace, scoped to this fixture JVM. */
	private static final Path FD_ROOT=Paths.get("/proc/self/fd");
	/** Literal score trace for an empty candidate population. */
	private static final String SCORE_HEADER="event\tk\tscore_bits\tscore\taccepted\n";
	/** Immutable regular file that makes creating a child output impossible. */
	private static final String BLOCKER="unchanged model-gate output blocker\n";
	/** Exact artifact set produced by the previously accepted generator capture. */
	private static final String[] GENERATED={"constant.bbnet", "badscore.bbnet", "badid.tsv", "badwidth.tsv", "nonfinite.tsv", "badk.tsv", "missingtrailer.tsv"};
}
