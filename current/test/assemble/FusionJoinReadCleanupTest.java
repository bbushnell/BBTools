package test.assemble;

import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.lang.management.GarbageCollectorMXBean;
import java.lang.management.ManagementFactory;
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipException;

import assemble.FusionJoinDiagnostic;
import assemble.SameKJoinAudit;
import assemble.SameKTrainingSample;
import assemble.SameKModelEvaluate;
import assemble.LiveBranchModelEvaluate;
import assemble.LiveBranchTrainingCorpus;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.Parse;
import parse.PreParser;
import parse.LineParser1;
import shared.Shared;

/** Finite real-reader cleanup checks for trace parsing and paired input acquisition.
 * One case runs per fresh Linux JVM with BF1 forced temporarily and restored in finally.
 * Exact readable descriptors, unchanged input bytes and unchanged GC counters separate
 * parser cleanup from unrelated I/O or collection. Paired-input modes call the real
 * audit/sampler mains with one thread and Java gzip, including normal output workers.
 * Planner cases separately exercise public plan dispatch and real output acquisition.
 * plan-close-* cases inject an unchecked diagnostic failure after real input closure,
 * distinguish EOF from processing-finally closure, and rescue owned output workers.
 * Evaluator cases run a pinned tiny analytic model; no assembler or training runs.
 * Corrupt gzip is real decoder failure, not injection.
 * training-pair-* cases use static readable directories to test real input-open failures.
 * The independent post-call descriptor expectation is zero on every source revision.
 * @author Brian Bushnell, Sucrose
 */
public final class FusionJoinReadCleanupTest{

	/** Accepts original parser cases or audit-valid, audit-gzip-failure, sample-valid,
	 * sample-gzip-failure, or isolated plan-valid/plan-output-failure, with outdir=NEW.
	 * plan-close-valid, plan-close-malformed, plan-close-eof-fault and
	 * plan-close-finally-fault also use outdir=NEW.
	 * training-pair-valid, training-pair-first-directory, training-pair-second-directory
	 * and training-pair-malformed use one literal validation genome and outdir=NEW.
	 * samek-eval-valid/samek-eval-gzip-failure and live-eval-valid/live-eval-gzip-failure
	 * additionally require inputs=PINNED_DIRECTORY.
	 * samek-eval-output-failure and live-eval-output-failure use their respective pinned
	 * inputs and an isolated output-open observer. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, FusionJoinReadCleanupTest.class, false);
		try{
			if(evaluatorOutputCase(pp.args)){
				runEvaluatorOutputCase(pp.args);
				return;
			}
			if(evaluatorCase(pp.args)){
				runEvaluatorCase(pp.args);
				return;
			}
			String name=null;
			Path root=null;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(pair.length==2 && !pair[1].isEmpty(), "Expected nonempty flag=value: "+arg);
				if(pair[0].equalsIgnoreCase("case")){name=pair[1];}
				else if(pair[0].equalsIgnoreCase("outdir")){root=Paths.get(pair[1]).toAbsolutePath().normalize();}
				else{throw new IllegalArgumentException("Unknown fixture option: "+arg);}
			}
			if("plan-valid".equals(name) || "plan-output-failure".equals(name)){
				runPlannerCase(name, root);
				return;
			}
			if(plannerCloseCase(name)){
				runPlannerCloseCase(name, root);
				return;
			}
			if(trainingPairCase(name)){
				runTrainingPairCase(name, root);
				return;
			}
			check(root!=null && ("valid".equals(name) || "malformed".equals(name) ||
					"empty-allowed".equals(name) || "empty-rejected".equals(name) || pairedCase(name)),
					"Required a supported parser or paired-reader case and outdir=NEW.");
			Files.createDirectory(root);
			final boolean[] original=backendFlags();
			final int originalThreads=Shared.threads();
			final boolean originalRaw=ReadWrite.RAWMODE, paired=pairedCase(name);
			final int remaining;
			try{
				ByteFile.FORCE_MODE_BF1=true;
				ByteFile.FORCE_MODE_BF2=ByteFile.FORCE_MODE_BF3=ByteFile.FORCE_MODE_BF4=false;
				if(paired){Shared.setThreads(1); ReadWrite.RAWMODE=false;}
				final FusionJoinReadCleanupTest test=new FusionJoinReadCleanupTest();
				remaining=paired ? test.runPaired(name, root) : test.run(name, root);
			}finally{
				ByteFile.FORCE_MODE_BF1=original[0]; ByteFile.FORCE_MODE_BF2=original[1];
				ByteFile.FORCE_MODE_BF3=original[2]; ByteFile.FORCE_MODE_BF4=original[3];
				if(paired){Shared.setThreads(originalThreads); ReadWrite.RAWMODE=originalRaw;}
			}
			final boolean[] restored=backendFlags();
			check(Arrays.equals(original, restored), "Forced reader flags were not restored.");
			final StringBuilder flags=new StringBuilder("flag\tbefore\tafter\n");
			for(int i=0; i<original.length; i++){
				flags.append("BF").append(i+1).append('\t').append(original[i]).append('\t').append(restored[i]).append('\n');
			}
			write(root.resolve("backend-flags.tsv"), flags.toString());
			if(paired){
				check(Shared.threads()==originalThreads && ReadWrite.RAWMODE==originalRaw, "Paired-reader settings were not restored.");
				write(root.resolve("settings.tsv"), "setting\tbefore\tduring\tafter\nthreads\t"+originalThreads+
						"\t1\t"+Shared.threads()+"\nrawmode\t"+originalRaw+"\tfalse\t"+ReadWrite.RAWMODE+"\n");
			}
			final boolean passed=remaining==0;
			final String report="case\texpected_reader_fds\tobserved_reader_fds\tpassed\n"+
					name+"\t0\t"+remaining+"\t"+passed+"\n";
			write(root.resolve("checks.tsv"), report);
			System.out.print(report);
			System.out.println("FUSION_JOIN_READ case="+name+" assertions="+
					FusionJoinReadCleanupTest.class.desiredAssertionStatus()+" setup=true flags_restored=true passed="+passed);
			if(!passed){throw new AssertionError((paired ? "Paired-reader setup left an input open after rejection: " :
					"Fusion trace parser left an input open after rejection: ")+remaining);}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Selects the isolated output-open witness; old evaluator cases retain their own dispatch. */
	private static boolean evaluatorOutputCase(final String[] args){
		String name=null;
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length==2 && pair[0].equalsIgnoreCase("case")){name=pair[1];}
		}
		return "samek-eval-output-failure".equals(name) || "live-eval-output-failure".equals(name);
	}

	/** One healthy warmup and one real output-open failure, with desired zero readers on every revision. */
	private static void runEvaluatorOutputCase(final String[] args) throws Exception{
		Path root=null, inputs=null;
		String name=null;
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			check(pair.length==2 && !pair[1].isEmpty(), "Expected nonempty flag=value: "+arg);
			if(pair[0].equalsIgnoreCase("outdir")){root=Paths.get(pair[1]).toAbsolutePath().normalize();}
			else if(pair[0].equalsIgnoreCase("inputs")){inputs=Paths.get(pair[1]).toAbsolutePath().normalize();}
			else if(pair[0].equalsIgnoreCase("case")){name=pair[1];}
			else{throw new IllegalArgumentException("Unknown output-open fixture option: "+arg);}
		}
		final boolean live="live-eval-output-failure".equals(name);
		check(root!=null && inputs!=null && Files.isDirectory(inputs) && Files.isDirectory(FD_ROOT) &&
				FusionJoinReadCleanupTest.class.desiredAssertionStatus(), "Output-open witness requires Linux, -ea, inputs=PINNED_DIRECTORY and outdir=NEW.");
		Files.createDirectory(root);
		final boolean[] modes=backendFlags();
		final int threads=Shared.threads(), buffers=Shared.numBuffers();
		final boolean raw=ReadWrite.RAWMODE, fma=Shared.SIMD_FMA, verbose=ReadWrite.verbose;
		final PrintStream stderr=System.err;
		final int remaining;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=ByteFile.FORCE_MODE_BF3=ByteFile.FORCE_MODE_BF4=false;
			Shared.setThreads(1);ReadWrite.RAWMODE=false;Shared.SIMD_FMA=false;ReadWrite.verbose=false;
			final FusionJoinReadCleanupTest test=new FusionJoinReadCleanupTest();
			final Path warmup=root.resolve("eval-warmup");Files.createDirectory(warmup);
			check(test.observeEvaluator(warmup, inputs, false, false, live)==0, "Healthy output-open control leaked a reader.");
			remaining=test.observeEvaluatorOutput(root, inputs, live);
		}finally{
			System.setErr(stderr);ReadWrite.verbose=verbose;ReadWrite.RAWMODE=raw;Shared.SIMD_FMA=fma;
			ByteFile.FORCE_MODE_BF1=modes[0];ByteFile.FORCE_MODE_BF2=modes[1];
			ByteFile.FORCE_MODE_BF3=modes[2];ByteFile.FORCE_MODE_BF4=modes[3];
			Shared.setThreads(threads);Shared.setBuffers(buffers);
		}
		check(System.err==stderr && ReadWrite.verbose==verbose && ReadWrite.RAWMODE==raw && Shared.SIMD_FMA==fma &&
				Shared.threads()==threads && Shared.numBuffers()==buffers && Arrays.equals(modes, backendFlags()),
				"Output-open witness must restore every changed global.");
		final StringBuilder flags=new StringBuilder("flag\tbefore\tafter\n");
		for(int i=0; i<4; i++){flags.append("BF").append(i+1).append('\t').append(modes[i]).append('\t').append(backendFlags()[i]).append('\n');}
		write(root.resolve("backend-flags.tsv"), flags.toString());
		write(root.resolve("settings.tsv"), "setting\tbefore\tafter\nthreads\t"+threads+"\t"+Shared.threads()+
				"\nbuffers\t"+buffers+"\t"+Shared.numBuffers()+"\nrawmode\t"+raw+"\t"+ReadWrite.RAWMODE+
				"\nsimdfma\t"+fma+"\t"+Shared.SIMD_FMA+"\nverbose\t"+verbose+"\t"+ReadWrite.verbose+"\nstderr_restored\ttrue\ttrue\n");
		final boolean passed=remaining==0;
		final String report="case\texpected_reader_fds\tobserved_reader_fds\tpassed\n"+name+"\t0\t"+remaining+"\t"+passed+"\n";
		write(root.resolve("checks.tsv"), report);System.out.print(report);
		System.out.println((live ? "LIVE_EVAL_OUTPUT" : "SAMEK_EVAL_OUTPUT")+" assertions=true setup=true globals_restored=true passed="+passed);
		if(!passed){throw new AssertionError("Evaluator left inputs open after prediction output rejection: "+remaining);}
	}

	/** Blocks only scores.tsv at the synchronous open diagnostic, after observing both real input descriptors. */
	private int observeEvaluatorOutput(final Path root, final Path inputs, final boolean live) throws Exception{
		for(String file : new String[]{"rows.tsv", "validation.tsv", "model.bbnet"}){Files.copy(inputs.resolve(file), root.resolve(file));}
		final Path meta=root.resolve("rows.tsv"), data=root.resolve("validation.tsv.gz"), model=root.resolve("model.bbnet");
		final Path output=root.resolve("evaluation"), scores=output.resolve("scores.tsv");
		try(OutputStream stream=new GZIPOutputStream(Files.newOutputStream(data, StandardOpenOption.CREATE_NEW))){
			stream.write(Files.readAllBytes(root.resolve("validation.tsv")));
		}
		final byte[] compressed=Files.readAllBytes(data);
		check(countReaderFds(meta)==0 && countReaderFds(data)==0 && !Files.exists(output), "Output-open paths must start unopened.");
		final String caller=live ? "assemble.LiveBranchModelEvaluate" : "assemble.SameKModelEvaluate";
		final String success=live ? "LIVE_BRANCH_EVALUATE_PASS " : "SAME_K_EVALUATE_PASS ";
		final String[] args={"in="+data, "rows="+meta, "net="+model, live ? "partition=validation" : "fold=validation", "out="+output,
				"simdff=false", "simdfma=false", "t=1", "bf1=t"};
		write(root.resolve("call.txt"), caller+"\n"+String.join("\n", args)+"\n");
		final ByteArrayOutputStream diagnostic=new ByteArrayOutputStream();
		final EvaluatorOutputObserver observer=new EvaluatorOutputObserver(diagnostic, meta, data, scores, caller, success);
		final PrintStream stderr=System.err;
		final boolean verbose=ReadWrite.verbose;
		Throwable failure=null;
		final int metaFds, dataFds, outputFds, workers;
		final long[] before=gcCounts(), after;
		try{
			System.setErr(observer);ReadWrite.verbose=true;
			try{if(live){LiveBranchModelEvaluate.main(args);}else{SameKModelEvaluate.main(args);}}
			catch(RuntimeException|Error caught){failure=caught;}
			metaFds=countReaderFds(meta);dataFds=countReaderFds(data);
			outputFds=countPlannerWriterFds(scores);workers=evaluatorOutputWorkers(scores);after=gcCounts();
		}finally{System.setErr(stderr);ReadWrite.verbose=verbose;}
		write(root.resolve("evaluation.stderr"), diagnostic.toString("UTF-8"));
		write(root.resolve("open-hook.tsv"), "opens\tmetadata_fds\tfeature_fds\twriter_fds\tworkers\tblock_created\n"+
				observer.opens+"\t"+observer.metaFds+"\t"+observer.dataFds+"\t"+observer.writerFds+"\t"+observer.workers+"\t"+observer.blocked+"\n");
		write(root.resolve("output-observation.tsv"), "metadata_fds\tfeature_fds\twriter_fds\tworkers\tfailure_class\tcause_class\n"+
				metaFds+"\t"+dataFds+"\t"+outputFds+"\t"+workers+"\t"+(failure==null ? "none" : failure.getClass().getName())+
				"\t"+(failure==null || failure.getCause()==null ? "none" : failure.getCause().getClass().getName())+"\n");
		final StringBuilder gc=new StringBuilder("collector\tbefore\tafter\n");
		for(int i=0; i<before.length; i++){gc.append(collectors.get(i).getName()).append('\t').append(before[i]).append('\t').append(after[i]).append('\n');}
		write(root.resolve("gc.tsv"), gc.toString());
		if(failure!=null){
			final StringWriter trace=new StringWriter();final PrintWriter writer=new PrintWriter(trace);
			failure.printStackTrace(writer);writer.flush();write(root.resolve("output-error.txt"), trace.toString());
		}
		check(observer.setupFailure==null && observer.opens==1 && observer.blocked && observer.successes==0,
				"Output-open witness must create one blocker after real input acquisition, without a setup failure or success report.");
		check(Arrays.equals(before, after), "GC could conceal missing evaluator input cleanup.");
		check(failure!=null && failure.getClass()==RuntimeException.class && failure.getCause() instanceof FileNotFoundException &&
				failure.getCause().getMessage().equals(scores+" (Is a directory)") && failure.getSuppressed().length==0 &&
				failure.getCause().getSuppressed().length==0 && plannerFrame(failure, "fileIO.ReadWrite", "getRawOutputStream") &&
				plannerFrame(failure, "fileIO.ByteStreamWriter", "<init>") && plannerFrame(failure, caller, "open") &&
				plannerFrame(failure, caller, "run") && plannerFrame(failure, caller, "main"),
				"Expected the original real output-open failure after both inputs, not an earlier rejection: "+failure);
		check(outputFds==0 && workers==0 && Files.isDirectory(scores), "Failed output acquisition must not leave an output stream/worker.");
		try(DirectoryStream<Path> files=Files.newDirectoryStream(scores)){check(!files.iterator().hasNext(), "Blocker directory changed.");}
		try(DirectoryStream<Path> files=Files.newDirectoryStream(output)){
			int count=0;for(Path path : files){check(path.equals(scores), "Failed output acquisition published another output: "+path);count++;}
			check(count==1, "Only the deliberate scores.tsv blocker may exist after failure.");
		}
		for(String file : new String[]{"rows.tsv", "validation.tsv", "model.bbnet"}){
			check(Arrays.equals(Files.readAllBytes(inputs.resolve(file)), Files.readAllBytes(root.resolve(file))), "Evaluator changed input: "+file);
		}
		check(Arrays.equals(compressed, Files.readAllBytes(data)), "Evaluator changed compressed input bytes.");
		return metaFds+dataFds;
	}

	/** Counts only live output workers for the single path this witness can acquire. */
	private static int evaluatorOutputWorkers(final Path scores){
		int count=0;
		for(Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread instanceof ByteStreamWriter && thread.isAlive() && ((ByteStreamWriter)thread).fname.equals(scores.toString())){count++;}
		}
		return count;
	}

	/** Arranges an ordinary filesystem open failure; never substitutes a backend exception. */
	private static final class EvaluatorOutputObserver extends PrintStream{
		/** Borrows diagnostics and exact paths; never takes ownership of the original stderr. */
		EvaluatorOutputObserver(final ByteArrayOutputStream sink, final Path meta_, final Path data_, final Path scores_,
				final String caller_, final String success_) throws Exception{
			super(sink, true, "UTF-8");meta=meta_;data=data_;scores=scores_;caller=caller_;success=success_;
		}
		/** Confirms live input ownership before creating the one filesystem blocker. */
		@Override
		public void println(final String text){
			super.println(text);
			if(text.startsWith(success)){successes++;}
			if(!text.equals("getOutputStream("+scores+", false, true, false)")){return;}
			try{
				metaFds=countReaderFds(meta);dataFds=countReaderFds(data);
				writerFds=countPlannerWriterFds(scores);workers=evaluatorOutputWorkers(scores);
				check(++opens==1 && metaFds==1 && dataFds==1 && writerFds==0 && workers==0 && !Files.exists(scores),
						"Prediction output must follow two live input descriptors and precede writer acquisition.");
				final Throwable location=new Throwable();
				check(ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && !ByteFile.FORCE_MODE_BF3 && !ByteFile.FORCE_MODE_BF4 &&
						Shared.threads()==1 && !ReadWrite.RAWMODE && plannerFrame(location, "fileIO.ReadWrite", "getOutputStream") &&
						plannerFrame(location, caller, "open") && plannerFrame(location, caller, "run"),
						"Block only the real serial evaluator output acquisition with BF1 inputs.");
				Files.createDirectory(scores);blocked=true;
			}catch(Exception|Error failure){setupFailure=failure;throw new IllegalStateException("Evaluator output-open setup failed.", failure);}
		}
		/** Synchronous observations only; no worker thread performs fixture mutation. */
		private final Path meta, data, scores;
		/** Explicit selected evaluator identity and its success prefix, never inferred from the failure. */
		private final String caller, success;
		private int opens, successes, metaFds=-1, dataFds=-1, writerFds=-1, workers=-1;
		private boolean blocked;
		private Throwable setupFailure;
	}

	/** Four production converter cases; directory failures are malformed filesystem inputs, not injected exceptions. */
	private static boolean trainingPairCase(final String name){
		return "training-pair-valid".equals(name) || "training-pair-first-directory".equals(name) ||
				"training-pair-second-directory".equals(name) || "training-pair-malformed".equals(name);
	}

	/** Warms the healthy path, measures one call without GC, and restores every changed global. */
	private static void runTrainingPairCase(final String name, final Path root) throws Exception{
		check(root!=null && Files.isDirectory(FD_ROOT) && FusionJoinReadCleanupTest.class.desiredAssertionStatus(),
				"Training input panel requires Linux, -ea and outdir=NEW.");
		Files.createDirectory(root);
		final boolean[] modes=backendFlags();
		final int threads=Shared.threads(), buffers=Shared.numBuffers();
		final boolean raw=ReadWrite.RAWMODE;
		final int remaining;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=ByteFile.FORCE_MODE_BF3=ByteFile.FORCE_MODE_BF4=false;
			Shared.setThreads(1);ReadWrite.RAWMODE=false;
			final FusionJoinReadCleanupTest test=new FusionJoinReadCleanupTest();
			final Path warmup=root.resolve("warmup");Files.createDirectory(warmup);
			check(test.observeTrainingPair(warmup, "training-pair-valid", false)==0, "Healthy training warmup leaked an input.");
			remaining=test.observeTrainingPair(root, name, true);
		}finally{
			ByteFile.FORCE_MODE_BF1=modes[0];ByteFile.FORCE_MODE_BF2=modes[1];
			ByteFile.FORCE_MODE_BF3=modes[2];ByteFile.FORCE_MODE_BF4=modes[3];
			Shared.setThreads(threads);Shared.setBuffers(buffers);ReadWrite.RAWMODE=raw;
		}
		check(Arrays.equals(modes, backendFlags()) && Shared.threads()==threads && Shared.numBuffers()==buffers && ReadWrite.RAWMODE==raw,
				"Training fixture must restore reader flags, threads, buffers and raw mode.");
		write(root.resolve("settings.tsv"), "setting\tbefore\tafter\nthreads\t"+threads+"\t"+Shared.threads()+
				"\nbuffers\t"+buffers+"\t"+Shared.numBuffers()+"\nrawmode\t"+raw+"\t"+ReadWrite.RAWMODE+"\n");
		final StringBuilder flags=new StringBuilder("flag\tbefore\tafter\n");
		for(int i=0; i<4; i++){flags.append("BF").append(i+1).append('\t').append(modes[i]).append('\t').append(backendFlags()[i]).append('\n');}
		write(root.resolve("backend-flags.tsv"), flags.toString());
		final boolean passed=remaining==0;
		final String report="case\texpected_reader_fds\tobserved_reader_fds\tpassed\n"+name+"\t0\t"+remaining+"\t"+passed+"\n";
		write(root.resolve("checks.tsv"), report);System.out.print(report);
		System.out.println("LIVE_TRAINING_READ case="+name+" assertions=true setup=true globals_restored=true passed="+passed);
		if(!passed){throw new AssertionError("Training converter left vectors open after metadata acquisition failed.");}
	}

	/** Runs the public converter with static inputs and observes resources before writer rescue or assertions. */
	private int observeTrainingPair(final Path root, final String name, final boolean measured) throws Exception{
		final boolean firstDirectory=name.equals("training-pair-first-directory"), secondDirectory=name.equals("training-pair-second-directory");
		final boolean malformed=name.equals("training-pair-malformed"), healthy=name.equals("training-pair-valid");
		final Path labels=root.resolve("labels"), run=root.resolve("run"), output=root.resolve("result");
		Files.createDirectory(labels);Files.createDirectory(run);
		final Path data=labels.resolve("vectors.tsv"), meta=labels.resolve("rows.tsv"), manifest=root.resolve("selected.tsv");
		final String inputRows=TRAIN_PAIR_ROWS+(malformed ? "1\t"+TRAIN_PAIR_WORD+"\tgrey\n" : "1\t"+TRAIN_PAIR_WORD+"\tgrey\t1\n");
		if(firstDirectory){Files.createDirectory(data);}else{write(data, TRAIN_PAIR_VECTOR_HEADER+TRAIN_PAIR_VECTOR);}
		if(secondDirectory){Files.createDirectory(meta);}else{write(meta, inputRows);}
		write(run.resolve("exit_code.txt"), "0\n");write(labels.resolve("census.tsv"), TRAIN_PAIR_CENSUS);
		write(labels.resolve("complete.txt"), TRAIN_PAIR_INPUT_COMPLETE);
		final String selected=TRAIN_PAIR_MANIFEST+"g0001\tvalidation\t24\t30\t123\tnew\t"+run+"\t"+labels+"\tunused.fq\tunused.fa\n";
		write(manifest, selected);
		check(Files.isReadable(data) && Files.isReadable(meta), "Both static input paths must pass the converter readability preflight.");
		check(countReaderFds(data)==0 && countReaderFds(meta)==0, "Fresh training inputs already have readers.");
		if(!measured){
			final ByteFile probe=ByteFile.makeByteFile(data.toString(), true);
			try{check(probe.getClass()==ByteFile1.class && countReaderFds(data)==1, "Training observation requires synchronous BF1 with a visible descriptor.");}
			finally{check(!probe.close(), "Training BF1 probe failed to close.");}
			check(countReaderFds(data)==0, "Training BF1 probe descriptor remained open.");
			write(root.resolve("preflight.tsv"), "backend\treader_fds_live\treader_fds_closed\nByteFile1\t1\t0\n");
		}
		final String[] args={"in="+manifest, "expected=1", "out="+output, "t=1", "bf1=t"};
		write(root.resolve("call.txt"), "assemble.LiveBranchTrainingCorpus\n"+String.join("\n", args)+"\n");
		final long[] before=gcCounts();
		Throwable failure=null;
		int dataFds=-1, metaFds=-1, outputFds=-1, workers=-1;
		long[] after=null;
		try{
			try{LiveBranchTrainingCorpus.main(args);}catch(RuntimeException|Error caught){failure=caught;}
			dataFds=countReaderFds(data);metaFds=countReaderFds(meta);
			outputFds=trainingOutputFds(output);workers=trainingWriters(output).size();after=gcCounts();
		}finally{
			Throwable rescueFailure=null;
			for(ByteStreamWriter writer : trainingWriters(output)){
				try{check(!writer.poisonAndWait() && !writer.isAlive(), "Training output rescue failed: "+writer.fname);}
				catch(RuntimeException|Error cleanup){
					if(rescueFailure==null){rescueFailure=cleanup;}else if(rescueFailure!=cleanup){rescueFailure.addSuppressed(cleanup);}
				}
			}
			if(rescueFailure!=null){throw new IllegalStateException("Training fixture could not rescue owned writers.", rescueFailure);}
		}
		write(root.resolve("observation.tsv"), "data_fds\tmetadata_fds\toutput_fds\tworkers\tfailure_class\tcause_class\n"+
				dataFds+"\t"+metaFds+"\t"+outputFds+"\t"+workers+"\t"+(failure==null ? "none" : failure.getClass().getName())+
				"\t"+(failure==null || failure.getCause()==null ? "none" : failure.getCause().getClass().getName())+"\n");
		final StringBuilder gc=new StringBuilder("collector\tbefore\tafter\n");
		for(int i=0; i<before.length; i++){gc.append(collectors.get(i).getName()).append('\t').append(before[i]).append('\t').append(after[i]).append('\n');}
		write(root.resolve("gc.tsv"), gc.toString());
		if(failure!=null){
			final StringWriter trace=new StringWriter();final PrintWriter writer=new PrintWriter(trace);
			failure.printStackTrace(writer);writer.flush();write(root.resolve("input-error.txt"), trace.toString());
		}
		check(!measured || Arrays.equals(before, after), "GC changed during training input observation.");
		check(outputFds==0 && workers==0 && trainingOutputFds(output)==0 && trainingWriters(output).isEmpty(),
				"Converter output cleanup must finish before fixture rescue.");
		check(metaFds==0, "Failed metadata acquisition or processing left its descriptor open.");
		if(firstDirectory || secondDirectory){
			check(failure!=null && failure.getClass()==RuntimeException.class && failure.getCause() instanceof FileNotFoundException &&
					failure.getSuppressed().length==0 && failure.getCause().getSuppressed().length==0 &&
					failure.getCause().getMessage().contains((firstDirectory ? data : meta).toString()) &&
					plannerFrame(failure, "assemble.LiveBranchTrainingCorpus", "convert") &&
					plannerFrame(failure, "fileIO.ReadWrite", "getRawInputStream"), "Expected real directory-open failure inside convert: "+failure);
		}else if(malformed){
			check(failure!=null && failure.getClass()==IllegalArgumentException.class &&
					"Wrong provenance width/ordinal at g0001".equals(failure.getMessage()) && failure.getCause()==null &&
					failure.getSuppressed().length==0, "Wrong malformed metadata rejection: "+failure);
		}else{check(healthy && failure==null && dataFds==0, "Healthy training conversion failed or leaked: "+failure);}
		checkBytes(manifest, selected);checkBytes(run.resolve("exit_code.txt"), "0\n");
		checkBytes(labels.resolve("census.tsv"), TRAIN_PAIR_CENSUS);checkBytes(labels.resolve("complete.txt"), TRAIN_PAIR_INPUT_COMPLETE);
		if(firstDirectory){check(Files.isDirectory(data) && data.toFile().list().length==0, "First directory input changed.");}
		else{checkBytes(data, TRAIN_PAIR_VECTOR_HEADER+TRAIN_PAIR_VECTOR);}
		if(secondDirectory){check(Files.isDirectory(meta) && meta.toFile().list().length==0, "Second directory input changed.");}
		else{checkBytes(meta, inputRows);}
		checkBytes(output.resolve("train.tsv"), TRAIN_PAIR_OUTPUT_HEADER);
		checkBytes(output.resolve("validation_controls.tsv"), TRAIN_PAIR_OUTPUT_HEADER);
		checkBytes(output.resolve("validation.tsv"), TRAIN_PAIR_OUTPUT_HEADER+(healthy ? TRAIN_PAIR_VECTOR : ""));
		checkBytes(output.resolve("rows.tsv"), TRAIN_PAIR_OUTPUT_ROWS+(healthy ? "validation\t1\tg0001\t24\t"+TRAIN_PAIR_WORD+"\tgrey\t1\n" : ""));
		checkBytes(output.resolve("census.tsv"), TRAIN_PAIR_OUTPUT_CENSUS+(healthy ? "g0001\tvalidation\t24\t1\t0\t0\t0\t0\t0\t1\t0\n" : ""));
		if(healthy){checkBytes(output.resolve("complete.txt"), TRAIN_PAIR_OUTPUT_COMPLETE);}
		else{check(!Files.exists(output.resolve("complete.txt")), "Failed conversion published completion.");}
		return dataFds;
	}

	/** Collects only live workers for the six exact converter output paths. */
	private static ArrayList<ByteStreamWriter> trainingWriters(final Path output){
		final ArrayList<ByteStreamWriter> writers=new ArrayList<ByteStreamWriter>();
		for(Thread thread : Thread.getAllStackTraces().keySet()){
			if(!(thread instanceof ByteStreamWriter) || !thread.isAlive()){continue;}
			final ByteStreamWriter writer=(ByteStreamWriter)thread;
			for(String name : TRAIN_PAIR_OUTPUTS){if(output.resolve(name).toString().equals(writer.fname)){writers.add(writer);break;}}
		}
		return writers;
	}

	/** Counts writable descriptors for this call, including outputs that failed before creation. */
	private static int trainingOutputFds(final Path output) throws Exception{
		int count=0;
		for(String name : TRAIN_PAIR_OUTPUTS){count+=countPlannerWriterFds(output.resolve(name));}
		return count;
	}

	/** Selects only the last explicit evaluator case, leaving the original argument parser intact. */
	private static boolean evaluatorCase(final String[] args){
		String name=null;
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length==2 && pair[0].equalsIgnoreCase("case")){name=pair[1];}
		}
		return "samek-eval-valid".equals(name) || "samek-eval-gzip-failure".equals(name) ||
				"live-eval-valid".equals(name) || "live-eval-gzip-failure".equals(name);
	}

	/** Runs one evaluator case with inputs=PINNED_DIRECTORY and outdir=NEW; restores all changed globals. */
	private static void runEvaluatorCase(final String[] args) throws Exception{
		String name=null;
		Path root=null, inputs=null;
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			check(pair.length==2 && !pair[1].isEmpty(), "Expected nonempty evaluator flag=value: "+arg);
			if(pair[0].equalsIgnoreCase("case")){name=pair[1];}
			else if(pair[0].equalsIgnoreCase("outdir")){root=Paths.get(pair[1]).toAbsolutePath().normalize();}
			else if(pair[0].equalsIgnoreCase("inputs")){inputs=Paths.get(pair[1]).toAbsolutePath().normalize();}
			else{throw new IllegalArgumentException("Unknown evaluator fixture option: "+arg);}
		}
		check(root!=null && inputs!=null && Files.isDirectory(inputs), "Evaluator requires outdir=NEW and inputs=PINNED_DIRECTORY.");
		final boolean live=name.startsWith("live-");
		Files.createDirectory(root);
		final boolean[] original=backendFlags();
		final int threads=Shared.threads();
		final boolean raw=ReadWrite.RAWMODE, fma=Shared.SIMD_FMA;
		final int remaining;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=ByteFile.FORCE_MODE_BF3=ByteFile.FORCE_MODE_BF4=false;
			Shared.setThreads(1); ReadWrite.RAWMODE=false; Shared.SIMD_FMA=false;
			final FusionJoinReadCleanupTest test=new FusionJoinReadCleanupTest();
			final Path warmup=root.resolve("eval-warmup");
			Files.createDirectory(warmup);
			check(test.observeEvaluator(warmup, inputs, false, false, live)==0, "Healthy evaluator warmup leaked metadata.");
			remaining=test.observeEvaluator(root, inputs, name.equals(live ? "live-eval-gzip-failure" : "samek-eval-gzip-failure"), true, live);
		}finally{
			ByteFile.FORCE_MODE_BF1=original[0]; ByteFile.FORCE_MODE_BF2=original[1];
			ByteFile.FORCE_MODE_BF3=original[2]; ByteFile.FORCE_MODE_BF4=original[3];
			Shared.setThreads(threads); ReadWrite.RAWMODE=raw; Shared.SIMD_FMA=fma;
		}
		final boolean[] restored=backendFlags();
		check(Arrays.equals(original, restored) && Shared.threads()==threads && ReadWrite.RAWMODE==raw && Shared.SIMD_FMA==fma,
				"Evaluator must restore reader flags, threads, raw-mode and SIMD FMA.");
		final StringBuilder flags=new StringBuilder("flag\tbefore\tafter\n");
		for(int i=0; i<4; i++){flags.append("BF").append(i+1).append('\t').append(original[i]).append('\t').append(restored[i]).append('\n');}
		write(root.resolve("backend-flags.tsv"), flags.toString());
		write(root.resolve("settings.tsv"), "setting\tbefore\tduring\tafter\nthreads\t"+threads+"\t1\t"+Shared.threads()+
				"\nrawmode\t"+raw+"\tfalse\t"+ReadWrite.RAWMODE+"\nsimdfma\t"+fma+"\tfalse\t"+Shared.SIMD_FMA+"\n");
		final boolean passed=remaining==0;
		final String report="case\texpected_reader_fds\tobserved_reader_fds\tpassed\n"+name+"\t0\t"+remaining+"\t"+passed+"\n";
		write(root.resolve("checks.tsv"), report); System.out.print(report);
		System.out.println((live ? "LIVE_EVAL_READ case=" : "SAMEK_EVAL_READ case=")+name+" assertions="+
				FusionJoinReadCleanupTest.class.desiredAssertionStatus()+" setup=true flags_restored=true passed="+passed);
		if(!passed){throw new AssertionError("Evaluator left metadata open after feature input rejection: "+remaining);}
	}

	/** Observes real public evaluation after a normal or malformed second gzip acquisition. */
	private int observeEvaluator(final Path root, final Path inputs, final boolean corrupt, final boolean measured, final boolean live) throws Exception{
		for(String file : new String[]{"rows.tsv", "validation.tsv", "model.bbnet"}){Files.copy(inputs.resolve(file), root.resolve(file));}
		final Path meta=root.resolve("rows.tsv"), data=root.resolve("validation.tsv.gz"), model=root.resolve("model.bbnet");
		final Path output=root.resolve("evaluation");
		if(corrupt){write(data, "not gzip\n");}
		else{
			try(OutputStream stream=new GZIPOutputStream(Files.newOutputStream(data, StandardOpenOption.CREATE_NEW))){
				stream.write(Files.readAllBytes(root.resolve("validation.tsv")));
			}
		}
		final byte[] savedData=Files.readAllBytes(data);
		if(!measured){
			final ByteFile probe=ByteFile.makeByteFile(meta.toString(), true);
			try{check(probe.getClass()==ByteFile1.class && countReaderFds(meta)==1, "Evaluator observation requires a live BF1 metadata reader.");}
			finally{check(!probe.close(), "Evaluator metadata probe close failed.");}
			check(countReaderFds(meta)==0, "Evaluator probe remained open after close.");
			write(root.resolve("preflight.tsv"), "backend\treader_fds_live\treader_fds_closed\nByteFile1\t1\t0\n");
		}
		check(countReaderFds(meta)==0 && countReaderFds(data)==0 && !Files.exists(output), "Evaluator paths must be fresh and unopened.");
		final String[] args={"in="+data, "rows="+meta, "net="+model, live ? "partition=validation" : "fold=validation", "out="+output,
				"simdff=false", "simdfma=false", "t=1", "bf1=t"};
		final String caller=live ? "assemble.LiveBranchModelEvaluate" : "assemble.SameKModelEvaluate";
		write(root.resolve("call.txt"), caller+"\n"+String.join("\n", args)+"\n");
		final long[] before=gcCounts();
		RuntimeException failure=null;
		try{if(live){LiveBranchModelEvaluate.main(args);}else{SameKModelEvaluate.main(args);}}
		catch(RuntimeException e){failure=e;}
		final int metaFds=countReaderFds(meta), dataFds=countReaderFds(data);
		int writers=0;
		final String[] outputFiles=live ? new String[]{"scores.tsv", "summary.tsv", "grey.tsv", "control.tsv", "complete.txt"} :
				new String[]{"scores.tsv", "summary.tsv", "candidate.tsv", "control.tsv", "candidate_points.tsv", "control_points.tsv", "complete.txt"};
		for(String file : outputFiles){
			writers+=countPlannerWriterFds(output.resolve(file));
		}
		final long[] after=gcCounts();
		write(root.resolve("eval-observation.tsv"), "metadata_fds\tfeature_fds\twriter_fds\tfailure_class\tcause_class\n"+
				metaFds+"\t"+dataFds+"\t"+writers+"\t"+(failure==null ? "none" : failure.getClass().getName())+"\t"+
				(failure==null || failure.getCause()==null ? "none" : failure.getCause().getClass().getName())+"\n");
		final StringBuilder gc=new StringBuilder("collector\tbefore\tafter\n");
		for(int i=0; i<before.length; i++){gc.append(collectors.get(i).getName()).append('\t').append(before[i]).append('\t').append(after[i]).append('\n');}
		write(root.resolve("gc.tsv"), gc.toString());
		if(failure!=null){
			final StringWriter trace=new StringWriter();
			final PrintWriter writer=new PrintWriter(trace); failure.printStackTrace(writer); writer.flush();
			write(root.resolve("eval-error.txt"), trace.toString());
		}
		check(!measured || Arrays.equals(before, after), "GC during evaluator observation could conceal missing cleanup.");
		for(String file : new String[]{"rows.tsv", "validation.tsv", "model.bbnet"}){
			check(Arrays.equals(Files.readAllBytes(inputs.resolve(file)), Files.readAllBytes(root.resolve(file))), "Evaluator changed input: "+file);
		}
		check(Arrays.equals(savedData, Files.readAllBytes(data)), "Evaluator changed the compressed feature input.");
		check(dataFds==0 && writers==0, "Failed decoder or evaluator output retained an unexpected descriptor.");
		if(corrupt){
			check(failure!=null && failure.getClass()==RuntimeException.class && failure.getCause() instanceof ZipException &&
					plannerFrame(failure, "fileIO.ReadWrite", "getGZipInputStream") &&
					plannerFrame(failure, caller, "run") && plannerFrame(failure, caller, "main"),
					"Expected real feature gzip rejection after metadata acquisition: "+failure);
			check(Files.isDirectory(output), "Evaluator creates its output directory before input acquisition.");
			try(DirectoryStream<Path> files=Files.newDirectoryStream(output)){
				check(!files.iterator().hasNext(), "Second-input failure must leave the created output directory empty.");
			}
		}else{
			check(failure==null && metaFds==0, "Valid evaluator control failed or leaked metadata: "+failure);
			if(live){checkLiveEvaluatorPredictions(output.resolve("scores.tsv"));}
			else{checkEvaluatorPredictions(output.resolve("scores.tsv"));}
			checkBytes(output.resolve("complete.txt"), live ? "LIVE_BRANCH_EVALUATE_PASS partition=validation grey=6 controls=0\n" :
					"SAME_K_EVALUATE_PASS fold=validation candidates=6 controls=2 simdff=false\n");
		}
		return metaFds;
	}

	/** Checks all eight analytic predictions independently; controller also compares all seven output files. */
	private static void checkEvaluatorPredictions(final Path path) throws Exception{
		final List<String> rows=Files.readAllLines(path, StandardCharsets.US_ASCII);
		check(rows.size()==9 && rows.get(0).equals("fold\trow\tgenome\tkind\tevent\ttarget\tscore_bits\tscore"), "Expected eight predictions and exact schema.");
		final float[] scores={1, 0.75f, 0.75f, 0.5f, 0.25f, 0, 0.75f, 0.25f};
		final int[] labels={1, 1, 0, 1, 0, 1, 1, 1};
		final LineParser1 fields=new LineParser1('\t');
		for(int i=0; i<8; i++){
			fields.set(rows.get(i+1).getBytes(StandardCharsets.US_ASCII));
			check(fields.terms()==8 && fields.termEquals("validation", 0) && fields.parseInt(1)==i+1 && fields.termEquals("g0002", 2) &&
					fields.termEquals(i<6 ? "candidate" : "control", 3) && fields.parseInt(4)==(i<6 ? i+1 : i-5) &&
					fields.parseInt(5)==labels[i] && fields.parseInt(6)==Float.floatToRawIntBits(scores[i]) && fields.parseFloat(7)==scores[i],
					"Prediction differs from the independent first-feature linear model: "+rows.get(i+1));
		}
	}

	/** Checks six grey predictions, including K and the literal canonical-word provenance. */
	private static void checkLiveEvaluatorPredictions(final Path path) throws Exception{
		final List<String> rows=Files.readAllLines(path, StandardCharsets.US_ASCII);
		check(rows.size()==7 && rows.get(0).equals("partition\trow\tgenome\tk\tword\trole\ttarget\tscore_bits\tscore"), "Expected six live-branch predictions and exact schema.");
		final float[] scores={1, 0.75f, 0.75f, 0.5f, 0.25f, 0};
		final int[] labels={1, 1, 0, 1, 0, 1};
		final String[] suffixes={"AA", "AC", "AG", "AT", "CA", "CC"};
		final LineParser1 fields=new LineParser1('\t');
		for(int i=0; i<6; i++){
			fields.set(rows.get(i+1).getBytes(StandardCharsets.US_ASCII));
			check(fields.terms()==9 && fields.termEquals("validation", 0) && fields.parseInt(1)==i+1 && fields.termEquals("g0002", 2) &&
					fields.parseInt(3)==24 && fields.termEquals("AAAAAAAAAAAAAAAAAAAAAA"+suffixes[i], 4) && fields.termEquals("grey", 5) &&
					fields.parseInt(6)==labels[i] && fields.parseInt(7)==Float.floatToRawIntBits(scores[i]) && fields.parseFloat(8)==scores[i],
					"Live-branch prediction differs from the independent analytic rows: "+rows.get(i+1));
		}
	}

	/** Selects the isolated processing-finally panel without changing original planner cases. */
	private static boolean plannerCloseCase(final String name){
		return "plan-close-valid".equals(name) || "plan-close-malformed".equals(name) ||
				"plan-close-eof-fault".equals(name) || "plan-close-finally-fault".equals(name);
	}

	/** Captures a real planner call before rescue; only the named diagnostic fault is injected. */
	private static void runPlannerCloseCase(final String name, final Path root) throws Exception{
		check(root!=null && Files.isDirectory(FD_ROOT) && FusionJoinReadCleanupTest.class.desiredAssertionStatus(),
				"Planner close panel requires Linux, -ea and outdir=NEW.");
		Files.createDirectory(root);
		final boolean malformed=name.equals("plan-close-malformed") || name.equals("plan-close-finally-fault");
		final boolean inject=name.endsWith("-fault"), expectedEof=!malformed;
		final Path manifest=root.resolve("manifest.tsv"), output=root.resolve("samples.tsv");
		final Path distribution=root.resolve("samples.tsv.distribution.tsv");
		final String inputText=malformed ? "genomeA\tunused.fa\n" : PLAN_MANIFEST;
		write(manifest, inputText);
		final String[] args={"plan", "in="+manifest, "out="+output, "expected=1", "seed=350194", "t=1", "bf1=t"};
		write(root.resolve("call.txt"), "assemble.FusionJoinDiagnostic\n"+String.join("\n", args)+"\n");
		final ByteArrayOutputStream diagnostic=new ByteArrayOutputStream();
		final PlannerCloseObserver observer=new PlannerCloseObserver(diagnostic, manifest, output, distribution, inject, expectedEof);
		final PrintStream original=System.err;
		final boolean verbose=ReadWrite.verbose, raw=ReadWrite.RAWMODE;
		final boolean[] modes=backendFlags();
		final int threads=Shared.threads(), buffers=Shared.numBuffers();
		Throwable observed=null, problem=null;
		int readers=-1, writers=-1, summaries=-1, workers=-1;
		boolean alive=false, daemon=false, distributionBefore=false;
		String state="ABSENT";
		byte[] beforeOutput=null;
		try{
			System.setErr(observer);ReadWrite.verbose=true;ReadWrite.RAWMODE=false;
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=ByteFile.FORCE_MODE_BF3=ByteFile.FORCE_MODE_BF4=false;
			Shared.setThreads(1);
			try{FusionJoinDiagnostic.main(args);}catch(RuntimeException|Error failure){observed=failure;}
			workers=observer.retainWriters();
			readers=countReaderFds(manifest);writers=countPlannerWriterFds(output);summaries=countPlannerWriterFds(distribution);
			beforeOutput=Files.readAllBytes(output);distributionBefore=Files.exists(distribution);
			if(observer.first!=null){state=observer.first.getState().toString();alive=observer.first.isAlive();daemon=observer.first.isDaemon();}
		}catch(Exception|Error failure){problem=failure;}
		finally{
			try{
				System.setErr(original);ReadWrite.verbose=verbose;ReadWrite.RAWMODE=raw;
				ByteFile.FORCE_MODE_BF1=modes[0];ByteFile.FORCE_MODE_BF2=modes[1];
				ByteFile.FORCE_MODE_BF3=modes[2];ByteFile.FORCE_MODE_BF4=modes[3];
				Shared.setThreads(threads);Shared.setBuffers(buffers);
			}finally{
				try{observer.rescue();}
				catch(RuntimeException|Error cleanup){
					if(observed!=null && observed!=cleanup){cleanup.addSuppressed(observed);}
					if(problem==null){problem=cleanup;}else if(problem!=cleanup){problem.addSuppressed(cleanup);}
				}
			}
		}
		if(problem!=null){throw new IllegalStateException("Planner fixture observation/rescue failed.", problem);}
		check(System.err==original && ReadWrite.verbose==verbose && ReadWrite.RAWMODE==raw &&
				Shared.threads()==threads && Shared.numBuffers()==buffers && Arrays.equals(modes, backendFlags()),
				"Planner close fixture must restore every changed global.");
		check(observer.retainWriters()==0 && countReaderFds(manifest)==0 && countPlannerWriterFds(output)==0 &&
				countPlannerWriterFds(distribution)==0, "Fixture rescue must finish all owned workers and descriptors.");
		write(root.resolve("planner.stderr"), diagnostic.toString("UTF-8"));
		Files.write(root.resolve("output-before-rescue.tsv"), beforeOutput, StandardOpenOption.CREATE_NEW);
		write(root.resolve("close-observation.tsv"), "reader_fds\twriter_fds\tdistribution_fds\tworkers\twriter_state\twriter_alive\twriter_daemon\tdistribution_present\n"+
				readers+"\t"+writers+"\t"+summaries+"\t"+workers+"\t"+state+"\t"+alive+"\t"+daemon+"\t"+distributionBefore+"\n");
		write(root.resolve("close-hook.tsv"), "input_live_fds\tinput_closed_fds\toutput_live_fds\tpost_close_calls\teof_frame\tmarker_thrown\n"+
				observer.inputLive+"\t"+observer.inputClosed+"\t"+observer.outputLive+"\t"+observer.closes+"\t"+observer.eof+"\t"+(observer.marker!=null)+"\n");
		write(root.resolve("close-rescue.tsv"), "reader_fds\twriter_fds\tdistribution_fds\tworkers\tglobals_restored\n0\t0\t0\t0\ttrue\n");
		write(root.resolve("precedence.tsv"), "processing_rejection_expected\tfailure_class\tmarker_identity\tcause_class\tsuppressed_count\n"+
				malformed+"\t"+(observed==null ? "none" : observed.getClass().getName())+"\t"+
				(observed!=null && observed==observer.marker)+"\t"+
				(observed==null || observed.getCause()==null ? "none" : observed.getCause().getClass().getName())+"\t"+
				(observed==null ? 0 : observed.getSuppressed().length)+"\n");
		if(observed!=null){
			final StringWriter trace=new StringWriter();
			final PrintWriter writer=new PrintWriter(trace);observed.printStackTrace(writer);writer.flush();
			write(root.resolve("close-error.txt"), trace.toString());
		}
		check(observer.setupFailure==null, "Invalid planner close hook: "+observer.setupFailure);
		check(observer.closes==1 && observer.inputLive==1 && observer.inputClosed==0 && observer.outputLive==1 &&
				observer.eof==expectedEof, "Planner close must be reached once through the requested EOF/finally path.");
		if(inject){
			check(observed==observer.marker && observed!=null && observed.getClass()==RuntimeException.class &&
					observed.getCause()==null && observed.getSuppressed().length==0,
					"Current close-marker precedence must survive without inventing the earlier processing exception.");
		}else if(malformed){
			check(observed!=null && observed.getClass()==IllegalArgumentException.class &&
					"Corpus references must have seven columns and unique IDs/paths.".equals(observed.getMessage()) &&
					observed.getCause()==null && observed.getSuppressed().length==0,
					"Ordinary malformed-row control must retain its original processing failure.");
		}else{check(observed==null, "Healthy planner failed: "+observed);}
		checkBytes(manifest, inputText);
		final String expectedPlan=malformed ? PLAN_HEADER : PLAN_CLOSE_EXPECTED;
		checkBytes(output, expectedPlan);
		if(!malformed && !inject){checkPlannerOutput(output, distribution);}
		else{check(!distributionBefore && !Files.exists(distribution), "Failure must not produce a distribution or success summary.");}
		check(observer.successes==(!malformed && !inject ? 1 : 0), "Unexpected planner success report.");
		final boolean passed=readers==0 && writers==0 && summaries==0 && workers==0 && !alive &&
				Arrays.equals(beforeOutput, expectedPlan.getBytes(StandardCharsets.US_ASCII));
		final String report="case\treader_fds\twriter_fds\tdistribution_fds\tworkers\tpassed\n"+
				name+"\t"+readers+"\t"+writers+"\t"+summaries+"\t"+workers+"\t"+passed+"\n";
		write(root.resolve("checks.tsv"), report);System.out.print(report);
		System.out.println("FUSION_PLAN_CLOSE case="+name+" assertions=true setup=true rescued=true passed="+passed);
		if(!passed){throw new AssertionError("Planner skipped output cleanup after unchecked input close failure.");}
	}

	/** Synchronous diagnostic observer; faults only after the actual manifest descriptor has closed. */
	private static final class PlannerCloseObserver extends PrintStream{
		/** Borrows no original stderr ownership; all diagnostic bytes go to the supplied buffer. */
		PlannerCloseObserver(final ByteArrayOutputStream sink, final Path manifest_, final Path output_,
				final Path distribution_, final boolean inject_, final boolean expectedEof_) throws Exception{
			super(sink, true, "UTF-8");manifest=manifest_;output=output_;distribution=distribution_;
			inject=inject_;expectedEof=expectedEof_;
		}
		/** Retains the output before setup checks can throw, then raises a distinct marker outside the setup catch. */
		@Override
		public void println(final String text){
			super.println(text);
			if("FUSION_CORPUS_PLAN_PASS genomes=1 seed=350194".equals(text)){successes++;}
			final boolean opening=text.equals("getOutputStream("+output+", false, true, false)");
			final boolean closing=text.startsWith("finishReading(") && text.endsWith(", "+manifest+", false) returned false");
			if(!opening && !closing){return;}
			try{
				final int live=retainWriters();
				if(opening){
					check(++opens==1 && live==0, "Plan output must be newly acquired once.");
					inputLive=countReaderFds(manifest);
					check(inputLive==1, "Output acquisition must follow a real open manifest descriptor.");
					return;
				}
				inputClosed=countReaderFds(manifest);outputLive=countPlannerWriterFds(output);
				check(++closes==1 && inputClosed==0 && outputLive==1 && live==1 && owned.size()==1,
						"Close diagnostic must follow actual input closure while output remains live.");
				first=owned.get(0);
				check(first.fname.equals(output.toString()) && first.isAlive() && !first.isDaemon(),
						"Expected the single real started planner output worker.");
				final Throwable location=new Throwable();
				eof=plannerFrame(location, "fileIO.ByteFile1", "nextLine") || plannerFrame(location, "fileIO.ByteFile1", "nextLineList");
				check(plannerFrame(location, "fileIO.ReadWrite", "finishReading") && plannerFrame(location, "fileIO.ByteFile1", "close") &&
						plannerFrame(location, "assemble.FusionCorpusPlan", "run") && plannerFrame(location, "assemble.FusionJoinDiagnostic", "main") &&
						eof==expectedEof, "Close fault must originate in the exact public planner EOF/finally path.");
			}catch(Exception|Error failure){
				setupFailure=failure;throw new IllegalStateException("Planner close diagnostic setup failed.", failure);
			}
			if(inject){marker=new RuntimeException("Injected planner input-close failure after real closure.");throw marker;}
		}
		/** Finds only this fixture's output/distribution workers; strong references prevent lost ownership. */
		int retainWriters(){
			int count=0;
			for(Thread thread : Thread.getAllStackTraces().keySet()){
				if(!(thread instanceof ByteStreamWriter) || !thread.isAlive()){continue;}
				final ByteStreamWriter writer=(ByteStreamWriter)thread;
				if(writer.fname.equals(output.toString()) || writer.fname.equals(distribution.toString())){
					if(!owned.contains(writer)){owned.add(writer);}count++;
				}
			}
			return count;
		}
		/** Tries every retained writer even after an observation or prior cleanup error. */
		void rescue(){
			Throwable failure=null;
			try{retainWriters();}catch(RuntimeException|Error e){failure=e;}
			for(ByteStreamWriter writer : owned){
				try{check(!writer.poisonAndWait() && !writer.isAlive(), "Planner fixture could not rescue "+writer.fname);}
				catch(RuntimeException|Error e){if(failure==null){failure=e;}else if(failure!=e){failure.addSuppressed(e);}}
			}
			if(failure!=null){throw new IllegalStateException("Planner output rescue failed.", failure);}
		}
		/** Exact fixture paths and selected injection policy; accessed synchronously by the caller. */
		private final Path manifest, output, distribution;
		private final boolean inject, expectedEof;
		private final ArrayList<ByteStreamWriter> owned=new ArrayList<ByteStreamWriter>();
		private ByteStreamWriter first;
		private RuntimeException marker;
		private Throwable setupFailure;
		private int opens, closes, successes, inputLive=-1, inputClosed=-1, outputLive=-1;
		private boolean eof;
	}

	/** Runs only the two planner cases; leaves all original case dispatch and settings intact. */
	private static void runPlannerCase(final String name, final Path root) throws Exception{
		check(root!=null, "Planner case requires outdir=NEW.");
		Files.createDirectory(root);
		final boolean[] original=backendFlags();
		final int threads=Shared.threads();
		final boolean raw=ReadWrite.RAWMODE;
		final int remaining;
		try{
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=ByteFile.FORCE_MODE_BF3=ByteFile.FORCE_MODE_BF4=false;
			Shared.setThreads(1); ReadWrite.RAWMODE=false;
			final FusionJoinReadCleanupTest test=new FusionJoinReadCleanupTest();
			final Path warmup=root.resolve("plan-warmup");
			Files.createDirectory(warmup);
			check(test.observePlanner(warmup, false, false)==0, "Healthy planner warmup leaked its manifest.");
			remaining=test.observePlanner(root, name.equals("plan-output-failure"), true);
		}finally{
			ByteFile.FORCE_MODE_BF1=original[0]; ByteFile.FORCE_MODE_BF2=original[1];
			ByteFile.FORCE_MODE_BF3=original[2]; ByteFile.FORCE_MODE_BF4=original[3];
			Shared.setThreads(threads); ReadWrite.RAWMODE=raw;
		}
		final boolean[] restored=backendFlags();
		check(Arrays.equals(original, restored) && Shared.threads()==threads && ReadWrite.RAWMODE==raw,
				"Planner fixture must restore reader flags, threads and raw-mode setting.");
		final StringBuilder flags=new StringBuilder("flag\tbefore\tafter\n");
		for(int i=0; i<4; i++){flags.append("BF").append(i+1).append('\t').append(original[i]).append('\t').append(restored[i]).append('\n');}
		write(root.resolve("backend-flags.tsv"), flags.toString());
		write(root.resolve("settings.tsv"), "setting\tbefore\tduring\tafter\nthreads\t"+threads+"\t1\t"+Shared.threads()+
				"\nrawmode\t"+raw+"\tfalse\t"+ReadWrite.RAWMODE+"\n");
		final boolean passed=remaining==0;
		final String report="case\texpected_reader_fds\tobserved_reader_fds\tpassed\n"+name+"\t0\t"+remaining+"\t"+passed+"\n";
		write(root.resolve("checks.tsv"), report); System.out.print(report);
		System.out.println("FUSION_PLAN_READ case="+name+" assertions="+
				FusionJoinReadCleanupTest.class.desiredAssertionStatus()+" setup=true flags_restored=true passed="+passed);
		if(!passed){throw new AssertionError("Corpus planner left its manifest open after output rejection: "+remaining);}
	}

	/** Calls the public planner route, proving backend reachability rather than accepting a preflight refusal. */
	private int observePlanner(final Path root, final boolean blocked, final boolean measured) throws Exception{
		final Path manifest=root.resolve("manifest.tsv"), blocker=root.resolve("blocker");
		final Path output=blocked ? blocker.resolve("samples.tsv") : root.resolve("samples.tsv");
		final Path distribution=Paths.get(output+".distribution.tsv");
		write(manifest, PLAN_MANIFEST);
		if(blocked){write(blocker, PLAN_BLOCKER);}
		if(!measured){
			final ByteFile probe=ByteFile.makeByteFile(manifest.toString(), true);
			try{
				check(probe.getClass()==ByteFile1.class && countReaderFds(manifest)==1,
						"Planner observation must see a live synchronous BF1 manifest descriptor.");
			}finally{check(!probe.close(), "Planner probe close failed.");}
			check(countReaderFds(manifest)==0, "Planner probe remained open after close.");
			write(root.resolve("preflight.tsv"), "backend\treader_fds_live\treader_fds_closed\nByteFile1\t1\t0\n");
		}
		check(countReaderFds(manifest)==0 && countPlannerWriterFds(output)==0 && countPlannerWriterFds(distribution)==0,
				"Planner paths must be unopened before the observed call.");
		final String[] args={"plan", "in="+manifest, "out="+output, "expected=1", "seed=350194", "t=1", "bf1=t"};
		write(root.resolve("call.txt"), "assemble.FusionJoinDiagnostic\n"+String.join("\n", args)+"\n");
		final long[] before=gcCounts();
		RuntimeException failure=null;
		try{FusionJoinDiagnostic.main(args);}
		catch(RuntimeException e){failure=e;}
		final int readers=countReaderFds(manifest), writers=countPlannerWriterFds(output), summaries=countPlannerWriterFds(distribution);
		final long[] after=gcCounts();
		write(root.resolve("plan-observation.tsv"), "reader_fds\twriter_fds\tsummary_fds\tfailure_class\tcause_class\n"+
				readers+"\t"+writers+"\t"+summaries+"\t"+(failure==null ? "none" : failure.getClass().getName())+"\t"+
				(failure==null || failure.getCause()==null ? "none" : failure.getCause().getClass().getName())+"\n");
		final StringBuilder gc=new StringBuilder("collector\tbefore\tafter\n");
		for(int i=0; i<before.length; i++){
			gc.append(collectors.get(i).getName()).append('\t').append(before[i]).append('\t').append(after[i]).append('\n');
		}
		write(root.resolve("gc.tsv"), gc.toString());
		if(failure!=null){
			final StringWriter trace=new StringWriter();
			final PrintWriter writer=new PrintWriter(trace);
			failure.printStackTrace(writer); writer.flush();
			write(root.resolve("plan-error.txt"), trace.toString());
		}
		check(!measured || Arrays.equals(before, after), "GC during planner observation could conceal missing cleanup.");
		checkBytes(manifest, PLAN_MANIFEST);
		check(writers==0 && summaries==0, "Planner must leave neither output writer open.");
		if(blocked){
			check(failure!=null && failure.getClass()==RuntimeException.class && failure.getCause() instanceof FileNotFoundException &&
					failure.getCause().getMessage()!=null && failure.getCause().getMessage().contains(output.toString()) &&
					plannerFrame(failure, "fileIO.ReadWrite", "getRawOutputStream") &&
					plannerFrame(failure, "fileIO.ByteStreamWriter", "<init>") &&
					plannerFrame(failure, "assemble.FusionCorpusPlan", "run") &&
					plannerFrame(failure, "assemble.FusionJoinDiagnostic", "main"),
					"Expected real output-open failure after planner input acquisition, not path preflight rejection: "+failure);
			checkBytes(blocker, PLAN_BLOCKER);
			check(!Files.exists(output) && !Files.exists(distribution), "Failed output acquisition created planner output.");
		}else{
			check(failure==null && readers==0, "Valid planner control failed or leaked its input: "+failure);
			checkPlannerOutput(output, distribution);
		}
		return readers;
	}

	/** Requires an exact owner/method frame so an earlier rejection cannot mimic the target failure. */
	private static boolean plannerFrame(final Throwable error, final String owner, final String method){
		for(StackTraceElement frame : error.getStackTrace()){
			if(frame.getClassName().equals(owner) && frame.getMethodName().equals(method)){return true;}
		}
		return false;
	}

	/** Checks schema, preserved metadata, schedule domains and all251 histogram bins without regenerating RNG draws. */
	private static void checkPlannerOutput(final Path output, final Path distribution) throws Exception{
		final List<String> lines=Files.readAllLines(output, StandardCharsets.US_ASCII);
		check(lines.size()==2 && (lines.get(0)+"\n").equals(PLAN_HEADER), "Expected one planner row and exact header.");
		final String row=lines.get(1);
		final LineParser1 fields=new LineParser1('\t'), schedule=new LineParser1(',');
		fields.set(row.getBytes(StandardCharsets.US_ASCII));
		check(fields.terms()==12 && row.startsWith(PLAN_ROW+"\t"), "Planner changed supplied manifest metadata or output width.");
		final int depth=fields.parseInt(7), quality=fields.parseInt(8);
		check(depth>=15 && depth<=150 && quality>=25 && quality<=35 && fields.parseLong(10)>=0 && fields.parseLong(11)==350194,
				"Planner random values must stay within documented domains and preserve the requested seed.");
		schedule.set(fields.parseByteArray(9));
		final int n=schedule.terms();
		check(n==3 || n==4, "Schedule must contain one low, one or two middle and one high K.");
		final int[] counts=new int[128];
		int previous=-1;
		for(int i=0; i<n; i++){
			final int k=schedule.parseInt(i), low=i==0 ? 24 : i==n-1 ? 97 : 49, high=i==0 ? 48 : i==n-1 ? 127 : 96;
			check(k>=low && k<=high && k>previous, "K schedule is not sorted, unique or within its documented bands.");
			counts[k]++; previous=k;
		}
		final List<String> histogram=Files.readAllLines(distribution, StandardCharsets.US_ASCII);
		check(histogram.size()==252 && histogram.get(0).equals("parameter\tvalue\tcount"), "Expected all251 distribution bins and header.");
		int index=1;
		for(int d=15; d<=150; d++){checkPlannerBin(fields, histogram.get(index++), "depth", d, d==depth ? 1 : 0);}
		for(int q=25; q<=35; q++){checkPlannerBin(fields, histogram.get(index++), "quality", q, q==quality ? 1 : 0);}
		for(int k=24; k<=127; k++){checkPlannerBin(fields, histogram.get(index++), "k", k, counts[k]);}
	}

	/** Checks one complete histogram row against the emitted plan, including zero bins. */
	private static void checkPlannerBin(final LineParser1 fields, final String row, final String parameter, final int value, final int expected){
		fields.set(row.getBytes(StandardCharsets.US_ASCII));
		check(fields.terms()==3 && fields.termEquals(parameter, 0) && fields.parseInt(1)==value && fields.parseInt(2)==expected,
				"Planner histogram disagrees with its single emitted row: "+row);
	}

	/** Counts exact writable paths, including a not-yet-created blocked destination. */
	private static int countPlannerWriterFds(final Path path) throws Exception{
		final Path real=path.toFile().getCanonicalFile().toPath();
		int count=0;
		try(DirectoryStream<Path> descriptors=Files.newDirectoryStream(FD_ROOT)){
			for(Path fd : descriptors){
				final Path target;
				try{target=Files.readSymbolicLink(fd);}
				catch(NoSuchFileException disappeared){continue;}
				if(target.equals(real)){
					final int access=fdFlags(fd)&3;
					check(access!=3, "Invalid planner output descriptor access mode.");
					if(access==1 || access==2){count++;}
				}
			}
		}
		return count;
	}

	/** Resolves the existing parser without adding a production access hook. */
	private FusionJoinReadCleanupTest() throws Exception{
		check(Files.isDirectory(FD_ROOT), "Fixture requires Linux /proc/self/fd.");
		parser=FusionJoinDiagnostic.class.getDeclaredMethod("readJoins", String.class, boolean.class);
		parser.setAccessible(true);
		collectors=ManagementFactory.getGarbageCollectorMXBeans();
		check(!collectors.isEmpty(), "GC counters are required to exclude collection-masked cleanup.");
		gcCounts();
	}

	/** Only these four additional cases invoke the audit/sampler entry points. */
	private static boolean pairedCase(final String name){
		return "audit-valid".equals(name) || "audit-gzip-failure".equals(name) ||
				"sample-valid".equals(name) || "sample-gzip-failure".equals(name);
	}

	/** Warms the same production path with valid gzip before the one measured call. */
	private int runPaired(final String name, final Path root) throws Exception{
		preflight(root);
		final boolean audit=name.startsWith("audit-");
		final Path warmup=root.resolve("paired-warmup");
		Files.createDirectory(warmup);
		check(observePair(audit, warmup, false, false)==0, "Healthy paired-reader warmup leaked an input.");
		return observePair(audit, root, name.endsWith("gzip-failure"), true);
	}

	/** Creates real inputs and records descriptors immediately after the public main returns or throws.
	 * The warmup's collector counts are informational; measured calls require no collection.
	 * The failed decoder must close its own input. An unreachable leaked BF1 has no worker;
	 * the fresh process exits after reporting the independent desired zero-FD contract. */
	private int observePair(final boolean audit, final Path root, final boolean corrupt, final boolean measured) throws Exception{
		final Path first=root.resolve("first.tsv"), second=root.resolve("second.tsv.gz");
		final Path output=root.resolve(audit ? "audit.tsv" : "sample");
		final String firstText=audit ? AUDIT_CANDIDATES : SAMPLE_ROWS;
		final String secondText=audit ? AUDIT_LABELS : SAMPLE_DATA;
		write(first, firstText);
		if(corrupt){write(second, "not gzip\n");}
		else{
			try(OutputStream stream=new GZIPOutputStream(Files.newOutputStream(second, StandardOpenOption.CREATE_NEW))){
				stream.write(secondText.getBytes(StandardCharsets.US_ASCII));
			}
		}
		final byte[] savedSecond=Files.readAllBytes(second);
		check(Files.isReadable(first) && Files.isReadable(second), "Both paired inputs must exist and pass readability checks.");
		check(countReaderFds(first)==0 && countReaderFds(second)==0, "New paired inputs already have readers.");
		check(Shared.threads()==1 && !ReadWrite.RAWMODE && ByteFile.FORCE_MODE_BF1 &&
				!ByteFile.FORCE_MODE_BF2 && !ByteFile.FORCE_MODE_BF3 && !ByteFile.FORCE_MODE_BF4,
				"One-thread BF1 with decompression enabled is required for synchronous Java gzip rejection.");
		final String[] args=audit ? new String[]{"in="+first, "labels="+second, "out="+output, "t=1", "bf1=t"} :
				new String[]{"in="+second, "rows="+first, "out="+output, "controls=0", "seed=1", "t=1", "bf1=t"};
		write(root.resolve("call.txt"), (audit ? "assemble.SameKJoinAudit" : "assemble.SameKTrainingSample")+"\n"+
				String.join("\n", args)+"\n");
		final long[] before=gcCounts();
		RuntimeException failure=null;
		try{
			if(audit){SameKJoinAudit.main(args);}else{SameKTrainingSample.main(args);}
		}catch(RuntimeException e){failure=e;}
		final int firstFds=countReaderFds(first), secondFds=countReaderFds(second);
		final long[] after=gcCounts();
		write(root.resolve("paired-observation.tsv"), "first_reader_fds\tsecond_reader_fds\tfailure_class\tcause_class\n"+
				firstFds+"\t"+secondFds+"\t"+(failure==null ? "none" : failure.getClass().getName())+"\t"+
				(failure==null || failure.getCause()==null ? "none" : failure.getCause().getClass().getName())+"\n");
		final StringBuilder gc=new StringBuilder("collector\tbefore\tafter\n");
		for(int i=0; i<before.length; i++){
			gc.append(collectors.get(i).getName()).append('\t').append(before[i]).append('\t').append(after[i]).append('\n');
		}
		write(root.resolve("gc.tsv"), gc.toString());
		if(failure!=null){
			final StringWriter stack=new StringWriter();
			final PrintWriter writer=new PrintWriter(stack);
			failure.printStackTrace(writer); writer.flush();
			write(root.resolve("parse-error.txt"), stack.toString());
		}
		check(!measured || Arrays.equals(before, after), "GC changed during paired-reader observation.");
		checkBytes(first, firstText);
		check(Arrays.equals(Files.readAllBytes(second), savedSecond), "Paired gzip bytes changed.");
		check(secondFds==0, "Second decoder did not close its own input.");
		if(corrupt){
			check(failure!=null && failure.getClass()==RuntimeException.class && failure.getCause() instanceof ZipException,
					"Expected real Java gzip header rejection after the first reader opened: "+failure);
			check(!Files.exists(output), "Second-input rejection must precede output creation.");
		}else{
			check(failure==null && firstFds==0, "Valid paired inputs failed or leaked: "+failure);
			if(audit){checkBytes(output, "k\tsource_code\tdest_code\tsource_junction\tdest_junction\tlabel\tcount\n");}
			else{
				checkBytes(output.resolve("train.tsv"), "#dims\t27\t1\n#schema=samek27-v1; training controls sampled uniformly without replacement\n");
				checkBytes(output.resolve("rows.tsv"), "row\tsource_row\tgenome\tkind\tevent\ttarget\n");
				checkBytes(output.resolve("complete.txt"), "SAME_K_TRAINING_SAMPLE_PASS candidates=0 controls=0 available_controls=0 emitted=0 seed=1\n");
			}
		}
		return firstFds;
	}

	/** Records a real parser outcome before checking the independent zero-reader expectation. */
	private int run(final String name, final Path root) throws Exception{
		preflight(root);
		final String inputText=name.equals("valid") ? VALID_TRACE : name.equals("malformed") ? MALFORMED_TRACE : EMPTY_TRACE;
		final Path input=root.resolve("trace.tsv");
		write(input, inputText);
		check(countReaderFds(input)==0, "New trace already has a readable descriptor.");
		final boolean allowEmpty=name.equals("empty-allowed");
		final long[] before=gcCounts();
		ArrayList<?> joins=null;
		RuntimeException failure=null;
		try{joins=parse(input, allowEmpty);}
		catch(RuntimeException e){failure=e;}
		final int remaining=countReaderFds(input);
		final long[] after=gcCounts();
		final StringBuilder gc=new StringBuilder("collector\tbefore\tafter\n");
		for(int i=0; i<before.length; i++){
			gc.append(collectors.get(i).getName()).append('\t').append(before[i]).append('\t').append(after[i]).append('\n');
		}
		write(root.resolve("gc.tsv"), gc.toString());
		write(root.resolve("observation.tsv"), "allow_empty\treturned_rows\treader_fds_after\n"+
				allowEmpty+"\t"+(joins==null ? -1 : joins.size())+"\t"+remaining+"\n");
		if(failure!=null){
			final StringWriter stack=new StringWriter();
			final PrintWriter writer=new PrintWriter(stack);
			failure.printStackTrace(writer); writer.flush();
			write(root.resolve("parse-error.txt"), stack.toString());
		}
		check(Arrays.equals(before, after), "GC changed during the parser observation; cleanup result is invalid.");
		checkBytes(input, inputText);
		if(name.equals("valid")){
			check(failure==null, "Valid trace unexpectedly failed: "+failure);
			write(root.resolve("parsed.tsv"), verifiedJoin(joins));
		}else if(allowEmpty){
			check(failure==null && joins!=null && joins.isEmpty(), "Allowed empty trace did not return an empty list.");
			write(root.resolve("parsed.tsv"), JOIN_HEADER);
		}else{
			final String expected=name.equals("malformed") ? "Unrecognized fusion trace record." : "Unreadable/empty fusion trace: "+input;
			check(joins==null && failure!=null && failure.getClass()==IllegalArgumentException.class &&
					expected.equals(failure.getMessage()), "Unexpected parser rejection: "+failure);
		}
		return remaining;
	}

	/** Proves BF1 selection, visible first-line ownership, normal close and valid parsing. */
	private void preflight(final Path root) throws Exception{
		final Path input=root.resolve("probe.tsv");
		write(input, VALID_TRACE);
		final ByteFile reader=ByteFile.makeByteFile(input.toString(), true);
		try{
			check(reader.getClass()==ByteFile1.class, "Fixture did not select the synchronous BF1 backend.");
			check(Arrays.equals(reader.nextLine(), "FusionPass\t3".getBytes(StandardCharsets.US_ASCII)), "Unexpected probe first line.");
			check(reader.isOpen() && countReaderFds(input)==1, "Probe must still own one readable FD before its next EOF read.");
		}finally{check(!reader.close(), "Probe reader close reported an I/O failure.");}
		check(countReaderFds(input)==0, "Probe close did not release its descriptor.");
		final String parsed=verifiedJoin(parse(input, false));
		check(countReaderFds(input)==0, "Valid warmup parse did not release its descriptor.");
		checkBytes(input, VALID_TRACE);
		write(root.resolve("probe-parsed.tsv"), parsed);
		write(root.resolve("preflight.tsv"), "backend\treader_fds_live\treader_fds_closed\nByteFile1\t1\t0\n");
	}

	/** Calls the real parser and preserves its unchecked exception identity and cause. */
	private ArrayList<?> parse(final Path input, final boolean allowEmpty) throws Exception{
		try{return (ArrayList<?>)parser.invoke(null, input.toString(), allowEmpty);}
		catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw e;
		}
	}

	/** Checks every retained field against literal independent geometry and sequence expectations. */
	private static String verifiedJoin(final ArrayList<?> joins) throws Exception{
		check(joins!=null && joins.size()==1, "Expected exactly one valid join.");
		final Object join=joins.get(0);
		final String[] integers={"event", "phase", "source", "dest", "overlap", "sourceTrim", "destTrim"};
		final int[] expected={1, 3, 0, 1, 3, 0, 0};
		final StringBuilder row=new StringBuilder(JOIN_HEADER);
		for(int i=0; i<integers.length; i++){
			final Field field=join.getClass().getDeclaredField(integers[i]); field.setAccessible(true);
			final int value=field.getInt(join);
			check(value==expected[i], "Wrong parsed join field: "+integers[i]+"="+value);
			if(i>0){row.append('\t');}row.append(value);
		}
		final String[] arrays={"left", "right", "product"}, sequences={"AACGTAA", "TAACTGC", "AACGTAACTGC"};
		for(int i=0; i<arrays.length; i++){
			final Field field=join.getClass().getDeclaredField(arrays[i]); field.setAccessible(true);
			final byte[] actual=(byte[])field.get(join);
			check(Arrays.equals(actual, sequences[i].getBytes(StandardCharsets.US_ASCII)), "Wrong parsed sequence: "+arrays[i]);
			row.append('\t').append(new String(actual, StandardCharsets.US_ASCII));
		}
		return row.append('\n').toString();
	}

	/** Snapshots all collector counts; missing support cannot establish a no-GC observation. */
	private long[] gcCounts(){
		final long[] counts=new long[collectors.size()];
		for(int i=0; i<counts.length; i++){
			counts[i]=collectors.get(i).getCollectionCount();
			check(counts[i]>=0, "Unavailable GC counter: "+collectors.get(i).getName());
		}
		return counts;
	}

	/** Counts only descriptors readable from this exact input; ignores only vanished readlink entries. */
	private static int countReaderFds(final Path input) throws Exception{
		final Path real=input.toRealPath();
		int count=0;
		try(DirectoryStream<Path> descriptors=Files.newDirectoryStream(FD_ROOT)){
			for(Path fd : descriptors){
				final Path target;
				try{target=Files.readSymbolicLink(fd);}
				catch(NoSuchFileException disappeared){continue;}
				if(target.equals(real)){
					final int access=fdFlags(fd)&3;
					check(access!=3, "Invalid descriptor access mode: "+fd);
					if(access==0 || access==2){count++;}
				}
			}
		}
		return count;
	}

	/** Requires one actual octal access-flags field; all fdinfo errors propagate as invalid setup. */
	private static int fdFlags(final Path fd) throws Exception{
		int flags=-1;
		final Path info=Paths.get("/proc/self/fdinfo", fd.getFileName().toString());
		for(String line : Files.readAllLines(info, StandardCharsets.US_ASCII)){
			if(line.startsWith("flags:")){
				check(flags==-1, "Duplicate descriptor flags: "+info);
				flags=Integer.parseInt(line.substring(6).trim(), 8);
			}
		}
		check(flags>=0, "Missing descriptor flags: "+info);
		return flags;
	}

	/** Captures the complete forced-reader configuration in BF1/BF2/BF3/BF4 order. */
	private static boolean[] backendFlags(){
		return new boolean[]{ByteFile.FORCE_MODE_BF1, ByteFile.FORCE_MODE_BF2, ByteFile.FORCE_MODE_BF3, ByteFile.FORCE_MODE_BF4};
	}

	/** Compares complete tiny input files after parser execution. */
	private static void checkBytes(final Path path, final String expected) throws Exception{
		check(Arrays.equals(Files.readAllBytes(path), expected.getBytes(StandardCharsets.US_ASCII)), "Input bytes changed: "+path);
	}

	/** Saves new evidence without overwriting an earlier observation. */
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}

	/** Setup and result checks remain active with assertions disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Existing static parser entry point, without modified production visibility. */
	private final Method parser;
	/** Collector roster resolved before observing a parser call. */
	private final List<GarbageCollectorMXBean> collectors;
	private static final Path FD_ROOT=Paths.get("/proc/self/fd");
	/** Literal empty-census schemas; expectations do not call production header helpers. */
	private static final String AUDIT_CANDIDATES="event\tk\tsource\tsource_end\tdest\tdest_end\tsource_len\tdest_len\tsource_cov\tdest_cov\tsource_code\tdest_code\tsource_word\tdest_word\tsource_junction\tdest_junction"+
			"\tsource_depth\tsource_rA\tsource_rC\tsource_rG\tsource_rT\tsource_lA\tsource_lC\tsource_lG\tsource_lT"+
			"\tdest_depth\tdest_rA\tdest_rC\tdest_rG\tdest_rT\tdest_lA\tdest_lC\tdest_lG\tdest_lT\tleft_context\tright_context\n";
	private static final String AUDIT_LABELS="# paired-input cleanup control\n"+
			"event\tk\tlabel\ttarget\twhole_hits\tleft_hits\tright_hits\tsource_hits\tsource_unknown\tsource_left_mask\tsource_right_mask\tsource_pair_mask\tdest_hits\tdest_unknown\tdest_left_mask\tdest_right_mask\tdest_pair_mask\n";
	/** Literal zero-row training input; quota zero remains a valid ordinary control. */
	private static final String SAMPLE_ROWS="fold\trow\tgenome\tkind\tevent\ttarget\n";
	private static final String SAMPLE_DATA="#dims\t27\t1\n#schema=samek27-v1; unknown labels excluded; no class balancing\n";
	/** Native grammar with one exact three-base overlap and no trimming. */
	private static final String VALID_TRACE="FusionPass\t3\nOverlapPair\t0\t1\t7\t7\t3\t0\t0\t100\tAACGTAA\tTAACTGC\n";
	/** Recognized record type with too few fields; its error must precede the next EOF read. */
	private static final String MALFORMED_TRACE="FusionPass\t3\nOverlapPair\t0\n";
	/** Nonrecord line preserves the existing empty-trace policy in both modes. */
	private static final String EMPTY_TRACE="# empty trace\n";
	private static final String JOIN_HEADER="event\tphase\tsource\tdest\toverlap\tsource_trim\tdest_trim\tleft\tright\tproduct\n";
	/** Literal planner metadata; referenced genome bytes are deliberately not read by the planner. */
	private static final String PLAN_ROW="genomeA\tunused.fa\t1\tgroupA\t100\t1000\t1";
	private static final String PLAN_MANIFEST="# one-row planner fixture\n\n"+PLAN_ROW+"\n";
	private static final String PLAN_BLOCKER="unchanged planner output blocker\n";
	private static final String PLAN_HEADER="#genome\tsource\ttaxid\tgroup\tcontig_length50\tbases\tscaffolds\tdepth\tquality\tk_schedule\tread_seed\tplan_seed\n";
	/** Saved one-row seed350194 result, independent of the production random generator. */
	private static final String PLAN_CLOSE_EXPECTED=PLAN_HEADER+PLAN_ROW+"\t59\t30\t43,65,84,108\t2287708796041542693\t350194\n";
	/** Literal one-genome validation corpus; no generator or learned model supplies the oracle. */
	private static final String TRAIN_PAIR_WORD="AAAAAAAAAAAAAAAAAAAAAAAA";
	private static final String TRAIN_PAIR_MANIFEST="genome\tfold\tk\tdepth\tseed\taction\trun\tlabels\treads\treference\n";
	private static final String TRAIN_PAIR_ROWS="row\tword\trole\ttarget\n";
	private static final String TRAIN_PAIR_VECTOR_HEADER="#dims\t12\t1\n#schema=livebranch12-v1; one row per known grey-zone/control node; no balancing\n";
	private static final String TRAIN_PAIR_VECTOR="0\t0\t0\t0\t0\t0\t0\t0\t0\t0.5731203\t0\t0\t1\n";
	private static final String TRAIN_PAIR_CENSUS="role\tsafe\tunsafe\tunknown\ngrey\t1\t0\t0\nhard\t0\t0\t0\ncontrol\t0\t0\t0\nexcluded\t0\t0\t0\n";
	private static final String TRAIN_PAIR_INPUT_COMPLETE="LIVE_BRANCH_CORPUS_PASS encounters=1 nodes=1 safe=1 unsafe=0 unknown=0 training=1 grey=1 controls=0\n";
	private static final String TRAIN_PAIR_OUTPUT_HEADER="#dims\t12\t1\n#schema=livebranch12-v1; genome folds fixed; validation candidates and controls separated\n";
	private static final String TRAIN_PAIR_OUTPUT_ROWS="partition\trow\tgenome\tk\tword\trole\ttarget\n";
	private static final String TRAIN_PAIR_OUTPUT_CENSUS="genome\tfold\tk\tgrey_safe\tgrey_unsafe\tgrey_unknown\tcontrol_safe_available\tcontrol_unsafe_available\tcontrol_unknown\tgrey_emitted\tcontrols_emitted\n";
	private static final String TRAIN_PAIR_OUTPUT_COMPLETE="LIVE_BRANCH_TRAINING_PASS genomes=1 train=0 validation=1 validation_controls=0\n";
	private static final String[] TRAIN_PAIR_OUTPUTS={"train.tsv", "validation.tsv", "validation_controls.tsv", "rows.tsv", "census.tsv", "complete.txt"};
}
