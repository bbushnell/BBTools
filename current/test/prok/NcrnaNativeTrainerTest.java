package test.prok;

import java.io.File;
import java.lang.management.ManagementFactory;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

import ml.CellNet;
import ml.CellNetParser;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaNativeTrainer;
import prok.NcrnaNativeTrainer.Method;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Real native trainers plus owned-group failure and interruption cleanup.
 * Synthetic matrices and tiny test networks remain under a fresh out= directory.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaNativeTrainerTest{

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, null, false);
		try(AutoCloseable log=()->Shared.closeStream(pp.outstream)){
			run(pp.args);
		}
	}

	/** Main owns the selected log through parsing, training and failure controls. */
	private static void run(final String[] args) throws Exception{
		final Parser parser=new Parser();
		String installation=null;
		for(String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(), b=split.length>1 ? split[1] : null;
			if(a.equals("installation")){installation=b;}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown fixture argument: "+arg);}
		}
		if(parser.out1==null || installation==null){throw new IllegalArgumentException("Fixture needs out= and installation=");}
		final Path root=Paths.get(parser.out1);
		Files.createDirectory(root);
		final String first=root.resolve("first.tsv.gz").toString(), second=root.resolve("second.tsv.gz").toString();
		final String validation=root.resolve("validation.tsv.gz").toString(), all=root.resolve("all.tsv.gz").toString();
		write(first, matrix(0, 4, false));
		write(second, matrix(4, 8, true));
		write(validation, matrix(8, 10, false));
		write(all, matrix(0, 10, false));
		final String classifier=root.resolve("classifier.config").toString();
		final String regression=root.resolve("regression.config").toString();
		final String classifierRecipe="dims=10,4,1\nseed=-1\nnetseed=3\nbatches=2\nnets=1\ncycles=1\nscan=0\nvfraction=0\nbalance=0\nfpb=1\ntriage=0\nfinal=rslog\nsimd=f\nt=2\nprintinterval=1\now=f\n";
		final String regressionRecipe="dims=10,4,1\nseed=3\nepochs=2\nvfraction=0\nfinal=linear\nsimd=f\nt=1\nprintevery=1\now=t\n";
		write(classifier, classifierRecipe);
		write(regression, regressionRecipe);
		final NcrnaNativeTrainer trainer=new NcrnaNativeTrainer(new File(installation), 2);
		final String[] shards={first, second};
		final File classifierNet=root.resolve("classifier.bbnet").toFile(), regressionNet=root.resolve("regression.bbnet").toFile();
		trainer.train(Method.CLASSIFIER, classifier, shards, validation, classifierNet);
		trainer.train(Method.REGRESSION, regression, shards, validation, regressionNet);
		checkNet(classifierNet);
		checkNet(regressionNet);
		final File refitClassifier=root.resolve("refit_classifier.bbnet").toFile();
		final File refitRegression=root.resolve("refit_regression.bbnet").toFile();
		trainer.train(Method.CLASSIFIER, classifier, new String[]{all}, all, refitClassifier);
		trainer.train(Method.REGRESSION, regression, new String[]{all}, all, refitRegression);
		checkNet(refitClassifier);
		checkNet(refitRegression);

		final String invalid=root.resolve("invalid.config").toString();
		for(final String option:new String[]{"dims", "seed", "epochs", "vfraction", "t"}){
			for(final String value:new String[]{"", "not_a_number"}){
				write(invalid, regressionRecipe+option+"="+value+"\n");
				fails(()->trainer.train(Method.REGRESSION, invalid, shards, validation,
					root.resolve("bad_number.bbnet").toFile()), "Invalid numeric option "+option+"=");
			}
		}
		write(invalid, regressionRecipe.replace("dims=10,4,1", "dims=9,4,1"));
		fails(()->trainer.train(Method.REGRESSION, invalid, shards, validation, root.resolve("bad_dims.bbnet").toFile()), "recipe requires");
		write(invalid, regressionRecipe.replace("epochs=2", "epochs=0"));
		fails(()->trainer.train(Method.REGRESSION, invalid, shards, validation, root.resolve("bad_epochs.bbnet").toFile()), "recipe requires");
		write(invalid, regressionRecipe.replace("vfraction=0", "vfraction=0.1"));
		fails(()->trainer.train(Method.REGRESSION, invalid, shards, validation, root.resolve("bad_split.bbnet").toFile()), "recipe requires");
		write(invalid, regressionRecipe+"in=another_population.tsv\n");
		fails(()->trainer.train(Method.REGRESSION, invalid, shards, validation, root.resolve("bad_input.bbnet").toFile()), "builder-owned");
		write(invalid, classifierRecipe+"netoutbest="+root.resolve("escaped.bbnet")+"\n");
		fails(()->trainer.train(Method.CLASSIFIER, invalid, shards, validation, root.resolve("bad_route.bbnet").toFile()), "builder-owned");
		write(invalid, classifierRecipe+"netseed=-1\n");
		fails(()->trainer.train(Method.CLASSIFIER, invalid, shards, validation, root.resolve("bad_seed.bbnet").toFile()), "recipe requires");
		assert(!Files.exists(root.resolve("escaped.bbnet"))) : "Rejected recipe output aliases must never launch a writer outside the selected product path";

		final Path failed=mock(root, "failed", "cat >/dev/null\nexit 7\n");
		fails(()->new NcrnaNativeTrainer(failed.toFile(), 1).train(Method.REGRESSION, regression,
			new String[]{first}, validation, failed.resolve("net.bbnet").toFile()), "exited 7");
		gone(failed);
		checkStackEnvironment(failed);
		final Path absent=mock(root, "absent", "cat >/dev/null\nexit 0\n");
		fails(()->new NcrnaNativeTrainer(absent.toFile(), 1).train(Method.REGRESSION, regression,
			new String[]{first}, validation, absent.resolve("net.bbnet").toFile()), "without producing");
		gone(absent);

		// More than pipe capacity guarantees the closed input is actually encountered.
		final String large=root.resolve("large.tsv").toString();
		final byte[] block=new byte[65536];
		Arrays.fill(block, (byte)'x');
		try(Output output=new Output(large)){
			for(int i=0; i<32; i++){output.writer.print(block);}
		}
		final Path feed=mock(root, "feed", "exec 0<&-\nexec sleep 60\n");
		fails(()->new NcrnaNativeTrainer(feed.toFile(), 1).train(Method.REGRESSION, regression,
			new String[]{large, first}, validation, feed.resolve("net.bbnet").toFile()), "input stream failed");
		gone(feed);
		// Self-SIGKILL and feeder-triggered SIGKILL produce the same exit status,
		// so this combined-failure case does not race a shell's normal exit cleanup.
		final Path early=mock(root, "early_exit", "kill -KILL \"$$\"\n");
		boolean combined=false;
		try{
			new NcrnaNativeTrainer(early.toFile(), 1).train(Method.REGRESSION, regression,
				new String[]{large, first}, validation, early.resolve("net.bbnet").toFile());
		}catch(IllegalStateException caught){
			combined=caught.getMessage().contains("exited 137") && caught.getMessage().contains("input stream failed")
				&& caught.getCause()!=null;
		}
		assert(combined) : "Early native exit and a broken input pipe must both remain visible";
		gone(early);
		interrupt(root, regression, new String[]{large, first}, validation);
		System.out.println("Native trainers passed: classifier/regression, all-shard stdin, single-file all-data refits, loaded test networks, recipe rejection and owned-group cleanup including a TERM-ignoring descendant");
	}

	/** A literal native matrix; the second shard deliberately omits its final newline. */
	private static String matrix(final int from, final int to, final boolean noNewline){
		assert(to>from) : "Each trainer fixture contains both target classes";
		final ByteBuilder b=new ByteBuilder().append("#dims\t10\t1\n");
		for(int row=from; row<to; row++){
			for(int col=0; col<10; col++){b.append(((row+col)%7)/7.0, 6).tab();}
			b.append(row&1).nl();
		}
		if(noNewline){b.length--;}
		return b.toString();
	}

	/** Loading is a test of the public artifact consumer, not a production self-readback step. */
	private static void checkNet(final File file){
		final CellNet net=CellNetParser.load(file.getPath(), false);
		assert(net.dims[0]==10 && net.dims[net.dims.length-1]==1) : "Every test product must retain the declared feature/target dimensions";
		net.applyInput(new float[10]);
		assert(Float.isFinite(net.feedForward())) : "The actual network consumer must produce a finite prediction";
	}

	private static Path mock(final Path root, final String name, final String body) throws Exception{
		final Path directory=root.resolve(name);
		Files.createDirectories(directory.resolve("dev"));
		write(directory.resolve("dev/regressiontrainer.sh").toString(),
			"#!/usr/bin/env bash\nprintf '%s\\n' \"$$\" > "+quote(directory.resolve("pid.txt").toString())+"\n"+
			"printf '%s\\n' \"$@\" > "+quote(directory.resolve("args.txt").toString())+"\n"+
			"for option in ${JAVA_TOOL_OPTIONS-}; do case \"$option\" in -Xss*) printf '%s\\n' \"$option\";; esac; done > "+
			quote(directory.resolve("stack.txt").toString())+"\n"+body);
		return directory;
	}

	/** Checks only the explicit stack token, never records unrelated inherited JVM options. */
	private static void checkStackEnvironment(final Path directory) throws Exception{
		String expected=null;
		for(String option:ManagementFactory.getRuntimeMXBean().getInputArguments()){
			if(option.startsWith("-Xss")){expected=option;}
		}
		if(expected!=null){
			final java.util.List<String> options=Files.readAllLines(directory.resolve("stack.txt"), StandardCharsets.UTF_8);
			assert(!options.isEmpty() && options.get(options.size()-1).equals(expected))
				: "The parent's explicit stack setting must be the last child JAVA_TOOL_OPTIONS stack token: "+expected;
			System.out.println("Trainer stack environment forwarded: "+expected);
		}
	}

	/** Verifies actual process disappearance, rather than merely observing a thrown exception. */
	private static void gone(final Path directory) throws Exception{
		final String pid;
		try(Input input=new Input(directory.resolve("pid.txt").toString())){
			pid=new String(input.reader.nextLine(), StandardCharsets.US_ASCII);
		}
		assert(Long.parseLong(pid)>0) : "Fixture must record the exact directly owned child PID";
		final Process probe=new ProcessBuilder("kill", "-0", pid)
			.redirectError(directory.resolve("probe.err").toFile()).start();
		assert(probe.waitFor()!=0) : "Failed or interrupted native child must no longer exist: "+pid;
	}

	private static void interrupt(final Path root, final String recipe, final String[] inputs, final String validation) throws Exception{
		final Path descendant=root.resolve("interrupted/descendant.pid");
		final String body="trap '' TERM\nprintf '%s\\n' \"$$\" > "+quote(descendant.toString())+"\nexec sleep 60\n";
		final Path directory=mock(root, "interrupted", "bash -c "+quote(body)+" &\nwait\n");
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final AtomicBoolean interruptPreserved=new AtomicBoolean();
		final Thread worker=new Thread(()->{
			try{
				new NcrnaNativeTrainer(directory.toFile(), 1).train(Method.REGRESSION, recipe,
					inputs, validation, directory.resolve("net.bbnet").toFile());
			}catch(Throwable caught){
				failure.set(caught);
				interruptPreserved.set(Thread.currentThread().isInterrupted());
			}
		});
		worker.start();
		try{
			final long deadline=System.nanoTime()+10_000_000_000L;
			while(!Files.exists(descendant) && worker.isAlive() && System.nanoTime()<deadline){Thread.sleep(10);}
			assert(Files.exists(descendant)) : "The TERM-ignoring descendant must start before cancellation is tested";
			worker.interrupt();
			worker.join(5000);
			assert(!worker.isAlive() && failure.get()!=null && interruptPreserved.get())
				: "Interruption during pipe backpressure must stop the owned group and preserve the caller's interrupt";
			gone(directory);
			descendantStopped(descendant);
		}finally{
			worker.interrupt();
			worker.join(5000);
		}
	}

	/** Orphan zombies are already terminated; the operating system owns their final reaping. */
	private static void descendantStopped(final Path file) throws Exception{
		final String pid=new String(Files.readAllBytes(file), StandardCharsets.US_ASCII).trim();
		assert(pid.matches("[1-9][0-9]*")) : "Probe only the descendant recorded by this fixture";
		final File status=file.resolveSibling("descendant.status").toFile();
		String state="";
		for(int attempt=0; attempt<50; attempt++){
			final Process probe=new ProcessBuilder("ps", "-p", pid, "-o", "stat=").redirectOutput(status).start();
			final int exit=probe.waitFor();
			state=new String(Files.readAllBytes(status.toPath()), StandardCharsets.US_ASCII).trim();
			if(exit==1 && state.isEmpty() || exit==0 && state.startsWith("Z")){return;}
			Thread.sleep(20);
		}
		throw new AssertionError("Owned-group cleanup left the TERM-ignoring descendant alive: "+state);
	}

	private static void fails(final Runnable action, final String message){
		try{action.run();}
		catch(RuntimeException expected){
			assert(expected.getMessage().contains(message)) : "Unexpected failure: "+expected;
			return;
		}
		throw new AssertionError("Expected native-training failure containing: "+message);
	}

	private static void write(final String path, final String text){NcrnaLearningPopulationTest.write(path, text);}
	private static String quote(final String text){return "'"+text.replace("'", "'\\''")+"'";}
}
