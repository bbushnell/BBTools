package test.assemble;

import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;

import assemble.SameKTrainingSample;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.Parse;
import parse.PreParser;
import shared.Shared;

/** Checks paired output acquisition through the real sampler in one fresh Linux JVM.
 * The existing ReadWrite diagnostic synchronously creates a directory at a selected
 * output path after canWrite but before FileOutputStream. This is an injected
 * blocked-path witness, not a spontaneous input failure. All observations precede
 * explicit fixture rescue; a failed resource contract is reported only after rescue.
 * No production hook, polling race, timeout or copied sampling algorithm is used.
 * @author Brian Bushnell, Sucrose
 */
public final class SameKTrainingSampleConstructionTest{

	/** Requires case=healthy|first-output-blocked|second-output-blocked, inputs=DIR,
	 * expected=DIR and outdir=NEW. Inputs and expected outputs must be saved tiny fixtures. */
	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, SameKTrainingSampleConstructionTest.class, false);
		try{
			String name=null;
			Path inputs=null, expected=null, root=null;
			for(String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(pair.length==2 && !pair[1].isEmpty(), "Expected nonempty flag=value: "+arg);
				if(pair[0].equalsIgnoreCase("case")){name=pair[1];}
				else if(pair[0].equalsIgnoreCase("inputs")){inputs=path(pair[1]);}
				else if(pair[0].equalsIgnoreCase("expected")){expected=path(pair[1]);}
				else if(pair[0].equalsIgnoreCase("outdir")){root=path(pair[1]);}
				else{throw new IllegalArgumentException("Unknown fixture option: "+arg);}
			}
			check(inputs!=null && expected!=null && root!=null &&
					("healthy".equals(name) || "first-output-blocked".equals(name) || "second-output-blocked".equals(name)),
					"Required case, inputs, expected and new outdir.");
			check(SameKTrainingSampleConstructionTest.class.desiredAssertionStatus(), "This baseline requires -ea.");
			check(Files.isDirectory(FD_ROOT), "Exact writer ownership requires Linux procfs.");
			Files.createDirectory(root);
			new SameKTrainingSampleConstructionTest().run(name, inputs, expected, root);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Runs production once, restores diagnostic/backend settings and rescues before any result report. */
	private void run(final String name, final Path inputs, final Path expected, final Path root) throws Exception{
		final Path output=root.resolve("sample"), train=output.resolve("train.tsv"), rows=output.resolve("rows.tsv");
		final Path complete=output.resolve("complete.txt");
		final Path blocked=name.equals("first-output-blocked") ? train : name.equals("second-output-blocked") ? rows : null;
		final ByteArrayOutputStream diagnostic=new ByteArrayOutputStream();
		final OpenObserver observer=new OpenObserver(diagnostic, output, blocked);
		final PrintStream original=System.err;
		final boolean verbose=ReadWrite.verbose;
		final boolean[] modes={ByteFile.FORCE_MODE_BF1, ByteFile.FORCE_MODE_BF2, ByteFile.FORCE_MODE_BF3, ByteFile.FORCE_MODE_BF4};
		final int threads=Shared.threads(), buffers=Shared.numBuffers();
		Throwable failure=null, problem=null;
		int trainFds=-1, rowsFds=-1, completeFds=-1, workers=-1;
		String state="ABSENT";
		boolean alive=false, daemon=false;
		boolean trainEmpty=false, blockedDirectory=false, completeAbsent=false, rowsAbsent=false;
		try{
			System.setErr(observer);
			ReadWrite.verbose=true;
			ByteFile.FORCE_MODE_BF1=true;
			ByteFile.FORCE_MODE_BF2=ByteFile.FORCE_MODE_BF3=ByteFile.FORCE_MODE_BF4=false;
			Shared.setThreads(1);
			try{
				SameKTrainingSample.main(new String[]{"in="+inputs.resolve("train.tsv"), "rows="+inputs.resolve("rows.tsv"),
						"out="+output, "controls=40", "seed=1", "t=1", "bf1=t"});
			}catch(RuntimeException|Error e){failure=e;}
			workers=observer.retainWorkers();
			trainFds=countWriterFds(train);rowsFds=countWriterFds(rows);completeFds=countWriterFds(complete);
			trainEmpty=Files.isRegularFile(train) && Files.size(train)==0;
			blockedDirectory=blocked!=null && Files.isDirectory(blocked);
			completeAbsent=!Files.exists(complete);rowsAbsent=!Files.exists(rows);
			if(observer.first!=null){
				state=observer.first.getState().toString();
				alive=observer.first.isAlive();daemon=observer.first.isDaemon();
			}
		}catch(Exception|Error e){problem=e;}
		finally{
			// Sampler.main has unwound before restoring stderr: Shared.closeStream excludes the current stream.
			try{
				System.setErr(original);ReadWrite.verbose=verbose;
				ByteFile.FORCE_MODE_BF1=modes[0];ByteFile.FORCE_MODE_BF2=modes[1];
				ByteFile.FORCE_MODE_BF3=modes[2];ByteFile.FORCE_MODE_BF4=modes[3];
				Shared.setThreads(threads);Shared.setBuffers(buffers);
			}finally{
				try{observer.rescue();}
				catch(RuntimeException|Error cleanup){
					if(failure!=null && failure!=cleanup){cleanup.addSuppressed(failure);}
					if(problem==null){problem=cleanup;}else if(problem!=cleanup){problem.addSuppressed(cleanup);}
				}
			}
		}
		if(problem!=null){throw new IllegalStateException("Fixture observation or rescue failed.", problem);}
		check(System.err==original && ReadWrite.verbose==verbose && Shared.threads()==threads && Shared.numBuffers()==buffers &&
				Arrays.equals(modes, new boolean[]{ByteFile.FORCE_MODE_BF1, ByteFile.FORCE_MODE_BF2,
						ByteFile.FORCE_MODE_BF3, ByteFile.FORCE_MODE_BF4}), "Fixture globals were not restored.");
		final int rescued=observer.retainWorkers(), rescuedTrain=countWriterFds(train), rescuedRows=countWriterFds(rows);
		check(rescued==0 && rescuedTrain==0 && rescuedRows==0 && countWriterFds(complete)==0,
				"Fixture rescue must leave no owned workers or writable descriptors.");
		write(root.resolve("sampler.stderr"), diagnostic.toString("UTF-8"));
		write(root.resolve("observation.tsv"), "train_fds\trows_fds\tcomplete_fds\tworkers\tfirst_state\tfirst_alive\tfirst_daemon\n"+
				trainFds+"\t"+rowsFds+"\t"+completeFds+"\t"+workers+"\t"+state+"\t"+alive+"\t"+daemon+"\n");
		write(root.resolve("rescue.tsv"), "train_fds\trows_fds\tworkers\tglobals_restored\n0\t0\t0\ttrue\n");
		write(root.resolve("artifacts.tsv"), "train_empty_before_rescue\tblocker_directory_before_rescue\tcomplete_absent_before_rescue\trows_absent_before_rescue\n"+
				trainEmpty+"\t"+blockedDirectory+"\t"+completeAbsent+"\t"+rowsAbsent+"\n");
		write(root.resolve("hook.tsv"), "train_opens\trows_opens\tcomplete_opens\tblocks\tfirst_live_fds\n"+
				observer.opens[0]+"\t"+observer.opens[1]+"\t"+observer.opens[2]+"\t"+observer.blocks+"\t"+observer.firstLiveFds+"\n");
		if(failure!=null){
			final StringWriter stack=new StringWriter();
			final PrintWriter writer=new PrintWriter(stack);
			failure.printStackTrace(writer);writer.flush();
			write(root.resolve("construction-error.txt"), stack.toString());
		}
		check(observer.setupFailure==null, "Diagnostic hook setup failed: "+observer.setupFailure);
		check(observer.opens[0]==1 && observer.opens[1]==(blocked==train ? 0 : 1) &&
				observer.opens[2]==(blocked==null ? 1 : 0) && observer.blocks==(blocked==null ? 0 : 1),
				"Unexpected diagnostic path; case did not reach the intended backend boundary.");
		if(blocked==null){
			check(failure==null, "Healthy sampler unexpectedly failed: "+failure);
			for(String file : OUTPUTS){
				check(Arrays.equals(Files.readAllBytes(output.resolve(file)), Files.readAllBytes(expected.resolve(file))),
						"Healthy sampler changed saved bytes: "+file);
			}
			check(observer.summaries==1, "Healthy sampler must report its exact 52-row summary once.");
		}else{
			check(failure!=null && failure.getClass()==RuntimeException.class && failure.getCause() instanceof FileNotFoundException &&
					failure.getCause().getMessage().equals(blocked+" (Is a directory)") && failure.getSuppressed().length==0,
					"Expected real blocked-path FileNotFoundException without cleanup failure, observed: "+failure);
			check(blockedDirectory && completeAbsent && Files.isDirectory(blocked) && !Files.exists(complete) && observer.summaries==0,
					"Blocked acquisition must not announce a completed sample.");
			try(DirectoryStream<Path> entries=Files.newDirectoryStream(blocked)){
				check(!entries.iterator().hasNext(), "Injected blocker must remain empty.");
			}
			if(blocked==train){check(rowsAbsent && !Files.exists(rows), "First-open refusal must not attempt the second writer.");}
			else{check(trainEmpty && Files.isRegularFile(train) && Files.size(train)==0, "Acquisition failure precedes all train-row emission.");}
		}
		final boolean passed=trainFds==0 && rowsFds==0 && completeFds==0 && workers==0 && !alive;
		final String report="case\ttrain_fds\trows_fds\tcomplete_fds\tworkers\tpassed\n"+
				name+"\t"+trainFds+"\t"+rowsFds+"\t"+completeFds+"\t"+workers+"\t"+passed+"\n";
		write(root.resolve("checks.tsv"), report);System.out.print(report);
		System.out.println("SAMPLE_CONSTRUCTION case="+name+" assertions=true setup=true rescued=true passed="+passed);
		if(!passed){throw new AssertionError("Sampler left a writer open after acquisition failure.");}
	}

	/** Observes only this sample's three output paths; holds first worker strongly across observation. */
	private static final class OpenObserver extends PrintStream{
		/** Keeps verbose output in a tiny in-memory transcript without owning the original stderr. */
		OpenObserver(final ByteArrayOutputStream sink, final Path root_, final Path blocked_) throws Exception{
			super(sink, true, "UTF-8");root=root_;blocked=blocked_;
		}
		/** This synchronous call precedes the real backend open, after its canWrite assertion. */
		@Override
		public void println(final String text){
			super.println(text);
			if(SUMMARY.equals(text)){summaries++;}
			for(int i=0; i<OUTPUTS.length; i++){
				final Path target=root.resolve(OUTPUTS[i]);
				if(!text.equals("getOutputStream("+target+", false, true, false)")){continue;}
				try{
					final int live=retainWorkers();
					check(++opens[i]==1 && Files.isDirectory(root) && !Files.exists(target),
							"Hook must observe one backend open in a newly created output directory.");
					if(i==1){
						check(live==1 && owned.size()==1 && owned.get(0).fname.equals(root.resolve("train.tsv").toString()),
								"Second open must follow exactly one started training writer.");
						first=owned.get(0);firstLiveFds=countWriterFds(root.resolve("train.tsv"));
						check(first.isAlive() && !first.isDaemon() && firstLiveFds==1,
								"Observer must see the first live non-daemon worker and its real writable descriptor.");
					}else{check(live==0, "First/completion open must not overlap an earlier live writer.");}
					if(target.equals(blocked)){Files.createDirectory(target);blocks++;}
				}catch(Exception|Error e){
					setupFailure=e;
					throw new IllegalStateException("Sampler diagnostic hook setup failed.", e);
				}
			}
		}
		/** Retains exact-path live writers before any checks that could throw, also covering unexpected failure paths. */
		int retainWorkers(){
			int count=0;
			for(Thread thread : Thread.getAllStackTraces().keySet()){
				if(!(thread instanceof ByteStreamWriter) || !thread.isAlive()){continue;}
				final ByteStreamWriter writer=(ByteStreamWriter)thread;
				for(String file : OUTPUTS){
					if(writer.fname.equals(root.resolve(file).toString())){
						if(!owned.contains(writer)){owned.add(writer);}count++;break;
					}
				}
			}
			return count;
		}
		/** Attempts every retained worker even if another reports cleanup failure; never interrupts a queue waiter. */
		void rescue(){
			Throwable failure=null;
			try{retainWorkers();}catch(RuntimeException|Error e){failure=e;}
			for(ByteStreamWriter writer : owned){
				try{check(!writer.poisonAndWait() && !writer.isAlive(), "Could not rescue owned writer: "+writer.fname);}
				catch(RuntimeException|Error e){
					if(failure==null){failure=e;}else if(failure!=e){failure.addSuppressed(e);}
				}
			}
			if(failure!=null){throw new IllegalStateException("Owned writer rescue failed.", failure);}
		}
		/** Exact paths, retained worker references, and synchronous producer-owned observations. */
		private final Path root, blocked;
		private final ArrayList<ByteStreamWriter> owned=new ArrayList<ByteStreamWriter>();
		private final int[] opens=new int[3];
		private ByteStreamWriter first;
		private Throwable setupFailure;
		private int blocks, summaries, firstLiveFds;
	}

	/** Counts only writable descriptors for an exact canonical path; disappearing enumeration entries are ignored. */
	private static int countWriterFds(final Path file) throws Exception{
		final Path real=file.toFile().getCanonicalFile().toPath();
		int count=0;
		try(DirectoryStream<Path> fds=Files.newDirectoryStream(FD_ROOT)){
			for(Path fd : fds){
				final Path target;
				try{target=Files.readSymbolicLink(fd);}catch(NoSuchFileException disappeared){continue;}
				if(!target.equals(real)){continue;}
				int flags=-1;
				for(String line : Files.readAllLines(Paths.get("/proc/self/fdinfo", fd.getFileName().toString()), StandardCharsets.US_ASCII)){
					if(line.startsWith("flags:")){
						check(flags==-1, "Duplicate procfs flags: "+fd);flags=Integer.parseInt(line.substring(6).trim(), 8);
					}
				}
				check(flags>=0 && (flags&3)!=3, "Missing or invalid procfs flags: "+fd);
				if((flags&3)!=0){count++;}
			}
		}
		return count;
	}
	/** Normalizes arguments before exact-path comparisons. */
	private static Path path(final String value){return Paths.get(value).toAbsolutePath().normalize();}
	/** Saves tiny reports only after resource rescue, without replacing earlier evidence. */
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
	}
	/** Setup checks stay active regardless of VM assertion settings. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}
	/** Only these output names are owned by this finite fixture. */
	private static final String[] OUTPUTS={"train.tsv", "rows.tsv", "complete.txt"};
	private static final Path FD_ROOT=Paths.get("/proc/self/fd");
	/** Independent saved fixture has twelve candidates and forty controls, all retained. */
	private static final String SUMMARY="SAME_K_TRAINING_SAMPLE_PASS candidates=12 controls=40 available_controls=40 emitted=52 seed=1";
}
