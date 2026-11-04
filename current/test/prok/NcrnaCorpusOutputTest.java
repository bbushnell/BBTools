package test.prok;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.function.Function;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.PreParser;
import prok.NcrnaAnnotatedCorpus;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;

/** Controlled output-acquisition failures through a private per-instance factory.
 * Prior writers are real native workers; the fixture owns and closes baseline leaks.
 * Raw output makes exact owned Linux descriptor checks possible without subprocesses.
 * This does not certify shared native writer cancellation or JVM-exit recovery.
 * @author Raiden */
public final class NcrnaCorpusOutputTest{

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaCorpusOutputTest.class, false);
		String root=null, roster=null, expect=null;
		for(String arg:pp.args){
			if(arg.startsWith("outdir=")){root=arg.substring(7);}
			else if(arg.startsWith("roster=")){roster=arg.substring(7);}
			else if(arg.startsWith("expect=")){expect=arg.substring(7);}
			else{throw new IllegalArgumentException("Unknown output fixture argument: "+arg);}
		}
		check(root!=null && roster!=null && ("before".equals(expect) || "after".equals(expect)), "Declare fresh outdir, three-genome GC roster and before/after expectation");
		check(new File(root).mkdir(), "Fixture owns a fresh output root");
		new NcrnaCorpusOutputTest(root, roster, expect.equals("after")).run();
	}

	private NcrnaCorpusOutputTest(String root_, String roster_, boolean after_){root=root_;roster=roster_;after=after_;}
	private void run() throws Exception{
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;
		ReadWrite.setZipThreads(1);ReadWrite.USE_UNPIGZ=false;
		report.append("case\topen_attempts\tconstructed\talive_after_call\topen_descriptors_after_call\tprimary\tsuppressed\n");
		for(int i=1; i<=6; i++){one("open_"+i, i, 0, false, 0, false);}
		one("native_open_4", 4, 0, true, 0, false);
		one("start_3", 0, 3, false, 0, false);
		one("open_and_close_flags", 4, 0, false, 1, false);
		one("open_and_close_exceptions", 4, 0, false, 2, false);
		one("cleanup_only", 0, 0, false, 1, false);
		one("body_and_close_exception", 0, 0, false, 2, true);
		one("ordinary_success", 0, 0, false, 0, false);
		check(cases==13, "All declared ownership controls must execute");
		try(Output out=new Output(FileFormat.testOutput(root+"/observations.tsv", FileFormat.TEXT, null, true, false, false, false))){
			out.writer.print(report);
		}
		System.out.println("NCRNA_CORPUS_OUTPUT_PASS cases="+cases+" expectation="+(after ? "after" : "before"));
	}

	/** Tests the actual private write path after normal native loading/extraction. */
	private void one(String label, int failOpen, int failStart, boolean nativeOpen, int closeMode, boolean body) throws Exception{
		final Factory factory=new Factory(failOpen, failStart, nativeOpen, closeMode);
		final Constructor<NcrnaAnnotatedCorpus> constructor=NcrnaAnnotatedCorpus.class.getDeclaredConstructor(Function.class);
		constructor.setAccessible(true);
		final NcrnaAnnotatedCorpus tool=constructor.newInstance(factory);
		invoke(tool, "load", new Class<?>[]{String.class, int.class, int.class}, roster, 3, 0);
		invoke(tool, "extract", new Class<?>[]{String.class, int.class, int.class, boolean.class}, "RF00174", 2, 1, true);
		final PrintStream prior=System.err;
		final ByteArrayOutputStream diagnostics=new ByteArrayOutputStream();
		Throwable problem=null;int alive=-1, descriptors=-1;
		final boolean raw=ReadWrite.RAWMODE, bgzf=ReadWrite.USE_BGZF, force=ReadWrite.FORCE_BGZIP;
		try(PrintStream capture=new PrintStream(diagnostics, true, "UTF-8")){
			ReadWrite.RAWMODE=true;ReadWrite.USE_BGZF=false;ReadWrite.FORCE_BGZIP=false;System.setErr(capture);
			try{invoke(tool, "write", new Class<?>[]{String.class, String.class, int.class, int.class, boolean.class}, root+"/"+label, "RF00174", body ? 14 : 15, 2, true);}
			catch(Throwable caught){problem=caught;}
			alive=0;for(TrackedWriter writer:factory.writers){if(writer.isAlive()){alive++;}}
			descriptors=openDescriptors(factory.writers);
			final boolean failure=failOpen>0 || failStart>0 || closeMode>0 || body;
			check((problem!=null)==failure, "Invocation outcome must match the controlled failure: "+label+" "+problem);
			if(failOpen>0 || failStart>0){check(problem==factory.primary, "Original acquisition failure identity must survive: "+label);}
			if(after){
				check(alive==0 && descriptors==0, "Every acquired output must close before the call returns: "+label);
				if(closeMode>0 && failOpen>0){check(problem.getSuppressed().length==2, "Both independent cleanup failures must remain under the open error");}
				if(body){check(problem instanceof IllegalArgumentException && problem.getMessage().contains("Extracted raw loci"), "Body failure must remain primary");check(problem.getSuppressed().length==2, "Body failure retains both cleanup exceptions");}
			}else{
				if(failOpen>0 || failStart>0){
					final int started=failOpen>0 ? failOpen-1 : failStart-1;
					check(alive==started && descriptors==factory.writers.size(), "Baseline must expose its known acquired writers/backends: "+label);
				}else if(closeMode==2){check(problem==factory.writers.get(0).cleanup, "Baseline cleanup exception masks the body failure");check(alive==5 && descriptors==5, "Baseline stops closing after its first thrown cleanup");}
			}
			check(factory.attempts==(failOpen>0 ? failOpen : failStart>0 ? failStart : 6), "No later output may open after acquisition fails");
			check(diagnostics.toString("UTF-8").contains("ANNOTATED_CORPUS_PASS")==!failure, "Success marker requires complete cleanup: "+label);
		}finally{
			System.setErr(prior);
			for(TrackedWriter writer:factory.writers){
				try{if(!writer.closed){if(writer.getState()==Thread.State.NEW){writer.closeBeforeStart();}else{writer.poisonAndWait();}}}
				catch(RuntimeException | Error cleanup){check(writer.closed, "Fixture cleanup must complete actual closure before its injected error");}
			}
			ReadWrite.RAWMODE=raw;ReadWrite.USE_BGZF=bgzf;ReadWrite.FORCE_BGZIP=force;
		}
		for(TrackedWriter writer:factory.writers){check(!writer.isAlive() && writer.closed, "Fixture must leave no intentionally leaked worker");}
		check(openDescriptors(factory.writers)==0, "Fixture must leave no owned descriptor");
		report.append(label).tab().append(factory.attempts).tab().append(factory.writers.size()).tab().append(alive).tab().append(descriptors)
			.tab().append(problem==null ? "none" : problem.getClass().getName()).tab().append(problem==null ? 0 : problem.getSuppressed().length).nl();
		cases++;
	}
	private static Object invoke(Object tool, String name, Class<?>[] types, Object... values) throws Exception{
		final Method method=NcrnaAnnotatedCorpus.class.getDeclaredMethod(name, types);method.setAccessible(true);
		try{return method.invoke(tool, values);}
		catch(InvocationTargetException wrapped){
			final Throwable cause=wrapped.getCause();if(cause instanceof Error){throw (Error)cause;}
			if(cause instanceof Exception){throw (Exception)cause;}throw new AssertionError(cause);
		}
	}
	/** Only checks descriptors naming the explicitly constructed writers' exact files. */
	private static int openDescriptors(ArrayList<TrackedWriter> writers) throws Exception{
		int count=0;
		try(DirectoryStream<Path> entries=Files.newDirectoryStream(new File("/proc/self/fd").toPath())){
			for(Path entry:entries){
				try{
					final String target=Files.readSymbolicLink(entry).toString();
					for(TrackedWriter writer:writers){if(target.equals(writer.path)){count++;break;}}
				}catch(NoSuchFileException closedConcurrently){}
			}
		}
		return count;
	}

	private static final class Factory implements Function<String, ByteStreamWriter>{
		Factory(int open, int start, boolean native_, int close){failOpen=open;failStart=start;nativeOpen=native_;closeMode=close;}
		@Override public ByteStreamWriter apply(String path){
			final int index=++attempts;
			if(index==failOpen && !nativeOpen){primary=new IllegalArgumentException("declared open failure");throw (RuntimeException)primary;}
			if(index==failOpen){check(new File(path).mkdir(), "Native rejection must be created after freshness validation");}
			try{
				final TrackedWriter writer=new TrackedWriter(path, index==failStart, index<=2 ? closeMode : 0);
				writers.add(writer);if(index==failStart){primary=writer.startup;}return writer;
			}catch(RuntimeException | Error error){primary=error;throw error;}
		}
		final int failOpen, failStart, closeMode;final boolean nativeOpen;
		final ArrayList<TrackedWriter> writers=new ArrayList<TrackedWriter>();int attempts;Throwable primary;
	}
	/** Injected cleanup outcomes occur AFTER native terminal delivery and real closure. */
	private static final class TrackedWriter extends ByteStreamWriter{
		TrackedWriter(String path_, boolean fail, int mode){
			super(FileFormat.testOutput(path_, FileFormat.TEXT, null, true, false, false, false));
			path=path_;failStart=fail;closeMode=mode;
		}
		@Override public synchronized void start(){if(failStart){throw startup;}super.start();}
		@Override public synchronized boolean closeBeforeStart(){final boolean error=super.closeBeforeStart();closed=true;return finish(error);}
		@Override public boolean poisonAndWait(){final boolean error=super.poisonAndWait();closed=true;return finish(error);}
		private boolean finish(boolean error){if(closeMode==2){throw cleanup;}return error || closeMode==1;}
		final String path;final boolean failStart;final int closeMode;boolean closed;
		final AssertionError startup=new AssertionError("declared startup failure");
		final IllegalStateException cleanup=new IllegalStateException("declared cleanup failure");
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private final String root, roster;private final boolean after;
	private final ByteBuilder report=new ByteBuilder();private int cases;
}
