package test.prok;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Function;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.PreParser;
import parse.LineParser1;
import prok.NcrnaFamilySurvey;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Native survey-output failures with exact owned-thread/descriptor observations.
 * The baseline's intentional leaks are closed by the test, including NEW writers.
 * The optional detail stream is tested before workers and during worker failure.
 * @author Raiden */
public final class NcrnaSurveyOutputTest{

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaSurveyOutputTest.class, false);
		try{runFixture(pp.args);}finally{Shared.closeStream(pp.outstream);}
	}
	private static void runFixture(String[] args) throws Exception{
		String root=null, input=null, expectation=null;
		for(String arg:args){
			if(arg.startsWith("outdir=")){root=arg.substring(7);}
			else if(arg.startsWith("inputroot=")){input=arg.substring(10);}
			else if(arg.startsWith("expect=")){expectation=arg.substring(7);}
			else{throw new IllegalArgumentException("Unknown survey fixture option: "+arg);}
		}
		check(root!=null && new File(root).mkdir() && ("before".equals(expectation) || "after".equals(expectation)), "Declare a fresh output root and before/after expectation");
		if(input==null){input=root+"/input";NcrnaFamilySurveyTest.main(new String[]{input});}
		check(new File(input, "map.tsv").isFile(), "Use the existing three-genome survey fixture map");
		new NcrnaSurveyOutputTest(root, input, expectation.equals("after")).run();
	}
	private NcrnaSurveyOutputTest(String root_, String input_, boolean after_){root=root_;input=input_;after=after_;}
	private void run() throws Exception{
		Shared.setThreads(1);ReadWrite.setZipThreads(1);ReadWrite.USE_UNPIGZ=false;
		report.append("case\tattempts\tconstructed\talive_after_call\topen_descriptors_after_call\tprimary\tsuppressed\n");
		for(int n=1; n<=6; n++){one("open_"+n, n, 0, false, 0, false, false, false, false);}
		one("native_open_4", 4, 0, true, 0, false, false, false, false);
		one("start_4", 0, 4, false, 0, false, false, false, false);
		one("open_and_close_flags", 5, 0, false, 1, false, false, false, false);
		one("open_and_close_exceptions", 5, 0, false, 2, false, false, false, false);
		one("cleanup_only", 0, 0, false, 2, false, false, false, false);
		one("body_and_cleanup", 0, 0, false, 2, true, false, false, false);
		one("ordinary_success", 0, 0, false, 0, false, false, false, false);
		one("detail_start", 0, 1, false, 0, false, true, false, false);
		one("detail_header", 0, 0, false, 0, false, true, true, false);
		one("detail_worker_and_cleanup", 0, 0, false, 2, false, true, false, true);
		one("detail_success", 0, 0, false, 0, false, true, false, false);
		check(cases==17, "All declared survey ownership cases must execute");
		try(Output out=new Output(FileFormat.testOutput(root+"/observations.tsv", FileFormat.TEXT, null, true, false, false, false))){
			out.writer.print(report);
		}
		System.out.println("NCRNA_SURVEY_OUTPUT_PASS cases="+cases+" expectation="+(after ? "after" : "before"));
	}

	private void one(String label, int open, int start, boolean nativeOpen, int close, boolean body, boolean detail, boolean header, boolean worker) throws Exception{
		final Factory factory=new Factory(open, start, nativeOpen, close, body, detail, header);
		final Constructor<NcrnaFamilySurvey> constructor=NcrnaFamilySurvey.class.getDeclaredConstructor(Function.class);constructor.setAccessible(true);
		final NcrnaFamilySurvey tool=constructor.newInstance(factory);factory.owner=tool;
		invoke(tool, "load", new Class<?>[]{String.class, String.class, int.class, int.class}, input+"/map.tsv", worker ? root+"/absent_input" : input, 3, 0);
		final PrintStream prior=System.err;final ByteArrayOutputStream diagnostics=new ByteArrayOutputStream();
		Throwable problem=null;int alive=-1, descriptors=-1;
		try(PrintStream capture=new PrintStream(diagnostics, true, "UTF-8")){
			System.setErr(capture);
			try{invoke(tool, "process", new Class<?>[]{String.class, int.class, boolean.class}, root+"/"+label, 1, detail);}
			catch(Throwable caught){problem=caught;}
			alive=0;for(TrackedWriter writer:factory.writers){if(writer.isAlive()){alive++;}}
			descriptors=openDescriptors(factory.writers);
			final boolean failure=open>0 || start>0 || close>0 || body || header || worker;
			check((problem!=null)==failure, "Invocation must match its controlled failure: "+label+" "+problem);
			if(open>0 || start>0 || header){check(problem==factory.primary, "Original acquisition/header error identity must survive: "+label);}
			if(after){
				check(alive==0 && descriptors==0, "All acquired survey outputs must close: "+label);
				if(open>0 && close>0){check(problem.getSuppressed().length==2, "Both cleanup failures must remain under acquisition failure");}
				if(body){check(problem instanceof IllegalArgumentException && problem.getMessage().contains("Every declared genome"), "Body conservation error must remain primary");check(problem.getSuppressed().length==2, "Both table-cleanup failures must remain visible");}
				if(worker){check(problem instanceof IllegalArgumentException && problem.getMessage().contains("GFF worker failed"), "Worker failure must precede detail cleanup");check(problem.getSuppressed().length==1, "Detail cleanup failure must remain visible");}
			}else{
				if(open>0){check(alive==Math.max(0, open-3) && descriptors==alive, "Earlier alias/exclusion outputs close; only table group leaks: "+label);}
				else if(start>0){check(alive==(detail ? 0 : 1) && descriptors==alive+1, "Constructed unstarted writer and earlier active table must remain visible");}
				else if(header){check(alive==1 && descriptors==1, "Baseline detail header fails before cleanup protection");}
				else if(close>0){
					final TrackedWriter first=factory.writers.get(detail ? 0 : 2);
					check(problem==first.cleanup, "Baseline close masks any body/worker error");
					check(alive==(detail ? 0 : 3) && descriptors==alive, "Baseline stops after its first throwing close");
				}
			}
			final int attempts=open>0 ? open : start>0 ? start : header || worker ? 1 : detail ? 7 : 6;
			check(factory.attempts==attempts, "No output opens after the controlled acquisition failure");
			check(diagnostics.toString("UTF-8").contains("NCRNA_SURVEY_PASS")==!failure, "Success requires complete output ownership");
		}finally{
			System.setErr(prior);
			for(TrackedWriter writer:factory.writers){
				try{if(!writer.closed){if(writer.getState()==Thread.State.NEW){writer.closeBeforeStart();}else{writer.poisonAndWait();}}}
				catch(RuntimeException | Error cleanup){check(writer.closed, "Injected cleanup error must follow actual native closure");}
			}
		}
		for(TrackedWriter writer:factory.writers){check(writer.closed && !writer.isAlive(), "Fixture must clean every intentional baseline leak");}
		check(openDescriptors(factory.writers)==0, "Fixture must leave no owned descriptor");
		if(after && label.equals("detail_success")){verifyDuplicateLine(root+"/"+label+".loci.tsv");}
		report.append(label).tab().append(factory.attempts).tab().append(factory.writers.size()).tab().append(alive).tab().append(descriptors)
			.tab().append(problem==null ? "none" : problem.getClass().getName()).tab().append(problem==null ? 0 : problem.getSuppressed().length).nl();cases++;
	}
	/** Literal input has duplicate components at lines4 and7; provenance must use4. */
	private static void verifyDuplicateLine(String path) throws Exception{
		final LineParser1 parser=new LineParser1('\t');int found=0;
		try(Input input=new Input(path)){
			check(input.reader.nextLine()!=null, "Detail output requires its header");
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				parser.set(row);check(parser.terms()==9, "Detail provenance retains the nine-column schema");
				if(parser.termEquals("g1", 0) && parser.termEquals("RF00504", 1) && parser.termEquals("1-40,61-100", 4)){
					check(parser.parseLong(7)==4, "Both duplicate components support this locus; the earliest source is literal GFF line4, not line7");found++;
				}
			}
		}
		check(found==1, "Duplicate collapse must retain exactly one glycine locus with earliest provenance");
	}
	/** Test-owned native input; resource scope retains a parsing/assertion failure. */
	private static final class Input implements AutoCloseable{
		Input(String path){reader=ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.TEXT, null, false, true));}
		@Override public void close(){check(!reader.close(), "Detail verification input must close without a native I/O error");}
		final ByteFile reader;
	}
	private static Object invoke(Object tool, String name, Class<?>[] types, Object... values) throws Exception{
		final Method method=NcrnaFamilySurvey.class.getDeclaredMethod(name, types);method.setAccessible(true);
		try{return method.invoke(tool, values);}
		catch(InvocationTargetException wrapped){final Throwable cause=wrapped.getCause();if(cause instanceof Error){throw (Error)cause;}if(cause instanceof Exception){throw (Exception)cause;}throw new AssertionError(cause);}
	}
	private static int openDescriptors(ArrayList<TrackedWriter> writers) throws Exception{
		assert(writers!=null):"Descriptor counts must be restricted to this fixture's acquired native writers";
		int count=0;
		try(DirectoryStream<Path> entries=Files.newDirectoryStream(new File("/proc/self/fd").toPath())){
			for(Path entry:entries){try{
				final String target=Files.readSymbolicLink(entry).toString();
				for(TrackedWriter writer:writers){if(target.equals(writer.path)){count++;break;}}
			}catch(NoSuchFileException closedConcurrently){}}
		}return count;
	}
	private static final class Factory implements Function<String, ByteStreamWriter>{
		Factory(int open, int start, boolean native_, int close, boolean body_, boolean detail_, boolean header_){failOpen=open;failStart=start;nativeOpen=native_;closeMode=close;body=body_;detail=detail_;header=header_;}
		@Override public ByteStreamWriter apply(String path){
			final int n=++attempts;
			if(n==failOpen && !nativeOpen){primary=new IllegalArgumentException("declared open failure");throw (RuntimeException)primary;}
			if(n==failOpen){check(new File(path).mkdir(), "Real open rejection occurs after freshness checks");}
			if(body && n==6){
				try{final Field field=NcrnaFamilySurvey.class.getDeclaredField("results");field.setAccessible(true);Arrays.fill((Object[])field.get(owner), null);}
				catch(ReflectiveOperationException error){throw new IllegalStateException("Fixture must exercise actual missing-result validation", error);}
			}
			try{
				final int cleanup=(detail ? n==1 : n==3 || n==4) ? closeMode : 0;
				final TrackedWriter writer=new TrackedWriter(path, n==failStart, header && n==1, cleanup);writers.add(writer);
				if(n==failStart){primary=writer.startup;}else if(header && n==1){primary=writer.headerFailure;}return writer;
			}catch(RuntimeException | Error error){primary=error;throw error;}
		}
		final int failOpen, failStart, closeMode;final boolean nativeOpen, body, detail, header;
		final ArrayList<TrackedWriter> writers=new ArrayList<TrackedWriter>();int attempts;Throwable primary;NcrnaFamilySurvey owner;
	}
	private static final class TrackedWriter extends ByteStreamWriter{
		TrackedWriter(String path_, boolean start, boolean header, int close){super(FileFormat.testOutput(path_, FileFormat.TEXT, null, true, false, false, false));path=path_;failStart=start;failHeader=header;closeMode=close;}
		@Override public synchronized void start(){if(failStart){throw startup;}super.start();}
		@Override public ByteStreamWriter println(String line){if(failHeader){throw headerFailure;}return super.println(line);}
		@Override public synchronized boolean closeBeforeStart(){final boolean error=super.closeBeforeStart();closed=true;return finish(error);}
		@Override public boolean poisonAndWait(){final boolean error=super.poisonAndWait();closed=true;return finish(error);}
		private boolean finish(boolean error){if(closeMode==2){throw cleanup;}return error || closeMode==1;}
		final String path;final boolean failStart, failHeader;final int closeMode;boolean closed;
		final AssertionError startup=new AssertionError("declared startup failure");
		final IllegalArgumentException headerFailure=new IllegalArgumentException("declared detail header failure");
		final IllegalStateException cleanup=new IllegalStateException("declared cleanup failure");
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private final String root, input;private final boolean after;private int cases;
	private final ByteBuilder report=new ByteBuilder();
}
