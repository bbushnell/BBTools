package test.template;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.FileFormat;
import shared.Timer;
import stream.Read;
import stream.SamLine;
import stream.Streamer;
import stream.Writer;
import structures.ListNum;
import template.A_SampleStreamer;

/** Completion-status checks for Streamer template reader and writer finalization.
 * Prepare mode writes an instrumented A_SampleStreamerMT copy that redirects only
 * concrete factory calls; run mode checks precompiled copies with finite fakes.
 * @author Jean
 */
public final class TemplateStreamerCompletionCase{

	/** Prepares the MT source copy or runs finite completion-status cases. */
	public static void main(String[] args) throws Exception{
		final HashMap<String, String> map=parseArgs(args);
		final String mode=map.get("mode");
		if("prepare".equals(mode)){
			prepare(Paths.get(required(map, "source")), Paths.get(required(map, "out")));
		}else if("run".equals(mode)){
			run(Paths.get(required(map, "classes")), Paths.get(required(map, "work")),
					required(map, "expect"));
		}else{
			throw new IllegalArgumentException("Expected mode=prepare or mode=run");
		}
	}

	/** Writes a source-instrumented A_SampleStreamerMT copy under out/template. */
	private static void prepare(final Path source, final Path out) throws Exception{
		String text=new String(Files.readAllBytes(
				source.resolve("template/A_SampleStreamerMT.java")), StandardCharsets.UTF_8);
		text=text.replace("A_SampleStreamerMT", "TemplateStreamerCompletionMT");
		text=replaceOnce(text, "StreamerFactory.makeStreamer",
				"test.template.TemplateStreamerCompletionCase.makeStreamer");
		text=replaceOnce(text, "WriterFactory.makeWriter",
				"test.template.TemplateStreamerCompletionCase.makeWriter");

		final Path javaFile=out.resolve("template/TemplateStreamerCompletionMT.java");
		Files.createDirectories(javaFile.getParent());
		check(!Files.exists(javaFile), "prepared source already exists: "+javaFile);
		Files.write(javaFile, text.getBytes(StandardCharsets.UTF_8));
		System.out.println("prepared\t"+javaFile);
	}

	/** Runs completion-status cases against precompiled sample classes. */
	private static void run(final Path classes, final Path work, final String expect)
			throws Exception{
		check("baseline".equals(expect) || "fixed".equals(expect),
				"expect=baseline or expect=fixed is required");
		check(!Files.exists(work), "Fresh work directory required: "+work);
		Files.createDirectories(work);
		try(ChildFirstLoader loader=new ChildFirstLoader(classes.toUri().toURL(),
				TemplateStreamerCompletionCase.class.desiredAssertionStatus())){
			final Class<?> mtClass=Class.forName("template.TemplateStreamerCompletionMT",
					true, loader);
			single(work.resolve("single-normal"), expect, false, false);
			single(work.resolve("single-writer-close-error"), expect, true, false);
			single(work.resolve("single-reader-close-error"), expect, false, true);
			mt(mtClass, work.resolve("mt-normal"), expect, false, false);
			mt(mtClass, work.resolve("mt-writer-close-error"), expect, true, false);
			mt(mtClass, work.resolve("mt-reader-close-error"), expect, false, true);
		}
	}

	/** Checks the direct single-consumer Streamer seam. */
	private static void single(final Path work, final String expect, final boolean writerError,
			final boolean readerError) throws Exception{
		Files.createDirectories(work);
		final String name=work.getFileName().toString();
		final A_SampleStreamer sample=new A_SampleStreamer(new String[] {
				"in="+writeFastq(work.resolve(name+".fq")), "t=1"});
		final FakeStreamer streamer=new FakeStreamer(readerError);
		final FakeWriter writer=new FakeWriter(writerError);
		final Method process=A_SampleStreamer.class.getDeclaredMethod("process",
				Streamer.class, Writer.class, Timer.class, Boolean.TYPE);
		process.setAccessible(true);
		final Throwable thrown=invoke(process, sample, streamer, writer, new Timer(), Boolean.TRUE);
		checkExpect(name, expect, writerError || readerError, thrown);
		check(streamer.closed, name+" should close the reader");
		check(writer.started, name+" should start the writer");
		check(writer.finished, name+" should finish the writer");
		check(writer.readsWritten==1, name+" should write 1 read, saw "+writer.readsWritten);
		check(writer.basesWritten==4, name+" should write 4 bases, saw "+writer.basesWritten);
		pass(name, thrown);
	}

	/** Checks A_SampleStreamerMT with only its factory calls redirected to fakes. */
	private static void mt(final Class<?> type, final Path work, final String expect,
			final boolean writerError, final boolean readerError) throws Exception{
		Files.createDirectories(work);
		configureFakes(writerError, readerError);
		final String name=work.getFileName().toString();
		final Object sample=type.getConstructor(String[].class).newInstance((Object)new String[] {
				"in="+writeFastq(work.resolve(name+".fq")), "out="+work.resolve(name+".out.fq"),
				"t=1"});
		final Method process=type.getDeclaredMethod("process", Timer.class);
		process.setAccessible(true);
		final Throwable thrown=invoke(process, sample, new Timer());
		final FakeStreamer streamer=lastStreamer.get();
		final FakeWriter writer=lastWriter.get();
		checkExpect(name, expect, writerError || readerError, thrown);
		check(streamer!=null && streamer.closed, name+" should close the reader");
		check(writer!=null && writer.started, name+" should start the writer");
		check(writer!=null && writer.finished, name+" should finish the writer");
		check(writer!=null && writer.readsWritten==1,
				name+" should write 1 read, saw "+(writer==null ? -1 : writer.readsWritten));
		check(writer!=null && writer.basesWritten==4,
				name+" should write 4 bases, saw "+(writer==null ? -1 : writer.basesWritten));
		pass(name, thrown);
	}

	/** Redirect target matching StreamerFactory.makeStreamer for A_SampleStreamerMT. */
	public static Streamer makeStreamer(final FileFormat ff1, final FileFormat ff2,
			final boolean ordered, final long maxReads, final boolean saveHeader,
			final boolean makeReads, final int threads){
		final FakeStreamer fake=new FakeStreamer(configuredReaderError);
		lastStreamer.set(fake);
		return fake;
	}

	/** Redirect target matching WriterFactory.makeWriter for A_SampleStreamerMT. */
	public static Writer makeWriter(final FileFormat ffout1, final FileFormat ffout2,
			final int threads, final ArrayList<byte[]> header, final boolean useSharedHeader){
		final FakeWriter fake=new FakeWriter(configuredWriterError);
		lastWriter.set(fake);
		return fake;
	}

	/** Configures fake factories for one instrumented run. */
	private static void configureFakes(final boolean writerError, final boolean readerError){
		configuredWriterError=writerError;
		configuredReaderError=readerError;
		lastStreamer.set(null);
		lastWriter.set(null);
	}

	/** Invokes a reflected method and unwraps the real throwable. */
	private static Throwable invoke(final Method method, final Object target, final Object... args)
			throws Exception{
		try{
			method.invoke(target, args);
			return null;
		}catch(InvocationTargetException e){
			return e.getCause();
		}
	}

	/** Verifies baseline/fixed thrown state for one row. */
	private static void checkExpect(final String name, final String expect,
			final boolean error, final Throwable thrown){
		if("baseline".equals(expect)){
			check(thrown==null, "baseline should miss "+name+", got "+thrown);
		}else if(error){
			check(thrown!=null && thrown.getClass()==RuntimeException.class,
					"fixed should throw RuntimeException for "+name+", got "+thrown);
			check(thrown.getMessage()!=null
					&& thrown.getMessage().contains("terminated in an error state"),
					"fixed should report termination for "+name+", got "+thrown.getMessage());
		}else{
			check(thrown==null, "fixed normal should not throw "+thrown);
		}
	}

	/** Parses key=value arguments while preserving later '=' characters. */
	private static HashMap<String, String> parseArgs(final String[] args){
		final HashMap<String, String> map=new HashMap<String, String>();
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			check(split.length==2, "Arguments require key=value: "+arg);
			map.put(split[0], split[1]);
		}
		return map;
	}

	/** Returns a required key from the argument map. */
	private static String required(final HashMap<String, String> map, final String key){
		final String value=map.get(key);
		check(value!=null && value.length()>0, "Missing "+key);
		return value;
	}

	/** Replaces exactly one source seam. */
	private static String replaceOnce(final String text, final String target,
			final String replacement){
		final int first=text.indexOf(target);
		check(first>=0, "Missing target: "+target);
		check(text.indexOf(target, first+target.length())<0, "Duplicate target: "+target);
		return text.substring(0, first)+replacement+text.substring(first+target.length());
	}

	/** Writes a tiny valid FASTQ input and returns its path. */
	private static Path writeFastq(final Path path) throws Exception{
		Files.write(path, FASTQ);
		return path;
	}

	/** Creates one tiny read batch. */
	private static ListNum<Read> batch(){
		final ArrayList<Read> list=new ArrayList<Read>(1);
		list.add(new Read("ACGT".getBytes(StandardCharsets.UTF_8),
				"!!!!".getBytes(StandardCharsets.UTF_8), "r0", 0));
		return new ListNum<Read>(list, 0);
	}

	/** Prints one successful case row. */
	private static void pass(final String name, final Throwable thrown){
		System.out.println("PASS\t"+name+"\t"+thrownName(thrown));
	}

	/** Formats a captured throwable type. */
	private static String thrownName(final Throwable thrown){
		return thrown==null ? "none" : thrown.getClass().getName();
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Loads the prepared MT sample before any same-named class from the parent loader. */
	private static final class ChildFirstLoader extends URLClassLoader{
		ChildFirstLoader(final URL url, final boolean assertions){
			super(new URL[] {url}, TemplateStreamerCompletionCase.class.getClassLoader());
			setDefaultAssertionStatus(assertions);
			setPackageAssertionStatus("template", assertions);
		}

		@Override
		protected Class<?> loadClass(final String name, final boolean resolve)
				throws ClassNotFoundException{
			if(name.equals("template.TemplateStreamerCompletionMT")
					|| name.startsWith("template.TemplateStreamerCompletionMT$")){
				synchronized(getClassLoadingLock(name)){
					Class<?> clazz=findLoadedClass(name);
					if(clazz==null){clazz=findClass(name);}
					if(resolve){resolveClass(clazz);}
					return clazz;
				}
			}
			return super.loadClass(name, resolve);
		}
	}

	/** Finite Read streamer with a controlled post-close error flag. */
	public static final class FakeStreamer implements Streamer{

		/** Configures close status and one finite batch. */
		FakeStreamer(final boolean closeError_){
			closeError=closeError_;
			next=batch();
		}

		@Override
		public String fname(){return "fake.fq";}
		@Override
		public void start(){started=true;}
		@Override
		public void close(){closed=true;}
		@Override
		public boolean paired(){return false;}
		@Override
		public int pairnum(){return 0;}
		@Override
		public long readsProcessed(){return 1;}
		@Override
		public long basesProcessed(){return 4;}
		@Override
		public void setSampleRate(final float rate, final long seed){}
		@Override
		public synchronized ListNum<Read> nextList(){
			final ListNum<Read> ln=next;
			next=null;
			return ln;
		}
		@Override
		public ListNum<SamLine> nextLines(){return null;}
		@Override
		public boolean hasMore(){return next!=null;}
		@Override
		public boolean errorState(){return closed && closeError;}

		/** Controlled close status. */
		final boolean closeError;
		/** Next finite batch, then null. */
		private ListNum<Read> next;
		/** Whether start was called. */
		boolean started=false;
		/** Whether close was called. */
		boolean closed=false;
	}

	/** Finite Writer with a controlled poisonAndWait status. */
	public static final class FakeWriter implements Writer{

		/** Configures completion status. */
		FakeWriter(final boolean finishError_){finishError=finishError_;}

		@Override
		public void start(){started=true;}
		@Override
		public long readsWritten(){return readsWritten;}
		@Override
		public long basesWritten(){return basesWritten;}
		@Override
		public synchronized void add(final ArrayList<Read> reads, final long id){
			count(reads);
		}
		@Override
		public void addReads(final ListNum<Read> reads){count(reads.list);}
		@Override
		public void addLines(final ListNum<SamLine> lines){}
		@Override
		public void poison(){}
		@Override
		public boolean waitForFinish(){return poisonAndWait();}
		@Override
		public boolean poisonAndWait(){
			finished=true;
			return finishError;
		}
		@Override
		public void finishError(){finished=true;}
		@Override
		public String fname(){return "fake.out";}
		@Override
		public boolean errorState(){return finished && finishError;}
		@Override
		public boolean finishedSuccessfully(){return finished && !finishError;}

		/** Counts submitted reads and bases like the template read hook. */
		private void count(final ArrayList<Read> reads){
			for(Read r : reads){
				if(r!=null){
					readsWritten+=r.pairCount();
					basesWritten+=r.pairLength();
				}
			}
		}

		/** Controlled finish status. */
		final boolean finishError;
		/** Count of written individual reads. */
		long readsWritten=0;
		/** Count of written bases. */
		long basesWritten=0;
		/** Whether start was called. */
		boolean started=false;
		/** Whether poisonAndWait was called. */
		boolean finished=false;
	}

	/** Tiny valid FASTQ payload used only so constructor paths exist. */
	private static final byte[] FASTQ="@r0\nACGT\n+\n!!!!\n".getBytes(StandardCharsets.UTF_8);
	/** Writer error configured for the next instrumented factory call. */
	private static volatile boolean configuredWriterError;
	/** Reader error configured for the next instrumented factory call. */
	private static volatile boolean configuredReaderError;
	/** Last streamer built by the instrumented factory. */
	private static final AtomicReference<FakeStreamer> lastStreamer=
			new AtomicReference<FakeStreamer>();
	/** Last writer built by the instrumented factory. */
	private static final AtomicReference<FakeWriter> lastWriter=
			new AtomicReference<FakeWriter>();
}
