package test.template;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Timer;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/** Completion-status checks for A_SampleSamStreamer input and output finalization.
 * Prepare mode writes a source-instrumented copy with constructor placeholders
 * satisfied and real finite streams; run mode checks controlled close statuses.
 * @author Jean
 */
public final class TemplateSamStreamerCompletionCase{

	/** Prepares the SamStreamer source copy or runs finite completion-status cases. */
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

	/** Writes a source-instrumented A_SampleSamStreamer copy under out/template. */
	private static void prepare(final Path source, final Path out) throws Exception{
		String text=new String(Files.readAllBytes(
				source.resolve("template/A_SampleSamStreamer.java")), StandardCharsets.UTF_8);
		text=text.replace("A_SampleSamStreamer", "TemplateSamStreamerCompletionSample");
		text=replaceOnce(text, "assert(false) : \"TODO: Default output format might be fasta.\";",
				"// TemplateSamStreamerCompletionCase supplies a concrete SAM output.");
		text=replaceOnce(text, "assert(false) : \"TODO: Check.\";",
				"// TemplateSamStreamerCompletionCase supplies a concrete reference.");
		text=replaceOnce(text, "\tprivate boolean validateParams(){\n"+
				"//\t\tassert(minfoo>0 && minfoo<=maxfoo) : minfoo+\", \"+maxfoo;\n"+
				"\t\tassert(false) : \"TODO\";\n"+
				"\t\treturn true;\n"+
				"\t}",
				"\tprivate boolean validateParams(){\n"+
				"//\t\tassert(minfoo>0 && minfoo<=maxfoo) : minfoo+\", \"+maxfoo;\n"+
				"\t\t// TemplateSamStreamerCompletionCase uses finite constructor settings.\n"+
				"\t\treturn true;\n"+
				"\t}");
		text=replaceOnce(text, "StreamerFactory.makeSamOrBamStreamer",
				"test.template.TemplateSamStreamerCompletionCase.makeSamOrBamStreamer");
		text=replaceOnce(text, "ReadWrite.closeStream(ros)",
				"test.template.TemplateSamStreamerCompletionCase.closeOutput(ros)");

		final Path javaFile=out.resolve("template/TemplateSamStreamerCompletionSample.java");
		Files.createDirectories(javaFile.getParent());
		check(!Files.exists(javaFile), "prepared source already exists: "+javaFile);
		Files.write(javaFile, text.getBytes(StandardCharsets.UTF_8));
		System.out.println("prepared\t"+javaFile);
	}

	/** Runs completion-status cases against a precompiled SamStreamer sample. */
	private static void run(final Path classes, final Path work, final String expect)
			throws Exception{
		check("baseline".equals(expect) || "fixed".equals(expect),
				"expect=baseline or expect=fixed is required");
		check(!Files.exists(work), "Fresh work directory required: "+work);
		Files.createDirectories(work);

		try(ChildFirstLoader loader=new ChildFirstLoader(classes.toUri().toURL(),
				TemplateSamStreamerCompletionCase.class.desiredAssertionStatus())){
			final Class<?> type=Class.forName("template.TemplateSamStreamerCompletionSample",
					true, loader);
			case1(type, work.resolve("normal"), expect, false, false);
			case1(type, work.resolve("input-close-error"), expect, true, false);
			case1(type, work.resolve("output-close-error"), expect, false, true);
		}
	}

	/** Checks one real empty-SAM input/output terminal status pair. */
	private static void case1(final Class<?> type, final Path work, final String expect,
			final boolean inputError, final boolean outputError) throws Exception{
		Files.createDirectories(work);
		configure(inputError, outputError);
		final String name=work.getFileName().toString();
		final Path out=work.resolve(name+".sam");
		CountingStreamer streamer=null;
		Throwable thrown=null;
		int productInputCloseCount=-1;
		try{
			final Object sample=type.getConstructor(String[].class).newInstance((Object)new String[] {
					"in="+writeEmptySam(work.resolve(name+".input.sam")),
					"ref="+writeReference(work.resolve(name+".ref.fa")), "out="+out, "t=1"});
			final Method process=type.getDeclaredMethod("process", Timer.class);
			process.setAccessible(true);
			thrown=invoke(process, sample, new Timer());
			streamer=lastStreamer;
			check(streamer!=null, name+" should create a Streamer");
			productInputCloseCount=streamer.closeCount;
		}finally{
			if(streamer==null){streamer=lastStreamer;}
			if(streamer!=null){streamer.cleanup();}
		}
		checkExpect(name, expect, inputError, outputError, thrown);
		check(productInputCloseCount==expectedInputCloseCount(expect), name+
				" product input close count "+productInputCloseCount);
		check(outputCloseCount==1, name+" should close one output, saw "+outputCloseCount);
		checkBytes(name, EXPECTED_EMPTY_OUTPUT, Files.readAllBytes(out));
		pass(name, thrown);
	}

	/** Redirect target matching StreamerFactory.makeSamOrBamStreamer. */
	public static Streamer makeSamOrBamStreamer(final FileFormat ffin, final int threads,
			final boolean saveHeader, final boolean ordered, final long maxReads,
			final boolean makeReads){
		final CountingStreamer streamer=new CountingStreamer(
				StreamerFactory.makeSamOrBamStreamer(ffin, threads, saveHeader, ordered,
						maxReads, makeReads));
		lastStreamer=streamer;
		return streamer;
	}

	/** Closes the output stream and injects a controlled terminal status. */
	public static boolean closeOutput(final ConcurrentReadOutputStream ros){
		outputCloseCount++;
		return ReadWrite.closeStream(ros) | outputCloseError;
	}

	/** Configures controlled close status for one case. */
	private static void configure(final boolean inputError, final boolean outputError){
		inputCloseError=inputError;
		outputCloseError=outputError;
		outputCloseCount=0;
		lastStreamer=null;
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
			final boolean inputError, final boolean outputError, final Throwable thrown){
		final boolean shouldThrow=outputError || (inputError && "fixed".equals(expect));
		if(shouldThrow){
			check(thrown!=null && thrown.getClass()==RuntimeException.class,
					expect+" should throw RuntimeException for "+name+", got "+thrown);
			check(thrown.getMessage()!=null
					&& thrown.getMessage().contains("terminated in an error state"),
					expect+" should report terminal error state for "+name+", got "+
					thrown.getMessage());
		}else{
			check(thrown==null, expect+" should miss "+name+", got "+thrown);
		}
	}

	/** Expected product-side input close count for this source revision family. */
	private static int expectedInputCloseCount(final String expect){
		return "fixed".equals(expect) ? 1 : 0;
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

	/** Writes a tiny FASTA reference and returns its path. */
	private static Path writeReference(final Path path) throws Exception{
		Files.write(path, FASTA);
		return path;
	}

	/** Writes a valid empty SAM input and returns its path. */
	private static Path writeEmptySam(final Path path) throws Exception{
		Files.write(path, EMPTY_SAM);
		return path;
	}

	/** Verifies exact output bytes. */
	private static void checkBytes(final String name, final byte[] expected,
			final byte[] observed){
		check(observed.length==expected.length,
				name+" output length "+observed.length+" != "+expected.length);
		for(int i=0; i<expected.length; i++){
			check(observed[i]==expected[i], name+" output differs at "+i+": "+
					observed[i]+" != "+expected[i]);
		}
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

	/** Loads the prepared SamStreamer sample before any same-named parent class. */
	private static final class ChildFirstLoader extends URLClassLoader{
		ChildFirstLoader(final URL url, final boolean assertions){
			super(new URL[] {url}, TemplateSamStreamerCompletionCase.class.getClassLoader());
			setDefaultAssertionStatus(assertions);
			setPackageAssertionStatus("template", assertions);
		}

		@Override
		protected Class<?> loadClass(final String name, final boolean resolve)
				throws ClassNotFoundException{
			if(name.equals("template.TemplateSamStreamerCompletionSample")
					|| name.startsWith("template.TemplateSamStreamerCompletionSample$")){
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

	/** Streamer wrapper that records product close calls and injects errorState. */
	private static final class CountingStreamer implements Streamer{
		CountingStreamer(final Streamer delegate_){delegate=delegate_;}

		@Override
		public String fname(){return delegate.fname();}

		@Override
		public void start(){delegate.start();}

		@Override
		public void close(){closeCount++; delegate.close();}

		@Override
		public boolean paired(){return delegate.paired();}

		@Override
		public int pairnum(){return delegate.pairnum();}

		@Override
		public long readsProcessed(){return delegate.readsProcessed();}

		@Override
		public long basesProcessed(){return delegate.basesProcessed();}

		@Override
		public void setSampleRate(final float rate, final long seed){
			delegate.setSampleRate(rate, seed);
		}

		@Override
		public ListNum<Read> nextList(){return delegate.nextList();}

		@Override
		public ListNum<SamLine> nextLines(){return delegate.nextLines();}

		@Override
		public boolean hasMore(){return delegate.hasMore();}

		@Override
		public boolean errorState(){return delegate.errorState() | inputCloseError;}

		@Override
		public void returnList(final ListNum<Read> ln){delegate.returnList(ln);}

		@Override
		public void returnList(final long id, final boolean b){delegate.returnList(id, b);}

		/** Closes any reader still owned by the fixture after product behavior is observed. */
		void cleanup(){
			if(closeCount==0){delegate.close();}
		}

		/** Wrapped real SAM Streamer. */
		private final Streamer delegate;
		/** Product close calls observed before fixture cleanup. */
		int closeCount=0;
	}

	/** Tiny valid SAM with no alignment records. */
	private static final byte[] EMPTY_SAM=("@HD\tVN:1.6\tSO:unsorted\n"+
			"@SQ\tSN:ref\tLN:4\n").getBytes(StandardCharsets.UTF_8);

	/** Tiny reference satisfying the SamStreamer template's scaffold map. */
	private static final byte[] FASTA=">ref\nACGT\n".getBytes(StandardCharsets.UTF_8);

	/** Expected empty SAM from makeCros's useSharedHeader=false output path. */
	private static final byte[] EXPECTED_EMPTY_OUTPUT=("@HD\tVN:1.4\tSO:unsorted\n"+
			"@PG\tID:BBMap\tPN:BBMap\tVN:40.04\n").getBytes(StandardCharsets.UTF_8);

	/** Last constructed Streamer wrapper, or null before makeStreamer runs. */
	private static CountingStreamer lastStreamer;
	/** Controlled input terminal status ORed into CountingStreamer.errorState(). */
	private static boolean inputCloseError;
	/** Controlled output terminal status ORed into closeOutput. */
	private static boolean outputCloseError;
	/** Number of product output close calls observed. */
	private static int outputCloseCount;
}
