package test.template;

import java.lang.reflect.Field;
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
import parse.Parser;
import shared.Shared;
import shared.Timer;
import stream.ConcurrentReadOutputStream;
import stream.ConcurrentReadStreamInterface;

/** Completion-status and output-policy checks for A_Sample2 and A_SampleD.
 * @author Jean
 */
public final class TemplateSample2DCompletionPolicyCase{

	/** Prepares instrumented source copies or runs finite policy cases. */
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

	/** Writes source-instrumented A_Sample2 and A_SampleD copies under out/template. */
	private static void prepare(final Path source, final Path out) throws Exception{
		String sample2=new String(Files.readAllBytes(source.resolve("template/A_Sample2.java")),
				StandardCharsets.UTF_8);
		sample2=sample2.replace("A_Sample2", "TemplateSample2CompletionPolicySample");
		sample2=replaceOnce(sample2, "ReadWrite.closeStreams(cris, ros)",
				"test.template.TemplateSample2DCompletionPolicyCase.closeSample2(cris, ros)");
		writePrepared(out, "TemplateSample2CompletionPolicySample.java", sample2);

		String sampleD=new String(Files.readAllBytes(source.resolve("template/A_SampleD.java")),
				StandardCharsets.UTF_8);
		sampleD=sampleD.replace("A_SampleD", "TemplateSampleDCompletionPolicySample");
		sampleD=replaceOnce(sampleD, "ReadWrite.closeStreams(crisD, ros)",
				"test.template.TemplateSample2DCompletionPolicyCase.closeSampleD(crisD, ros)");
		writePrepared(out, "TemplateSampleDCompletionPolicySample.java", sampleD);
	}

	/** Writes one prepared sample source. */
	private static void writePrepared(final Path out, final String name, final String text)
			throws Exception{
		final Path javaFile=out.resolve("template").resolve(name);
		Files.createDirectories(javaFile.getParent());
		check(!Files.exists(javaFile), "prepared source already exists: "+javaFile);
		Files.write(javaFile, text.getBytes(StandardCharsets.UTF_8));
		System.out.println("prepared\t"+javaFile);
	}

	/** Runs descriptor and close-status cases against precompiled prepared samples. */
	private static void run(final Path classes, final Path work, final String expect)
			throws Exception{
		check("baseline".equals(expect) || "fixed".equals(expect),
				"expect=baseline or expect=fixed is required");
		check(!Files.exists(work), "Fresh work directory required: "+work);
		Files.createDirectories(work);

		try(ChildFirstLoader loader=new ChildFirstLoader(classes.toUri().toURL(),
				TemplateSample2DCompletionPolicyCase.class.desiredAssertionStatus())){
			final Class<?> sample2=Class.forName(
					"template.TemplateSample2CompletionPolicySample", true, loader);
			final Class<?> sampleD=Class.forName(
					"template.TemplateSampleDCompletionPolicySample", true, loader);
			checkOrigin(sample2, classes);
			checkOrigin(sampleD, classes);

			descriptorCases(sample2, work.resolve("sample2-descriptor"), expect);
			descriptorCases(sampleD, work.resolve("sampleD-descriptor"), expect);
			processCase(sample2, work.resolve("sample2-normal"), expect, false);
			processCase(sample2, work.resolve("sample2-close-error"), expect, true);
			processCase(sampleD, work.resolve("sampleD-normal"), expect, false);
			processCase(sampleD, work.resolve("sampleD-close-error"), expect, true);
		}
	}

	/** Runs the constructor-only output-policy rows for one prepared sample. */
	private static void descriptorCases(final Class<?> type, final Path work,
			final String expect) throws Exception{
		Files.createDirectories(work);
		normalOutput(type, work.resolve("normal"), expect);
		overwriteFalse(type, work.resolve("overwrite-false"), expect);
		appendTrue(type, work.resolve("append-true"), expect);
	}

	/** Verifies the ordinary constructor descriptor policy. */
	private static void normalOutput(final Class<?> type, final Path work,
			final String expect) throws Exception{
		Files.createDirectories(work);
		final Path in=writeFastq(work.resolve("normal.fq"));
		final Path out=work.resolve("normal.out.fq");
		final Constructed observed=construct(type, "in="+in, "out="+out);
		check(observed.thrown==null, type.getSimpleName()+" normal threw "+observed.thrown);
		checkFlags(observed.output, true, false, type.getSimpleName()+" normal");
		check(!Files.exists(out), "constructor-only normal check must not write output");
		pass(expect, type, "descriptor-normal", flags(observed.output));
	}

	/** Captures descriptor overwrite policy separately from actual output writes. */
	private static void overwriteFalse(final Class<?> type, final Path work,
			final String expect) throws Exception{
		Files.createDirectories(work);
		final Path in=writeFastq(work.resolve("ow-input.fq"));
		final Path out=work.resolve("existing.fq");
		Files.write(out, SENTINEL);
		final Constructed observed=construct(type, "in="+in, "out="+out, "ow=f");
		check(observed.thrown==null, type.getSimpleName()+" ow=f threw "+observed.thrown);
		checkFlags(parsedOutput(out, "ow=f"), false, false,
				"FileFormat-normalized Parser overwrite=false");
		if("baseline".equals(expect)){
			checkFlags(observed.output, true, false, type.getSimpleName()+" baseline ow=f");
		}else{
			checkFlags(observed.output, false, false, type.getSimpleName()+" fixed ow=f");
		}
		checkBytes("existing output sentinel", SENTINEL, Files.readAllBytes(out));
		pass(expect, type, "descriptor-overwrite-false", flags(observed.output));
	}

	/** Captures descriptor append policy separately from actual output writes. */
	private static void appendTrue(final Class<?> type, final Path work, final String expect)
			throws Exception{
		Files.createDirectories(work);
		final Path in=writeFastq(work.resolve("append-input.fq"));
		final Path out=work.resolve("append.fq");
		final Constructed observed=construct(type, "in="+in, "out="+out, "append=t");
		check(observed.thrown==null, type.getSimpleName()+" append=t threw "+observed.thrown);
		checkFlags(parsedOutput(out, "append=t"), false, true,
				"FileFormat-normalized Parser append=true");
		if("baseline".equals(expect)){
			checkFlags(observed.output, true, false,
					type.getSimpleName()+" baseline append=t");
		}else{
			checkFlags(observed.output, false, true, type.getSimpleName()+" fixed append=t");
		}
		check(!Files.exists(out), "constructor-only append check must not write output");
		pass(expect, type, "descriptor-append-true", flags(observed.output));
	}

	/** Runs one finite process row through the prepared sample. */
	private static void processCase(final Class<?> type, final Path work,
			final String expect, final boolean error) throws Exception{
		Files.createDirectories(work);
		closeError=error;
		closeSample2Count=0;
		closeSampleDCount=0;

		final String name=work.getFileName().toString();
		final Path out=work.resolve(name+".out.fq");
		final boolean oldUseMpi=Shared.USE_MPI;
		Throwable thrown;
		try{
			Shared.USE_MPI=false;
			final Object sample=type.getConstructor(String[].class).newInstance((Object)
					new String[] {"in="+writeFastq(work.resolve(name+".fq")),
					"out="+out, "t=1"});
			final Method process=type.getDeclaredMethod("process", Timer.class);
			process.setAccessible(true);
			thrown=invoke(process, sample, new Timer());
		}finally{
			Shared.USE_MPI=oldUseMpi;
		}
		final int expectedSample2=sample2(type) ? 1 : 0;
		final int expectedSampleD=sample2(type) ? 0 : 1;
		check(closeSample2Count==expectedSample2, name+" Sample2 close count "+
				closeSample2Count+" != "+expectedSample2);
		check(closeSampleDCount==expectedSampleD, name+" SampleD close count "+
				closeSampleDCount+" != "+expectedSampleD);
		checkBytes(name, EXPECTED_OUTPUT, Files.readAllBytes(out));
		checkExpect(name, expect, error, thrown);
		pass(expect, type, error ? "process-close-error" : "process-normal",
				thrownName(thrown));
	}

	/** Closes A_Sample2's streams and injects a controlled combined status. */
	public static boolean closeSample2(final ConcurrentReadStreamInterface cris,
			final ConcurrentReadOutputStream ros){
		closeSample2Count++;
		return ReadWrite.closeStreams(cris, ros) | closeError;
	}

	/** Closes A_SampleD's streams and injects a controlled combined status. */
	public static boolean closeSampleD(final ConcurrentReadStreamInterface crisD,
			final ConcurrentReadOutputStream ros){
		closeSampleDCount++;
		return ReadWrite.closeStreams(crisD, ros) | closeError;
	}

	/** Constructs a prepared sample and returns its private output descriptor. */
	private static Constructed construct(final Class<?> type, final String... args)
			throws Exception{
		try{
			final Object sample=type.getConstructor(String[].class).newInstance(
					(Object)args);
			return new Constructed(output(type, sample), null);
		}catch(InvocationTargetException e){
			return new Constructed(null, e.getCause());
		}
	}

	/** Builds an output descriptor from a parsed ordinary output-policy argument. */
	private static FileFormat parsedOutput(final Path out, final String arg){
		final Parser parser=new Parser();
		final String[] split=arg.split("=", 2);
		check(split.length==2 && parser.parse(arg, split[0], split[1]),
				"Parser should accept "+arg);
		return FileFormat.testOutput(out.toString(), FileFormat.FASTQ, null, true,
				parser.overwrite, parser.append, false);
	}

	/** Returns the private output descriptor. */
	private static FileFormat output(final Class<?> type, final Object sample)
			throws Exception{
		final Field field=type.getDeclaredField("ffout1");
		field.setAccessible(true);
		return (FileFormat)field.get(sample);
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

	/** Verifies baseline/fixed thrown state for one close-status row. */
	private static void checkExpect(final String name, final String expect,
			final boolean error, final Throwable thrown){
		if("baseline".equals(expect) || !error){
			check(thrown==null, expect+" should miss "+name+", got "+thrown);
		}else{
			check(thrown!=null && thrown.getClass()==RuntimeException.class,
					"fixed should throw RuntimeException for "+name+", got "+thrown);
			check(thrown.getMessage()!=null &&
					thrown.getMessage().contains("terminated in an error state"),
					"fixed should report terminal error state for "+name+", got "+
					thrown.getMessage());
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

	/** Verifies that the instrumented sample came from the prepared classes root. */
	private static void checkOrigin(final Class<?> clazz, final Path classes) throws Exception{
		final Path origin=Paths.get(clazz.getProtectionDomain().getCodeSource().
				getLocation().toURI()).toAbsolutePath().normalize();
		check(origin.equals(classes.toAbsolutePath().normalize()),
				"Wrong sample origin: "+origin+" != "+classes);
	}

	/** Writes a tiny valid FASTQ input and returns its path. */
	private static Path writeFastq(final Path path) throws Exception{
		Files.write(path, FASTQ);
		return path;
	}

	/** Checks descriptor flags. */
	private static void checkFlags(final FileFormat ff, final boolean overwrite,
			final boolean append, final String name){
		check(ff.overwrite()==overwrite, name+" overwrite flag "+flags(ff));
		check(ff.append()==append, name+" append flag "+flags(ff));
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

	/** Formats a descriptor's policy flags. */
	private static String flags(final FileFormat ff){
		return "ow="+ff.overwrite()+",append="+ff.append();
	}

	/** Formats a captured throwable type. */
	private static String thrownName(final Throwable thrown){
		return thrown==null ? "none" : thrown.getClass().getName();
	}

	/** True for the prepared A_Sample2 copy. */
	private static boolean sample2(final Class<?> type){
		return type.getName().endsWith("TemplateSample2CompletionPolicySample");
	}

	/** Prints one successful case row. */
	private static void pass(final String expect, final Class<?> type,
			final String name, final String observed){
		System.out.println("PASS\t"+expect+"\t"+type.getSimpleName()+"\t"+
				name+"\t"+observed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Loads the prepared samples before any same-named parent classes. */
	private static final class ChildFirstLoader extends URLClassLoader{
		ChildFirstLoader(final URL url, final boolean assertions){
			super(new URL[] {url}, TemplateSample2DCompletionPolicyCase.class.getClassLoader());
			setDefaultAssertionStatus(assertions);
			setPackageAssertionStatus("template", assertions);
		}

		@Override
		protected Class<?> loadClass(final String name, final boolean resolve)
				throws ClassNotFoundException{
			if(name.equals("template.TemplateSample2CompletionPolicySample") ||
					name.equals("template.TemplateSampleDCompletionPolicySample")){
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

	/** Constructor result or thrown exception. */
	private static final class Constructed{

		/** Stores constructor observations. */
		Constructed(final FileFormat output_, final Throwable thrown_){
			output=output_;
			thrown=thrown_;
		}

		/** Output descriptor, or null when construction failed or output is disabled. */
		final FileFormat output;
		/** Constructor throwable, or null on success. */
		final Throwable thrown;
	}

	/** Tiny valid FASTQ payload used for real finite input. */
	private static final byte[] FASTQ="@r0\nACGT\n+\nIIII\n".getBytes(StandardCharsets.UTF_8);
	/** Exact output from the templates' one real FASTQ read. */
	private static final byte[] EXPECTED_OUTPUT=
			"@r0\nACGT\n+\nIIII\n".getBytes(StandardCharsets.UTF_8);
	/** Bytes used to prove constructor-only checks do not write output. */
	private static final byte[] SENTINEL="sentinel\n".getBytes(StandardCharsets.UTF_8);
	/** Error status injected after real combined stream closure. */
	private static volatile boolean closeError;
	/** Count of A_Sample2 close wrapper calls for the current case. */
	private static volatile int closeSample2Count;
	/** Count of A_SampleD close wrapper calls for the current case. */
	private static volatile int closeSampleDCount;
}
