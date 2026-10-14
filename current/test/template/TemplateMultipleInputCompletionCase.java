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

import fileIO.ReadWrite;
import shared.Timer;
import stream.ConcurrentReadOutputStream;
import stream.ConcurrentReadStreamInterface;

/** Completion-status checks for A_SampleMultipleInput reader and writer finalization.
 * @author Sucrose
 */
public final class TemplateMultipleInputCompletionCase{

	/** Prepares the source copy or runs finite completion-status cases. */
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

	/** Writes a source-instrumented A_SampleMultipleInput copy under out/template. */
	private static void prepare(final Path source, final Path out) throws Exception{
		String text=new String(Files.readAllBytes(
				source.resolve("template/A_SampleMultipleInput.java")), StandardCharsets.UTF_8);
		text=text.replace("A_SampleMultipleInput", "TemplateMultipleInputCompletionSample");
		text=replaceOnce(text, "ReadWrite.closeStream(ros)",
				"test.template.TemplateMultipleInputCompletionCase.closeOutput(ros)");
		text=replaceOnce(text, "ReadWrite.closeStream(cris)",
				"test.template.TemplateMultipleInputCompletionCase.closeInput(cris)");
		final Path javaFile=out.resolve("template/TemplateMultipleInputCompletionSample.java");
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
				TemplateMultipleInputCompletionCase.class.desiredAssertionStatus())){
			final Class<?> type=Class.forName("template.TemplateMultipleInputCompletionSample",
					true, loader);
			checkOrigin(type, classes);
			case1(type, work.resolve("normal"), expect, false, false);
			case1(type, work.resolve("input-close-error"), expect, true, false);
			case1(type, work.resolve("output-close-error"), expect, false, true);
		}
	}

	/** Checks one controlled input/output terminal status pair. */
	private static void case1(final Class<?> type, final Path work, final String expect,
			final boolean inputError, final boolean outputError) throws Exception{
		Files.createDirectories(work);
		inputCloseError=inputError;
		outputCloseError=outputError;
		inputCloseCount=0;
		outputCloseCount=0;
		final String name=work.getFileName().toString();
		final Path out=work.resolve(name+".out.fq");
		final Object sample=type.getConstructor(String[].class).newInstance((Object)new String[] {
			"in="+writeFastq(work.resolve(name+".fq")), "out="+out, "t=1"});
		final Method process=type.getDeclaredMethod("process", Timer.class);
		process.setAccessible(true);
		final Throwable thrown=invoke(process, sample, new Timer());
		check(inputCloseCount==1, name+" should close one input, saw "+inputCloseCount);
		check(outputCloseCount==1, name+" should close one output, saw "+outputCloseCount);
		checkBytes(name, EXPECTED_OUTPUT, Files.readAllBytes(out));
		checkExpect(name, expect, inputError || outputError, thrown);
		pass(name, thrown);
	}

	/** Closes the input stream and injects a controlled terminal status. */
	public static boolean closeInput(final ConcurrentReadStreamInterface cris){
		inputCloseCount++;
		return ReadWrite.closeStream(cris) | inputCloseError;
	}

	/** Closes the output stream and injects a controlled terminal status. */
	public static boolean closeOutput(final ConcurrentReadOutputStream ros){
		outputCloseCount++;
		return ReadWrite.closeStream(ros) | outputCloseError;
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
			check(thrown.getMessage()!=null &&
					thrown.getMessage().contains("terminated in an error state"),
					"fixed should report terminal error state for "+name+", got "+
					thrown.getMessage());
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

	/** Verifies that the instrumented sample came from the prepared classes root. */
	private static void checkOrigin(final Class<?> clazz, final Path classes) throws Exception{
		final Path origin=Paths.get(clazz.getProtectionDomain().getCodeSource().
				getLocation().toURI()).toAbsolutePath().normalize();
		check(origin.equals(classes.toAbsolutePath().normalize()),
				"Wrong sample origin: "+origin+" != "+classes);
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

	/** Writes a tiny valid FASTQ input and returns its path. */
	private static Path writeFastq(final Path path) throws Exception{
		Files.write(path, FASTQ);
		return path;
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

	/** Loads the prepared MultipleInput sample before any same-named parent class. */
	private static final class ChildFirstLoader extends URLClassLoader{
		ChildFirstLoader(final URL url, final boolean assertions){
			super(new URL[] {url}, TemplateMultipleInputCompletionCase.class.getClassLoader());
			setDefaultAssertionStatus(assertions);
			setPackageAssertionStatus("template", assertions);
		}

		@Override
		protected Class<?> loadClass(final String name, final boolean resolve)
				throws ClassNotFoundException{
			if(name.equals("template.TemplateMultipleInputCompletionSample")
					|| name.startsWith("template.TemplateMultipleInputCompletionSample$")){
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

	/** Tiny valid FASTQ payload used only so constructor paths exist. */
	private static final byte[] FASTQ="@r0\nACGT\n+\nIIII\n".getBytes(StandardCharsets.UTF_8);
	/** Exact output from the template's one real FASTQ read. */
	private static final byte[] EXPECTED_OUTPUT=
			"@r0\nACGT\n+\nIIII\n".getBytes(StandardCharsets.UTF_8);
	/** Error status injected after real input closure. */
	private static volatile boolean inputCloseError;
	/** Error status injected after real output closure. */
	private static volatile boolean outputCloseError;
	/** Count of input wrapper calls for the current case. */
	private static volatile int inputCloseCount;
	/** Count of output wrapper calls for the current case. */
	private static volatile int outputCloseCount;
}
