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
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/** Source-derived check for StreamerMT root filtering through real writers.
 * Prepare mode replaces only the documented read-pair hook in an owned source
 * copy; run mode loads the guard-compiled copy and checks that a rejected
 * nonnull root is removed from FASTQ/SAM output.
 * @author Jean
 */
public final class TemplateStreamerFilterCase{

	/** Prepares patched source or runs a finite FASTQ/SAM filter case. */
	public static void main(String[] args) throws Exception{
		final HashMap<String, String> map=parseArgs(args);
		final String mode=map.get("mode");
		if("prepare".equals(mode)){
			prepare(Paths.get(required(map, "source")), Paths.get(required(map, "out")));
		}else if("run".equals(mode)){
			run(Paths.get(required(map, "classes")), Paths.get(required(map, "work")),
					required(map, "format"), required(map, "scenario"),
					required(map, "expect"));
		}else{
			throw new IllegalArgumentException("Expected mode=prepare or mode=run");
		}
	}

	/** Writes the patched sample source under out/template/A_SampleStreamerMT.java. */
	private static void prepare(final Path source, final Path out) throws Exception{
		final String original=new String(Files.readAllBytes(source), StandardCharsets.UTF_8);
		String patched=original;
		int occurrences=0;
		for(String[] pair : PROCESS_READ_PAIR_REPLACEMENTS){
			final int count=countOccurrences(original, pair[0]);
			if(count>0){patched=original.replace(pair[0], pair[1]);}
			occurrences+=count;
		}
		check(occurrences==1, "processReadPair anchor occurrence");

		final Path sourceFile=out.resolve("template/A_SampleStreamerMT.java");
		Files.createDirectories(sourceFile.getParent());
		check(!Files.exists(sourceFile), "prepared source already exists: "+sourceFile);
		Files.write(sourceFile, patched.getBytes(StandardCharsets.UTF_8));
		System.out.println("prepared\t"+sourceFile);
	}

	/** Runs one finite main-program case and checks the exact output records. */
	private static void run(final Path classes, final Path work, final String format,
			final String scenario, final String expect) throws Exception{
		check("fastq".equals(format) || "sam".equals(format),
				"Expected format=fastq or format=sam");
		check("allkeep".equals(scenario) || "drop".equals(scenario),
				"Expected scenario=allkeep or scenario=drop");
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");

		check(!Files.exists(work), "work directory already exists: "+work);
		Files.createDirectories(work);
		final Path in=work.resolve("in."+format);
		final Path out=work.resolve("out."+format);
		writeInput(in, format, inputNames(scenario));

		try(ChildFirstLoader loader=new ChildFirstLoader(classes.toUri().toURL(),
				TemplateStreamerFilterCase.class.desiredAssertionStatus())){
			final Class<?> sample=Class.forName("template.A_SampleStreamerMT", true, loader);
			final Method main=sample.getMethod("main", String[].class);
			final String[] args=new String[] {"in="+in, "out="+out, "ow=t", "t=1",
					"ordered=t", "qin=33", "qout=33"};
			try{
				main.invoke(null, (Object)args);
			}catch(InvocationTargetException ite){
				final Throwable cause=ite.getCause();
				if(cause instanceof Exception){throw (Exception)cause;}
				if(cause instanceof Error){throw (Error)cause;}
				throw new RuntimeException(cause);
			}
		}

		final List<String> actual=outputRecords(out, format);
		final List<String> expected=expectedRecords(format, scenario, expect);
		check(actual.equals(expected), "output records: "+actual+" != "+expected);
		System.out.println("PASS\t"+format+"\t"+scenario+"\t"+expect+"\t"+actual);
	}

	/** Returns input names for the ordinary control or filtering case. */
	private static List<String> inputNames(final String scenario){
		if("allkeep".equals(scenario)){return Arrays.asList("keep1", "keep2");}
		return Arrays.asList("keep", "drop");
	}

	/** Returns the exact FASTQ records or SAM @SQ/body lines expected. */
	private static List<String> expectedRecords(final String format, final String scenario,
			final String expect){
		final List<String> names=("allkeep".equals(scenario) || "baseline".equals(expect)) ?
				inputNames(scenario) : Arrays.asList("keep");
		return records(format, names);
	}

	/** Writes a two-record FASTQ or SAM input. */
	private static void writeInput(final Path in, final String format, final List<String> names)
			throws Exception{
		final List<String> lines=records(format, names);
		Files.write(in, String.join("\n", lines).concat("\n").getBytes(StandardCharsets.UTF_8));
	}

	/** Returns literal FASTQ records or a SAM @SQ row plus literal body rows. */
	private static List<String> records(final String format, final List<String> names){
		final ArrayList<String> lines=new ArrayList<String>();
		if("sam".equals(format)){
			lines.add("@SQ\tSN:ref\tLN:100");
			for(int i=0; i<names.size(); i++){
				lines.add(names.get(i)+"\t0\tref\t"+(i+1)
						+"\t255\t4M\t*\t0\t0\t"+BASES[i]+"\tIIII\tMD:Z:4");
			}
		}else{
			for(int i=0; i<names.size(); i++){
				lines.add("@"+names.get(i));
				lines.add(BASES[i]);
				lines.add("+");
				lines.add("IIII");
			}
		}
		return lines;
	}

	/** Extracts full FASTQ records or SAM @SQ/body lines. */
	private static ArrayList<String> outputRecords(final Path out, final String format)
			throws Exception{
		final ArrayList<String> lines=new ArrayList<String>();
		for(String line : Files.readAllLines(out, StandardCharsets.UTF_8)){
			if("fastq".equals(format) || line.startsWith("@SQ\t")
					|| (line.length()>0 && line.charAt(0)!='@')){
				lines.add(line);
			}
		}
		return lines;
	}

	/** Parses key=value arguments while preserving later '=' characters. */
	private static HashMap<String, String> parseArgs(final String[] args){
		final HashMap<String, String> map=new HashMap<String, String>();
		for(String arg : args){
			final String[] pair=arg.split("=", 2);
			check(pair.length==2, "Arguments require key=value");
			map.put(pair[0], pair[1]);
		}
		return map;
	}

	/** Fetches a required key. */
	private static String required(final HashMap<String, String> map, final String key){
		final String value=map.get(key);
		check(value!=null && value.length()>0, "Missing "+key);
		return value;
	}

	/** Counts literal substring occurrences. */
	private static int countOccurrences(final String text, final String query){
		int count=0;
		for(int index=text.indexOf(query); index>=0;
				index=text.indexOf(query, index+query.length())){
			count++;
		}
		return count;
	}

	/** Remains enabled under both -ea and -da. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Loads the patched sample before any same-named class from the parent loader. */
	private static final class ChildFirstLoader extends URLClassLoader{
		ChildFirstLoader(final URL url, final boolean assertions){
			super(new URL[] {url}, TemplateStreamerFilterCase.class.getClassLoader());
			setDefaultAssertionStatus(assertions);
			setPackageAssertionStatus("template", assertions);
		}

		@Override
		protected Class<?> loadClass(final String name, final boolean resolve)
				throws ClassNotFoundException{
			if(name.startsWith("template.A_SampleStreamerMT")){
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

	private static final String[][] PROCESS_READ_PAIR_REPLACEMENTS={{
		"\t\tprivate boolean processReadPair(Read r1, Read r2) {\n"
		+"\t\t\treturn true;\n"
		+"\t\t}\n",
		"\t\tprivate boolean processReadPair(Read r1, Read r2) {\n"
		+"\t\t\t// Simulate a documented custom filter hook for a concrete tool.\n"
		+"\t\t\treturn !\"drop\".equals(r1.id);\n"
		+"\t\t}\n"}, {
		"\t\tprivate boolean processReadPair(Read r1, Read r2){\n"
		+"\t\t\treturn true;\n"
		+"\t\t}\n",
		"\t\tprivate boolean processReadPair(Read r1, Read r2){\n"
		+"\t\t\t// Simulate a documented custom filter hook for a concrete tool.\n"
		+"\t\t\treturn !\"drop\".equals(r1.id);\n"
		+"\t\t}\n"}
	};

	private static final String[] BASES={"ACGT", "TGCA"};
}
