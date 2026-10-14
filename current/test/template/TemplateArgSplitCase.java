package test.template;

import java.lang.reflect.Constructor;
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

/** Source-instrumented checks for template-local key=value decoding loops.
 * Prepare mode adds local decode probes to owned source copies; run mode checks
 * the compiled copies preserve values containing '=' without changing existing
 * empty and null handling.
 * @author Jean
 */
public final class TemplateArgSplitCase{

	/** Prepares patched sources or checks the compiled probes. */
	public static void main(String[] args) throws Exception{
		final HashMap<String, String> map=parseArgs(args);
		final String mode=map.get("mode");
		final String group=map.containsKey("group") ? map.get("group") : "default";
		if("prepare".equals(mode)){
			prepare(Paths.get(required(map, "source")), Paths.get(required(map, "out")),
					group);
		}else if("run".equals(mode)){
			run(Paths.get(required(map, "classes")), Paths.get(required(map, "work")),
					required(map, "expect"), group);
		}else{
			throw new IllegalArgumentException("Expected mode=prepare or mode=run");
		}
	}

	/** Writes patched concrete samples and generated BBTool_ST probes under out/template. */
	private static void prepare(final Path source, final Path out, final String group)
			throws Exception{
		final String[][] loops=concreteLoops(group);
		final Path targetDir=out.resolve("template");
		Files.createDirectories(targetDir);
		for(String[] loop : loops){
			final String file=loop[0];
			final String loopName=loop[1];
			final Path sourceFile=source.resolve(file);
			final Path targetFile=targetDir.resolve(file);
			check(!Files.exists(targetFile), "prepared source already exists: "+targetFile);
			final String original=new String(Files.readAllBytes(sourceFile), StandardCharsets.UTF_8);
			final String patched=insertProbe(original, loopName, file);
			Files.write(targetFile, patched.getBytes(StandardCharsets.UTF_8));
			System.out.println("prepared\t"+loopName+"\t"+targetFile);
		}
		writeGenerated(targetDir.resolve("TemplateArgSplitProbe.java"), PROBE_SOURCE);
		if(hasBBToolLoops(group)){
			writeGenerated(targetDir.resolve("TemplateArgSplitBBToolProbe.java"),
					BBTOOL_PROBE_SOURCE);
		}
	}

	/** Inserts one source-instrumented observation branch after local normalization. */
	private static String insertProbe(final String original, final String loopName,
			final String file){
		final StringBuilder sb=new StringBuilder(original.length()+200);
		int inserted=0;
		int start=0;
		int index=original.indexOf(NULL_NORMALIZATION);
		while(index>=0){
			inserted++;
			sb.append(original, start, index+NULL_NORMALIZATION.length());
			final int lineStart=original.lastIndexOf('\n', index)+1;
			final String indent=original.substring(lineStart, index);
			sb.append('\n').append(indent).append("if(a.equals(\"jean_probe\")){"
					+"throw new TemplateArgSplitProbe(\"").append(loopName)
					.append("\", arg, a, b);}");
			start=index+NULL_NORMALIZATION.length();
			index=original.indexOf(NULL_NORMALIZATION, start);
		}
		check(inserted==1, file+" probe insertion count "+inserted);
		sb.append(original, start, original.length());
		return sb.toString();
	}

	/** Writes one generated helper source file. */
	private static void writeGenerated(final Path path, final String text) throws Exception{
		check(!Files.exists(path), "prepared source already exists: "+path);
		Files.write(path, text.getBytes(StandardCharsets.UTF_8));
		System.out.println("prepared\tgenerated\t"+path);
	}

	/** Runs all local-decode observations against precompiled classes. */
	private static void run(final Path classes, final Path work, final String expect,
			final String group) throws Exception{
		final String[][] loops=concreteLoops(group);
		check("baseline".equals(expect) || "fixed".equals(expect),
				"Expected expect=baseline or expect=fixed");
		check(!Files.exists(work), "work directory already exists: "+work);
		Files.createDirectories(work);
		try(ChildFirstLoader loader=new ChildFirstLoader(classes.toUri().toURL(),
				TemplateArgSplitCase.class.desiredAssertionStatus())){
			for(String[] loop : loops){
				for(String[] valueCase : VALUE_CASES){
					checkConcrete(loader, loop[1], valueCase, expect);
				}
			}
			if(hasBBToolLoops(group)){
				for(String[] valueCase : VALUE_CASES){
					checkBBToolConstructor(loader, valueCase, expect);
					checkBBToolReparse(loader, work, valueCase, expect);
				}
			}
		}
	}

	/** Checks one constructor-local or private-parse loop. */
	private static void checkConcrete(final ClassLoader loader, final String name,
			final String[] valueCase, final String expect) throws Exception{
		final Class<?> sample=Class.forName("template."+name, true, loader);
		final Constructor<?> constructor=sample.getConstructor(String[].class);
		final TemplateArgSplitObservation obs=invokeProbe(constructor, new String[] {arg(valueCase)});
		checkObservation(obs, name, valueCase, expect);
	}

	/** Checks BBTool_ST's constructor loop through a generated subclass. */
	private static void checkBBToolConstructor(final ClassLoader loader, final String[] valueCase,
			final String expect) throws Exception{
		final Class<?> sample=Class.forName("template.TemplateArgSplitBBToolProbe", true, loader);
		final Constructor<?> constructor=sample.getConstructor(String[].class);
		final TemplateArgSplitObservation obs=invokeProbe(constructor, new String[] {arg(valueCase)});
		checkObservation(obs, "BBTool_ST.constructor", valueCase, expect);
	}

	/** Checks BBTool_ST's protected reparse loop through a generated subclass. */
	private static void checkBBToolReparse(final ClassLoader loader, final Path work,
			final String[] valueCase, final String expect) throws Exception{
		final Path in=work.resolve("bbtool.fq");
		Files.write(in, "@r\nACGT\n+\nIIII\n".getBytes(StandardCharsets.UTF_8));
		final Class<?> sample=Class.forName("template.TemplateArgSplitBBToolProbe", true, loader);
		final Constructor<?> constructor=sample.getConstructor(String[].class);
		final Object instance=constructor.newInstance((Object)new String[] {
				"in="+in, "out="+work.resolve("bbtool.out"), "ow=t"});
		final Method reparse=sample.getMethod("reparsePublic", String[].class);
		final TemplateArgSplitObservation obs=invokeProbe(reparse, instance,
				new String[] {arg(valueCase)});
		checkObservation(obs, "BBTool_ST.reparse", valueCase, expect);
	}

	/** Builds one raw probe argument. */
	private static String arg(final String[] valueCase){
		return "jean_probe="+valueCase[1];
	}

	/** Returns the loop-local probe argument after existing preprocessing. */
	private static String localArg(final String loop, final String[] valueCase){
		if(!"BBTool_ST.reparse".equals(loop) && "null".equals(valueCase[0])){
			return "jean_probe=";
		}
		return arg(valueCase);
	}

	/** Invokes a constructor expected to throw one generated probe exception. */
	private static TemplateArgSplitObservation invokeProbe(final Constructor<?> constructor,
			final String[] args) throws Exception{
		try{
			constructor.newInstance((Object)args);
		}catch(InvocationTargetException ite){
			return observation(ite.getCause());
		}
		throw new IllegalStateException("Expected TemplateArgSplitProbe from constructor");
	}

	/** Invokes a method expected to throw one generated probe exception. */
	private static TemplateArgSplitObservation invokeProbe(final Method method,
			final Object instance, final String[] args) throws Exception{
		try{
			method.invoke(instance, (Object)args);
		}catch(InvocationTargetException ite){
			return observation(ite.getCause());
		}
		throw new IllegalStateException("Expected TemplateArgSplitProbe from method");
	}

	/** Converts a generated probe exception into parent-loader data. */
	private static TemplateArgSplitObservation observation(final Throwable cause)
			throws Exception{
		check(cause!=null && cause.getClass().getName().equals("template.TemplateArgSplitProbe"),
				"Unexpected exception: "+cause);
		final Class<?> clazz=cause.getClass();
		return new TemplateArgSplitObservation(
				stringField(clazz, cause, "loop"),
				stringField(clazz, cause, "arg"),
				stringField(clazz, cause, "a"),
				stringField(clazz, cause, "b"));
	}

	/** Reads a String field from a generated probe exception. */
	private static String stringField(final Class<?> clazz, final Object target,
			final String name) throws Exception{
		final Field field=clazz.getDeclaredField(name);
		field.setAccessible(true);
		return (String)field.get(target);
	}

	/** Checks one observation against the baseline or repaired expectation. */
	private static void checkObservation(final TemplateArgSplitObservation obs,
			final String loop, final String[] valueCase, final String expect){
		final String expectedB=expectedB(loop, valueCase, expect);
		check(loop.equals(obs.loop), "loop", obs.loop, loop);
		check(localArg(loop, valueCase).equals(obs.arg), "arg", obs.arg,
				localArg(loop, valueCase));
		check("jean_probe".equals(obs.a), "a", obs.a, "jean_probe");
		check(equals(expectedB, obs.b), "b", obs.b, expectedB);
		System.out.println("PASS\t"+expect+"\t"+loop+"\t"+valueCase[0]
				+"\t"+(obs.b==null)+"\t"+String.valueOf(obs.b));
	}

	/** Returns the expected right hand side for one case. */
	private static String expectedB(final String loop, final String[] valueCase,
			final String expect){
		if("equals".equals(valueCase[0]) && "baseline".equals(expect)){return "left";}
		if("BBTool_ST.reparse".equals(loop) && "null".equals(valueCase[0])){
			return "null";
		}
		return valueCase[2];
	}

	/** Null-safe String equality. */
	private static boolean equals(final String a, final String b){
		return a==null ? b==null : a.equals(b);
	}

	/** Selects the default repaired loops or the remaining concrete-template loops. */
	private static String[][] concreteLoops(final String group){
		if("default".equals(group)){return CONCRETE_LOOPS;}
		if("remaining".equals(group)){return REMAINING_LOOPS;}
		throw new IllegalArgumentException("Expected group=default or group=remaining");
	}

	/** Returns true when the selected group should also probe BBTool_ST. */
	private static boolean hasBBToolLoops(final String group){
		return "default".equals(group);
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

	/** Remains enabled under both -ea and -da. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Remains enabled under both -ea and -da. */
	private static void check(final boolean condition, final String name, final String actual,
			final String expected){
		if(!condition){throw new IllegalStateException(
				name+": "+String.valueOf(actual)+" != "+String.valueOf(expected));}
	}

	/** Loads prepared sample classes before any same-named class from the parent loader. */
	private static final class ChildFirstLoader extends URLClassLoader{
		ChildFirstLoader(final URL url, final boolean assertions){
			super(new URL[] {url}, TemplateArgSplitCase.class.getClassLoader());
			setDefaultAssertionStatus(assertions);
			setPackageAssertionStatus("template", assertions);
		}

		@Override
		protected Class<?> loadClass(final String name, final boolean resolve)
				throws ClassNotFoundException{
			if(name.startsWith("template.A_Sample")
					|| name.startsWith("template.TemplateArgSplitBBToolProbe")
					|| name.startsWith("template.TemplateArgSplitProbe")){
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

	/** Parent-loader copy of one local-decode observation. */
	private static final class TemplateArgSplitObservation{
		TemplateArgSplitObservation(final String loop_, final String arg_, final String a_,
				final String b_){
			loop=loop_;
			arg=arg_;
			a=a_;
			b=b_;
		}

		final String loop;
		final String arg;
		final String a;
		final String b;
	}

	private static final String NULL_NORMALIZATION=
			"if(b!=null && b.equalsIgnoreCase(\"null\")){b=null;}";

	private static final String[][] CONCRETE_LOOPS={{
		"A_SampleBasic.java", "A_SampleBasic"
	}, {
		"A_SampleBasic2.java", "A_SampleBasic2"
	}, {
		"A_SampleD.java", "A_SampleD"
	}, {
		"A_Sample2.java", "A_Sample2"
	}, {
		"A_SampleSummary.java", "A_SampleSummary"
	}, {
		"A_Sample_Generator.java", "A_Sample_Generator"
	}};

	private static final String[][] REMAINING_LOOPS={{
		"A_SampleByteFile.java", "A_SampleByteFile"
	}, {
		"A_SampleByteFileMT.java", "A_SampleByteFileMT"
	}, {
		"A_SampleMT.java", "A_SampleMT"
	}, {
		"A_SampleMultipleInput.java", "A_SampleMultipleInput"
	}, {
		"A_Sample_Unpaired.java", "A_Sample_Unpaired"
	}, {
		"A_Sample.java", "A_Sample"
	}, {
		"A_SampleStreamerMT.java", "A_SampleStreamerMT"
	}, {
		"A_SampleStreamer.java", "A_SampleStreamer"
	}, {
		"A_SampleSamStreamer.java", "A_SampleSamStreamer"
	}};

	private static final String[][] VALUE_CASES={{
		"equals", "left=right", "left=right"
	}, {
		"ordinary", "ordinary", "ordinary"
	}, {
		"empty", "", null
	}, {
		"null", "null", null
	}};

	private static final String PROBE_SOURCE=
			"package template;\n"
			+"\n"
			+"final class TemplateArgSplitProbe extends RuntimeException{\n"
			+"\tTemplateArgSplitProbe(final String loop_, final String arg_, "
			+"final String a_, final String b_){\n"
			+"\t\tloop=loop_;\n"
			+"\t\targ=arg_;\n"
			+"\t\ta=a_;\n"
			+"\t\tb=b_;\n"
			+"\t}\n"
			+"\n"
			+"\tfinal String loop;\n"
			+"\tfinal String arg;\n"
			+"\tfinal String a;\n"
			+"\tfinal String b;\n"
			+"}\n";

	private static final String BBTOOL_PROBE_SOURCE=
			"package template;\n"
			+"\n"
			+"import shared.Timer;\n"
			+"import stream.Read;\n"
			+"\n"
			+"public final class TemplateArgSplitBBToolProbe extends BBTool_ST{\n"
			+"\tpublic TemplateArgSplitBBToolProbe(String[] args){super(args);}\n"
			+"\n"
			+"\tpublic void reparsePublic(String[] args){\n"
			+"\t\tloop=\"BBTool_ST.reparse\";\n"
			+"\t\ttry{\n"
			+"\t\t\treparse(args);\n"
			+"\t\t}finally{\n"
			+"\t\t\tloop=\"BBTool_ST.constructor\";\n"
			+"\t\t}\n"
			+"\t}\n"
			+"\n"
			+"\t@Override\n"
			+"\tprotected void setDefaults(){}\n"
			+"\n"
			+"\t@Override\n"
			+"\tpublic boolean parseArgument(String arg, String a, String b){\n"
			+"\t\tif(a.equals(\"jean_probe\")){"
			+"throw new TemplateArgSplitProbe(loop, arg, a, b);}\n"
			+"\t\treturn false;\n"
			+"\t}\n"
			+"\n"
			+"\t@Override\n"
			+"\tprotected void startupSubclass(){}\n"
			+"\n"
			+"\t@Override\n"
			+"\tprotected void shutdownSubclass(){}\n"
			+"\n"
			+"\t@Override\n"
			+"\tprotected void showStatsSubclass(final Timer t, long readsIn, "
			+"long basesIn){}\n"
			+"\n"
			+"\t@Override\n"
			+"\tprotected boolean processReadPair(final Read r1, final Read r2){\n"
			+"\t\treturn true;\n"
			+"\t}\n"
			+"\n"
			+"\t@Override\n"
			+"\tprotected boolean useSharedHeader(){return false;}\n"
			+"\n"
			+"\tprivate static String loop=\"BBTool_ST.constructor\";\n"
			+"}\n";
}
