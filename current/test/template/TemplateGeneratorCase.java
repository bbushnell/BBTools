package test.template;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.net.URL;
import java.net.URLClassLoader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicLong;

import stream.ConcurrentReadOutputStream;

/** Source-derived check for the generator sample's paired-read output accounting.
 * Prepare mode fills only the documented generation hook in an owned source copy;
 * run mode loads the guard-compiled copy and exercises the original worker
 * counter path without crossing explicit scaffold assertions.
 * @author Jean
 */
public final class TemplateGeneratorCase{

	/** Prepares patched source or checks a compiled patched generator. */
	public static void main(String[] args) throws Exception{
		final HashMap<String, String> map=parseArgs(args);
		final String mode=map.get("mode");
		if("prepare".equals(mode)){
			prepare(Paths.get(required(map, "source")), Paths.get(required(map, "out")));
		}else if("run".equals(mode)){
			run(Paths.get(required(map, "classes")), required(map, "expect"),
					fragments(map));
		}else{
			throw new IllegalArgumentException("Expected mode=prepare or mode=run");
		}
	}

	/** Writes the patched sample source under out/template/A_Sample_Generator.java. */
	private static void prepare(final Path source, final Path out) throws Exception{
		final String original=new String(Files.readAllBytes(source), StandardCharsets.UTF_8);
		final String patched=original.replace(GENERATE_READ_STUB, GENERATE_LITERAL_READS);
		check(!patched.equals(original), "generateRead replacement");

		final Path sourceFile=out.resolve("template/A_Sample_Generator.java");
		Files.createDirectories(sourceFile.getParent());
		check(!Files.exists(sourceFile), "prepared source already exists: "+sourceFile);
		Files.write(sourceFile, patched.getBytes(StandardCharsets.UTF_8));
		System.out.println("prepared\t"+sourceFile);
	}

	/** Runs the patched worker and checks baseline or fixed counter expectations. */
	private static void run(final Path classes, final String expect, final int fragments)
			throws Exception{
		if(!"baseline".equals(expect) && !"fixed".equals(expect)){
			throw new IllegalArgumentException("Expected expect=baseline or expect=fixed");
		}
		check(fragments==1 || fragments==2, "Expected fragments=1 or fragments=2");

		try(ChildFirstLoader loader=new ChildFirstLoader(classes.toUri().toURL(),
				TemplateGeneratorCase.class.desiredAssertionStatus())){
			final Class<?> workerClass=Class.forName(
					"template.A_Sample_Generator$ProcessThread", true, loader);
			final Constructor<?> constructor=workerClass.getDeclaredConstructor(
					ConcurrentReadOutputStream.class, int.class, long.class, AtomicLong.class);
			constructor.setAccessible(true);
			final Object worker=constructor.newInstance(null, 0, (long)fragments,
					new AtomicLong(0));

			final Method processInner=workerClass.getDeclaredMethod("processInner");
			processInner.setAccessible(true);
			processInner.invoke(worker);

			final long reads=longField(workerClass, worker, "readsOutT");
			final long bases=longField(workerClass, worker, "basesOutT");
			final long expectedReads=(fragments==1 ? 1 : 3);
			final long expectedBases=(fragments==1 ? 5 : ("fixed".equals(expect) ? 17 : 10));
			check(reads==expectedReads, "readsOutT", reads, expectedReads);
			check(bases==expectedBases, "basesOutT", bases, expectedBases);
			System.out.println("PASS\t"+expect+"\t"+fragments+"\t"+reads+"\t"+bases);
		}
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

	/** Parses the number of generated fragments for the worker run. */
	private static int fragments(final HashMap<String, String> map){
		final int fragments=Integer.parseInt(required(map, "fragments"));
		check(fragments==1 || fragments==2, "Expected fragments=1 or fragments=2");
		return fragments;
	}

	/** Reads a long counter from the otherwise private sample worker. */
	private static long longField(final Class<?> owner, final Object target, final String name)
			throws Exception{
		final Field field=owner.getDeclaredField(name);
		field.setAccessible(true);
		return field.getLong(target);
	}

	/** Remains enabled under both -ea and -da. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Remains enabled under both -ea and -da. */
	private static void check(final boolean condition, final String name, final long actual,
			final long expected){
		if(!condition){throw new IllegalStateException(name+": "+actual+" != "+expected);}
		System.out.println(name+"\t"+actual);
	}

	/** Loads the patched generator before any same-named class from the parent loader. */
	private static final class ChildFirstLoader extends URLClassLoader{
		ChildFirstLoader(final URL url, final boolean assertions){
			super(new URL[] {url}, TemplateGeneratorCase.class.getClassLoader());
			setDefaultAssertionStatus(assertions);
			setPackageAssertionStatus("template", assertions);
		}

		@Override
		protected Class<?> loadClass(final String name, final boolean resolve)
				throws ClassNotFoundException{
			if(name.startsWith("template.A_Sample_Generator")){
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

	private static final String GENERATE_READ_STUB=
			"\t\tprivate Read generateRead(final long nextID){\n"
			+"//\t\t\tRead r=new Read(null, null, nextID);\n"
			+"\t\t\tthrow new RuntimeException(\"TODO: Implement this method.\"); //TODO\n"
			+"\t\t}\n";

	private static final String GENERATE_LITERAL_READS=
			"\t\tprivate Read generateRead(final long nextID){\n"
			+"\t\t\tif(nextID==0){return new Read(new byte[] {'A', 'C', 'G', 'T', 'A'}, null, \"single\", nextID);}\n"
			+"\t\t\tRead r1=new Read(new byte[] {'A', 'C', 'G', 'T', 'A'}, null, \"literal/1\", nextID);\n"
			+"\t\t\tRead r2=new Read(new byte[] {'T', 'G', 'C', 'A', 'T', 'G', 'C'}, null, \"literal/2\", nextID);\n"
			+"\t\t\tr1.mate=r2;\n"
			+"\t\t\tr2.mate=r1;\n"
			+"\t\t\treturn r1;\n"
			+"\t\t}\n";
}
