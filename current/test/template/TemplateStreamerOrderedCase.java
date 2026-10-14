package test.template;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import fileIO.FileFormat;
import template.A_SampleStreamer;
import template.A_SampleStreamerMT;

/** Constructor-local checks for Streamer template output ordering descriptors.
 * @author Jean
 */
public final class TemplateStreamerOrderedCase{

	/** Runs finite constructor cases without invoking process or writing output. */
	public static void main(String[] args) throws Exception{
		Path work=null;
		String expect=null;
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			check(split.length==2, "Arguments require key=value: "+arg);
			if(split[0].equals("work")){work=Paths.get(split[1]).toAbsolutePath().normalize();}
			else if(split[0].equals("expect")){expect=split[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		check(work!=null, "work= is required");
		check("baseline".equals(expect) || "fixed".equals(expect),
				"expect=baseline or expect=fixed is required");
		check(!Files.exists(work), "Fresh work directory required: "+work);
		Files.createDirectories(work);

		noOutput(work, A_SampleStreamerMT.class, "streamer-mt-no-out", expect);
		singleOutput(work, A_SampleStreamerMT.class, "streamer-mt-ordered-f", false, expect);
		singleOutput(work, A_SampleStreamerMT.class, "streamer-mt-ordered-t", true, expect);
		pairedOutput(work, A_SampleStreamerMT.class, "streamer-mt-paired", expect);

		noOutput(work, A_SampleStreamer.class, "streamer-no-out", expect);
		singleOutput(work, A_SampleStreamer.class, "streamer-ordered-f", false, expect);
		singleOutput(work, A_SampleStreamer.class, "streamer-ordered-t", true, expect);
		pairedOutput(work, A_SampleStreamer.class, "streamer-paired", expect);
	}

	/** Verifies that a false ordered request is parsed when output is absent. */
	private static void noOutput(final Path work, final Class<?> type, final String name,
			final String expect) throws Exception{
		final Constructed observed=construct(type, work, name, false, false);
		check(observed.ordered==false, name+" should parse ordered=false");
		check(observed.out1==null, name+" should leave ffout1 null");
		check(observed.out2==null, name+" should leave ffout2 null");
		pass(name, "ffout1=null,ffout2=null", expect);
	}

	/** Verifies single-output ordering metadata for true and false requests. */
	private static void singleOutput(final Path work, final Class<?> type, final String name,
			final boolean ordered, final String expect) throws Exception{
		final Constructed observed=construct(type, work, name, true, false, ordered);
		final boolean expected=ordered || "baseline".equals(expect);
		check(observed.ordered==ordered, name+" parsed ordered "+observed.ordered);
		check(observed.out1!=null, name+" should set ffout1");
		check(observed.out1.ordered()==expected, name+" ffout1 "+flags(observed.out1));
		check(observed.out2==null, name+" should leave ffout2 null");
		pass(name, flags(observed.out1)+",ffout2=null", expect);
	}

	/** Preserves the current paired-output ordered descriptors for WriterFactory. */
	private static void pairedOutput(final Path work, final Class<?> type, final String name,
			final String expect) throws Exception{
		final Constructed observed=construct(type, work, name, true, true, false);
		check(observed.ordered==false, name+" should parse ordered=false");
		check(observed.out1!=null, name+" should set ffout1");
		check(observed.out2!=null, name+" should set ffout2");
		check(observed.out1.ordered(), name+" ffout1 "+flags(observed.out1));
		check(observed.out2.ordered(), name+" ffout2 "+flags(observed.out2));
		pass(name, flags(observed.out1)+","+flags(observed.out2), expect);
	}

	/** Constructs one template and returns its private ordered output state. */
	private static Constructed construct(final Class<?> type, final Path work, final String name,
			final boolean output, final boolean paired) throws Exception{
		return construct(type, work, name, output, paired, false);
	}

	/** Constructs one template and returns its private ordered output state. */
	private static Constructed construct(final Class<?> type, final Path work, final String name,
			final boolean output, final boolean paired, final boolean ordered) throws Exception{
		final Path in=writeFastq(work.resolve(name+".fq"));
		final String[] args;
		if(paired){
			args=new String[] {"in="+in, "out="+work.resolve(name+"-1.fq"),
				"out2="+work.resolve(name+"-2.fq"), "ordered="+ordered};
		}else if(output){
			args=new String[] {"in="+in, "out="+work.resolve(name+".out.fq"),
				"ordered="+ordered};
		}else{
			args=new String[] {"in="+in, "ordered="+ordered};
		}
		final Object sample=type.getConstructor(String[].class).newInstance((Object)args);
		return new Constructed(ordered(sample, type), output(sample, type, "ffout1"),
				output(sample, type, "ffout2"));
	}

	/** Returns a private output descriptor. */
	private static FileFormat output(final Object sample, final Class<?> type,
			final String fieldName) throws Exception{
		final Field field=type.getDeclaredField(fieldName);
		field.setAccessible(true);
		return (FileFormat)field.get(sample);
	}

	/** Returns the private parsed ordered flag. */
	private static boolean ordered(final Object sample, final Class<?> type) throws Exception{
		final Field field=type.getDeclaredField("ordered");
		field.setAccessible(true);
		return field.getBoolean(sample);
	}

	/** Writes a tiny valid FASTQ input and returns its path. */
	private static Path writeFastq(final Path path) throws Exception{
		Files.write(path, FASTQ);
		return path;
	}

	/** Formats the ordered flag under test. */
	private static String flags(final FileFormat ff){return "ordered="+ff.ordered();}

	/** Prints one successful case row. */
	private static void pass(final String name, final String observed, final String expect){
		System.out.println("PASS\t"+expect+"\t"+name+"\t"+observed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Constructor-local observations. */
	private static final class Constructed{

		/** Stores constructor observations. */
		Constructed(final boolean ordered_, final FileFormat out1_, final FileFormat out2_){
			ordered=ordered_;
			out1=out1_;
			out2=out2_;
		}

		/** Parsed ordered flag. */
		final boolean ordered;
		/** Primary output descriptor, or null when output is disabled. */
		final FileFormat out1;
		/** Secondary output descriptor, or null when absent. */
		final FileFormat out2;
	}

	/** Small valid FASTQ payload used only so paths exist. */
	private static final byte[] FASTQ="@r0\nACGT\n+\n!!!!\n".getBytes(StandardCharsets.UTF_8);
}
