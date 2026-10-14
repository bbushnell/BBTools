package test.template;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import fileIO.FileFormat;
import parse.Parser;
import template.A_SampleMultipleInput;

/** Constructor-local checks for A_SampleMultipleInput input and output parsing.
 * @author Sucrose
 */
public final class TemplateMultipleInputCase{

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

		namedInput(work, expect);
		positionalInput(work, expect);
		nullInput(work, expect);
		emptyInput(work, expect);
		existingOutputOverwriteFalse(work, expect);
		appendTrue(work, expect);
	}

	/** Verifies an explicit in= path is retained as the sole input. */
	private static void namedInput(final Path work, final String expect) throws Exception{
		final Path in=writeFastq(work.resolve("named.fq"));
		final Constructed observed=construct("in="+in);
		check(observed.thrown==null, "named input threw "+observed.thrown);
		checkList("named", observed.inputs, in.toString());
		pass("named", "ok", expect);
	}

	/** Verifies a positional existing file is retained as the sole input. */
	private static void positionalInput(final Path work, final String expect) throws Exception{
		final Path in=writeFastq(work.resolve("positional.fq"));
		final Constructed observed=construct(in.toString());
		check(observed.thrown==null, "positional input threw "+observed.thrown);
		checkList("positional", observed.inputs, in.toString());
		pass("positional", "ok", expect);
	}

	/** Captures the current in=null constructor failure. */
	private static void nullInput(final Path work, final String expect) throws Exception{
		final Constructed observed=construct("in=null");
		if("baseline".equals(expect)){
			check(observed.thrown instanceof NullPointerException,
					"baseline in=null should throw NullPointerException, got "+observed.thrown);
		}else{
			checkFailure(observed, IllegalArgumentException.class,
					"in= requires a non-empty input path", "fixed in=null");
		}
		pass("null-input", thrownName(observed), expect);
	}

	/** Captures the current in= constructor failure. */
	private static void emptyInput(final Path work, final String expect) throws Exception{
		final Constructed observed=construct("in=");
		if("baseline".equals(expect)){
			check(observed.thrown instanceof NullPointerException,
					"baseline in= should throw NullPointerException, got "+observed.thrown);
		}else{
			checkFailure(observed, IllegalArgumentException.class,
					"in= requires a non-empty input path", "fixed in=");
		}
		pass("empty-input", thrownName(observed), expect);
	}

	/** Captures descriptor overwrite policy separately from actual output writes. */
	private static void existingOutputOverwriteFalse(final Path work, final String expect)
			throws Exception{
		final Path in=writeFastq(work.resolve("ow-input.fq"));
		final Path out=work.resolve("existing.fq");
		Files.write(out, SENTINEL);
		final Constructed observed=construct("in="+in, "out="+out, "ow=f");
		check(observed.thrown==null, "ow=f constructor threw "+observed.thrown);
		checkList("ow=false", observed.inputs, in.toString());
		checkFlags(parsedOutput(out, "ow=f"), false, false,
				"Parser overwrite=false passed to FileFormat");
		if("baseline".equals(expect)){
			check(observed.output.overwrite(), "baseline descriptor should force overwrite");
			check(!observed.output.append(), "baseline descriptor should force append=false");
		}else{
			check(!observed.output.overwrite(), "fixed descriptor should honor overwrite=false");
			check(!observed.output.append(), "fixed descriptor should keep append=false");
		}
		checkBytes("existing output sentinel", SENTINEL, Files.readAllBytes(out));
		pass("overwrite-false", flags(observed.output), expect);
	}

	/** Captures descriptor append policy separately from actual output writes. */
	private static void appendTrue(final Path work, final String expect) throws Exception{
		final Path in=writeFastq(work.resolve("append-input.fq"));
		final Path out=work.resolve("append.fq");
		final Constructed observed=construct("in="+in, "out="+out, "append=t");
		check(observed.thrown==null, "append=t constructor threw "+observed.thrown);
		checkList("append", observed.inputs, in.toString());
		checkFlags(parsedOutput(out, "append=t"), false, true,
				"FileFormat-normalized Parser append=true");
		if("baseline".equals(expect)){
			check(observed.output.overwrite(), "baseline descriptor should force overwrite");
			check(!observed.output.append(), "baseline descriptor should force append=false");
		}else{
			check(!observed.output.overwrite(), "fixed descriptor should clear overwrite");
			check(observed.output.append(), "fixed descriptor should honor append=true");
		}
		check(!Files.exists(out), "constructor-only descriptor check must not write output");
		pass("append-true", flags(observed.output), expect);
	}

	/** Constructs A_SampleMultipleInput and returns private constructor state. */
	private static Constructed construct(final String... args) throws Exception{
		try{
			final A_SampleMultipleInput sample=new A_SampleMultipleInput(args);
			return new Constructed(inputs(sample), output(sample), null);
		}catch(Throwable t){return new Constructed(null, null, t);}
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

	/** Returns the private input list. */
	@SuppressWarnings("unchecked")
	private static ArrayList<String> inputs(final A_SampleMultipleInput sample)
			throws Exception{
		final Field field=A_SampleMultipleInput.class.getDeclaredField("in");
		field.setAccessible(true);
		return (ArrayList<String>)field.get(sample);
	}

	/** Returns the private output descriptor. */
	private static FileFormat output(final A_SampleMultipleInput sample) throws Exception{
		final Field field=A_SampleMultipleInput.class.getDeclaredField("ffout1");
		field.setAccessible(true);
		return (FileFormat)field.get(sample);
	}

	/** Writes a tiny valid FASTQ input and returns its path. */
	private static Path writeFastq(final Path path) throws Exception{
		Files.write(path, FASTQ);
		return path;
	}

	/** Checks exact input-list contents. */
	private static void checkList(final String name, final ArrayList<String> observed,
			final String... expected){
		check(observed!=null, name+" input list is null");
		check(observed.size()==expected.length, name+" input count "+observed.size());
		for(int i=0; i<expected.length; i++){
			check(expected[i].equals(observed.get(i)), name+" input "+i+" = "+observed.get(i));
		}
	}

	/** Checks exact bytes. */
	private static void checkBytes(final String name, final byte[] expected,
			final byte[] observed){
		check(expected.length==observed.length, name+" length");
		for(int i=0; i<expected.length; i++){
			check(expected[i]==observed[i], name+" byte "+i);
		}
	}

	/** Checks a captured failure's exact type and diagnostic. */
	private static void checkFailure(final Constructed observed,
			final Class<? extends Throwable> expectedType, final String expectedMessage,
			final String name){
		check(observed.thrown!=null, name+" should throw "+expectedType.getName());
		check(observed.thrown.getClass()==expectedType,
				name+" threw "+observed.thrown.getClass().getName());
		check(expectedMessage.equals(observed.thrown.getMessage()),
				name+" message = "+observed.thrown.getMessage());
	}

	/** Checks descriptor flags. */
	private static void checkFlags(final FileFormat ff, final boolean overwrite,
			final boolean append, final String name){
		check(ff.overwrite()==overwrite, name+" overwrite flag "+flags(ff));
		check(ff.append()==append, name+" append flag "+flags(ff));
	}

	/** Formats the descriptor flags under test. */
	private static String flags(final FileFormat ff){
		return "ow="+ff.overwrite()+",append="+ff.append();
	}

	/** Formats a captured throwable type. */
	private static String thrownName(final Constructed observed){
		return observed.thrown==null ? "none" : observed.thrown.getClass().getName();
	}

	/** Prints one successful case row. */
	private static void pass(final String name, final String observed, final String expect){
		System.out.println("PASS\t"+expect+"\t"+name+"\t"+observed);
	}

	/** Lightweight assertion independent of the JVM assertion flag. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Constructor result or thrown exception. */
	private static final class Constructed{

		/** Stores constructor observations. */
		Constructed(final ArrayList<String> inputs_, final FileFormat output_,
				final Throwable thrown_){
			inputs=inputs_;
			output=output_;
			thrown=thrown_;
		}

		/** Parsed input list, or null when construction failed. */
		final ArrayList<String> inputs;
		/** Output descriptor, or null when construction failed or output is disabled. */
		final FileFormat output;
		/** Constructor throwable, or null on success. */
		final Throwable thrown;
	}

	/** Small valid FASTQ payload used only so paths exist. */
	private static final byte[] FASTQ="@r0\nACGT\n+\n!!!!\n".getBytes(StandardCharsets.UTF_8);
	/** Bytes used to prove constructor-only checks do not write output. */
	private static final byte[] SENTINEL="sentinel\n".getBytes(StandardCharsets.UTF_8);
}
