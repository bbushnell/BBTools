package test.jgi;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Objects;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.AdjustHomopolymers;
import parse.Parse;
import parse.Parser;
import shared.Shared;
import stream.FASTQ;
import stream.Read;
import tracker.ReadStats;

/**
 * Constructor-only rate/parser regressions plus two ordinary direct-edit controls.
 * Invalid rates NEVER reach processRead or the driver. Private valid full/prefix
 * paths keep the equals regression away from fatal missing-input handling.
 * Baseline accepts four invalid rates and the placeholder, truncates equals paths,
 * and reports a null-value failure through Float.parseFloat rather than a rate diagnostic.
 * No driver processing, large arithmetic, invalid-rate editing or golden-test changes.
 * @author Sucrose
 */
public final class AdjustHomopolymersArgumentsTest{

	/** Requires phase=baseline|fixed and a new evidence directory without '='. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2, "Expected flag=value: "+arg);
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown argument fixture option: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null, "Set a new outdir.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		require(plainPath(root), "Only the explicit equals-path scenario may contain '='; evidence paths must fit one row.");
		require(AdjustHomopolymers.class.desiredAssertionStatus(), "The baseline zero-rate and unknown-option contracts require -ea.");
		Files.createDirectory(root);
		final AdjustHomopolymersArgumentsTest test=new AdjustHomopolymersArgumentsTest(root, phase.equals("baseline"));
		final Globals inherited=new Globals();
		try{
			for(final InvalidRate row : INVALID_RATES){test.invalidRate(row);}
			test.paths();
			test.lastOption();
			test.placeholder();
			test.direct("minus-one", "-1", "AACN", new byte[]{10, 20, 30, 0}, "N", new byte[]{0});
			test.direct("plus-two", "2", "A", new byte[]{7}, "AAA", new byte[]{7, 7, 7});
		}finally{inherited.restore();}
		test.check(inherited.matches(), "Restore all22 inherited process settings at fixture exit.");
		test.check(test.cases==11, "Frozen matrix requires11 scenarios: "+test.cases);
		test.check(test.constructorCalls==11, "Each scenario constructs exactly once: "+test.constructorCalls);
		test.check(test.readCalls==2, "Only ordinary -1/+2 controls may call processRead: "+test.readCalls);
		test.check(test.invalidAccepted==(test.baseline ? 4 : 0), "Expected4 baseline invalid-rate accepts or0 fixed: "+test.invalidAccepted);
		test.check(test.pathMismatches==(test.baseline ? 1 : 0), "Expected1 baseline equals-path mismatch or0 fixed: "+test.pathMismatches);
		test.check(test.placeholderAccepted==(test.baseline ? 1 : 0), "Expected1 baseline placeholder accept or0 fixed: "+test.placeholderAccepted);
		if(test.failures>0){throw new AssertionError("Argument fixture failed "+test.failures+" of "+test.checks+" checks; inspect "+root);}
		System.out.println("ADJUST_HOMOPOLYMERS_ARGUMENTS_PASS phase="+phase+" cases="+test.cases+
			" constructorCalls="+test.constructorCalls+" readCalls="+test.readCalls+" invalidAccepted="+test.invalidAccepted+
			" pathMismatches="+test.pathMismatches+" placeholderAccepted="+test.placeholderAccepted+" checks="+test.checks);
	}

	private AdjustHomopolymersArgumentsTest(final Path root_, final boolean baseline_) throws Exception{
		root=root_; baseline=baseline_;
		edit=AdjustHomopolymers.class.getDeclaredMethod("processRead", Read.class); edit.setAccessible(true);
		write(root.resolve("cases.tsv"), "case\tphaseExpected\tactual\n");
	}

	/** Constructor only, even when the baseline unexpectedly accepts a dangerous value. */
	private void invalidRate(final InvalidRate row) throws Exception{
		try(final Work work=new Work(row.name, false)){
			final Outcome result=construct(work, "rate="+row.value);
			final Class<?> expected=baseline ? row.baselineFailure : IllegalArgumentException.class;
			boolean correct=outcome(result, expected, expected==NullPointerException.class ? null : "rate");
			if(baseline && expected==NullPointerException.class){
				correct&=frame(result.error, Float.class.getName(), "parseFloat");
			}
			check(correct, row.name+": expected "+typeName(expected)+" with the rate-specific phase contract.");
			if(result.app!=null){invalidAccepted++;}
			record(work, result, typeName(expected));
		}
	}

	/** Both valid input names and both output sentinels exist; construction must preserve all four. */
	private void paths() throws Exception{
		try(final Work work=new Work("equals-paths", true)){
			final Outcome result=construct(work, "rate=0.5");
			check(outcome(result, null, null), "Equals paths must construct without an acquisition failure.");
			require(result.app!=null, "Cannot inspect paths after a failed constructor.");
			final String[] intended=pathValues(work.input, work.output), expected=baseline ? pathValues(work.inputPrefix, work.outputPrefix) : intended;
			final String[] actual=pathValues(result.app);
			check(Arrays.equals(actual, expected), "Parsed paths/descriptors must match the literal phase oracle.");
			if(!Arrays.equals(actual, intended)){pathMismatches++;}
			write(work.directory.resolve("paths.tsv"), pathTable(intended, expected, actual));
			record(work, result, "returned");
		}
	}

	/** Validate the final selected rate, preserving an earlier zero overridden by a valid value. */
	private void lastOption() throws Exception{
		try(final Work work=new Work("last-option", false)){
			final Outcome result=construct(work, "rate=0", "rate=0.5");
			check(outcome(result, null, null), "Last valid rate must override earlier zero before final validation.");
			require(result.app!=null, "Cannot inspect settings after a failed constructor.");
			final float actual=(Float)value(result.app, "rate");
			check(Float.floatToIntBits(actual)==Float.floatToIntBits(0.5f), "Final rate must be exactly0.5, observed="+actual);
			check(Arrays.equals(pathValues(result.app), pathValues(work.input, work.output)), "Ordinary input/output names and descriptors must be preserved.");
			write(work.directory.resolve("rate.tsv"), "expected\tactual\n0.5\t"+actual+"\n");
			record(work, result, "returned");
		}
	}

	private void placeholder() throws Exception{
		try(final Work work=new Work("placeholder", false)){
			final Outcome result=construct(work, "rate=0.5", "parse_flag_goes_here=1");
			check(outcome(result, baseline ? null : AssertionError.class, "parse_flag_goes_here"),
				"Placeholder must follow the existing unknown-option rejection under -ea after its branch is removed.");
			if(result.app!=null){placeholderAccepted++;}
			record(work, result, baseline ? "returned" : AssertionError.class.getName());
		}
	}

	/** Finite ordinary endpoint/control only; invalid-rate cases have no path into this method. */
	private void direct(final String name, final String rate, final String input, final byte[] quality,
		final String expected, final byte[] expectedQuality) throws Exception{
		try(final Work work=new Work(name, false)){
			final Outcome result=construct(work, "rate="+rate);
			check(outcome(result, null, null), name+": ordinary control must construct.");
			require(result.app!=null, "No editing call is allowed after constructor failure.");
			final Read read=new Read(bytes(input), quality.clone(), name, 42, false);
			final Metadata saved=new Metadata(read);
			try{edit.invoke(result.app, read); readCalls++;}
			catch(InvocationTargetException e){
				final Throwable cause=e.getCause();
				if(cause instanceof Exception){throw (Exception)cause;}
				if(cause instanceof Error){throw (Error)cause;}
				throw e;
			}
			check(Arrays.equals(read.bases, bytes(expected)), name+": exact literal sequence oracle.");
			check(Arrays.equals(read.quality, expectedQuality), name+": exact literal Phred-quality oracle.");
			check(saved.unchanged(), name+": metadata and mate identities must remain unchanged.");
			check(read.bases.length==read.quality.length, name+": output must retain one quality per base.");
			write(work.directory.resolve("observed.tsv"), "input\tinputQualities\texpected\texpectedQualities\tactual\tactualQualities\n"+
				input+"\t"+Arrays.toString(quality)+"\t"+expected+"\t"+Arrays.toString(expectedQuality)+"\t"+
				new String(read.bases, StandardCharsets.US_ASCII)+"\t"+Arrays.toString(read.quality)+"\n");
			record(work, result, "returned");
		}
	}

	private Outcome construct(final Work work, final String... options) throws Exception{
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+work.input, "out="+work.output,
			"int=f", "t=1", "simd=f", "qin=33", "qout=33", "ow=t", "app=f"));
		args.addAll(Arrays.asList(options));
		Files.write(work.directory.resolve("constructor.argv"), args, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		constructorCalls++;
		try{return new Outcome(new AdjustHomopolymers(args.toArray(new String[args.size()])), null);}
		catch(RuntimeException e){return new Outcome(null, e);}
		catch(AssertionError e){return new Outcome(null, e);}
	}

	private static boolean outcome(final Outcome result, final Class<?> expected, final String fragment){
		if(expected==null){return result.app!=null && result.error==null;}
		return result.app==null && result.error!=null && result.error.getClass()==expected &&
			(frame(result.error, AdjustHomopolymers.class.getName(), "<init>") || frame(result.error, AdjustHomopolymers.class.getName(), "parse")) &&
			(fragment==null || String.valueOf(result.error.getMessage()).toLowerCase(Locale.ROOT).contains(fragment));
	}
	private void record(final Work work, final Outcome result, final String expected) throws Exception{
		final String actual=result.error==null ? "returned" : result.error.getClass().getName();
		write(work.directory.resolve("outcome.tsv"), "phaseExpected\tactual\tmessage\n"+expected+"\t"+actual+"\t"+
			cell(result.error==null ? "" : String.valueOf(result.error.getMessage()))+"\n");
		if(result.error!=null){
			final StringWriter text=new StringWriter(); final PrintWriter writer=new PrintWriter(text);
			result.error.printStackTrace(writer); writer.flush();
			write(work.directory.resolve("failure.txt"), text.toString());
		}
		append(root.resolve("cases.tsv"), work.name+"\t"+expected+"\t"+actual+"\n");
	}
	private static boolean frame(final Throwable error, final String owner, final String method){
		if(error==null){return false;}
		for(final StackTraceElement row : error.getStackTrace()){if(row.getClassName().equals(owner) && row.getMethodName().equals(method)){return true;}}
		return false;
	}
	private static Object value(final Object object, final String name) throws Exception{
		final Field field=AdjustHomopolymers.class.getDeclaredField(name); field.setAccessible(true); return field.get(object);
	}
	private static String[] pathValues(final AdjustHomopolymers app) throws Exception{
		return new String[]{(String)value(app, "in1"), (String)value(app, "out1"),
			((FileFormat)value(app, "ffin1")).name(), ((FileFormat)value(app, "ffout1")).name()};
	}
	private static String[] pathValues(final Path input, final Path output){
		return new String[]{input.toString(), output.toString(), input.toString(), output.toString()};
	}
	private static String pathTable(final String[] intended, final String[] expected, final String[] actual){
		final StringBuilder out=new StringBuilder("field\tintended\tphaseExpected\tactual\n");
		final String[] fields={"in1", "out1", "ffin1.name", "ffout1.name"};
		for(int i=0; i<fields.length; i++){out.append(fields[i]).append('\t').append(intended[i]).append('\t').append(expected[i]).append('\t').append(actual[i]).append('\n');}
		return out.toString();
	}

	/** Owns private preflight files and restores all touched process settings, including failed constructors. */
	private final class Work implements AutoCloseable{
		Work(final String name_, final boolean equals) throws Exception{
			name=name_; directory=Files.createDirectory(root.resolve(name)); globals=new Globals();
			inputPrefix=directory.resolve("input.fq"); outputPrefix=directory.resolve("output.fq");
			input=equals ? directory.resolve("input.fq=one=two.fq") : inputPrefix;
			output=equals ? directory.resolve("output.fq=one=two.fq") : outputPrefix;
			fullPaths=equals;
			write(inputPrefix, INPUT_PREFIX);
			if(fullPaths){write(input, INPUT_FULL); write(outputPrefix, OUTPUT_PREFIX); write(output, OUTPUT_FULL);}
			write(directory.resolve("globals-before.tsv"), globals.table());
			write(directory.resolve("paths-before.tsv"), "role\tpath\tpreserved\n");
			require(capture("before", true) && capture("before", false), "Private literal file setup must be intact before construction.");
		}
		@Override
		public void close() throws Exception{
			try{
				write(directory.resolve("paths-after.tsv"), "role\tpath\tpreserved\n");
				check(capture("after", true), name+": preserve both original/prefix input byte strings.");
				check(capture("after", false), name+": preserve output sentinels or leave ordinary output absent.");
			}finally{
				globals.restore();
				write(directory.resolve("globals-restored.tsv"), new Globals().table());
				check(globals.matches(), name+": restore all22 captured settings after terminal observation.");
				cases++;
			}
		}
		private boolean capture(final String stage, final boolean inputs) throws Exception{
			boolean same=observe(stage, inputs ? "input-prefix" : "output-prefix",
				inputs ? inputPrefix : outputPrefix, inputs ? INPUT_PREFIX : (fullPaths ? OUTPUT_PREFIX : null));
			if(fullPaths){same&=observe(stage, inputs ? "input-full" : "output-full", inputs ? input : output, inputs ? INPUT_FULL : OUTPUT_FULL);}
			return same;
		}
		private boolean observe(final String stage, final String role, final Path path, final String expected) throws Exception{
			final boolean regular=Files.isRegularFile(path);
			final byte[] actual=regular ? Files.readAllBytes(path) : null;
			final boolean same=expected==null ? !Files.exists(path) : regular && Arrays.equals(actual, bytes(expected));
			if(regular){Files.write(directory.resolve(role+"."+stage+".bin"), actual, StandardOpenOption.CREATE_NEW);}
			append(directory.resolve("paths-"+stage+".tsv"), role+"\t"+path+"\t"+same+"\n");
			return same;
		}
		final String name;
		final Path directory, input, output, inputPrefix, outputPrefix;
		final boolean fullPaths;
		final Globals globals;
	}

	private static final class Outcome{
		Outcome(final AdjustHomopolymers app_, final Throwable error_){app=app_; error=error_;}
		final AdjustHomopolymers app;
		final Throwable error;
	}
	private static final class InvalidRate{
		InvalidRate(final String name_, final String value_, final Class<?> baselineFailure_){name=name_; value=value_; baselineFailure=baselineFailure_;}
		final String name, value;
		final Class<?> baselineFailure;
	}
	private static final class Metadata{
		Metadata(final Read read_){
			read=read_; id=read.id; mate=read.mate; numericID=read.numericID; flags=read.flags;
			chrom=read.chrom; start=read.start; stop=read.stop; insert=read.insert();
		}
		boolean unchanged(){return read.id==id && read.mate==mate && read.numericID==numericID && read.flags==flags &&
			read.chrom==chrom && read.start==start && read.stop==stop && read.insert()==insert;}
		final Read read, mate;
		final String id;
		final long numericID;
		final int flags, chrom, start, stop, insert;
	}

	/** Explicit snapshot matches the existing golden fixture's22 named constructor settings. */
	private static final class Globals{
		Globals() throws Exception{
			add(FASTQ.class, "TEST_INTERLEAVED", "FORCE_INTERLEAVED", "PARSE_CUSTOM", "ASCII_OFFSET", "ASCII_OFFSET_OUT",
				"DETECT_QUALITY", "DETECT_QUALITY_OUT", "SET_QIN");
			add(Shared.class, "THREADS", "READ_BUFFER_NUM_BUFFERS", "SIMD", "COMMAND_LINE", "mainClass");
			add(ReadWrite.class, "USE_PIGZ", "USE_UNPIGZ", "MAX_ZIP_THREADS", "SET_ZIP_THREADS");
			add(Parser.class, "qin", "qout", "parsedQuality");
			add(ReadStats.class, "overwrite", "append");
		}
		private void add(final Class<?> type, final String... names) throws Exception{
			for(final String name : names){
				final Field field=type.getDeclaredField(name); field.setAccessible(true);
				fields.add(field); values.add(field.get(null));
			}
		}
		void restore() throws Exception{for(int i=0; i<fields.size(); i++){fields.get(i).set(null, values.get(i));}}
		boolean matches() throws Exception{
			for(int i=0; i<fields.size(); i++){if(!Objects.equals(values.get(i), fields.get(i).get(null))){return false;}}
			return fields.size()==22;
		}
		String table(){
			final StringBuilder out=new StringBuilder("setting\tvalue\n");
			for(int i=0; i<fields.size(); i++){
				final Object value=values.get(i);
				out.append(fields.get(i).getDeclaringClass().getSimpleName()).append('.').append(fields.get(i).getName()).append('\t')
					.append(cell(value instanceof String[] ? Arrays.toString((String[])value) : String.valueOf(value))).append('\n');
			}
			return out.toString();
		}
		private final ArrayList<Field> fields=new ArrayList<Field>();
		private final ArrayList<Object> values=new ArrayList<Object>();
	}

	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static String cell(final String text){return text.replace("\t", "\\t").replace("\n", "\\n").replace("\r", "\\r");}
	private static String typeName(final Class<?> type){return type==null ? "returned" : type.getName();}
	private static boolean plainPath(final Path path){
		final String text=path.toString(); return text.indexOf('=')<0 && text.indexOf('\t')<0 && text.indexOf('\n')<0 && text.indexOf('\r')<0;
	}
	private static void write(final Path path, final String text) throws Exception{Files.write(path, bytes(text), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path path, final String text) throws Exception{Files.write(path, bytes(text), StandardOpenOption.APPEND);}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}
	private void check(final boolean condition, final String message){checks++; if(!condition){failures++; System.err.println("CHECK_FAILED "+message);}}

	private static final InvalidRate[] INVALID_RATES={
		new InvalidRate("zero", "0", AssertionError.class),
		new InvalidRate("nan", "NaN", null),
		new InvalidRate("positive-infinity", "Infinity", null),
		new InvalidRate("negative-infinity", "-Infinity", null),
		new InvalidRate("below-minus-one", "-1.01", null),
		new InvalidRate("missing-value", "", NullPointerException.class)
	};
	private static final String INPUT_PREFIX="@prefix\nA\n+\nI\n", INPUT_FULL="@full\nC\n+\nJ\n";
	private static final String OUTPUT_PREFIX="prefix-output-sentinel\n", OUTPUT_FULL="full-output-sentinel\n";
	private final Path root;
	private final boolean baseline;
	private final Method edit;
	private int cases, constructorCalls, readCalls, invalidAccepted, pathMismatches, placeholderAccepted, checks, failures;
}
