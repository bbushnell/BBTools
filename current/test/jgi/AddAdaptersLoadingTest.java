package test.jgi;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Objects;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import jgi.AddAdapters;
import parse.Parse;
import parse.Parser;
import shared.Shared;
import stream.FASTQ;
import stream.FastaReadInputStream;
import tracker.ReadStats;

/**
 * Twelve constructor-only adapter-loading and literal equals-path baseline/fixed cases.
 * Main input/output names remain absent; only two private FASTA references are read.
 * An existing directory causes a catchable FileInputStream failure; a missing reference
 * is forbidden because ReadWrite would terminate the VM through KillSwitch.
 * Each case restores the named process settings even when construction fails.
 * @author Sucrose
 */
public final class AddAdaptersLoadingTest{

	/** Requires phase=baseline|fixed and a new outdir; never calls the driver's read/write methods. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown constructor fixture argument: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null, "Set a new outdir for private constructor evidence.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		final String name=root.toString();
		require(name.indexOf('=')<0 && name.indexOf('\t')<0 && name.indexOf('\n')<0 && name.indexOf('\r')<0,
			"Only the explicit equals-path cases may contain '='; evidence paths must fit one row.");
		require(FastaReadInputStream.MIN_READ_LEN==1 && !FastaReadInputStream.SPLIT_READS,
			"Two-base FASTA oracles require MIN_READ_LEN=1 and SPLIT_READS=false in a fresh JVM.");
		Files.createDirectory(root);
		final AddAdaptersLoadingTest test=new AddAdaptersLoadingTest(root, phase.equals("baseline"));
		final Globals inherited=new Globals();
		try{
			for(int i=0; i<4; i++){test.flags((i&2)!=0, (i&1)!=0, false);}
			for(int i=0; i<4; i++){test.flags((i&2)!=0, (i&1)!=0, true);}
			for(int i=0; i<4; i++){test.paths(i);}
		}finally{inherited.restore();}
		test.check(inherited.matches(), "Fixture must restore all21 captured global settings at exit.");
		test.check(test.cases==12, "Frozen matrix requires12 constructor cases, observed="+test.cases);
		test.check(test.flagMismatches==(test.baseline ? 3 : 0), "Expected3 baseline flag-case mismatches or0 fixed, observed="+test.flagMismatches);
		test.check(test.pathMismatches==(test.baseline ? 2 : 0), "Expected2 baseline path-case mismatches or0 fixed, observed="+test.pathMismatches);
		if(test.failures>0){throw new AssertionError("AddAdapters loading failed "+test.failures+" of "+test.checks+" checks; inspect "+root);}
		System.out.println("ADD_ADAPTERS_LOADING_PASS phase="+phase+" cases="+test.cases+" flagMismatches="+
			test.flagMismatches+" pathMismatches="+test.pathMismatches+" checks="+test.checks);
	}

	private AddAdaptersLoadingTest(final Path root_, final boolean baseline_) throws Exception{
		root=root_; baseline=baseline_;
		input=root.resolve("primary.fq"); output=root.resolve("output.fq");
		equalInput=root.resolve("primary.fq=one=two.fq"); equalOutput=root.resolve("output.fq=one=two.fq");
		badReference=Files.createDirectory(root.resolve("bad-adapter.fa"));
		reference=root.resolve("adapter.fa"); equalReference=root.resolve("adapter.fa=kept.fa");
		write(reference, REF_AC); write(equalReference, REF_TG);
		write(root.resolve("reference-prefix-before.fa"), REF_AC); write(root.resolve("reference-full-before.fa"), REF_TG);
		write(root.resolve("cases.tsv"), "case\tintended\tphaseExpected\tactual\tmatchesIntended\n");
		require(mainPathsAbsent(), "All full and truncated main paths must be absent before construction.");
	}

	/** The failed constructor is synchronous; observe flags after its caught terminal exception. */
	private void flags(final boolean ti, final boolean fi, final boolean reject) throws Exception{
		final String label=(reject ? "directory-" : "literal-")+ti+"-"+fi;
		final Path work=Files.createDirectory(root.resolve(label));
		final Globals saved=new Globals();
		try{
			write(work.resolve("globals-before.tsv"), saved.table());
			FASTQ.TEST_INTERLEAVED=ti; FASTQ.FORCE_INTERLEAVED=fi;
			require(!reject || Files.isDirectory(badReference), "The failing reference must remain an existing directory.");
			final String[] args=arguments(input, output, reject ? badReference : null, "AC");
			Files.write(work.resolve("constructor.argv"), Arrays.asList(args), StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
			AddAdapters app=null;
			RuntimeException failure=null;
			try{app=new AddAdapters(args);}
			catch(RuntimeException e){
				if(!reject || !acquisitionFailure(e, badReference.toString())){throw e;}
				failure=e;
			}
			final boolean actualTI=FASTQ.TEST_INTERLEAVED, actualFI=FASTQ.FORCE_INTERLEAVED;
			final boolean expectedTI=baseline && reject ? false : ti, expectedFI=baseline && reject ? false : fi;
			final boolean intended=actualTI==ti && actualFI==fi;
			if(!intended){flagMismatches++;}
			write(work.resolve("flags.tsv"), "state\tTEST_INTERLEAVED\tFORCE_INTERLEAVED\nentry\t"+ti+"\t"+fi+
				"\nintended\t"+ti+"\t"+fi+"\nphaseExpected\t"+expectedTI+"\t"+expectedFI+
				"\nactual\t"+actualTI+"\t"+actualFI+"\n");
			check(reject ? failure!=null && app==null : failure==null && app!=null, label+": wrong constructor outcome.");
			check(actualTI==expectedTI && actualFI==expectedFI, label+": pairing flags differ from literal phase oracle.");
			if(reject){
				check(acquisitionFailure(failure, badReference.toString()), label+": expected directory FileInputStream acquisition failure.");
				final StringWriter stack=new StringWriter();
				if(failure!=null){failure.printStackTrace(new PrintWriter(stack));}
				else{stack.write("Constructor unexpectedly returned instead of rejecting the adapter directory.\n");}
				write(work.resolve("failure.txt"), stack.toString());
			}else{
				check(adapterText(app).equals("[AC]"), label+": literal-only adapter must be exactly[AC].");
				write(work.resolve("adapters.txt"), adapterText(app)+"\n");
			}
			check(mainPathsAbsent(), label+": constructor must leave all main input/output names absent.");
			write(work.resolve("terminal.tsv"), terminalTable());
			append(root.resolve("cases.tsv"), label+"\t"+ti+","+fi+"\t"+expectedTI+","+expectedFI+
				"\t"+actualTI+","+actualFI+"\t"+intended+"\n");
			cases++;
		}finally{
			saved.restore();
			write(work.resolve("globals-restored.tsv"), new Globals().table());
			check(saved.matches(), label+": restore every captured process setting after observation.");
		}
	}

	/** Literal full/prefix references both exist so truncation never enters the fatal missing-file branch. */
	private void paths(final int kind) throws Exception{
		final String label=PATH_CASES[kind];
		final Path work=Files.createDirectory(root.resolve(label));
		final Globals saved=new Globals();
		try{
			write(work.resolve("globals-before.tsv"), saved.table());
			FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
			final Path in=kind==1 ? equalInput : input, out=kind==1 ? equalOutput : output;
			final Path ref=kind==2 ? reference : kind==3 ? equalReference : null;
			final String[] args=arguments(in, out, ref, kind<2 ? "AC" : "GG");
			Files.write(work.resolve("constructor.argv"), Arrays.asList(args), StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
			final AddAdapters app=new AddAdapters(args);
			final String[] intended=pathOracle(kind, false), expected=pathOracle(kind, baseline), actual=pathValues(app);
			final boolean same=Arrays.equals(actual, intended);
			if(!same){pathMismatches++;}
			write(work.resolve("paths.tsv"), pathTable(intended, expected, actual));
			check(Arrays.equals(actual, expected), label+": parsed paths/descriptors/adapters differ from literal phase oracle.");
			check(mainPathsAbsent(), label+": no full or truncated main path may be opened/created.");
			check(referenceUnchanged(), label+": both private adapter-reference byte strings must remain unchanged.");
			check(!FASTQ.TEST_INTERLEAVED && !FASTQ.FORCE_INTERLEAVED, label+": normal loading must preserve false pairing flags.");
			write(work.resolve("terminal.tsv"), terminalTable());
			append(root.resolve("cases.tsv"), label+"\t"+rowText(intended)+"\t"+rowText(expected)+"\t"+rowText(actual)+"\t"+same+"\n");
			cases++;
		}finally{
			saved.restore();
			write(work.resolve("globals-restored.tsv"), new Globals().table());
			check(saved.matches(), label+": restore every captured process setting after observation.");
		}
	}

	private static String[] arguments(final Path in, final Path out, final Path ref, final String literal){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+in, "out="+out, "mode=write", "literal="+literal,
			"arc=f", "adderrors=f", "simd=f", "t=1", "qin=33", "qout=33", "ow=t", "app=f"));
		if(ref!=null){args.add("ref="+ref);}
		return args.toArray(new String[args.size()]);
	}
	private String[] pathOracle(final int kind, final boolean beforeRepair){
		final Path in=kind==1 && !beforeRepair ? equalInput : input, out=kind==1 && !beforeRepair ? equalOutput : output;
		final Path ref=kind==2 || kind==3 && beforeRepair ? reference : kind==3 ? equalReference : null;
		final String adapters=kind<2 ? "[AC]" : kind==3 && !beforeRepair ? "[TG, GG]" : "[AC, GG]";
		return new String[]{in.toString(), out.toString(), in.toString(), out.toString(), ref==null ? "<null>" : ref.toString(),
			ref==null ? "<null>" : ref.toString(), adapters};
	}
	private static String[] pathValues(final AddAdapters app) throws Exception{
		return new String[]{fieldText(app, "in1"), fieldText(app, "out1"), descriptor(app, "ffin1"), descriptor(app, "ffout1"),
			fieldText(app, "adapterFile"), descriptor(app, "ffa"), adapterText(app)};
	}
	private static String pathTable(final String[] intended, final String[] expected, final String[] actual){
		final StringBuilder text=new StringBuilder("field\tintended\tphaseExpected\tactual\n");
		for(int i=0; i<PATH_FIELDS.length; i++){text.append(PATH_FIELDS[i]).append('\t').append(intended[i]).append('\t').append(expected[i]).append('\t').append(actual[i]).append('\n');}
		return text.toString();
	}
	private static String rowText(final String[] values){return Arrays.toString(values);}
	private static Object value(final Object object, final String name) throws Exception{
		final Field field=object.getClass().getDeclaredField(name); field.setAccessible(true); return field.get(object);
	}
	private static String fieldText(final AddAdapters app, final String name) throws Exception{
		final Object field=value(app, name); return field==null ? "<null>" : field.toString();
	}
	private static String descriptor(final AddAdapters app, final String name) throws Exception{
		final FileFormat ff=(FileFormat)value(app, name); return ff==null ? "<null>" : ff.name();
	}
	private static String adapterText(final AddAdapters app) throws Exception{
		final Object list=value(app, "adapters");
		if(!(list instanceof ArrayList<?>)){return "<null>";}
		final ArrayList<String> text=new ArrayList<String>();
		for(final Object bases : (ArrayList<?>)list){text.add(new String((byte[])bases, StandardCharsets.US_ASCII));}
		return text.toString();
	}

	/** Accepts only the documented synchronous directory failure, never a different runtime error. */
	private static boolean acquisitionFailure(final RuntimeException failure, final String directory){
		if(failure==null || failure.getClass()!=RuntimeException.class || !(failure.getCause() instanceof FileNotFoundException)){return false;}
		final Throwable cause=failure.getCause();
		return cause.getMessage()!=null && cause.getMessage().contains(directory) && frame(cause, "java.io.FileInputStream", "<init>") &&
			frame(failure, "fileIO.ReadWrite", "getRawInputStream") && frame(failure, "stream.FastaReadInputStream", "open") &&
			frame(failure, "jgi.AddAdapters", "makeAdapterList");
	}
	private static boolean frame(final Throwable error, final String owner, final String method){
		for(final StackTraceElement frame : error.getStackTrace()){
			if(frame.getClassName().equals(owner) && frame.getMethodName().equals(method)){return true;}
		}
		return false;
	}
	private boolean mainPathsAbsent(){
		return !Files.exists(input) && !Files.exists(output) && !Files.exists(equalInput) && !Files.exists(equalOutput);
	}
	private boolean referenceUnchanged() throws Exception{
		return Arrays.equals(Files.readAllBytes(reference), REF_AC.getBytes(StandardCharsets.US_ASCII)) &&
			Arrays.equals(Files.readAllBytes(equalReference), REF_TG.getBytes(StandardCharsets.US_ASCII));
	}
	private String terminalTable() throws Exception{
		return "mainPathsAbsent\tadapterDirectoryExists\treferencePrefixUnchanged\treferenceFullUnchanged\n"+
			mainPathsAbsent()+"\t"+Files.isDirectory(badReference)+"\t"+
			Arrays.equals(Files.readAllBytes(reference), REF_AC.getBytes(StandardCharsets.US_ASCII))+"\t"+
			Arrays.equals(Files.readAllBytes(equalReference), REF_TG.getBytes(StandardCharsets.US_ASCII))+"\n";
	}
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}
	private static void append(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);
	}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}
	private void check(final boolean condition, final String message){checks++; if(!condition){failures++; System.err.println("CHECK_FAILED "+message);}}

	/** Explicit bounded snapshot of the21 mutable settings touched by these arguments and constructors. */
	private static final class Globals{
		Globals() throws Exception{
			add(FASTQ.class, "TEST_INTERLEAVED", "FORCE_INTERLEAVED", "PARSE_CUSTOM", "ASCII_OFFSET", "ASCII_OFFSET_OUT",
				"DETECT_QUALITY", "DETECT_QUALITY_OUT", "SET_QIN");
			add(Shared.class, "THREADS", "READ_BUFFER_NUM_BUFFERS", "SIMD", "COMMAND_LINE", "mainClass");
			add(ReadWrite.class, "USE_PIGZ", "USE_UNPIGZ", "MAX_ZIP_THREADS", "SET_ZIP_THREADS");
			add(Parser.class, "qin", "qout", "parsedQuality");
			add(ReadStats.class, "append");
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
			return fields.size()==21;
		}
		String table(){
			final StringBuilder text=new StringBuilder("setting\tvalue\n");
			for(int i=0; i<fields.size(); i++){
				final Object value=values.get(i);
				final String rendered=value instanceof String[] ? Arrays.toString((String[])value) : String.valueOf(value);
				text.append(fields.get(i).getDeclaringClass().getSimpleName()).append('.').append(fields.get(i).getName()).append('\t')
					.append(rendered.replace("\t", "\\t").replace("\n", "\\n").replace("\r", "\\r")).append('\n');
			}
			return text.toString();
		}
		private final ArrayList<Field> fields=new ArrayList<Field>();
		private final ArrayList<Object> values=new ArrayList<Object>();
	}

	private static final String[] PATH_CASES={"plain-paths", "equals-main-paths", "plain-reference", "equals-reference"};
	private static final String[] PATH_FIELDS={"in1", "out1", "ffin1.name", "ffout1.name", "adapterFile", "ffa.name", "adapters"};
	private static final String REF_AC=">a\nAC\n", REF_TG=">b\nTG\n";
	private final Path root, input, output, equalInput, equalOutput, badReference, reference, equalReference;
	private final boolean baseline;
	private int cases, checks, failures, flagMismatches, pathMismatches;
}
