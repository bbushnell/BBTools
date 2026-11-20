package test.jgi;

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
import java.util.Objects;

import fileIO.ReadWrite;
import jgi.AdjustHomopolymers;
import parse.Parse;
import parse.Parser;
import shared.Shared;
import stream.FASTQ;
import stream.Read;
import tracker.ReadStats;

/**
 * Eight finite direct-method baseline/fixed scenarios for homopolymer run editing.
 * Literal sequence and Phred-byte quality oracles cover boundaries, final runs,
 * null qualities, unchanged controls and a linked pair. A private valid input is
 * used only for constructor preflight; process and all read/write streams are excluded.
 * Read constructor validation is disabled through its existing overload so it cannot
 * normalize the literal N qualities before the editing helper is tested.
 * @author Sucrose
 */
public final class AdjustHomopolymersTest{

	/** Requires phase=baseline|fixed and a new outdir; no extreme rates are accepted by the matrix. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown direct-method fixture argument: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null, "Set a new outdir for literal oracles and observed values.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		final String name=root.toString();
		require(name.indexOf('=')<0 && name.indexOf('\t')<0 && name.indexOf('\n')<0 && name.indexOf('\r')<0,
			"Current constructor splits equals signs; evidence paths must fit one row.");
		final Method single=AdjustHomopolymers.class.getDeclaredMethod("processRead", Read.class);
		final Method paired=AdjustHomopolymers.class.getDeclaredMethod("processReadPair", Read.class, Read.class);
		single.setAccessible(true); paired.setAccessible(true);
		Files.createDirectory(root);
		final AdjustHomopolymersTest test=new AdjustHomopolymersTest(root, phase.equals("baseline"), single, paired);
		final Globals inherited=new Globals();
		try{for(final Scenario row : SCENARIOS){test.scenario(row);}}
		finally{inherited.restore();}
		test.check(inherited.matches(), "Restore all22 captured process settings at fixture exit.");
		test.check(test.cases==8, "Frozen matrix requires8 scenarios, observed="+test.cases);
		test.check(test.readCalls==9, "Seven individual reads and one linked pair require9 editing calls: "+test.readCalls);
		test.check(test.pairCalls==1, "Exactly one scenario must exercise processReadPair: "+test.pairCalls);
		test.check(test.mismatchCases==(test.baseline ? 6 : 0), "Expected6 baseline mismatch scenarios or0 fixed: "+test.mismatchCases);
		if(test.failures>0){throw new AssertionError("AdjustHomopolymers failed "+test.failures+" of "+test.checks+" checks; inspect "+root);}
		System.out.println("ADJUST_HOMOPOLYMERS_PASS phase="+phase+" cases="+test.cases+" readCalls="+test.readCalls+
			" pairCalls="+test.pairCalls+" mismatchCases="+test.mismatchCases+" checks="+test.checks);
	}

	private AdjustHomopolymersTest(final Path root_, final boolean baseline_, final Method single_, final Method paired_) throws Exception{
		root=root_; baseline=baseline_; single=single_; paired=paired_;
		write(root.resolve("cases.tsv"), "case\trate\treads\tmatchesIntended\n");
	}

	/** Fresh instance per scenario; completed method calls precede terminal input/global observations. */
	private void scenario(final Scenario row) throws Exception{
		final Path work=Files.createDirectory(root.resolve(row.name));
		final Path input=work.resolve("input.fq"), output=work.resolve("unused-output.fq");
		write(input, INPUT); write(work.resolve("input-before.fq"), INPUT);
		final Globals saved=new Globals();
		try{
			write(work.resolve("globals-before.tsv"), saved.table());
			final String[] args={"in="+input, "out="+output, "rate="+row.rate, "int=f", "t=1", "simd=f",
				"qin=33", "qout=33", "ow=t", "app=f"};
			Files.write(work.resolve("constructor.argv"), Arrays.asList(args), StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
			final AdjustHomopolymers app=new AdjustHomopolymers(args);
			final Read[] reads=new Read[row.reads.length];
			for(int i=0; i<reads.length; i++){
				final ReadCase expected=row.reads[i];
				reads[i]=new Read(bytes(expected.input), copy(expected.inputQual), row.name+"/"+(i+1), 700+cases, false);
			}
			if(reads.length==2){
				reads[0].mate=reads[1]; reads[1].mate=reads[0];
				reads[0].setPairnum(0); reads[1].setPairnum(1);
			}
			final Metadata[] before=new Metadata[reads.length];
			for(int i=0; i<reads.length; i++){before[i]=new Metadata(reads[i]);}
			write(work.resolve("reads-before.tsv"), snapshots(before));
			if(reads.length==2){
				final Object retained=invoke(paired, app, reads[0], reads[1]); pairCalls++; readCalls+=2;
				write(work.resolve("pair-result.txt"), String.valueOf(retained)+"\n");
				check(Boolean.TRUE.equals(retained), row.name+": paired wrapper must retain the pair.");
			}else{
				invoke(single, app, reads[0]); readCalls++;
			}
			boolean intended=true;
			final StringBuilder observed=new StringBuilder("read\tinputBases\tinputQualities\tbaselineBases\tbaselineQualities\tintendedBases\tintendedQualities\tphaseExpectedBases\tphaseExpectedQualities\tactualBases\tactualQualities\n");
			final Metadata[] after=new Metadata[reads.length];
			for(int i=0; i<reads.length; i++){
				final Read read=reads[i]; final ReadCase expected=row.reads[i];
				final String expectedBases=baseline ? expected.oldBases : expected.newBases;
				final byte[] expectedQual=baseline ? expected.oldQual : expected.newQual;
				check(Arrays.equals(read.bases, bytes(expectedBases)), row.name+"/"+(i+1)+": sequence differs from literal phase oracle.");
				check(Arrays.equals(read.quality, expectedQual), row.name+"/"+(i+1)+": Phred qualities/nullness differ from literal phase oracle.");
				check(before[i].unchanged(), row.name+"/"+(i+1)+": editing must preserve borrowed metadata and mate links.");
				intended&=Arrays.equals(read.bases, bytes(expected.newBases)) && Arrays.equals(read.quality, expected.newQual);
				observed.append(i+1).append('\t').append(expected.input).append('\t').append(Arrays.toString(expected.inputQual)).append('\t')
					.append(expected.oldBases).append('\t').append(Arrays.toString(expected.oldQual)).append('\t')
					.append(expected.newBases).append('\t').append(Arrays.toString(expected.newQual)).append('\t')
					.append(expectedBases).append('\t').append(Arrays.toString(expectedQual)).append('\t')
					.append(new String(read.bases, StandardCharsets.US_ASCII)).append('\t').append(Arrays.toString(read.quality)).append('\n');
				after[i]=new Metadata(read);
			}
			if(!intended){mismatchCases++;}
			write(work.resolve("observed.tsv"), observed.toString());
			write(work.resolve("reads-after.tsv"), snapshots(after));
			Files.write(work.resolve("input-after.fq"), Files.readAllBytes(input), StandardOpenOption.CREATE_NEW);
			check(matches(input, INPUT), row.name+": constructor/direct editing changed its private preflight input.");
			check(!Files.exists(output), row.name+": direct-method fixture must not create a main data output.");
			append(root.resolve("cases.tsv"), row.name+"\t"+row.rate+"\t"+reads.length+"\t"+intended+"\n");
			cases++;
		}finally{
			saved.restore();
			write(work.resolve("globals-restored.tsv"), new Globals().table());
			check(saved.matches(), row.name+": restore all captured process settings after observation.");
		}
	}

	/** Unwraps reflection without accepting a target failure as successful test behavior. */
	private static Object invoke(final Method method, final Object target, final Object... args) throws Exception{
		try{return method.invoke(target, args);}
		catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof Exception){throw (Exception)cause;}
			if(cause instanceof Error){throw (Error)cause;}
			throw e;
		}
	}
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static byte[] copy(final byte[] values){return values==null ? null : values.clone();}
	private static byte[] q(final int... values){
		final byte[] out=new byte[values.length];
		for(int i=0; i<values.length; i++){out[i]=(byte)values[i];}
		return out;
	}
	private static boolean matches(final Path path, final String text) throws Exception{
		return Files.isRegularFile(path) && Arrays.equals(Files.readAllBytes(path), bytes(text));
	}
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, bytes(text), StandardOpenOption.CREATE_NEW);
	}
	private static void append(final Path path, final String text) throws Exception{
		Files.write(path, bytes(text), StandardOpenOption.APPEND);
	}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}
	private void check(final boolean condition, final String message){checks++; if(!condition){failures++; System.err.println("CHECK_FAILED "+message);}}

	/** Copies metadata values, allowing legitimate replacement/resizing of sequence and quality arrays. */
	private static final class Metadata{
		Metadata(final Read read_){
			read=read_; id=read.id; numericID=read.numericID; mate=read.mate;
			flags=read.flags; chrom=read.chrom; start=read.start; stop=read.stop; insert=read.insert();
		}
		boolean unchanged(){
			return read.id==id && read.numericID==numericID && read.mate==mate && read.flags==flags &&
				read.chrom==chrom && read.start==start && read.stop==stop && read.insert()==insert;
		}
		String row(){return id+"\t"+numericID+"\t"+flags+"\t"+chrom+"\t"+start+"\t"+stop+"\t"+insert+"\t"+(mate==null ? "none" : mate.id)+"\n";}
		final Read read, mate;
		final String id;
		final long numericID;
		final int flags, chrom, start, stop, insert;
	}
	private static String snapshots(final Metadata[] values){
		final StringBuilder out=new StringBuilder("id\tnumericID\tflags\tchrom\tstart\tstop\tinsert\tmateId\n");
		for(final Metadata value : values){out.append(value.row());}
		return out.toString();
	}

	/** Explicit snapshot of the22 named mutable settings touched by these constructors/arguments. */
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
				final String text=value instanceof String[] ? Arrays.toString((String[])value) : String.valueOf(value);
				out.append(fields.get(i).getDeclaringClass().getSimpleName()).append('.').append(fields.get(i).getName()).append('\t')
					.append(text.replace("\t", "\\t").replace("\n", "\\n").replace("\r", "\\r")).append('\n');
			}
			return out.toString();
		}
		private final ArrayList<Field> fields=new ArrayList<Field>();
		private final ArrayList<Object> values=new ArrayList<Object>();
	}

	private static final class ReadCase{
		ReadCase(final String input_, final byte[] inputQual_, final String oldBases_, final byte[] oldQual_, final String newBases_, final byte[] newQual_){
			input=input_; inputQual=inputQual_; oldBases=oldBases_; oldQual=oldQual_; newBases=newBases_; newQual=newQual_;
		}
		final String input, oldBases, newBases;
		final byte[] inputQual, oldQual, newQual;
	}
	private static final class Scenario{
		Scenario(final String name_, final String rate_, final ReadCase... reads_){name=name_; rate=rate_; reads=reads_;}
		final String name, rate;
		final ReadCase[] reads;
	}

	private static final Scenario[] SCENARIOS={
		new Scenario("expand-boundary", "0.5", new ReadCase("AAC", q(10, 20, 30), "AACA", q(10, 20, 30, 20), "AAAC", q(10, 20, 20, 30))),
		new Scenario("shrink-boundary", "-0.5", new ReadCase("AAC", q(10, 20, 30), "AA", q(10, 20), "AC", q(10, 30))),
		new Scenario("expand-final", "0.5", new ReadCase("AAAA", q(10, 20, 30, 40), "AAAA", q(10, 20, 30, 40), "AAAAAA", q(10, 20, 30, 40, 40, 40))),
		new Scenario("shrink-final", "-0.5", new ReadCase("CAAAA", q(5, 10, 20, 30, 40), "CAAAA", q(5, 10, 20, 30, 40), "CAA", q(5, 10, 20))),
		new Scenario("singleton-control", "0.5", new ReadCase("ACGT", q(10, 20, 30, 40), "ACGT", q(10, 20, 30, 40), "ACGT", q(10, 20, 30, 40))),
		new Scenario("undefined-control", "0.5", new ReadCase("NNNN", q(7, 8, 9, 10), "NNNN", q(7, 8, 9, 10), "NNNN", q(7, 8, 9, 10))),
		new Scenario("null-quality", "0.5", new ReadCase("AAAA", null, "AAAA", null, "AAAAAA", null)),
		new Scenario("paired-shrink", "-0.5", new ReadCase("AAC", q(10, 20, 30), "AA", q(10, 20), "AC", q(10, 30)),
			new ReadCase("CAAAA", q(5, 10, 20, 30, 40), "CAAAA", q(5, 10, 20, 30, 40), "CAA", q(5, 10, 20)))
	};
	private static final String INPUT="@constructor-control\nACGT\n+\nIIII\n";
	private final Path root;
	private final boolean baseline;
	private final Method single, paired;
	private int cases, readCalls, pairCalls, checks, failures, mismatchCases;
}
