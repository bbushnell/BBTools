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

import jgi.AddAdapters;
import parse.Parse;
import stream.FASTQ;
import stream.Read;

/**
 * Five direct-grader scenarios with literal baseline and intended counters.
 * Calls only the actual private grade(Read), six times in total; never starts
 * read/write, FASTQ parsing or CRIS. Each scenario constructs a fresh app using
 * generated private input and a literal adapter. Modern controls have explicit
 * valid metadata; legacy reads retain ordinary unknown coordinates/insert.
 * No decoder implementation, insert-aware generation or end-to-end claim.
 * @author Sucrose
 */
public final class AddAdaptersGradingTest{

	/** Requires {@code phase=baseline|fixed} and a new {@code outdir=...}. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else{throw new IllegalArgumentException("Unknown direct-grader argument: "+arg);}
		}
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		require(directory!=null, "Set a new outdir for direct-grader observations.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Direct grade baseline must preserve production assertion behavior; use -ea.");
		final Path root=Paths.get(directory).toAbsolutePath().normalize();
		final String path=root.toString();
		require(path.indexOf('=')<0 && path.indexOf('\t')<0 && path.indexOf('\n')<0 && path.indexOf('\r')<0,
			"Constructor paths must survive AddAdapters argument parsing and fit one table row.");
		Files.createDirectory(root);
		final Method grade=AddAdapters.class.getDeclaredMethod("grade", Read.class);
		grade.setAccessible(true);
		final Field[] fields=new Field[COUNTERS.length];
		for(int i=0; i<fields.length; i++){
			fields[i]=AddAdapters.class.getDeclaredField(COUNTERS[i]);
			fields[i].setAccessible(true);
		}
		final boolean originalParser=FASTQ.PARSE_CUSTOM;
		try{
			for(int i=0; i<NAMES.length; i++){scenario(root, grade, fields, i, phase.equals("baseline"));}
		}finally{
			FASTQ.PARSE_CUSTOM=originalParser;
		}
		check(FASTQ.PARSE_CUSTOM==originalParser, "Fixture must restore its inherited custom-parser flag.");
		check(cases==5 && gradeCalls==6, "Frozen semantic baseline requires five scenarios/six direct calls: "+cases+"/"+gradeCalls);
		if(failures>0){throw new AssertionError("Direct grading failed "+failures+" of "+checks+" checks; inspect "+root);}
		System.out.println("ADD_ADAPTERS_GRADING_PASS phase="+phase+" cases="+cases+" gradeCalls="+gradeCalls+" checks="+checks);
	}

	/** Retains all oracle/observation rows before checking any counter differences. */
	private static void scenario(final Path root, final Method grade, final Field[] fields,
		final int index, final boolean baseline) throws Exception{
		assert(index>=0 && index<NAMES.length) : "Each scenario has one literal read setup and two independent counter rows: "+index;
		final Path work=Files.createDirectory(root.resolve(NAMES[index]));
		final Path input=work.resolve("input.fq"), output=work.resolve("unused.fq");
		write(input, INPUT); write(work.resolve("input-before.fq"), INPUT);
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+input, "out="+output,
			"mode=grade", "literal=AC", "minlen=1", "parsecustom=f", "int=f", "qin=33", "qout=33",
			"simd=f", "t=1", "ow=t", "app=f", "adderrors=f", "arc=f"));
		Files.write(work.resolve("constructor.argv"), args, StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW);
		final AddAdapters app=new AddAdapters(args.toArray(new String[args.size()]));
		final Read[] reads=reads(index);
		final ReadSnapshot[] before=new ReadSnapshot[reads.length];
		for(int i=0; i<reads.length; i++){before[i]=new ReadSnapshot(reads[i]);}
		write(work.resolve("reads-before.tsv"), snapshots(before));
		final boolean parser=FASTQ.PARSE_CUSTOM;
		final long[] expectedBaseline=BASELINE[index], intended=INTENDED[index];
		write(work.resolve("counters.tsv"), counterHeader()+counterRow("baseline", expectedBaseline)+counterRow("intended", intended));
		for(final Read read : reads){
			try{grade.invoke(app, read);}
			catch(InvocationTargetException e){throw new AssertionError("Actual direct grade failed: "+NAMES[index], e.getCause());}
			gradeCalls++;
			boolean unchanged=true;
			for(final ReadSnapshot snapshot : before){unchanged&=snapshot.unchanged();}
			check(unchanged, NAMES[index]+": grading must preserve both Read objects, mate links, arrays and metadata after each call.");
		}
		final long[] actual=new long[fields.length];
		for(int i=0; i<actual.length; i++){actual[i]=fields[i].getLong(app);}
		append(work.resolve("counters.tsv"), counterRow("actual", actual));
		final ReadSnapshot[] after=new ReadSnapshot[reads.length];
		for(int i=0; i<reads.length; i++){after[i]=new ReadSnapshot(reads[i]);}
		write(work.resolve("reads-after.tsv"), snapshots(after));
		final long[] expected=baseline ? expectedBaseline : intended;
		for(int i=0; i<actual.length; i++){
			check(actual[i]==expected[i], NAMES[index]+"/"+COUNTERS[i]+": expected="+expected[i]+", observed="+actual[i]);
		}
		check(FASTQ.PARSE_CUSTOM==parser, NAMES[index]+": grade must not change parser globals.");
		check(matches(input, INPUT), NAMES[index]+": generated constructor input changed.");
		check(!Files.exists(output), NAMES[index]+": direct grading must not open the unused output.");
		Files.write(work.resolve("input-after.fq"), Files.readAllBytes(input), StandardOpenOption.CREATE_NEW);
		cases++;
	}

	/** Modern metadata is baseline setup only, not a representation for repaired legacy length tags. */
	private static Read[] reads(final int index){
		if(index==0){return new Read[]{read("7_6_2", 2)};}
		if(index==1){return new Read[]{read("7_6_0", 1)};}
		if(index==2){
			final Read first=read("7_6_2_8_5 /1", 2), second=read("7_6_2_8_5 /2", 5);
			first.mate=second; second.mate=first;
			first.setPairnum(0); second.setPairnum(1);
			return new Read[]{first, second};
		}
		assert(index==3 || index==4) : "Only the final two scenarios carry explicit modern metadata: "+index;
		final int insert=index==3 ? 3 : 9;
		final Read modern=read("SYN_7_0_5_"+insert+"_+_100_1_._ref 1:", index==3 ? 3 : 6);
		modern.start=100; modern.stop=105; modern.chrom=1;
		modern.setInsert(insert);
		return new Read[]{modern};
	}
	private static Read read(final String id, final int length){
		final byte[] bases=new byte[length], quality=new byte[length];
		Arrays.fill(bases, (byte)'G'); Arrays.fill(quality, (byte)30);
		return new Read(bases, quality, id, 7);
	}

	/** Snapshots the state a length decoder/grader must borrow without changing, including both linked mates. */
	private static final class ReadSnapshot{
		ReadSnapshot(final Read read_){
			read=read_; bases=read.bases; quality=read.quality;
			basesCopy=bases.clone(); qualityCopy=quality.clone();
			id=read.id; numericID=read.numericID; mate=read.mate;
			start=read.start; stop=read.stop; chrom=read.chrom; flags=read.flags; insert=read.insert();
		}
		boolean unchanged(){
			return read.bases==bases && read.quality==quality && Arrays.equals(read.bases, basesCopy) &&
				Arrays.equals(read.quality, qualityCopy) && read.id==id && read.numericID==numericID &&
				read.mate==mate && read.start==start && read.stop==stop && read.chrom==chrom &&
				read.flags==flags && read.insert()==insert;
		}
		String row(){
			return id+"\t"+numericID+"\t"+basesCopy.length+"\t"+new String(basesCopy, StandardCharsets.US_ASCII)+
				"\t"+Arrays.toString(qualityCopy)+"\t"+start+"\t"+stop+"\t"+chrom+
				"\t"+flags+"\t"+insert+"\t"+(mate==null ? "none" : mate.id)+"\n";
		}
		final Read read, mate;
		final byte[] bases, quality, basesCopy, qualityCopy;
		final String id;
		final long numericID;
		final int start, stop, chrom, flags, insert;
	}
	private static String snapshots(final ReadSnapshot[] snapshots){
		final StringBuilder out=new StringBuilder("id\tnumericID\tlength\tbases\tquality\tstart\tstop\tchrom\tflags\tinsert\tmateId\n");
		for(final ReadSnapshot snapshot : snapshots){out.append(snapshot.row());}
		return out.toString();
	}
	private static String counterHeader(){
		final StringBuilder out=new StringBuilder("row");
		for(final String name : COUNTERS){out.append('\t').append(name);}
		return out.append('\n').toString();
	}
	private static String counterRow(final String label, final long[] values){
		final StringBuilder out=new StringBuilder(label);
		for(final long value : values){out.append('\t').append(value);}
		return out.append('\n').toString();
	}
	/** Literal (column,value) pairs leave all unlisted counters at zero. */
	private static long[] counters(final long... pairs){
		assert((pairs.length&1)==0) : "Counter oracles are named column/value pairs, never an odd list.";
		final long[] out=new long[COUNTERS.length];
		for(int i=0; i<pairs.length; i+=2){
			assert(pairs[i]>=0 && pairs[i]<out.length) : "Oracle column must name a reflected grading counter: "+pairs[i];
			out[(int)pairs[i]]=pairs[i+1];
		}
		return out;
	}
	private static boolean matches(final Path path, final String text) throws Exception{
		return Files.isRegularFile(path) && Arrays.equals(Files.readAllBytes(path), text.getBytes(StandardCharsets.US_ASCII));
	}
	private static void write(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}
	private static void append(final Path path, final String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);
	}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}
	private static void check(final boolean condition, final String message){
		checks++;
		if(!condition){failures++; System.err.println("CHECK_FAILED "+message);}
	}

	private static final String[] COUNTERS={"readsProcessed", "basesProcessed", "correct", "correctBases",
		"validBasesExpected", "validBasesCounted", "adapterReadsTotal", "adapterBasesTotal", "truePos", "trueNeg",
		"falsePos", "falseNeg", "tooShort", "tooShortBases", "tooShortReadBases", "tooLong", "tooLongBases",
		"tooLongReadBases", "adapterReadsRemaining", "adapterBasesRemaining", "broken", "mispaired", "fullyRemoved"};
	private static final int READS=0, BASES=1, CORRECT=2, CORRECT_BASES=3, EXPECTED=4, COUNTED=5,
		ADAPTER_READS=6, ADAPTER_BASES=7, TP=8, TN=9, FN=11, TOO_LONG=15, TOO_LONG_BASES=16,
		TOO_LONG_READ_BASES=17, ADAPTER_READS_LEFT=18, ADAPTER_BASES_LEFT=19, BROKEN=20;
	private static final long[][] BASELINE={
		counters(READS, 1, BASES, 2, EXPECTED, -1, COUNTED, -1, ADAPTER_READS, 1, ADAPTER_BASES, 2,
			FN, 1, TOO_LONG, 1, TOO_LONG_BASES, 3, TOO_LONG_READ_BASES, 2, ADAPTER_READS_LEFT, 1, ADAPTER_BASES_LEFT, 3, BROKEN, 1),
		counters(READS, 1, BASES, 1, CORRECT, 1, CORRECT_BASES, -1, EXPECTED, -1, COUNTED, -1,
			ADAPTER_READS, 1, ADAPTER_BASES, 2, TP, 1),
		counters(READS, 2, BASES, 7, EXPECTED, -2, COUNTED, -2, ADAPTER_READS, 2, ADAPTER_BASES, 4,
			FN, 2, TOO_LONG, 2, TOO_LONG_BASES, 9, TOO_LONG_READ_BASES, 7, ADAPTER_READS_LEFT, 2, ADAPTER_BASES_LEFT, 9, BROKEN, 2),
		counters(READS, 1, BASES, 3, CORRECT, 1, CORRECT_BASES, 3, EXPECTED, 3, COUNTED, 3, ADAPTER_READS, 1, ADAPTER_BASES, 3, TP, 1),
		counters(READS, 1, BASES, 6, CORRECT, 1, CORRECT_BASES, 6, EXPECTED, 6, COUNTED, 6, TN, 1)
	};
	private static final long[][] INTENDED={
		counters(READS, 1, BASES, 2, CORRECT, 1, CORRECT_BASES, 2, EXPECTED, 2, COUNTED, 2, ADAPTER_READS, 1, ADAPTER_BASES, 4, TP, 1),
		counters(READS, 1, BASES, 1, CORRECT, 1, EXPECTED, 0, COUNTED, 0, ADAPTER_READS, 1, ADAPTER_BASES, 6, TP, 1),
		counters(READS, 2, BASES, 7, CORRECT, 2, CORRECT_BASES, 7, EXPECTED, 7, COUNTED, 7, ADAPTER_READS, 2, ADAPTER_BASES, 7, TP, 2),
		counters(READS, 1, BASES, 3, CORRECT, 1, CORRECT_BASES, 3, EXPECTED, 3, COUNTED, 3, ADAPTER_READS, 1, ADAPTER_BASES, 3, TP, 1),
		counters(READS, 1, BASES, 6, CORRECT, 1, CORRECT_BASES, 6, EXPECTED, 6, COUNTED, 6, TN, 1)
	};
	private static final String[] NAMES={"legacy-trimmed", "legacy-zero-placeholder", "legacy-pair", "modern-trimmed", "modern-untrimmed"};
	private static final String INPUT="@ordinary\nG\n+\nI\n";
	private static int checks, failures, cases, gradeCalls;
}
