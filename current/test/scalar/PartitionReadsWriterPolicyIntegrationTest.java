package test.scalar;

import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import fileIO.ByteFile;
import shared.Timer;
import stream.MultiWriterPolicy;

/** Integration checks for PartitionReads bounded writer-policy routes. */
public class PartitionReadsWriterPolicyIntegrationTest{

	private static final String FASTQ1="@p0\nACGT\n+\nIIII\n@p1\nTGCA\n+\nJJJJ\n";
	private static final String FASTQ2="@p0\nGGGG\n+\nHHHH\n@p1\nCCCC\n+\nFFFF\n";
	private static final String SAM_TEXT="@HD\tVN:1.6\tSO:unsorted\n"
			+"@SQ\tSN:ref\tLN:100\n"
			+"s0\t0\tref\t1\t60\t4=\t*\t0\t0\tACGT\tIIII\n"
			+"s1\t0\tref\t41\t60\t4=\t*\t0\t0\tTGCA\tJJJJ\n";
	private static final String AUTO_SAM_TEXT="@HD\tVN:1.6\tSO:unsorted\n"
			+"@SQ\tSN:ref\tLN:100\n"
			+"low0\t0\tref\t1\t60\t16=\t*\t0\t0\tAAAACCCCAAAAAAAA\tIIIIIIIIIIIIIIII\n"
			+"low1\t0\tref\t1\t60\t16=\t*\t0\t0\tAAAACCCCAAAAAAAA\tIIIIIIIIIIIIIIII\n"
			+"low2\t0\tref\t1\t60\t16=\t*\t0\t0\tAAAACCCCAAAAAAAA\tIIIIIIIIIIIIIIII\n"
			+"low3\t0\tref\t1\t60\t16=\t*\t0\t0\tAAAACCCCAAAAAAAA\tIIIIIIIIIIIIIIII\n"
			+"high0\t0\tref\t41\t60\t16=\t*\t0\t0\tGGGGCCCCGGGGAAAA\tIIIIIIIIIIIIIIII\n"
			+"high1\t0\tref\t41\t60\t16=\t*\t0\t0\tGGGGCCCCGGGGAAAA\tIIIIIIIIIIIIIIII\n"
			+"high2\t0\tref\t41\t60\t16=\t*\t0\t0\tGGGGCCCCGGGGAAAA\tIIIIIIIIIIIIIIII\n"
			+"high3\t0\tref\t41\t60\t16=\t*\t0\t0\tGGGGCCCCGGGGAAAA\tIIIIIIIIIIIIIIII\n";

	/** Runs all tiny writer-policy integration cases for one partitioner variant. */
	public static void main(String[] args) throws Exception{
		final String variant=args[0];
		final Path dir=Paths.get(args[1]);
		Files.createDirectories(dir);

		requireMode(constructFastq(variant, dir, "forced_st", "multiwritermode=st",
				"ways=3"), MultiWriterPolicy.ST, "forced ST");
		requireMode(constructFastq(variant, dir, "forced_zt", "multiwritermode=zt",
				"ways=3"), MultiWriterPolicy.ZT, "forced ZT");
		requireMode(constructFastq(variant, dir, "auto_fastq", "multiwritermode=auto",
				"multiwriterstmax=4", "multiwriterstmemory=1", "t=5", "ways=2"),
				MultiWriterPolicy.ST, "AUTO paired FASTQ");
		requireMode(constructFastaQual(variant, dir, "auto_fastaqual", "t=5"),
				MultiWriterPolicy.ZT, "AUTO FASTA+QUAL physical output count");

		runPairedFastq(variant, dir.resolve("paired"));
		runFastaQual(variant, dir.resolve("fasta_qual"));
		if(!variant.equals("1")){runAutoSam(variant, dir.resolve("auto_sam"));}
		runSam(variant, dir.resolve("sam"), ".sam");
		runSam(variant, dir.resolve("bam"), ".bam");
		runZip(variant, dir.resolve("zip"));
		runUnsupported(variant, dir.resolve("unsupported"));

		System.out.println("PASS "+variant+" writer-policy-integration");
	}

	/** Constructs a paired FASTQ case and snapshots its selected writer policy. */
	private static Object constructFastq(String variant, Path dir, String label,
			String mode, String... extra) throws Exception{
		Files.createDirectories(dir);
		final Path in1=dir.resolve(label+"_1.fq");
		final Path in2=dir.resolve(label+"_2.fq");
		write(in1, FASTQ1);
		write(in2, FASTQ2);
		final ArrayList<String> args=new ArrayList<String>();
		args.addAll(Arrays.asList("in="+in1, "in2="+in2,
				"out="+dir.resolve(label+"_%_1.fq"), "out2="+dir.resolve(label+"_%_2.fq"),
				mode, "ow=t"));
		args.addAll(Arrays.asList(extra));
		return construct(variant, args);
	}

	/** Constructs a paired FASTA+QUAL case whose eight physical outputs force AUTO ZT at t=5. */
	private static Object constructFastaQual(String variant, Path dir, String label,
			String threadsArg) throws Exception{
		Files.createDirectories(dir);
		final Path in1=dir.resolve(label+"_1.fq");
		final Path in2=dir.resolve(label+"_2.fq");
		write(in1, FASTQ1);
		write(in2, FASTQ2);
		return construct(variant, Arrays.asList("in="+in1, "in2="+in2,
				"out="+dir.resolve(label+"_%_1.fa"),
				"out2="+dir.resolve(label+"_%_2.fa"),
				"qfout="+dir.resolve(label+"_%_1.qual"),
				"qfout2="+dir.resolve(label+"_%_2.qual"),
				"multiwritermode=auto", "multiwriterstmax=8", "multiwriterstmemory=1",
				threadsArg, "ways=2", "ow=t"));
	}

	/** Confirms twin FASTQ primary and mate records land in the matching shards. */
	private static void runPairedFastq(String variant, Path dir) throws Exception{
		final Object partition=constructFastq(variant, dir, "pairs", "multiwritermode=zt",
				"ways=2");
		process(partition);
		requireFastq(dir.resolve("pairs_0_1.fq"), rec("p0", "ACGT", "IIII"), "R1 part0");
		requireFastq(dir.resolve("pairs_0_2.fq"), rec("p0", "GGGG", "HHHH"), "R2 part0");
		requireFastq(dir.resolve("pairs_1_1.fq"), rec("p1", "TGCA", "JJJJ"), "R1 part1");
		requireFastq(dir.resolve("pairs_1_2.fq"), rec("p1", "CCCC", "FFFF"), "R2 part1");
	}

	/** Confirms FASTA+QUAL output writes both bases and numeric qualities per shard. */
	private static void runFastaQual(String variant, Path dir) throws Exception{
		final Object partition=constructFastaQual(variant, dir, "fq", "t=9");
		requireMode(partition, MultiWriterPolicy.ST, "AUTO FASTA+QUAL ST");
		process(partition);
		requireFasta(dir.resolve("fq_0_1.fa"), rec("p0", "ACGT", ""), "FASTA R1 part0");
		requireFasta(dir.resolve("fq_0_2.fa"), rec("p0", "GGGG", ""), "FASTA R2 part0");
		requireQual(dir.resolve("fq_0_1.qual"), rec("p0", "", "40 40 40 40"), "QUAL R1 part0");
		requireQual(dir.resolve("fq_0_2.qual"), rec("p0", "", "39 39 39 39"), "QUAL R2 part0");
		requireFasta(dir.resolve("fq_1_1.fa"), rec("p1", "TGCA", ""), "FASTA R1 part1");
		requireFasta(dir.resolve("fq_1_2.fa"), rec("p1", "CCCC", ""), "FASTA R2 part1");
		requireQual(dir.resolve("fq_1_1.qual"), rec("p1", "", "41 41 41 41"), "QUAL R1 part1");
		requireQual(dir.resolve("fq_1_2.qual"), rec("p1", "", "37 37 37 37"), "QUAL R2 part1");
	}

	/** Confirms SAM or native BAM output can be decoded with exact mapped records. */
	private static void runSam(String variant, Path dir, String extension) throws Exception{
		Files.createDirectories(dir);
		final Path input=dir.resolve("reads.sam");
		write(input, SAM_TEXT);
		final Object partition=construct(variant, Arrays.asList("in="+input,
				"out="+dir.resolve("part%"+extension), "multiwritermode=zt",
				"ways=2", "ow=t"));
		process(partition);
		requireSamRecords(readSamLines(dir.resolve("part0"+extension)),
				sam("s0", 1, "ACGT", "IIII"));
		requireSamRecords(readSamLines(dir.resolve("part1"+extension)),
				sam("s1", 41, "TGCA", "JJJJ"));
	}

	/** Confirms ZIP-compressed FASTQ shards decode to their assigned records. */
	private static void runZip(String variant, Path dir) throws Exception{
		final Object partition=constructFastq(variant, dir, "zip", "multiwritermode=zt",
				"ways=2", "out="+dir.resolve("zip%_1.fq.zip"),
				"out2="+dir.resolve("zip%_2.fq.zip"));
		process(partition);
		requireFastq(zipText(dir.resolve("zip0_1.fq.zip")), rec("p0", "ACGT", "IIII"), "ZIP R1 part0");
		requireFastq(zipText(dir.resolve("zip0_2.fq.zip")), rec("p0", "GGGG", "HHHH"), "ZIP R2 part0");
		requireFastq(zipText(dir.resolve("zip1_1.fq.zip")), rec("p1", "TGCA", "JJJJ"), "ZIP R1 part1");
		requireFastq(zipText(dir.resolve("zip1_2.fq.zip")), rec("p1", "CCCC", "FFFF"), "ZIP R2 part1");
	}

	/** Confirms unsupported output compression is rejected before processing starts. */
	private static void runUnsupported(String variant, Path dir) throws Exception{
		Files.createDirectories(dir);
		final Path in1=dir.resolve("bad_1.fq");
		final Path in2=dir.resolve("bad_2.fq");
		write(in1, FASTQ1);
		write(in2, FASTQ2);
		requireConstructUnsupported(variant, Arrays.asList("in="+in1, "in2="+in2,
				"out="+dir.resolve("bad%.fq.bz2"), "multiwritermode=zt", "ways=2", "ow=t"),
				"bz2 FASTQ output");
		requireConstructUnsupported(variant, Arrays.asList("in="+in1, "in2="+in2,
				"out="+dir.resolve("bad%_1.fa"), "out2="+dir.resolve("bad%_2.fa"),
				"qfout="+dir.resolve("bad%_1.qual.bz2"), "multiwritermode=zt",
				"ways=2", "ow=t"), "bz2 QUAL output");
		if(!variant.equals("1")){
			requireConstructUnsupported(variant, Arrays.asList("in="+in1, "in2="+in2,
					"out="+dir.resolve("bad_auto%.fq.bz2"), "mode=gc", "auto=t",
					"maxp=2", "multiwritermode=auto", "ow=t"),
					"auto-partition bz2 FASTQ output");
		}
		requireConstructAllowed(variant, Arrays.asList("in="+in1, "in2="+in2,
				"qfout="+dir.resolve("ignored%.qual.bz2"), "multiwritermode=zt",
				"ways=2", "ow=t"), "qfout-only no-output control");
	}

	/** Requires no-output construction to preserve ignored QUAL template compatibility. */
	private static void requireConstructAllowed(String variant, List<String> args,
			String label) throws Exception{
		try{
			construct(variant, args);
		}catch(IllegalArgumentException unexpected){
			throw new AssertionError(label+" was rejected", unexpected);
		}
	}

	/** Requires deterministic lightweight rejection during construction itself. */
	private static void requireConstructUnsupported(String variant, List<String> args,
			String label) throws Exception{
		try{
			construct(variant, args);
		}catch(IllegalArgumentException expected){
			require(causeContains(expected, "No lightweight compression backend"),
					"unexpected unsupported-compression failure: "+expected);
			return;
		}
		throw new AssertionError(label+" was accepted");
	}

	/** Confirms the auto-partition SAM path creates headers and records after detecting ways. */
	private static void runAutoSam(String variant, Path dir) throws Exception{
		Files.createDirectories(dir);
		final Path input=dir.resolve("reads.sam");
		write(input, AUTO_SAM_TEXT);
		final Object partition=construct(variant, Arrays.asList("in="+input,
				"out="+dir.resolve("auto%.sam"), "mode=gc", "auto=t", "maxp=2",
				"multiwritermode=auto", "multiwriterstmax=2", "multiwriterstmemory=1",
				"t=3", "ow=t"));
		process(partition);
		requireMode(partition, MultiWriterPolicy.ST, "AUTO SAM after partition detection");
		requireSamRecords(readSamLines(dir.resolve("auto0.sam")),
				sam("low0", 1, "AAAACCCCAAAAAAAA", "IIIIIIIIIIIIIIII"),
				sam("low1", 1, "AAAACCCCAAAAAAAA", "IIIIIIIIIIIIIIII"),
				sam("low2", 1, "AAAACCCCAAAAAAAA", "IIIIIIIIIIIIIIII"),
				sam("low3", 1, "AAAACCCCAAAAAAAA", "IIIIIIIIIIIIIIII"));
		requireSamRecords(readSamLines(dir.resolve("auto1.sam")),
				sam("high0", 41, "GGGGCCCCGGGGAAAA", "IIIIIIIIIIIIIIII"),
				sam("high1", 41, "GGGGCCCCGGGGAAAA", "IIIIIIIIIIIIIIII"),
				sam("high2", 41, "GGGGCCCCGGGGAAAA", "IIIIIIIIIIIIIIII"),
				sam("high3", 41, "GGGGCCCCGGGGAAAA", "IIIIIIIIIIIIIIII"));
	}

	/** Reflectively constructs the selected historical or current partitioner. */
	private static Object construct(String variant, List<String> args) throws Exception{
		final Class<?> clazz=Class.forName(className(variant));
		final Constructor<?> c=clazz.getConstructor(String[].class);
		try{
			return c.newInstance((Object)args.toArray(new String[args.size()]));
		}catch(InvocationTargetException e){
			throw rethrow(e.getCause());
		}
	}

	/** Maps fixture variant numbers to production classes. */
	private static String className(String variant){
		if(variant.equals("1")){return "scalar.PartitionReads";}
		if(variant.equals("2")){return "scalar.PartitionReads2";}
		if(variant.equals("3")){return "scalar.PartitionReads3";}
		throw new IllegalArgumentException("Unknown PartitionReads variant: "+variant);
	}

	/** Invokes package-private processing without changing production visibility. */
	private static void process(Object target) throws Exception{
		final Method m=target.getClass().getDeclaredMethod("process", Timer.class);
		m.setAccessible(true);
		try{
			m.invoke(target, new Timer());
		}catch(InvocationTargetException e){
			throw rethrow(e.getCause());
		}
	}

	/** Reads the private writerMode selected by a partitioner instance. */
	private static int writerMode(Object target) throws Exception{
		final Field f=target.getClass().getDeclaredField("writerMode");
		f.setAccessible(true);
		return f.getInt(target);
	}

	/** Requires one selected writer policy. */
	private static void requireMode(Object target, int expected, String label) throws Exception{
		final int actual=writerMode(target);
		require(actual==expected, label+" selected "+actual+" instead of "+expected);
	}

	/** Reads SAM text or ByteFile-decoded BAM and rejects deferred read or close errors. */
	private static List<String> readSamLines(Path path) throws Exception{
		final ArrayList<String> lines=new ArrayList<String>();
		final ByteFile bf=ByteFile.makeByteFile(path.toString(), true);
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				lines.add(new String(line, StandardCharsets.UTF_8));
			}
		}catch(Exception|Error e){
			if(bf.close()){e.addSuppressed(new RuntimeException("Error closing "+path));}
			throw e;
		}
		require(!bf.close(), "Error reading or closing "+path);
		return lines;
	}

	/** Requires one decoded SAM/BAM shard to retain its dictionary and exact records. */
	private static void requireSamRecords(List<String> lines, SamRec... expected){
		require(lines.contains("@SQ\tSN:ref\tLN:100"), "SAM/BAM output lost @SQ ref header: "+lines);
		final ArrayList<String> records=new ArrayList<String>();
		for(String line : lines){
			if(!line.startsWith("@") && line.length()>0){records.add(line);}
		}
		require(records.size()==expected.length,
				"record count changed: "+records.size()+" != "+expected.length+" in "+records);
		for(int i=0; i<expected.length; i++){
			final String[] fields=records.get(i).split("\t", -1);
			final SamRec e=expected[i];
			require(fields.length>=11, "record "+i+" has too few SAM fields: "+records.get(i));
			require(fields[0].equals(e.qname), "record "+i+" qname changed: "+records.get(i));
			require(fields[1].equals("0"), "record "+i+" flag changed: "+records.get(i));
			require(fields[2].equals("ref"), "record "+i+" rname changed: "+records.get(i));
			require(fields[3].equals(Integer.toString(e.pos)), "record "+i+" pos changed: "+records.get(i));
			require(fields[4].equals("60"), "record "+i+" mapq changed: "+records.get(i));
			require(fields[5].equals(e.seq.length()+"="), "record "+i+" cigar changed: "+records.get(i));
			require(fields[6].equals("*"), "record "+i+" mate rname changed: "+records.get(i));
			require(fields[7].equals("0"), "record "+i+" mate pos changed: "+records.get(i));
			require(fields[8].equals("0"), "record "+i+" template length changed: "+records.get(i));
			require(fields[9].equals(e.seq), "record "+i+" sequence changed: "+records.get(i));
			require(fields[10].equals(e.qual), "record "+i+" quality changed: "+records.get(i));
		}
	}

	/** Reads and decodes a single-entry ZIP shard. */
	private static String zipText(Path path) throws Exception{
		try(InputStream raw=Files.newInputStream(path);
				ZipInputStream zip=new ZipInputStream(raw);
				ByteArrayOutputStream bytes=new ByteArrayOutputStream()){
			final ZipEntry entry=zip.getNextEntry();
			require(entry!=null, "zip shard had no entry: "+path);
			final byte[] buffer=new byte[1024];
			for(int len=zip.read(buffer); len>=0; len=zip.read(buffer)){bytes.write(buffer, 0, len);}
			require(zip.getNextEntry()==null, "zip shard had multiple entries: "+path);
			return new String(bytes.toByteArray(), StandardCharsets.UTF_8);
		}
	}

	/** Requires a file to contain one exact FASTQ record. */
	private static void requireFastq(Path path, Rec expected, String label) throws Exception{
		requireFastq(read(path), expected, label);
	}

	/** Requires text to contain one exact FASTQ record. */
	private static void requireFastq(String text, Rec expected, String label){
		final String[] lines=text.trim().split("\n");
		require(lines.length==4, label+" FASTQ line count changed: "+lines.length);
		require(lines[0].equals("@"+expected.id), label+" id changed: "+lines[0]);
		require(lines[1].equals(expected.seq), label+" sequence changed: "+lines[1]);
		require(lines[2].equals("+"), label+" plus line changed: "+lines[2]);
		require(lines[3].equals(expected.qual), label+" quality changed: "+lines[3]);
	}

	/** Requires a FASTA file to contain one exact sequence. */
	private static void requireFasta(Path path, Rec expected, String label) throws Exception{
		final String[] lines=read(path).trim().split("\n");
		require(lines.length==2, label+" FASTA line count changed: "+lines.length);
		require(lines[0].equals(">"+expected.id), label+" id changed: "+lines[0]);
		require(lines[1].equals(expected.seq), label+" sequence changed: "+lines[1]);
	}

	/** Requires a QUAL file to contain one exact numeric quality vector. */
	private static void requireQual(Path path, Rec expected, String label) throws Exception{
		final String[] lines=read(path).trim().split("\n");
		require(lines.length==2, label+" QUAL line count changed: "+lines.length);
		require(lines[0].equals(">"+expected.id), label+" id changed: "+lines[0]);
		require(lines[1].equals(expected.qual), label+" quality changed: "+lines[1]);
	}

	/** Reads UTF-8 text. */
	private static String read(Path path) throws Exception{
		return new String(Files.readAllBytes(path), StandardCharsets.UTF_8);
	}

	/** Writes UTF-8 text. */
	private static void write(Path path, String text) throws Exception{
		Files.write(path, text.getBytes(StandardCharsets.UTF_8));
	}

	/** Builds an expected decoded record. */
	private static Rec rec(String id, String seq, String qual){return new Rec(id, seq, qual);}

	/** Builds an expected decoded SAM/BAM record. */
	private static SamRec sam(String qname, int pos, String seq, String qual){
		return new SamRec(qname, pos, seq, qual);
	}

	/** Unwraps reflective target exceptions. */
	private static Exception rethrow(Throwable cause) throws Exception{
		if(cause instanceof Error){throw (Error)cause;}
		if(cause instanceof Exception){throw (Exception)cause;}
		throw new RuntimeException(cause);
	}

	/** Returns whether a throwable chain mentions the expected text. */
	private static boolean causeContains(Throwable t, String text){
		for(Throwable x=t; x!=null; x=x.getCause()){
			if(String.valueOf(x.getMessage()).contains(text)){return true;}
		}
		return false;
	}

	/** These checks remain active under both -ea and -da. */
	private static void require(boolean valid, String message){
		if(!valid){throw new AssertionError(message);}
	}

	/** Expected decoded sequence record. */
	private static final class Rec{
		/** Constructs an expected decoded sequence record. */
		Rec(String id_, String seq_, String qual_){
			id=id_;
			seq=seq_;
			qual=qual_;
		}

		/** Read identifier without leading FASTA/FASTQ sigil. */
		final String id;
		/** Base string; empty for QUAL-only expectations. */
		final String seq;
		/** ASCII or numeric quality string; empty for FASTA-only expectations. */
		final String qual;
	}

	/** Expected mapped SAM/BAM alignment fields. */
	private static final class SamRec{
		/** Constructs expected mapped SAM/BAM alignment fields. */
		SamRec(String qname_, int pos_, String seq_, String qual_){
			qname=qname_;
			pos=pos_;
			seq=seq_;
			qual=qual_;
		}

		/** Expected QNAME. */
		final String qname;
		/** Expected 1-based POS. */
		final int pos;
		/** Expected SEQ. */
		final String seq;
		/** Expected QUAL text. */
		final String qual;
	}
}
