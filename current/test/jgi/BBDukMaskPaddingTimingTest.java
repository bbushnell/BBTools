package test.jgi;

import java.io.IOException;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import jgi.BBDuk;
import structures.ByteBuilder;

/** Bounded normal-process padding timing observer; Bash times only the sample JVM.
 * Prepare generates immutable input; sample constructs one ordinary owner and calls
 * public process() once; verify checks every literal input/output byte in bounded
 * memory. Preparation and verification are separate JVMs outside sample timing.
 * GNU process CPU includes startup/JIT, normal loading, allocation, masking, I/O,
 * closure and small settings/counter/storage observations. Main-thread CPU would
 * miss BBDuk worker work. No fields/tables/counters/release flags are seeded.
 * Baseline masks2051 positions but counts2052 per read; intended fixed masks3077.
 * Thus timing includes the extra corrected base/quality writes, not just arithmetic.
 * The literal workload is513 forward k4 hits, including a terminal hit, per4102 bases.
 * Root owns execution; R is capped at262144 and bulk scratch at20GiB. Fully-covered,
 * negative-padding and other modes are not repaired or certified by this observer.
 * @author Sucrose
 */
public final class BBDukMaskPaddingTimingTest{

	/** Use exactly mode, phase, reads, input, output and outdir as flag=value arguments. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final int split=arg.indexOf('=');
			require(split>0 && split<arg.length()-1, "Use nonempty flag=value arguments: "+arg);
			require(options.put(arg.substring(0, split), arg.substring(split+1))==null, "Duplicate option: "+arg);
		}
		require(options.size()==6 && options.keySet().containsAll(Arrays.asList("mode", "phase", "reads", "input", "output", "outdir")), "Use exactly mode, phase, reads, input, output and outdir.");
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbdukmaskpaddingtimingtest.sh with assertions enabled.");
		final BBDukMaskPaddingTimingTest test=new BBDukMaskPaddingTimingTest(options);
		if(test.mode.equals("prepare")){test.prepare();}
		else if(test.mode.equals("sample")){test.sample();}
		else{test.verify();}
	}

	private BBDukMaskPaddingTimingTest(final Map<String, String> options) throws Exception{
		mode=options.get("mode"); phase=options.get("phase");
		require(mode.equals("prepare") || mode.equals("sample") || mode.equals("verify"), "mode=prepare|sample|verify.");
		require(phase.equals("baseline") || phase.equals("fixed"), "phase=baseline|fixed.");
		baseline=phase.equals("baseline");
		final String value=options.get("reads");
		require(value.matches("[1-9][0-9]*"), "reads must be a positive decimal integer: "+value);
		reads=Integer.parseInt(value); require(reads<=262144, "Approved finite workload caps reads at262144: "+reads);
		input=bulkPath(options.get("input")); output=bulkPath(options.get("output")); out=path(options.get("outdir"));
		require(!input.equals(output) && !out.startsWith(SCRATCH), "Separate immutable input/output; small records belong outside bulk scratch.");
		inputBases=bytes(repeat("GGACGTGG", 512)+"GGACGT"); inputQuality=bytes(repeat("I", LENGTH));
		outputBases=bytes(repeat(baseline ? "GGANNNNG" : "GNNNNNNG", 512)+(baseline ? "GGANNN" : "GNNNNN"));
		outputQuality=bytes(repeat(baseline ? "III!!!!I" : "I!!!!!!I", 512)+(baseline ? "III!!!" : "I!!!!!"));
		require(inputBases.length==LENGTH && inputQuality.length==LENGTH && outputBases.length==LENGTH && outputQuality.length==LENGTH, "Literal block construction must retain4102 bases/qualities.");
		if(!mode.equals("verify")){Files.createDirectory(out);}
		else{require(Files.isDirectory(out) && Files.isRegularFile(out.resolve("sample.tsv")), "Verify the existing sample records, not a new/unrelated directory.");}
		write(out.resolve(mode+"-args.txt"), arguments(options));
		write(out.resolve(mode+"-scope.txt"), "GNU process CPU=user+system and wall; sample JVM includes startup/JIT, real loading, allocation, masking, I/O and small observations. Preparation and full verification are outside timing. Corrected masking performs3077 instead of2051 in-read base/quality writes per read.\n");
	}

	/** Literal generation through the native buffered writer; no BBDuk execution. */
	private void prepare() throws Exception{
		require(!Files.exists(input), "Never overwrite immutable timing input: "+input);
		budget(expectedBytes());
		final ByteStreamWriter writer=new ByteStreamWriter(input.toString(), false, false, false);
		writer.start(); final ByteBuilder buffer=new ByteBuilder(2*LENGTH+64);
		boolean writerError;
		try{
			for(int i=0; i<reads; i++){
				buffer.clear().append(ID_PREFIX).append(i).append('\n').append(inputBases).append('\n').append('+').append('\n').append(inputQuality).append('\n');
				writer.print(buffer);// Copies bytes; this builder remains caller-owned and reusable.
			}
		}finally{writerError=writer.poisonAndWait();}
		check(!writerError, "Native input writer must finish without error before any timing sample.");
		check(Files.size(input)==expectedBytes(), "Generated input must contain the exact finite byte count.");
		verifyFile(input, inputBases, inputQuality);
		require(checks==reads+4, "Prepare checks include every input record, EOF, close, writer and size; actual="+checks);
		write(out.resolve("prepare.tsv"), "reads\tlength\tbytes\thits_per_read\twindows_per_read\tchecks\n"+reads+"\t"+LENGTH+"\t"+Files.size(input)+"\t513\t4099\t"+checks+"\n");
		finish("prepare", checks);
	}

	/** Exactly one ordinary owner/process call; GNU time includes this complete JVM. */
	private void sample() throws Exception{
		require(Files.isRegularFile(input) && !Files.exists(output), "Sample needs existing immutable input and a new output path.");
		require(Files.size(input)==expectedBytes(), "Wrong finite input length for requested read count."); budget(expectedBytes());
		write(out.resolve("settings.tsv"), String.join("\n", SETTINGS)+"\n");
		write(out.resolve("unrelated-counter-schema.tsv"), String.join("\n", ZERO_COUNTERS)+"\n");
		write(out.resolve("observations.tsv"), "stage\treads_in\tbases_in\treads_out\tbases_out\tstored_kmers\tref_reads\tref_bases\tref_kmers\treads_ktrimmed\tbases_ktrimmed\tunrelated_counts\n");
		write(out.resolve("reflection.tsv"), "operation\tmember\nORDINARY_PUBLIC_CALL\tnew BBDuk(public arguments).process()\nREAD_ONLY\tBBDuk parsed settings, reference scalars and statistics\nCOUNT_OBSERVATION\tbasesKTrimmed via normal xsum/join aggregation; NOT direct private return\nFIELD_WRITES\tNONE\n");
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(SETTINGS));
		args.add("reads="+reads); args.add("trimpad=1"); args.add("literal=ACGT"); args.add("in="+input); args.add("out="+output);
		write(out.resolve("constructor-args.txt"), String.join("\n", args)+"\n");
		final BBDuk owner=new BBDuk(args.toArray(new String[args.size()]));
		verifySettings(owner); observe(owner, "before", false);
		owner.process();
		observe(owner, "after", true);
		final long callerCount=number(owner, "basesKTrimmed").longValue();
		check(callerCount==(long)reads*(baseline ? 2052 : 3077) && number(owner, "readsKTrimmed").longValue()==reads, "Normal worker aggregation must report the phase-specific count for every hit-bearing read; actual="+callerCount);
		check(Files.size(output)==expectedBytes(), "Ordinary output must retain every header and4102-byte base/quality record.");
		require(checks==14, "Six settings groups, six before/after groups, mask count and output size give14 sample checks; actual="+checks);
		write(out.resolve("sample.tsv"), "phase\treads\tlength\tbases\tcaller_count\tinput_bytes\toutput_bytes\tchecks\n"+phase+"\t"+reads+"\t"+LENGTH+"\t"+((long)reads*LENGTH)+"\t"+callerCount+"\t"+Files.size(input)+"\t"+Files.size(output)+"\t"+checks+"\n");
		finish("sample", checks);
	}

	/** Reads all literal FASTQ records after timing, with bounded byte-line storage. */
	private void verify() throws Exception{
		final List<String> sample=Files.readAllLines(out.resolve("sample.tsv"), StandardCharsets.US_ASCII);
		require(sample.size()==2 && sample.get(1).equals(phase+"\t"+reads+"\t"+LENGTH+"\t"+((long)reads*LENGTH)+"\t"+((long)reads*(baseline ? 2052 : 3077))+"\t"+expectedBytes()+"\t"+expectedBytes()+"\t14"), "Verify must match the exact completed sample metadata/phase/read count.");
		verifyFile(input, inputBases, inputQuality); verifyFile(output, outputBases, outputQuality);
		check(Files.size(input)==expectedBytes(), "Immutable input byte count must remain exact after timing.");
		check(Files.size(output)==expectedBytes(), "Complete verified output byte count must remain exact.");
		require(checks==2*reads+6, "Verify checks every input/output record plus EOF/close/byte counts; actual="+checks);
		write(out.resolve("verified.tsv"), "phase\treads\tlength\tinput_records\toutput_records\tmasked_per_read\tcaller_per_read\tchecks\n"+phase+"\t"+reads+"\t"+LENGTH+"\t"+reads+"\t"+reads+"\t"+(baseline ? 2051 : 3077)+"\t"+(baseline ? 2052 : 3077)+"\t"+checks+"\n");
		finish("verify", checks);
	}

	/** Raw ASCII literal oracle, not a general FASTQ parser; exact order and four lines matter. */
	private void verifyFile(final Path file, final byte[] bases, final byte[] quality){
		assert(bases.length==LENGTH && quality.length==LENGTH) : "Literal oracle arrays must match the finite record length.";
		final ByteFile reader=ByteFile.makeByteFile1(file.toString(), false); boolean closeError;
		try{
			for(int i=0; i<reads; i++){
				final byte[] id=reader.nextLine(), actualBases=reader.nextLine(), plus=reader.nextLine(), actualQuality=reader.nextLine();
				check(Arrays.equals(id, bytes(ID_PREFIX+i)) && Arrays.equals(actualBases, bases) && Arrays.equals(plus, PLUS) && Arrays.equals(actualQuality, quality), "Every literal record/ID/base/quality must match; file="+file+" record="+i);
			}
			check(reader.nextLine()==null, "No extra FASTQ record or trailing byte line is permitted: "+file);
		}finally{closeError=reader.close();}
		check(!closeError, "Byte-line verifier must close without an input error: "+file);
	}
	/** Read effective settings after public parsing, without seeding any owner state. */
	private void verifySettings(final BBDuk owner) throws Exception{
		check(flag(owner, "ktrimN") && !flag(owner, "ktrimLeft") && !flag(owner, "ktrimRight") && !flag(owner, "ksplit") && !flag(owner, "kmaskLowercase") && !flag(owner, "kmaskFullyCovered") && number(owner, "trimSymbol").intValue()=='N' && number(owner, "trimPad").intValue()==1, "Use only ordinary kmask=N with the publicly requested padding.");
		check(number(owner, "k").intValue()==4 && number(owner, "mink").intValue()==4 && !flag(owner, "useShortKmers") && !flag(owner, "rcomp") && !flag(owner, "maskMiddle") && number(owner, "hammingDistance").intValue()==0 && number(owner, "editDistance").intValue()==0 && number(owner, "qHammingDistance").intValue()==0 && number(owner, "qSkip").intValue()==1 && number(owner, "minSkip").intValue()==1 && number(owner, "maxSkip").intValue()==1 && number(owner, "speed").intValue()==0, "Requested mink=-1 disables short kmers; BBDuk constructor normalizes effective mink to min(6,k)=4. Use exact forward k4 keys without middle-mask, speed or skipping.");
		check(number(owner, "THREADS").intValue()==1 && number(owner, "maxReads").longValue()==reads && number(owner, "minReadLength").intValue()==1 && number(owner, "maxReadLength").intValue()==Integer.MAX_VALUE && number(owner, "minLenFraction").doubleValue()==0 && !flag(owner, "rename") && !flag(owner, "skipR1") && !flag(owner, "skipR2") && number(owner, "restrictLeft").intValue()==0 && number(owner, "restrictRight").intValue()==0, "Process the finite unpaired workload without renaming, length rejection or restricted matching.");
		check(!flag(owner, "qtrimLeft") && !flag(owner, "qtrimRight") && !flag(owner, "trimClip") && number(owner, "minAvgQuality").doubleValue()==0 && number(owner, "minBaseQuality").intValue()==0 && number(owner, "maxNs").intValue()==-1 && number(owner, "maxNRate").doubleValue()==1 && number(owner, "minConsecutiveBases").intValue()<=0 && number(owner, "minBaseFrequency").doubleValue()==0, "Retain masked N/quality0 output; disable quality and undefined-base filtering/trimming.");
		check(number(owner, "forceTrimLeft").intValue()<=0 && number(owner, "forceTrimRight").intValue()<=0 && number(owner, "forceTrimRight2").intValue()<=0 && number(owner, "forceTrimModulo").intValue()<=0 && number(owner, "trimPolyA").intValue()<=0 && number(owner, "trimPolyGLeft").intValue()<=0 && number(owner, "trimPolyGRight").intValue()<=0 && number(owner, "trimPolyCLeft").intValue()<=0 && number(owner, "trimPolyCRight").intValue()<=0 && number(owner, "filterPolyG").intValue()<=0 && number(owner, "filterPolyC").intValue()<=0 && !flag(owner, "trimByOverlap") && !flag(owner, "swift") && !flag(owner, "ecc") && !flag(owner, "calcEntropy") && !flag(owner, "entropyMask") && !flag(owner, "entropyMark") && number(owner, "entropyTrim").intValue()==0 && !flag(owner, "filterGC"), "No force/polymer/overlap/Swift/ECC/entropy/GC operation may confound mask output or counters.");
		check(field("ref").get(owner)==null && Arrays.equals((String[])field("literal").get(owner), new String[]{"ACGT"}) && ((Object[])field("keySets").get(owner)).length==7 && ((List<?>)field("scaffoldNames").get(owner)).size()==1, "The ordinary constructor allocates seven tables and dummy scaffold0; normal process() must load the single literal as scaffold1.");
		final StringBuilder parsed=new StringBuilder("field\tvalue\n");
		for(final String name : PARSED_FIELDS){final Object value=field(name).get(owner); parsed.append(name).append('\t').append(value instanceof String[] ? Arrays.toString((String[])value) : String.valueOf(value)).append('\n');}
		write(out.resolve("parsed-settings.tsv"), parsed.toString());
	}

	private void observe(final BBDuk owner, final String stage, final boolean processed) throws Exception{
		final long count=processed ? reads : 0, bases=count*LENGTH;
		check(number(owner, "readsIn").longValue()==count && number(owner, "readsOut").longValue()==count && number(owner, "basesIn").longValue()==bases && number(owner, "basesOut").longValue()==bases, "Normal processing must emit every unchanged-length read; stage="+stage);
		final long stored=number(owner, "storedKmers").longValue(), refReads=number(owner, "refReads").longValue(), refBases=number(owner, "refBases").longValue(), refKmers=number(owner, "refKmers").longValue();
		check(stored==(processed ? 1 : 0) && refReads==(processed ? 1 : 0) && refBases==(processed ? 4 : 0) && refKmers==(processed ? 1 : 0), "Only real public-process loading may produce the single reference key; stage="+stage);
		final long[] unrelated=new long[ZERO_COUNTERS.length]; boolean zero=true;
		for(int i=0; i<unrelated.length; i++){unrelated[i]=number(owner, ZERO_COUNTERS[i]).longValue(); zero&=unrelated[i]==0;}
		check(zero && (processed || (number(owner, "readsKTrimmed").longValue()==0 && number(owner, "basesKTrimmed").longValue()==0)), "All unrelated trim/filter counters and initial mask counts must stay zero; stage="+stage+" values="+Arrays.toString(unrelated));
		append(out.resolve("observations.tsv"), stage+"\t"+number(owner, "readsIn")+"\t"+number(owner, "basesIn")+"\t"+number(owner, "readsOut")+"\t"+number(owner, "basesOut")+"\t"+stored+"\t"+refReads+"\t"+refBases+"\t"+refKmers+"\t"+number(owner, "readsKTrimmed")+"\t"+number(owner, "basesKTrimmed")+"\t"+Arrays.toString(unrelated)+"\n");
	}

	/** Flat owned scratch only; prospective writes must fit the approved logical-byte peak. */
	private void budget(final long additional) throws IOException{
		assert(additional>=0) : "A prospective FASTQ write cannot reduce the20GiB scratch budget.";
		long used=0;
		try(DirectoryStream<Path> paths=Files.newDirectoryStream(SCRATCH)){
			for(final Path path : paths){require(!Files.isSymbolicLink(path) && Files.isRegularFile(path), "Timing scratch is a flat owned root with no symlink/subdirectory: "+path); used=Math.addExact(used, Files.size(path));}
		}
		require(Math.addExact(used, additional)<=CAP, "Scratch peak would exceed20GiB; used="+used+" additional="+additional+". Report actual cost and redesign; do not expand silently.");
		require(Files.getFileStore(SCRATCH).getUsableSpace()>additional, "Scratch filesystem lacks space for the finite output.");
	}
	private long expectedBytes(){long result=(long)reads*(2*LENGTH+ID_PREFIX.length()+5); for(int i=0; i<reads; i++){result+=Integer.toString(i).length();} return result;}
	private void finish(final String operation, final int count) throws IOException{final String summary="BBDUK_MASK_PADDING_TIMING_PASS version=1 mode="+operation+" phase="+phase+" reads="+reads+" length="+LENGTH+" checks="+count+" failures=0\n"; write(out.resolve(operation+"-summary.txt"), summary); System.out.print(summary);}
	private static Path bulkPath(final String value) throws IOException{final Path result=path(value); require(result.getParent().equals(SCRATCH) && result.toString().endsWith(".fq") && !Files.isSymbolicLink(result) && result.getParent().toRealPath().equals(SCRATCH), "Bulk FASTQ must be directly inside the actual approved scratch root: "+result); return result;}
	private static Path path(final String value){require(value!=null && !value.isEmpty() && value.indexOf('=')<0 && value.indexOf('\t')<0 && value.indexOf('\n')<0, "Use a nonempty path without equals,tabs,newlines."); return Paths.get(value).toAbsolutePath().normalize();}
	private static String arguments(final Map<String, String> options){final StringBuilder result=new StringBuilder(); for(final Map.Entry<String, String> entry : options.entrySet()){result.append(entry.getKey()).append('=').append(entry.getValue()).append('\n');} return result.toString();}
	private static String repeat(final String value, final int count){assert(count>0 && value.length()>0) : "Literal oracle repetition needs a positive finite block count."; final StringBuilder result=new StringBuilder(value.length()*count); for(int i=0; i<count; i++){result.append(value);} return result.toString();}
	private static Field field(final String name) throws NoSuchFieldException{final Field result=BBDuk.class.getDeclaredField(name); result.setAccessible(true); return result;}
	private static Number number(final BBDuk owner, final String name) throws Exception{return (Number)field(name).get(owner);}
	private static boolean flag(final BBDuk owner, final String name) throws Exception{return field(name).getBoolean(owner);}
	private void check(final boolean value, final String why){checks++; require(value, why);}
	private static void require(final boolean value, final String why){if(!value){throw new IllegalArgumentException(why);}}
	private static byte[] bytes(final String value){return value.getBytes(StandardCharsets.US_ASCII);}
	private static void write(final Path file, final String value) throws IOException{Files.write(file, bytes(value), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path file, final String value) throws IOException{Files.write(file, bytes(value), StandardOpenOption.APPEND);}

	private final String mode, phase;
	private final boolean baseline;
	private final int reads;
	private final Path input, output, out;
	private final byte[] inputBases, inputQuality, outputBases, outputQuality;
	private int checks;
	private static final int LENGTH=4102;
	private static final long CAP=20L*1024*1024*1024;
	private static final Path SCRATCH=Paths.get("/clusterfs/jgi/scratch/gentech/genome_analysis/bbushnell/Shinobu/jgi-mask-padding-timing-20261118");
	private static final String ID_PREFIX="@mask-padding-";
	private static final byte[] PLUS={(byte)'+'};
	private static final String[] SETTINGS={"k=4", "mink=-1", "kmask=N", "maskfullycovered=f", "rcomp=f", "maskmiddle=f", "hdist=0", "edist=0", "qhdist=0", "qskip=1", "rskip=1", "speed=0", "rename=f", "minlen=1", "qtrim=f", "maxns=-1", "maq=0", "t=1", "simd=f", "qin=33", "qout=33", "array=t", "prealloc=f", "ow=f", "app=f"};
	private static final String[] ZERO_COUNTERS={"readsKFiltered", "basesKFiltered", "readsQTrimmed", "basesQTrimmed", "readsFTrimmed", "basesFTrimmed", "readsQFiltered", "basesQFiltered", "readsNFiltered", "basesNFiltered", "readsEFiltered", "basesEFiltered", "readsPolyTrimmed", "basesPolyTrimmed", "readsTrimmedByOverlap", "basesTrimmedByOverlap", "readsTrimmedBySwift", "basesTrimmedBySwift", "badGcReads", "badGcBases", "badHeaderReads", "badHeaderBases"};
	private static final String[] PARSED_FIELDS={"k", "mink", "ktrimN", "ktrimLeft", "ktrimRight", "ksplit", "trimSymbol", "trimPad", "kmaskLowercase", "kmaskFullyCovered", "useShortKmers", "rcomp", "maskMiddle", "hammingDistance", "editDistance", "qHammingDistance", "qSkip", "minSkip", "maxSkip", "speed", "THREADS", "maxReads", "minReadLength", "maxReadLength", "minLenFraction", "rename", "skipR1", "skipR2", "restrictLeft", "restrictRight", "qtrimLeft", "qtrimRight", "trimClip", "minAvgQuality", "minBaseQuality", "maxNs", "maxNRate", "minConsecutiveBases", "minBaseFrequency", "forceTrimLeft", "forceTrimRight", "forceTrimRight2", "forceTrimModulo", "trimPolyA", "trimPolyGLeft", "trimPolyGRight", "trimPolyCLeft", "trimPolyCRight", "filterPolyG", "filterPolyC", "trimByOverlap", "swift", "ecc", "calcEntropy", "entropyMask", "entropyMark", "entropyTrim", "filterGC", "literal"};
}
