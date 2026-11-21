package test.jgi;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import fileIO.ReadWrite;
import jgi.BBDuk;
import kmer.AbstractKmerTable;
import shared.KillSwitch;
import structures.ByteBuilder;
import structures.IntList;

/** Streaming X-mask observer and literal verifier, each in a separate fresh JVM.
 * Producer stdout contains only ordinary BBDuk FASTQ. Small reports go to outdir.
 * No bulk input/output is generated. The immutable input is either read completely
 * or consumed through the declared512-read prefix by normal reads= handling.
 * Verifier checks raw bytes, including line endings, without reparsing X as N.
 * Content mismatches drain through EOF before exit2. The finite readfatal injection
 * deliberately exits74, so the owning controller must stop its producer session.
 * Baseline/fixed support freezes verifier bytes before a separately released repair;
 * the baseline controller does not release large samples or the fixed phase.
 * @author Sucrose
 */
public final class BBDukCustomMaskPaddingTimingTest{

	/** mode=producer|verify|input phase=baseline|fixed reads=512|262144 fault=NAME outdir=NEW. */
	public static void main(final String[] args){
		try{
			final Map<String, String> options=new LinkedHashMap<String, String>();
			for(final String arg : args){final int split=arg.indexOf('='); require(split>0 && split<arg.length()-1, "Use nonempty flag=value arguments."); require(options.put(arg.substring(0, split), arg.substring(split+1))==null, "Duplicate option.");}
			require(options.size()==5 && options.keySet().containsAll(Arrays.asList("mode", "phase", "reads", "fault", "outdir")), "Use exactly mode,phase,reads,fault,outdir.");
			boolean assertions=false; assert(assertions=true); require(assertions, "Use the assertion-enabled launcher.");
			final BBDukCustomMaskPaddingTimingTest test=new BBDukCustomMaskPaddingTimingTest(options);
			if(test.mode.equals("producer")){test.produce();}
			else if(!test.verify()){System.exit(2);}
		}catch(final Throwable failure){failure.printStackTrace(System.err); KillSwitch.kill("Custom-mask stream fixture failed; aborting its dedicated JVM instead of abandoning workers.");}
	}

	private BBDukCustomMaskPaddingTimingTest(final Map<String, String> options) throws Exception{
		mode=options.get("mode"); phase=options.get("phase"); fault=options.get("fault");
		require(mode.equals("producer") || mode.equals("verify") || mode.equals("input"), "Unknown mode.");
		require(phase.equals("baseline") || phase.equals("fixed"), "Unknown phase."); baseline=phase.equals("baseline");
		require(options.get("reads").equals("512") || options.get("reads").equals("262144"), "Only declared finite prefix512 or full262144 reads are supported."); reads=Integer.parseInt(options.get("reads"));
		final boolean finiteFault=(reads==512 && ((mode.equals("producer") && (fault.equals("corrupt") || fault.equals("truncate") || fault.equals("extra"))) || (mode.equals("verify") && fault.equals("readfatal"))));
		require(fault.equals("none") || finiteFault, "Fault injection is restricted to the declared512-read design checks.");
		require(!mode.equals("input") || (reads==262144 && fault.equals("none")), "Input preflight must check the complete immutable input.");
		out=Paths.get(options.get("outdir")).toAbsolutePath().normalize(); require(!out.startsWith(BULK) && !Files.exists(out), "Use a new small-record directory outside bulk."); Files.createDirectory(out);
		require(Files.isRegularFile(INPUT) && !Files.isSymbolicLink(INPUT) && Files.size(INPUT)==INPUT_BYTES, "Immutable full input path/size must match; no input generation is allowed.");
		write("args.txt", lines(options));
		write("scope.txt", "Producer CPU includes complete ordinary process/startup/JIT/loading/masking/formatting/pipe I/O/cleanup/observations. Verifier CPU is separate. Wall is coupled, not isolated kernel speed. Faults are finite design checks only. No bulk generation or BBDuk field/counter writes.\n");
	}

	/** One ordinary process call; only finite negative cases replace fixture stdout. */
	private void produce() throws Exception{
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(SETTINGS));
		args.add("reads="+reads); args.add("literal=ACGT"); args.add("in="+INPUT); args.add("out=stdout.fq");
		write("constructor-args.txt", String.join("\n", args)+"\n");
		write("reflection.tsv", "operation\tmember\nCALL\tordinary public constructor/process once\nREAD_ONLY\tsettings/counters/reference/scaffold release\nBBDUK_FIELD_WRITES\tNONE\nFAULT_STDOUT_WRAPPER\t"+(fault.equals("corrupt") || fault.equals("truncate"))+"\n");
		final PrintStream original=System.out;
		if(fault.equals("corrupt") || fault.equals("truncate")){System.setOut(new PrintStream(new FaultOutput(original, fault, expectedBytes(reads)), false, "US-ASCII"));}
		final BBDuk owner=new BBDuk(args.toArray(new String[args.size()]));
		settings(owner); observe(owner, false); owner.process(); observe(owner, true);
		final long count=number(owner, "basesKTrimmed");
		require(count==(long)reads*(baseline ? 2052 : 3077) && number(owner, "readsKTrimmed")==reads, "Ordinary caller-aggregated mask counts must match the phase literal.");
		if(fault.equals("extra")){System.out.write('!');}
		System.out.flush(); require(!System.out.checkError() && !original.checkError(), "FASTQ stdout must finish without a retained write error.");
		write("producer.tsv", "phase\treads\tlength\tbases\tcaller_count\texpected_output_bytes\tfault\n"+phase+"\t"+reads+"\t4102\t"+((long)reads*LENGTH)+"\t"+count+"\t"+expectedBytes(reads)+"\t"+fault+"\n");
		write("summary.txt", "CUSTOM_MASK_PRODUCER_PASS phase="+phase+" reads="+reads+" processCalls=1 callerCount="+count+" fault="+fault+"\n");
	}

	/** Fixed buffers compare raw records and hash all received bytes; mismatches drain. */
	private boolean verify() throws Exception{
		final boolean inputMode=mode.equals("input");
		final byte[] bases=bytes(inputMode ? repeat("GGACGTGG", 512)+"GGACGT" : repeat(baseline ? "GGAXXXXG" : "GXXXXXXG", 512)+(baseline ? "GGAXXX" : "GXXXXX"));
		final byte[] quality=bytes(repeat("I", LENGTH)), buffer=new byte[65536];
		final ByteBuilder expected=new ByteBuilder(2*LENGTH+64);
		final MessageDigest digest=MessageDigest.getInstance("SHA-256");
		final InputStream stream=ReadWrite.getRawInputStream(inputMode ? INPUT.toString() : "stdin.fq", false);
		long physical=0, compared=0, first=-1; int wanted=-2, actual=-2, available=0, position=0, completed=0;
		boolean eof=false;
		try{
			for(int id=0; id<reads && !eof; id++){
				expected.clear().append("@mask-padding-").append(id).nl().append(bases).nl().append('+').nl().append(quality).nl();
				int offset=0;
				while(offset<expected.length()){
					if(position==available){
						available=stream.read(buffer); position=0;
						if(available<0){if(first<0){first=compared; wanted=expected.array[offset]&255; actual=-1;} eof=true; break;}
						if(available==0){continue;}
						digest.update(buffer, 0, available); physical+=available;
						if(fault.equals("readfatal")){write("injected-fatal.txt", "READFATAL after first actual buffer; controller must stop its owned producer.\n"); System.exit(74);}
					}
					final int count=Math.min(available-position, expected.length()-offset);
					for(int i=0; i<count; i++){if(buffer[position+i]!=expected.array[offset+i] && first<0){first=compared+i; wanted=expected.array[offset+i]&255; actual=buffer[position+i]&255;}}
					position+=count; offset+=count; compared+=count;
				}
				if(!eof){completed++; if((completed&32767)==0){Files.write(out.resolve("progress.tsv"), bytes(completed+"\t"+physical+"\n"), StandardOpenOption.CREATE, StandardOpenOption.APPEND);}}
			}
			if(position<available && first<0){first=compared; wanted=-1; actual=buffer[position]&255;}
			if(!eof){
				int count;
				while((count=stream.read(buffer))>=0){
					if(count==0){continue;}
					if(first<0){first=compared; wanted=-1; actual=buffer[0]&255;}
					digest.update(buffer, 0, count); physical+=count;
				}
			}
		}finally{if(stream!=System.in){stream.close();}}
		final String hash=tailHex(digest.digest());
		final boolean valid=first<0 && completed==reads && physical==expectedBytes(reads);
		if(inputMode){require(valid && hash.equals("f5e6a195c3e9b40de93d"), "Complete immutable input must match both raw literals and known sha80.");}
		write("stream.tsv", "phase\treads\tcomplete_records\texpected_bytes\tactual_bytes\tquality\tfirst_mismatch\texpected_byte\tactual_byte\tsha80\tvalid\n"+phase+"\t"+reads+"\t"+completed+"\t"+expectedBytes(reads)+"\t"+physical+"\t40\t"+first+"\t"+wanted+"\t"+actual+"\t"+hash+"\t"+valid+"\n");
		write("summary.txt", "CUSTOM_MASK_STREAM_"+(valid ? "PASS" : "MISMATCH")+" mode="+mode+" phase="+phase+" reads="+reads+" records="+completed+" bytes="+physical+" drained=true fault="+fault+"\n");
		return valid;
	}

	/** Read only the exact public configuration; X never implies N quality policy. */
	private void settings(final BBDuk owner) throws Exception{
		require(flag(owner, "ktrimN") && !flag(owner, "ktrimLeft") && !flag(owner, "ktrimRight") && !flag(owner, "ksplit") && !flag(owner, "kmaskLowercase") && !flag(owner, "kmaskFullyCovered") && number(owner, "trimSymbol")==88 && number(owner, "trimPad")==1, "Only ordinary positive X replacement masking is timed.");
		require(number(owner, "k")==4 && number(owner, "mink")==4 && !flag(owner, "useShortKmers") && !flag(owner, "rcomp") && !flag(owner, "maskMiddle") && number(owner, "hammingDistance")==0 && number(owner, "editDistance")==0 && number(owner, "qHammingDistance")==0 && number(owner, "qSkip")==1 && number(owner, "minSkip")==1 && number(owner, "maxSkip")==1 && number(owner, "speed")==0, "Use every exact forward full-kmer window; no short/middle/approximate/skipped matches.");
		require(number(owner, "THREADS")==1 && number(owner, "maxReads")==reads && number(owner, "minReadLength")==1 && number(owner, "maxReadLength")==Integer.MAX_VALUE && decimal(owner, "minLenFraction")==0 && !flag(owner, "rename") && !flag(owner, "skipR1") && !flag(owner, "skipR2") && number(owner, "restrictLeft")==0 && number(owner, "restrictRight")==0, "Only the declared unpaired input prefix/full set is requested.");
		require(!flag(owner, "qtrimLeft") && !flag(owner, "qtrimRight") && !flag(owner, "trimClip") && decimal(owner, "minAvgQuality")==0 && number(owner, "minBaseQuality")==0 && number(owner, "maxNs")==-1 && decimal(owner, "maxNRate")==1 && number(owner, "minConsecutiveBases")<=0 && decimal(owner, "minBaseFrequency")==0, "No quality/undefined-base filtering or trimming confounds X output.");
		for(final String name : NONPOS){require(number(owner, name)<=0, "Unrelated operation must remain disabled: "+name);}
		for(final String name : FALSE_FIELDS){require(!flag(owner, name), "Unrelated operation must remain disabled: "+name);}
		require(number(owner, "entropyTrim")==0 && field("ref").get(owner)==null && Arrays.equals((String[])field("literal").get(owner), new String[]{"ACGT"}) && INPUT.toString().equals(field("in1").get(owner)) && "stdout.fq".equals(field("out1").get(owner)), "Use one literal, immutable input and native stdout FASTQ.");
		final ByteBuilder report=new ByteBuilder(); report.append("field\tvalue\n");
		for(final String name : SETTINGS_FIELDS){final Object value=field(name).get(owner); report.append(name).append('\t').append(String.valueOf(value)).nl();}
		write("parsed-settings.tsv", report.toString());
	}

	/** Observe constructor/load/aggregation and normal release without retaining tables. */
	private void observe(final BBDuk owner, final boolean done) throws Exception{
		final long n=done ? reads : 0;
		require(number(owner, "readsIn")==n && number(owner, "readsOut")==n && number(owner, "basesIn")==n*LENGTH && number(owner, "basesOut")==n*LENGTH && number(owner, "storedKmers")==(done ? 1 : 0) && number(owner, "refReads")==(done ? 1 : 0) && number(owner, "refBases")==(done ? 4 : 0) && number(owner, "refKmers")==(done ? 1 : 0), "One real literal load and every requested read must survive normal processing.");
		for(final String name : ZERO_COUNTERS){require(number(owner, name)==0, "Unrelated counter changed: "+name);}
		if(!done){require(number(owner, "readsKTrimmed")==0 && number(owner, "basesKTrimmed")==0, "Fresh mask counters start at zero without seeding.");}
		final AbstractKmerTable[] tables=(AbstractKmerTable[])field("keySets").get(owner); require(tables.length==7, "The ordinary constructor owns seven tables.");
		for(final AbstractKmerTable table : tables){require(done ? table==null : table!=null && table.size()==0, "BBDuk process2 releases tables at1264; constructor tables are empty.");}
		final List<?> names=(List<?>)field("scaffoldNames").get(owner); final IntList lengths=(IntList)field("scaffoldLengths").get(owner);
		require(flag(owner, "RELEASE_TABLES") && (done ? names.isEmpty() && lengths==null : names.size()==1 && "".equals(names.get(0)) && lengths.size==1 && lengths.get(0)==0) && Arrays.equals((int[])field("refScafCounts").get(owner), new int[]{done ? 1 : 0}) && field("scaffoldReadCounts").get(owner)==null && field("scaffoldBaseCounts").get(owner)==null && !owner.errorState, "Normal scaffold cleanup at1276 must run; never disable release to observe retained data.");
		write(done ? "after.tsv" : "before.tsv", "reads_in\tbases_in\treads_out\tbases_out\tstored\tref_reads\tref_bases\tref_kmers\treads_masked\tbases_masked\tcleaned\tunrelated_zero\n"+n+"\t"+(n*LENGTH)+"\t"+n+"\t"+(n*LENGTH)+"\t"+number(owner, "storedKmers")+"\t"+number(owner, "refReads")+"\t"+number(owner, "refBases")+"\t"+number(owner, "refKmers")+"\t"+number(owner, "readsKTrimmed")+"\t"+number(owner, "basesKTrimmed")+"\t"+done+"\ttrue\n");
	}

	/** Finite byte perturbations only; healthy timing never installs this wrapper. */
	private static final class FaultOutput extends OutputStream{
		FaultOutput(final OutputStream out_, final String fault_, final long bytes_){out=out_; fault=fault_; expected=bytes_; assert(expected>0) : "Finite fault boundary must be positive.";}
		@Override
		public void write(final int value) throws IOException{if(!fault.equals("truncate") || seen<expected-1){out.write(fault.equals("corrupt") && seen==0 ? '!' : value);} seen++;}
		@Override
		public void write(final byte[] bytes, final int offset, final int length) throws IOException{
			assert(offset>=0 && length>=0 && offset+length<=bytes.length) : "Caller-owned stdout slice must be valid.";
			if(fault.equals("corrupt") && seen==0 && length>0){out.write('!'); out.write(bytes, offset+1, length-1);}
			else{final int count=fault.equals("truncate") ? (int)Math.max(0, Math.min(length, expected-1-seen)) : length; if(count>0){out.write(bytes, offset, count);}}
			seen+=length;
		}
		@Override
		public void flush() throws IOException{out.flush();}
		private final OutputStream out; private final String fault; private final long expected; private long seen;
	}

	private static long expectedBytes(final int reads){long result=0; for(int i=0; i<reads; i++){int digits=1; for(int x=i; x>=10; x/=10){digits++;} result+=8223+digits;} return result;}
	private static String tailHex(final byte[] digest){assert(digest.length==32) : "SHA-256 produces32bytes; publish only the final80bits."; final char[] hex="0123456789abcdef".toCharArray(), result=new char[20]; for(int i=0; i<10; i++){final int b=digest[digest.length-10+i]&255; result[2*i]=hex[b>>>4]; result[2*i+1]=hex[b&15];} return new String(result);}
	private static String repeat(final String s, final int count){assert(count>0 && !s.isEmpty()) : "Literal body requires positive repetitions."; final StringBuilder b=new StringBuilder(count*s.length()); for(int i=0; i<count; i++){b.append(s);} return b.toString();}
	private static byte[] bytes(final String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static Field field(final String name) throws Exception{final Field f=BBDuk.class.getDeclaredField(name); f.setAccessible(true); return f;}
	private static long number(final BBDuk owner, final String name) throws Exception{return ((Number)field(name).get(owner)).longValue();}
	private static double decimal(final BBDuk owner, final String name) throws Exception{return ((Number)field(name).get(owner)).doubleValue();}
	private static boolean flag(final BBDuk owner, final String name) throws Exception{return field(name).getBoolean(owner);}
	private static void require(final boolean value, final String why){if(!value){throw new IllegalArgumentException(why);}}
	private void write(final String name, final String text) throws IOException{Files.write(out.resolve(name), bytes(text), StandardOpenOption.CREATE_NEW);}
	private static String lines(final Map<String, String> options){final StringBuilder b=new StringBuilder(); for(final Map.Entry<String, String> entry : options.entrySet()){b.append(entry.getKey()).append('=').append(entry.getValue()).append('\n');} return b.toString();}

	private final String mode, phase, fault; private final boolean baseline; private final int reads; private final Path out;
	private static final int LENGTH=4102;
	private static final long INPUT_BYTES=2157071866L;
	private static final Path BULK=Paths.get("/clusterfs/jgi/scratch/gentech/genome_analysis/bbushnell/Shinobu/jgi-mask-padding-timing-20261118"), INPUT=BULK.resolve("input-r262144.fq");
	private static final String[] SETTINGS={"k=4", "mink=-1", "kmask=X", "trimpad=1", "maskfullycovered=f", "rcomp=f", "maskmiddle=f", "hdist=0", "edist=0", "qhdist=0", "qskip=1", "rskip=1", "speed=0", "rename=f", "minlen=1", "qtrim=f", "maxns=-1", "maq=0", "t=1", "simd=f", "qin=33", "qout=33", "array=t", "prealloc=f", "ow=f", "app=f"};
	private static final String[] NONPOS={"forceTrimLeft", "forceTrimRight", "forceTrimRight2", "forceTrimModulo", "trimPolyA", "trimPolyGLeft", "trimPolyGRight", "trimPolyCLeft", "trimPolyCRight", "filterPolyG", "filterPolyC"};
	private static final String[] FALSE_FIELDS={"trimByOverlap", "swift", "ecc", "calcEntropy", "entropyMask", "entropyMark", "filterGC"};
	private static final String[] SETTINGS_FIELDS={"k", "mink", "trimSymbol", "trimPad", "ktrimN", "ktrimLeft", "ktrimRight", "ksplit", "kmaskLowercase", "kmaskFullyCovered", "useShortKmers", "rcomp", "maskMiddle", "hammingDistance", "editDistance", "qHammingDistance", "qSkip", "minSkip", "maxSkip", "speed", "THREADS", "maxReads", "minReadLength", "maxReadLength", "minLenFraction", "qtrimLeft", "qtrimRight", "maxNs", "maxNRate", "minAvgQuality", "minBaseQuality", "minConsecutiveBases", "minBaseFrequency", "rename", "skipR1", "skipR2", "restrictLeft", "restrictRight", "trimClip", "entropyTrim", "RELEASE_TABLES"};
	private static final String[] ZERO_COUNTERS={"readsKFiltered", "basesKFiltered", "readsQTrimmed", "basesQTrimmed", "readsFTrimmed", "basesFTrimmed", "readsQFiltered", "basesQFiltered", "readsNFiltered", "basesNFiltered", "readsEFiltered", "basesEFiltered", "readsPolyTrimmed", "basesPolyTrimmed", "readsTrimmedByOverlap", "basesTrimmedByOverlap", "readsTrimmedBySwift", "basesTrimmedBySwift", "badGcReads", "badGcBases", "badHeaderReads", "badHeaderBases"};
}
