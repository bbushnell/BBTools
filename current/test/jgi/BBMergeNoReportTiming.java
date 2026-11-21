package test.jgi;

import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.lang.reflect.Constructor;
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
import java.util.Map;

import fileIO.ByteStreamWriter;
import jgi.BBMerge;
import parse.Parse;
import shared.Shared;
import stream.Read;
import stream.SamLine;
import stream.Streamer;
import stream.Writer;
import structures.ListNum;

/** Matched in-memory timing of the real no-report BBMerge batch and pair paths.
 * Calls a NEW worker's run method on the measuring thread; no worker is started.
 * Finite input resets owned metadata and allocates one list wrapper per64pairs.
 * Those costs are timed, as are both changed disabled-report conditions. No parser,
 * output writer, graph, training, neural accuracy or whole-CLI claim is made.
 * @author Sucrose
 */
public final class BBMergeNoReportTiming{

	/** mode=pilot|measure, phase=baseline|fixed, root, new outdir, measure-only counts. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("mode", "phase", "root", "outdir", "counts").contains(pair[0]), "Expected documented flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate option: "+pair[0]);
		}
		boolean assertions=false;assert(assertions=true);
		final String mode=options.get("mode"), phase=options.get("phase");
		final boolean pilot="pilot".equals(mode), baseline="baseline".equals(phase);
		require(assertions && (pilot || "measure".equals(mode)) && (baseline || "fixed".equals(phase)), "Use launcher -ea, mode=pilot|measure and phase=baseline|fixed.");
		require(!pilot || baseline, "Only the baseline calibrates matched counts.");
		require(pilot ? !options.containsKey("counts") : options.containsKey("counts"), "Only measure requires counts=BASELINE_PILOT/counts.tsv.");
		final BBMergeNoReportTiming test=new BBMergeNoReportTiming(path(options.get("root")), path(options.get("outdir")), phase);
		try{test.run(pilot, pilot ? null : path(options.get("counts")));}
		finally{test.input.close();}
	}

	private BBMergeNoReportTiming(final Path root_, final Path out_, final String phase_) throws Exception{
		root=root_;out=out_;phase=phase_;baseline=phase.equals("baseline");
		Files.createDirectory(out);
		require(Runtime.getRuntime().availableProcessors()==1 && !Shared.SIMD, "Use AP1 and simd=f for the frozen single-thread workload.");
		require(Paths.get(BBMerge.class.getProtectionDomain().getCodeSource().getLocation().toURI()).equals(root.resolve("current")), "Load production from the single selected runtime, not an overlay.");
		clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported(), "Current-thread CPU time is required; wall time is not a compute substitute.");
		if(!clock.isThreadCpuTimeEnabled()){clock.setThreadCpuTimeEnabled(true);}
		sourcePin=baseline ? BASELINE_PIN : FIXED_PIN;
		require(pin(root.resolve("current/jgi/BBMerge.java")).equals(sourcePin) && pin(root.resolve("networks/bbmerge.bbnet.gz")).equals(MODEL_PIN), "Preserve the reviewed phase source and existing model.");
		for(int i=0; i<RING; i++){ring[i]=new Pair(i);}
		write(out.resolve("A.fq"), fastq(ring[0].a));write(out.resolve("B.fq"), fastq(ring[0].b));
		final BBMerge merge=new BBMerge(arguments().toArray(new String[0]));
		write(out.resolve("settings.tsv"), settings(merge));
		input=new FiniteInput(ring, clock);
		final Class<?> type=Class.forName("jgi.BBMerge$MateThread");
		final Constructor<?> constructor=type.getDeclaredConstructor(BBMerge.class, Streamer.class, Writer.class, Writer.class, baseline ? Writer.class : ByteStreamWriter.class, boolean.class, boolean.class);
		constructor.setAccessible(true);worker=constructor.newInstance(merge, input, null, null, null, false, false);thread=(Thread)worker;
		require(thread.getState()==Thread.State.NEW, "The timing caller must never start the worker thread.");
		for(int i=0; i<COUNTERS.length; i++){counterFields[i]=field(type, COUNTERS[i]);}
		hist=(long[])value(worker, "hist");require(hist.length==2000, "Literal histogram oracle requires the reviewed2000bins.");
		require(value(worker, "rosi")==null && value(worker, "rosgood")==null && value(worker, "rosbad")==null && value(worker, "loglogT")==null && value(worker, "bridgeT")==null, "No report/read writer, cardinality or graph helper may participate.");
		//The existing constructor may copy net0 even with MAKE_VECTOR=false; no inference runs in this workload.
		input.start();
		write(out.resolve("measurements.tsv"), "stage\tcalls\tbatches\twallNanos\tcpuNanos\tpairs\tbases\thist100\tincorrect\tinsertSum\tchecksum\tworkerNEW\tfailureNull\n");
		write(out.resolve("counters.tsv"), "stage\tfield\tdelta\texpected\n");
		write(out.resolve("histograms.tsv"), "stage\tbin\tdelta\texpected\n");
		write(out.resolve("ring.tsv"), ringText());
		write(out.resolve("oracles.tsv"), "recipe\tring\tperfect\timperfect\tr1Length\tr2Length\tinsert\tperfectErrors\timperfectErrors\tbaseCountPerPair\n"+RECIPE+"\t64\t32\t32\t160\t160\t100\t0\t1\t320\n");
		write(out.resolve("runtime.txt"), "phase="+phase+"\nsource.sha80="+sourcePin+"\nmodel.sha80="+MODEL_PIN+"\njava.version="+System.getProperty("java.version")+"\nactiveProcessors=1\nsimd=false\nthreadStarts=0\ninput=owned finite ring, no parser\ntimed=worker.run+processReads+findOverlapInThread+processReadPair+metadata-reset+batch-wrapper\npilotBudget=90CPU seconds checked at batch boundaries; one64-pair batch can cross the boundary\n");
	}

	private void run(final boolean pilot, final Path counts) throws Exception{
		final byte[][] reference=pilot ? null : reference(counts);
		final long saved=pilot ? INITIAL_CALLS : savedCalls(reference[0]);
		if(reference!=null){sameReference(reference);}
		if(pilot){input.deadline=clock.getCurrentThreadCpuTime()+PILOT_BUDGET;}
		final String before=validate("validate-before");write(out.resolve("validation.tsv"), before);
		if(reference!=null){require(Arrays.equals(bytes(before), reference[4]), "Literal validation differs from baseline pilot.");}
		sample("warmup", INITIAL_CALLS);
		long calls=saved;long[] selected=null;
		for(int attempt=1; attempt<=(pilot ? MAX_ATTEMPTS : 1); attempt++){
			selected=sample(pilot ? "pilot-"+attempt : "measure", calls);
			if(!pilot || selected[1]>=PILOT_MIN){break;}
			require(attempt<MAX_ATTEMPTS, "Pilot attempt cap reached; review sizing before another invocation.");
			final long scaled=(long)Math.min(calls*16.0, Math.ceil(calls*(TARGET_CPU/(double)selected[1])));
			calls=((Math.max(calls+RING, scaled)+RING-1)/RING)*RING;
			require(calls<=MAX_CALLS, "Pilot call cap exceeded; review sizing instead of silently expanding.");
		}
		require(selected!=null && selected[1]>=(pilot ? PILOT_MIN : MIN_CPU), "Compute sample too short; recalibrate baseline counts rather than resizing fixed only.");
		final String after=validate("validate-after");write(out.resolve("validation-after.tsv"), after);require(before.equals(after), "Literal ring outputs changed after timing.");
		if(reference!=null){final byte[][] again=reference(counts);for(int i=0; i<again.length; i++){require(Arrays.equals(reference[i], again[i]), "Saved baseline reference changed: "+i);}}
		require(pin(root.resolve("current/jgi/BBMerge.java")).equals(sourcePin) && pin(root.resolve("networks/bbmerge.bbnet.gz")).equals(MODEL_PIN), "Source/model changed during measurement.");
		if(pilot){write(out.resolve("counts.tsv"), COUNTS_HEADER+RECIPE+"\tbaseline\tpilot\t64\t"+calls+"\t"+selected[1]+"\n");}
		final String summary="BBMERGE_NO_REPORT_TIMING_PASS mode="+(pilot ? "pilot" : "measure")+" phase="+phase+" ring=64 calls="+calls+" batches="+(calls/RING)+" cpuNanos="+selected[1]+" wallNanos="+selected[2]+" pairDelta="+calls+" baseDelta="+(320*calls)+" hist100="+calls+" workerNEW=true failureNull=true literalRing=true\n";
		write(out.resolve("summary.txt"), summary);System.out.print(summary);
	}

	/** Includes real worker dispatch, owned metadata reset and wrapper allocation. Checks/I/O follow both clocks. */
	private long[] sample(final String stage, final long calls) throws Exception{
		require(calls>=RING && calls<=MAX_CALLS && calls%RING==0 && thread.getState()==Thread.State.NEW && value(worker, "failure")==null, "Use whole bounded ring cycles and a never-started, successful worker.");
		final long[] before=snapshot();final long[] oldHist=hist.clone();input.configure(calls/RING);
		final long cpuStart=clock.getCurrentThreadCpuTime(), wallStart=System.nanoTime();
		thread.run();
		final long wall=System.nanoTime()-wallStart, cpu=clock.getCurrentThreadCpuTime()-cpuStart;
		require(value(worker, "failure")==null && thread.getState()==Thread.State.NEW, "The real worker captured a failure or was started: "+value(worker, "failure"));
		require(input.delivered==calls/RING && input.terminalCalls==1 && input.owner==Thread.currentThread(), "Finite input did not deliver exactly the planned batches and one terminal call.");
		final long[] after=snapshot();
		final StringBuilder counterTable=new StringBuilder(), histogramTable=new StringBuilder();
		for(int i=0; i<COUNTERS.length; i++){
			final long delta=after[i]-before[i], expected=expectedCounter(COUNTERS[i], calls);
			counterTable.append(stage).append('\t').append(COUNTERS[i]).append('\t').append(delta).append('\t').append(expected).append('\n');require(delta==expected, "Wrong actual counter delta: "+stage+"/"+COUNTERS[i]+"/"+delta+"/"+expected);
		}
		for(int i=0; i<hist.length; i++){
			final long delta=hist[i]-oldHist[i], expected=i==100 ? calls : 0;
			histogramTable.append(stage).append('\t').append(i).append('\t').append(delta).append('\t').append(expected).append('\n');require(delta==expected, "Wrong insert histogram bin: "+stage+"/"+i+"/"+delta);
		}
		append(out.resolve("counters.tsv"), counterTable.toString());append(out.resolve("histograms.tsv"), histogramTable.toString());
		require(((Number)value(worker, "insertMin")).intValue()==100 && ((Number)value(worker, "insertMax")).intValue()==100, "Every accepted candidate must have insert100.");
		long checksum=0;for(final Pair pair : ring){pair.check();checksum+=pair.a.insert()+pair.b.insert()+pair.a.errors+pair.b.errors;}
		require(checksum==13408 && cpu>0 && wall>0, "Literal final ring or CPU clock failed: "+checksum);
		append(out.resolve("measurements.tsv"), stage+"\t"+calls+"\t"+(calls/RING)+"\t"+wall+"\t"+cpu+"\t"+calls+"\t"+(320*calls)+"\t"+calls+"\t"+calls+"\t"+(100*calls)+"\t"+checksum+"\ttrue\ttrue\n");
		return new long[]{calls, cpu, wall};
	}

	/** Full literal state is serialized only outside the timed worker call. */
	private String validate(final String stage) throws Exception{
		sample(stage, RING);final StringBuilder text=new StringBuilder("slot\terrors\tinsert1\tinsert2\tA\tB\tAquality\tBquality\tAflags\tBflags\n");
		for(final Pair p : ring){p.check();text.append(p.index).append('\t').append(p.a.errors).append("\t100\t100\t").append(new String(p.a.bases, StandardCharsets.US_ASCII)).append('\t').append(new String(p.b.bases, StandardCharsets.US_ASCII)).append('\t').append(Arrays.toString(p.a.quality)).append('\t').append(Arrays.toString(p.b.quality)).append('\t').append(p.a.flags).append('\t').append(p.b.flags).append('\n');}
		return text.toString();
	}
	private long[] snapshot() throws IllegalAccessException{final long[] values=new long[counterFields.length];for(int i=0; i<values.length; i++){values[i]=counterFields[i].getLong(worker);}return values;}
	private static long expectedCounter(final String name, final long calls){return name.equals("pairsProcessed") || name.equals("incorrectCount") ? calls : name.equals("basesProcessed") ? 320*calls : name.equals("insertSumIncorrect") ? 100*calls : 0;}
	private String ringText(){final StringBuilder s=new StringBuilder("slot\tA\tB\tAquality\tBquality\n");for(final Pair p : ring){s.append(p.index).append('\t').append(new String(p.ab, StandardCharsets.US_ASCII)).append('\t').append(new String(p.bb, StandardCharsets.US_ASCII)).append('\t').append(Arrays.toString(p.aq)).append('\t').append(Arrays.toString(p.bq)).append('\n');}return s.toString();}
	private void sameReference(final byte[][] reference) throws IOException{final String[] files={"ring.tsv", "oracles.tsv", "settings.tsv"};for(int i=0; i<files.length; i++){require(Arrays.equals(Files.readAllBytes(out.resolve(files[i])), reference[i+1]), "Frozen workload differs from baseline: "+files[i]);}}
	private static byte[][] reference(final Path counts) throws IOException{final Path dir=counts.getParent();return new byte[][]{Files.readAllBytes(counts), Files.readAllBytes(dir.resolve("ring.tsv")), Files.readAllBytes(dir.resolve("oracles.tsv")), Files.readAllBytes(dir.resolve("settings.tsv")), Files.readAllBytes(dir.resolve("validation.tsv"))};}
	private static long savedCalls(final byte[] data){final String[] rows=new String(data, StandardCharsets.US_ASCII).split("\n");require(rows.length==2 && (rows[0]+"\n").equals(COUNTS_HEADER), "Expected the exact saved count schema.");final String[] f=rows[1].split("\t", -1);require(f.length==6 && f[0].equals(RECIPE) && f[1].equals("baseline") && f[2].equals("pilot") && f[3].equals("64"), "Saved pilot recipe differs.");final long calls=Long.parseLong(f[4]);require(calls>=RING && calls<=MAX_CALLS && calls%RING==0 && Long.parseLong(f[5])>=PILOT_MIN, "Saved pilot must qualify with bounded whole-ring counts.");return calls;}

	/** Keeps every workload selection explicit; constructor reads only the tiny descriptor mates. */
	private ArrayList<String> arguments(){return new ArrayList<String>(Arrays.asList("in1="+out.resolve("A.fq"), "in2="+out.resolve("B.fq"), "t=1", "simd=f", "qin=33", "qout=33", "ow=f", "app=f", "merge=f", "ecco=f", "changequality=f", "qtrim=f", "maq=0", "minlength=1", "usequality=t", "ouq=f", "owq=t", "ratiomode=t", "flatmode=f", "entropy=f", "minoverlap=20", "minoverlap0=20", "ratiominoverlapreduction=0", "mininsert=1", "mininsert0=1", "efilter=f", "pfilter=0", "mee=0", "tof=0", "extend=0", "extend2=0", "kfilter=0", "nn=f", "maxratio=0.05", "maxmismatches=5"));}
	private String settings(final BBMerge merge) throws Exception{
		final String[] yes={"useRatioMode", "useQuality", "overlapWithoutQuality"};
		final String[] no={"MAKE_VECTOR", "TRAINING", "TAG_CUSTOM", "parseCustom", "qtrim1", "qtrim2", "qtrimLeft", "qtrimRight", "ecco", "join", "changeQuality", "useFlatMode", "useEntropy", "useEfilter", "useMEEfilter", "overlapUsingQuality", "kmerBridge", "eccTadpole", "eccBloom", "verifyAdapters", "findAdapterSequence", "recalibrateQuality", "iupacToN"};
		final StringBuilder s=new StringBuilder("field\tvalue\n");for(final String f : yes){require(Boolean.TRUE.equals(value(merge, f)), "Expected enabled setting: "+f);s.append(f).append("\ttrue\n");}for(final String f : no){require(Boolean.FALSE.equals(value(merge, f)), "Expected disabled setting: "+f);s.append(f).append("\tfalse\n");}
		final String[] numbers={"THREADS", "minReadLength", "minInsert", "minInsert0", "MIN_OVERLAPPING_BASES", "MIN_OVERLAPPING_BASES_0", "MIN_OVERLAPPING_BASES_RATIO_REDUCTION", "MAX_MISMATCHES_R", "pfilterRatio", "extendRight1", "extendRight2", "forceTrimLeft", "forceTrimRight", "forceTrimRight2", "forceTrimModulo"};
		final int[] expected={1, 1, 20, 1, 20, 20, 0, 5, 0, 0, 0, -1, -1, -1, -1};for(int i=0; i<numbers.length; i++){require(((Number)value(merge, numbers[i])).doubleValue()==expected[i], "Wrong literal numeric setting: "+numbers[i]);s.append(numbers[i]).append('\t').append(expected[i]).append('\n');}
		require(((Number)value(merge, "MAX_RATIO")).floatValue()==.05f && value(merge, "tadpole")==null && value(merge, "bloomFilter")==null && value(merge, "outinsert")==null && value(merge, "out1")==null && value(merge, "outb1")==null, "The bounded search/report/output selection changed.");
		return s.append("MAX_RATIO\t0.05\ntadpole\tnull\nbloomFilter\tnull\noutinsert\tnull\nout1\tnull\noutb1\tnull\n").toString();
	}

	/** Reuses bytes only after the synchronous worker finishes the previous batch. */
	private static final class FiniteInput implements Streamer{
		FiniteInput(final Pair[] ring_, final ThreadMXBean clock_){ring=ring_;clock=clock_;for(final Pair p : ring){reads.add(p.a);}}
		void configure(final long count){require(started && !closed && count>0 && owner==Thread.currentThread(), "Configure only the owned idle finite input.");limit=count;delivered=terminalCalls=0;}
		@Override public ListNum<Read> nextList(){
			require(started && !closed && owner==Thread.currentThread(), "Finite input is confined to its single measuring thread.");
			if(deadline!=Long.MAX_VALUE){require(clock.getCurrentThreadCpuTime()<deadline, "Aggregate pilot90CPU-second budget reached at a batch boundary.");}
			if(delivered==limit){terminalCalls++;return null;}
			for(final Pair p : ring){p.restore();}return new ListNum<Read>(reads, delivered++);
		}
		@Override public void start(){require(!started, "Finite input starts once.");started=true;owner=Thread.currentThread();}
		@Override public void close(){closed=true;}
		@Override public String fname(){return "owned-no-report-ring";}
		@Override public boolean paired(){return true;}
		@Override public int pairnum(){return 0;}
		@Override public long readsProcessed(){return delivered*RING;}
		@Override public long basesProcessed(){return delivered*RING*320L;}
		@Override public void setSampleRate(final float rate, final long seed){require(rate==1, "The finite timing ring is never sampled.");}
		@Override public ListNum<SamLine> nextLines(){throw new UnsupportedOperationException("Finite timing supplies reads only.");}
		@Override public boolean hasMore(){return delivered<limit;}
		@Override public boolean errorState(){return false;}
		final Pair[] ring;final ThreadMXBean clock;final ArrayList<Read> reads=new ArrayList<Read>(RING);
		Thread owner;long limit, delivered, terminalCalls, deadline=Long.MAX_VALUE;boolean started, closed;
	}
	/** Literal insert100 pair; only flags, insert and error scratch need resetting. */
	private static final class Pair{
		Pair(final int index_){
			index=index_;final byte[] aBases=bytes(F100+ADAPTER_A), bBases=bytes(RAW_F100+ADAPTER_B);
			final byte[] aQual=new byte[160], bQual=new byte[160];Arrays.fill(aQual, (byte)20);Arrays.fill(bQual, (byte)20);Arrays.fill(aQual, 100, 160, (byte)40);Arrays.fill(bQual, 100, 160, (byte)40);
			if((index&1)!=0){bBases[49]='A';aQual[50]=30;}
			a=new Read(aBases, aQual, "ring-"+index+"/1", index, false);b=new Read(bBases, bQual, "ring-"+index+"/2", index, false);a.mate=b;b.mate=a;b.setPairnum(1);a.validate(true);b.validate(true);
			ab=a.bases;bb=b.bases;aq=a.quality;bq=b.quality;abSaved=ab.clone();bbSaved=bb.clone();aqSaved=aq.clone();bqSaved=bq.clone();af=a.flags;bf=b.flags;
			require(ab.length==160 && bb.length==160 && F100.length()==100 && RAW_F100.length()==100 && F100.charAt(50)=='G' && RAW_F100.charAt(49)=='C', "Frozen literal geometry requires insert100 and one raw mate substitution.");
		}
		void restore(){a.setInsert(-1);a.flags=af;b.flags=bf;a.errors=7;b.errors=9;}
		void check(){
			require(a.bases==ab && b.bases==bb && a.quality==aq && b.quality==bq && Arrays.equals(ab, abSaved) && Arrays.equals(bb, bbSaved) && Arrays.equals(aq, aqSaved) && Arrays.equals(bq, bqSaved), "Actual worker changed literal arrays/references: "+index);
			require(a.insert()==100 && b.insert()==100 && a.insertvalid() && b.insertvalid() && a.flags==(af|Read.INSERTMASK) && b.flags==(bf|Read.INSERTMASK) && a.errors==(index&1) && b.errors==9, "Actual no-report pair state differs: "+index);
			require(a.mate==b && b.mate==a && a.numericID==index && b.numericID==index && a.id.equals("ring-"+index+"/1") && b.id.equals("ring-"+index+"/2") && a.obj==null && b.obj==null, "Pair identity or attachment changed: "+index);
			require(a.chrom==-1 && b.chrom==-1 && a.start==-1 && b.start==-1 && a.stop==-1 && b.stop==-1 && a.mapScore==0 && b.mapScore==0 && a.copies==1 && b.copies==1, "Unrelated pair metadata changed: "+index);
		}
		final int index, af, bf;final Read a, b;final byte[] ab, bb, aq, bq, abSaved, bbSaved, aqSaved, bqSaved;
	}

	private static Field field(final Class<?> type, final String name) throws NoSuchFieldException{final Field f=type.getDeclaredField(name);f.setAccessible(true);return f;}
	private static Object value(final Object owner, final String name) throws Exception{return field(owner.getClass(), name).get(owner);}
	private static Path path(final String s){require(s!=null && !s.isEmpty(), "Missing required path.");return Paths.get(s).toAbsolutePath().normalize();}
	private static byte[] bytes(final String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static void require(final boolean okay, final String why){if(!okay){throw new IllegalArgumentException(why);}}
	private static void write(final Path p, final String s) throws IOException{Files.write(p, bytes(s), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path p, final String s) throws IOException{Files.write(p, bytes(s), StandardOpenOption.APPEND);}
	private static String fastq(final Read r){final StringBuilder s=new StringBuilder("@").append(r.id).append('\n').append(new String(r.bases, StandardCharsets.US_ASCII)).append("\n+\n");for(final byte q : r.quality){s.append((char)(q+33));}return s.append('\n').toString();}
	private static String pin(final Path p) throws Exception{final byte[] digest=MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(p));final StringBuilder s=new StringBuilder();for(int i=digest.length-10; i<digest.length; i++){s.append(Character.forDigit((digest[i]>>>4)&15, 16)).append(Character.forDigit(digest[i]&15, 16));}return s.toString();}
	private final Path root, out;private final String phase, sourcePin;private final boolean baseline;
	private final ThreadMXBean clock;private final Pair[] ring=new Pair[RING];private final FiniteInput input;private final Object worker;private final Thread thread;private final long[] hist;
	private final Field[] counterFields=new Field[COUNTERS.length];
	private static final String[] COUNTERS={"pairsProcessed", "basesProcessed", "matedCount", "correctCount", "unmergedOverlappingCount", "ambiguousCount", "tooShortCount", "tooLongCount", "incorrectCount", "noSolutionCount", "insertSumCorrect", "insertSumIncorrect", "fullyExtendedT", "partlyExtendedT", "notExtendedT", "extensionsAttemptedT", "bridgeAttemptsT", "bridgeSuccessT", "bridgeNoPathT", "bridgeNoTargetT", "bridgeAmbiguousT", "bridgeCycleT", "bridgeNonreciprocalT", "bridgeRejectedT", "adaptersExpectedT", "adaptersFoundT", "errorsCorrectedT"};
	private static final int RING=64, MAX_ATTEMPTS=6;private static final long INITIAL_CALLS=32768, MAX_CALLS=536870912;
	private static final long MIN_CPU=10000000000L, PILOT_MIN=15000000000L, TARGET_CPU=18000000000L, PILOT_BUDGET=90000000000L;
	private static final String RECIPE="bbmerge-no-report-ring-v1", COUNTS_HEADER="recipe\tphase\tmode\tring\tcalls\tpilotCpuNanos\n";
	private static final String BASELINE_PIN="989c3ef47aebab252185", FIXED_PIN="5a4527a67db65afb8ab1", MODEL_PIN="24998a8aa257552d4a22";
	private static final String F100="ACGTTGCATGTCGACATGCTAGCTTCGATCGATGACCTGATCCGATAGCAGTACCGTATCATGCTAGGTCCGATCTACGATACGGTTCAGCAGTACGATG";
	private static final String RAW_F100="CATCGTACTGCTGAACCGTATCGTAGATCGGACCTAGCATGATACGGTACTGCTATCGGATCAGGTCATCGATCGAAGCTAGCATGTCGACATGCAACGT";
	private static final String ADAPTER_A="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTG";
	private static final String ADAPTER_B="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGC";
}
