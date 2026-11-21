package test.jgi;

import java.io.IOException;
import java.lang.invoke.MethodHandle;
import java.lang.invoke.MethodHandles;
import java.lang.invoke.MethodType;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
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

import jgi.BBMerge;
import parse.Parse;
import stream.Read;
import stream.Streamer;
import stream.Writer;

/** Times the actual worker search/restoration path with 64 owned nonnull-quality pairs.
 * The literal 40/56 pair is the public null-quality fixture's passing FASTQ control.
 * Each call starts in raw sequencer orientation with Q30 qualities, usequality=f,
 * no neural model evaluation, trimming, correction or extension. The measured region
 * includes Read reference/metadata reset, ring lookup, typed worker dispatch and a
 * rotating output checksum. It excludes setup reflection, full validation and I/O.
 * This is neither isolated restoreQualities cost nor joined-output/CLI throughput.
 * The external runner owns build/runtime provenance and three fresh JVMs per phase.
 * Setup and pilot/count handling follow ErrorFilterTiming and RetainedQualityTiming.
 * @author Sucrose
 */
public final class BBMergeNullQualityTiming{

	/** mode=pilot|measure, phase=baseline|fixed, outdir=NEW; measure needs counts=BASELINE/counts.tsv. */
	public static void main(final String[] args) throws Throwable{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("mode", "phase", "outdir", "counts").contains(pair[0]), "Expected documented flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate option: "+pair[0]);
		}
		final String mode=options.get("mode"), phase=options.get("phase");
		require("pilot".equals(mode) || "measure".equals(mode), "Set mode=pilot|measure.");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		final boolean pilot=mode.equals("pilot");
		require(!pilot || phase.equals("baseline"), "Only the baseline pilot selects shared call counts.");
		require(pilot ? !options.containsKey("counts") : options.containsKey("counts"), "Only measure requires counts=BASELINE_PILOT/counts.tsv.");
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbmergenullqualitytiming.sh with assertions enabled.");
		final Path saved=pilot ? null : path(options.get("counts"));
		final byte[][] reference=saved==null ? null : readReference(saved);
		final BBMergeNullQualityTiming test=new BBMergeNullQualityTiming(path(options.get("outdir")), phase);
		final long calls=pilot ? INITIAL_CALLS : savedCalls(reference[0]);
		final String validation=test.validateRing();
		write(test.out.resolve("validation.tsv"), validation);
		final String pins=test.inputPins(); write(test.out.resolve("input-pins.tsv"), pins);
		if(reference!=null){test.compareReference(reference);}
		test.sample("warmup", 0, INITIAL_CALLS);
		Result selected=null; long candidate=calls;
		for(int attempt=1; attempt<=(pilot ? MAX_ATTEMPTS : 1); attempt++){
			selected=test.sample(mode, attempt, candidate);
			if(!pilot || selected.cpu>=MIN_PILOT_CPU){break;}
			if(attempt<MAX_ATTEMPTS){candidate=nextCalls(candidate, selected.cpu);}
		}
		require(selected!=null, "Every invocation must retain its measured result.");
		final String after=test.validateRing(); write(test.out.resolve("validation-after.tsv"), after);
		require(validation.equals(after), "Full ring outputs or state changed after timing.");
		require(pins.equals(test.inputPins()), "Frozen workload files changed during timing.");
		if(reference!=null){
			final byte[][] afterReference=readReference(saved);
			for(int i=0; i<reference.length; i++){require(Arrays.equals(reference[i], afterReference[i]), "Saved pilot input changed during measurement: "+i);}
		}
		final boolean enough=selected.cpu>=(pilot ? MIN_PILOT_CPU : MIN_CPU);
		if(pilot && enough){
			write(test.out.resolve("counts.tsv"), COUNTS_HEADER+"\n"+RECIPE+"\tbaseline\tpilot\t64\t"+selected.calls+"\t40\t56\t40\t"+selected.cpu+"\n");
		}
		final String summary="BBMERGE_NULL_QUALITY_TIMING_"+(enough ? "PASS" : "SHORT")+" mode="+mode+" phase="+phase+
			" ring=64 validationCalls=128 calls="+selected.calls+" r1Length=40 r2Length=56 insert=40"+
			" wallNanos="+selected.wall+" cpuNanos="+selected.cpu+" checksum="+selected.checksum+
			" pairDelta="+selected.pairs+" baseDelta="+selected.bases+" arraysUnchanged=true stateMatches=true\n";
		write(test.out.resolve("summary.txt"), summary); System.out.print(summary);
		require(enough, "Compute region too short; recalibrate matched baseline counts, never increase fixed alone.");
	}

	/** Creates tiny constructor inputs and one unstarted worker; no streams are started. */
	private BBMergeNullQualityTiming(final Path out_, final String phase_) throws Exception{
		out=out_; phase=phase_; Files.createDirectory(out);
		clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported(), "Current measuring-thread CPU time is required; wall time is not a substitute.");
		if(!clock.isThreadCpuTimeEnabled()){clock.setThreadCpuTimeEnabled(true);}
		write(out.resolve("A.fq"), fastq(FRAGMENT)); write(out.resolve("B.fq"), fastq(RAW_INPUT));
		final ArrayList<String> arguments=new ArrayList<String>(Arrays.asList(COMMON));
		arguments.add("in1="+out.resolve("A.fq")); arguments.add("in2="+out.resolve("B.fq"));
		final StringBuilder settings=new StringBuilder("index\targument\n");
		for(int i=0; i<COMMON.length; i++){settings.append(i).append('\t').append(COMMON[i]).append('\n');}
		// Paths vary by fresh JVM; literal input bytes are pinned separately.
		write(out.resolve("settings.tsv"), settings.toString());
		final BBMerge merge=new BBMerge(arguments.toArray(new String[0]));
		require(Runtime.getRuntime().availableProcessors()==1 && !shared.Shared.SIMD, "Use one active processor and SIMD disabled.");
		for(final String flag : new String[]{"useQuality", "MAKE_VECTOR", "TAG_CUSTOM", "parseCustom", "ONLY_OUTPUT_FALSE_POSITIVE", "useFlatMode", "useEntropy", "useEfilter", "useMEEfilter", "qtrim1", "qtrim2", "qtrimLeft", "qtrimRight", "recalibrateQuality", "iupacToN", "verifyAdapters", "useKFilter", "USE_MAPPING", "kmerBridge"}){
			require(Boolean.FALSE.equals(value(merge, flag)), "Frozen ordinary quality-disabled worker requires "+flag+"=false.");
		}
		require(Boolean.TRUE.equals(value(merge, "useRatioMode")) && Boolean.TRUE.equals(value(merge, "MATE_BY_OVERLAP")), "The frozen worker must actually search for its overlap.");
		for(final String count : new String[]{"extendRight1", "extendRight2", "forceTrimLeft", "forceTrimRight", "forceTrimRight2", "forceTrimModulo", "minAvgQuality", "pfilterRatio"}){
			require(((Number)value(merge, count)).doubleValue()<=0, "Frozen worker excludes positive optional preprocessing/filtering settings: "+count);
		}
		require(value(merge, "tadpole")==null && value(merge, "bloomFilter")==null, "No graph or Bloom correction/extension belongs in this workload.");
		final Class<?> type=Class.forName("jgi.BBMerge$MateThread");
		final Constructor<?> constructor=type.getDeclaredConstructor(BBMerge.class, Streamer.class, Writer.class, Writer.class, Writer.class, boolean.class, boolean.class);
		constructor.setAccessible(true); worker=constructor.newInstance(merge, null, null, null, null, true, false);
		require(((Thread)worker).getState()==Thread.State.NEW, "The worker is invoked on the measuring thread, never started.");
		pairsField=field(type, "pairsProcessed"); basesField=field(type, "basesProcessed");
		for(int i=0; i<RING; i++){ring[i]=new Pair(i);}
		write(out.resolve("ring.tsv"), inputs());
		write(out.resolve("oracles.tsv"), oracles());
		write(out.resolve("measurements.tsv"), "stage\tattempt\tphase\tcalls\tchecksum\texpectedChecksum\twallNanos\tcpuNanos\tpairDelta\tbaseDelta\tarraysUnchanged\tstateMatches\n");
		write(out.resolve("runtime.txt"), "java.version="+System.getProperty("java.version")+"\njava.vendor="+System.getProperty("java.vendor")+
			"\njava.vm.name="+System.getProperty("java.vm.name")+"\nactiveProcessors=1\nsimd=false\nrecipe="+RECIPE+
			"\nkernel=processReadPair+Read-reset+typed-dispatch+checksum\nparsecustom=false: outer truth parsing is outside processReadPair\n");
	}

	/** Full array/reference/metadata checks and actual counter deltas stay outside clocks. */
	private String validateRing() throws Throwable{
		final long beforePairs=pairsField.getLong(worker), beforeBases=basesField.getLong(worker);
		final StringBuilder table=new StringBuilder("ring\treturn\tA\tB\tAquality\tBquality\tAflags\tBflags\tAerrors\tBerrors\tchecksum\n");
		cycleChecksum=0;
		for(int i=0; i<RING; i++){
			final Pair p=ring[i]; require(p.arraysUnchanged(), "Golden input arrays changed before validation: "+i); p.restore();
			require(p.stateMatches(false), "Reset failed to restore complete input state: "+i);
			final int returned=(int)PROCESS_PAIR.invokeExact(worker, p.a, p.b);
			require(returned==40 && p.arraysUnchanged() && p.stateMatches(true), "Actual worker did not return insert40 and restore full raw pair state: "+i);
			final long actual=p.checksum(returned), expected=p.expectedChecksum();
			require(actual==expected, "Output checksum does not consume the literal raw pair/retained quality identities: "+i);
			cycleChecksum=Math.addExact(cycleChecksum, expected);
			table.append(i).append('\t').append(returned).append('\t').append(new String(p.a.bases, StandardCharsets.US_ASCII)).append('\t')
				.append(new String(p.b.bases, StandardCharsets.US_ASCII)).append('\t').append(Arrays.toString(p.a.quality)).append('\t')
				.append(Arrays.toString(p.b.quality)).append('\t').append(p.a.flags).append('\t').append(p.b.flags).append('\t')
				.append(p.a.errors).append('\t').append(p.b.errors).append('\t').append(actual).append('\n');
		}
		require(pairsField.getLong(worker)-beforePairs==RING && basesField.getLong(worker)-beforeBases==96L*RING, "preprocess must account for every actual 40+56-base pair.");
		return table.toString();
	}

	/** No reflective calls, full-state walks, array reconstruction or I/O inside either clock. */
	private Result sample(final String stage, final int attempt, final long calls) throws Throwable{
		require(calls>0 && calls<=MAX_CALLS && calls%RING==0, "Checksum and sizing require bounded whole-ring cycles: "+calls);
		require(validOutputs(), "Full state changed before sample.");
		final long beforePairs=pairsField.getLong(worker), beforeBases=basesField.getLong(worker);
		long checksum=0;
		final long cpuStart=clock.getCurrentThreadCpuTime(), wallStart=System.nanoTime();
		for(long i=0; i<calls; i++){
			final Pair p=ring[(int)i&(RING-1)]; p.restore();
			final int returned=(int)PROCESS_PAIR.invokeExact(worker, p.a, p.b);
			checksum+=p.checksum(returned);
		}
		final long wall=System.nanoTime()-wallStart, cpu=clock.getCurrentThreadCpuTime()-cpuStart;
		final long pairDelta=pairsField.getLong(worker)-beforePairs, baseDelta=basesField.getLong(worker)-beforeBases;
		final long expected=Math.multiplyExact(cycleChecksum, calls/RING);
		boolean arrays=true, state=true;
		for(final Pair p : ring){arrays&=p.arraysUnchanged(); state&=p.stateMatches(true);}
		Files.write(out.resolve("measurements.tsv"), bytes(stage+"\t"+attempt+"\t"+phase+"\t"+calls+"\t"+checksum+"\t"+expected+"\t"+wall+"\t"+cpu+
			"\t"+pairDelta+"\t"+baseDelta+"\t"+arrays+"\t"+state+"\n"), StandardOpenOption.APPEND);
		require(checksum==expected && arrays && state && pairDelta==calls && baseDelta==96L*calls && cpu>0 && wall>0,
			"Timed checksum/state/arrays/accounting/clock failure; inspect measurements.tsv.");
		return new Result(calls, checksum, wall, cpu, pairDelta, baseDelta);
	}

	private boolean validOutputs(){for(final Pair p : ring){if(!p.arraysUnchanged() || !p.stateMatches(true)){return false;}} return true;}
	private String inputs(){
		final StringBuilder table=new StringBuilder("ring\tA\tB\tAquality\tBquality\tAflags\tBflags\tAerrors\tBerrors\n");
		for(int i=0; i<RING; i++){
			final Pair p=ring[i]; table.append(i).append('\t').append(FRAGMENT).append('\t').append(RAW_INPUT).append('\t')
				.append(Arrays.toString(p.aqGolden)).append('\t').append(Arrays.toString(p.bqGolden)).append('\t')
				.append(p.aFlags).append('\t').append(p.bFlags).append("\t7\t9\n");
		}
		return table.toString();
	}
	private String oracles(){
		final StringBuilder table=new StringBuilder("ring\tinsert\tAerrors\tBerrors\tqualityIdentity\tbaseIdentity\tstate\tchecksum\n");
		for(int i=0; i<RING; i++){table.append(i).append("\t40\t0\t9\toriginal\toriginal\traw-validated-pair\t").append(ring[i].expectedChecksum()).append('\n');}
		return table.toString();
	}
	private String inputPins() throws Exception{
		final StringBuilder table=new StringBuilder("file\tsha80\n");
		for(final String name : REFERENCE_FILES){table.append(name).append('\t').append(pin(Files.readAllBytes(out.resolve(name)))).append('\n');}
		return table.toString();
	}
	private void compareReference(final byte[][] reference) throws IOException{
		for(int i=0; i<REFERENCE_FILES.length; i++){
			require(Arrays.equals(reference[i+1], Files.readAllBytes(out.resolve(REFERENCE_FILES[i]))), "Frozen baseline input/oracle/state differs: "+REFERENCE_FILES[i]);
		}
	}
	private static byte[][] readReference(final Path counts) throws IOException{
		final byte[][] reference=new byte[REFERENCE_FILES.length+1][]; reference[0]=Files.readAllBytes(counts);
		for(int i=0; i<REFERENCE_FILES.length; i++){reference[i+1]=Files.readAllBytes(counts.getParent().resolve(REFERENCE_FILES[i]));}
		return reference;
	}
	private static long savedCalls(final byte[] data){
		final String[] rows=new String(data, StandardCharsets.US_ASCII).split("\n");
		require(rows.length==2 && rows[0].equals(COUNTS_HEADER), "Saved baseline counts have the wrong schema.");
		final String[] f=rows[1].split("\t", -1);
		require(f.length==9 && f[0].equals(RECIPE) && f[1].equals("baseline") && f[2].equals("pilot") && f[3].equals("64") &&
			f[5].equals("40") && f[6].equals("56") && f[7].equals("40"), "Saved counts do not describe this frozen workload.");
		final long calls=Long.parseLong(f[4]), cpu=Long.parseLong(f[8]);
		require(calls>0 && calls<=MAX_CALLS && calls%RING==0 && cpu>=MIN_PILOT_CPU, "Saved pilot must meet whole-ring bounds and the CPU minimum."); return calls;
	}
	private static long nextCalls(final long calls, final long cpu){
		require(cpu>0, "Calibration requires measured current-thread CPU time.");
		final long scaled=(long)Math.min(calls*16.0, Math.ceil(calls*(TARGET_CPU/(double)cpu)));
		final long next=((Math.max(calls+RING, scaled)+RING-1)/RING)*RING;
		require(next<=MAX_CALLS, "Pilot call cap exceeded; review sizing instead of silently expanding work."); return next;
	}
	/** Resolves a single unbound handle whose receiver has a stable Object type at invokeExact. */
	private static MethodHandle processHandle(){
		try{
			final Method method=Class.forName("jgi.BBMerge$MateThread").getDeclaredMethod("processReadPair", Read.class, Read.class);
			method.setAccessible(true);
			return MethodHandles.lookup().unreflect(method).asType(MethodType.methodType(int.class, Object.class, Read.class, Read.class));
		}catch(final ReflectiveOperationException e){throw new ExceptionInInitializerError(e);}
	}
	private static Field field(final Class<?> type, final String name) throws NoSuchFieldException{final Field f=type.getDeclaredField(name); f.setAccessible(true); return f;}
	private static Object value(final Object owner, final String name) throws ReflectiveOperationException{return field(BBMerge.class, name).get(owner);}
	private static String pin(final byte[] bytes) throws Exception{
		final byte[] digest=MessageDigest.getInstance("SHA-256").digest(bytes); final StringBuilder s=new StringBuilder(20);
		for(int i=digest.length-10; i<digest.length; i++){s.append(Character.forDigit((digest[i]>>>4)&15, 16)).append(Character.forDigit(digest[i]&15, 16));} return s.toString();
	}
	private static void require(final boolean value, final String why){if(!value){throw new IllegalArgumentException(why);}}
	private static Path path(final String s){require(s!=null && !s.isEmpty(), "Missing required path."); return Paths.get(s).toAbsolutePath().normalize();}
	private static byte[] bytes(final String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static void write(final Path p, final String s) throws IOException{Files.write(p, bytes(s), StandardOpenOption.CREATE_NEW);}
	private static String fastq(final String bases){final char[] q=new char[bases.length()]; Arrays.fill(q, '?'); return "@insert=40\n"+bases+"\n+\n"+new String(q)+"\n";}

	private static final class Pair{
		Pair(final int index){
			require(index>=0 && index<RING, "Only frozen ring slots are valid: "+index);
			numericID=index; positionA=index%40; positionB=index%56;
			abGolden=bytes(FRAGMENT); bbGolden=bytes(RAW_INPUT); aqGolden=new byte[40]; bqGolden=new byte[56];
			Arrays.fill(aqGolden, (byte)30); Arrays.fill(bqGolden, (byte)30);
			ab=abGolden.clone(); bb=bbGolden.clone(); aq=aqGolden.clone(); bq=bqGolden.clone();
			a=new Read(ab, aq, ID, index, false); b=new Read(bb, bq, ID, index, false); a.mate=b; b.mate=a; b.setPairnum(1);
			a.validate(true); b.validate(true); aFlags=a.flags; bFlags=b.flags; restore();
			require(stateMatches(false) && arraysUnchanged(), "Constructor/validation changed the literal pair's initial contract: "+index);
		}
		/** Reinstates references and scalar/auxiliary state without reconstructing arrays. */
		void restore(){reset(a, b, ab, aq, aFlags, 7); reset(b, a, bb, bq, bFlags, 9);}
		private void reset(final Read r, final Read mate, final byte[] bases, final byte[] quality, final int flags, final int errors){
			r.bases=bases; r.quality=quality; r.mate=mate; r.id=ID; r.numericID=numericID;
			r.chrom=-1; r.start=-1; r.stop=-1; r.copies=1; r.errors=errors; r.mapScore=0; r.neuralMapqs=0;
			r.match=null; r.gaps=null; r.sites=null; r.originalSite=null; r.obj=null; r.samline=null; r.rand=-1;
			r.setInsert(-1); r.flags=flags;
		}
		boolean arraysUnchanged(){return Arrays.equals(ab, abGolden) && Arrays.equals(bb, bbGolden) && Arrays.equals(aq, aqGolden) && Arrays.equals(bq, bqGolden);}
		boolean stateMatches(final boolean output){return state(a, b, ab, aq, aFlags, output ? 0 : 7) && state(b, a, bb, bq, bFlags, 9);}
		private boolean state(final Read r, final Read mate, final byte[] bases, final byte[] quality, final int flags, final int errors){
			return r.bases==bases && r.quality==quality && r.mate==mate && r.id==ID && r.numericID==numericID && r.flags==flags &&
				r.chrom==-1 && r.start==-1 && r.stop==-1 && r.copies==1 && r.errors==errors && r.mapScore==0 && r.neuralMapqs==0 &&
				r.match==null && r.gaps==null && r.sites==null && r.originalSite==null && r.obj==null && r.samline==null && r.rand==-1 && r.insert()==-1;
		}
		long checksum(final int returned){return returned+17L*a.length()+19L*b.length()+31L*a.bases[positionA]+131L*b.bases[positionB]+
			257L*a.quality[positionA]+65537L*b.quality[positionB]+a.errors+b.errors+a.flags+b.flags+
			((a.quality==aq && b.quality==bq) ? 104729L : 0L)+((a.bases==ab && b.bases==bb) ? 130363L : 0L);}
		long expectedChecksum(){return 40+17L*40+19L*56+31L*abGolden[positionA]+131L*bbGolden[positionB]+257L*30+65537L*30+9+aFlags+bFlags+104729L+130363L;}
		final Read a, b;
		final byte[] ab, bb, aq, bq, abGolden, bbGolden, aqGolden, bqGolden;
		final int numericID, positionA, positionB, aFlags, bFlags;
	}
	private static final class Result{
		Result(final long n, final long sum, final long w, final long c, final long p, final long b){calls=n; checksum=sum; wall=w; cpu=c; pairs=p; bases=b;}
		final long calls, checksum, wall, cpu, pairs, bases;
	}
	private final Path out;
	private final String phase;
	private final ThreadMXBean clock;
	private final Object worker;
	private final Field pairsField, basesField;
	private final Pair[] ring=new Pair[RING];
	private long cycleChecksum;
	private static final int RING=64, MAX_ATTEMPTS=6;
	private static final long INITIAL_CALLS=32768, MAX_CALLS=536870912;
	private static final long MIN_CPU=10000000000L, MIN_PILOT_CPU=15000000000L, TARGET_CPU=18000000000L;
	private static final String RECIPE="null-quality-worker-ring-v1", ID="insert=40";
	private static final String COUNTS_HEADER="recipe\tphase\tmode\tring\tcalls\tr1Length\tr2Length\tinsert\tpilotCpuNanos";
	private static final String[] REFERENCE_FILES={"ring.tsv", "oracles.tsv", "settings.tsv", "A.fq", "B.fq", "validation.tsv"};
	private static final String[] COMMON={"t=1", "simd=f", "bf1=t", "nn=f", "usequality=f", "qtrim=f", "maq=0", "minlength=1", "ordered=t", "ow=f", "app=f", "qin=33", "qout=33",
		"ratiomode=t", "flatmode=f", "entropy=f", "minoverlap=20", "minoverlap0=20", "ratiominoverlapreduction=0", "mininsert=1", "mininsert0=1",
		"ouq=f", "owq=t", "efilter=f", "pfilter=0", "mee=0", "tof=0", "extend=0", "extend2=0", "kfilter=0", "ecc=f", "join=t", "mix=f",
		"parsecustom=f", "samplerate=1", "sampleseed=1", "reads=-1", "recalibrate=f", "quantize=f", "train=f"};
	private static final String FRAGMENT="ACGTTGCATGTCAGTACGATCGTACCTGAGCTAGTCGATC";
	private static final String RAW_INPUT="GATCGACTAGCTCAGGTACGATCGTACTGACATGCAACGTAGATCGGAAGAGCACA";
	private static final MethodHandle PROCESS_PAIR=processHandle();
}
