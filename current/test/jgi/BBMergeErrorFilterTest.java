package test.jgi;

import java.io.IOException;
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
import java.util.List;
import java.util.Map;

import fileIO.FileFormat;
import jgi.BBMerge;
import jgi.BBMergeOverlapper;
import ml.CellNet;
import parse.Parse;
import stream.FASTQ;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import structures.FloatList;

/** Five independent worker-error-filter cases for JGI086, after the helper repair.
 * Baseline means the worker still ignores classical error-filter ambiguity; fixed
 * means only the strict non-neural mismatch is rejected. Every probe and CLI uses
 * a fresh JVM. Literal raw mates, fragment and quality oracles do not use joinRead
 * or the production reverse-complement routine to construct expected output.
 * The NN case deliberately uses test-only cutoff=0 and the unchanged real model:
 * it checks network reachability, not default-cutoff accuracy. Custom/training
 * semantics remain a source-review obligation and are not exercised here.
 * @author Sucrose
 */
public final class BBMergeErrorFilterTest{

	/** mode=suite|probe, phase=baseline|fixed, root, new outdir; fixed reference is
	 * the baseline FIXTURE output directory, containing summary.txt and case folders. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("mode", "phase", "root", "outdir", "reference", "case").contains(pair[0]), "Expected documented flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate fixture option: "+pair[0]);
		}
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use dev/bbmergeerrorfiltertest.sh with assertions enabled.");
		final String mode=options.containsKey("mode") ? options.get("mode") : "suite", phase=options.get("phase");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed for the worker-filter repair.");
		final Path root=path(options.get("root")), out=path(options.get("outdir"));
		if(mode.equals("probe")){
			require(!options.containsKey("reference") && Files.isDirectory(out), "Probe uses an existing case directory and no reference.");
			probe(root, out, spec(options.get("case")), phase);
			return;
		}
		require(mode.equals("suite") && !options.containsKey("case"), "Use mode=suite or an internal mode=probe case=name.");
		final Path reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require(phase.equals("baseline") ? reference==null : reference!=null, "Only fixed phase requires reference=<baseline fixture outdir>.");
		if(reference!=null){require(read(reference.resolve("summary.txt")).startsWith("BBMERGE_ERROR_FILTER_PASS phase=baseline cases=5 "), "Reference is not a successful five-case baseline fixture directory.");}
		final BBMergeErrorFilterTest test=new BBMergeErrorFilterTest(root, out, reference, phase);
		for(final String name : NAMES){test.run(spec(name));}
		test.check(modelPin(root).equals(MODEL_PIN), "Model bytes changed during the suite.");
		final String summary="BBMERGE_ERROR_FILTER_"+(test.failures==0 ? "PASS" : "FAIL")+" phase="+phase+" cases=5 children="+test.children+
			" checks="+test.checks+" failures="+test.failures+" customTraining=source-review-only defaultNnAccuracy=not-tested\n";
		write(out.resolve("summary.txt"), summary);
		System.out.print(summary);
		if(test.failures>0){throw new AssertionError(summary+"Inspect retained evidence; never weaken a failed reachability or output oracle.");}
	}

	/** Creates new evidence only; model and source/runtime dependencies are read-only. */
	private BBMergeErrorFilterTest(final Path root_, final Path out_, final Path reference_, final String phase_) throws Exception{
		root=root_; out=out_; reference=reference_; phase=phase_;
		Files.createDirectory(out);
		require(modelPin(root).equals(MODEL_PIN), "Unexpected model bytes; no model substitution or retraining.");
		write(out.resolve("model.sha80"), MODEL_PIN+"\n");
		write(out.resolve("cases.tsv"), "case\tprobeExit\tcliExit\texpectedMerged\tmerged\tunmerged1\tunmerged2\toracleMatches\n");
		if(reference!=null){same(out.resolve("model.sha80"), reference.resolve("model.sha80"), "Model bytes differ from baseline.");}
	}

	/** Retains every child's command, exit and output; root owns execution monitoring. */
	private void run(final Spec spec) throws Exception{
		final Path work=Files.createDirectory(out.resolve(spec.name));
		writeFastq(work.resolve("A.fq"), spec.first()); writeFastq(work.resolve("B.fq"), spec.second());
		final byte[] inputA=Files.readAllBytes(work.resolve("A.fq")), inputB=Files.readAllBytes(work.resolve("B.fq"));
		if(reference!=null){
			same(work.resolve("A.fq"), reference.resolve(spec.name).resolve("A.fq"), "Raw mate1 changed: "+spec.name);
			same(work.resolve("B.fq"), reference.resolve(spec.name).resolve("B.fq"), "Raw mate2 changed: "+spec.name);
		}
		final int probeExit=child(new ArrayList<String>(Arrays.asList("bash", root.resolve("dev/bbmergeerrorfiltertest.sh").toString(),
			"mode=probe", "phase="+phase, "root="+root, "outdir="+work, "case="+spec.name)), work, "probe");
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("bash", root.resolve("bbmerge.sh").toString(), "-Xmx128m", "-Xms32m"));
		command.addAll(arguments(root, work, spec));
		final int cliExit=child(command, work, "cli");
		check(probeExit==0, spec.name+": actual candidate/gate/worker prerequisite failed.");
		check(cliExit==0, spec.name+": real CLI failed.");
		check(Arrays.equals(inputA, Files.readAllBytes(work.resolve("A.fq"))) && Arrays.equals(inputB, Files.readAllBytes(work.resolve("B.fq"))), spec.name+": child changed input FASTQs.");
		if(probeExit!=0 || cliExit!=0){append(out.resolve("cases.tsv"), spec.name+"\t"+probeExit+"\t"+cliExit+"\tNA\tNA\tNA\tNA\tfalse\n"); return;}
		final List<Read> merged=reads(work.resolve("merged.fq")), first=reads(work.resolve("unmerged1.fq")), second=reads(work.resolve("unmerged2.fq"));
		final boolean merge=spec.merge(phase), matches=merge ? one(merged, spec.joined()) && first.isEmpty() && second.isEmpty() :
			merged.isEmpty() && one(first, spec.first()) && one(second, spec.second());
		check(matches, spec.name+": routing, IDs, length, full sequence or numeric qualities differ from independent oracle.");
		write(work.resolve("outputs.tsv"), records("merged", merged)+records("unmerged1", first)+records("unmerged2", second));
		append(out.resolve("cases.tsv"), spec.name+"\t"+probeExit+"\t"+cliExit+"\t"+(merge ? 1 : 0)+"\t"+merged.size()+"\t"+first.size()+"\t"+second.size()+"\t"+matches+"\n");
		if(reference!=null){
			final Path before=reference.resolve(spec.name);
			if(!spec.name.equals("strict")){
				for(final String file : new String[]{"merged.fq", "unmerged1.fq", "unmerged2.fq", "outputs.tsv"}){same(work.resolve(file), before.resolve(file), "Unchanged control output differs: "+spec.name+"/"+file);}
			}
			final Map<String, String> old=table(before.resolve("probe.tsv")), now=table(work.resolve("probe.tsv"));
			for(final String key : old.keySet()){
				if(!key.equals("phase") && !(spec.name.equals("strict") && key.equals("workerReturned"))){check(old.get(key).equals(now.get(key)), spec.name+": unexpected probe-state change: "+key);}
			}
			check(old.size()==now.size(), spec.name+": probe field inventory changed.");
			if(spec.nn){same(work.resolve("vector.tsv"), before.resolve("vector.tsv"), "All23 NN inputs/raw bits must remain identical.");}
		}
	}

	/** Calls the real worker once, never starting a thread or adding a production hook.
	 * For non-NN cases a separate ratio-only call first exposes the known candidate;
	 * NN uses the actual worker-produced vector and reads its cached activation once. */
	private static void probe(final Path root, final Path work, final Spec spec, final String phase) throws Exception{
		require(modelPin(root).equals(MODEL_PIN), "Unexpected probe model bytes.");
		final BBMerge merge=new BBMerge(arguments(root, work, spec).toArray(new String[0]));
		require(bool(merge, "MAKE_VECTOR")==spec.nn && !bool(merge, "TAG_CUSTOM") && !bool(merge, "parseCustom"), "NN/custom/training mode differs from this bounded fixture.");
		require(bool(merge, "useQuality") && !bool(merge, "overlapUsingQuality") && bool(merge, "overlapWithoutQuality") && bool(merge, "useRatioMode") && !bool(merge, "useFlatMode"), "Expected quality-filtered unweighted ratio candidate path.");
		require(!bool(merge, "useMEEfilter") && !bool(merge, "qtrim1") && !bool(merge, "qtrim2") && !bool(merge, "useEntropy"), "A prefilter or alternate path hides the error-filter gate.");
		require(bool(merge, "useEfilter")==spec.enabled && number(merge, "efilterRatio")==spec.ratio && number(merge, "efilterOffset")==.05f,
			"Actual error-filter boolean/ratio/offset differs from planned case.");
		require(number(merge, "pfilterRatio")==spec.pfilter, "Probability filter differs from the isolated case.");
		final Class<?> type=Class.forName("jgi.BBMerge$MateThread");
		final Constructor<?> constructor=type.getDeclaredConstructor(BBMerge.class, Streamer.class, Writer.class, Writer.class, fileIO.ByteStreamWriter.class, boolean.class, boolean.class);
		constructor.setAccessible(true);
		final Object worker=constructor.newInstance(merge, null, null, null, null, true, false);
		require(((Thread)worker).getState()==Thread.State.NEW, "A probe must not start its worker.");
		int candidate=-1, bad=-1, ambiguous=-1;
		if(!spec.nn){
			final Read a=spec.first(), b=spec.alignedSecond();
			final int min=(Integer)invoke(worker, "calcMinOverlapFromEntropy", new Class<?>[]{Read.class, Read.class}, a, b);
			candidate=(Integer)invoke(worker, "mateByOverlap_ratioMode", new Class<?>[]{Read.class, Read.class, int.class}, a, b, min);
			final int[] rv=(int[])value(worker, "rvector"); bad=rv[2]; ambiguous=rv[4];
		}
		final Read a=spec.first(), b=spec.alignedSecond();
		a.mate=b; b.mate=a; b.setPairnum(1);
		final byte[] ab=a.bases.clone(), aq=a.quality.clone(), bb=b.bases.clone(), bq=b.quality.clone();
		final int returned=(Integer)invoke(worker, "mateByOverlap", new Class<?>[]{Read.class, Read.class}, a, b);
		final float expected=BBMergeOverlapper.expectedMismatches(a, b, 100), probability=BBMergeOverlapper.probability(a, b, 100);
		final boolean classicalReject=spec.enabled && (expected+.05f)*spec.ratio<spec.bad;
		final Map<String, String> evidence=new LinkedHashMap<String, String>();
		evidence.put("case", spec.name); evidence.put("phase", phase);
		evidence.put("useEfilter", Boolean.toString(spec.enabled)); evidence.put("efilterRatio", Float.toString(spec.ratio));
		evidence.put("efilterOffset", "0.05"); evidence.put("pfilter", Float.toString(spec.pfilter));
		evidence.put("helperExpected", Float.toString(expected)); evidence.put("probability", Float.toString(probability));
		evidence.put("classicalWouldReject", Boolean.toString(classicalReject)); evidence.put("workerReturned", Integer.toString(returned));
		evidence.put("probeScope", "one-mateByOverlap-call; separate-CLI-retries-not-instrumented");
		if(spec.nn){
			final FloatList vector=b.fetchVector();
			require(vector!=null && vector.size==23, "REACHABILITY_FAIL: NN case must produce all23 actual inputs.");
			final StringBuilder text=new StringBuilder("input\tfloat\trawBitsHex\n");
			for(int i=0; i<23; i++){
				require(Float.isFinite(vector.get(i)), "NN input is non-finite at"+i);
				text.append(i).append('\t').append(vector.get(i)).append('\t').append(Integer.toHexString(Float.floatToRawIntBits(vector.get(i)))).append('\n');
			}
			write(work.resolve("vector.tsv"), text.toString());
			candidate=Math.round(vector.get(5)/.004f);
			require(Float.floatToRawIntBits(vector.get(5))==Float.floatToRawIntBits(100*.004f), "NN input5 must encode actual candidate100, not a supplied result.");
			require(Float.floatToRawIntBits(vector.get(21))==Float.floatToRawIntBits(expected) &&
				Float.floatToRawIntBits(vector.get(22))==Float.floatToRawIntBits((float)(Math.sqrt(probability)+.0000015)), "Actual NN helper/probability inputs differ from their documented units.");
			final int[] rv=(int[])value(worker, "rvector"); bad=rv[2]; ambiguous=rv[4];
			final CellNet net=(CellNet)value(worker, "net"); require(net!=null, "Real worker network is missing.");
			final float score=net.getOutput(0), cutoff=number(merge, "netCutoff");
			require(Float.isFinite(score) && cutoff==0 && score>=cutoff, "Actual cached NN output must reach the explicit test-only zero cutoff.");
			evidence.put("score", Float.toString(score)); evidence.put("scoreBits", Integer.toHexString(Float.floatToRawIntBits(score)));
			evidence.put("cutoff", Float.toString(cutoff)); evidence.put("defaultAccuracy", "not-tested; explicit-test-only-cutoff0");
		}
		evidence.put("candidateInsert", Integer.toString(candidate)); evidence.put("candidateBad", Integer.toString(bad)); evidence.put("candidateAmbiguous", Integer.toString(ambiguous));
		final StringBuilder text=new StringBuilder();
		for(final Map.Entry<String, String> entry : evidence.entrySet()){text.append(entry.getKey()).append('\t').append(entry.getValue()).append('\n');}
		write(work.resolve("probe.tsv"), text.toString());
		require(candidate==100 && bad==spec.bad && ambiguous==0, "REACHABILITY_FAIL: ratio candidate must be insert100 with the literal mismatch count, before error filtering.");
		// Perfect:100*(1-.99*.99). Mismatch case replaces one term with1-.999*.99;
		// changing the base does not change this quality-only expectation.
		require(Float.isFinite(expected) && Math.abs(expected-(spec.mismatch ? 1.98109 : 1.99))<.0005, "Helper differs from independent decimal quality oracle.");
		require(classicalReject==(spec.enabled && spec.ratio==.06f && spec.mismatch), "Independent strict/perfect/permissive gate inequality differs.");
		require(probability>=spec.pfilter, "Probability filter hides the intended classical threshold.");
		require(returned==(spec.merge(phase) ? 100 : BBMerge.RET_AMBIG), "Worker return differs from phase oracle; an unconditional ambiguity return must fail the NN control.");
		require(Arrays.equals(ab, a.bases) && Arrays.equals(aq, a.quality) && Arrays.equals(bb, b.bases) && Arrays.equals(bq, b.quality), "Worker probe changed input bases or qualities.");
	}

	/** Shared CLI/constructor options isolate efilter; NN alone keeps positive pfilter
	 * so vector22 contains actual probability, and uses an explicitly test-only cutoff. */
	private static ArrayList<String> arguments(final Path root, final Path work, final Spec spec){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in1="+work.resolve("A.fq"), "in2="+work.resolve("B.fq"),
			"out="+work.resolve("merged.fq"), "outu1="+work.resolve("unmerged1.fq"), "outu2="+work.resolve("unmerged2.fq"),
			"t=1", "simd=f", "qin=33", "qout=33", "ow=t", "qtrim=f", "maq=0", "minlength=1", "usequality=t", "ouq=f", "owq=t", "ratiomode=t", "flatmode=f",
			"entropy=f", "minoverlap=20", "minoverlap0=20", "ratiominoverlapreduction=0", "mininsert=1", "mininsert0=1",
			"efilter="+(spec.enabled ? Float.toString(spec.ratio) : "f"), "efilteroffset=0.05", "pfilter="+spec.pfilter,
			"mee=0", "tof=0", "extend=0", "extend2=0", "kfilter=0", "nn="+spec.nn));
		if(spec.nn){args.add("net="+root.resolve("networks/bbmerge.bbnet.gz")); args.add("cutoff=0");}
		return args;
	}

	/** Owns exactly the returned Process; no process-name searches or signaling. */
	private int child(final ArrayList<String> command, final Path work, final String name) throws Exception{
		children++;
		final StringBuilder saved=new StringBuilder();
		for(int i=0; i<command.size(); i++){saved.append(i).append('\t').append(command.get(i)).append('\n');}
		write(work.resolve(name+".command.tsv"), saved.toString());
		final ProcessBuilder builder=new ProcessBuilder(command).directory(root.toFile());
		builder.redirectOutput(work.resolve(name+".stdout").toFile()); builder.redirectError(work.resolve(name+".stderr").toFile());
		builder.environment().remove("_JAVA_OPTIONS"); builder.environment().remove("JDK_JAVA_OPTIONS");
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final int code=builder.start().waitFor();
		write(work.resolve(name+".exit"), code+"\n");
		return code;
	}

	/** Native FASTQ reader; single files cannot infer interleaving or ASCII offset. */
	private static List<Read> reads(final Path path){
		require(Files.isRegularFile(path), "Missing output file: "+path);
		FASTQ.ASCII_OFFSET=33; FASTQ.DETECT_QUALITY=false;
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		return StreamerFactory.getReads(-1, false, FileFormat.testInput(path.toString(), FileFormat.FASTQ, null, true, false), null, null, null);
	}
	/** Exact cardinality, ID, orientation, sequence and every numeric quality. */
	private static boolean one(final List<Read> actual, final Read expected){
		if(actual.size()!=1){return false;}
		final Read read=actual.get(0);
		return read.mate==null && read.id.equals(expected.id) && Arrays.equals(read.bases, expected.bases) && Arrays.equals(read.quality, expected.quality);
	}
	private static String records(final String channel, final List<Read> reads){
		final StringBuilder text=new StringBuilder();
		for(final Read read : reads){text.append(channel).append('\t').append(read.id).append('\t').append(new String(read.bases, StandardCharsets.US_ASCII)).append('\t').append(Arrays.toString(read.quality)).append('\n');}
		return text.toString();
	}
	/** Literal ASCII33 writer, independent of the output formatter under test. */
	private static void writeFastq(final Path path, final Read read) throws IOException{
		final StringBuilder text=new StringBuilder("@").append(read.id).append('\n').append(new String(read.bases, StandardCharsets.US_ASCII)).append("\n+\n");
		for(final byte q : read.quality){text.append((char)(q+33));}
		write(path, text.append('\n').toString());
	}

	/** Five fixed cases; only strict non-NN changes its expected routing by phase. */
	private static Spec spec(final String name){
		require(name!=null && Arrays.asList(NAMES).contains(name), "Unknown worker error-filter case: "+name);
		return new Spec(name);
	}
	private static final class Spec{
		Spec(final String name_){
			name=name_; mismatch=!name.equals("perfect"); nn=name.equals("nn"); enabled=!name.equals("disabled");
			ratio=enabled ? (name.equals("permissive") ? 6f : .06f) : -1f;
			pfilter=nn ? .00004f : 0f; bad=mismatch ? 1 : 0;
			require(F100.length()==100 && RAW_F100.length()==100 && ADAPTER_A.length()==60 && ADAPTER_B_RAW.length()==60 && ADAPTER_B_RC.length()==60,
				"Literal fragment/raw mate/adapter lengths must preserve insert100 and read160 geometry.");
			require(F100.charAt(50)=='G' && RAW_F100.charAt(49)=='C', "Mismatch must be aligned fragment50 G/T, raw mate2 site49 C/A.");
		}
		boolean merge(final String phase){return !phase.equals("fixed") || !name.equals("strict");}
		Read first(){
			final byte[] q=qualities(160, 20); Arrays.fill(q, 100, 160, (byte)40); if(mismatch){q[50]=30;}
			return new Read((F100+ADAPTER_A).getBytes(StandardCharsets.US_ASCII), q, name+"/1", 0, false);
		}
		Read second(){
			final char[] raw=(RAW_F100+ADAPTER_B_RAW).toCharArray(); if(mismatch){raw[49]='A';}
			final byte[] q=qualities(160, 20); Arrays.fill(q, 100, 160, (byte)40);
			return new Read(new String(raw).getBytes(StandardCharsets.US_ASCII), q, name+"/2", 0, false);
		}
		/** Supplies aligned probe input from separate literals, not production reversal. */
		Read alignedSecond(){
			final char[] aligned=(ADAPTER_B_RC+F100).toCharArray(); if(mismatch){aligned[110]='T';}
			final byte[] q=qualities(160, 20); Arrays.fill(q, 0, 60, (byte)40);
			return new Read(new String(aligned).getBytes(StandardCharsets.US_ASCII), q, name+"/2", 0, false);
		}
		Read joined(){final byte[] q=qualities(100, 25); if(mismatch){q[50]=10;} return new Read(F100.getBytes(StandardCharsets.US_ASCII), q, name+"/1", 0, false);}
		final String name;
		final boolean mismatch, nn, enabled;
		final float ratio, pfilter;
		final int bad;
	}

	private static byte[] qualities(final int length, final int quality){final byte[] q=new byte[length]; Arrays.fill(q, (byte)quality); return q;}
	/** Reflection is confined to test observation/invocation; no state is overwritten. */
	private static Object value(final Object owner, final String name) throws Exception{final Field f=owner.getClass().getDeclaredField(name); f.setAccessible(true); return f.get(owner);}
	private static boolean bool(final Object owner, final String name) throws Exception{return (Boolean)value(owner, name);}
	private static float number(final Object owner, final String name) throws Exception{return ((Number)value(owner, name)).floatValue();}
	private static Object invoke(final Object owner, final String name, final Class<?>[] types, final Object... args) throws Exception{final Method m=owner.getClass().getDeclaredMethod(name, types); m.setAccessible(true); return m.invoke(owner, args);}
	private void same(final Path a, final Path b, final String message) throws IOException{check(Arrays.equals(Files.readAllBytes(a), Files.readAllBytes(b)), message);}
	private void check(final boolean condition, final String message){checks++; if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}}
	private static void require(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
	private static Path path(final String value){require(value!=null && !value.isEmpty(), "Missing required path."); return Paths.get(value).toAbsolutePath().normalize();}
	private static void write(final Path path, final String text) throws IOException{Files.write(path, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path path, final String text) throws IOException{Files.write(path, text.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);}
	private static String read(final Path path) throws IOException{return new String(Files.readAllBytes(path), StandardCharsets.UTF_8);}
	private static Map<String, String> table(final Path path) throws IOException{
		final Map<String, String> map=new LinkedHashMap<String, String>();
		for(final String line : read(path).split("\n")){final String[] pair=Parse.splitOnFirst(line, '\t'); require(pair.length==2 && map.put(pair[0], pair[1])==null, "Malformed or duplicate probe field.");}
		return map;
	}
	/** Final80bits of the compressed model file; never exposes the full digest. */
	private static String modelPin(final Path root) throws Exception{
		final byte[] bytes=MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(root.resolve("networks/bbmerge.bbnet.gz")));
		final StringBuilder pin=new StringBuilder();
		for(int i=bytes.length-10; i<bytes.length; i++){pin.append(Character.forDigit((bytes[i]>>>4)&15, 16)).append(Character.forDigit(bytes[i]&15, 16));}
		return pin.toString();
	}

	private final Path root, out, reference;
	private final String phase;
	private int checks, failures, children;
	private static final String[] NAMES={"strict", "perfect", "disabled", "permissive", "nn"};
	private static final String MODEL_PIN="24998a8aa257552d4a22";
	private static final String F100="ACGTTGCATGTCGACATGCTAGCTTCGATCGATGACCTGATCCGATAGCA"+"GTACCGTATCATGCTAGGTCCGATCTACGATACGGTTCAGCAGTACGATG";
	private static final String RAW_F100="CATCGTACTGCTGAACCGTATCGTAGATCGGACCTAGCATGATACGGTACTGCTATCGGATCAGGTCATCGATCGAAGCTAGCATGTCGACATGCAACGT";
	/* Recognized60-base tails from the accepted equal-length witness. Raw mate2 ends
	 * in ADAPTER_B_RAW; its aligned probe begins with independently literal B_RC. */
	private static final String ADAPTER_A="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTG";
	private static final String ADAPTER_B_RAW="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGC";
	private static final String ADAPTER_B_RC="GCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC";
}
