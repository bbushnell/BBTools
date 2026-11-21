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

import bbduk.BBDukS;
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

/** Fifteen finite probability caller families using four predetermined paired inputs.
 * Reuses the retained expected-error caller fixture's fresh-process, observation and
 * native FASTQ comparison pattern without invoking its earlier ten-family suite.
 * Covers current BBDukS/ProcessorS, classical BBMerge worker/static calls, and three
 * real-model packets. Ordinary and disabled-filter controls retain baseline bits;
 * changed neural scores/decisions are reported without direction or accuracy claims.
 * Probes invoke actual code on owned arrays; reflection never writes production state.
 * Probes start no read-processing worker or I/O pipeline. Real index construction
 * may start and join allocation helpers. CLI children run independently.
 * @author Sucrose
 */
public final class BBMergeProbabilityCallerTest{

	/** mode=suite|probe; phase=baseline|fixed; root; new outdir; fixed-only reference. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("mode", "phase", "root", "outdir", "reference", "family").contains(pair[0]), "Expected documented flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate fixture option: "+pair[0]);
		}
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbmergeprobabilitycallertest.sh with assertions enabled.");
		final String mode=options.containsKey("mode") ? options.get("mode") : "suite", phase=options.get("phase");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		final Path root=path(options.get("root")), out=path(options.get("outdir"));
		require(Files.isRegularFile(root.resolve("bbmerge.sh")) && Files.isRegularFile(root.resolve("bbduk.sh")), "Set root to the single runtime.");
		if(mode.equals("probe")){
			require(!options.containsKey("reference") && Files.isDirectory(out), "Probe uses an existing family directory and no reference.");
			probe(root, out, family(options.get("family")), phase); return;
		}
		require(mode.equals("suite") && !options.containsKey("family"), "Use suite or an internal probe family.");
		final Path reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require(phase.equals("baseline") ? reference==null : reference!=null, "Only fixed requires reference=BASELINE_CALLER_FIXTURE_DIRECTORY.");
		if(reference!=null){require(read(reference.resolve("summary.txt")).matches("BBMERGE_PROBABILITY_CALLER_PASS phase=baseline version=1 families=15 probes=15 cli=11 children=26 checks=[0-9]+ failures=0\\n"), "Reference must be a successful fifteen-family probability baseline.");}
		final BBMergeProbabilityCallerTest test=new BBMergeProbabilityCallerTest(root, out, reference, phase);
		for(final String name : NAMES){test.run(family(name));}
		test.check(test.families==15 && test.probes==15 && test.cli==11 && test.children==26, "Expected15 probes and11 CLIs, without replaying old families.");
		test.check(modelPin(root).equals(MODEL_PIN), "Model bytes changed during suite.");
		final String summary="BBMERGE_PROBABILITY_CALLER_"+(test.failures==0 ? "PASS" : "FAIL")+" phase="+phase+
			" version=1 families="+test.families+" probes="+test.probes+" cli="+test.cli+" children="+test.children+" checks="+test.checks+" failures="+test.failures+"\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
		if(test.failures!=0){throw new AssertionError(summary+"Retain prerequisite failures; do not fit inputs or oracles to results.");}
	}

	private BBMergeProbabilityCallerTest(final Path root_, final Path out_, final Path reference_, final String phase_) throws Exception{
		root=root_; out=out_; reference=reference_; phase=phase_;
		require(modelPin(root).equals(MODEL_PIN), "Model substitution/retraining is outside this fixture.");
		Files.createDirectory(out);
		write(out.resolve("model.sha80"), MODEL_PIN+"\n");
		write(out.resolve("fragment1024.txt"), fragment(1024)+"\n");
		write(out.resolve("cases.tsv"), "family\tprobeExit\tcliExit\taccepted\toracleMatches\n");
		write(out.resolve("decisions.tsv"), "family\tbaselineAccepted\tfixedAccepted\tbaselineScoreBits\tfixedScoreBits\tfirstCallVsCliScope\n");
		if(reference!=null){same(out.resolve("model.sha80"), reference.resolve("model.sha80"), "Model changed."); same(out.resolve("fragment1024.txt"), reference.resolve("fragment1024.txt"), "Predetermined fragment changed.");}
	}

	/** A failed probe skips its dependent CLI; independent families still run. */
	private void run(final Family f) throws Exception{
		families++;
		final Path work=Files.createDirectory(out.resolve(f.name));
		writeFastq(work.resolve("A.fq"), f.spec.first()); writeFastq(work.resolve("B.fq"), f.spec.second());
		write(work.resolve("inputs.tsv"), "recipe\t"+RECIPE+"\n"+records("A", Arrays.asList(f.spec.first()))+records("B-raw", Arrays.asList(f.spec.second()))+records("B-aligned", Arrays.asList(f.spec.alignedSecond())));
		verifyAdapters(work, f.spec);
		final byte[] inputA=Files.readAllBytes(work.resolve("A.fq")), inputB=Files.readAllBytes(work.resolve("B.fq"));
		if(reference!=null){for(final String file : new String[]{"A.fq", "B.fq", "inputs.tsv"}){same(work.resolve(file), reference.resolve(f.name).resolve(file), f.name+": input changed: "+file);}}
		final int probeExit=child(new ArrayList<String>(Arrays.asList("bash", root.resolve("dev/bbmergeprobabilitycallertest.sh").toString(),
			"mode=probe", "phase="+phase, "root="+root, "outdir="+work, "family="+f.name)), work, "probe");
		probes++;
		check(probeExit==0, f.name+": prerequisite/probe failed; inspect probe.stderr.");
		int cliExit=-1;
		if(!f.isStatic() && probeExit==0){
			final ArrayList<String> command=new ArrayList<String>(Arrays.asList("bash", root.resolve(f.duk() ? "bbduk.sh" : "bbmerge.sh").toString(), "-Xmx128m", "-Xms32m"));
			command.addAll(arguments(root, work, f)); cliExit=child(command, work, "cli");
			cli++;
			check(cliExit==0, f.name+": actual paired CLI failed.");
		}
		check(Arrays.equals(inputA, Files.readAllBytes(work.resolve("A.fq"))) && Arrays.equals(inputB, Files.readAllBytes(work.resolve("B.fq"))), f.name+": probe/CLI changed input files.");
		final String cliStatus=f.isStatic() ? "NA" : probeExit!=0 ? "SKIPPED" : Integer.toString(cliExit);
		if(probeExit!=0 || (!f.isStatic() && cliExit!=0)){
			append(out.resolve("cases.tsv"), f.name+"\t"+probeExit+"\t"+cliStatus+"\tNA\tfalse\n"); return;
		}
		final Map<String, String> packet=table(work.resolve("probe.tsv"));
		boolean accepted=Boolean.parseBoolean(packet.get("accepted")), truth=true;
		if(f.duk()){
			final List<Read> a=reads(work.resolve("trimmed1.fq")), b=reads(work.resolve("trimmed2.fq"));
			final int length=f.accept(phase) ? f.spec.length : f.spec.length+60;
			truth=one(a, f.spec.first(), length) && one(b, f.spec.second(), length);
			check(truth, f.name+": trimming must preserve literal raw prefixes/IDs/qualities.");
			write(work.resolve("outputs.tsv"), records("mate1", a)+records("mate2", b));
		}else if(!f.isStatic()){
			final List<Read> merged=reads(work.resolve("merged.fq")), a=reads(work.resolve("unmerged1.fq")), b=reads(work.resolve("unmerged2.fq"));
			check(merged.size()<=1 && a.size()<=1 && b.size()<=1 && merged.size()+a.size()==1 && a.size()==b.size(), f.name+": duplicated/missing pair routing.");
			accepted=merged.size()==1;
			truth=accepted ? one(merged, f.spec.joined(), f.spec.length) && a.isEmpty() && b.isEmpty() : merged.isEmpty() && one(a, f.spec.first(), f.spec.length+60) && one(b, f.spec.second(), f.spec.length+60);
			check(truth, f.name+": CLI records differ from the phase-local literal consensus/raw-pair oracle.");
			if(!f.neural()){check(accepted==f.accept(phase), f.name+": classical CLI cutoff routing differs from oracle.");}
			write(work.resolve("outputs.tsv"), records("merged", merged)+records("mate1", a)+records("mate2", b));
		}
		append(out.resolve("cases.tsv"), f.name+"\t0\t"+cliStatus+"\t"+accepted+"\t"+truth+"\n");
		if(reference!=null && f.neural()){comparePacket(f, work, packet, accepted);}
	}

	/** Only probability input22 may change for active last; controls preserve all bits. */
	private void comparePacket(final Family f, final Path work, final Map<String, String> now, final boolean accepted) throws Exception{
		final Path before=reference.resolve(f.name);
		final String[] a=read(before.resolve("vector.tsv")).split("\n"), b=read(work.resolve("vector.tsv")).split("\n");
		check(a.length==24 && b.length==24, f.name+": expected header and23 vector rows.");
		final boolean control=f.disabled() || f.spec.name.equals("ordinary");
		if(a.length==24 && b.length==24){for(int i=0; i<23; i++){if(control || i!=22){check(a[i+1].equals(b[i+1]), f.name+": unexpected vector raw-bit change at input"+i);}}}
		final Map<String, String> old=table(before.resolve("probe.tsv"));
		check(old.get("cutoff").equals(now.get("cutoff")), f.name+": default network cutoff changed.");
		if(control){
			check(old.get("scoreBits").equals(now.get("scoreBits")) && old.get("returned").equals(now.get("returned")), f.name+": control network score/first-call return changed.");
			for(final String file : new String[]{"merged.fq", "unmerged1.fq", "unmerged2.fq"}){same(work.resolve(file), before.resolve(file), f.name+": control output bytes changed.");}
		}
		append(out.resolve("decisions.tsv"), f.name+"\t"+(!reads(before.resolve("merged.fq")).isEmpty())+"\t"+accepted+"\t"+old.get("scoreBits")+"\t"+now.get("scoreBits")+"\tfirst-candidate packet; CLI retry count not instrumented\n");
	}

	/** Actual settings and candidates precede any asserted probability-filter outcome. */
	private static void probe(final Path root, final Path work, final Family f, final String phase) throws Exception{
		require(modelPin(root).equals(MODEL_PIN), "Unexpected compressed model bytes.");
		final Map<String, String> evidence=new LinkedHashMap<String, String>();
		evidence.put("family", f.name); evidence.put("phase", phase); evidence.put("recipe", RECIPE);
		final Read a=f.spec.first(), b=f.spec.alignedSecond();
		link(a, b);
		final Snapshot snapshot=new Snapshot(a, b);
		final float probability=BBMergeOverlapper.probability(a, b, f.spec.length);
		require(close(probability, f.spec.probability(phase)), "Direct helper score differs from independent phase oracle: "+f.name+" "+probability);
		evidence.put("probability", Float.toString(probability)); evidence.put("probabilityBits", bits(probability));
		evidence.put("probabilityOracle", Double.toString(f.spec.probability(phase)));
		if(f.duk()){probeDuk(root, work, f, phase, evidence);}
		else if(f.isStatic()){probeStatic(f, phase, a, b, snapshot, evidence);}
		else{probeWorker(root, work, f, phase, a, b, snapshot, evidence);}
		final StringBuilder text=new StringBuilder();
		for(final Map.Entry<String, String> entry : evidence.entrySet()){text.append(entry.getKey()).append('\t').append(entry.getValue()).append('\n');}
		write(work.resolve("probe.tsv"), text.toString());
	}

	/** Current ProcessorS processes an owned pair without starting its stream loop. */
	private static void probeDuk(final Path root, final Path work, final Family f, final String phase, final Map<String, String> evidence) throws Exception{
		final BBDukS driver=new BBDukS(arguments(root, work, f).toArray(new String[0]));
		final Object parser=value(driver, "parser"), index=value(driver, "index");
		final Class<?> parserType=Class.forName("bbduk.BBDukParser"), indexType=Class.forName("bbduk.BBDukIndex"), type=Class.forName("bbduk.BBDukProcessorS");
		final Constructor<?> ctor=type.getDeclaredConstructor(parserType, indexType, Streamer.class, Writer.class, Writer.class, Writer.class);
		ctor.setAccessible(true);
		final Object processor=ctor.newInstance(parser, index, null, null, null, null);
		require(bool(processor, "trimByOverlap") && bool(processor, "useQualityForOverlap") && number(processor, "pfilterRatio")==.00005f &&
			number(processor, "efilterRatio")==6 && number(processor, "efilterOffset")==.05f && number(processor, "meeFilter")==999999999f,
			"Current permissive TBO/quality/filter settings were not reached.");
		require(value(processor, "ref")==null && value(processor, "literal")==null && longNumber(processor, "storedKmers")==0 && value(processor, "sidechannel")==null,
			"Processor probe must use its real reference-free index and no mapping sidechannel.");
		require(!bool(processor, "ecc") && !bool(processor, "quantizeQuality") && !bool(processor, "qtrimLeft") && !bool(processor, "qtrimRight") && number(processor, "minAvgQuality")==0,
			"Unexpected correction, quantization or quality trimming/filtering.");
		final Read ca=f.spec.first(), cb=f.spec.alignedSecond();
		final int[] rv=new int[5];
		final int candidate=BBMergeOverlapper.mateByOverlapRatio(ca, cb, new float[ca.length()], new float[cb.length()], rv,
			integer(processor, "minOverlap0"), integer(processor, "minOverlap"), integer(processor, "minInsert0"), integer(processor, "minInsert"),
			number(processor, "maxRatio"), .12f, number(processor, "ratioMargin"), number(processor, "ratioOffset"), .95f, .95f, true);
		candidate(f, candidate, rv);
		final Read a=f.spec.first(), b=f.spec.second(); link(a, b);
		final float pre=Math.max(a.expectedErrors(false, -1), b.expectedErrors(false, -1));
		final float expected=BBMergeOverlapper.expectedMismatches(ca, cb, candidate);
		require(pre<number(processor, "meeFilter") && expected<=number(processor, "meeFilter") && (expected+.05f)*6>=f.spec.mismatches,
			"REACHABILITY_FAIL: MEE or efilter hides the probability gate.");
		final ArrayList<Read> list=new ArrayList<Read>(); list.add(a);
		final int removed=(Integer)invoke(processor, "processList", new Class<?>[]{ArrayList.class, ArrayList.class, ArrayList.class, ArrayList.class},
			list, new ArrayList<Read>(), new ArrayList<Read>(), new ArrayList<Read>());
		final int[] actual=(int[])value(processor, "overlapVector");
		require(removed==0 && list.size()==1 && list.get(0)==a && actual[2]==f.spec.mismatches && actual[4]==0, "REACHABILITY_FAIL: actual ProcessorS rejected/changed the candidate state.");
		final boolean accepted=a.length()==f.spec.length && b.length()==f.spec.length;
		require(accepted==f.accept(phase), "Actual ProcessorS trim decision differs from phase oracle.");
		final int length=accepted ? f.spec.length : f.spec.length+60;
		require(rawPrefix(a, f.spec.first(), length) && rawPrefix(b, f.spec.second(), length), "Actual ProcessorS trim result differs from literal raw-prefix oracle.");
		require(longNumber(processor, "readsTrimmedByOverlap")==(accepted ? 2 : 0) && longNumber(processor, "basesTrimmedByOverlap")==(accepted ? 120 : 0), "ProcessorS overlap trim counters disagree.");
		evidence.put("candidateInsert", Integer.toString(candidate)); evidence.put("candidateSource", "standalone actual finder; actual processor rv[0] is not an insert observation");
		evidence.put("actualBad", Integer.toString(actual[2])); evidence.put("actualAmbiguous", Integer.toString(actual[4]));
		evidence.put("preMaximumFullReadErrors", Float.toString(pre)); evidence.put("helperExpected", Float.toString(expected));
		evidence.put("pfilter", Float.toString(number(processor, "pfilterRatio"))); evidence.put("efilter", Float.toString(number(processor, "efilterRatio")));
		evidence.put("mee", Float.toString(number(processor, "meeFilter"))); evidence.put("expectedAccepted", Boolean.toString(f.accept(phase)));
		evidence.put("accepted", Boolean.toString(accepted)); evidence.put("trimmedBases", Long.toString(longNumber(processor, "basesTrimmedByOverlap")));
	}

	/** Static API receives fresh raw mates and restores their original orientation. */
	private static void probeStatic(final Family f, final String phase, final Read alignedA, final Read alignedB, final Snapshot snapshot, final Map<String, String> evidence){
		final int[] rv=new int[5];
		final int candidate=BBMergeOverlapper.mateByOverlapRatio(alignedA, alignedB, null, null, rv, 20, 20, 1, 1, .09f, .1f, 5.5f, .55f, .95f, .95f, false);
		candidate(f, candidate, rv);
		final Read a=f.spec.first(), b=f.spec.second(); final Snapshot raw=new Snapshot(a, b);
		final int permissive=staticCall(a, b, 0);
		require(permissive==f.spec.length, "REACHABILITY_FAIL: static API cannot pass adapters/efilter without probability filter.");
		require(raw.unchanged(), "Static permissive call changed raw bases/qualities/orientation.");
		final Read filteredA=f.spec.first(), filteredB=f.spec.second();
		final Snapshot filtered=new Snapshot(filteredA, filteredB);
		final int returned=staticCall(filteredA, filteredB, f.cutoff());
		require(returned==(f.accept(phase) ? f.spec.length : -1), "Static probability cutoff differs from phase oracle.");
		require(raw.unchanged() && filtered.unchanged() && snapshot.unchanged(), "Static/helper calls changed input arrays.");
		evidence.put("candidateInsert", Integer.toString(candidate)); evidence.put("candidateBad", Integer.toString(rv[2])); evidence.put("candidateAmbiguous", Integer.toString(rv[4]));
		evidence.put("candidateSource", "standalone ratio finder with minOverlap20, not the static API internal vector");
		evidence.put("staticMinOverlapPolicy", "entropy=0 replaces supplied minOverlap with BBMerge.MIN_OVERLAPPING_BASES");
		evidence.put("staticPermissiveReturned", Integer.toString(permissive)); evidence.put("returned", Integer.toString(returned));
		evidence.put("pfilter", Float.toString(f.cutoff())); evidence.put("accepted", Boolean.toString(returned>0));
	}
	private static int staticCall(final Read a, final Read b, final float cutoff){return BBMerge.findOverlap(a, b, false, 20, 20, 1, 1, 0, .09f, .1f, 5.5f, .55f, 6, .05f, cutoff);}

	/** One real NEW-worker call, optionally capturing all23 actual network inputs. */
	private static void probeWorker(final Path root, final Path work, final Family f, final String phase, final Read a, final Read b, final Snapshot snapshot, final Map<String, String> evidence) throws Exception{
		final BBMerge merge=new BBMerge(arguments(root, work, f).toArray(new String[0]));
		require(!bool(merge, "overlapUsingQuality") && bool(merge, "overlapWithoutQuality") && bool(merge, "useRatioMode") && !bool(merge, "useFlatMode") &&
			!bool(merge, "useMEEfilter") && !bool(merge, "qtrim1") && !bool(merge, "qtrim2") && bool(merge, "useQuality") && bool(merge, "useEfilter"), "Unexpected BBMerge candidate/prefilter settings.");
		require(number(merge, "efilterRatio")==6 && number(merge, "pfilterRatio")==f.cutoff() && bool(merge, "MAKE_VECTOR")==f.neural() &&
			!bool(merge, "TAG_CUSTOM") && !bool(merge, "ONLY_OUTPUT_FALSE_POSITIVE") && integer(merge, "MAX_MISMATCHES_R")>=12,
			"Actual filter/vector/mismatch settings differ from the case.");
		final Class<?> type=Class.forName("jgi.BBMerge$MateThread");
		final Constructor<?> ctor=type.getDeclaredConstructor(BBMerge.class, Streamer.class, Writer.class, Writer.class, fileIO.ByteStreamWriter.class, boolean.class, boolean.class);
		ctor.setAccessible(true); final Object worker=ctor.newInstance(merge, null, null, null, null, true, false);
		require(((Thread)worker).getState()==Thread.State.NEW, "Probe must not start a worker or production stream.");
		int candidate=-1;
		if(!f.neural()){
			final int min=(Integer)invoke(worker, "calcMinOverlapFromEntropy", new Class<?>[]{Read.class, Read.class}, a, b);
			candidate=(Integer)invoke(worker, "mateByOverlap_ratioMode", new Class<?>[]{Read.class, Read.class, int.class}, a, b, min);
			candidate(f, candidate, (int[])value(worker, "rvector"));
		}
		final int returned=(Integer)invoke(worker, "mateByOverlap", new Class<?>[]{Read.class, Read.class}, a, b);
		final int[] rv=(int[])value(worker, "rvector");
		if(f.neural()){
			final FloatList vector=b.fetchVector(); require(vector!=null && vector.size==23, "REACHABILITY_FAIL: missing complete actual network packet.");
			candidate=Math.round(vector.get(5)/.004f);
			require(Float.floatToRawIntBits(vector.get(5))==Float.floatToRawIntBits(f.spec.length*.004f), "Network candidate insert differs from known fragment.");
			final float probability=BBMergeOverlapper.probability(a, b, candidate), expected=BBMergeOverlapper.expectedMismatches(a, b, candidate);
			final float feature=f.disabled() ? .1f : (float)(Math.sqrt(probability)+.0000015);
			require(bits(vector.get(22)).equals(bits(feature)) && bits(vector.get(21)).equals(bits(expected)), "Actual helper feature units changed; preserve the existing extra sqrt and disabled fallback.");
			final StringBuilder values=new StringBuilder("input\tfloat\trawBitsHex\n");
			for(int i=0; i<23; i++){require(Float.isFinite(vector.get(i)), "Nonfinite vector input"+i); values.append(i).append('\t').append(vector.get(i)).append('\t').append(bits(vector.get(i))).append('\n');}
			write(work.resolve("vector.tsv"), values.toString());
			final CellNet net=(CellNet)value(worker, "net"); require(net!=null, "Real worker network missing.");
			final float score=net.getOutput(0), cutoff=number(merge, "netCutoff");
			require(Float.isFinite(score) && cutoff==.872857f && returned==(score<cutoff ? BBMerge.RET_AMBIG : f.spec.length), "Cached default-cutoff decision differs from first-call return.");
			evidence.put("score", Float.toString(score)); evidence.put("scoreBits", bits(score)); evidence.put("cutoff", Float.toString(cutoff));
			evidence.put("packetScope", "one mateByOverlap call; separate CLI retries not instrumented; accuracy not tested");
		}else{require(returned==(f.accept(phase) ? f.spec.length : BBMerge.RET_NO_SOLUTION), "Classical worker probability cutoff differs from phase oracle.");}
		candidate(f, candidate, rv);
		final float expected=BBMergeOverlapper.expectedMismatches(a, b, candidate);
		require((expected+.05f)*6>=f.spec.mismatches && snapshot.unchanged(), "REACHABILITY_FAIL: efilter intervened or worker mutated input arrays.");
		evidence.put("candidateInsert", Integer.toString(candidate)); evidence.put("candidateBad", Integer.toString(rv[2])); evidence.put("candidateAmbiguous", Integer.toString(rv[4]));
		evidence.put("helperExpected", Float.toString(expected)); evidence.put("pfilter", Float.toString(f.cutoff())); evidence.put("returned", Integer.toString(returned)); evidence.put("accepted", Boolean.toString(returned>0));
	}

	private static void candidate(final Family f, final int insert, final int[] rv){require(insert==f.spec.length && rv[2]==f.spec.mismatches && rv[4]==0, "REACHABILITY_FAIL: expected known insert/mismatch count/nonambiguity in "+f.name+", got "+insert+"/"+Arrays.toString(rv));}

	/** Shared explicit constructor/CLI settings; ordered paired output is required. */
	private static ArrayList<String> arguments(final Path root, final Path work, final Family f){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in1="+work.resolve("A.fq"), "in2="+work.resolve("B.fq"), "ordered=t",
			"t=1", "simd=f", "qin=33", "qout=33", "ow=t", "qtrim=f", "maq=0", "minlength=1"));
		if(f.duk()){
			args.addAll(Arrays.asList("out1="+work.resolve("trimmed1.fq"), "out2="+work.resolve("trimmed2.fq"), "tbo=t", "usequality=t", "strictoverlap=f", "minoverlap=20", "mininsert=1", "tpe=f"));
		}else{
			args.addAll(Arrays.asList("out="+work.resolve("merged.fq"), "outu1="+work.resolve("unmerged1.fq"), "outu2="+work.resolve("unmerged2.fq"),
				"usequality=t", "ouq=f", "owq=t", "ratiomode=t", "flatmode=f", "entropy=f", "minoverlap=20", "minoverlap0=20", "ratiominoverlapreduction=0",
				"mininsert=1", "mininsert0=1", "efilter=6", "efilteroffset=.05", "pfilter="+f.cutoff(), "mee=0", "tof=0", "extend=0", "extend2=0", "kfilter=0", "nn="+f.neural()));
			if(f.neural()){args.add("net="+root.resolve("networks/bbmerge.bbnet.gz")); args.add("cutoff=.872857");}
		}
		return args;
	}

	/** Owns only its returned Process; root monitors the finite child JVMs. */
	private int child(final ArrayList<String> command, final Path work, final String name) throws Exception{
		final StringBuilder saved=new StringBuilder();
		for(int i=0; i<command.size(); i++){saved.append(i).append('\t').append(command.get(i)).append('\n');}
		write(work.resolve(name+".command.tsv"), saved.toString());
		final ProcessBuilder builder=new ProcessBuilder(command).directory(root.toFile());
		builder.redirectOutput(work.resolve(name+".stdout").toFile()); builder.redirectError(work.resolve(name+".stderr").toFile());
		builder.environment().remove("_JAVA_OPTIONS"); builder.environment().remove("JDK_JAVA_OPTIONS");
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final Process process=builder.start(); children++;
		final int code=process.waitFor(); write(work.resolve(name+".exit"), code+"\n"); return code;
	}

	private static List<Read> reads(final Path path){
		require(Files.isRegularFile(path), "Missing output file: "+path);
		FASTQ.ASCII_OFFSET=33; FASTQ.DETECT_QUALITY=false; FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		return StreamerFactory.getReads(-1, false, FileFormat.testInput(path.toString(), FileFormat.FASTQ, null, true, false), null, null, null);
	}
	private static boolean one(final List<Read> actual, final Read expected, final int length){return actual.size()==1 && actual.get(0).mate==null && rawPrefix(actual.get(0), expected, length);}
	private static boolean rawPrefix(final Read actual, final Read expected, final int length){return actual.id.equals(expected.id) && Arrays.equals(actual.bases, Arrays.copyOf(expected.bases, length)) && Arrays.equals(actual.quality, Arrays.copyOf(expected.quality, length));}
	private static String records(final String channel, final List<Read> reads){final StringBuilder text=new StringBuilder(); for(final Read r : reads){text.append(channel).append('\t').append(r.id).append('\t').append(new String(r.bases, StandardCharsets.US_ASCII)).append('\t').append(Arrays.toString(r.quality)).append('\n');} return text.toString();}
	private static void writeFastq(final Path path, final Read r) throws IOException{
		final StringBuilder text=new StringBuilder("@").append(r.id).append('\n').append(new String(r.bases, StandardCharsets.US_ASCII)).append("\n+\n");
		for(final byte q : r.quality){text.append((char)(q+33));} write(path, text.append('\n').toString());
	}
	private static void verifyAdapters(final Path work, final Spec spec) throws IOException{
		final List<String> a=Files.readAllLines(work.resolve("A.fq"), StandardCharsets.US_ASCII), b=Files.readAllLines(work.resolve("B.fq"), StandardCharsets.US_ASCII);
		require(a.size()==4 && b.size()==4 && a.get(1).substring(spec.length).equals(ADAPTER_A60) && b.get(1).substring(spec.length).equals(ADAPTER_B_RAW60), "Literal raw adapter wire format differs.");
		require(new String(reverseComplement(spec.alignedSecond().bases), StandardCharsets.US_ASCII).equals(b.get(1)), "Independent raw/aligned B literals disagree.");
	}
	private static void link(final Read a, final Read b){a.mate=b; b.mate=a; b.setPairnum(1);}
	private static boolean close(final float observed, final double expected){return Float.isFinite(observed) && observed>=0 && observed<=1 && (expected==0 ? Float.floatToRawIntBits(observed)==0 : observed>0 && Math.abs(observed-expected)<=Math.max(8*Math.ulp((float)expected), .000002*expected));}
	private static Object value(final Object owner, final String name) throws Exception{final Field f=owner.getClass().getDeclaredField(name); f.setAccessible(true); return f.get(owner);}
	private static Object invoke(final Object owner, final String name, final Class<?>[] types, final Object... args) throws Exception{final Method m=owner.getClass().getDeclaredMethod(name, types); m.setAccessible(true); return m.invoke(owner, args);}
	private static boolean bool(final Object owner, final String name) throws Exception{return (Boolean)value(owner, name);}
	private static int integer(final Object owner, final String name) throws Exception{return ((Number)value(owner, name)).intValue();}
	private static long longNumber(final Object owner, final String name) throws Exception{return ((Number)value(owner, name)).longValue();}
	private static float number(final Object owner, final String name) throws Exception{return ((Number)value(owner, name)).floatValue();}
	private static String modelPin(final Path root) throws Exception{
		final byte[] digest=MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(root.resolve("networks/bbmerge.bbnet.gz"))); final StringBuilder text=new StringBuilder();
		for(int i=digest.length-10; i<digest.length; i++){text.append(Character.forDigit((digest[i]>>>4)&15, 16)).append(Character.forDigit(digest[i]&15, 16));} return text.toString();
	}
	private void same(final Path a, final Path b, final String why) throws IOException{check(Arrays.equals(Files.readAllBytes(a), Files.readAllBytes(b)), why);}
	private void check(final boolean condition, final String why){checks++; if(!condition){failures++; System.err.println("CHECK_FAIL\t"+why);}}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}
	private static Path path(final String value){require(value!=null && !value.isEmpty(), "Missing required path option."); return Paths.get(value).toAbsolutePath().normalize();}
	private static String bits(final float f){return String.format(java.util.Locale.ROOT, "%08x", Float.floatToRawIntBits(f));}
	private static void write(final Path p, final String text) throws IOException{Files.write(p, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path p, final String text) throws IOException{Files.write(p, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);}
	private static String read(final Path p) throws IOException{return new String(Files.readAllBytes(p), StandardCharsets.US_ASCII);}
	private static Map<String, String> table(final Path p) throws IOException{
		final Map<String, String> map=new LinkedHashMap<String, String>();
		for(final String line : read(p).split("\n")){final String[] pair=Parse.splitOnFirst(line, '\t'); require(pair.length==2 && map.put(pair[0], pair[1])==null, "Malformed/duplicate evidence row.");} return map;
	}

	/** Input recipe is fixed in advance; failure never triggers a seed or sequence search. */
	private static String fragment(final int length){
		int state=0x6d2b79f5; final StringBuilder text=new StringBuilder();
		for(int i=0; i<length; i++){state^=state<<13; state^=state>>>17; state^=state<<5; text.append("ACGT".charAt((state>>>30)&3));} return text.toString();
	}
	private static byte[] reverseComplement(final byte[] bases){
		final byte[] out=new byte[bases.length];
		for(int i=0; i<bases.length; i++){final byte b=bases[i]; require(b=='A' || b=='C' || b=='G' || b=='T', "Independent complement only accepts ACGT."); out[bases.length-1-i]=(byte)(b=='A' ? 'T' : b=='T' ? 'A' : b=='C' ? 'G' : 'C');} return out;
	}
	private static byte[] quality(final int n, final int q){final byte[] out=new byte[n]; Arrays.fill(out, (byte)q); return out;}
	private static Family family(final String name){require(name!=null && Arrays.asList(NAMES).contains(name), "Unknown probability family: "+name); final int dash=name.indexOf('-'); return new Family(name, name.substring(0, dash), new Spec(name.substring(dash+1)));}

	private static final class Family{
		Family(final String name_, final String route_, final Spec spec_){name=name_; route=route_; spec=spec_;}
		boolean duk(){return route.equals("duk");} boolean isStatic(){return route.equals("static");} boolean disabled(){return route.equals("disabled");} boolean neural(){return route.equals("neural") || disabled();}
		float cutoff(){return neural() ? (disabled() ? 0 : .00004f) : spec.name.equals("first") ? .98f : spec.name.equals("block") ? .00005f : .5f;}
		boolean accept(final String phase){
			if(spec.name.equals("ordinary")){return true;}
			if(duk()){return !phase.equals("fixed") ? !spec.name.equals("last") : !spec.name.equals("block");}
			return phase.equals("baseline") ? !spec.name.equals("last") : spec.name.equals("last");
		}
		final String name, route; final Spec spec;
	}
	private static final class Spec{
		Spec(final String name_){
			name=name_; require(Arrays.asList("ordinary", "first", "last", "block").contains(name), "Unknown input recipe.");
			length=name.equals("ordinary") ? 100 : 1024; q=name.equals("ordinary") ? 20 : name.equals("block") ? 10 : 4;
			mismatches=name.equals("ordinary") ? 0 : name.equals("block") ? 12 : 1;
			fragment=name.equals("ordinary") ? F100 : fragment(length);
			require(fragment.length()==length, "Predetermined fragment must match the insert/quality geometry: "+name);
			final byte[] changed=fragment.getBytes(StandardCharsets.US_ASCII);
			for(int i=0; i<mismatches; i++){final int at=name.equals("last") ? length-1 : i; changed[at]=(byte)"CGTA".charAt("ACGT".indexOf(changed[at]));}
			mutated=new String(changed, StandardCharsets.US_ASCII);
			require(ADAPTER_A60.length()==60 && ADAPTER_B_RAW60.length()==60 && ADAPTER_B_RC60.length()==60, "Keep recognized adapter geometry60.");
		}
		Read first(){return new Read((fragment+ADAPTER_A60).getBytes(StandardCharsets.US_ASCII), qualities(false), name+"/1", 0, false);}
		Read second(){return new Read((new String(reverseComplement(mutated.getBytes(StandardCharsets.US_ASCII)), StandardCharsets.US_ASCII)+ADAPTER_B_RAW60).getBytes(StandardCharsets.US_ASCII), qualities(false), name+"/2", 1, false);}
		Read alignedSecond(){return new Read((ADAPTER_B_RC60+mutated).getBytes(StandardCharsets.US_ASCII), qualities(true), name+"/2", 1, false);}
		byte[] qualities(final boolean aligned){final byte[] out=quality(length+60, q); if(aligned){Arrays.fill(out, 0, 60, (byte)40);}else{Arrays.fill(out, length, length+60, (byte)40);} return out;}
		Read joined(){
			final byte[] bases=fragment.getBytes(StandardCharsets.US_ASCII), quals=quality(length, q+q/4);
			for(int i=0; i<mismatches; i++){final int at=name.equals("last") ? length-1 : i; bases[at]='N'; quals[at]=0;}
			return new Read(bases, quals, name+"/1", 0, false);
		}
		double probability(final String phase){
			if(name.equals("ordinary")){return 1;}if(phase.equals("baseline")){return name.equals("last") ? 0 : 1;}
			final float correct=q==4 ? .6019f : .9000f, pc=correct*correct, match=pc+(1-pc)*.25f, error=1-match;
			return Math.pow((double)error/match, mismatches*.5);
		}
		final String name, fragment, mutated; final int length, q, mismatches;
	}
	private static final class Snapshot{
		Snapshot(final Read a_, final Read b_){a=a_; b=b_; ab=a.bases; bb=b.bases; aq=a.quality; bq=b.quality; ac=ab.clone(); bc=bb.clone(); qc=aq.clone(); rc=bq.clone();}
		boolean unchanged(){return a.bases==ab && b.bases==bb && a.quality==aq && b.quality==bq && Arrays.equals(ab, ac) && Arrays.equals(bb, bc) && Arrays.equals(aq, qc) && Arrays.equals(bq, rc);}
		final Read a, b; final byte[] ab, bb, aq, bq, ac, bc, qc, rc;
	}

	private final Path root, out, reference; private final String phase;
	private int families, probes, cli, children, checks, failures;
	private static final String MODEL_PIN="24998a8aa257552d4a22", RECIPE="probability-callers-v1;xorshift32=6d2b79f5;adapters60;equal-qualities";
	private static final String[] NAMES={"duk-ordinary", "duk-first", "duk-last", "duk-block", "worker-ordinary", "worker-first", "worker-last", "worker-block", "static-ordinary", "static-first", "static-last", "static-block", "neural-ordinary", "neural-last", "disabled-last"};
	private static final String ADAPTER_A60="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTG";
	private static final String ADAPTER_B_RAW60="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGC";
	private static final String ADAPTER_B_RC60="GCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC";
	private static final String F100="ACGTTGCATGTCGACATGCTAGCTTCGATCGATGACCTGATCCGATAGCAGTACCGTATCATGCTAGGTCCGATCTACGATACGGTTCAGCAGTACGATG";
}
