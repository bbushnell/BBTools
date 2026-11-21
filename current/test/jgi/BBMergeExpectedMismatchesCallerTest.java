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
import jgi.BBDuk;
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

/** Ten literal caller families for the unequal-length expected-error repair (JGI-083).
 * Four current-model packets, three static efilter boundaries and three strict BBDuk
 * trimming boundaries share independent fragment/quality oracles. Every CLI and probe
 * uses a fresh JVM. A probe starts no BBMerge/BBDuk stream or worker and adds no hook.
 * Fixed phase requires saved baseline evidence; model bytes, cutoff, other22 input
 * columns and unchanged controls must agree. Changed NN decisions are reported for
 * review, without assuming monotonic scores or general model accuracy.
 * The separate equal-length witness preserves the historical ignored filter in
 * baseline; fixed requires the later JGI086 worker-filter repair to reject it.
 * @author Sucrose
 */
public final class BBMergeExpectedMismatchesCallerTest{

	/** mode=suite|probe; phase=baseline|fixed; root=BBTools; outdir=new; fixed reference=baseline. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("mode", "phase", "root", "outdir", "reference", "profile", "geometry").contains(pair[0]),
				"Expected a documented fixture flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate fixture option: "+pair[0]);
		}
		final String mode=options.containsKey("mode") ? options.get("mode") : "suite";
		require(mode.equals("suite") || mode.equals("probe"), "Set mode=suite|probe.");
		final String phase=options.get("phase");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Use the dev launcher with assertions enabled.");
		final Path root=path(options.get("root")), out=path(options.get("outdir"));
		require(Files.isRegularFile(root.resolve("bbmerge.sh")) && Files.isRegularFile(root.resolve("bbduk.sh")), "Set root to the single BBTools runtime.");
		if(mode.equals("probe")){
			final Spec spec=new Spec(options.get("profile"), Integer.parseInt(options.get("geometry")));
			probe(root, out, spec, phase);
			return;
		}
		require(!options.containsKey("profile") && !options.containsKey("geometry"), "Profile/geometry are internal probe arguments.");
		final Path reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require(phase.equals("baseline") ? reference==null : reference!=null, "Only fixed phase requires reference=<successful baseline outdir>.");
		if(reference!=null){
			require(read(reference.resolve("summary.txt")).startsWith("BBMERGE_CALLER_PASS phase=baseline families=10 "), "Reference is not a successful ten-family baseline.");
		}
		final BBMergeExpectedMismatchesCallerTest test=new BBMergeExpectedMismatchesCallerTest(root, out, reference, phase);
		for(final String profile : new String[]{"N", "E", "M"}){
			for(int geometry=1; geometry<=(profile.equals("N") ? 4 : 3); geometry++){
				test.run(new Spec(profile, geometry));
			}
		}
		test.check(test.families==10, "Expected exactly ten case/profile families.");
		test.check(test.ignoredFilterWitnesses==1, "Expected one separate equal-length JGI086 witness.");
		test.check(modelPin(root).equals(MODEL_PIN), "Model bytes changed during this phase.");
		test.check(test.modelRegressions==0, "A known true fragment accepted at baseline was rejected after the helper repair; policy review is required.");
		final String summary="BBMERGE_CALLER_"+(test.failures==0 ? "PASS" : "FAIL")+" phase="+phase+
			" families="+test.families+" ignoredFilterWitnesses="+test.ignoredFilterWitnesses+" modelRegressions="+test.modelRegressions+" checks="+test.checks+" failures="+test.failures+"\n";
		write(out.resolve("summary.txt"), summary);
		System.out.print(summary);
		if(test.failures>0){throw new AssertionError(summary+"Inspect retained CLI/probe evidence; never fit the oracle to a failed candidate.");}
	}

	/** Creates only new evidence; source/model files are read-only dependencies. */
	private BBMergeExpectedMismatchesCallerTest(final Path root_, final Path out_, final Path reference_, final String phase_) throws Exception{
		root=root_; out=out_; reference=reference_; phase=phase_;
		Files.createDirectory(out);
		require(modelPin(root).equals(MODEL_PIN), "Unexpected model bytes; do not substitute or retrain the model.");
		write(out.resolve("model.sha80"), modelPin(root)+"\n");
		write(out.resolve("cases.tsv"), "family\tprobeExit\tcliExit\tmerged\tunmerged1\tunmerged2\toutputMatchesOracle\n");
		write(out.resolve("decisions.tsv"), "family\tbaselineMerged\tfixedMerged\tdecisionChanged\toutputMatchesOracle\tknownTruth\tbaselineOutcome\tfixedOutcome\tregression\n");
		if(reference!=null){same(out.resolve("model.sha80"), reference.resolve("model.sha80"), "Model bytes changed.");}
	}

	/** Preserves all ten outputs, including a failed reachability prerequisite in another family. */
	private void run(final Spec spec) throws Exception{
		families++;
		final Path work=Files.createDirectory(out.resolve(spec.name()));
		writeFastq(work.resolve("A.fq"), spec.first());
		writeFastq(work.resolve("B.fq"), spec.second());
		verifyWrittenAdapters(work, spec);
		final byte[] inputA=Files.readAllBytes(work.resolve("A.fq")), inputB=Files.readAllBytes(work.resolve("B.fq"));
		if(reference!=null){
			same(work.resolve("A.fq"), reference.resolve(spec.name()).resolve("A.fq"), "A input changed.");
			same(work.resolve("B.fq"), reference.resolve(spec.name()).resolve("B.fq"), "Raw B input changed.");
		}
		final ArrayList<String> probe=new ArrayList<String>(Arrays.asList("bash", root.resolve("dev/bbmergeexpectedmismatchescallertest.sh").toString(),
			"mode=probe", "phase="+phase, "root="+root, "outdir="+work, "profile="+spec.profile, "geometry="+spec.geometry));
		final int probeExit=child(probe, work, "probe");
		if(spec.profile.equals("E")){
			check(probeExit==0, spec.name()+": static caller/candidate prerequisite failed.");
			check(Arrays.equals(inputA, Files.readAllBytes(work.resolve("A.fq"))) && Arrays.equals(inputB, Files.readAllBytes(work.resolve("B.fq"))), spec.name()+": static probe changed input files.");
			append(out.resolve("cases.tsv"), spec.name()+"\t"+probeExit+"\tNA\t0\t0\t0\t"+(probeExit==0)+"\n");
			if(spec.geometry==3){ignoredFilterWitness(work);}
			return;
		}
		final ArrayList<String> cli=new ArrayList<String>();
		cli.add("bash"); cli.add(root.resolve(spec.profile.equals("M") ? "bbduk.sh" : "bbmerge.sh").toString());
		cli.add("-Xmx128m"); cli.add("-Xms32m"); cli.addAll(arguments(root, work, spec));
		final int cliExit=child(cli, work, "cli");
		check(probeExit==0, spec.name()+": probe prerequisite failed; inspect probe.stderr.");
		check(cliExit==0, spec.name()+": real CLI failed; inspect cli.stderr.");
		check(Arrays.equals(inputA, Files.readAllBytes(work.resolve("A.fq"))) && Arrays.equals(inputB, Files.readAllBytes(work.resolve("B.fq"))), spec.name()+": CLI/probe changed input files.");
		if(probeExit!=0 || cliExit!=0){return;}
		final Map<String, String> packet=table(work.resolve("probe.tsv"));
		final List<Read> merged=reads(work.resolve(spec.profile.equals("M") ? "trimmed1.fq" : "merged.fq"));
		final List<Read> first=spec.profile.equals("M") ? merged : reads(work.resolve("unmerged1.fq"));
		final List<Read> second=reads(work.resolve(spec.profile.equals("M") ? "trimmed2.fq" : "unmerged2.fq"));
		boolean truth=true;
		if(spec.profile.equals("M")){
			final boolean trim=spec.geometry==3 || phase.equals("baseline");
			truth=one(first, spec.first(), trim ? 100 : spec.a.length()) && one(second, spec.second(), trim ? 100 : spec.b.length());
			check(truth, spec.name()+": strict-MEE trimming, IDs, bases or qualities differ from literal oracle.");
		}else{
			check(merged.size()<=1 && first.size()<=1 && second.size()<=1, spec.name()+": duplicated output records.");
			check(merged.size()+first.size()==1 && first.size()==second.size(), spec.name()+": missing/duplicated pair routing.");
			if(!merged.isEmpty()){
				truth=one(merged, spec.joined(), spec.insert);
				check(truth, spec.name()+": merged ID/sequence/quality differs from known fragment.");
			}else{
				truth=one(first, spec.first(), spec.a.length()) && one(second, spec.second(), spec.b.length());
				check(truth, spec.name()+": unmerged raw mates changed.");
			}
		}
		append(out.resolve("cases.tsv"), spec.name()+"\t"+probeExit+"\t"+cliExit+"\t"+(spec.profile.equals("M") ? 0 : merged.size())+
			"\t"+first.size()+"\t"+second.size()+"\t"+truth+"\n");
		write(work.resolve("outputs.tsv"), records("primary", merged)+records("mate1", spec.profile.equals("M") ? new ArrayList<Read>() : first)+records("mate2", second));
		if(spec.profile.equals("N")){
			write(work.resolve("routing.tsv"), "firstCandidateReturned\tcliMerged\tfirstCallDecisionAgrees\tcliRetryCount\tknownTruth\tcliOutcome\n"+
				packet.get("returned")+"\t"+merged.size()+"\t"+((Integer.parseInt(packet.get("returned"))>0)==(merged.size()==1))+
				"\tnot-instrumented\ttrue-overlap\t"+(merged.size()==1 ? "true-positive" : "false-negative")+"\n");
		}
		if(reference!=null && spec.profile.equals("N")){comparePacket(spec, work, merged.size(), truth);}
	}

	/** Equal lengths isolate JGI086 from the geometry repair: baseline merges, fixed rejects.
	 * The historical witness remains counted within E-G3 in both phases, not as an eleventh family. */
	private void ignoredFilterWitness(final Path control) throws Exception{
		ignoredFilterWitnesses++;
		final Spec spec=new Spec("W", 3);
		final Path work=Files.createDirectory(control.resolve("ignored-filter-witness"));
		writeFastq(work.resolve("A.fq"), spec.first()); writeFastq(work.resolve("B.fq"), spec.second());
		verifyWrittenAdapters(work, spec);
		final byte[] inputA=Files.readAllBytes(work.resolve("A.fq")), inputB=Files.readAllBytes(work.resolve("B.fq"));
		final int probeExit=child(new ArrayList<String>(Arrays.asList("bash", root.resolve("dev/bbmergeexpectedmismatchescallertest.sh").toString(),
			"mode=probe", "phase="+phase, "root="+root, "outdir="+work, "profile=W", "geometry=3")), work, "probe");
		final ArrayList<String> command=new ArrayList<String>(Arrays.asList("bash", root.resolve("bbmerge.sh").toString(), "-Xmx128m", "-Xms32m"));
		command.addAll(arguments(root, work, spec));
		final int cliExit=child(command, work, "cli");
		check(probeExit==0 && cliExit==0, "JGI086 equal-length witness did not complete.");
		check(Arrays.equals(inputA, Files.readAllBytes(work.resolve("A.fq"))) && Arrays.equals(inputB, Files.readAllBytes(work.resolve("B.fq"))), "JGI086 witness changed input files.");
		if(probeExit!=0 || cliExit!=0){return;}
		final List<Read> merged=reads(work.resolve("merged.fq"));
		final List<Read> first=reads(work.resolve("unmerged1.fq")), second=reads(work.resolve("unmerged2.fq"));
		final boolean baseline=phase.equals("baseline");
		check(baseline ? one(merged, spec.joined(), 100) && first.isEmpty() && second.isEmpty()
				: merged.isEmpty() && one(first, spec.first(), spec.a.length()) && one(second, spec.second(), spec.b.length()),
				"JGI086 equal-length routing/IDs/bases/qualities differ from the "+phase+" oracle.");
		write(work.resolve("outputs.tsv"), records("merged", merged)+records("unmerged1", first)+records("unmerged2", second));
		final String policyHeader="sourceContract\tintendedPolicy\tgeometryChanged\n";
		write(work.resolve("policy.tsv"), policyHeader+(baseline ? "merge" : "reject")+"\treject\tfalse\n");
		if(reference!=null){
			final Path before=reference.resolve("E-G3/ignored-filter-witness");
			for(final String file : new String[]{"A.fq", "B.fq"}){same(work.resolve(file), before.resolve(file), "JGI086 equal-length witness input changed: "+file);}
			check(read(before.resolve("policy.tsv")).equals(policyHeader+"merge\treject\tfalse\n"),
					"JGI086 reference must retain the historical ignored-filter policy witness.");
			final Map<String, String> oldProbe=table(before.resolve("probe.tsv"));
			check("true".equals(oldProbe.get("efilterWouldReject")) && "false".equals(oldProbe.get("intendedCliMerge"))
					&& "true".equals(oldProbe.get("ignoredFilterWitness")) && Integer.toString(spec.insert).equals(oldProbe.get("returned")),
					"JGI086 baseline probe must record intended rejection and its actual positive return.");
			write(work.resolve("policy-comparison.tsv"), "phase\t"+policyHeader+
					"baseline\tmerge\treject\tfalse\nfixed\treject\treject\tfalse\n");
		}
	}

	/** Compares raw float bits and complete control FASTQs, without assuming a score direction. */
	private void comparePacket(final Spec spec, final Path work, final int merged, final boolean truth) throws Exception{
		final Path before=reference.resolve(spec.name());
		final String[] a=lines(before.resolve("vector.tsv")), b=lines(work.resolve("vector.tsv"));
		check(a.length==24 && b.length==24, spec.name()+": expected header plus23 inputs.");
		if(a.length==24 && b.length==24){
			for(int i=0; i<23; i++){
				if(i!=21 || spec.geometry>=3){check(a[i+1].equals(b[i+1]), spec.name()+": unexpected vector change at input"+i);}
			}
		}
		final Map<String, String> old=table(before.resolve("probe.tsv")), now=table(work.resolve("probe.tsv"));
		check(old.get("cutoff").equals(now.get("cutoff")), spec.name()+": model cutoff changed.");
		if(spec.geometry>=3){
			check(old.get("scoreBits").equals(now.get("scoreBits")), spec.name()+": unchanged control score changed.");
			for(final String file : new String[]{"merged.fq", "unmerged1.fq", "unmerged2.fq"}){same(work.resolve(file), before.resolve(file), spec.name()+": control output changed: "+file);}
		}
		final int oldMerged=reads(before.resolve("merged.fq")).size();
		final boolean regression=oldMerged==1 && merged==0;
		if(regression){modelRegressions++;}
		append(out.resolve("decisions.tsv"), spec.name()+"\t"+oldMerged+"\t"+merged+"\t"+(oldMerged!=merged)+"\t"+truth+
			"\ttrue-overlap\t"+(oldMerged==1 ? "true-positive" : "false-negative")+"\t"+(merged==1 ? "true-positive" : "false-negative")+"\t"+regression+"\n");
	}

	/** Reads actual constructor settings and one real worker candidate, not a surrogate NN. */
	private static void probe(final Path root, final Path work, final Spec spec, final String phase) throws Exception{
		require(modelPin(root).equals(MODEL_PIN), "Probe model fingerprint changed.");
		final Map<String, String> evidence=new LinkedHashMap<String, String>();
		evidence.put("family", spec.name()); evidence.put("phase", phase);
		final String[] args=arguments(root, work, spec).toArray(new String[0]);
		final Read a=spec.first(), rawB=spec.second();
		a.mate=rawB; rawB.mate=a; rawB.setPairnum(1);
		final float pre=Math.max(a.expectedErrors(false, -1), rawB.expectedErrors(false, -1));
		evidence.put("preMaximumFullReadErrors", Float.toString(pre));
		rawB.reverseComplementFast();
		final Read b=rawB;
		int candidate, bad, ambiguous;
		float expected, probability;
		if(spec.profile.equals("M")){
			final BBDuk duk=new BBDuk(args);
			require(bool(duk, "trimByOverlap") && bool(duk, "useQualityForOverlap") && bool(duk, "strictOverlap"), "BBDuk strict quality-TBO configuration was not reached.");
			require(number(duk, "meeFilter")==15 && number(duk, "efilterRatio")==3.5f && number(duk, "pfilterRatio")==.001f, "Actual strict BBDuk gates differ from reviewed constants.");
			final int[] rv=new int[5];
			candidate=BBMergeOverlapper.mateByOverlapRatio(a, b, new float[a.length()], new float[b.length()], rv,
				integer(duk, "minOverlap0"), integer(duk, "minOverlap"), integer(duk, "minInsert0"), integer(duk, "minInsert"),
				number(duk, "maxRatio"), .12f, number(duk, "ratioMargin"), number(duk, "ratioOffset"), .95f, .95f, true);
			bad=rv[2]; ambiguous=rv[4];
			require(pre<number(duk, "meeFilter"), "REACHABILITY_FAIL: pre-overlap per-read maximum blocks TBO.");
			require(Math.abs(pre-(spec.geometry==3 ? 1.006 : 10.006))<.0005, "Full-read expected errors differ from literal qualities.");
			evidence.put("preGatePassed", "true");
		}else{
			final BBMerge merge=new BBMerge(args);
			require(!bool(merge, "overlapUsingQuality") && bool(merge, "overlapWithoutQuality") && bool(merge, "useRatioMode") && !bool(merge, "useFlatMode"), "Expected unweighted ratio candidate path.");
			require(!bool(merge, "useMEEfilter") && !bool(merge, "qtrim1") && !bool(merge, "qtrim2"), "Unexpected prefilter or trim retry; do not bypass it in the probe.");
			final Class<?> type=Class.forName("jgi.BBMerge$MateThread");
			final Constructor<?> constructor=type.getDeclaredConstructor(BBMerge.class, Streamer.class, Writer.class, Writer.class, fileIO.ByteStreamWriter.class, boolean.class, boolean.class);
			constructor.setAccessible(true);
			final Object worker=constructor.newInstance(merge, null, null, null, null, true, false);
			require(((Thread)worker).getState()==Thread.State.NEW, "Probe must not start a worker or any production stream.");
			if(spec.profile.equals("N")){
				final int returned=(Integer)invoke(worker, "mateByOverlap", new Class<?>[]{Read.class, Read.class}, a, b);
				final FloatList vector=b.fetchVector();
				require(vector!=null && vector.size==23, "REACHABILITY_FAIL: actual first candidate lacks23vector inputs.");
				final StringBuilder values=new StringBuilder("input\tfloat\trawBitsHex\n");
				for(int i=0; i<23; i++){
					require(Float.isFinite(vector.get(i)), "Non-finite actual vector input"+i);
					values.append(i).append('\t').append(vector.get(i)).append('\t').append(Integer.toHexString(Float.floatToRawIntBits(vector.get(i)))).append('\n');
				}
				write(work.resolve("vector.tsv"), values.toString());
				// Overlapper's MAKE_VECTOR block: column5 is INSERT*.004, not overlap or raw bases.
				candidate=Math.round(vector.get(5)/.004f);
				require(Float.floatToRawIntBits(vector.get(5))==Float.floatToRawIntBits(spec.insert*.004f), "REACHABILITY_FAIL: vector candidate is not known insert"+spec.insert);
				final int[] rv=(int[])value(worker, "rvector"); bad=rv[2]; ambiguous=rv[4];
				final CellNet net=(CellNet)value(worker, "net");
				require(net!=null, "Actual candidate was not scored by the current network.");
				// Reads cached activation only; never calls feedForward a second time.
				final float score=net.getOutput(0), cutoff=number(merge, "netCutoff");
				require(Float.isFinite(score) && cutoff==.872857f, "Invalid cached score or changed cutoff.");
				evidence.put("score", Float.toString(score)); evidence.put("scoreBits", Integer.toHexString(Float.floatToRawIntBits(score)));
				evidence.put("cutoff", Float.toString(cutoff)); evidence.put("returned", Integer.toString(returned));
				evidence.put("returnCategory", returned>0 ? "INSERT" : returned==BBMerge.RET_AMBIG ? "AMBIG" : "OTHER");
				evidence.put("packetScope", "one-mateByOverlap-call; separate-CLI-retries-not-instrumented");
				require(returned==(score<cutoff ? BBMerge.RET_AMBIG : spec.insert), "First-call return does not match cached network decision.");
				expected=vector.get(21); probability=BBMergeOverlapper.probability(a, b, spec.insert);
				require(Float.floatToRawIntBits(vector.get(22))==Float.floatToRawIntBits((float)(Math.sqrt(probability)+.0000015)), "Vector22 must retain sqrt(probability)+.0000015 units.");
			}else{
				final int min=(Integer)invoke(worker, "calcMinOverlapFromEntropy", new Class<?>[]{Read.class, Read.class}, a, b);
				candidate=(Integer)invoke(worker, "mateByOverlap_ratioMode", new Class<?>[]{Read.class, Read.class, int.class}, a, b, min);
				final int[] rv=(int[])value(worker, "rvector"); bad=rv[2]; ambiguous=rv[4];
				require(number(merge, "efilterRatio")==.06f && number(merge, "efilterOffset")==.05f && number(merge, "pfilterRatio")==0, "E profile did not isolate the actual error filter.");
				expected=BBMergeOverlapper.expectedMismatches(a, b, spec.insert);
				probability=BBMergeOverlapper.probability(a, b, spec.insert);
				final boolean intendedReject=(expected+.05f)*.06f<bad;
				evidence.put("efilterWouldReject", Boolean.toString(intendedReject));
				evidence.put("intendedCliMerge", Boolean.toString(!intendedReject));
				if(spec.profile.equals("W")){
					final Read freshA=spec.first(), freshB=spec.second(); freshB.reverseComplementFast();
					final int returned=(Integer)invoke(worker, "mateByOverlap", new Class<?>[]{Read.class, Read.class}, freshA, freshB);
					evidence.put("returned", Integer.toString(returned));
					evidence.put("ignoredFilterWitness", Boolean.toString(intendedReject && returned>0));
					final int expectedReturn=phase.equals("baseline") ? spec.insert : BBMerge.RET_AMBIG;
					require(intendedReject && returned==expectedReturn, "Equal-length worker-filter return differs for "+phase+
							"; expected="+expectedReturn+", actual="+returned);
				}
				// Unlike the worker, public findOverlap1659 propagates ambig as-1.
				// It receives NEW raw-orientation mates and reverses B itself.
				final int permissive=BBMerge.findOverlap(spec.first(), spec.second(), false, 20, 20, 1, 1, 0,
					number(merge, "MAX_RATIO"), number(merge, "MIN_SECOND_RATIO"), number(merge, "RATIO_MARGIN"), number(merge, "RATIO_OFFSET"), 6f, .05f, 0);
				evidence.put("staticPermissiveReturned", Integer.toString(permissive));
				require(permissive==spec.insert, "REACHABILITY_FAIL: static path cannot find the known insert before tightening only efilter.");
				final int staticReturned=BBMerge.findOverlap(spec.first(), spec.second(), false, 20, 20, 1, 1, 0,
					number(merge, "MAX_RATIO"), number(merge, "MIN_SECOND_RATIO"), number(merge, "RATIO_MARGIN"), number(merge, "RATIO_OFFSET"), .06f, .05f, 0);
				evidence.put("staticReturned", Integer.toString(staticReturned));
				evidence.put("staticExpected", Integer.toString(intendedReject ? -1 : spec.insert));
				require(staticReturned==(intendedReject ? -1 : spec.insert), "REACHABILITY_FAIL: live static efilter boundary differs from independent inequality.");
			}
			finishProbe(work, spec, phase, evidence, candidate, bad, ambiguous, expected, probability);
			return;
		}
		expected=BBMergeOverlapper.expectedMismatches(a, b, spec.insert);
		probability=BBMergeOverlapper.probability(a, b, spec.insert);
		require(bad==0 && probability>=.001f && (expected+.05f)*3.5f>=bad, "REACHABILITY_FAIL: another strict-TBO filter hides the MEE boundary.");
		evidence.put("postMeeWouldReject", Boolean.toString(expected>15));
		finishProbe(work, spec, phase, evidence, candidate, bad, ambiguous, expected, probability);
	}

	/** A missing/ambiguous candidate is a failure, never an observed insert supplied by the test. */
	private static void finishProbe(final Path work, final Spec spec, final String phase, final Map<String, String> evidence,
			final int candidate, final int bad, final int ambiguous, final float expected, final float probability) throws Exception{
		evidence.put("candidateInsert", Integer.toString(candidate)); evidence.put("candidateBad", Integer.toString(bad));
		evidence.put("candidateAmbiguous", Integer.toString(ambiguous)); evidence.put("helperExpected", Float.toString(expected));
		evidence.put("probability", Float.toString(probability));
		final StringBuilder text=new StringBuilder();
		for(final Map.Entry<String, String> entry : evidence.entrySet()){text.append(entry.getKey()).append('\t').append(entry.getValue()).append('\n');}
		write(work.resolve("probe.tsv"), text.toString());
		require(candidate==spec.insert && ambiguous==0 && bad==((spec.profile.equals("E") && spec.geometry<3) || spec.profile.equals("W") ? 1 : 0), "REACHABILITY_FAIL: actual insert/bad/ambiguity differs from known candidate in "+spec.name());
		final double oracle=phase.equals("baseline") ? spec.sourceExpected : spec.fixedExpected;
		require(Float.isFinite(expected) && Math.abs(expected-oracle)<.0005, "Expected-error sum differs from independent decimal oracle: "+expected+" versus "+oracle);
	}

	/** Identical real CLI options are also used by the separate constructor probe. */
	private static ArrayList<String> arguments(final Path root, final Path work, final Spec spec){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in1="+work.resolve("A.fq"), "in2="+work.resolve("B.fq"),
			"t=1", "simd=f", "qin=33", "qout=33", "ow=t", "qtrim=f", "maq=0", "minlength=1"));
		if(spec.profile.equals("M")){
			args.addAll(Arrays.asList("out1="+work.resolve("trimmed1.fq"), "out2="+work.resolve("trimmed2.fq"),
				"tbo=t", "usequality=t", "strictoverlap=t", "minoverlap=20", "mininsert=1", "tpe=f"));
		}else{
			args.addAll(Arrays.asList("out="+work.resolve("merged.fq"), "outu1="+work.resolve("unmerged1.fq"), "outu2="+work.resolve("unmerged2.fq")));
			if(spec.profile.equals("N")){args.addAll(Arrays.asList("nn=t", "net="+root.resolve("networks/bbmerge.bbnet.gz"), "cutoff=0.872857"));}
			else{args.addAll(Arrays.asList("nn=f", "usequality=t", "ouq=f", "entropy=f", "minoverlap=20", "minoverlap0=20", "ratiominoverlapreduction=0",
				"mininsert=1", "mininsert0=1", "efilter=0.06", "efilteroffset=0.05", "pfilter=0", "mee=0", "tof=0", "extend=0", "extend2=0", "kfilter=0"));}
		}
		return args;
	}

	/** Starts exactly the requested shell launcher; root monitors these finite tiny inputs.
	 * Waits on the owned Process; never searches for or signals processes by name. */
	private int child(final ArrayList<String> command, final Path work, final String name) throws Exception{
		final StringBuilder saved=new StringBuilder();
		for(int i=0; i<command.size(); i++){saved.append(i).append('\t').append(command.get(i)).append('\n');}
		write(work.resolve(name+".command.tsv"), saved.toString());
		final ProcessBuilder builder=new ProcessBuilder(command).directory(root.toFile());
		builder.redirectOutput(work.resolve(name+".stdout").toFile()); builder.redirectError(work.resolve(name+".stderr").toFile());
		builder.environment().remove("_JAVA_OPTIONS"); builder.environment().remove("JDK_JAVA_OPTIONS");
		builder.environment().put("JAVA_TOOL_OPTIONS", "-XX:ActiveProcessorCount=1 -Dfile.encoding=UTF-8");
		final Process process=builder.start();
		final int code=process.waitFor();
		write(work.resolve(name+".exit"), code+"\n");
		return code;
	}

	/** Uses the native FASTQ reader with explicit33 and no interleaving inference. */
	private static List<Read> reads(final Path path){
		require(Files.isRegularFile(path), "Missing expected output file: "+path);
		FASTQ.ASCII_OFFSET=33; FASTQ.DETECT_QUALITY=false;
		FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=false;
		final FileFormat ff=FileFormat.testInput(path.toString(), FileFormat.FASTQ, null, true, false);
		return StreamerFactory.getReads(-1, false, ff, null, null, null);
	}

	/** Exact record oracle; trimming removes the raw3-prime tail from either mate. */
	private static boolean one(final List<Read> actual, final Read expected, final int length){
		if(actual.size()!=1){return false;}
		final Read r=actual.get(0);
		return r.mate==null && expected.id.equals(r.id) && Arrays.equals(r.bases, Arrays.copyOf(expected.bases, length)) &&
			Arrays.equals(r.quality, Arrays.copyOf(expected.quality, length));
	}

	/** Saves every parsed output ID, sequence and numeric quality in ordinary text. */
	private static String records(final String channel, final List<Read> reads){
		final StringBuilder text=new StringBuilder();
		for(final Read r : reads){text.append(channel).append('\t').append(r.id).append('\t').append(new String(r.bases, StandardCharsets.US_ASCII)).append('\t').append(Arrays.toString(r.quality)).append('\n');}
		return text.toString();
	}

	/** Writes literal ASCII33 FASTQ; no tested formatter participates in the input oracle. */
	private static void writeFastq(final Path path, final Read read) throws IOException{
		final StringBuilder text=new StringBuilder("@").append(read.id).append('\n').append(new String(read.bases, StandardCharsets.US_ASCII)).append("\n+\n");
		for(final byte q : read.quality){text.append((char)(q+33));}
		write(path, text.append('\n').toString());
	}

	/** Checks the literal input wire format, not output parsing. E/W alone need the
	 * static API's default-adapter prerequisite; N/M keep their accepted C/G tails. */
	private static void verifyWrittenAdapters(final Path work, final Spec spec) throws IOException{
		if(!spec.recognizedAdapters){return;}
		if(spec.geometry==2 || spec.geometry==3){
			final List<String> a=Files.readAllLines(work.resolve("A.fq"), StandardCharsets.US_ASCII);
			require(a.size()==4 && a.get(1).endsWith(ADAPTER_A60), "Written E/W mate1 must retain its recognized60-base adapter.");
		}
		if(spec.geometry==1 || spec.geometry==3){
			final List<String> b=Files.readAllLines(work.resolve("B.fq"), StandardCharsets.US_ASCII);
			require(b.size()==4 && b.get(1).endsWith(ADAPTER_B_RAW60), "Written raw mate2 must end in the intended adapter, not its reverse complement.");
		}
	}

	/** Reflection stays in tests; only actual state is read, never overwritten. */
	private static Object value(final Object owner, final String name) throws Exception{
		final Field field=owner.getClass().getDeclaredField(name); field.setAccessible(true); return field.get(owner);
	}
	private static boolean bool(final Object owner, final String name) throws Exception{return (Boolean)value(owner, name);}
	private static int integer(final Object owner, final String name) throws Exception{return ((Number)value(owner, name)).intValue();}
	private static float number(final Object owner, final String name) throws Exception{return ((Number)value(owner, name)).floatValue();}
	private static Object invoke(final Object owner, final String name, final Class<?>[] types, final Object... args) throws Exception{
		final Method method=owner.getClass().getDeclaredMethod(name, types); method.setAccessible(true); return method.invoke(owner, args);
	}

	/** Final80bits only, including the compressed model bytes exactly as supplied. */
	private static String modelPin(final Path root) throws Exception{
		final byte[] digest=MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(root.resolve("networks/bbmerge.bbnet.gz")));
		final StringBuilder pin=new StringBuilder();
		for(int i=digest.length-10; i<digest.length; i++){pin.append(Character.forDigit((digest[i]>>>4)&15, 16)).append(Character.forDigit(digest[i]&15, 16));}
		return pin.toString();
	}

	private void same(final Path a, final Path b, final String message) throws IOException{check(Arrays.equals(Files.readAllBytes(a), Files.readAllBytes(b)), message);}
	private void check(final boolean condition, final String message){checks++; if(!condition){failures++; System.err.println("CHECK_FAIL\t"+message);}}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}
	private static Path path(final String value){require(value!=null && !value.isEmpty(), "Missing required path option."); return Paths.get(value).toAbsolutePath().normalize();}
	private static void write(final Path path, final String value) throws IOException{Files.write(path, value.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path path, final String value) throws IOException{Files.write(path, value.getBytes(StandardCharsets.UTF_8), StandardOpenOption.APPEND);}
	private static String read(final Path path) throws IOException{return new String(Files.readAllBytes(path), StandardCharsets.UTF_8);}
	private static String[] lines(final Path path) throws IOException{return read(path).split("\n");}
	private static Map<String, String> table(final Path path) throws IOException{
		final Map<String, String> map=new LinkedHashMap<String, String>();
		for(final String line : lines(path)){final String[] pair=Parse.splitOnFirst(line, '\t'); require(pair.length==2, "Malformed probe evidence."); require(map.put(pair[0], pair[1])==null, "Duplicate probe evidence key.");}
		return map;
	}

	/** Four fixed geometries. B is aligned/RC; second() restores sequencer orientation. */
	private static final class Spec{
		Spec(final String profile_, final int geometry_){
			profile=profile_; geometry=geometry_;
			require(Arrays.asList("N", "E", "M", "W").contains(profile) && geometry>=1 && geometry<=(profile.equals("N") ? 4 : 3) && (!profile.equals("W") || geometry==3), "Unknown literal family.");
			recognizedAdapters=profile.equals("E") || profile.equals("W");
			if(recognizedAdapters){
				require(ADAPTER_A60.length()==60 && ADAPTER_B_RAW60.length()==60 && ADAPTER_B_RC60.length()==60, "Keep E/W adapter geometry exactly60bases.");
				final Read adapter=new Read(ADAPTER_B_RC60.getBytes(StandardCharsets.US_ASCII), null, "adapter-rc", 0, false);
				adapter.reverseComplementFast();
				require(ADAPTER_B_RAW60.equals(new String(adapter.bases, StandardCharsets.US_ASCII)), "Explicit rawB60 and aligned-prefix RC literals disagree.");
			}
			insert=geometry==4 ? 200 : 100;
			a=geometry==2 || geometry==3 ? F100+(recognizedAdapters ? ADAPTER_A60 : repeat('G', 60)) : F100;
			String aligned=geometry==1 || geometry==3 ? (recognizedAdapters ? ADAPTER_B_RC60 : repeat('C', 60))+F100 : geometry==2 ? F100 : F200.substring(40);
			aq=qualities(a.length(), geometry>=3 ? 20 : 10);
			bq=qualities(aligned.length(), geometry>=3 ? 20 : 10);
			if(geometry==2 || geometry==3){Arrays.fill(aq, 100, 160, (byte)40);}
			if(geometry==1 || geometry==3){Arrays.fill(bq, 0, 60, (byte)40);}
			if(profile.equals("E") && geometry<3){
				aq[50]=20;
				final char[] changed=aligned.toCharArray(); final int index=geometry==1 ? 110 : 50;
				require(changed[index]=='G', "Literal mismatch site must be fragment base50=G."); changed[index]='T'; aligned=new String(changed);
			}
			if(profile.equals("W")){
				aq[50]=30;
				final char[] changed=aligned.toCharArray(); require(changed[110]=='G', "Equal-length witness mismatch must be G."); changed[110]='T'; aligned=new String(changed);
			}
			b=aligned;
			sourceExpected=profile.equals("W") ? 1.98109 : geometry==1 ? (profile.equals("E") ? 13.515409 : 13.6054) : geometry==2 ? 7.6 : geometry==3 ? 1.99 : 1.194;
			fixedExpected=geometry<3 ? (profile.equals("E") ? 18.919 : 19) : sourceExpected;
			require(F200.length()==200 && aq.length==a.length() && bq.length==b.length(), "Literal sequence/quality length mismatch.");
		}
		String name(){return profile+"-G"+geometry;}
		Read first(){return new Read(a.getBytes(StandardCharsets.US_ASCII), aq.clone(), "G"+geometry+"/1", 0, false);}
		Read second(){final Read read=new Read(b.getBytes(StandardCharsets.US_ASCII), bq.clone(), "G"+geometry+"/2", 0, false); read.reverseComplementFast(); return read;}
		Read joined(){
			final byte[] q=qualities(insert, geometry>=3 ? 25 : 12);
			if(geometry==4){Arrays.fill(q, (byte)20); Arrays.fill(q, 40, 100, (byte)25);}
			if((profile.equals("E") && geometry<3) || profile.equals("W")){q[50]=10;}
			return new Read((geometry==4 ? F200 : F100).getBytes(StandardCharsets.US_ASCII), q, "G"+geometry+"/1", 0, false);
		}
		final String profile, a, b;
		final boolean recognizedAdapters;
		final int geometry, insert;
		final byte[] aq, bq;
		final double sourceExpected, fixedExpected;
	}

	private static byte[] qualities(final int length, final int q){final byte[] values=new byte[length]; Arrays.fill(values, (byte)q); return values;}
	private static String repeat(final char base, final int count){final char[] values=new char[count]; Arrays.fill(values, base); return new String(values);}
	private final Path root, out, reference;
	private final String phase;
	private int families, checks, failures, modelRegressions, ignoredFilterWitnesses;
	private static final String MODEL_PIN="24998a8aa257552d4a22";
	/* First60 bases of resources/adapters.fa Reverse_adapter and TruSeq_Adapter_Index_2.
	 * Default static verification compares their exact16-base prefixes. These affect
	 * E/W input reachability only; no quality, helper oracle, N/M input or gate changes. */
	private static final String ADAPTER_A60="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTG";
	private static final String ADAPTER_B_RAW60="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGC";
	private static final String ADAPTER_B_RC60="GCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC";
	private static final String F200="ACGTTGCATGTCGACATGCTAGCTTCGATCGATGACCTGATCCGATAGCA"+
		"GTACCGTATCATGCTAGGTCCGATCTACGATACGGTTCAGCAGTACGATG"+
		"GCATTCAGCTTGCAGATCGACTAGGCTTACGTCATCGAGTAGCATGTCAC"+
		"CTGACTAGCGAATCGGTACAGACTTCAGGTTCGTAGACCTGCTACATGGA";
	private static final String F100=F200.substring(0, 100);
}
