package test.jgi;

import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import jgi.BBMerge;
import parse.Parse;

/** Two real BBMerge main phases in one fresh JVM expose stale false-negative totals.
 * The same truthful insert=40 pair first fails a length filter, then must merge at 40.
 * Snapshots only read production counters; no field is seeded and no reset is inserted
 * between phases. One final public reset supplies controls from actual populated state.
 * This does not model BBMergeMulti's normal fresh-JVM commands or measure performance.
 * @author Sucrose
 */
public final class BBMergeCounterResetTest{

	/** phase=baseline|fixed, root=RUNTIME, outdir=NEW; fixed requires reference=BASELINE. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("phase", "root", "outdir", "reference").contains(pair[0]), "Expected documented flag=value: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate option: "+pair[0]);
		}
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbmergecounterresettest.sh with assertions enabled.");
		final String phase=options.get("phase");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		final Path root=path(options.get("root")), out=path(options.get("outdir"));
		final Path reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require(phase.equals("baseline") ? reference==null : reference!=null, "Only fixed requires a saved baseline reference.");
		require(Files.isRegularFile(root.resolve("bbmerge.sh")), "Set root to the single BBTools runtime.");
		if(reference!=null){
			require(read(reference.resolve("summary.txt")).matches("BBMERGE_COUNTER_RESET_PASS phase=baseline version=1 mainCalls=2 phases=2 directResets=1 discrepancies=3 checks=[0-9]+ failures=0\\n"), "Reference must be a successful two-phase baseline.");
		}
		final BBMergeCounterResetTest test=new BBMergeCounterResetTest(out, reference, phase);
		try{test.exercise();}
		finally{test.saveDispatch();}
		test.finish();
	}

	private BBMergeCounterResetTest(final Path out_, final Path reference_, final String phase_) throws IOException{
		out=out_; reference=reference_; phase=phase_; baseline=phase.equals("baseline");
		Files.createDirectory(out);
		write(out.resolve("inputs.tsv"), "mate\tbases\tinsert\nR1\t"+FRAGMENT+"\t40\nR2\t"+RAW_INPUT+"\t40\n");
		write(out.resolve("oracles.tsv"), "stage\tpairs\tbases\tcorrect\tincorrect\thist40\thistSum\tbaselineFN\tfixedFN\tcorrectFNR\n"+
			"phase-1\t1\t96\t0\t0\t0\t0\t1\t1\t100.0000%\n"+
			"phase-2\t1\t96\t1\t0\t1\t1\t1\t0\t0.0000%\n"+
			"direct-reset\t0\t0\t0\t0\t0\t0\t1\t0\t-\n");
		write(out.resolve("settings.tsv"), "scope\targument\n"+settings());
		write(out.resolve("controls.tsv"), "stage\tfield\tvalue\n");
	}

	/** Every later operation depends on the preceding real phase and literal oracle. */
	private void exercise() throws Exception{
		final Snapshot initial=snapshot("initial");
		for(final String name : COUNTERS){check(initial.value(name)==0, "Fresh JVM required; nonzero initial "+name);}
		check(sum(initial.hist)==0, "Fresh JVM histogram must be empty; do not manufacture it with a reset.");
		check(!field("errorState").getBoolean(null), "Fresh runtime already reports an I/O error.");

		final Snapshot first=runMain(1, 100);
		check(first.value("readsProcessedTotal")==1 && first.value("basesProcessedTotal")==96, "REACHABILITY_FAIL: phase 1 must process the one literal 40/56 pair.");
		check(first.value("correctCountTotal")==0 && first.value("incorrectCountTotal")==0 && sum(first.hist)==0, "REACHABILITY_FAIL: the minimum-length filter must prevent every phase-1 merge.");
		check(first.value(FN)==1 && fnr(out.resolve("phase-1")).equals("100.0000%"), "REACHABILITY_FAIL: truthful insert=40 rejection must produce exactly one false negative and FNR 100%.");
		outputs(1, false);
		phase1="PASSED";

		final Snapshot second=runMain(2, 1);
		check(second.value("readsProcessedTotal")==1 && second.value("basesProcessedTotal")==96, "REACHABILITY_FAIL: phase 2 must process the same single pair.");
		check(second.value("correctCountTotal")==1 && second.value("incorrectCountTotal")==0 && second.value("insertSumCorrectTotal")==40, "REACHABILITY_FAIL: phase 2 must accept the literal insert 40, without fitting inputs or thresholds.");
		check(second.hist.length>40 && second.hist[40]==1 && sum(second.hist)==1 && second.value("insertMinTotal")==40 && second.value("insertMaxTotal")==40, "REACHABILITY_FAIL: phase-2 histogram and extrema must contain exactly insert 40.");
		outputs(2, true);
		final long expected=baseline ? 1 : 0;
		check(second.value(FN)==expected, "Wrong phase-2 false-negative total for "+phase+": "+second.value(FN));
		if(second.value(FN)!=0){discrepancies++;}
		final String rate=fnr(out.resolve("phase-2"));
		check(rate.equals(baseline ? "100.0000%" : "0.0000%"), "Actual phase-2 FNR differs from its independent oracle: "+rate);
		if(!rate.equals("0.0000%")){discrepancies++;}
		phase2="PASSED";

		reset="STARTED"; directResets++;
		BBMerge.resetCounters();
		final Snapshot cleared=snapshot("direct-reset");
		check(cleared.originalHist==second.originalHist && sum(cleared.hist)==0, "Direct reset must clear the existing populated histogram in place.");
		for(final String name : COUNTERS){
			check(cleared.value(name)==(name.equals(FN) ? expected : 0), "Direct-reset control differs for "+name+": "+cleared.value(name));
		}
		check(cleared.value("insertMinTotal")==999999999 && cleared.value("insertMaxTotal")==0, "Direct reset must preserve the established extrema sentinels.");
		if(cleared.value(FN)!=0){discrepancies++;}
		reset="PASSED";
		check(mainCalls==2 && phases==2 && directResets==1, "Keep two real main phases and exactly one final direct reset.");
		check(discrepancies==(baseline ? 3 : 0), "Expected three baseline stale-count/log/reset witnesses, and none after repair.");
		if(reference!=null){compareReference();}
	}

	/** Captures actual main diagnostics, restoring both borrowed process streams on every exit. */
	private Snapshot runMain(final int number, final int minimum) throws Exception{
		final Path work=out.resolve("phase-"+number); Files.createDirectory(work);
		write(work.resolve("input1.fa"), R1); write(work.resolve("input2.fa"), R2);
		final ArrayList<String> arguments=new ArrayList<String>(Arrays.asList(COMMON));
		arguments.addAll(Arrays.asList("in1="+work.resolve("input1.fa"), "in2="+work.resolve("input2.fa"),
			"out="+work.resolve("merged.fa"), "outu1="+work.resolve("unmerged1.fa"), "outu2="+work.resolve("unmerged2.fa"), "minlength="+minimum));
		final StringBuilder command=new StringBuilder("entry\tBBMerge.main in this same JVM\n");
		for(int i=0; i<arguments.size(); i++){command.append(i).append('\t').append(arguments.get(i)).append('\n');}
		write(work.resolve("command.tsv"), command.toString());
		final PrintStream previousOut=System.out, previousErr=System.err;
		try(PrintStream capturedOut=new PrintStream(Files.newOutputStream(work.resolve("main.stdout"), StandardOpenOption.CREATE_NEW), true, "UTF-8");
			PrintStream capturedErr=new PrintStream(Files.newOutputStream(work.resolve("main.stderr"), StandardOpenOption.CREATE_NEW), true, "UTF-8")){
			try{
				System.setOut(capturedOut); System.setErr(capturedErr);
				if(number==1){phase1="STARTED";}else{phase2="STARTED";}
				mainCalls++;
				BBMerge.main(arguments.toArray(new String[arguments.size()]));
				phases++;
				if(number==1){phase1="RETURNED";}else{phase2="RETURNED";}
			}finally{
				System.setOut(previousOut); System.setErr(previousErr);
			}
		}
		check(System.out==previousOut && System.err==previousErr, "Borrowed capture streams must be restored after actual main.");
		check(!field("errorState").getBoolean(null), "Actual main reported a stream error.");
		check(read(work.resolve("input1.fa")).equals(R1) && read(work.resolve("input2.fa")).equals(R2), "Actual main changed literal input bytes.");
		return snapshot("phase-"+number);
	}

	private Snapshot snapshot(final String stage) throws Exception{
		final LinkedHashMap<String, Long> values=new LinkedHashMap<String, Long>();
		for(final String name : COUNTERS){values.put(name, ((Number)field(name).get(null)).longValue());}
		for(final String name : new String[]{"insertMinTotal", "insertMaxTotal"}){values.put(name, ((Number)field(name).get(null)).longValue());}
		final long[] original=(long[])field("histTotal").get(null);
		check(original!=null, "Production histogram must already exist; snapshots never install one.");
		final Snapshot result=new Snapshot(values, original);
		final StringBuilder counters=new StringBuilder("field\tvalue\n"), controls=new StringBuilder();
		for(final Map.Entry<String, Long> entry : values.entrySet()){
			counters.append(entry.getKey()).append('\t').append(entry.getValue()).append('\n');
			if(!entry.getKey().equals(FN)){controls.append(stage).append('\t').append(entry.getKey()).append('\t').append(entry.getValue()).append('\n');}
		}
		write(out.resolve(stage+".counters.tsv"), counters.toString());
		append(out.resolve("controls.tsv"), controls.toString());
		final StringBuilder histogram=new StringBuilder("insert\tcount\n");
		for(int i=0; i<result.hist.length; i++){histogram.append(i).append('\t').append(result.hist[i]).append('\n');}
		write(out.resolve(stage+".histogram.tsv"), histogram.toString());
		return result;
	}

	private void outputs(final int number, final boolean merged) throws IOException{
		final Path work=out.resolve("phase-"+number);
		check(read(work.resolve("merged.fa")).equals(merged ? R1 : ""), "Native merged FASTA must be the literal 40-base fragment, or empty in phase 1.");
		check(read(work.resolve("unmerged1.fa")).equals(merged ? "" : R1), "Native unmerged R1 differs from the literal oracle.");
		check(read(work.resolve("unmerged2.fa")).equals(merged ? "" : R2), "Native unmerged R2 differs from the literal oracle.");
	}

	/** Extracts one printed metric, not sequence records or a reconstructed rate. */
	private String fnr(final Path work) throws IOException{
		String result=null;
		for(final String line : read(work.resolve("main.stderr")).split("\\r?\\n")){
			if(line.startsWith("FNR:")){check(result==null, "Expected exactly one actual FNR line per main phase."); result=line.substring(4).trim();}
		}
		check(result!=null, "Actual custom-header FNR line is missing; do not substitute a calculated value.");
		return result;
	}

	private void compareReference() throws IOException{
		for(final String file : new String[]{"inputs.tsv", "oracles.tsv", "settings.tsv", "controls.tsv"}){same(out.resolve(file), reference.resolve(file));}
		for(final String stage : new String[]{"initial", "phase-1", "phase-2", "direct-reset"}){same(out.resolve(stage+".histogram.tsv"), reference.resolve(stage+".histogram.tsv"));}
		for(final String stage : new String[]{"phase-1", "phase-2"}){
			for(final String file : new String[]{"input1.fa", "input2.fa", "merged.fa", "unmerged1.fa", "unmerged2.fa"}){same(out.resolve(stage).resolve(file), reference.resolve(stage).resolve(file));}
		}
	}

	private void saveDispatch() throws IOException{
		write(out.resolve("dispatch.tsv"), "operation\tstatus\nphase-1\t"+phase1+"\nphase-2\t"+phase2+"\ndirect-reset\t"+reset+"\n");
	}

	private void finish() throws IOException{
		final String summary="BBMERGE_COUNTER_RESET_PASS phase="+phase+" version=1 mainCalls="+mainCalls+" phases="+phases+
			" directResets="+directResets+" discrepancies="+discrepancies+" checks="+checks+" failures=0\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}

	/** Reflection only observes existing production fields; there are no field writes. */
	private static Field field(final String name) throws ReflectiveOperationException{
		final Field field=BBMerge.class.getDeclaredField(name);
		require(Modifier.isStatic(field.getModifiers()), "Snapshot field must be static: "+name);
		field.setAccessible(true); return field;
	}
	private static long sum(final long[] values){long sum=0; for(final long value : values){require(value>=0 && value<=Long.MAX_VALUE-sum, "Finite histogram oracle requires nonnegative counts and a fitting sum."); sum+=value;} return sum;}
	private static String settings(){final StringBuilder b=new StringBuilder(); for(final String arg : COMMON){b.append("both\t").append(arg).append('\n');} return b.append("phase-1\tminlength=100\nphase-2\tminlength=1\n").toString();}
	private void same(final Path a, final Path b) throws IOException{check(Arrays.equals(Files.readAllBytes(a), Files.readAllBytes(b)), "Saved baseline control/input/output changed: "+a.getFileName());}
	private void check(final boolean condition, final String why){checks++; if(!condition){throw new IllegalStateException(why);}}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}
	private static Path path(final String s){require(s!=null && !s.isEmpty(), "Missing required path."); return Paths.get(s).toAbsolutePath().normalize();}
	private static String read(final Path file) throws IOException{return new String(Files.readAllBytes(file), StandardCharsets.US_ASCII);}
	private static void write(final Path file, final String text) throws IOException{Files.write(file, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path file, final String text) throws IOException{Files.write(file, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.APPEND);}

	private static final class Snapshot{
		Snapshot(final Map<String, Long> values_, final long[] hist_){values=values_; originalHist=hist_; hist=hist_.clone();}
		long value(final String name){return values.get(name);}
		final Map<String, Long> values;
		final long[] originalHist, hist;
	}

	private final Path out, reference;
	private final String phase;
	private final boolean baseline;
	private int checks, mainCalls, phases, directResets, discrepancies;
	private String phase1="SKIPPED_PREREQUISITE", phase2="SKIPPED_PREREQUISITE", reset="SKIPPED_PREREQUISITE";
	private static final String FN="unmergedOverlappingTotal";
	private static final String[] COUNTERS={"readsProcessedTotal", "basesProcessedTotal", "matedCountTotal", "correctCountTotal", FN,
		"ambiguousCountTotal", "tooShortCountTotal", "tooLongCountTotal", "incorrectCountTotal", "noSolutionCountTotal",
		"insertSumCorrectTotal", "insertSumIncorrectTotal", "fullyExtendedTotal", "partlyExtendedTotal", "notExtendedTotal", "extensionsAttempted",
		"bridgeAttemptsTotal", "bridgeSuccessTotal", "bridgeNoPathTotal", "bridgeNoTargetTotal", "bridgeAmbiguousTotal", "bridgeCycleTotal",
		"bridgeNonreciprocalTotal", "bridgeRejectedTotal", "adaptersExpected", "adaptersFound"};
	private static final String[] COMMON={"t=1", "simd=f", "bf1=t", "nn=f", "qtrim=f", "maq=0", "ordered=t", "ow=f", "app=f",
		"ratiomode=t", "flatmode=f", "entropy=f", "minoverlap=20", "minoverlap0=20", "ratiominoverlapreduction=0", "mininsert=1", "mininsert0=1",
		"ouq=f", "owq=t", "usequality=t", "efilter=f", "pfilter=0", "mee=0", "tof=0", "extend=0", "extend2=0", "kfilter=0",
		"ecc=f", "join=t", "mix=f", "parsecustom=t", "samplerate=1", "sampleseed=1", "reads=-1", "recalibrate=f", "quantize=f", "train=f"};
	private static final String FRAGMENT="ACGTTGCATGTCAGTACGATCGTACCTGAGCTAGTCGATC";
	private static final String RAW_INPUT="GATCGACTAGCTCAGGTACGATCGTACTGACATGCAACGTAGATCGGAAGAGCACA";
	private static final String R1=">insert=40\n"+FRAGMENT+"\n", R2=">insert=40\n"+RAW_INPUT+"\n";
}
