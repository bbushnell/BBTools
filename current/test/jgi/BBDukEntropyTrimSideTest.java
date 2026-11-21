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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import jgi.BBDuk;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import tracker.EntropyTracker;

/** Five finite checks of legacy BBDuk entropy trimming's parsed left/right/both choice.
 * Each case creates an ordinary parsed owner and one NEW ProcessThread, reads its
 * actual tracker, and calls the real private trim helper once with a fresh unmapped
 * read and null BitSet. At k1/window4/cutoff0.1, only the two AAAA windows in
 * AAAACGTAAAA fail; every mixed window has entropy above0.4. Expected bases and
 * numeric quality slices are independent literals, not another entropy algorithm.
 * No configuration/counter/final fields are seeded; no worker, queue, pipeline,
 * stream or report runs. Normal constructor table allocators may start and join.
 * Two all-low AAAA cases distinguish the retained quality for left-only trimming
 * from the existing both-end minimum-length fallback. Baseline version2 requires
 * the accepted version1 directory and compares every original case record.
 * The helper's existing diagnostic print is retained once per case. This is correctness
 * verification only: no timing loop and no printing-policy or production repair.
 * @author Sucrose
 */
public final class BBDukEntropyTrimSideTest{

	/** phase=baseline|fixed outdir=NEW reference=V1_BASELINE|V2_BASELINE, respectively. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final int split=arg.indexOf('=');
			require(split>0 && split<arg.length()-1, "Use nonempty flag=value arguments: "+arg);
			final String name=arg.substring(0, split), value=arg.substring(split+1);
			require(Arrays.asList("phase", "outdir", "reference").contains(name) && options.put(name, value)==null, "Unknown or duplicate option: "+name);
		}
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbdukentropytrimsidetest.sh with assertions enabled.");
		final String phase=options.get("phase");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		final Path reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require(reference!=null, "Version2 requires the accepted version1 baseline for baseline, or the expanded version2 baseline for fixed.");
		final String expectedReference="baseline".equals(phase) ? "version=1 scenarios=3 owners=3 workers=3 calls=3 controls=1 discrepancies=2 checks=37" : "version=2 scenarios=5 owners=5 workers=5 calls=5 controls=2 discrepancies=3 checks=[0-9]+";
		require(text(reference.resolve("summary.txt")).matches("BBDUK_ENTROPY_TRIM_SIDE_PASS phase=baseline "+expectedReference+" failures=0\\n"), "Reference must match the phase's accepted baseline version and case contract.");
		final BBDukEntropyTrimSideTest test=new BBDukEntropyTrimSideTest(path(options.get("outdir")), reference, phase);
		try{test.exercise();}finally{test.saveDispatch();}
		test.verifyReference();
		test.finish();
	}

	private BBDukEntropyTrimSideTest(final Path out_, final Path reference_, final String phase_) throws Exception{
		out=out_; reference=reference_; phase=phase_; baseline=phase.equals("baseline"); Files.createDirectory(out);
		final Class<?> type=Class.forName("jgi.BBDuk$ProcessThread");
		constructor=type.getDeclaredConstructor(BBDuk.class, ConcurrentReadInputStream.class, ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class, boolean.class); constructor.setAccessible(true);
		trim=type.getDeclaredMethod("trimLowEntropy", Read.class, java.util.BitSet.class, EntropyTracker.class); trim.setAccessible(true);
		trackerField=field(type, "eTrackerT"); workerReads=field(type, "readsEFilteredT"); workerBases=field(type, "basesEFilteredT");
		ownerReads=field(BBDuk.class, "readsEFiltered"); ownerBases=field(BBDuk.class, "basesEFiltered");
		write(out.resolve("settings.tsv"), String.join("\n", SETTINGS)+"\n");
		write(out.resolve("reflection.tsv"), "operation\tmember\nCONSTRUCT_ONLY\tBBDuk.ProcessThread(null streams,localArrays=false)\nREAD_ONLY\tBBDuk.entropyTrim,entropyCutoff,entropyHighpass,calcEntropy,entropyMask,entropyMark,readsEFiltered,basesEFiltered\nREAD_ONLY\tProcessThread.eTrackerT,readsEFilteredT,basesEFilteredT;EntropyTracker.highPass\nINVOKE_SYNC_ONCE_PER_CASE\tProcessThread.trimLowEntropy(Read,null,actual worker tracker)\nFIELD_WRITES\tNONE\n");
		write(out.resolve("setups.tsv"), "scenario\tparsed_side\tcutoff\tcalc_entropy\towner_highpass\ttracker_k\ttracker_window\ttracker_cutoff\ttracker_highpass\n");
		write(out.resolve("observations.tsv"), "scenario\tstage\tworker_state\tworker_reads\tworker_bases\towner_reads\towner_bases\ttracker_left\ttracker_right\ttracker_ns\n");
		write(out.resolve("cases.tsv"), "scenario\tmode\tparsed_side\tremoved\tbases\tqualities\tintended_removed\tintended_bases\tintended_qualities\tworker_state\n");
		write(out.resolve("controls.tsv"), "scenario\tremoved\tbases\tqualities\tworker_state\n");
		final StringBuilder inputs=new StringBuilder("scenario\tmode\tparsed_side\tinput_bases\tinput_qualities\tbaseline_removed\tbaseline_bases\tbaseline_qualities\tfixed_removed\tfixed_bases\tfixed_qualities\n");
		for(int i=0; i<NAMES.length; i++){statuses[i]="SKIPPED"; inputs.append(NAMES[i]).append('\t').append(MODES[i]).append('\t').append(SIDES[i]).append('\t').append(INPUT_BASES[i]).append('\t').append(Arrays.toString(INPUT_QUALITY[i])).append('\t').append(BASELINE_REMOVED[i]).append('\t').append(BASELINE_BASES[i]).append('\t').append(Arrays.toString(BASELINE_QUALITY[i])).append('\t').append(FIXED_REMOVED[i]).append('\t').append(FIXED_BASES[i]).append('\t').append(Arrays.toString(FIXED_QUALITY[i])).append('\n');}
		write(out.resolve("inputs.tsv"), inputs.toString());
	}

	private void exercise() throws Exception{
		for(int i=0; i<NAMES.length; i++){statuses[i]="FAILED"; exercise(i); statuses[i]="PASSED";}
		check(owners==5 && workers==5 && calls==5 && controls==2 && discrepancies==(baseline ? 3 : 0), "Require five separate owners/NEW workers/direct calls, two both-end controls and the phase's three side/quality discrepancies.");
	}

	/** Compare original records during expansion, then the whole expanded setup after repair. */
	private void verifyReference() throws IOException{
		if(baseline){
			for(final String file : new String[]{"settings.tsv", "reflection.tsv"}){same(file);}
			for(final String file : new String[]{"inputs.tsv", "setups.tsv", "cases.tsv", "observations.tsv", "controls.tsv", "dispatch.tsv"}){samePrefix(file);}
			for(int i=0; i<3; i++){same(NAMES[i]+"/constructor-input.fq"); sameConstructorArgs(NAMES[i]);}
		}else{
			for(final String file : new String[]{"settings.tsv", "reflection.tsv", "inputs.tsv", "setups.tsv", "observations.tsv", "controls.tsv", "dispatch.tsv"}){same(file);}
			for(final String name : NAMES){same(name+"/constructor-input.fq"); sameConstructorArgs(name);}
		}
	}

	/** Keep each owner/worker local to its one case; never invoke its processing loop. */
	private void exercise(final int ordinal) throws Exception{
		final String name=NAMES[ordinal]; final Path directory=out.resolve(name); Files.createDirectory(directory);
		write(directory.resolve("constructor-input.fq"), "@constructor_input\n"+INPUT_BASES[ordinal]+"\n+\n"+qualityText(INPUT_QUALITY[ordinal])+"\n");
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(SETTINGS));
		args.add("entropytrim="+MODES[ordinal]); args.add("literal=ACGA"); args.add("in="+directory.resolve("constructor-input.fq"));
		write(directory.resolve("constructor-args.txt"), String.join("\n", args)+"\n");
		final BBDuk owner=new BBDuk(args.toArray(new String[args.size()])); owners++;
		final int side=field(BBDuk.class, "entropyTrim").getInt(owner);
		check(side==SIDES[ordinal] && field(BBDuk.class, "entropyCutoff").getFloat(owner)==0.1f && field(BBDuk.class, "calcEntropy").getBoolean(owner) && field(BBDuk.class, "entropyHighpass").getBoolean(owner) && !field(BBDuk.class, "entropyMask").getBoolean(owner) && !field(BBDuk.class, "entropyMark").getBoolean(owner), "Normal public parsing must retain the requested side and active high-pass cutoff, without masking or marking.");
		final Thread worker=(Thread)constructor.newInstance(owner, null, null, null, null, Boolean.FALSE); workers++;
		final EntropyTracker tracker=(EntropyTracker)trackerField.get(worker);
		check(tracker!=null && tracker.k()==1 && tracker.windowBases()==4 && tracker.cutoff()==0.1f && field(EntropyTracker.class, "highPass").getBoolean(tracker), "Use this normal worker's actual k1/window4/high-pass0.1 tracker; never construct or seed a substitute.");
		append(out.resolve("setups.tsv"), name+"\t"+side+"\t0.1\ttrue\ttrue\t"+tracker.k()+"\t"+tracker.windowBases()+"\t"+tracker.cutoff()+"\ttrue\n");
		observe(name, "before", owner, worker, tracker, -4, -1);
		final Read read=new Read(bytes(INPUT_BASES[ordinal]), INPUT_QUALITY[ordinal].clone(), name, ordinal);
		final int removed=((Integer)trim.invoke(worker, read, null, tracker)).intValue(); calls++;
		final String expectedBases=baseline ? BASELINE_BASES[ordinal] : FIXED_BASES[ordinal];
		final byte[] expectedQuality=baseline ? BASELINE_QUALITY[ordinal] : FIXED_QUALITY[ordinal];
		check(removed==(baseline ? BASELINE_REMOVED[ordinal] : FIXED_REMOVED[ordinal]), "Removed count must honor the phase's literal side oracle; case="+name+" removed="+removed);
		check(Arrays.equals(read.bases, bytes(expectedBases)), "Actual retained bases must match the independent literal slice; case="+name);
		check(Arrays.equals(read.quality, expectedQuality), "Actual numeric qualities must retain exactly the same slice as bases; case="+name);
		check(name.equals(read.id) && read.numericID==ordinal && read.mate==null && read.obj==null && read.match==null && read.samline==null && read.sites==null && !read.discarded(), "The direct helper preserves this fresh unmapped read's identity, null mate/metadata and non-discarded state.");
		observe(name, "after", owner, worker, tracker, INPUT_BASES[ordinal].length()-4, INPUT_BASES[ordinal].length()-1);
		if(removed!=FIXED_REMOVED[ordinal] || !Arrays.equals(read.bases, bytes(FIXED_BASES[ordinal])) || !Arrays.equals(read.quality, FIXED_QUALITY[ordinal])){discrepancies++;}
		append(out.resolve("cases.tsv"), name+"\t"+MODES[ordinal]+"\t"+side+"\t"+removed+"\t"+new String(read.bases, StandardCharsets.US_ASCII)+"\t"+Arrays.toString(read.quality)+"\t"+FIXED_REMOVED[ordinal]+"\t"+FIXED_BASES[ordinal]+"\t"+Arrays.toString(FIXED_QUALITY[ordinal])+"\t"+worker.getState()+"\n");
		if(ordinal==2 || ordinal==4){controls++; append(out.resolve("controls.tsv"), name+"\t"+removed+"\t"+new String(read.bases, StandardCharsets.US_ASCII)+"\t"+Arrays.toString(read.quality)+"\t"+worker.getState()+"\n");}
	}

	private void observe(final String name, final String stage, final BBDuk owner, final Thread worker, final EntropyTracker tracker, final int left, final int right) throws Exception{
		check(worker.getState()==Thread.State.NEW && !worker.isAlive(), "The fixture constructs the worker only; do not start/run/process reads or streams.");
		check(workerReads.getLong(worker)==0 && workerBases.getLong(worker)==0 && ownerReads.getLong(owner)==0 && ownerBases.getLong(owner)==0, "Only the surrounding processing loop aggregates entropy statistics; direct trim calls must leave all four counters zero.");
		check(tracker.leftPos()==left && tracker.rightPos()==right && tracker.ns()==0, "The normal tracker must start empty and then scan every original defined base; case="+name+" stage="+stage);
		append(out.resolve("observations.tsv"), name+"\t"+stage+"\t"+worker.getState()+"\t"+workerReads.getLong(worker)+"\t"+workerBases.getLong(worker)+"\t"+ownerReads.getLong(owner)+"\t"+ownerBases.getLong(owner)+"\t"+tracker.leftPos()+"\t"+tracker.rightPos()+"\t"+tracker.ns()+"\n");
	}
	private void same(final String relative) throws IOException{check(Arrays.equals(Files.readAllBytes(out.resolve(relative)), Files.readAllBytes(reference.resolve(relative))), "Frozen baseline setup/control changed: "+relative);}
	private void samePrefix(final String relative) throws IOException{final String actual=text(out.resolve(relative)), expected=text(reference.resolve(relative)); check(!expected.isEmpty() && actual.length()>expected.length() && actual.startsWith(expected), "The original three-case records must remain an exact prefix: "+relative);}
	private void sameConstructorArgs(final String name) throws IOException{final String actual=text(out.resolve(name+"/constructor-args.txt")).replace("in="+out.resolve(name+"/constructor-input.fq"), "in=constructor-input.fq"), expected=text(reference.resolve(name+"/constructor-args.txt")).replace("in="+reference.resolve(name+"/constructor-input.fq"), "in=constructor-input.fq"); check(actual.equals(expected), "Preserve every normal constructor argument, normalizing only its input path: "+name);}
	/** Encode the literal numeric qualities for the otherwise unused normal constructor FASTQ. */
	private static String qualityText(final byte[] values){final StringBuilder result=new StringBuilder(); for(final byte value : values){assert(value>=0 && value<=93) : "Constructor FASTQ uses Phred+33; quality="+value; result.append((char)(value+33));} return result.toString();}
	private void saveDispatch() throws IOException{final StringBuilder value=new StringBuilder("scenario\tstatus\n"); for(int i=0; i<NAMES.length; i++){value.append(NAMES[i]).append('\t').append(statuses[i]).append('\n');} write(out.resolve("dispatch.tsv"), value.toString());}
	private void finish() throws IOException{final String summary="BBDUK_ENTROPY_TRIM_SIDE_PASS phase="+phase+" version=2 scenarios=5 owners="+owners+" workers="+workers+" calls="+calls+" controls="+controls+" discrepancies="+discrepancies+" checks="+checks+" failures=0\n"; write(out.resolve("summary.txt"), summary); System.out.print(summary);}
	private void check(final boolean condition, final String why){checks++; require(condition, why);}
	private static Field field(final Class<?> type, final String name) throws NoSuchFieldException{final Field result=type.getDeclaredField(name); result.setAccessible(true); return result;}
	private static Path path(final String value){require(value!=null && !value.isEmpty() && value.indexOf('=')<0 && value.indexOf('\t')<0 && value.indexOf('\n')<0, "Use a nonempty path without equals,tabs,newlines for the legacy constructor."); return Paths.get(value).toAbsolutePath().normalize();}
	private static byte[] bytes(final String value){return value.getBytes(StandardCharsets.US_ASCII);}
	private static String text(final Path file) throws IOException{return new String(Files.readAllBytes(file), StandardCharsets.US_ASCII);}
	private static void write(final Path file, final String value) throws IOException{Files.write(file, bytes(value), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path file, final String value) throws IOException{Files.write(file, bytes(value), StandardOpenOption.APPEND);}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}

	private final Path out, reference;
	private final String phase;
	private final boolean baseline;
	private final Constructor<?> constructor;
	private final Method trim;
	private final Field trackerField, workerReads, workerBases, ownerReads, ownerBases;
	private final String[] statuses=new String[5];
	private int owners, workers, calls, controls, discrepancies, checks;
	private static final byte[] QUALITY={30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40};
	private static final byte[] LOW_QUALITY={30, 31, 32, 33};
	private static final String[] INPUT_BASES={"AAAACGTAAAA", "AAAACGTAAAA", "AAAACGTAAAA", "AAAA", "AAAA"};
	private static final byte[][] INPUT_QUALITY={QUALITY, QUALITY, QUALITY, LOW_QUALITY, LOW_QUALITY};
	private static final String[] NAMES={"left", "right", "both", "all_low_left", "all_low_both"}, MODES={"l", "r", "rl", "l", "rl"};
	private static final int[] SIDES={2, 1, 3, 2, 3}, BASELINE_REMOVED={8, 8, 8, 3, 3}, FIXED_REMOVED={4, 4, 8, 3, 3};
	private static final String[] BASELINE_BASES={"CGT", "CGT", "CGT", "A", "A"}, FIXED_BASES={"CGTAAAA", "AAAACGT", "CGT", "A", "A"};
	private static final byte[][] BASELINE_QUALITY={{34, 35, 36}, {34, 35, 36}, {34, 35, 36}, {30}, {30}};
	private static final byte[][] FIXED_QUALITY={{34, 35, 36, 37, 38, 39, 40}, {30, 31, 32, 33, 34, 35, 36}, {34, 35, 36}, {33}, {30}};
	private static final String[] SETTINGS={"k=4", "mink=-1", "edist=0", "hdist=0", "ktrim=f", "maskmiddle=f", "rcomp=f", "findbestmatch=f", "rename=f", "minlen=1", "t=1", "simd=f", "qin=33", "qout=33", "array=t", "prealloc=f", "ow=f", "app=f", "entropy=0.1", "entropyk=1", "entropywindow=4", "entropymask=f", "entropymark=f"};
}
