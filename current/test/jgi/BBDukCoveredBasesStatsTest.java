package test.jgi;

import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
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
import java.util.concurrent.atomic.AtomicLongArray;

import jgi.BBDuk;
import kmer.AbstractKmerTable;
import parse.Parse;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;

/** Finite direct-helper check of legacy BBDuk's covered-base scaffold statistics.
 * Histogram-enabled AA scans return full coverage two, while the no-histogram AA
 * control intentionally returns one at its first threshold crossing. Both should
 * record one read and two input bases. A miss and two single-hit calls on a reused
 * worker are unaffected controls. Four processing workers remain NEW for all five
 * synchronous calls; ordinary constructor table allocators may start and join.
 * Histogram selection uses ordinary outduk arguments. No final or histogram field
 * is assigned reflectively; no processing/loader/stream/report writer runs.
 * @author Sucrose
 */
public final class BBDukCoveredBasesStatsTest{

	/** phase=baseline|fixed outdir=NEW; fixed requires reference=ORIGINAL_BASELINE_CASES. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && !pair[1].isEmpty() && Arrays.asList("phase", "outdir", "reference").contains(pair[0]), "Use documented nonempty flag=value arguments: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate option: "+pair[0]);
		}
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbdukcoveredbasesstatstest.sh with assertions enabled.");
		final String phase=options.get("phase");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		final Path out=path(options.get("outdir")), reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require("baseline".equals(phase) ? reference==null : reference!=null, "Only fixed requires an accepted baseline fixture directory.");
		if(reference!=null){
			require(read(reference.resolve("summary.txt")).matches("BBDUK_COVERED_BASES_STATS_PASS phase=baseline version=1 owners=2 scenarios=4 workers=4 calls=5 controls=3 discrepancies=1 checks=[0-9]+ failures=0\\n"), "Reference must be the successful original covered-base statistics baseline.");
		}
		final BBDukCoveredBasesStatsTest test=new BBDukCoveredBasesStatsTest(out, reference, phase);
		try{test.exercise();}
		finally{test.saveDispatch();}
		test.finish();
	}

	private BBDukCoveredBasesStatsTest(final Path out_, final Path reference_, final String phase_) throws Exception{
		out=out_; reference=reference_; phase=phase_; baseline=phase.equals("baseline");
		Files.createDirectory(out);
		write(out.resolve("settings.tsv"), String.join("\n", SETTINGS)+"\n");
		write(out.resolve("inputs.tsv"), "scenario\tcall\tid\tnumeric_id\tbases\tqualities\tthreshold\n"+
			"histogram_two_hits\t1\tfirst\t0\tAA\tII\t1\nno_histogram_early_return\t1\tfirst\t0\tAA\tII\t1\nhistogram_no_hit\t1\tfirst\t0\tC\tI\t1\n"+
			"histogram_reused_single_hit\t1\tfirst\t0\tA\tI\t1\nhistogram_reused_single_hit\t2\tsecond\t1\tA\tI\t1\n");
		write(out.resolve("oracles.tsv"), "scenario\treturns\tbaseline_reads\tbaseline_bases\tfixed_reads\tfixed_bases\thistogram\n"+
			"histogram_two_hits\t2\t2\t4\t1\t2\tbin2=1,others0\nno_histogram_early_return\t1\t1\t2\t1\t2\tnull\nhistogram_no_hit\t0\t0\t0\t0\t0\tbin0=1,others0\n"+
			"histogram_reused_single_hit\t1,1\t2\t2\t2\t2\tbin1=2,others0\n");
		write(out.resolve("reflection.tsv"), "operation\tmember\nREAD_CONTENTS\tBBDuk.keySets,BBDuk.scaffoldNames\nWRITE_NONFINAL\tBBDuk.storedKmers,BBDuk.scaffoldReadCounts,BBDuk.scaffoldBaseCounts\nCONSTRUCT_ONLY\tBBDuk.ProcessThread\nINVOKE_SYNC\tProcessThread.countCoveredBases\nOBSERVE\tProcessThread.scaffoldReadCountsT,scaffoldBaseCountsT,hitCountsT,BBDuk.hitCounts\n");
		final Class<?> type=Class.forName("jgi.BBDuk$ProcessThread");
		constructor=type.getDeclaredConstructor(BBDuk.class, ConcurrentReadInputStream.class, ConcurrentReadOutputStream.class,
			ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class, boolean.class); constructor.setAccessible(true);
		covered=type.getDeclaredMethod("countCoveredBases", Read.class, AbstractKmerTable[].class, int.class); covered.setAccessible(true);
		readCounts=field(type, "scaffoldReadCountsT"); baseCounts=field(type, "scaffoldBaseCountsT"); histogram=field(type, "hitCountsT"); ownerHistogram=field(BBDuk.class, "hitCounts");
		histogramOn=newOwner(true); histogramOff=newOwner(false);
		write(out.resolve("observations.tsv"), "scenario\tcall\tid\treturn\treads_before\tbases_before\thistogram_before\treads_after\tbases_after\thistogram_after\tstate_before\tstate_after\n");
		final String header="scenario\treturns\tscaffold_reads\tscaffold_bases\thistogram\tpayloads\tworker_state\n";
		write(out.resolve("cases.tsv"), header); write(out.resolve("controls.tsv"), header);
	}

	/** Ordinary options choose histogram allocation; reflection seeds only existing table contents and three nonfinal fields. */
	private Owner newOwner(final boolean enabled) throws Exception{
		final Path directory=out.resolve(enabled ? "histogram_on" : "histogram_off"); Files.createDirectory(directory);
		write(directory.resolve("constructor-input.fq"), "@constructor_input\nAA\n+\nII\n");
		final Path report=directory.resolve("unused-duk.txt");
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(SETTINGS)); args.add("in="+directory.resolve("constructor-input.fq"));
		if(enabled){args.add("outduk="+report);}
		write(directory.resolve("constructor-args.txt"), String.join("\n", args)+"\n");
		final BBDuk tool=new BBDuk(args.toArray(new String[args.size()]));
		final AbstractKmerTable[] sets=(AbstractKmerTable[])field(BBDuk.class, "keySets").get(tool);
		@SuppressWarnings("unchecked")
		final ArrayList<String> names=(ArrayList<String>)field(BBDuk.class, "scaffoldNames").get(tool);
		check(sets.length==7 && names.size()==1 && "".equals(names.get(0)), "Retain normal reserved scaffold ID0 and seven empty tables.");
		long size=0; for(final AbstractKmerTable set : sets){size+=set.size();}
		check(size==0, "No loader may populate reference tables in this direct-helper fixture.");
		names.add("target"); sets[4].set(4L, 1);
		check(sets[4].getValue(4L)==1, "Forward unmasked k1 DNA A is key4 in shard4, scaffold1.");
		final Owner result=new Owner(tool, sets, enabled, report);
		setNonFinal(tool, field(BBDuk.class, "storedKmers"), Long.valueOf(1));
		setNonFinal(tool, field(BBDuk.class, "scaffoldReadCounts"), result.reads);
		setNonFinal(tool, field(BBDuk.class, "scaffoldBaseCounts"), result.bases);
		final long[] counts=(long[])ownerHistogram.get(tool);
		check(enabled ? counts!=null && counts.length==1001 && allZero(counts) : counts==null, "Normal outduk selection alone must allocate or omit the owner histogram.");
		check(!Files.exists(report), "Construction must not create or run the DUK report writer.");
		write(directory.resolve("seed.tsv"), "key\t4\nshard\t4\nscaffold_id\t1\nstored_kmers\t1\nscaffold_name_0\t"+names.get(0)+"\nscaffold_name_1\t"+names.get(1)+"\nscaffold_counter_length\t2\nhistogram_enabled\t"+enabled+"\nhistogram_length\t"+(counts==null ? 0 : counts.length)+"\nfinal_or_histogram_assignments\t0\n");
		return result;
	}

	private void exercise() throws Exception{
		for(int scenario=0; scenario<NAMES.length; scenario++){run(scenario);}
		check(workers==4 && calls==5 && controls==3 && discrepancies==(baseline ? 1 : 0), "Keep four NEW workers, five calls, three unaffected controls and only the expected baseline statistics discrepancy.");
		for(final Owner owner : new Owner[]{histogramOn, histogramOff}){ownerUnchanged(owner);}
		if(reference!=null){
			for(final String file : new String[]{"settings.tsv", "inputs.tsv", "oracles.tsv", "reflection.tsv", "controls.tsv",
				"histogram_on/seed.tsv", "histogram_off/seed.tsv", "histogram_on/constructor-input.fq", "histogram_off/constructor-input.fq"}){
				check(Arrays.equals(Files.readAllBytes(out.resolve(file)), Files.readAllBytes(reference.resolve(file))), "Baseline literal/setup/control changed: "+file);
			}
		}
	}

	/** Five calls total; histogram and no-histogram coverage returns intentionally differ. */
	private void run(final int scenario) throws Exception{
		statuses[scenario]="FAILED";
		final Owner owner=scenario==1 ? histogramOff : histogramOn;
		final Thread worker=(Thread)constructor.newInstance(owner.tool, null, null, null, null, Boolean.TRUE); workers++;
		check(worker.getState()==Thread.State.NEW && !worker.isAlive(), "Construct processing workers without starting them.");
		final StringBuilder returns=new StringBuilder(), payloads=new StringBuilder();
		Snapshot after=snapshot(worker, owner.enabled);
		check(allZero(after.reads) && allZero(after.bases) && (after.histogram==null || allZero(after.histogram)), "Every new worker must start with zero local statistics and histogram.");
		for(int call=0; call<SEQUENCES[scenario].length; call++){
			final String sequence=SEQUENCES[scenario][call], id=call==0 ? "first" : "second";
			final byte[] expectedBases=bytes(sequence), qualities=new byte[sequence.length()]; Arrays.fill(qualities, (byte)40);
			final Read read=new Read(expectedBases.clone(), qualities.clone(), id, call);
			payload(read, id, call, expectedBases, qualities);
			check((int)Math.ceil(.5*sequence.length())==1, "Literal threshold one must agree with the ordinary mcf=.5 caller.");
			final Snapshot before=snapshot(worker, owner.enabled);
			final int result;
			try{result=((Integer)covered.invoke(worker, read, owner.sets, Integer.valueOf(1))).intValue();}
			catch(final InvocationTargetException e){throw new IllegalStateException("Unexpected covered-base helper failure in "+NAMES[scenario], e.getCause());}
			calls++; after=snapshot(worker, owner.enabled);
			append(out.resolve("observations.tsv"), NAMES[scenario]+"\t"+(call+1)+"\t"+id+"\t"+result+"\t"+Arrays.toString(before.reads)+"\t"+Arrays.toString(before.bases)+"\t"+Arrays.toString(before.histogram)+
				"\t"+Arrays.toString(after.reads)+"\t"+Arrays.toString(after.bases)+"\t"+Arrays.toString(after.histogram)+"\t"+before.state+"\t"+after.state+"\n");
			payload(read, id, call, expectedBases, qualities);
			check(result==COVERAGE[scenario][call], "Preserve full histogram-enabled coverage and intentional no-histogram early return.");
			final int fixedReads=READS[scenario][call], fixedBases=BASES[scenario][call];
			final boolean witness=baseline && scenario==0;
			check(Arrays.equals(after.reads, new long[]{0, witness ? 2 : fixedReads}) && Arrays.equals(after.bases, new long[]{0, witness ? 4 : fixedBases}), "Literal per-read scaffold statistics differ from the selected baseline/fixed expectation.");
			if(!Arrays.equals(after.reads, new long[]{0, fixedReads}) || !Arrays.equals(after.bases, new long[]{0, fixedBases})){discrepancies++;}
			final long[] expectedHistogram=owner.enabled ? new long[1001] : null;
			if(expectedHistogram!=null){expectedHistogram[HISTOGRAM_BIN[scenario]]=call+1;}
			check(Arrays.equals(after.histogram, expectedHistogram), "All 1001 histogram bins must count the full covered-base result exactly once per completed call.");
			ownerUnchanged(owner);
			if(call>0){returns.append(','); payloads.append('|');}
			returns.append(result); payloads.append(id).append(':').append(sequence).append(':'); for(int i=0; i<qualities.length; i++){payloads.append('I');}
		}
		final String row=NAMES[scenario]+"\t"+returns+"\t"+Arrays.toString(after.reads)+"\t"+Arrays.toString(after.bases)+"\t"+Arrays.toString(after.histogram)+"\t"+payloads+"\t"+after.state+"\n";
		append(out.resolve("cases.tsv"), row); if(scenario>0){append(out.resolve("controls.tsv"), row); controls++;} statuses[scenario]="PASSED";
	}

	private Snapshot snapshot(final Thread worker, final boolean enabled) throws IllegalAccessException{
		final long[] reads=(long[])readCounts.get(worker), bases=(long[])baseCounts.get(worker), counts=(long[])histogram.get(worker);
		check(worker.getState()==Thread.State.NEW && !worker.isAlive(), "The helper must execute synchronously while its processing worker remains NEW.");
		check(reads!=null && reads.length==2 && bases!=null && bases.length==2 && (enabled ? counts!=null && counts.length==1001 : counts==null), "Use two-scaffold local arrays and the normally selected 1001-bin histogram or null.");
		return new Snapshot(reads.clone(), bases.clone(), counts==null ? null : counts.clone(), worker.getState());
	}
	private void ownerUnchanged(final Owner owner) throws IllegalAccessException{
		check(owner.reads.get(0)==0 && owner.reads.get(1)==0 && owner.bases.get(0)==0 && owner.bases.get(1)==0, "Never aggregate local worker statistics into owner counters or run a pipeline.");
		final long[] counts=(long[])ownerHistogram.get(owner.tool);
		check(owner.enabled ? counts!=null && counts.length==1001 && allZero(counts) : counts==null, "Owner histogram must remain untouched by direct worker-local calls.");
		check(!Files.exists(owner.report), "The finite helper fixture must not write a DUK report.");
	}
	private void payload(final Read read, final String id, final long numeric, final byte[] bases, final byte[] qualities){check(id.equals(read.id) && read.numericID==numeric && read.mate==null && Arrays.equals(read.bases, bases) && Arrays.equals(read.quality, qualities), "Preserve literal read identity, pairing, bases and Q40 qualities.");}
	private void setNonFinal(final BBDuk tool, final Field field, final Object value) throws IllegalAccessException{check(!Modifier.isFinal(field.getModifiers()) && Arrays.asList("storedKmers", "scaffoldReadCounts", "scaffoldBaseCounts").contains(field.getName()), "Only three named NONFINAL owner fields may be seeded: "+field.getName()); field.set(tool, value);}
	private void saveDispatch() throws IOException{final StringBuilder text=new StringBuilder("scenario\tstatus\n"); for(int i=0; i<NAMES.length; i++){text.append(NAMES[i]).append('\t').append(statuses[i]).append('\n');} write(out.resolve("dispatch.tsv"), text.toString());}
	private void finish() throws IOException{final String summary="BBDUK_COVERED_BASES_STATS_PASS phase="+phase+" version=1 owners=2 scenarios=4 workers="+workers+" calls="+calls+" controls="+controls+" discrepancies="+discrepancies+" checks="+checks+" failures=0\n"; write(out.resolve("summary.txt"), summary); System.out.print(summary);}
	private void check(final boolean condition, final String why){checks++; if(!condition){throw new IllegalStateException(why);}}
	private static boolean allZero(final long[] values){for(final long value : values){if(value!=0){return false;}}return true;}
	private static Field field(final Class<?> type, final String name) throws NoSuchFieldException{final Field field=type.getDeclaredField(name); field.setAccessible(true); return field;}
	private static Path path(final String text){require(text!=null && !text.isEmpty() && text.indexOf('=')<0 && text.indexOf('\n')<0 && text.indexOf('\t')<0, "Use a nonempty path without equals, tabs or newlines for the legacy constructor."); return Paths.get(text).toAbsolutePath().normalize();}
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static String read(final Path file) throws IOException{return new String(Files.readAllBytes(file), StandardCharsets.US_ASCII);}
	private static void write(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.APPEND);}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}

	private static final class Owner{
		Owner(final BBDuk tool_, final AbstractKmerTable[] sets_, final boolean enabled_, final Path report_){tool=tool_; sets=sets_; enabled=enabled_; report=report_;}
		final BBDuk tool; final AbstractKmerTable[] sets; final boolean enabled; final Path report;
		final AtomicLongArray reads=new AtomicLongArray(2), bases=new AtomicLongArray(2);
	}
	private static final class Snapshot{
		Snapshot(final long[] reads_, final long[] bases_, final long[] histogram_, final Thread.State state_){reads=reads_; bases=bases_; histogram=histogram_; state=state_;}
		final long[] reads, bases, histogram; final Thread.State state;
	}
	private final Path out, reference;
	private final String phase;
	private final boolean baseline;
	private final Constructor<?> constructor;
	private final Method covered;
	private final Field readCounts, baseCounts, histogram, ownerHistogram;
	private final Owner histogramOn, histogramOff;
	private int workers, calls, controls, discrepancies, checks;
	private final String[] statuses={"SKIPPED_PREREQUISITE", "SKIPPED_PREREQUISITE", "SKIPPED_PREREQUISITE", "SKIPPED_PREREQUISITE"};
	private static final String[] NAMES={"histogram_two_hits", "no_histogram_early_return", "histogram_no_hit", "histogram_reused_single_hit"};
	private static final String[][] SEQUENCES={{"AA"}, {"AA"}, {"C"}, {"A", "A"}};
	private static final int[][] COVERAGE={{2}, {1}, {0}, {1, 1}}, READS={{1}, {1}, {0}, {1, 2}}, BASES={{2}, {2}, {0}, {1, 2}};
	private static final int[] HISTOGRAM_BIN={2, -1, 0, 1};
	private static final String[] SETTINGS={"k=1", "maskmiddle=f", "rcomp=f", "qskip=1", "speed=0", "qhdist=0", "hdist=0", "edist=0", "useshortkmers=f", "mcf=0.5",
		"findbestmatch=f", "rename=f", "skipr1=f", "skipr2=f", "restrictleft=0", "restrictright=0", "literal=A", "minlen=1", "t=1", "simd=f", "qin=33", "qout=33", "prealloc=f", "ow=f", "app=f"};
}
