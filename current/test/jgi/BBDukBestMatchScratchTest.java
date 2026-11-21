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
import structures.IntList;

/** Finite direct-helper check of legacy BBDuk's per-worker best-match scratch lifetime.
 * Three processing workers are constructed but never started: five synchronous helper
 * calls cover A then AA at threshold one, fresh AA, and A then AA at threshold zero.
 * Reflection is limited to the named existing constructor/helper and seed/observation
 * fields. Final options use the ordinary constructor; no final field is assigned.
 * Normal constructor table allocators may start and join. No processing/loader/stream
 * pipeline is created or run. Scratch observations are not fixed-phase cleanup rules.
 * @author Sucrose
 */
public final class BBDukBestMatchScratchTest{

	/** phase=baseline|fixed outdir=NEW; fixed requires reference=ACCEPTED_BASELINE. */
	public static void main(final String[] args) throws Exception{
		final Map<String, String> options=new LinkedHashMap<String, String>();
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && Arrays.asList("phase", "outdir", "reference").contains(pair[0]), "Use documented flag=value arguments: "+arg);
			require(options.put(pair[0], pair[1])==null, "Duplicate option: "+pair[0]);
		}
		boolean assertions=false; assert(assertions=true);
		require(assertions, "Use dev/bbdukbestmatchscratchtest.sh with assertions enabled.");
		final String phase=options.get("phase");
		require("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline|fixed.");
		final Path out=path(options.get("outdir")), reference=options.containsKey("reference") ? path(options.get("reference")) : null;
		require("baseline".equals(phase) ? reference==null : reference!=null, "Only fixed requires an accepted baseline reference.");
		if(reference!=null){
			require(read(reference.resolve("summary.txt")).matches("BBDUK_BEST_MATCH_SCRATCH_PASS phase=baseline version=1 scenarios=3 workers=3 calls=5 controls=2 discrepancies=1 expectedExceptions=1 checks=[0-9]+ failures=0\\n"), "Reference must be the successful finite scratch baseline.");
		}
		final BBDukBestMatchScratchTest test=new BBDukBestMatchScratchTest(out, reference, phase);
		try{test.exercise();}
		finally{test.saveDispatch();}
		test.finish();
	}

	private BBDukBestMatchScratchTest(final Path out_, final Path reference_, final String phase_) throws Exception{
		out=out_; reference=reference_; phase=phase_; baseline=phase.equals("baseline");
		Files.createDirectory(out);
		write(out.resolve("constructor-input.fq"), "@constructor_input\nA\n+\nI\n");
		write(out.resolve("inputs.tsv"), "id\tbases\tqualities\nfirst\tA\tI\nsecond\tAA\tII\n");
		write(out.resolve("oracles.tsv"), "scenario\tthreshold\tbaseline_returns\tfixed_returns\tbaseline_scaffold_reads\tbaseline_scaffold_bases\tfixed_scaffold_reads\tfixed_scaffold_bases\nreused_below\t1\t-1,EXCEPTION\t-1,1\t0\t0\t1\t2\nfresh_above\t1\t1\t1\t1\t2\t1\t2\nreused_above\t0\t1,1\t1,1\t2\t3\t2\t3\n");
		write(out.resolve("settings.tsv"), String.join("\n", SETTINGS)+"\n");
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList(SETTINGS));
		args.add("in="+out.resolve("constructor-input.fq"));
		write(out.resolve("constructor-args.txt"), String.join("\n", args)+"\n");
		owner=new BBDuk(args.toArray(new String[args.size()]));
		workerType=Class.forName("jgi.BBDuk$ProcessThread");
		constructor=workerType.getDeclaredConstructor(BBDuk.class, ConcurrentReadInputStream.class, ConcurrentReadOutputStream.class,
			ConcurrentReadOutputStream.class, ConcurrentReadOutputStream.class, boolean.class);
		constructor.setAccessible(true);
		find=workerType.getDeclaredMethod("findBestMatch", Read.class, AbstractKmerTable[].class, int.class); find.setAccessible(true);
		countArray=field(workerType, "countArray"); idList=field(workerType, "idList"); countList=field(workerType, "countList");
		readCounts=field(workerType, "scaffoldReadCountsT"); baseCounts=field(workerType, "scaffoldBaseCountsT");
		final Field tables=field(BBDuk.class, "keySets"), namesField=field(BBDuk.class, "scaffoldNames");
		sets=(AbstractKmerTable[])tables.get(owner);
		@SuppressWarnings("unchecked")
		final ArrayList<String> names=(ArrayList<String>)namesField.get(owner);
		check(sets.length==7 && names.size()==1 && "".equals(names.get(0)), "Retain the normal constructor's single empty scaffold name at reserved ID0.");
		long size=0; for(final AbstractKmerTable set : sets){size+=set.size();}
		check(size==0, "No reference loader may have populated tables before this direct-helper fixture.");
		names.add("target"); sets[4].set(4L, 1);
		check(sets[4].getValue(4L)==1, "Literal DNA A with k1/rcomp=f/maskmiddle=f is key4 in shard4, scaffold1.");
		outerReads=new AtomicLongArray(2); outerBases=new AtomicLongArray(2);
		setNonFinal(field(BBDuk.class, "storedKmers"), Long.valueOf(1));
		setNonFinal(field(BBDuk.class, "scaffoldReadCounts"), outerReads);
		setNonFinal(field(BBDuk.class, "scaffoldBaseCounts"), outerBases);
		write(out.resolve("seed.tsv"), "key\t4\nshard\t4\nscaffold_id\t1\nstored_kmers\t1\nscaffold_name_0\t"+names.get(0)+"\nscaffold_name_1\t"+names.get(1)+"\nscaffold_counter_length\t2\nlocal_arrays\ttrue\nfinal_field_assignments\t0\n");
		write(out.resolve("reflection.tsv"), "operation\tmember\nREAD_CONTENTS\tBBDuk.keySets\nREAD_CONTENTS\tBBDuk.scaffoldNames\nWRITE_NONFINAL\tBBDuk.storedKmers\nWRITE_NONFINAL\tBBDuk.scaffoldReadCounts\nWRITE_NONFINAL\tBBDuk.scaffoldBaseCounts\nCONSTRUCT_ONLY\tBBDuk.ProcessThread\nINVOKE_SYNC\tProcessThread.findBestMatch\nOBSERVE\tProcessThread.countArray,idList,countList,scaffoldReadCountsT,scaffoldBaseCountsT\n");
		write(out.resolve("observations.tsv"), "scenario\tcall\tid\tthreshold\tresult\texception\tcounts_before\tids_before\thits_before\tcounts_after\tids_after\thits_after\tscaffold_reads\tscaffold_bases\tstate_before\tstate_after\n");
		write(out.resolve("cases.tsv"), "scenario\treturns\tscaffold_reads\tscaffold_bases\tpayloads\tworker_state\n");
		write(out.resolve("controls.tsv"), "scenario\treturns\tscaffold_reads\tscaffold_bases\tpayloads\tworker_state\n");
	}

	private void exercise() throws Exception{
		for(int scenario=0; scenario<3; scenario++){run(scenario);}
		check(workers==3 && calls==5 && controls==2, "Keep three nonstarted workers, five calls and two unaffected control scenarios.");
		check(expectedExceptions==(baseline ? 1 : 0) && discrepancies==(baseline ? 1 : 0), "Only the second reused-below-threshold baseline call may fail its independent oracle.");
		check(outerReads.get(0)==0 && outerReads.get(1)==0 && outerBases.get(0)==0 && outerBases.get(1)==0, "Direct local-array calls must not aggregate into owner counters or run the processing pipeline.");
		if(reference!=null){
			for(final String file : new String[]{"constructor-input.fq", "inputs.tsv", "oracles.tsv", "settings.tsv", "seed.tsv", "reflection.tsv", "controls.tsv"}){
				check(Arrays.equals(Files.readAllBytes(out.resolve(file)), Files.readAllBytes(reference.resolve(file))), "Baseline literal/setup/control changed: "+file);
			}
		}
	}

	/** Each call executes only the existing helper on the fixture's main thread. */
	private void run(final int scenario) throws Exception{
		final Thread worker=(Thread)constructor.newInstance(owner, null, null, null, null, Boolean.TRUE); workers++;
		check(worker.getState()==Thread.State.NEW && !worker.isAlive(), "A constructed processing worker must remain NEW, never started.");
		final StringBuilder returns=new StringBuilder(), payloads=new StringBuilder();
		Snapshot after=null;
		for(int call=0; call<SEQUENCES[scenario].length; call++){
			final String bases=SEQUENCES[scenario][call], id=bases.length()==1 ? "first" : "second";
			final byte[] qualities=new byte[bases.length()]; Arrays.fill(qualities, (byte)40);
			final Read read=new Read(bytes(bases), qualities.clone(), id, bases.length()-1L);
			final Snapshot before=snapshot(worker);
			check(before.state==Thread.State.NEW && !worker.isAlive(), "The actual helper must be called before any processing thread starts.");
			Integer result=null; Throwable failure=null;
			try{result=(Integer)find.invoke(worker, read, sets, Integer.valueOf(THRESHOLDS[scenario]));}
			catch(final InvocationTargetException e){failure=e.getCause();}
			calls++; after=snapshot(worker);
			final String returned=failure==null ? String.valueOf(result) : "EXCEPTION";
			append(out.resolve("observations.tsv"), NAMES[scenario]+"\t"+(call+1)+"\t"+id+"\t"+THRESHOLDS[scenario]+"\t"+returned+"\t"+(failure==null ? "-" : failure.getClass().getName())+
				"\t"+Arrays.toString(before.counts)+"\t"+before.ids+"\t"+before.hits+"\t"+Arrays.toString(after.counts)+"\t"+after.ids+"\t"+after.hits+
				"\t"+Arrays.toString(after.reads)+"\t"+Arrays.toString(after.bases)+"\t"+before.state+"\t"+after.state+"\n");
			check(after.state==Thread.State.NEW && !worker.isAlive(), "Synchronous helper return or exception must leave the worker NEW.");
			check(id.equals(read.id) && read.numericID==bases.length()-1L && read.mate==null, "Helper changed the literal read identity or pairing.");
			check(Arrays.equals(read.bases, bytes(bases)) && Arrays.equals(read.quality, qualities), "Helper changed literal bases or Q40 qualities with rename=f.");
			final boolean witness=(scenario==0 && call==1);
			if(failure!=null){
				if(!(baseline && witness && failure instanceof ArrayIndexOutOfBoundsException)){
					throw new IllegalStateException("Unexpected direct-helper failure in "+NAMES[scenario]+" call "+(call+1), failure);
				}
				final StackTraceElement[] trace=failure.getStackTrace();
				check(trace.length>0 && trace[0].getClassName().equals(workerType.getName()) && trace[0].getMethodName().equals("findBestMatch"), "Expected failure must originate in the actual helper's scaffold update, not reflection/setup.");
				final StringBuilder text=new StringBuilder(failure.toString()+"\n"); for(final StackTraceElement frame : trace){text.append(frame).append('\n');}
				write(out.resolve("expected-exception.txt"), text.toString()); expectedExceptions++;
				check(before.counts[1]==1 && before.ids.equals("[1]") && after.counts[1]==3 && after.ids.equals("[]") && after.hits.equals("[]"), "Baseline must expose dropped touched-ID bookkeeping after the one-hit read.");
			}else{
				check(!(baseline && witness), "Baseline must reproduce the two-read stale-scratch witness before any repair.");
				check(result!=null && result.intValue()==IDS[scenario][call], "Literal per-read best-scaffold oracle mismatch.");
			}
			if(failure!=null || result.intValue()!=IDS[scenario][call]){discrepancies++;}
			final int readsExpected=(baseline && witness ? 0 : READ_COUNTS[scenario][call]);
			final int basesExpected=(baseline && witness ? 0 : BASE_COUNTS[scenario][call]);
			check(Arrays.equals(after.reads, new long[]{0, readsExpected}) && Arrays.equals(after.bases, new long[]{0, basesExpected}), "Literal local scaffold read/base counts differ.");
			if(call>0){returns.append(','); payloads.append('|');}
			returns.append(returned); payloads.append(id).append(':').append(bases).append(':');
			for(int i=0; i<qualities.length; i++){payloads.append('I');}
		}
		final String row=NAMES[scenario]+"\t"+returns+"\t"+Arrays.toString(after.reads)+"\t"+Arrays.toString(after.bases)+"\t"+payloads+"\t"+after.state+"\n";
		append(out.resolve("cases.tsv"), row);
		if(scenario>0){append(out.resolve("controls.tsv"), row); controls++;}
		statuses[scenario]="PASSED";
	}

	/** Reflection assigns only the three named nonfinal owner fields, never a final option/reference. */
	private void setNonFinal(final Field field, final Object value) throws IllegalAccessException{
		check(!Modifier.isFinal(field.getModifiers()) && Arrays.asList("storedKmers", "scaffoldReadCounts", "scaffoldBaseCounts").contains(field.getName()), "Refusing an unreleased or final field assignment: "+field.getName());
		field.set(owner, value);
	}
	private Snapshot snapshot(final Thread worker) throws IllegalAccessException{
		final int[] counts=(int[])countArray.get(worker);
		final long[] reads=(long[])readCounts.get(worker), bases=(long[])baseCounts.get(worker);
		check(counts!=null && counts.length==2 && reads!=null && reads.length==2 && bases!=null && bases.length==2, "Use only two-scaffold local scratch and counter arrays.");
		return new Snapshot(counts.clone(), list((IntList)idList.get(worker)), list((IntList)countList.get(worker)), reads.clone(), bases.clone(), worker.getState());
	}
	private static Field field(final Class<?> type, final String name) throws NoSuchFieldException{
		final Field field=type.getDeclaredField(name); field.setAccessible(true); return field;
	}
	private static String list(final IntList list){final StringBuilder text=new StringBuilder("["); for(int i=0; i<list.size; i++){if(i>0){text.append(',');}text.append(list.get(i));} return text.append(']').toString();}
	private void saveDispatch() throws IOException{final StringBuilder text=new StringBuilder("scenario\tstatus\n"); for(int i=0; i<NAMES.length; i++){text.append(NAMES[i]).append('\t').append(statuses[i]).append('\n');} write(out.resolve("dispatch.tsv"), text.toString());}
	private void finish() throws IOException{
		final String summary="BBDUK_BEST_MATCH_SCRATCH_PASS phase="+phase+" version=1 scenarios=3 workers="+workers+" calls="+calls+" controls="+controls+
			" discrepancies="+discrepancies+" expectedExceptions="+expectedExceptions+" checks="+checks+" failures=0\n";
		write(out.resolve("summary.txt"), summary); System.out.print(summary);
	}
	private void check(final boolean condition, final String why){checks++; if(!condition){throw new IllegalStateException(why);}}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}
	private static Path path(final String text){require(text!=null && !text.isEmpty() && text.indexOf('=')<0, "Set nonempty paths without '=' for the legacy constructor."); return Paths.get(text).toAbsolutePath().normalize();}
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.US_ASCII);}
	private static String read(final Path file) throws IOException{return new String(Files.readAllBytes(file), StandardCharsets.US_ASCII);}
	private static void write(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.CREATE_NEW);}
	private static void append(final Path file, final String text) throws IOException{Files.write(file, bytes(text), StandardOpenOption.APPEND);}

	private static final class Snapshot{
		Snapshot(final int[] counts_, final String ids_, final String hits_, final long[] reads_, final long[] bases_, final Thread.State state_){counts=counts_; ids=ids_; hits=hits_; reads=reads_; bases=bases_; state=state_;}
		final int[] counts; final String ids, hits; final long[] reads, bases; final Thread.State state;
	}
	private final Path out, reference;
	private final String phase;
	private final boolean baseline;
	private final BBDuk owner;
	private final Class<?> workerType;
	private final Constructor<?> constructor;
	private final Method find;
	private final Field countArray, idList, countList, readCounts, baseCounts;
	private final AbstractKmerTable[] sets;
	private final AtomicLongArray outerReads, outerBases;
	private int checks, workers, calls, controls, discrepancies, expectedExceptions;
	private final String[] statuses={"SKIPPED_PREREQUISITE", "SKIPPED_PREREQUISITE", "SKIPPED_PREREQUISITE"};
	private static final String[] NAMES={"reused_below", "fresh_above", "reused_above"};
	private static final String[][] SEQUENCES={{"A", "AA"}, {"AA"}, {"A", "AA"}};
	private static final int[] THRESHOLDS={1, 1, 0};
	private static final int[][] IDS={{-1, 1}, {1}, {1, 1}}, READ_COUNTS={{0, 1}, {1}, {1, 2}}, BASE_COUNTS={{0, 2}, {2}, {1, 3}};
	private static final String[] SETTINGS={"k=1", "maskmiddle=f", "rcomp=f", "qskip=1", "speed=0", "qhdist=0", "hdist=0", "edist=0", "useshortkmers=f",
		"findbestmatch=t", "rename=f", "skipr1=f", "skipr2=f", "restrictleft=0", "restrictright=0", "maxbadkmers=1", "literal=A", "minlen=1", "t=1", "simd=f", "qin=33", "qout=33", "prealloc=f"};
}
