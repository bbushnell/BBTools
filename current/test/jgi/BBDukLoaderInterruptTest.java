package test.jgi;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.atomic.AtomicLongArray;

import jgi.BBDuk;
import kmer.AbstractKmerTable;
import shared.KillSwitch;
import structures.IntList;

/** Two literal-only loader baselines, one fresh JVM per case.
 * Ordinary construction precedes one reflective spawnLoadThreads invocation.
 * The interrupted owner fails at its first literal put; finally always discovers
 * exact class/owner-matched live loaders, verifies empty queues, offers the existing
 * identity POISON and joins those exact handles. No queues are cleared and no worker
 * is interrupted. Observation/reporting failures cannot bypass that cleanup.
 * Unexpected fixture or cleanup Throwable prints and immediately halts this dedicated
 * JVM through KillSwitch; it never returns while abandoning its live workers.
 * Normal completion reports no live owned loaders, not retained terminated handles.
 * No reference input, process pipeline, field seeding, timeout or hang replay.
 * @author Sucrose
 */
public final class BBDukLoaderInterruptTest{

	/** Use scenario=normal|owner_interrupted and a new outdir via the dev launcher. */
	public static void main(final String[] args){
		try{
			String scenario=null, directory=null;
			for(final String arg : args){
				final int at=arg.indexOf('='); require(at>0 && at<arg.length()-1, "Use nonempty flag=value: "+arg);
				final String key=arg.substring(0, at), value=arg.substring(at+1);
				if(key.equals("scenario")){require(scenario==null, "Duplicate scenario."); scenario=value;}
				else if(key.equals("outdir")){require(directory==null, "Duplicate outdir."); directory=value;}
				else{throw new IllegalArgumentException("Unknown option: "+arg);}
			}
			require(args.length==2 && directory!=null && ("normal".equals(scenario) || "owner_interrupted".equals(scenario)), "Use exactly scenario and outdir.");
			boolean assertions=false; assert(assertions=true); require(assertions, "The dev launcher must enable assertions.");
			new BBDukLoaderInterruptTest(Paths.get(directory).toAbsolutePath().normalize(), scenario).run();
		}catch(final Throwable failure){abort("Loader-interruption baseline fixture failed; aborting its dedicated JVM.", failure);}
	}

	private BBDukLoaderInterruptTest(final Path out_, final String name_) throws Exception{
		out=out_; name=name_; interrupted=name.equals("owner_interrupted"); input=out.resolve("input.fq");
		require(out.toString().indexOf('=')<0 && out.toString().indexOf('\t')<0 && out.toString().indexOf('\n')<0, "Use a plain reproducible directory.");
		Files.createDirectory(out); write("input.fq", INPUT); write("settings.tsv", lines(SETTINGS)); write("counter-schema.tsv", lines(COUNTERS));
		loadType=Class.forName("jgi.BBDuk$LoadThread"); enclosing=member(loadType, "this$0");
		queueField=member(loadType, "queue"); partitionField=member(loadType, "tnum");
		write("reflection.tsv", "operation\tmember\nCONSTRUCTION\tordinary public BBDuk arguments\nCALL\tprivate spawnLoadThreads exactly once\nOWNERSHIP\texact runtime class and this$0 identity equal to fresh owner\nEVENT\tonly calling fixture thread interrupted after construction\nCLEANUP\tfinally; empty owned queues; identity POISON offer; exact-handle join\nWORKER_INTERRUPTS\t0\nQUEUE_CLEARS\t0\nBBDUK_FIELD_WRITES\t0\nPROCESS_CALLS\t0\nINPUT_FACTORIES\t0; literal only\nFAILURE_POLICY\tunexpected fixture or cleanup Throwable prints and immediately halts own JVM with exit1\n");
	}

	private void run() throws Exception{
		final ArrayList<String> options=new ArrayList<String>(Arrays.asList(SETTINGS)); options.add("literal=ACGT"); options.add("in="+input);
		final String[] raw=options.toArray(new String[options.size()]), passed=raw.clone(); write("constructor-args.txt", lines(raw));
		Throwable caught=null; Long returned=null; boolean flagBefore=false, flagAfter=false;
		String afterState=null;
		try{
			owner=new BBDuk(passed);
			check(!Thread.currentThread().isInterrupted() && number("k")==4 && number("mink")==4 && !bool("useShortKmers") && !bool("maskMiddle") && !bool("rcomp") && number("THREADS")==1 && field("ref").get(owner)==null && Arrays.equals((String[])field("literal").get(owner), new String[]{"ACGT"}), "Use an ordinary fresh literal-only owner before delivering the interrupt.");
			check(state(0, 0, false), "Constructor state must have seven empty tables, dummy scaffold and zero counters.");
			write("counters-before.tsv", counters());
			check(discover().isEmpty(), "Fresh owner must have no live LoadThreads before the one call.");
			final Method loader=BBDuk.class.getDeclaredMethod("spawnLoadThreads"); loader.setAccessible(true);
			if(interrupted){Thread.currentThread().interrupt();}
			flagBefore=Thread.currentThread().isInterrupted(); calls++;
			try{returned=(Long)loader.invoke(owner);}catch(final InvocationTargetException invocation){caught=invocation.getCause();}
			flagAfter=Thread.currentThread().isInterrupted();
			check(flagBefore==interrupted && !flagAfter, "The first put consumes only the owner's delivered interruption; before="+flagBefore+" after="+flagAfter);
			check(interrupted ? exactException(caught) && returned==null : caught==null && Long.valueOf(1).equals(returned), "Require the exact first literal-put exception or normal returned count1, never a seeded storedKmers value.");
			captured.putAll(discover()); observed=captured.size();
			check(observed==(interrupted ? 7 : 0) && emptyOwned(captured), "Interrupted first put must leave seven exact owner/class-matched empty queues; normal completion exposes no live handles.");
			check(state(interrupted ? 0 : 1, interrupted ? 0 : 1, true), "Require exact tables/ref counters/scaffolds/atomic arrays after the real loader call.");
			afterState=counters()+tablesAndScaffolds(); write("counters-after.tsv", counters()); write("state.tsv", tablesAndScaffolds());
			check(Arrays.equals(raw, passed) && Arrays.equals(Files.readAllBytes(input), INPUT.getBytes(StandardCharsets.US_ASCII)) && input.toString().equals(field("in1").get(owner)) && field("in2").get(owner)==null && field("out1").get(owner)==null && field("out2").get(owner)==null && field("sidechannel").get(owner)==null, "Loader-only observation preserves input, raw arguments and absent processing outputs.");
			write("constructor-args-after.txt", lines(passed)); write("exception.txt", caught==null ? "NONE\n" : trace(caught));
		}finally{cleanup();}
		check(cleanupCalls==1 && offered==(interrupted ? 7 : 0) && joined==offered && discover().isEmpty(), "Unconditional cleanup must offer/join exactly the observed owned handles and leave no live owned loader.");
		check(!interrupted || emptyCompleted(), "All seven interrupted-case loaders must consume identity POISON, finish successfully and retain zero per-worker load counters.");
		check(afterState!=null && afterState.equals(counters()+tablesAndScaffolds()) && !Thread.currentThread().isInterrupted(), "Fixture teardown must not change owner data or leave an interrupt pending.");
		check(calls==1 && !owner.errorState && field("ref").get(owner)==null, "Exactly one loader call and no reference input/process pipeline are allowed.");
		require(checks==12, "Expected twelve contracts, observed "+checks);
		write("workers.tsv", workers());
		write("case.tsv", "scenario\toutcome\treturned\terror\tcause\tflag_before\tflag_after\tlive_before\tcaptured_handles\tempty_queues\tpoison_offered\tjoined\tlive_after\tchecks\n"+name+"\t"+(interrupted ? "THREW" : "RETURNED")+"\t"+(returned==null ? "UNAVAILABLE" : returned)+"\t"+(caught==null ? "NONE\tNONE" : "RuntimeException\tInterruptedException")+"\t"+flagBefore+"\t"+flagAfter+"\t0\t"+observed+"\t"+observed+"\t"+offered+"\t"+joined+"\t0\t"+checks+"\n");
		final String summary="BBDUK_LOADER_INTERRUPT_BASELINE_PASS version=1 scenario="+name+" loaderCalls=1 outcome="+(interrupted ? "THREW" : "RETURNED")+" capturedHandles="+observed+" cleanupOffered="+offered+" cleanupJoined="+joined+" liveOwnedAfter=0 processCalls=0 factoryCalls=0 checks="+checks+" failures=0\n";
		write("summary.txt", summary); System.out.print(summary);
	}

	/** Always invoked, including while an observation/assertion/report write unwinds. */
	private void cleanup(){
		cleanupCalls++;
		try{
			if(owner==null){return;}
			for(final Map.Entry<Integer, Thread> entry : discover().entrySet()){
				final Thread previous=captured.put(entry.getKey(), entry.getValue()); require(previous==null || previous==entry.getValue(), "Conflicting loader handle for one partition.");
			}
			require(emptyOwned(captured), "Unexpected owned queue payload/identity during cleanup; cannot safely guess or clear it.");
			final Object poison=field("POISON").get(null);
			for(final Thread loader : captured.values()){
				require(loader.getClass()==loadType && enclosing.get(loader)==owner, "Recheck exact ownership immediately before offering termination.");
				require(queue(loader).offer(poison), "Empty owned queue refused identity POISON; abort rather than block."); offered++;
			}
			for(final Thread loader : captured.values()){loader.join(); require(!loader.isAlive(), "Exact joined loader must be terminated."); joined++;}
			require(discover().isEmpty(), "Live owned loaders remain after cleanup.");
		}catch(final Throwable failure){abort("Loader-interruption fixture cleanup failed; aborting its dedicated JVM without abandoning workers.", failure);}
	}
	/** Only live threads in this JVM with the exact private class and owner identity. */
	private TreeMap<Integer, Thread> discover() throws Exception{
		final TreeMap<Integer, Thread> result=new TreeMap<Integer, Thread>();
		for(final Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread.isAlive() && thread.getClass()==loadType && enclosing.get(thread)==owner){
				final int partition=partitionField.getInt(thread); require(partition>=0 && partition<7 && result.put(partition, thread)==null, "Owned loader partitions must be unique in [0,7).");
			}
		}
		return result;
	}
	private boolean emptyOwned(final TreeMap<Integer, Thread> threads) throws Exception{
		for(final Map.Entry<Integer, Thread> entry : threads.entrySet()){
			final Thread thread=entry.getValue(); if(thread.getClass()!=loadType || enclosing.get(thread)!=owner || partitionField.getInt(thread)!=entry.getKey() || !thread.isAlive() || queue(thread).size()!=0 || queue(thread).peek()!=null){return false;}
			if(number(thread, "addedT")!=0 || number(thread, "refReadsT")!=0 || number(thread, "refBasesT")!=0 || number(thread, "refKmersT")!=0 || member(loadType, "success").getBoolean(thread)){return false;}
		}
		return true;
	}
	private boolean emptyCompleted() throws Exception{
		if(captured.size()!=7){return false;}
		for(final Thread thread : captured.values()){
			if(thread.isAlive() || !queue(thread).isEmpty() || !member(loadType, "success").getBoolean(thread) || number(thread, "addedT")!=0 || number(thread, "refReadsT")!=0 || number(thread, "refBasesT")!=0 || number(thread, "refKmersT")!=0){return false;}
		}
		return true;
	}
	private boolean state(final long loaded, final long refReads, final boolean attempted) throws Exception{
		for(final String counter : COUNTERS){final long expected=counter.equals("refReads") || counter.equals("refKmers") ? refReads : counter.equals("refBases") ? 4*refReads : 0;if(number(counter)!=expected){return false;}}
		final AbstractKmerTable[] tables=(AbstractKmerTable[])field("keySets").get(owner); if(tables.length!=7){return false;} long sum=0;
		for(final AbstractKmerTable table : tables){if(table==null){return false;} sum+=table.size();} if(sum!=loaded){return false;}
		final List<?> names=(List<?>)field("scaffoldNames").get(owner); final IntList lengths=(IntList)field("scaffoldLengths").get(owner);
		final int size=attempted ? 2 : 1; if(names.size()!=size || lengths.size!=size || !"".equals(names.get(0)) || lengths.get(0)!=0 || (attempted && (!"1".equals(names.get(1)) || lengths.get(1)!=4))){return false;}
		if(!Arrays.equals((int[])field("refScafCounts").get(owner), new int[]{attempted ? 1 : 0})){return false;}
		for(final String key : new String[]{"scaffoldReadCounts", "scaffoldBaseCounts"}){
			final AtomicLongArray array=(AtomicLongArray)field(key).get(owner);
			if(refReads==0){if(array!=null){return false;}}else{if(array==null || array.length()!=2 || array.get(0)!=0 || array.get(1)!=0){return false;}}
		}
		return true;
	}
	private boolean exactException(final Throwable caught){
		if(caught==null || caught.getClass()!=RuntimeException.class || caught.getCause()==null || caught.getCause().getClass()!=InterruptedException.class){return false;}
		final StackTraceElement[] frames=caught.getStackTrace(); return frames.length>0 && frames[0].getClassName().equals("jgi.BBDuk") && frames[0].getMethodName().equals("spawnLoadThreads") && frames[0].getLineNumber()==1926;
	}
	private String counters() throws Exception{final StringBuilder b=new StringBuilder("field\tvalue\n");for(final String name : COUNTERS){b.append(name).append('\t').append(number(name)).append('\n');}return b.toString();}
	private String tablesAndScaffolds() throws Exception{
		long sum=0; final AbstractKmerTable[] tables=(AbstractKmerTable[])field("keySets").get(owner);for(final AbstractKmerTable table : tables){sum+=table.size();}
		return "tables\tstored\tscaffold_names\tscaffold_lengths\tref_scaffolds\tatomic_arrays\n7\t"+sum+"\t"+((List<?>)field("scaffoldNames").get(owner)).size()+"\t"+((IntList)field("scaffoldLengths").get(owner)).size+"\t"+Arrays.toString((int[])field("refScafCounts").get(owner))+"\t"+(field("scaffoldReadCounts").get(owner)==null ? "NULL" : "2x2_ZERO")+"\n";
	}
	private String workers() throws Exception{
		final StringBuilder b=new StringBuilder("partition\tthread_id\texact_class\towner_identity\tqueue_before\tpoison_offered\tjoined\tterminated\tsuccess\tadded\tref_reads\tref_bases\tref_kmers\n");
		for(final Map.Entry<Integer, Thread> entry : captured.entrySet()){final Thread t=entry.getValue();b.append(entry.getKey()).append('\t').append(t.getId()).append("\ttrue\ttrue\t0\ttrue\ttrue\ttrue\ttrue\t0\t0\t0\t0\n");}return b.toString();
	}
	@SuppressWarnings("unchecked")
	private ArrayBlockingQueue<Object> queue(final Thread thread) throws IllegalAccessException{return (ArrayBlockingQueue<Object>)queueField.get(thread);}
	private long number(final String name) throws Exception{return ((Number)field(name).get(owner)).longValue();}
	private long number(final Thread thread, final String name) throws Exception{return ((Number)member(loadType, name).get(thread)).longValue();}
	private boolean bool(final String name) throws Exception{return field(name).getBoolean(owner);}
	private static Field field(final String name) throws Exception{return member(BBDuk.class, name);}
	private static Field member(final Class<?> type, final String name) throws Exception{final Field field=type.getDeclaredField(name);field.setAccessible(true);return field;}
	private static String trace(final Throwable failure){final StringWriter text=new StringWriter();failure.printStackTrace(new PrintWriter(text));return text.toString();}
	private static String lines(final String[] values){return String.join("\n", values)+"\n";}
	private void write(final String name, final String value) throws Exception{Files.write(out.resolve(name), value.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);}
	private void check(final boolean condition, final String why){checks++;require(condition, why);}
	private static void require(final boolean condition, final String why){if(!condition){throw new IllegalArgumentException(why);}}
	private static void abort(final String why, final Throwable failure){failure.printStackTrace(System.err);KillSwitch.kill(why);}

	private final Path out, input;
	private final String name;
	private final boolean interrupted;
	private final Class<?> loadType;
	private final Field enclosing, queueField, partitionField;
	private final TreeMap<Integer, Thread> captured=new TreeMap<Integer, Thread>();
	private BBDuk owner;
	private int calls, checks, observed, cleanupCalls, offered, joined;
	private static final String INPUT="@loader-unused\nGGACGTGG\n+\nIIIIIIII\n";
	private static final String[] SETTINGS={"k=4", "mink=-1", "maskmiddle=f", "rcomp=f", "hdist=0", "edist=0", "qhdist=0", "qskip=1", "rskip=1", "speed=0", "rename=f", "minlen=1", "qtrim=f", "maxns=-1", "maq=0", "reads=1", "t=1", "simd=f", "qin=33", "array=t", "prealloc=f", "ow=f", "app=f", "json=t"};
	private static final String[] COUNTERS={"readsIn", "basesIn", "readsOut", "basesOut", "storedKmers", "refReads", "refBases", "refKmers", "readsKFiltered", "basesKFiltered", "readsKTrimmed", "basesKTrimmed", "readsQTrimmed", "basesQTrimmed", "readsFTrimmed", "basesFTrimmed", "readsQFiltered", "basesQFiltered", "readsNFiltered", "basesNFiltered", "readsEFiltered", "basesEFiltered", "readsPolyTrimmed", "basesPolyTrimmed", "readsTrimmedByOverlap", "basesTrimmedByOverlap", "readsTrimmedBySwift", "basesTrimmedBySwift", "badGcReads", "badGcBases", "badHeaderReads", "badHeaderBases"};
}
