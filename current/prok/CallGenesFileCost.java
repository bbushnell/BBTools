package prok;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.ArrayList;
import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Optional shared-load spans and post-join per-input costs; never attributes
 * overlapping process CPU to concurrent genomes. Enable with the JVM property
 * bbtools.diagnostic.callGenesFileCostOut=PREFIX (for example through JDK_JAVA_OPTIONS).
 * Writes PREFIX.startup.tsv, PREFIX.inputs.tsv and PREFIX.families.tsv only.
 * Add the existing bbtools.diagnostic.ncrnaCpuOut property to observe lazy INIT;
 * otherwise its column is explicitly NA. No clocks are read on the disabled path.
 * @author Ganyu */
public final class CallGenesFileCost implements AutoCloseable {

	static final String PROPERTY="bbtools.diagnostic.callGenesFileCostOut";
	public static CallGenesFileCost start(){
		final String prefix=System.getProperty(PROPERTY);
		return prefix==null ? null : new CallGenesFileCost(prefix, systemClock());
	}
	private static Clock systemClock(){
		final ThreadMXBean threads=ManagementFactory.getThreadMXBean();
		require(threads.isCurrentThreadCpuTimeSupported(), "File-cost diagnostics require current-thread CPU time");
		if(!threads.isThreadCpuTimeEnabled()){threads.setThreadCpuTimeEnabled(true);}
		final java.lang.management.OperatingSystemMXBean base=ManagementFactory.getOperatingSystemMXBean();
		require(base instanceof com.sun.management.OperatingSystemMXBean, "Shared-load diagnostics require process CPU time");
		final com.sun.management.OperatingSystemMXBean os=(com.sun.management.OperatingSystemMXBean)base;
		return new Clock(){
			public long wall(){return System.nanoTime();}
			public long threadCpu(){final long n=threads.getCurrentThreadCpuTime();require(n>=0, "Worker CPU clock unavailable");return n;}
			public long processCpu(){final long n=os.getProcessCpuTime();require(n>=0, "Startup process CPU clock unavailable");return n;}
		};
	}
	public interface Clock {long wall();long threadCpu();long processCpu();}
	public CallGenesFileCost(String prefix, Clock clock_){
		require(prefix!=null && !prefix.isEmpty() && clock_!=null, "File-cost output needs a fresh prefix and clock");
		RrnaResourceIO.fresh(prefix, ".startup.tsv", ".inputs.tsv", ".families.tsv");
		clock=clock_;owner=Thread.currentThread().getId();
		startup=RrnaResourceIO.writer(prefix+".startup.tsv");
		inputs=RrnaResourceIO.writer(prefix+".inputs.tsv");
		families=RrnaResourceIO.writer(prefix+".families.tsv");
		startup.println("# Selected shared-load spans only; process CPU includes JVM/helper threads. No file workers run in these spans.");
		startup.println("stage\tstatus\telapsed_ns\tprocess_cpu_ns");
		inputs.println("# Worker CPU excludes main/reader/writer/GC/JIT threads. Worker-span elapsed includes input/output waits. processReadInner elapsed is summed across contigs/workers and includes selected sequence-output waits; GFF enqueue is outside that span. Preparation before worker start and final output flush are excluded.");
		inputs.println("input_index\tinput\tworkers\treads\tbases\tworker_cpu_ns\tworker_span_elapsed_ns_including_waits\tprocess_read_inner_elapsed_ns_sum");
		families.println("# Actual post-join family counters. Lazy INIT CPU requires bbtools.diagnostic.ncrnaCpuOut; otherwise NA. INIT belongs to each GeneCaller, not the shared resource-load span. Counters include boundary-trim alignments.");
		families.println("input_index\tinput\tfamily\talignments\taligned_bases\thbm_calls\thbm_bases\tseed_hits\twindows\tlazy_init_cpu_ns");
	}
	public Span beginStartup(String stage){
		mainThread();cell(stage);require(files==null, "Shared resource spans must precede file-worker preparation");
		return new Span(stage, clock.wall(), clock.processCpu());
	}
	public void endStartup(Span span, boolean success){
		mainThread();require(span!=null && !span.done && files==null, "Complete each startup span once before dispatch");
		final long elapsed=delta(clock.wall(), span.wall), cpu=delta(clock.processCpu(), span.cpu);span.done=true;
		startup.print(buffer.clear().append(span.stage).tab().append(success?"COMPLETE":"FAILED").tab().append(elapsed).tab().append(cpu).nl());
	}
	public void skippedStartup(String stage){
		mainThread();cell(stage);require(files==null, "Disabled startup must be recorded before dispatch");
		startup.print(buffer.clear().append(stage).append("\tSKIPPED_DISABLED\t0\t0\n"));
	}
	public static final class Span{
		private Span(String s, long w, long c){stage=s;wall=w;cpu=c;}
		final String stage;final long wall, cpu;boolean done;
	}
	void prepare(ArrayList<String> paths, ArrayList<NcrnaFamily> registry){
		final String[] names=new String[registry.size()];for(int i=0; i<names.length; i++){names[i]=registry.get(i).name;}
		prepare(paths.toArray(new String[paths.size()]), names, NcrnaCpuProfile.enabled());
	}
	public void prepare(String[] paths, String[] names, boolean observeInit){
		mainThread();require(files==null && paths!=null && paths.length>0 && names!=null, "Initialize file identities once after shared resources load");
		familyNames=names.clone();for(String name:familyNames){cell(name);}
		initObserved=observeInit;files=new Input[paths.length];
		for(int i=0; i<files.length; i++){cell(paths[i]);files[i]=new Input(paths[i], names.length);}
	}
	public Worker worker(int input){
		mainThread();require(files!=null && input>=0 && input<files.length && !finished, "Worker must bind a prepared input ordinal");
		files[input].issued++;return new Worker(input, clock, this);
	}
	/** Call only on the coordinator after this worker has joined, before CPU-profile merge. */
	public void accept(Worker worker, long reads, long bases, long[] alignments, long[] alignedBases, long[] hbmCalls, long[] hbmBases,
			long[] hits, long[] windows, NcrnaCpuProfile cpu){
		mainThread();require(worker!=null && worker.parent==this && worker.stopped && worker.success && !worker.accepted && !finished,
			"Only successful joined workers may be recorded, exactly once");
		require(reads>=0 && bases>=0, "Joined input denominators cannot be negative");
		final long[][] values={alignments, alignedBases, hbmCalls, hbmBases, hits, windows};
		for(long[] row:values){require(row!=null && row.length==familyNames.length, "Family registry width must remain fixed until join");}
		final Input f=files[worker.input];worker.accepted=true;f.joined++;
		f.reads=Math.addExact(f.reads, reads);f.bases=Math.addExact(f.bases, bases);
		f.cpu=Math.addExact(f.cpu, worker.cpu);f.calling=Math.addExact(f.calling, worker.calling);
		f.first=Math.min(f.first, worker.first);f.last=Math.max(f.last, worker.last);
		for(int i=0; i<familyNames.length; i++){
			for(int j=0; j<values.length; j++){require(values[j][i]>=0, "Actual family activity cannot be negative");f.counts[i*7+j]=Math.addExact(f.counts[i*7+j], values[j][i]);}
			if(initObserved && cpu!=null){f.counts[i*7+6]=Math.addExact(f.counts[i*7+6], cpu.initializationCpu(i));}
		}
	}
	public void finish(long[] expectedAlignments){
		mainThread();require(files!=null && !finished && expectedAlignments!=null && expectedAlignments.length==familyNames.length,
			"Final file-cost report must bind the caller's completed global alignment counters");
		final long[] totals=new long[familyNames.length];
		for(Input f:files){
			require(f.issued>0 && f.issued==f.joined, "Every input's issued workers must join successfully before publication: "+f.path);
			for(int i=0; i<totals.length; i++){totals[i]=Math.addExact(totals[i], f.counts[i*7]);}
		}
		for(int i=0; i<totals.length; i++){require(totals[i]==expectedAlignments[i], "Per-input alignment counts must conserve the existing total for "+familyNames[i]);}
		for(int index=0; index<files.length; index++){
			final Input f=files[index];
			inputs.print(buffer.clear().append(index).tab().append(f.path).tab().append(f.joined).tab().append(f.reads).tab().append(f.bases)
				.tab().append(f.cpu).tab().append(delta(f.last, f.first)).tab().append(f.calling).nl());
			for(int i=0; i<familyNames.length; i++){
				buffer.clear().append(index).tab().append(f.path).tab().append(familyNames[i]);
				for(int j=0; j<6; j++){buffer.tab().append(f.counts[i*7+j]);}
				buffer.tab();if(initObserved){buffer.append(f.counts[i*7+6]);}else{buffer.append("NA");}families.print(buffer.nl());
			}
		}
		finished=true;
	}
	@Override public void close(){mainThread();closed=true;RrnaResourceIO.finish(startup, inputs, families);}
	public static final class Worker{
		private Worker(int input_, Clock clock_, CallGenesFileCost parent_){input=input_;clock=clock_;parent=parent_;}
		public void start(){require(!started, "Start one worker-cost span once");owner=Thread.currentThread().getId();first=clock.wall();cpuStart=clock.threadCpu();started=true;}
		public long beginCalling(){require(started && !stopped && owner==Thread.currentThread().getId(), "Calling time belongs to the active worker");return clock.wall();}
		public void endCalling(long begin){require(started && !stopped && owner==Thread.currentThread().getId(), "Finish calling time on its owning worker");calling=Math.addExact(calling, delta(clock.wall(), begin));}
		public void stop(boolean success_){require(started && !stopped && owner==Thread.currentThread().getId(), "Stop worker CPU on the same live thread, not after join");last=clock.wall();cpu=delta(clock.threadCpu(), cpuStart);delta(last, first);success=success_;stopped=true;}
		final int input;final Clock clock;final CallGenesFileCost parent;long owner, first, last, cpuStart, cpu, calling;boolean started, stopped, success, accepted;
	}
	private static final class Input{
		Input(String p, int families){path=p;counts=new long[families*7];}
		final String path;final long[] counts;int issued, joined;long reads, bases, cpu, calling, first=Long.MAX_VALUE, last=Long.MIN_VALUE;
	}
	private void mainThread(){require(!closed && Thread.currentThread().getId()==owner, "Only the coordinator writes file-cost records or reads joined workers");}
	private static void cell(String value){require(value!=null && !value.isEmpty() && value.indexOf('\t')<0 && value.indexOf('\n')<0 && value.indexOf('\r')<0, "Diagnostic identities must fit one TSV cell");}
	private static long delta(long end, long begin){final long n=end-begin;require(n>=0, "Diagnostic clocks must be monotonic within a measured span");return n;}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private final Clock clock;private final long owner;private final ByteStreamWriter startup, inputs, families;
	private final ByteBuilder buffer=new ByteBuilder();private Input[] files;private String[] familyNames;private boolean initObserved, finished, closed;
}
