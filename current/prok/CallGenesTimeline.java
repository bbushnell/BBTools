package prok;

import java.io.IOException;
import java.io.OutputStream;
import java.lang.management.ManagementFactory;
import java.lang.management.RuntimeMXBean;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import structures.ByteBuilder;

/** Optional single-input wall-clock milestones, with actual processRead occupancy.
 * This observes scheduling without modifying the caller or batch sizes.
 * @author Ganyu */
public final class CallGenesTimeline {
	static final String PROPERTY="bbtools.diagnostic.callGenesTimelineOut";
	static final String HEADER="event\tepoch_ms\tprocess_cpu_ns\tactive_workers\tpeak_workers\tcompleted_contigs\tpid";

	public static CallGenesTimeline start(){
		final String path=System.getProperty(PROPERTY);
		if(path==null){return null;}
		require(!path.isEmpty(), "Requested timeline requires a fresh output path");
		final RuntimeMXBean runtime=ManagementFactory.getRuntimeMXBean();
		final String name=runtime.getName();
		final int at=name.indexOf('@');
		require(at>0, "Java 8 runtime name must expose the JVM PID for the external sampler");
		final long pid=Long.parseLong(name.substring(0, at));
		final java.lang.management.OperatingSystemMXBean base=ManagementFactory.getOperatingSystemMXBean();
		require(base instanceof com.sun.management.OperatingSystemMXBean, "Process CPU clock unavailable on this JVM");
		final com.sun.management.OperatingSystemMXBean os=(com.sun.management.OperatingSystemMXBean)base;
		final Clock clock=new Clock(){
			public long epoch(){return System.currentTimeMillis();}
			public long cpu(){final long ns=os.getProcessCpuTime();require(ns>=0, "Process CPU clock became unavailable");return ns;}
		};
		final CallGenesTimeline result=new CallGenesTimeline(path, pid, runtime.getStartTime(), clock);
		Runtime.getRuntime().addShutdownHook(new Thread("CallGenes-timeline-exit"){
			@Override public void run(){
				try{result.shutdown();}catch(Throwable e){e.printStackTrace(System.err);Runtime.getRuntime().halt(1);}
			}
		});
		return result;
	}

	public interface Clock {long epoch();long cpu();}
	public CallGenesTimeline(String path, long pid_, long jvmStart, Clock clock_){
		require(pid_>0 && jvmStart>0 && clock_!=null, "Timeline requires a real PID, JVM start and clocks");
		pid=pid_;clock=clock_;
		RrnaResourceIO.fresh(path, "", ".pid");
		try{
			// Synchronous writes are intentional for these few milestones and prompt PID discovery.
			// No background writer is introduced whose own shutdown could race EXIT_HOOK.
			out=Files.newOutputStream(Paths.get(path), StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE);
			Files.write(Paths.get(path+".pid"), (pid+"\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
			out.write((HEADER+"\n").getBytes(StandardCharsets.US_ASCII));
		}catch(IOException e){throw new IllegalStateException("Cannot create timeline/PID files: "+path, e);}
		row("JVM_START", jvmStart, -1, 0, 0, 0);
		event("MAIN_START");
	}
	public void singleInput(int count){require(count==1, "Timeline diagnostic is for one concatenated input, not repeated per-file phases");}
	public synchronized void workers(int count){
		require(count>0 && target==0 && active==0 && !first, "Set the actual worker target once, before dispatch");
		target=count;
	}
	public synchronized void enter(){
		assert(target>0 && active<target && !finished) : shared.KillSwitch.assertDie("Timeline occupancy brackets actual processRead entries, bounded by launched workers");
		active++;peak=Math.max(peak, active);
		if(!first){first=true;event("FIRST_CONTIG_DISPATCHED");}
		if(!full && active==target){full=true;event("WORKERS_FULL_OCCUPANCY");}
	}
	public synchronized void leave(boolean success){
		assert(active>0 && !finished) : shared.KillSwitch.assertDie("Every processRead exit must balance exactly one timeline entry");
		active--;
		if(success){completed++;lastEpoch=clock.epoch();lastActive=active;}
	}
	public synchronized void finishContigs(){
		require(active==0 && !finished, "Join every worker before publishing the last actual contig completion");
		finished=true;
		// The wall snapshot was taken on completion; no historical process CPU was measured.
		if(completed>0){row("LAST_CONTIG_DONE", lastEpoch, -1, lastActive, peak, completed);}
	}
	public synchronized void event(String name){row(name, clock.epoch(), clock.cpu(), active, peak, completed);}
	public synchronized void shutdown(){
		require(!closed, "Timeline shutdown must run exactly once");
		event("EXIT_HOOK");
		try{out.close();closed=true;}catch(IOException e){throw new IllegalStateException("Timeline close failed", e);}
	}
	private void row(String name, long epoch, long cpu, int workers, int peakWorkers, long done){
		require(!closed && name.indexOf('\t')<0 && name.indexOf('\n')<0 && epoch>0 && cpu>=-1,
			"Milestone rows must retain literal names, valid time and explicit missing CPU");
		final ByteBuilder b=new ByteBuilder(160);
		b.append(name).tab().append(epoch).tab();
		if(cpu<0){b.append("NA");}else{b.append(cpu);}
		b.tab().append(workers).tab().append(peakWorkers).tab().append(done).tab().append(pid).nl();
		try{out.write(b.array, 0, b.length());out.flush();}
		catch(IOException e){
			// FIRST/FULL may be emitted by a worker; a swallowed failure would hang its consumers.
			shared.KillSwitch.kill("Timeline milestone write failed: "+name+"; "+e);
			throw new IllegalStateException(e);
		}
	}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private final Clock clock;
	private final long pid;
	private final OutputStream out;
	private int active,peak,target,lastActive;
	private long completed,lastEpoch;
	private boolean first,full,finished,closed;
}
