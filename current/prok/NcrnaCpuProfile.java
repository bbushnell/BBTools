package prok;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.ArrayList;
import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Optional worker CPU counters. Stores no callers, models or scavengers.
 * Enabled only by -Dbbtools.diagnostic.ncrnaCpuOut=FILE.
 * @author Ganyu */
public final class NcrnaCpuProfile {
	static final String PROPERTY="bbtools.diagnostic.ncrnaCpuOut";
	static final String OUTPUT=System.getProperty(PROPERTY);
	public static final int INIT=0, SCAVENGE=1, SHARED_SCAN=2, SSU_APPLY=3, SSU_FINISH=4, PATH=5, STAGE_COUNT=6;
	private static final Clock CLOCK=OUTPUT==null ? null : systemClock();
	public static boolean enabled(){return OUTPUT!=null;}
	public static NcrnaCpuProfile create(int families){return enabled() ? new NcrnaCpuProfile(families, CLOCK) : null;}
	static void checkOutput(){if(enabled()){require(!OUTPUT.isEmpty(), "CPU profile output path must not be empty"); RrnaResourceIO.fresh(OUTPUT, "");}}
	private static Clock systemClock(){
		final ThreadMXBean bean=ManagementFactory.getThreadMXBean();
		require(bean.isCurrentThreadCpuTimeSupported(), "Requested ncRNA CPU profiling is unsupported by this JVM");
		if(!bean.isThreadCpuTimeEnabled()){bean.setThreadCpuTimeEnabled(true);}
		require(bean.getCurrentThreadCpuTime()>=0, "CPU clock must be available before workers start");
		return new Clock(){public long cpu(){return bean.getCurrentThreadCpuTime();}};
	}
	public interface Clock {long cpu();}
	public NcrnaCpuProfile(int families_, Clock clock_){
		require(families_>=0 && clock_!=null, "CPU profile requires a family count and supported clock");
		families=families_; clock=clock_; times=new long[(families+1)*STAGE_COUNT*2]; counts=new long[families*6];
	}
	public void begin(int family, int stage){
		assert(!active && !merged && family>=0 && family<=families && stage>=0 && stage<STAGE_COUNT
			&& (family==families ? stage>=SHARED_SCAN : stage<SHARED_SCAN))
			: shared.KillSwitch.assertDie("CPU spans must be disjoint family INIT/SCAVENGE or GLOBAL scan/SSU/path stages");
		final long thread=Thread.currentThread().getId();
		assert(owner==0 || owner==thread) : shared.KillSwitch.assertDie("One CPU profile belongs to one worker thread");
		owner=thread; slot=(family*STAGE_COUNT+stage)*2; active=true;
		cpuStart=clock.cpu();
	}
	public void end(){
		final long cpu=clock.cpu();
		assert(active && owner==Thread.currentThread().getId() && cpuStart>=0 && cpu>=cpuStart)
			: shared.KillSwitch.assertDie("Current-thread CPU clocks must remain available and monotonic within a span");
		times[slot]+=cpu-cpuStart; times[slot+1]++; active=false;
	}
	public void activity(int family, long alignments, long bases, long hbmCalls, long hbmBases, long hits, long windows){
		require(!active && !merged && family>=0 && family<families, "Read family activity only after worker termination");
		require(alignments>=0 && bases>=0 && hbmCalls>=0 && hbmBases>=0 && hits>=0 && windows>=0, "Family workload counters must be nonnegative");
		final int i=family*6;
		counts[i]=alignments; counts[i+1]=bases; counts[i+2]=hbmCalls; counts[i+3]=hbmBases; counts[i+4]=hits; counts[i+5]=windows;
	}
	public void merge(NcrnaCpuProfile other){
		require(other!=null && other!=this && !active && !other.active && !other.merged && families==other.families,
			"Merge each completed worker exactly once, preserving the family registry order");
		for(int i=0; i<times.length; i++){times[i]=Math.addExact(times[i], other.times[i]);}
		for(int i=0; i<counts.length; i++){counts[i]=Math.addExact(counts[i], other.counts[i]);}
		other.merged=true;
	}
	void write(String path, ArrayList<NcrnaFamily> registry){
		require(registry.size()==families, "Family registry must not change before CPU report emission");
		final String[] names=new String[families];for(int i=0; i<families; i++){names[i]=registry.get(i).name;}
		write(path, names);
	}
	public void write(String path, String[] names){
		require(!active && !merged && names.length==families, "Only completed aggregate counters may be written");
		RrnaResourceIO.fresh(path, "");
		final ByteStreamWriter out=RrnaResourceIO.writer(path); final ByteBuilder b=new ByteBuilder();
		try{
			out.println("# Diagnostic current-thread total CPU nanoseconds, not process user time; excludes other threads (GC/JIT/I/O). Clock overhead remains included.");
			out.println("# Stages are disjoint. GLOBAL scan/SSU/path are not allocated to families. SSU_APPLY includes eligibility checks; SSU_FINISH is final-path checking.");
			out.println("# Family activity covers scavengers, excluding external SSU competition. Invocations count measured spans, including early returns.");
			out.println("family\tstage\tinvocations\tcpu_ns\talignments\taligned_bases\thbm_calls\thbm_bases\tseed_hits\twindows");
			for(int f=0; f<=families; f++){
				final String name=f==families ? "GLOBAL" : names[f];
				require(name!=null && name.indexOf('\t')<0 && name.indexOf('\n')<0, "Family name must fit one TSV cell");
				for(int stage=0; stage<STAGE_COUNT; stage++){
					if(f==families ? stage<SHARED_SCAN : stage>=SHARED_SCAN){continue;}
					final int i=(f*STAGE_COUNT+stage)*2;
					b.clear().append(name).tab().append(STAGES[stage]).tab().append(times[i+1]).tab().append(times[i]);
					for(int metric=0; metric<6; metric++){b.tab();if(stage==SCAVENGE){b.append(counts[f*6+metric]);}else{b.append("NA");}}
					out.print(b.nl());
				}
			}
		}finally{RrnaResourceIO.finish(out);}
	}
	static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}
	private static final String[] STAGES={"INIT", "SCAVENGE", "SHARED_SEED_SCAN", "SSU_APPLY", "SSU_FINISH", "FIND_PATH"};
	final int families;
	final long[] times, counts;
	private final Clock clock;
	private long owner, cpuStart;
	private int slot;
	private boolean active, merged;
}
