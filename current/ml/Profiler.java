package ml;

import java.util.Arrays;

/**
 * Manual profiler for caller-selected timing segments.
 * log() and printTimes() honor PROFILING; log(int) and accumulate(Profiler)
 * update caller-owned, unsynchronized instances directly.
 * @author Brian Bushnell
 */
public class Profiler{

	/**
	 * Creates a profiler with an output prefix and fixed number of timing buckets.
	 * @param prefix_ text prefix for output formatting
	 * @param len_ number of timing buckets to track
	 */
	public Profiler(String prefix_, int len_){
		prefix=prefix_;
		times=new long[len_];
		reset();
	}

	/** Records timing for the current profiling point and advances to the next point.
	 * Only performs timing if PROFILING is enabled. */
	void log(){
		if(!PROFILING){return;}
		log(idx);
		idx++;
	}

	/**
	 * Adds elapsed nanoseconds since the previous measurement to one timing bucket.
	 * @param idx valid index of the timing bucket to record
	 */
	void log(int idx){
		final long nanos=System.nanoTime();
		final long dif=nanos-start;
		times[idx]+=dif;
		start=nanos;
	}

	/**
	 * Formats all timing buckets as a tab-delimited line.
	 * @return prefix followed by each accumulated bucket in milliseconds
	 */
	public String toString(){
		final StringBuilder sb=new StringBuilder();
		sb.append(prefix);
		for(int i=0; i<times.length; i++){
			sb.append('\t').append(String.format("%d", times[i]/1000000));
		}
		return sb.toString();
	}

	/** Prints accumulated timing data to standard error.
	 * Only outputs if PROFILING is enabled. */
	void printTimes(){
		if(PROFILING){
			System.err.println(this);
		}
	}

	/**
	 * Adds another profiler's timing buckets into this profiler.
	 * @param p profiler whose timing buckets to add; must have at least this many
	 *          buckets
	 */
	void accumulate(Profiler p){
		for(int i=0; i<times.length; i++){
			times[i]+=p.times[i];
		}
	}

	/** Resets the start timestamp and sequential bucket index, preserving totals. */
	void reset(){
		start=System.nanoTime();
		idx=0;
	}

	/** Clears all accumulated timing buckets and resets profiler state. */
	void clear(){
		Arrays.fill(times, 0);
		reset();
	}

	/** Nanosecond timestamp of last measurement point */
	private long start;
	/** Current profiling point index for sequential logging */
	private int idx=0;
	/** Cumulative nanosecond timings for each bucket */
	final long[] times;
	/** Text prefix used in formatted output */
	final String prefix;
	
	/** Global flag gating log() and printTimes(). */
	public static boolean PROFILING=false;
}
