package ml;

import java.util.Arrays;

/**
 * Simple performance profiler for timing code segments.
 * Tracks cumulative execution times across multiple profiling points and provides
 * formatted output for performance analysis.
 * @author Brian Bushnell
 */
public class Profiler{

	/**
	 * Creates a new profiler with specified output prefix and number of timing points.
	 * @param prefix_ Text prefix for output formatting
	 * @param len_ Number of timing measurement points to track
	 */
	public Profiler(String prefix_, int len_){
		prefix=prefix_;
		times=new long[len_];
		reset();
	}

	/** Records timing for the current profiling point and advances to the next point.
	 * Only performs timing if PROFILING is enabled. */
	void log() {
		if(!PROFILING){return;}
		log(idx);
		idx++;
	}

	/**
	 * Records elapsed time since last measurement to specified timing point.
	 * Accumulates time difference from previous measurement into the specified index.
	 * @param idx Index of the timing point to record
	 */
	void log(int idx) {
		long nanos=System.nanoTime();
		long dif=nanos-start;
		times[idx]+=dif;
		start=nanos;
	}

	/**
	 * Formats all timing data as tab-delimited string.
	 * Outputs prefix followed by millisecond timings for each profiling point.
	 * @return Formatted timing data string
	 */
	public String toString() {
		StringBuilder sb=new StringBuilder();
		sb.append(prefix);
		for(int i=0; i<times.length; i++) {
//			sb.append('\t').append(String.format("%.4f", array[i]/1000000.0));
			sb.append('\t').append(String.format("%d", times[i]/1000000));
		}
		return sb.toString();
	}

	/** Prints accumulated timing data to standard error.
	 * Only outputs if PROFILING is enabled. */
	void printTimes() {
		if(PROFILING) {
			System.err.println(this);
		}
	}

	/**
	 * Adds timing data from another profiler to this profiler's accumulated times.
	 * Sums corresponding timing points across both profilers.
	 * @param p Profiler whose timing data to add to this profiler
	 */
	void accumulate(Profiler p) {
		for(int i=0; i<times.length; i++) {
			times[i]+=p.times[i];
		}
	}

	/** Resets timing state to begin new measurement cycle.
	 * Sets start time to current nanoseconds and resets index to zero. */
	void reset() {
		start=System.nanoTime();
		idx=0;
	}

	/** Clears all accumulated timing data and resets profiler state.
	 * Zeroes all timing arrays and calls reset(). */
	void clear() {
		Arrays.fill(times, 0);
		reset();
	}

	/** Nanosecond timestamp of last measurement point */
	private long start;
	/** Current profiling point index for sequential logging */
	private int idx=0;
	/** Array storing cumulative nanosecond timings for each profiling point */
	final long[] times;
	/** Text prefix used in formatted output */
	final String prefix;
	
	/** Global flag to enable or disable profiling operations */
	public static boolean PROFILING=false;//does not seem to impact speed
}