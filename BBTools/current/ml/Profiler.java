package ml;

import java.util.Arrays;

public class Profiler{

	public Profiler(String prefix_, int len_){
		prefix=prefix_;
		times=new long[len_];
		reset();
	}

	void log() {
		if(!PROFILING){return;}
		log(idx);
		idx++;
	}

	void log(int idx) {
		long nanos=System.nanoTime();
		long dif=nanos-start;
		times[idx]+=dif;
		start=nanos;
	}

	public String toString() {
		StringBuilder sb=new StringBuilder();
		sb.append(prefix);
		for(int i=0; i<times.length; i++) {
//			sb.append('\t').append(String.format("%.4f", array[i]/1000000.0));
			sb.append('\t').append(String.format("%d", times[i]/1000000));
		}
		return sb.toString();
	}

	void printTimes() {
		if(PROFILING) {
			System.err.println(this);
		}
	}

	void accumulate(Profiler p) {
		for(int i=0; i<times.length; i++) {
			times[i]+=p.times[i];
		}
	}

	void reset() {
		start=System.nanoTime();
		idx=0;
	}

	void clear() {
		Arrays.fill(times, 0);
		reset();
	}

	private long start;
	private int idx=0;
	final long[] times;
	final String prefix;
	
	public static boolean PROFILING=false;//does not seem to impact speed
}