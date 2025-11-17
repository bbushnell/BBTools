package shared;

import java.io.PrintStream;

/**
 * High-precision timer utility for measuring elapsed time in nanoseconds.
 * Provides convenient methods for starting, stopping, and reporting execution times
 * with configurable output formatting and automatic printing capabilities.
 *
 * @author Brian Bushnell
 * @date December 15, 2010
 */
public class Timer {
	
	/** Creates a Timer that outputs to System.err with tab formatting enabled. */
	public Timer(){this(System.err, true);}
	
	/** Creates a Timer and immediately prints the specified message.
	 * @param s Message to print upon creation */
	public Timer(String s){
		this(System.err, true);
		if(outstream!=null){outstream.println(s);}
	}
	
	/** Creates a Timer with custom output stream and tab formatting enabled.
	 * @param outstream_ Stream for timer output (null disables output) */
	public Timer(PrintStream outstream_){
		this(outstream_, true);
	}
	
	/**
	 * Creates a Timer with full configuration control and starts timing.
	 * @param outstream_ Stream for timer output (null disables output)
	 * @param addTab_ Whether to append tab characters to output messages
	 */
	public Timer(PrintStream outstream_, boolean addTab_){
		outstream=outstream_;
		addTab=addTab_;
		start();
	}
	
	/**
	 * Prints a message and starts timing.
	 * @param s Message to print before starting
	 * @return Start time in nanoseconds
	 */
	public long start(String s){
		if(outstream!=null){outstream.println(s);}
		return start();
	}
	
	/** Stops timing and prints the elapsed time.
	 * @return Stop time in nanoseconds */
	public long stopAndPrint(){
		long x=stop();
		if(outstream!=null){outstream.println(this);}
		return x;
	}
	
	/**
	 * Stops timing and prints message with elapsed time.
	 * @param s Message to print with timing result
	 * @return Stop time in nanoseconds
	 */
	public long stop(String s){
		long x=stop();
		if(addTab && s!=null && !s.endsWith("\t")){s=s+"\t";}
		if(outstream!=null){outstream.println(s+this);}
		return x;
	}
	
	/**
	 * Stops current timing, prints message with elapsed time, then restarts timing.
	 * @param s Message to print with timing result
	 * @return Stop time in nanoseconds
	 */
	public long stopAndStart(String s){
		long x=stop();
		if(addTab && s!=null && !s.endsWith("\t")){s=s+"\t";}
		if(outstream!=null){outstream.println(s+this);}
		start();
		return x;
	}
	
	/** Starts timing by capturing current nanosecond timestamp.
	 * @return Start time in nanoseconds */
	public long start(){
		time1=time2=System.nanoTime();
		elapsed=0;
		return time1;
	}
	
	/** Stops timing and calculates elapsed duration.
	 * @return Stop time in nanoseconds */
	public long stop(){
		time2=System.nanoTime();
		elapsed=time2-time1;
		return time2;
	}
	
	@Override
	public String toString(){
//		new Exception().printStackTrace();
		return timeInSeconds(3)+" seconds.";
	}
	
	/**
	 * Returns elapsed time as formatted string with specified decimal precision.
	 * @param decimals Number of decimal places to display
	 * @return Formatted time string
	 */
	public String timeInSeconds(int decimals) {
		return Tools.format("%."+decimals+"f", timeInSeconds());
	}
	
	/** Returns elapsed time converted from nanoseconds to seconds.
	 * @return Elapsed time in seconds as double */
	public double timeInSeconds() {
		return elapsed/1000000000d;
	}

	/** Start time in nanoseconds */
	public long time1;
	/** Stop time in nanoseconds */
	public long time2;
	/** in nanos */
	public long elapsed;
	
	/** Output stream for timer messages (null disables output) */
	public PrintStream outstream=System.err;
	/** Whether to append tab characters to output messages */
	public boolean addTab=true;
}
