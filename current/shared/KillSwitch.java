package shared;

import java.lang.management.ManagementFactory;
import java.lang.management.OperatingSystemMXBean;
import java.lang.reflect.Array;
import java.util.Arrays;
//import com.sun.management.OperatingSystemMXBean;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicIntegerArray;

/**
 * Monitors CPU utilization to determine if the program has crashed.
 * Also performs VM forced shutdowns and safe memory allocation.
 * @author Brian Bushnell
 * @date Feb 25, 2015
 *
 */
public final class KillSwitch extends Thread {
	
	/** Entry point for standalone execution that launches KillSwitch monitoring.
	 * @param args Command-line arguments: seconds (max duration) and load (minimum CPU load) */
	public static void main(String[] args){
		double seconds=Double.parseDouble(args[0]);
		double load=Double.parseDouble(args[1]);
		launch(seconds, load);
		if(args.length>2){
			
		}
	}
	
	/**
	 * Constructs a KillSwitch thread with specified monitoring parameters.
	 * @param seconds Maximum duration in seconds before shutdown if load is too low
	 * @param load Minimum system load threshold to remain active
	 */
	private KillSwitch(double seconds, double load) {
		maxSeconds=seconds;
		minLoad=load;
	}

	/** Launches KillSwitch with default timeout of 600 seconds and load threshold 0.002.
	 * @return false if already running, true if successfully started */
	public static boolean launch(){
		return launch(600);
	}

	/**
	 * Launches KillSwitch with specified timeout and default load threshold 0.002.
	 * @param seconds Maximum duration in seconds before shutdown if load is too low
	 * @return false if already running, true if successfully started
	 */
	public static boolean launch(double seconds){
		return launch(seconds, 0.002);
	}
	
	/**
	 * Launches KillSwitch thread with custom timeout and load parameters.
	 * Only one KillSwitch instance can run at a time.
	 * @param seconds Maximum duration in seconds before shutdown if load is too low
	 * @param load Minimum system load threshold to remain active
	 * @return false if already running, true if successfully started
	 */
	public static synchronized boolean launch(double seconds, double load){
		if(count>0){return false;}
		ks=new KillSwitch(seconds, load);
		ks.start();
		return true;
	}
	
	@Override
	public void run(){
		
		boolean success=monitor();
//		System.err.println("success: "+success);
		if(!success || killFlag.get()){
			if(!suppressMessages){
				System.err.println("Process has decided it has crashed, and will abort.\n" +
						"If this decision was incorrect, please re-run with the flag 'monitor=f'");
			}
			kill0();
		}
	}
	
	/**
	 * Monitors system load and returns false if load stays below threshold too long.
	 * Uses OperatingSystemMXBean to track load average every 500ms. Resets timeout
	 * whenever load exceeds minimum threshold. Returns false to trigger shutdown.
	 * @return false if timeout exceeded with low load, true for normal shutdown
	 */
	private boolean monitor(){
		
		final OperatingSystemMXBean bean=ManagementFactory.getOperatingSystemMXBean();
		if(bean.getSystemLoadAverage()<0){
			System.err.println("This OS does not support monitor, so monitoring was disabled.");
			return true;
		}
		
		final long start=System.currentTimeMillis();
		final long buffer=(long)(1+maxSeconds*1000);
		long stop=start+buffer;
//		System.err.println("start="+start+", stop="+stop+", buffer="+buffer);
//		System.err.println("shutdownFlag.get()="+shutdownFlag.get());
		while(!shutdownFlag.get()){
			try {
				sleep(500);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			final double load=bean.getSystemLoadAverage();
			final long time=System.currentTimeMillis();
			if(load>minLoad){stop=time+buffer;}
			if(time>stop){return false;}
//			System.err.println("stop-time="+(stop-time)+", load="+load);
		}
//		System.err.println("shutdownFlag.get()="+shutdownFlag.get());
		return true;
	}

	/*--------------------------------------------------------------*/
	
	/** Prints error message with stack trace and terminates the VM immediately.
	 * @param s Error message to print before termination */
	public static synchronized void kill(String s){
		ballast=null;
		Exception e=new Exception(s);
		e.printStackTrace();
		kill0();
	}
	
	/** Prints error message without stack trace and terminates the VM immediately.
	 * @param s Error message to print before termination */
	public static synchronized void killTraceless(String s){
		ballast=null;
		System.err.println(s);
		kill0();
	}
	
//	public static void kill(Throwable e){
//		e.printStackTrace();
//		kill0();
//	}
	
	/** Prints default "Aborting" message with stack trace and terminates the VM. */
	public static synchronized void kill(){
		ballast=null;
		Exception e=new Exception("Aborting.");
		e.printStackTrace();
		kill0();
	}
	
	/** Terminates the VM immediately without any output messages. */
	public static synchronized void killSilent(){
		ballast=null;
		kill0();
	}
	
	/**
	 * Internal method that releases ballast memory and halts the VM with exit code 1.
	 */
	private static void kill0(){
		ballast=null;
		Runtime.getRuntime().halt(1);
	}
	
	/** Sets shutdown flag to signal monitoring thread to stop gracefully. */
	public static void shutdown(){
		shutdownFlag.set(true);
	}
	
	/** Sets kill flag to trigger immediate termination on next monitoring cycle. */
	public static void setKillFlag(){
		killFlag.set(true);
	}
	
	/*--------------------------------------------------------------*/
	
	/** Prints Throwable stack trace and terminates VM with synchronized error handling.
	 * @param e Throwable that caused the termination */
	public static final void throwableKill(Throwable e){
		ballast=null;
		synchronized(MemKillMessage){
			e.printStackTrace();
			kill0();
		}
	}
	
	public static final void exceptionKill(Throwable e){
		ballast=null;
		synchronized(MemKillMessage){
			e.printStackTrace();
			kill0();
		}
	}
	
	/**
	 * Handles OutOfMemoryError by printing stack trace and memory guidance message.
	 * Releases ballast memory before termination to ensure error reporting succeeds.
	 * @param e OutOfMemoryError that triggered the termination
	 */
	public static final void memKill(OutOfMemoryError e){
		ballast=null;
		synchronized(MemKillMessage){
			e.printStackTrace();
			System.err.println(MemKillMessage);
//			Shared.printMemory();
//			killSilent();
			kill0();
		}
	}
	
	/** Handles AssertionError by printing stack trace and terminating the VM.
	 * @param e AssertionError that triggered the termination */
	public static final void assertionKill(AssertionError e){
		ballast=null;
		synchronized(MemKillMessage){
//			System.err.println(e);
			e.printStackTrace();
			kill0();
		}
	}
	
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Safely allocates AtomicIntegerArray with OutOfMemoryError handling.
	 * @param len Array length
	 * @return New AtomicIntegerArray or triggers termination on memory error
	 */
	public static final AtomicIntegerArray allocAtomicInt(int len){
		AtomicIntegerArray ret=null;
		try {
			ret=new AtomicIntegerArray(len);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	/**
	 * Safely allocates long array with OutOfMemoryError handling.
	 * @param len Array length
	 * @return New long array or triggers termination on memory error
	 */
	public static final long[] allocLong1D(int len){
		long[] ret=null;
		try {
			ret=new long[len];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	/**
	 * Safely allocates generic object array using reflection with memory error handling.
	 * @param <K> Array element type
	 * @param len Array length
	 * @param kc Class object for array element type
	 * @return New generic array or triggers termination on memory error
	 */
	@SuppressWarnings("unchecked")
	public static final <K> K[] allocObject1D(int len, Class<K> kc){
		K[] ret=null;
		try {
			ret=(K[]) Array.newInstance(kc, len);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	/**
	 * Safely allocates int array with OutOfMemoryError handling.
	 * @param len Array length
	 * @return New int array or triggers termination on memory error
	 */
	public static final int[] allocInt1D(int len){
		int[] ret=null;
		try {
			ret=new int[len];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}

	/**
	 * Safely allocates float array with OutOfMemoryError handling.
	 * @param len Array length
	 * @return New float array or triggers termination on memory error
	 */
	public static float[] allocFloat1D(int len) {
		float[] ret=null;
		try {
			ret=new float[len];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}

	/**
	 * Safely allocates double array with OutOfMemoryError handling.
	 * @param len Array length
	 * @return New double array or triggers termination on memory error
	 */
	public static double[] allocDouble1D(int len) {
		double[] ret=null;
		try {
			ret=new double[len];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	/**
	 * Safely allocates byte array with OutOfMemoryError handling.
	 * @param len Array length
	 * @return New byte array or triggers termination on memory error
	 */
	public static final byte[] allocByte1D(int len){
		byte[] ret=null;
		try {
			ret=new byte[len];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	/**
	 * Safely allocates char array with OutOfMemoryError handling.
	 * @param len Array length
	 * @return New char array or triggers termination on memory error
	 */
	public static final char[] allocChar1D(int len){
		char[] ret=null;
		try {
			ret=new char[len];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}

	/*--------------------------------------------------------------*/
	
	/**
	 * Safely allocates 2D int array with first dimension only, using memory error handling.
	 * @param x First dimension size
	 * @return New int[][] with null second dimension arrays
	 */
	public static final int[][] allocInt2D(int x){
		int[][] ret=null;
		try {
			ret=new int[x][];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	/**
	 * Safely allocates fully initialized 2D int array with memory error handling.
	 * Releases secondary ballast and prints dimensions on memory error.
	 * @param x First dimension size
	 * @param y Second dimension size
	 * @return New int[x][y] array or triggers termination on memory error
	 */
	public static final int[][] allocInt2D(int x, int y){
		int[][] ret=null;
		try {
			ret=new int[x][y];
		} catch (OutOfMemoryError e) {
			ballast2=null;
			System.err.print(x+","+y);
			memKill(e);
		}
		return ret;
	}
	
	/**
	 * Safely allocates 2D float array with first dimension only, using memory error handling.
	 * @param x First dimension size
	 * @return New float[][] with null second dimension arrays
	 */
	public static final float[][] allocFloat2D(int x){
		float[][] ret=null;
		try {
			ret=new float[x][];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	/**
	 * Safely allocates fully initialized 2D float array with memory error handling.
	 * Releases secondary ballast and prints dimensions on memory error.
	 * @param x First dimension size
	 * @param y Second dimension size
	 * @return New float[x][y] array or triggers termination on memory error
	 */
	public static final float[][] allocFloat2D(int x, int y){
		float[][] ret=null;
		try {
			ret=new float[x][y];
		} catch (OutOfMemoryError e) {
			ballast2=null;
			System.err.print(x+","+y);
			memKill(e);
		}
		return ret;
	}
	
	/**
	 * Safely allocates 2D long array with first dimension only, using memory error handling.
	 * @param x First dimension size
	 * @return New long[][] with null second dimension arrays
	 */
	public static final long[][] allocLong2D(int x){
		long[][] ret=null;
		try {
			ret=new long[x][];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	/**
	 * Safely allocates fully initialized 2D long array with memory error handling.
	 * @param x First dimension size
	 * @param y Second dimension size
	 * @return New long[x][y] array or triggers termination on memory error
	 */
	public static final long[][] allocLong2D(int x, int y){
		long[][] ret=null;
		try {
			ret=new long[x][y];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}

	/*--------------------------------------------------------------*/
	
	/**
	 * Safely allocates 3D int array with first dimension only, using memory error handling.
	 * @param x First dimension size
	 * @return New int[][][] with null second/third dimension arrays
	 */
	public static int[][][] allocInt3D(int x) {
		int[][][] ret=null;
		try {
			ret=new int[x][][];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	/**
	 * Safely allocates 3D int array with first two dimensions, using memory error handling.
	 * @param x First dimension size
	 * @param y Second dimension size
	 * @return New int[x][y][] with null third dimension arrays
	 */
	public static int[][][] allocInt3D(int x, int y) {
		int[][][] ret=null;
		try {
			ret=new int[x][y][];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	/**
	 * Safely allocates fully initialized 3D int array with memory error handling.
	 * @param x First dimension size
	 * @param y Second dimension size
	 * @param z Third dimension size
	 * @return New int[x][y][z] array or triggers termination on memory error
	 */
	public static int[][][] allocInt3D(int x, int y, int z) {
		int[][][] ret=null;
		try {
			ret=new int[x][y][z];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}

	/*--------------------------------------------------------------*/
	
	/**
	 * Creates copy of byte array with new length, enforcing array size limits.
	 * Terminates if requested length exceeds maximum array size.
	 * @param buffer Source array to copy
	 * @param newLength Desired length of new array
	 * @return New byte array copy or triggers termination on memory error
	 */
	public static byte[] copyOf(byte[] buffer, long newLength) {
		final int len=buffer.length;
		final int len2=(int)Tools.min(newLength, Shared.MAX_ARRAY_LEN);
		if(newLength>len2 && len2<=len){
			exceptionKill(new Exception("Tried to create an array above length limit: "+len+"," +newLength));
		}
		byte[] copy=null;
		try {
			copy=Arrays.copyOf(buffer, len2);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}
	
	/**
	 * Creates copy of float array with new length, enforcing array size limits.
	 * Terminates if requested length exceeds maximum array size.
	 * @param buffer Source array to copy
	 * @param newLength Desired length of new array
	 * @return New float array copy or triggers termination on memory error
	 */
	public static float[] copyOf(float[] buffer, long newLength) {
		final int len=buffer.length;
		final int len2=(int)Tools.min(newLength, Shared.MAX_ARRAY_LEN);
		if(newLength>len2 && len2<=len){
			exceptionKill(new Exception("Tried to create an array above length limit: "+len+"," +newLength));
		}
		float[] copy=null;
		try {
			copy=Arrays.copyOf(buffer, len2);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}
	
	/**
	 * Creates copy of double array with new length, enforcing array size limits.
	 * Terminates if requested length exceeds maximum array size.
	 * @param buffer Source array to copy
	 * @param newLength Desired length of new array
	 * @return New double array copy or triggers termination on memory error
	 */
	public static double[] copyOf(double[] buffer, long newLength) {
		final int len=buffer.length;
		final int len2=(int)Tools.min(newLength, Shared.MAX_ARRAY_LEN);
		if(newLength>len2 && len2<=len){
			exceptionKill(new Exception("Tried to create an array above length limit: "+len+"," +newLength));
		}
		double[] copy=null;
		try {
			copy=Arrays.copyOf(buffer, len2);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}
	
	/** 
	 * Copy the buffer into an array of size newLength.
	 * Fill extra cells with fillValue.
	 * @param buffer Old array
	 * @param newLength Length of new array
	 * @param fillValue Value to insert in extra cells
	 * @return New array
	 */
	public static int[] copyAndFill(int[] buffer, long newLength, int fillValue) {
		final int[] copy=copyOf(buffer, newLength);
		Arrays.fill(copy, buffer.length, copy.length, fillValue);
		return copy;
	}
	
	/**
	 * Creates copy of int array with new length, enforcing array size limits.
	 * Terminates if requested length exceeds maximum array size.
	 * @param buffer Source array to copy
	 * @param newLength Desired length of new array
	 * @return New int array copy or triggers termination on memory error
	 */
	public static int[] copyOf(int[] buffer, long newLength) {
		final int len=buffer.length;
		final int len2=(int)Tools.min(newLength, Shared.MAX_ARRAY_LEN);
		if(newLength>len2 && len2<=len){
			exceptionKill(new Exception("Tried to create an array above length limit: "+len+"," +newLength));
		}
		int[] copy=null;
		try {
			copy=Arrays.copyOf(buffer, len2);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}
	
	/**
	 * Creates copy of 2D int array with new length, enforcing array size limits.
	 * Terminates if requested length exceeds maximum array size.
	 * @param buffer Source 2D array to copy
	 * @param newLength Desired length of new array
	 * @return New int[][] copy or triggers termination on memory error
	 */
	public static int[][] copyOf(int[][] buffer, long newLength) {
		final int len=buffer.length;
		final int len2=(int)Tools.min(newLength, Shared.MAX_ARRAY_LEN);
		if(newLength>len2 && len2<=len){
			exceptionKill(new Exception("Tried to create an array above length limit: "+len+"," +newLength));
		}
		int[][] copy=null;
		try {
			copy=Arrays.copyOf(buffer, len2);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}
	
	/**
	 * Creates copy of long array with new length, enforcing array size limits.
	 * Terminates if requested length exceeds maximum array size.
	 * @param buffer Source array to copy
	 * @param newLength Desired length of new array
	 * @return New long array copy or triggers termination on memory error
	 */
	public static long[] copyOf(long[] buffer, long newLength) {
		final int len=buffer.length;
		final int len2=(int)Tools.min(newLength, Shared.MAX_ARRAY_LEN);
		if(newLength>len2 && len2<=len){
			exceptionKill(new Exception("Tried to create an array above length limit: "+len+"," +newLength));
		}
		long[] copy=null;
		try {
			copy=Arrays.copyOf(buffer, len2);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}
	
	/**
	 * Creates copy of long array with new length and fills extra cells with specified value.
	 * @param buffer Source array to copy
	 * @param newLength Desired length of new array
	 * @param fillValue Value to fill extra cells beyond original array length
	 * @return New long array copy with filled extra cells
	 */
	public static long[] copyAndFill(long[] buffer, long newLength, int fillValue) {
		final long[] copy=copyOf(buffer, newLength);
		Arrays.fill(copy, buffer.length, copy.length, fillValue);
		return copy;
	}
	
	/**
	 * Creates copy of 2D long array with new length, enforcing array size limits.
	 * Terminates if requested length exceeds maximum array size.
	 * @param buffer Source 2D array to copy
	 * @param newLength Desired length of new array
	 * @return New long[][] copy or triggers termination on memory error
	 */
	public static long[][] copyOf(long[][] buffer, long newLength) {
		final int len=buffer.length;
		final int len2=(int)Tools.min(newLength, Shared.MAX_ARRAY_LEN);
		if(newLength>len2 && len2<=len){
			exceptionKill(new Exception("Tried to create an array above length limit: "+len+"," +newLength));
		}
		long[][] copy=null;
		try {
			copy=Arrays.copyOf(buffer, len2);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}

	/*--------------------------------------------------------------*/
	
	/**
	 * Creates copy of byte array subrange with memory error handling.
	 * @param buffer Source array
	 * @param from Starting index (inclusive)
	 * @param to Ending index (exclusive)
	 * @return New byte array containing specified range
	 */
	public static byte[] copyOfRange(byte[] buffer, int from, int to) {
		byte[] copy=null;
		try {
			copy=Arrays.copyOfRange(buffer, from, to);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}

	/**
	 * Creates copy of int array subrange with memory error handling.
	 * @param buffer Source array
	 * @param start Starting index (inclusive)
	 * @param limit Ending index (exclusive)
	 * @return New int array containing specified range
	 */
	public static int[] copyOfRange(int[] buffer, int start, int limit) {
		int[] copy=null;
		try {
			copy=Arrays.copyOfRange(buffer, start, limit);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}

	/**
	 * Creates copy of long array subrange with memory error handling.
	 * @param buffer Source array
	 * @param start Starting index (inclusive)
	 * @param limit Ending index (exclusive)
	 * @return New long array containing specified range
	 */
	public static long[] copyOfRange(long[] buffer, int start, int limit) {
		long[] copy=null;
		try {
			copy=Arrays.copyOfRange(buffer, start, limit);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Maximum duration in seconds before shutdown if system load stays below threshold.
	 */
	private final double maxSeconds;
	/** Minimum system load average threshold to keep the program running. */
	private final double minLoad;
	
	/** Thread-safe flag to signal graceful shutdown of monitoring thread. */
	private static AtomicBoolean shutdownFlag=new AtomicBoolean(false);
	/**
	 * Thread-safe flag to trigger immediate termination on next monitoring cycle.
	 */
	private static AtomicBoolean killFlag=new AtomicBoolean(false);
	/** Counter to ensure only one KillSwitch instance runs at a time. */
	private static int count=0;
	/** Reference to the currently running KillSwitch instance. */
	private static KillSwitch ks;
	/** Flag to suppress crash detection messages when terminating. */
	private static boolean suppressMessages=false;
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Error message displayed when terminating due to OutOfMemoryError conditions.
	 */
	private static final String MemKillMessage=new String("\nThis program ran out of memory.\n"
			+ "Try increasing the -Xmx flag and using tool-specific memory-related parameters.");

	/** Allocates ballast memory arrays that can be released during out-of-memory conditions.
	 * Creates two 20000-element int arrays to provide emergency memory for error handling. */
	public static void addBallast() {
		synchronized(KillSwitch.class){
			ballast=new int[20000];
			ballast2=new int[20000];
			ballast[0]=1;
			assert(ballast[0]==1);
		}
	}

	/** Emergency memory that can be released during out-of-memory conditions. */
	private static int[] ballast;
	/** Secondary emergency memory for critical memory allocation failures. */
	private static int[] ballast2;
	
}
