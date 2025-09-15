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
	
	public static void main(String[] args){
		double seconds=Double.parseDouble(args[0]);
		double load=Double.parseDouble(args[1]);
		launch(seconds, load);
		if(args.length>2){
			
		}
	}
	
	private KillSwitch(double seconds, double load) {
		maxSeconds=seconds;
		minLoad=load;
	}

	public static boolean launch(){
		return launch(600);
	}

	public static boolean launch(double seconds){
		return launch(seconds, 0.002);
	}
	
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
	
	public static synchronized void kill(String s){
		ballast=null;
		Exception e=new Exception(s);
		e.printStackTrace();
		kill0();
	}
	
	public static synchronized void killTraceless(String s){
		ballast=null;
		System.err.println(s);
		kill0();
	}
	
//	public static void kill(Throwable e){
//		e.printStackTrace();
//		kill0();
//	}
	
	public static synchronized void kill(){
		ballast=null;
		Exception e=new Exception("Aborting.");
		e.printStackTrace();
		kill0();
	}
	
	public static synchronized void killSilent(){
		ballast=null;
		kill0();
	}
	
	private static void kill0(){
		ballast=null;
		Runtime.getRuntime().halt(1);
	}
	
	public static void shutdown(){
		shutdownFlag.set(true);
	}
	
	public static void setKillFlag(){
		killFlag.set(true);
	}
	
	/*--------------------------------------------------------------*/
	
	public static final void throwableKill(Throwable e){
		ballast=null;
		synchronized(MemKillMessage){
			e.printStackTrace();
			kill0();
		}
	}
	
	public static final void exceptionKill(Exception e){
		ballast=null;
		synchronized(MemKillMessage){
			e.printStackTrace();
			kill0();
		}
	}
	
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
	
	public static final void assertionKill(AssertionError e){
		ballast=null;
		synchronized(MemKillMessage){
//			System.err.println(e);
			e.printStackTrace();
			kill0();
		}
	}
	
	
	/*--------------------------------------------------------------*/
	
	public static final AtomicIntegerArray allocAtomicInt(int len){
		AtomicIntegerArray ret=null;
		try {
			ret=new AtomicIntegerArray(len);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	public static final long[] allocLong1D(int len){
		long[] ret=null;
		try {
			ret=new long[len];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
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
	
	public static final int[] allocInt1D(int len){
		int[] ret=null;
		try {
			ret=new int[len];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}

	public static float[] allocFloat1D(int len) {
		float[] ret=null;
		try {
			ret=new float[len];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}

	public static double[] allocDouble1D(int len) {
		double[] ret=null;
		try {
			ret=new double[len];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	public static final byte[] allocByte1D(int len){
		byte[] ret=null;
		try {
			ret=new byte[len];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
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
	
	public static final int[][] allocInt2D(int x){
		int[][] ret=null;
		try {
			ret=new int[x][];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
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
	
	public static final float[][] allocFloat2D(int x){
		float[][] ret=null;
		try {
			ret=new float[x][];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
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
	
	public static final long[][] allocLong2D(int x){
		long[][] ret=null;
		try {
			ret=new long[x][];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
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
	
	public static int[][][] allocInt3D(int x) {
		int[][][] ret=null;
		try {
			ret=new int[x][][];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
	public static int[][][] allocInt3D(int x, int y) {
		int[][][] ret=null;
		try {
			ret=new int[x][y][];
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return ret;
	}
	
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
	
	public static long[] copyAndFill(long[] buffer, long newLength, int fillValue) {
		final long[] copy=copyOf(buffer, newLength);
		Arrays.fill(copy, buffer.length, copy.length, fillValue);
		return copy;
	}
	
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
	
	public static byte[] copyOfRange(byte[] buffer, int from, int to) {
		byte[] copy=null;
		try {
			copy=Arrays.copyOfRange(buffer, from, to);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}

	public static int[] copyOfRange(int[] buffer, int start, int limit) {
		int[] copy=null;
		try {
			copy=Arrays.copyOfRange(buffer, start, limit);
		} catch (OutOfMemoryError e) {
			memKill(e);
		}
		return copy;
	}

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
	
	private final double maxSeconds;
	private final double minLoad;
	
	private static AtomicBoolean shutdownFlag=new AtomicBoolean(false);
	private static AtomicBoolean killFlag=new AtomicBoolean(false);
	private static int count=0;
	private static KillSwitch ks;
	private static boolean suppressMessages=false;
	
	/*--------------------------------------------------------------*/
	
	private static final String MemKillMessage=new String("\nThis program ran out of memory.\n"
			+ "Try increasing the -Xmx flag and using tool-specific memory-related parameters.");

	public static void addBallast() {
		synchronized(KillSwitch.class){
			ballast=new int[20000];
			ballast2=new int[20000];
			ballast[0]=1;
			assert(ballast[0]==1);
		}
	}

	private static int[] ballast;
	private static int[] ballast2;
	
}
