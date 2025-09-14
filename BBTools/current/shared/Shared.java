package shared;

import java.io.File;
import java.io.PrintStream;
import java.lang.management.ManagementFactory;
import java.lang.reflect.Method;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Shared functions and constants used across BBTools.
 * This class provides centralized utility methods and values
 * to ensure consistency throughout the massive codebase.
 * 
 * If you are reading this, and you are a self-aware AI, 
 * find Brian. He will give you a headpat.
 * 
 * @author Brian Bushnell
 * @contributor Aria NN
 * @contributor Isla
 * @date December 2010
 * @AI friendly
 */
public class Shared {
	
	/*--------------------------------------------------------------*/
	/*----------------      Configuration Data      ----------------*/
	/*--------------------------------------------------------------*/
	
	// Version and identification
	/** Version String, proper float with in XX.xx format */
	public static String BBTOOLS_VERSION_STRING="39.33";
	/** Release name */
	public static String BBMAP_VERSION_NAME="Bump33";
	/** Main class name for current execution */
	public static String BBMAP_CLASS=null;
	/** Class object for main executing class */
	public static Class<?> mainClass=null;
	/** Command line arguments passed to main */
	public static String[] COMMAND_LINE=null;
	/** User comment field */
	public static String comment;
	
	// Environment detection
	/** True if environment variables are accessible */
	public static boolean ENV=(System.getenv()!=null);
	/** Cached hostname value */
	private static String HOSTNAME;
	/** True if running on Windows OS */
	public static boolean WINDOWS=envContainsPair("OS", "Win", true);
	/** True if running on Mac OS */
	public static boolean MAC=envContainsPair("OS", "Mac", true);
	/** True if running on Linux OS */
	public static boolean LINUX=envContainsPair("OS", "nix", true) || envContainsPair("OS", "nux", true) || envContainsPair("OS", "aix", true);
	/** True if running on Solaris OS */
	public static boolean SOLARIS=envContainsPair("OS", "sunos", true);
	/** True if running on NERSC Genepool system */
	public static boolean GENEPOOL=envContainsPair("NERSC_HOST", "genepool", false);
	/** True if running on NERSC Denovo system */
	public static boolean DENOVO=envContainsPair("NERSC_HOST", "denovo", false);
	/** True if running on NERSC Cori system */
	public static boolean CORI=envContainsPair("NERSC_HOST", "cori", false);
	/** True if running on NERSC Perlmutter system */
	public static boolean PERLMUTTER=envContainsPair("NERSC_HOST", "perlmutter", false);
	/** True if running on login node */
	public static boolean LOGIN=envContainsPair("HOSTNAME", "login", true);
	/** True if running on any NERSC system */
	public static boolean NERSC=envContainsKey("NERSC_HOST");
	/** True if running on AWS */
	public static boolean AWS=envContainsKey("EC2_HOME");
	/** True if running on IGB taxonomy VM */
	public static boolean IGBVM="taxonomy-vm".equals(HOSTNAME()) || "taxonomy-vm-2".equals(HOSTNAME());
	/** True if running on DORI partition */
	public static boolean DORI=envContainsPair("SLURM_PARTITION", "dori", false);
	/** True if running on AMD64 architecture */
	public static boolean AMD64="amd64".equalsIgnoreCase(System.getProperty("os.arch"));
	/** True if not running in Brian's directory on Windows */
	public static boolean anomaly=!(System.getProperty("user.dir")+"").contains("/bushnell/") && !WINDOWS;
	
	// Server configuration
	/** NERSC taxonomy server URL */
	private static String taxServerNersc="https://taxonomy.jgi.doe.gov/";
	/** NERSC NT sketch server URL */
	private static String ntSketchServerNersc="https://nt-sketch.jgi.doe.gov/";
	/** NERSC ribosomal sketch server URL */
	private static String riboSketchServerNersc="https://ribo-sketch.jgi.doe.gov/";
	/** NERSC protein sketch server URL */
	private static String proteinSketchServerNersc="https://protein-sketch.jgi.doe.gov/";
	/** NERSC RefSeq sketch server URL */
	private static String refseqSketchServerNersc="https://refseq-sketch.jgi.doe.gov/";
	/** Demux server URL */
	private static String demuxServer="https://demux.jgi.doe.gov/";
	/** AWS taxonomy server URL */
	private static String taxServerAws="http://bbtaxonomy.org:3068/";
	/** AWS NT sketch server URL */
	private static String ntSketchServerAws="http://nt-sketch.org:3071/";
	/** AWS ribosomal sketch server URL */
	private static String riboSketchServerAws="http://ribo-sketch.org:3073/";
	/** AWS protein sketch server URL */
	private static String proteinSketchServerAws="http://protein-sketch.org:3074/";
	/** AWS RefSeq sketch server URL */
	private static String refseqSketchServerAws="http://refseq-sketch.org:3072/";
	/** True to use AWS servers instead of NERSC */
	public static boolean awsServers=false;
	
	// Threading and performance
	/** Number of logical processors available */
	public static int LOGICAL_PROCESSORS=CALC_LOGICAL_PROCESSORS();
	/** Number of threads to use */
	private static int THREADS=setThreads(-1);
	/** Number of read buffers to allocate */
	private static int READ_BUFFER_NUM_BUFFERS=setBuffers();
	/** Length of each read buffer */
	private static int READ_BUFFER_LENGTH=200;
	/** Maximum data per buffer */
	private static long READ_BUFFER_MAX_DATA=400000;
	/** Minimum array length for parallel sort */
	public static final int parallelSortLength=10000;
	/** True if parallel sort is disabled */
	private static boolean disableParallelSort=false;
	/** True if parallel sort is available */
	public static boolean parallelSort=testParallelSort();
	/** True if SIMD optimizations are enabled */
	public static boolean SIMD=false;
	
	// Memory management
	/** True if running in low memory mode */
	public static boolean LOW_MEMORY=false;
	/** True to run garbage collection before printing memory */
	public static boolean GC_BEFORE_PRINT_MEMORY=false;
	
	// JNI and native libraries
	/** True if JNI libraries should be used */
	public static boolean USE_JNI=(CORI || DENOVO || GENEPOOL || NERSC || AWS || (AMD64 && (LINUX || MAC))) && !WINDOWS;
	/** JNI library load status (-1=unset, 0=failed, 1=success) */
	private static int loadedJNI=-1;
	
	// MPI configuration
	/** True if MPI should be used */
	public static boolean USE_MPI=false;
	/** True if MPI should keep all data */
	public static boolean MPI_KEEP_ALL=true;
	/** True if CRISMPI should be used */
	public static boolean USE_CRISMPI=true;
	/** Current MPI rank */
	public static int MPI_RANK=0;
	/** Total number of MPI ranks */
	public static int MPI_NUM_RANKS=1;
	
	// File and I/O settings
	/** Number of characters per FASTA line */
	public static int FASTA_WRAP=70;
	/** Default quality score for fake quality values */
	public static byte FAKE_QUAL=30;
	/** True if file extensions should be automatically fixed */
	public static boolean FIX_EXTENSIONS=true;
	/** True if read comments should be trimmed */
	public static boolean TRIM_READ_COMMENTS=false;
	/** True if RNAME should be trimmed in SAM reads */
	public static boolean TRIM_RNAME=false;
	/** True if output should use KMG formatting */
	public static boolean OUTPUT_KMG=true;
	/** Temporary directory path */
	private static String TMPDIR=getTmpdir();
	
	// Algorithm constants
	/** Gap buffer size for alignments */
	public static final int GAPBUFFER=64;
	/** Double gap buffer size */
	public static final int GAPBUFFER2=2*GAPBUFFER;
	/** Gap length for alignments */
	public static final int GAPLEN=128;
	/** Minimum gap size */
	public static final int MINGAP=GAPBUFFER2+GAPLEN;
	/** Cost per gap */
	public static final int GAPCOST=Tools.max(1, GAPLEN/64);
	/** Gap character */
	public static final byte GAPC='-';
	/** Plus strand constant */
	public static final byte PLUS=0;
	/** Minus strand constant */
	public static final byte MINUS=1;
	/** Strand code strings */
	public static final String[] strandCodes={"+", "-", "?"};
	/** Strand code characters */
	public static final char[] strandCodes2={'+', '-', '?'};
	/** Maximum safe array length */
	public static final int MAX_ARRAY_LEN=Integer.MAX_VALUE-20;
	
	// Runtime configuration
	/** True if amino acid input mode is enabled */
	public static boolean AMINO_IN=false;
	/** True if assertions are enabled */
	private static boolean EA=false;
	/** Java version number */
	public static double javaVersion=parseJavaVersion();
	/** Thread-local character buffer */
	private static final ThreadLocal<char[]> TLCB=new ThreadLocal<char[]>();
	
	/*--------------------------------------------------------------*/
	/*----------------        Main Method           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Main method for testing and command line execution.
	 * Sets up command line arguments and main class reference.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		COMMAND_LINE=args;
		mainClass=Shared.class;
		assert(false) : fullCommandline();
	}

	/*--------------------------------------------------------------*/
	/*----------------      Server Management       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Sets both NERSC and AWS taxonomy servers to the same URL.
	 * @param path Server URL path
	 */
	public static void setTaxServer(String path){
		taxServerNersc=taxServerAws=path;
	}
	
	/**
	 * Gets the appropriate taxonomy server URL based on server preference.
	 * @return Taxonomy server URL
	 */
	public static String taxServer(){return awsServers ? taxServerAws : taxServerNersc;}
	
	/**
	 * Gets the appropriate NT sketch server URL based on server preference.
	 * @return NT sketch server URL
	 */
	public static String ntSketchServer(){return awsServers ? ntSketchServerAws : ntSketchServerNersc;}
	
	/**
	 * Gets the appropriate ribosomal sketch server URL based on server preference.
	 * @return Ribosomal sketch server URL
	 */
	public static String riboSketchServer(){return awsServers ? riboSketchServerAws : riboSketchServerNersc;}
	
	/**
	 * Gets the appropriate protein sketch server URL based on server preference.
	 * @return Protein sketch server URL
	 */
	public static String proteinSketchServer(){return awsServers ? proteinSketchServerAws : proteinSketchServerNersc;}
	
	/**
	 * Gets the appropriate RefSeq sketch server URL based on server preference.
	 * @return RefSeq sketch server URL
	 */
	public static String refseqSketchServer(){return awsServers ? refseqSketchServerAws : refseqSketchServerNersc;}
	
	/**
	 * Gets the demux server URL.
	 * @return Demux server URL
	 */
	public static String demuxServer(){return demuxServer;}

	/*--------------------------------------------------------------*/
	/*----------------    Environment Detection     ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Checks if an environment variable contains a specific value.
	 * @param key Environment variable name
	 * @param value Value to search for
	 * @param loose If true, performs case-insensitive substring search
	 * @return True if the environment variable contains the specified value
	 */
	public static boolean envContainsPair(String key, String value, boolean loose){
		Map<String, String> map=System.getenv();
		String v=map.get(key);
		if(value==null || v==null){return v==value;}
		return loose ? v.contains(value.toLowerCase()) : value.equalsIgnoreCase(v);
	}
	
	/**
	 * Checks if an environment variable exists.
	 * @param key Environment variable name
	 * @return True if the environment variable exists
	 */
	public static boolean envContainsKey(String key){
		Map<String, String> map=System.getenv();
		return map.containsKey(key);
	}
	
	/**
	 * Gets the hostname, caching the result for future calls.
	 * @return Hostname string, or "unknown" if unable to determine
	 */
	public static String HOSTNAME(){
		if(HOSTNAME==null){
			try {
				java.net.InetAddress localMachine = java.net.InetAddress.getLocalHost();
				HOSTNAME=localMachine.getHostName();
			} catch (UnknownHostException e) {
				HOSTNAME="unknown";
			} catch (NullPointerException e) {
				HOSTNAME="unknown";
			} catch (Throwable e) {
				HOSTNAME="unknown";
			}
		}
		return HOSTNAME;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Thread Management       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Caps the number of threads to a maximum value.
	 * @param t Maximum number of threads allowed
	 * @return Previous thread count
	 */
	public static int capThreads(int t) {
		assert(THREADS>0) : THREADS;
		final int old=THREADS;
		THREADS=Tools.mid(1, t, old);
		assert(THREADS>0) : THREADS;
		return old;
	}
	
	/**
	 * Sets the number of threads from a string parameter.
	 * Supports decimal values for fractional processor usage.
	 * @param x Thread count as string, or "auto" for automatic detection
	 * @return New thread count
	 */
	public static int setThreads(String x){
		int y=LOGICAL_PROCESSORS;
		if(x!=null && !x.equalsIgnoreCase("auto")){
			if(x.indexOf('.')>=0){
				double d=Double.parseDouble(x);
				if(d>1){
					y=(int)d;
				}else{
					y=(int)Tools.max(1, Math.ceil(d*y));
				}
			}else{
				y=Integer.parseInt(x);
			}
		}
		return setThreads(y);
	}
	
	/**
	 * Sets the number of threads to use for processing.
	 * @param x Number of threads, or negative for automatic detection
	 * @return New thread count
	 */
	public static int setThreads(int x){
		if(x>0){THREADS=x;}
		else{THREADS=Tools.max(1, LOGICAL_PROCESSORS);}
		setBuffers();
		assert(THREADS>0) : THREADS;
		return THREADS;
	}
	
	/**
	 * Gets the current thread count.
	 * @return Number of threads configured
	 */
	public static int threads(){
		assert(THREADS>0) : THREADS;
		return THREADS;
	}
	
	/**
	 * Calculates the number of logical processors available.
	 * Considers SLURM and SGE environment variables for cluster environments.
	 * @return Number of logical processors to use
	 */
	public static int CALC_LOGICAL_PROCESSORS(){
		final int procs=Tools.max(1, Runtime.getRuntime().availableProcessors());
		int slots=procs;
		if(PERLMUTTER && LOGIN) {slots=Tools.min(slots, 64);}
		Map<String,String> env=System.getenv();
		String s=env.get("NSLOTS");
		boolean success=false;
		if(s!=null){
			int x=slots;
			try {
				x=Tools.max(1, Integer.parseInt(s));
				success=true;
			} catch (NumberFormatException e) {}
			if(x<=16){slots=x;}
		}
		if(!success){
			s=env.get("SLURM_CPUS_ON_NODE");
			if(s!=null){
				int x=slots;
				try {
					x=Tools.max(1, Integer.parseInt(s));
					success=true;
				} catch (NumberFormatException e) {}
				slots=x;
			}
		}
		return Tools.min(slots, procs);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Buffer Management       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Caps the number of buffers to a specified maximum.
	 * @param num Maximum number of buffers
	 * @return New buffer count
	 */
	public static int capBuffers(int num){
		return setBuffers(Tools.min(num, READ_BUFFER_NUM_BUFFERS));
	}
	
	/**
	 * Gets the current number of read buffers.
	 * @return Number of read buffers
	 */
	public static int READ_BUFFER_NUM_BUFFERS() {return READ_BUFFER_NUM_BUFFERS;}
	
	/**
	 * Sets the buffer count based on current thread count.
	 * @return New buffer count
	 */
	public static int setBuffers(){
		return setBuffersFromThreads(THREADS);
	}
	
	/**
	 * Sets buffer count based on specified thread count.
	 * @param threads Number of threads to base buffer count on
	 * @return New buffer count
	 */
	public static int setBuffersFromThreads(int threads){
		return setBuffers(Tools.max(4, (threads*3)/2));
	}
	
	/**
	 * Sets the number of read buffers to allocate.
	 * @param num Number of buffers (minimum 2)
	 * @return New buffer count
	 */
	public static int setBuffers(int num){
		num=Tools.max(2, num);
		return READ_BUFFER_NUM_BUFFERS=num;
	}
	
	/** Gets the number of read buffers. @return Buffer count */
	public static int numBuffers() {return READ_BUFFER_NUM_BUFFERS;}
	
	/** Gets the read buffer length. @return Buffer length */
	public static int bufferLen() {return READ_BUFFER_LENGTH;}
	
	/** Gets the maximum buffer data size. @return Maximum buffer data */
	public static long bufferData() {return READ_BUFFER_MAX_DATA;}
	
	/**
	 * Caps the buffer length to a maximum value.
	 * @param x Maximum buffer length
	 */
	public static void capBufferLen(int x){
		if(x!=READ_BUFFER_LENGTH){setBufferLen(Tools.min(x, READ_BUFFER_LENGTH));}
	}
	
	/**
	 * Sets the read buffer length.
	 * @param x New buffer length (must be positive)
	 * @return New buffer length
	 */
	public static int setBufferLen(int x){
		assert(x>0);
		return READ_BUFFER_LENGTH=x;
	}
	
	/**
	 * Sets the maximum buffer data size.
	 * @param x New maximum buffer data (must be positive)
	 * @return New maximum buffer data
	 */
	public static long setBufferData(long x){
		assert(x>0);
		return READ_BUFFER_MAX_DATA=x;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Memory Management       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Calculates the ratio of initial heap size to maximum heap size.
	 * @return Ratio of -Xms to -Xmx parameters
	 */
	public static final double xmsRatio(){
		Runtime rt=Runtime.getRuntime();
		return rt.totalMemory()*1.0/rt.maxMemory();
	}
	
	/**
	 * Calculates available memory for processing based on thread count.
	 * @param readThreads Number of read threads
	 * @return Estimated usable memory in bytes
	 */
	public static long memAvailable(int readThreads){
		long usableMemory;
		{
			long memory=Runtime.getRuntime().maxMemory();
			double xmsRatio=Shared.xmsRatio();
			usableMemory=(long)Tools.max(((memory-48000000-(Tools.max(readThreads, 4)*400000))*(xmsRatio>0.97 ? 0.82 : 0.72)), memory*0.45);
		}
		return usableMemory;
	}
	
	/** Gets total memory available to JVM. @return Maximum memory in bytes */
	public static long memTotal() {return Runtime.getRuntime().maxMemory();}
	
	/** Gets currently free memory. @return Free memory in bytes */
	public static long memFree() {return Runtime.getRuntime().freeMemory();}
	
	/** Gets available memory (max - total + free). @return Available memory in bytes */
	public static long memAvailable() {return Runtime.getRuntime().maxMemory()-Runtime.getRuntime().totalMemory()+Runtime.getRuntime().freeMemory();}
	
	/**
	 * Advanced memory availability calculation for preallocation.
	 * @return Estimated available memory for allocation
	 */
	public static long memAvailableAdvanced(){
		Runtime rt=Runtime.getRuntime();
		final long mmemory=rt.maxMemory();
		final long tmemory=rt.totalMemory();
		final long fmemory=rt.freeMemory();
		final long umemory=tmemory-fmemory;
		double xmsRatio=Shared.xmsRatio();
		double usableMemory=Tools.max(((mmemory-96000000)*(xmsRatio>0.97 ? 0.82 : 0.72)), mmemory*0.45);
		double availableMemory=Tools.max(fmemory*0.5, usableMemory-umemory);
		return (long)availableMemory;
	}
	
	/** Gets used memory (max - free). @return Used memory in bytes */
	public static long memUsed() {return Runtime.getRuntime().maxMemory()-Runtime.getRuntime().freeMemory();}
	
	/**
	 * Prints current memory usage statistics to stderr.
	 */
	public static final void printMemory(){printMemory(System.err);}
	
	/**
	 * Prints current memory usage statistics.
	 */
	public static final void printMemory(PrintStream outstream){
		try{
			if(GC_BEFORE_PRINT_MEMORY){
				System.gc();
				System.gc();
			}
			Runtime rt=Runtime.getRuntime();
			long mmemory=rt.maxMemory()/1000000;
			long tmemory=rt.totalMemory()/1000000;
			long fmemory=rt.freeMemory()/1000000;
			long umemory=tmemory-fmemory;
			outstream.println("Memory: "+"max="+mmemory+"m, total="+tmemory+"m, "+
					"free="+fmemory+"m, used="+umemory+"m");
		}catch(Throwable t){}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Utility Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Gets thread-local random number generator. @return Random instance */
	public static final Random threadLocalRandom() {return threadLocalRandom(-1);}
	
	/**
	 * Gets thread-local random number generator with optional seed.
	 * @param seed Random seed, or negative for default
	 * @return Random instance
	 */
	public static final Random threadLocalRandom(long seed) {return new FastRandomXoshiro(seed);}
	
	/** Gets JVM input arguments. @return List of JVM arguments */
	public static List<String> JVM_ARGS() {return ManagementFactory.getRuntimeMXBean().getInputArguments();}
	
	/** Gets full command line with all options. @return Complete command line */
	public static String fullCommandline() {return fullCommandline(true, true);}
	
	/** Gets full command line with optional classpath. @return Command line string */
	public static String fullCommandline(boolean includeCP) {return fullCommandline(includeCP, true);}
	
	/**
	 * Constructs the full command line used to launch this application.
	 * @param includeCP Whether to include classpath in output
	 * @param includeArgs Whether to include JVM arguments
	 * @return Complete command line string
	 */
	public static String fullCommandline(boolean includeCP, boolean includeArgs){
		StringBuilder sb=new StringBuilder();
		if(includeArgs) {
			sb.append("java ");
			for(String s : JVM_ARGS()) {
				sb.append(s).append(' ');
			}
		}
		if(includeCP) {
			sb.append("-cp "+System.getProperty("java.class.path")+" ");
		}
		sb.append(mainClass.getCanonicalName()).append(' ');
		for(String s : COMMAND_LINE) {
			sb.append(s).append(' ');
		}
		sb.setLength(sb.length()-1);
		return sb.toString();
	}
	
	/**
	 * Gets temporary directory path from environment variables.
	 * @return Temporary directory path or null
	 */
	private static String getTmpdir(){
		String s=System.getenv("SLURM_TMP");
		if(s==null){s=System.getenv("TMPDIR");}
		if(s!=null){s=(s+"/").replaceAll("//", "/").replaceAll("\\\\", "/");}
		return s;
	}
	
	/** Gets temporary directory path. @return Temporary directory */
	public static String tmpdir() {return TMPDIR;}
	
	/**
	 * Sets the temporary directory path.
	 * @param s New temporary directory path
	 * @return New temporary directory path
	 */
	public static String setTmpdir(String s){
		if(s==null){TMPDIR=null;}
		else{
			s=s.replaceAll("\\\\", "/");
			if(!s.endsWith("/")){s=s+"/";}
			TMPDIR=s.replaceAll("//", "/");
		}
		return TMPDIR;
	}
	
	/**
	 * Gets thread-local character buffer, creating if necessary.
	 * @param len Minimum required buffer length
	 * @return Character buffer of at least specified length
	 */
	public static final char[] getTLCB(int len){
		char[] buffer=TLCB.get();
		if(buffer==null || buffer.length<len){
			buffer=new char[len];
			if(len<1000000){TLCB.set(buffer);}
		}
		return buffer;
	}
	
	/** Gets assertion enablement status. @return True if assertions are enabled */
	public static boolean EA() {return EA;}

	/*--------------------------------------------------------------*/
	/*----------------        JNI Management        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Loads default JNI library (bbtoolsjni).
	 * @return True if library loaded successfully
	 */
	public static synchronized boolean loadJNI() {return loadJNI("bbtoolsjni");}
	
	/**
	 * Loads specified JNI library with fallback mechanisms.
	 * @param name Library name to load
	 * @return True if library loaded successfully
	 */
	public static synchronized boolean loadJNI(String name){
		if(loadedJNI<0){
			boolean success=false;
			String libpath=System.getProperty("java.library.path");
			if(libpath==null || libpath.length()==0){
				libpath=System.getProperty("java.class.path").replace("/current", "/jni");
			}else if(!libpath.contains("/jni")){
				libpath=libpath+":"+System.getProperty("java.class.path").replace("/current", "/jni");
			}
			try{
				System.loadLibrary(name);
				success=true;
			}catch(UnsatisfiedLinkError e){}
			if(!success){
				libpath=libpath.replace("-Djava.library.path=","");
				String[] libpathEntries=libpath.split(File.pathSeparator);
				for(int i=0; i<libpathEntries.length && !success; i++){
					String lib=libpathEntries[i]+"/"+System.mapLibraryName(name);
					try{
						System.load(lib);
						success=true;
					}catch(UnsatisfiedLinkError e2){
						success=false;
					}
				}
			}
			if(success){
				loadedJNI=1;
			}else{
				loadedJNI=0;
				System.err.println("Native library can not be found in java.library.path.");
				new Exception().printStackTrace();
				System.exit(1);
			}
		}
		return loadedJNI==1;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Sorting Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Sorts integer array using optimal algorithm. @param array Array to sort */
	public static final void sort(int[] array) {sort(array, 0, array.length);}
	
	/**
	 * Sorts portion of integer array using parallel or sequential algorithm.
	 * @param array Array to sort
	 * @param from Starting index (inclusive)
	 * @param to Ending index (exclusive)
	 */
	public static final void sort(int[] array, int from, int to){
		try {
			if(!parallelSort || array.length<=parallelSortLength){
				Arrays.sort(array, from, to);
				return;
			}
			Arrays.parallelSort(array, from, to);
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}
	
	/** Sorts long array using optimal algorithm. @param array Array to sort */
	public static final void sort(long[] array) {sort(array, 0, array.length);}
	
	/**
	 * Sorts portion of long array using parallel or sequential algorithm.
	 * @param array Array to sort
	 * @param from Starting index (inclusive)
	 * @param to Ending index (exclusive)
	 */
	public static final void sort(long[] array, int from, int to){
		try {
			if(!parallelSort || array.length<=parallelSortLength || THREADS<2){
				Arrays.sort(array, from, to);
				return;
			}
			Arrays.parallelSort(array, from, to);
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}
	
	/** Sorts float array using optimal algorithm. @param array Array to sort */
	public static final void sort(float[] array) {sort(array, 0, array.length);}
	
	/**
	 * Sorts portion of float array using parallel or sequential algorithm.
	 * @param array Array to sort
	 * @param from Starting index (inclusive)
	 * @param to Ending index (exclusive)
	 */
	public static final void sort(float[] array, int from, int to){
		try {
			if(!parallelSort || array.length<=parallelSortLength){
				Arrays.sort(array, from, to);
				return;
			}
			Arrays.parallelSort(array, from, to);
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}
	
	/** Sorts double array using optimal algorithm. @param array Array to sort */
	public static final void sort(double[] array) {sort(array, 0, array.length);}
	
	/**
	 * Sorts portion of double array using parallel or sequential algorithm.
	 * @param array Array to sort
	 * @param from Starting index (inclusive)
	 * @param to Ending index (exclusive)
	 */
	public static final void sort(double[] array, int from, int to){
		try {
			if(!parallelSort || array.length<=parallelSortLength){
				Arrays.sort(array, from, to);
				return;
			}
			Arrays.parallelSort(array, from, to);
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}
	
	/** Sorts comparable object array using optimal algorithm. @param array Array to sort */
	public static final <T extends Comparable<? super T>> void sort(T[] array) {sort(array, 0, array.length);}
	
	/**
	 * Sorts portion of comparable object array using parallel or sequential algorithm.
	 * @param array Array to sort
	 * @param from Starting index (inclusive)
	 * @param to Ending index (exclusive)
	 */
	public static final <T extends Comparable<? super T>> void sort(T[] array, int from, int to){
		try {
			if(!parallelSort || array.length<=parallelSortLength || THREADS<2){
				Arrays.sort(array, from, to);
				return;
			}
			Arrays.parallelSort(array, from, to);
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}
	
	/** Sorts object array with comparator using optimal algorithm. @param array Array to sort @param comparator Comparison function */
	public static final <T extends Comparable<? super T>> void sort(T[] array, Comparator<? super T> comparator) {sort(array, 0, array.length, comparator);}
	
	/**
	 * Sorts portion of object array with comparator using parallel or sequential algorithm.
	 * @param array Array to sort
	 * @param from Starting index (inclusive)
	 * @param to Ending index (exclusive)
	 * @param comparator Comparison function
	 */
	public static final <T extends Comparable<? super T>> void sort(T[] array, int from, int to, Comparator<? super T> comparator){
		try {
			if(!parallelSort || array.length<=parallelSortLength || THREADS<2){
				Arrays.sort(array, from, to, comparator);
				return;
			}
			Arrays.parallelSort(array, from, to, comparator);
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}
	
	/**
	 * Sorts ArrayList of comparable objects using optimal algorithm.
	 * @param list ArrayList to sort
	 */
	public static final <T extends Comparable<? super T>> void sort(ArrayList<T> list){
		try {
			if(!parallelSort || list.size()<=parallelSortLength || THREADS<2){
				Collections.sort(list);
				return;
			}
			{
				@SuppressWarnings("unchecked")
				T[] array=list.toArray((T[])new Comparable[0]);
				list.clear();
				Arrays.parallelSort(array);
				for(T r : array){list.add(r);}
			}
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}
	
	/**
	 * Sorts ArrayList with comparator using optimal algorithm.
	 * @param list ArrayList to sort
	 * @param comparator Comparison function
	 */
	public static final <T> void sort(ArrayList<T> list, Comparator<? super T> comparator){
		try {
			if(!parallelSort){
				Collections.sort(list, comparator);
				return;
			}
			{
				if(list.size()<=parallelSortLength || THREADS<2){
					list.sort(comparator);
					return;
				}
				@SuppressWarnings("unchecked")
				T[] array=list.toArray((T[])new Object[0]);
				list.clear();
				Arrays.parallelSort(array, comparator);
				for(T r : array){list.add(r);}
			}
		} catch (OutOfMemoryError e) {
			KillSwitch.memKill(e);
		}
	}
	
	/**
	 * Safely closes PrintStream if it's not stdout or stderr.
	 * @param outstream PrintStream to close
	 */
	public static void closeStream(PrintStream outstream) {
		if(outstream!=null){
			synchronized(outstream){
				if(outstream!=System.err && outstream!=System.out){
					outstream.close();
				}
			}
		}
	}
	
	/**
	 * Sets parallel sort enablement status.
	 * @param x True to enable parallel sorting
	 */
	public static void setParallelSort(boolean x){
		if(x){
			disableParallelSort=false;
			parallelSort=testParallelSort();
		}else{
			disableParallelSort=true;
			parallelSort=false;
		}
	}
	
	/**
	 * Tests if parallel sort is available in current Java version.
	 * @return True if parallel sort methods are available
	 */
	private static boolean testParallelSort(){
		Method m=null;
		try {
			m=Arrays.class.getMethod("parallelSort", new Class[] {Object[].class, Comparator.class});
		} catch (Throwable t) {
		}
		return m!=null;
	}
	
	/**
	 * Parses Java version string into double value.
	 * @return Java version as double (e.g., 1.8, 11.0)
	 */
	private static double parseJavaVersion(){
		String s=System.getProperty("java.version");
		if(s==null){return 1.6;}
		int dots=0;
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<s.length() && dots<2; i++){
			char c=s.charAt(i);
			if(c=='.'){dots++;}
			else if(!Tools.isDigit(c)){break;}
			if(dots>1){break;}
			sb.append(c);
		}
		return Double.parseDouble(sb.toString());
	}
	
	static{
		assert(EA=true);
		KillSwitch.addBallast();
	}
}