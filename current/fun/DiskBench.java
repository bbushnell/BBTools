package fun;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.QuickFile;
import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FastaReadInputStream;
import structures.ByteBuilder;

/**
 * @author Brian Bushnell
 * @date December 6, 2017
 *
 */
public class DiskBench {
	
	/**
	 * Main entry point for disk benchmarking. Creates DiskBench instance,
	 * runs benchmarking process with timing, and closes output streams.
	 * @param args Command-line arguments for benchmark configuration
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		DiskBench x=new DiskBench(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Constructs DiskBench with command-line argument parsing. Sets up
	 * benchmark parameters including path, data size, passes, threading mode,
	 * and file access method. Configures compression settings and creates
	 * output directory if needed.
	 * @param args Command-line arguments for benchmark configuration
	 */
	public DiskBench(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=new Parser();
		parser.overwrite=true;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("path")){
				path=b;
				if(path==null){path="";}
				else if(!path.endsWith("/")){path=path+"/";}
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("data") || a.equals("size")){
				data=Parse.parseKMG(b);
			}else if(a.equals("passes")){
				passes=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("gzip")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("mode")){
				assert(b!=null) : "Bad parameter: "+arg;
				if(Tools.isDigit(b.charAt(0))){mode=Integer.parseInt(b);}
				else if("read".equalsIgnoreCase(b) || "r".equalsIgnoreCase(b)){mode=READ;}
				else if("write".equalsIgnoreCase(b) || "w".equalsIgnoreCase(b)){mode=WRITE;}
				else if("readwrite".equalsIgnoreCase(b) || "rw".equalsIgnoreCase(b)){mode=READWRITE;}
				else{assert(false) : "Bad mode: "+arg;}
			}else if(a.equals("read") || a.equals("r")){
				mode=READ;
			}else if(a.equals("write") || a.equals("w")){
				mode=WRITE;
			}else if(a.equals("readwrite") || a.equals("rw")){
				mode=READWRITE;
			}else if(a.equals("printtid")){
				printTid=Parse.parseBoolean(b);
			}else if(a.equals("processbis")){
				processBis=Parse.parseBoolean(b);
			}else if(a.equals("preread")){
				preRead=Parse.parseBoolean(b);
			}
			
			else if(a.equals("method")){
				assert(b!=null) : "Bad parameter: "+arg;
				if(Tools.isDigit(b.charAt(0))){method=Integer.parseInt(b);}
				else if("BYTEFILE".equalsIgnoreCase(b) || "bf".equalsIgnoreCase(b)){method=BYTEFILE;}
				else if("TEXTFILE".equalsIgnoreCase(b) || "tf".equalsIgnoreCase(b)){method=TEXTFILE;}
				else if("QUICKFILE".equalsIgnoreCase(b) || "qf".equalsIgnoreCase(b)){method=QUICKFILE;}
				else if("BUFFEREDINPUTSTREAM".equalsIgnoreCase(b) || "bis".equalsIgnoreCase(b)){method=BUFFEREDINPUTSTREAM;}
				else if("FILEINPUTSTREAM".equalsIgnoreCase(b) || "fis".equalsIgnoreCase(b)){method=FILEINPUTSTREAM;}
				else if("BUFFEREDINPUTSTREAM2".equalsIgnoreCase(b) || "bis2".equalsIgnoreCase(b)){method=BUFFEREDINPUTSTREAM2;}
				else if("FILEINPUTSTREAM2".equalsIgnoreCase(b) || "fis2".equalsIgnoreCase(b)){method=FILEINPUTSTREAM2;}
				else{assert(false) : "Bad mode: "+arg;}
			}
			else if("BYTEFILE".equalsIgnoreCase(a) || "bf".equalsIgnoreCase(a)){method=BYTEFILE;}
			else if("TEXTFILE".equalsIgnoreCase(a) || "tf".equalsIgnoreCase(a)){method=TEXTFILE;}
			else if("QUICKFILE".equalsIgnoreCase(a) || "qf".equalsIgnoreCase(a)){method=QUICKFILE;}
			else if("BUFFEREDINPUTSTREAM".equalsIgnoreCase(a) || "bis".equalsIgnoreCase(a)){method=BUFFEREDINPUTSTREAM;}
			else if("FILEINPUTSTREAM".equalsIgnoreCase(a) || "fis".equalsIgnoreCase(a)){method=FILEINPUTSTREAM;}
			else if("BUFFEREDINPUTSTREAM2".equalsIgnoreCase(a) || "bis2".equalsIgnoreCase(a)){method=BUFFEREDINPUTSTREAM2;}
			else if("FILEINPUTSTREAM2".equalsIgnoreCase(a) || "fis2".equalsIgnoreCase(a)){method=FILEINPUTSTREAM2;}
			else if("buffer".equalsIgnoreCase(a) || "bufferlen".equalsIgnoreCase(a)){
				bufferlen=(int)Parse.parseKMGBinary(b);
			}
			
			else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			overwrite=parser.overwrite;
			threads=Shared.threads();
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		File pfile=new File(path);
		if(!pfile.exists()){pfile.mkdirs();}
	}
	
	/** Worker thread for writing random data to disk during benchmarks.
	 * Generates specified amount of data and measures write performance. */
	class WriteThread extends Thread{
		
		/**
		 * Constructs write thread for specific file and data size.
		 * @param fname_ Output filename for random data
		 * @param size_ Number of bytes to write
		 */
		public WriteThread(String fname_, long size_){
			fname=fname_;
			size=size_;
		}
		
		@Override
		public void run(){
			t=new Timer();
			written=writeRandomData(fname, size, t, overwrite);
		}
		
		/** Target filename for write operations */
		String fname;
		/** Number of bytes to write to file */
		long size;
		/** Actual bytes written during operation */
		long written=0;
		/** Timer for measuring write operation duration */
		Timer t;
		
	}
	
	/**
	 * Generates and writes random ASCII data to file using ByteStreamWriter.
	 * Creates pseudo-random text with characters in range 33-96 (printable ASCII).
	 * Uses efficient bit-shifting to generate multiple characters per random int.
	 * @param fname Output filename
	 * @param size Target number of bytes to write
	 * @param t Timer for operation measurement (may be null)
	 * @param overwrite Whether to overwrite existing files
	 * @return Actual file size written to disk
	 */
	public static long writeRandomData(final String fname, final long size, final Timer t, final boolean overwrite){
		if(t!=null){t.start();}
		long written=0;
		final Random randy=Shared.threadLocalRandom();
		FileFormat ffout=FileFormat.testOutput(fname, FileFormat.TEXT, null, true, overwrite, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();
		final ByteBuilder bb=new ByteBuilder(66000);
		final int shift=6;
		final int shiftsPerRand=32/shift;
		assert(shiftsPerRand>0);
		final long limit=size-20-shiftsPerRand*1000;
		while(written<limit){
			for(int i=0; i<1000; i+=shiftsPerRand){
				int x=randy.nextInt();
				for(int j=0; j<shiftsPerRand; j++){
					byte b=(byte)(33+x&63);
					bb.append(b);
					x>>=shift;
				}
			}
//			for(int i=0; i<1000; i+=shiftsPerRand){
//				long x=randy.nextLong();
//				for(int j=0; j<shiftsPerRand; j++){
//					byte b=(byte)(33+x&63);
//					bb.append(b);
//					x>>=shift;
//				}
//			}
			bb.nl();
			written+=bb.length;
			bsw.print(bb);
			bb.clear();
		}
		while(written<size-1){
			bb.append((byte)(33+(randy.nextInt()&63)));
			written++;
		}
		bb.nl();
		written+=bb.length;
		bsw.print(bb);
		bb.clear();
		bsw.poisonAndWait();
		File f=new File(fname);
		long diskSize=(f.length());
		if(t!=null){t.stop();}
		return diskSize;
	}
	
	/** Worker thread for reading data from disk using specified access method.
	 * Supports multiple file reading strategies for performance comparison. */
	class ReadThread extends Thread{
		
		/**
		 * Constructs read thread for specific file and thread ID.
		 * @param fname_ Input filename to read
		 * @param tid_ Thread identifier for debugging output
		 */
		public ReadThread(String fname_, int tid_){
			fname=fname_;
			tid=tid_;
		}
		
		@Override
		public void run(){
			t=new Timer();
			FileFormat ffin=FileFormat.testInput(fname, FileFormat.TEXT, null, false, false, false);
			
			if(method==BYTEFILE){
				runBf(ffin);
			}else if(method==QUICKFILE){
				runQf(ffin);
			}else if(method==TEXTFILE){
				runTf(ffin);
			}else if(method==BUFFEREDINPUTSTREAM){
				runBis(ffin, true);
			}else if(method==FILEINPUTSTREAM){
				runBis(ffin, false);
			}else if(method==BUFFEREDINPUTSTREAM2){
				runBis2(ffin, true);
			}else if(method==FILEINPUTSTREAM2){
				runBis2(ffin, false);
			}
			if(printTid){System.err.print(tid+",");}
			t.stop();
		}
		
		/**
		 * Reads file using ByteFile line-by-line access.
		 * Counts bytes read including newline characters.
		 * @param ffin Input file format configuration
		 */
		private void runBf(FileFormat ffin){
			ByteFile bf=ByteFile.makeByteFile(ffin);
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				read+=line.length+1;
			}
		}
		
		/**
		 * Reads file using QuickFile line-by-line access.
		 * Counts bytes read including newline characters.
		 * @param ffin Input file format configuration
		 */
		private void runQf(FileFormat ffin){
			QuickFile qf=new QuickFile(ffin);
			for(byte[] line=qf.nextLine(); line!=null; line=qf.nextLine()){
				read+=line.length+1;
			}
		}
		
		/**
		 * Reads file using TextFile string-based line access.
		 * Counts characters read including newline characters.
		 * @param ffin Input file format configuration
		 */
		private void runTf(FileFormat ffin){
			TextFile tf=new TextFile(ffin);
			for(String line=tf.nextLine(); line!=null; line=tf.nextLine()){
				read+=line.length()+1;
			}
		}
		
		/**
		 * Reads file using direct InputStream access with configurable buffering.
		 * Optionally processes data to find line boundaries during reading.
		 * Tests raw stream performance vs. higher-level abstractions.
		 * @param ffin Input file format configuration
		 * @param bufferedStream Whether to use BufferedInputStream vs FileInputStream
		 */
		private void runBis(FileFormat ffin, boolean bufferedStream){
			
			final byte[] buffer=new byte[bufferlen];
			InputStream is=ReadWrite.getInputStream(ffin.name(), bufferedStream, false, false);
			
			for(int r=1; r>0; ){
				r=0;
				try {
					r=is.read(buffer);
					if(r>0){read+=r;}
					
					if(processBis){
						int last=0;
						for(int i=1; i<r; i++){
							byte b=buffer[i];
							if(b=='\n'){
								cache=Arrays.copyOfRange(buffer, last, i);
								last=i+1;
							}
						}
					}
					
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			
			try {
				is.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		/**
		 * Enhanced stream reading that stores line data in ArrayList for processing.
		 * Tests memory allocation overhead during high-throughput reading.
		 * Recreates ArrayList every 800 lines to simulate batch processing.
		 * @param ffin Input file format configuration
		 * @param bufferedStream Whether to use BufferedInputStream vs FileInputStream
		 */
		private void runBis2(FileFormat ffin, boolean bufferedStream){
			
			final byte[] buffer=new byte[bufferlen];
			InputStream is=ReadWrite.getInputStream(ffin.name(), bufferedStream, false, false);
			
			list=new ArrayList<byte[]>(800);
			for(int r=1; r>0; ){
				r=0;
				try {
					r=is.read(buffer);
					if(r>0){read+=r;}
					
					if(processBis){
						int last=0;
						for(int i=1; i<r; i++){
							byte b=buffer[i];
							if(b=='\n'){
								byte[] line=Arrays.copyOfRange(buffer, last, i);
								list.add(line);
								if(list.size()>=800){
									list=new ArrayList<byte[]>(800);
								}
								last=i+1;
							}
						}
					}
					
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			
			try {
				is.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		/** Temporary storage for line data during stream processing */
		byte[] cache;
		/** Collection for storing line data during enhanced stream processing */
		ArrayList<byte[]> list;
		/** Input filename for read operations */
		String fname;
		/** Total bytes read during operation */
		long read=0;
		/** Number of lines processed during read operation */
		long lines=0;
		/** Timer for measuring read operation duration */
		Timer t;
		/** Thread identifier for debugging and output correlation */
		final int tid;
		
	}
	
	/**
	 * Generates unique filenames for benchmark pass using timestamp and
	 * random components. Creates one filename per thread to avoid conflicts.
	 * @param pass Benchmark pass number for filename uniqueness
	 * @return Array of unique filenames for thread use
	 */
	String[] makeFnames(int pass){
		String[] fnames=new String[threads];
		Random randy=new Random();
		for(int i=0; i<threads; i++){
			fnames[i]=path+pass+"_"+i+"_"+(System.nanoTime()&0xFFFF)+"_"+randy.nextInt(4096);
		}
		return fnames;
	}
	
	/**
	 * Executes concurrent read/write benchmark pass. Starts all write threads
	 * followed by read threads, then waits for completion. Measures combined
	 * read/write performance under concurrent load.
	 * @param fnamesW Filenames for write operations
	 * @param fnamesR Filenames for read operations
	 * @return Timer with elapsed time for combined operations
	 */
	Timer readWrite(String[] fnamesW, String[] fnamesR){
		Timer t=new Timer();
		
		WriteThread[] wta=new WriteThread[threads];
		long size=data/threads;
		for(int i=0; i<threads; i++){
			wta[i]=new WriteThread(fnamesW[i], size);
		}
		for(int i=0; i<threads; i++){
			wta[i].start();
		}
		
		ReadThread[] rta=new ReadThread[threads];
		for(int i=0; i<threads; i++){
			rta[i]=new ReadThread(fnamesR[i], i);
		}
		for(int i=0; i<threads; i++){
			rta[i].start();
		}
		
		for(int i=0; i<threads; i++){
			while(wta[i].getState()!=Thread.State.TERMINATED){
				try {
					wta[i].join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		for(int i=0; i<threads; i++){
			while(rta[i].getState()!=Thread.State.TERMINATED){
				try {
					rta[i].join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		t.stop();
		return t;
	}
	
	/**
	 * Executes multi-threaded write benchmark. Distributes data size evenly
	 * across threads and measures time to completion for all writers.
	 * @param fnames Output filenames for each write thread
	 * @return Timer with elapsed time for write completion
	 */
	Timer write(String[] fnames){
		Timer t=new Timer();
		WriteThread[] wta=new WriteThread[threads];
		long size=data/threads;
		for(int i=0; i<threads; i++){
			wta[i]=new WriteThread(fnames[i], size);
		}
		for(int i=0; i<threads; i++){
			wta[i].start();
		}
		for(int i=0; i<threads; i++){
			while(wta[i].getState()!=Thread.State.TERMINATED){
				try {
					wta[i].join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		t.stop();
		return t;
	}
	
	/**
	 * Executes multi-threaded read benchmark with optional pre-read warmup.
	 * If preRead enabled, performs single-threaded read first to warm disk cache.
	 * Then runs concurrent readers and measures completion time.
	 * @param fnames Input filenames for each read thread
	 * @return Timer with elapsed time for read completion
	 */
	Timer read(String[] fnames){
		
		Timer t=new Timer();
		
		if(preRead){
			ReadThread rt=new ReadThread(fnames[0], 0);
			rt.start();
			while(rt.getState()!=Thread.State.TERMINATED){
				try {
					rt.join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		
		ReadThread[] rta=new ReadThread[threads];
		for(int i=0; i<threads; i++){
			rta[i]=new ReadThread(fnames[i], i);
		}
		for(int i=0; i<threads; i++){
			rta[i].start();
		}
		for(int i=0; i<threads; i++){
			while(rta[i].getState()!=Thread.State.TERMINATED){
				try {
					rta[i].join();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			linesInternal+=(rta[i].list==null ? 0 : rta[i].list.size());
		}
		t.stop();
		return t;
	}
	
	/**
	 * Cleans up benchmark files by deleting all specified filenames.
	 * Checks for file existence before deletion to avoid errors.
	 * @param fnames Array of filenames to delete
	 */
	void delete(String[] fnames){
		for(String s : fnames){
			File f=new File(s);
			if(f.exists()){
				f.delete();
			}
		}
	}
	
	/**
	 * Main benchmark execution logic. Performs initial write, then runs
	 * specified number of passes in chosen mode (read/write/readwrite).
	 * Reports throughput in MB/s for each pass and overall performance.
	 * Cleans up files after completion.
	 * @param t0 Overall timer started at program entry
	 */
	void process(Timer t0){
		
		t0.start();
		String[] fnamesW=makeFnames(0);
		
		Timer t=write(fnamesW);
		String[] fnamesR=fnamesW;
		fnamesW=null;
		
		final long initialWriteElapsed=t.elapsed;
		
		System.err.println("Initial write:   \t"+t.toString()+"  \t"+Tools.format("%.3f MB/s", (1000.0*data)/t.elapsed));
		
		for(int pass=0; pass<passes; pass++){
			if(mode==READWRITE){
				fnamesW=makeFnames(pass);
				t=readWrite(fnamesW, fnamesR);
				delete(fnamesR);
				fnamesR=fnamesW;
				fnamesW=null;
			}else if(mode==READ){
				t=read(fnamesR);
			}else{
				delete(fnamesR);
				fnamesW=makeFnames(pass);
				t=write(fnamesW);
				fnamesR=fnamesW;
				fnamesW=null;
			}
			System.err.println("Pass        "+pass+":   \t"+t.toString()+"  \t"+Tools.format("%.3f MB/s", (1000.0*data)/t.elapsed));
		}
		delete(fnamesR);
		
		t0.stop();
		System.err.println("Overall:         \t"+t0.toString()+"  \t"+Tools.format("%.3f MB/s", (1000.0*(data*passes))/(t0.elapsed-initialWriteElapsed)));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/
	
	/** Base directory path for benchmark file creation */
	private String path="";
	
	/*--------------------------------------------------------------*/
	
	/** Buffer size for stream-based reading operations */
	private int bufferlen=4096;
	/** Total bytes of data to write per benchmark pass */
	private long data=8000000000L;
	/** Number of benchmark passes to execute */
	private int passes=2;
	
	/** Internal line count for benchmark statistics */
	public int linesInternal;
	
	/** Number of concurrent threads for benchmark operations */
	private int threads;
	
	/** Maximum number of lines to process (currently unused) */
	private long maxLines=Long.MAX_VALUE;
	
	/** Benchmark mode: READWRITE, READ, or WRITE operations */
	int mode=READWRITE;
	/** Mode constant for write-only benchmarking */
	/** Mode constant for read-only benchmarking */
	/** Mode constant for combined read/write benchmarking */
	static final int READWRITE=1, READ=2, WRITE=3;

	/** Whether to print thread IDs during execution for debugging */
	boolean printTid=false;
	/** Whether to process line boundaries during stream reading */
	boolean processBis=false;
	/** Whether to perform warmup read before main benchmark */
	boolean preRead=false;
	
	/** File access method for read operations */
	int method=BYTEFILE;
	/** Method constant for ByteFile-based reading */
	static final int BYTEFILE=1;
	/** Method constant for TextFile-based reading */
	static final int TEXTFILE=2;
	/** Method constant for BufferedInputStream reading */
	static final int BUFFEREDINPUTSTREAM=3;
	/** Method constant for direct FileInputStream reading */
	static final int FILEINPUTSTREAM=4;
	/** Method constant for enhanced BufferedInputStream with line storage */
	static final int BUFFEREDINPUTSTREAM2=5;
	/** Method constant for enhanced FileInputStream with line storage */
	static final int FILEINPUTSTREAM2=6;
	/** Method constant for QuickFile-based reading */
	static final int QUICKFILE=7;
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	/** Output stream for benchmark results and messages */
	private PrintStream outstream=System.err;
	/** Whether to enable verbose output during benchmarking */
	public static boolean verbose=false;
	/** Tracks whether any errors occurred during benchmark execution */
	public boolean errorState=false;
	/** Whether to overwrite existing benchmark files */
	private boolean overwrite=true;
	
}
