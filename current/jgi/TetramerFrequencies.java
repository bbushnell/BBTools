/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jgi;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;

import bin.BinObject;
import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import shared.Vector;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;

/**
 *
 * @author syao
 * @contributor Brian Bushnell
 * Last updated : 01102018
 */
public class TetramerFrequencies {
	/** Program entry point for tetramer frequency analysis.
	 * @param args Command-line arguments */
	public static void main(String[] args){

		System.out.println("Start Tetramer Frequencies analysis ...");
		
		final int oldThreads=Shared.threads();
		Shared.capThreads(16);

		Timer t=new Timer();

		//Create an instance of this class
		TetramerFrequencies x=new TetramerFrequencies(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
		
		Shared.setThreads(oldThreads);
	}

	/**
	 * Constructs TetramerFrequencies analyzer with command-line arguments.
	 * Parses parameters for window size, step, output format, and k-mer length.
	 * Initializes k-mer remapping tables and output streams.
	 * @param args Command-line arguments including input files and parameters
	 */
	public TetramerFrequencies(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		int k_=4;
		Parser parser=new Parser();
		for (String arg : args) {
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if (a.equals("s") || a.equals("step")){
				step = Integer.parseInt(b);
			} else if (a.equals("w") || a.equals("window")){
				winSize = Integer.parseInt(b);
			} else if (a.equals("out") || a.equals("freq")){
				out1 = b;
			} else if (a.equals("dropshort")){
				keepShort=!Parse.parseBoolean(b);
			} else if (a.equals("keepshort") || a.equals("short")){
				keepShort=Parse.parseBoolean(b);
			} else if (a.equalsIgnoreCase("gc")){
				printGC=Parse.parseBoolean(b);
			} else if (a.equalsIgnoreCase("float") || a.equals("floats") || a.equals("fraction")){
				printFloats=Parse.parseBoolean(b);
			} else if (a.equalsIgnoreCase("gccomp") || a.equals("comp") || a.equals("compensate")
					|| a.equalsIgnoreCase("gccompensate")|| a.equalsIgnoreCase("gccompensated") 
					|| a.equals("compensated")){
				gcCompensate=Parse.parseBoolean(b);
			} else if (a.equals("k")){
				k_ = Integer.parseInt(b);
			} else if(parser.parse(arg, a, b)){
				//do nothing
			} else {
				throw new RuntimeException("Unknown argument "+arg);
			}
		}

		{//Process parser fields
			Parser.processQuality();

			maxReads=parser.maxReads;
			in1=parser.in1;
		}

		k=k_;
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
		assert(ffin1!=null) : "No input file.";
		assert(ffin1.exists() && ffin1.canRead()) : "Cannot read input file "+in1+".";


		// initialize the class member textstream writer here, so no need of keep results in memory
		if (out1==null || out1.equals("")){
			out1 = "stdout";
		}
		
		remap=BinObject.makeRemap(k);
		canonicalKmers=Tools.max(remap)+1;
		gcmap=BinObject.gcmap(k, remap);
		
		threads=Tools.mid(1, Shared.threads(), 32);
		inq=new ArrayBlockingQueue<Line>(threads+1);
		
		FileFormat ff=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, true, false, true);
		bsw = new ByteStreamWriter(ff);
		bsw.start();
		
		ByteBuilder header=header().nl();
		bsw.add(header, nextID);
		nextID++;
	}

	/**
	 * Main processing method that executes tetramer frequency analysis.
	 * Spawns worker threads, reads input sequences, and generates windowed
	 * k-mer profiles for each sequence.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		
		ArrayList<PrintThread> alpt=spawnThreads();

		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start();
		}

		long readsProcessed=0;
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);

			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}

				for (Read r1 : reads){
					windowedTetramerProfile(r1.bases, r1.id);
					readsProcessed++;
				}

				cris.returnList(ln);
				if(verbose){
					outstream.println("Returned a list.");
				}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}

			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		putLine(POISON_LINE);
		waitForFinish(alpt);

		// close the output file
		bsw.poisonAndWait();

		ReadWrite.closeStream(cris);
		if(verbose){outstream.println("Finished.");}

		t.stop();
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+Tools.format("%.2fk reads/sec", (readsProcessed/(double)(t.elapsed))*1000000));
	}

	/**
	 * Generates tetramer frequency profiles using sliding windows across sequence.
	 * Creates Line objects for each window position and queues them for processing.
	 * @param bases Sequence bases to analyze
	 * @param header Sequence identifier/header
	 */
	private void windowedTetramerProfile(byte[] bases, String header){
		int sidx = 0;
		int eidx = winSize<1 ? bases.length : Tools.min(bases.length, winSize);
		
		while (eidx <= bases.length){
			Line line=new Line(header, bases, sidx, eidx, nextID);
			putLine(line);
			sidx+=windowsPerLine*step;
			eidx+=windowsPerLine*step;
			nextID++;
		}
	}
	
	/**
	 * Formats and appends k-mer frequency data to output buffer.
	 * Handles both count and fractional output formats with optional GC compensation.
	 *
	 * @param line Line object containing sequence window information
	 * @param counts K-mer count array for the window
	 * @param gc GC content fraction for the window
	 * @param bb ByteBuilder for output formatting
	 */
	void append(Line line, int[] counts, float gc, ByteBuilder bb){
		bb.append(line.header);
		bb.tab();
		bb.append(line.sidx+1);
		bb.append('-');
		bb.append(line.eidx);
		if(printGC) {
			bb.tab().append(gc, 4);
		}
		if(printFloats) {
			if(gcCompensate) {
				float[] freqs=Vector.compensate(counts, k, gcmap);
				for (float f: freqs){
					bb.tab();
					bb.append(f, 5);
				}
			}else {
				long sum=Tools.sum(counts);
				float inv=1f/Math.max(1, sum);
				for (int cnt: counts){
					bb.tab();
					bb.append(cnt*inv, 5);
				}
			}
		}else {
			if(gcCompensate) {
				long sum=Tools.sum(counts);
				float[] freqs=Vector.compensate(counts, k, gcmap);
				for (float f: freqs){
					bb.tab();
					bb.append(Math.round(f*sum));
				}
			}else {
				for (int cnt: counts){
					bb.tab();
					bb.append(cnt);
				}
			}
		}
		bb.nl();
	}
	
	// factor this out so we can work on reads
	/**
	 * Counts k-mers in a sequence window and calculates GC content.
	 * Uses rolling hash approach with canonical k-mer representation.
	 * Skips k-mers containing ambiguous bases (N).
	 *
	 * @param bases Sequence bases array
	 * @param startidx Starting index of window (inclusive)
	 * @param endidx Ending index of window (exclusive)
	 * @param counts Output array for k-mer counts (modified in-place)
	 * @param gc Output array for GC content [0] (modified in-place)
	 * @return The counts array (same reference as input)
	 */
	public int[] countKmers(byte[] bases, int startidx, int endidx, int[] counts, float[] gc){
		final int mask=(1<<(2*k))-1;
		final int[] acgtn=new int[5];
		
		int len=0, kmer=0;
		for (int i=startidx; i<endidx; i++){
			int x = AminoAcid.baseToNumberACGTother[bases[i]];
			acgtn[x]++;
			kmer=((kmer<<2)|x)&mask;   // this can produce -1 vaue if any base in tetramer is N!
			if(x<4) {
				len++;
				if (len>=k){
					int idx = remap[kmer];
					counts[idx]++; 
				}
			}else {
				len=0;
				kmer=0;
			}
		}

		int gcCount=acgtn[1]+acgtn[2];
		int atCount=acgtn[0]+acgtn[3];
		gc[0]=gcCount/Math.max(1f, (gcCount+atCount));
		return counts;
	}
	
	/**
	 * Generates column header for output file.
	 * Includes scaffold, range, optional GC column, and canonical k-mer strings.
	 * @return ByteBuilder containing formatted header line
	 */
	public ByteBuilder header() {
		ByteBuilder bb=new ByteBuilder();

		bb.append("scaffold");
		bb.append("\trange");
		if(printGC) {bb.append("\tGC");}

		final int limit=(1<<(2*k));
		for(int kmer=0; kmer<limit; kmer++) {
			if(kmer<=AminoAcid.reverseComplementBinaryFast(kmer, k)) {
				bb.tab().append(AminoAcid.kmerToString(kmer, k));
			}
		}
		return bb;
	}

	/** Prints usage information and command-line parameter help. */
	public static void printHelp(){
		List<String> helplist = new ArrayList<String>();
		helplist.add("Program Name : TetramerFrequencies v1.1");
		helplist.add("Usage : ");
		helplist.add(" -h : this page");
		helplist.add(" -s : step [500]");
		helplist.add(" -w : window size [2000]. If set to 0 the whole sequence is processed");
		System.out.println(String.join("\n", helplist));
	}

	/*--------------------------------------------------------------*/
	
	/** Takes next Line from processing queue, blocking until available.
	 * @return Next Line object to process */
	final Line takeLine(){
		Line line=null;
		while(line==null){
			try {
				line=inq.take();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
//		System.err.println("takeLine("+line.id+")");
		return line;
	}
	
	/** Puts Line object into processing queue for worker threads.
	 * @param line Line object to queue for processing */
	final void putLine(Line line){
//		System.err.println("putLine("+line.id+")");
		while(line!=null){
			try {
				inq.put(line);
				line=null;
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	/** Spawn process threads */
	private ArrayList<PrintThread> spawnThreads(){
		
		//Do anything necessary prior to processing
		
		//Fill a list with PrintThreads
		ArrayList<PrintThread> alpt=new ArrayList<PrintThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new PrintThread());
		}
		if(verbose){outstream.println("Spawned threads.");}
		
		//Start the threads
		for(PrintThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
		
		//Do anything necessary after processing
		return alpt;
	}
	
	/** Waits for all worker threads to complete processing.
	 * @param alpt List of PrintThread objects to wait for */
	private void waitForFinish(ArrayList<PrintThread> alpt){
		//Wait for completion of all threads
		boolean allSuccess=true;
		for(PrintThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
		}
	}
	
	/** Worker thread for processing Line objects and generating k-mer frequency output.
	 * Each thread maintains its own count arrays and processes windows independently. */
	private class PrintThread extends Thread {
		
		/** Default constructor for PrintThread */
		PrintThread(){}
		
		@Override
		public void run(){
			Line line=takeLine();
			while(line!=null && line!=POISON_LINE){
				processLine(line);
				line=takeLine();
			}
			putLine(POISON_LINE);
		}
		
		/**
		 * Processes a single Line object to generate k-mer frequency output.
		 * Handles multiple windows per line and writes results to output stream.
		 * @param line Line object containing sequence window information
		 */
		private void processLine(Line line){
			ByteBuilder bb=new ByteBuilder(512);
			for(int i=0; i<windowsPerLine && line.eidx<=line.bases.length; i++){
				Arrays.fill(counts, 0);
				countKmers(line.bases, line.sidx, line.eidx, counts, gc);
				if(keepShort || line.length()>=winSize){
					append(line, counts, gc[0], bb);
				}
				line.sidx+=step;
				line.eidx+=step;
			}
			bsw.add(bb, line.id);
		}
		
		/** Thread-local k-mer count array for canonical k-mers */
		private final int[] counts=new int[canonicalKmers];
		/** Thread-local array for GC content calculation */
		private final float[] gc=new float[1];
	}
	
	/** Represents a sequence window for k-mer frequency analysis.
	 * Contains sequence data, window boundaries, and unique identifier. */
	private class Line {
		
		/**
		 * Constructs Line object with sequence window parameters.
		 *
		 * @param header_ Sequence identifier/header
		 * @param bases_ Sequence bases array
		 * @param sidx_ Starting index of window
		 * @param eidx_ Ending index of window
		 * @param id_ Unique identifier for this line
		 */
		Line(String header_, byte[] bases_, int sidx_, int eidx_, long id_){
			header=header_;
			bases=bases_;
			sidx=sidx_;
			eidx=eidx_;
			id=id_;
		}
		
		/** Calculates length of sequence window.
		 * @return Length of window in bases */
		public int length() {
			return eidx-sidx+1; //Possible bug: should this be eidx-sidx for exclusive end?
		}

		/** Sequence identifier/header string */
		final String header;
		/** Sequence bases array */
		final byte[] bases;
		/** Starting index of current window */
		int sidx;
		/** Ending index of current window */
		int eidx;
		/** Unique identifier for this line object */
		final long id;
		
	}
	
	/*--------------------------------------------------------------*/
	
	/** Sentinel object used to signal thread termination */
	final Line POISON_LINE=new Line(null, null, -1, -1, -1);
	/** Thread-safe queue for Line objects awaiting processing */
	private final ArrayBlockingQueue<Line> inq;
	/** Number of worker threads for parallel processing */
	private final int threads;
	/** Counter for assigning unique IDs to Line objects */
	private long nextID=0;
	
	/*--------------------------------------------------------------*/

	/** Input file path */
	private String in1 = null;
	/** Output file path */
	private String out1 = null;
	/** Output stream writer for results */
	private ByteStreamWriter bsw = null; // for output

	/** Input file format specification */
	private final FileFormat ffin1;

	/** Output stream for status messages */
	private java.io.PrintStream outstream=System.err;

	/*--------------------------------------------------------------*/

	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	/** Step size between consecutive windows */
	int step = 500;
	/** Size of sliding window for k-mer analysis */
	private int winSize = 2000;
	/** Whether to keep windows shorter than winSize */
	private boolean keepShort=false;
	/** Whether to include GC content in output */
	private boolean printGC=false;
	/** Whether to output frequencies as floats instead of counts */
	private boolean printFloats=false;
	/** Whether to apply GC compensation to k-mer frequencies */
	private boolean gcCompensate=false;
	/** K-mer length for frequency analysis */
	private final int k;
	/** Remapping array for canonical k-mer representation */
	private final int[] remap;
	/** GC content mapping array for k-mers */
	private final int[] gcmap;
	/** Number of canonical k-mers for the given k */
	private final int canonicalKmers;

	/*--------------------------------------------------------------*/

	/** Number of windows to process per Line object */
	private final static int windowsPerLine=8;
	/** Enable verbose logging output */
	public static boolean verbose=false;
}
