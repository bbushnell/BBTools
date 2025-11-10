package scalar;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;

import bin.AdjustEntropy;
import clade.Clade;
import clade.SendClade;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.DisplayParams;
import sketch.SendSketch;
import sketch.Sketch;
import sketch.SketchMakerMini;
import sketch.SketchObject;
import sketch.SketchTool;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.KmerTracker;

/**
 * Calculates compositional scalar metrics from sequencing data.
 * Computes GC-independent metrics (HH, CAGA, strandedness, etc.) either globally
 * or using a sliding window to characterize within-genome variance.
 * Outputs mean and standard deviation for each metric.
 *
 * @author Brian Bushnell
 * @date Oct 2, 2025
 */
public class ScalarIntervals {

	/**
	 * Main entry point for the Scalars program.
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		ScalarIntervals x=new ScalarIntervals(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructs a Scalars instance and parses command-line arguments.
	 * Supports windowed or global analysis of compositional metrics.
	 * @param args Command-line arguments
	 */
	public ScalarIntervals(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null/*getClass()*/, false);
			args=pp.args;
			outstream=pp.outstream;
		}

		Parser parser=new Parser();
		String stdout="stdout.txt";
		parser.out1=stdout;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("header") || a.equals("colheader") || a.equals("columnheader")){
				header=Parse.parseBoolean(b);
			}else if(a.equals("raw")){
				raw=Parse.parseBoolean(b);
			}else if(a.equals("window")){
				window=Parse.parseIntKMG(b);
			}else if(a.equals("interval")){
				interval=Parse.parseIntKMG(b);
			}else if(a.equals("shred")){
				interval=window=Parse.parseIntKMG(b);
			}else if(a.equals("minlen")){
				minlen=Parse.parseIntKMG(b);
			}else if(a.equals("break")){
				breakOnContig=Parse.parseBoolean(b);
			}else if(a.equals("sketch") | a.equals("bbsketch")){
				ScalarData.makeSketch=Parse.parseBoolean(b);
			}else if(a.equals("clade") || a.equals("quickclade")){
				ScalarData.makeClade=Parse.parseBoolean(b);
			}else if(a.equals("printname") || a.equals("printnames")){
				printName=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("printPos") || a.equals("pos")){
				printPos=Parse.parseBoolean(b);
			}else if(a.equals("parsetaxid") || a.equals("parsetax") || a.equals("parsetid")){
				ScalarData.parseTID=Parse.parseBoolean(b);
			}else if(a.equals("mt")){
				mt=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("sendInThread")){
				sendInThread=Parse.parseBoolean(b);
			}else if(a.equals("concurrency")){
				SendClade.maxConcurrency=SendSketch.maxConcurrency=Integer.parseInt(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("printTime") || a.equalsIgnoreCase("time")){
				printTime=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(new File(arg).exists() && FileFormat.isSequence(arg)) {
				in.add(arg);
			}else if(parser.out1==stdout && i>0 && Tools.looksLikeOutputStream(arg) && !FileFormat.isSequence(arg)){
				parser.out1=arg;
			}else{
				//				throw new RuntimeException("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				outstream.println("Unknown parameter "+args[i]);
			}
		}

		{//Process parser fields
			Parser.processQuality();

			maxReads=parser.maxReads;
			if(parser.in1!=null) {in.add(parser.in1);}
			out=parser.out1;
		}
		
		
		in=Tools.getFileOrFiles(in, true, false, false, false);
		ffout=FileFormat.testOutput(out, FileFormat.TXT, null, true, true, false, false);
	}

	/**
	 * Processes input reads and calculates compositional metrics.
	 * Either accumulates global dimer counts or builds histograms from sliding windows.
	 * @param t Timer for performance tracking
	 */
	void process(Timer t){

		readsProcessed=0;
		basesProcessed=0;
		if(verbose) {outstream.println("callingToIntervals");}
		
		if(verbose){outstream.println("Finished reading data; printing to "+out);}
		
		if(header && true) {System.err.print(ScalarData.header(true, printName, false));}
		
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout);
		for(int i=0; i<in.size(); i++) {
			FileFormat ffin=FileFormat.testInput(in.get(i), FileFormat.FASTA, null, true, true);
			ScalarData data=toIntervals(ffin, window, interval, minlen, breakOnContig, maxReads);
			data.print(bsw, printName, header && i==0, interval>0 && printPos, interval);
			if(true) {
				System.err.print(data.mean(true, ffin.name()));
				System.err.print(data.stdev(true, ffin.name()));
			}
		}
		if(bsw!=null) {bsw.poison();}

		t.stop();
		if(printTime) {
			outstream.println();
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		}
		assert(!errorState) : "An error was encountered.";
	}
	
	public static ScalarData toIntervals(String fname, int window, int interval, int minlen, boolean breakOnContig, long maxReads) {
		if(verbose) {System.err.println("callingToIntervals(string)");}
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		return toIntervals(ff, window, interval, minlen, breakOnContig, maxReads);
	}
	
	public static ScalarData toIntervals(FileFormat ff, int window, int interval, int minlen, boolean breakOnContig, long maxReads) {
		if(verbose) {System.err.println("callingToIntervals(ff)");}
		int tid=(ScalarData.parseTID ? bin.BinObject.parseTaxID(ff.name()) : -1);
		final ConcurrentReadInputStream cris;
		cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ff, null);
		cris.start();
		ScalarData data=toIntervals(cris, window, interval, minlen, tid, breakOnContig, maxReads);
		boolean errorState=ReadWrite.closeStreams(cris);
		if(verbose){System.err.println("Finished reading data.");}
		if(errorState){System.err.println("Something went wrong reading "+ff.name());}
		return data;
	}
	
	public static ScalarData toIntervals(ConcurrentReadInputStream cris, int window, int interval, int minlen, 
		int tid, boolean breakOnContig, long maxReads) {
		if(ScalarData.makeSketch || ScalarData.makeClade) {
			ArrayList<ScalarData> list=toIntervals_multi(
				cris, window, interval, minlen, tid, breakOnContig, maxReads);
			return collapse(list, mt);
		}
		if(verbose) {System.err.println("callingToIntervals(cris)");}
		ScalarData data=new ScalarData(printName, -1);
		final KmerTracker dimers=new KmerTracker(2, window);
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			while(ln!=null && reads!=null && reads.size()>0){
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx), r2=r1.mate;
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();

					data.add(r1, dimers, null, interval, minlen, tid, breakOnContig);
					if(r2!=null) {data.add(r2, dimers, null, interval, minlen, tid, breakOnContig);}
				}

				cris.returnList(ln);
				if(verbose){System.err.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		if(verbose) {System.err.println("finished ToIntervals(cris)");}
		return data;
	}
	
	public static ScalarData collapse(ArrayList<ScalarData> list, boolean sort) {
		if(list==null || list.isEmpty()) {return null;}
		if(sort) {Shared.sort(list);}
		ScalarData data=list.get(0);
		for(int i=1; i<list.size(); i++) {
			ScalarData sd=list.set(i, null);
			data.add(sd);
		}
		return data;
	}
	
	public static ArrayList<ScalarData> toIntervals_multi(ConcurrentReadInputStream cris, 
		int window, int interval, int minlen, int tid, boolean breakOnContig, long maxReads) {
		
		if(mt) {return toIntervalsMT(cris, window, interval, minlen, tid, breakOnContig, maxReads);}
		if(verbose) {System.err.println("callingToIntervals(cris)");}
		Timer t=new Timer();
		
		if(ScalarData.makeClade) {
			Clade.DELETE_COUNTS=false;
			AdjustEntropy.load();
		}
		SketchMakerMini smm=null;
		if(ScalarData.makeSketch) {
			SketchObject.processSSU=false;
			SketchObject.postParse();
			DisplayParams params=new DisplayParams();
			params.format=DisplayParams.FORMAT_JSON;
			params.taxLevel=TaxTree.GENUS;
			params.maxRecords=2;
			SketchTool tool=new SketchTool(2000, params);
			smm=new SketchMakerMini(tool, SketchObject.ONE_SKETCH, params);
		}
		final KmerTracker dimers=new KmerTracker(2, window);
		ArrayList<ScalarData> sdlist=new ArrayList<ScalarData>();
//		System.err.println("Made KmerTracker "+2+", "+window);
		{
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			while(ln!=null && reads!=null && reads.size()>0){
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();
					
					if(r1.length()<minlen) {continue;}
					ScalarData sd=new ScalarData(printName, -1);
					sd.add(r1, dimers, smm, interval, minlen, tid, breakOnContig);
					if(ScalarData.makeSketch) {
						if(r1.pairLength()>=ScalarData.minSketchSize) {
							sd.sketch=smm.toSketch(1);
						}else {smm.heap().clear(false);}
					}
					if(sd.clade!=null) {
						assert(!Clade.DELETE_COUNTS);
						sd.clade.finish();
					}
					sdlist.add(sd);
				}

				cris.returnList(ln);
				if(verbose){System.err.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		String fname=(cris.fname.replace(",null", ""));
		t.stopAndStart("Processed "+cris.readsIn()+" sequences from "+fname);
		if(verbose) {System.err.println("finished ToIntervals(cris)");}
		
		sendAndLabel(sdlist);
		
		return sdlist;
	}
	
	static boolean sendAndLabel(ArrayList<ScalarData> sdlist) {
		if(sdlist==null || sdlist.isEmpty()) {return true;}
		Timer t=new Timer();
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		ArrayList<Clade> clades=new ArrayList<Clade>();
		for(ScalarData sd : sdlist) {
			if(sd.sketch!=null) {sketches.add(sd.sketch);}
			if(sd.clade!=null) {clades.add(sd.clade);}
		}
		boolean cladeSuccess=true, sketchSuccess=true;
		if(!clades.isEmpty()) {
			cladeSuccess=SendClade.sendAndLabel(clades);
			if(printTime && (!sendInThread || !mt)){t.stopAndStart("Sent "+clades.size()+" clades.");}
			for(ScalarData sd : sdlist) {
				if(sd.clade!=null && sd.clade.taxID>0) {
					sd.taxIDs.fill(sd.clade.taxID);
				}
			}
		}
		if(!sketches.isEmpty()) {
			sketchSuccess=SendSketch.sendAndLabel(sketches, null);
			if(printTime && (!sendInThread || !mt)){t.stopAndStart("Sent "+sketches.size()+" sketches.");}
			for(ScalarData sd : sdlist) {
				if(sd.sketch!=null && sd.sketch.taxID>0) {
					sd.taxIDs.fill(sd.sketch.taxID);
				}
			}
		}
		return cladeSuccess && sketchSuccess;
	}
	
	public static ArrayList<ScalarData> toIntervalsMT(ConcurrentReadInputStream cris, 
		int window, int interval, int minlen, 
		int tid, boolean breakOnContig, long maxReads) {
		if(verbose) {System.err.println("callingToIntervals(cris)");}
		Timer t=new Timer();
		
		if(ScalarData.makeClade) {
			Clade.DELETE_COUNTS=false;
			AdjustEntropy.load();
		}
		DisplayParams params=null;
		if(ScalarData.makeSketch) {
			SketchObject.processSSU=false;
			SketchObject.postParse();
			params=new DisplayParams();
			params.format=DisplayParams.FORMAT_JSON;
			params.taxLevel=TaxTree.GENUS;
			params.maxRecords=2;
		}

		ArrayList<ScalarData> sdlist=spawnThreads(cris, params,  window, 
			interval, minlen, tid, breakOnContig);
		
		if(sendInThread && sdlist!=null) {
			int s=0, c=0;
			for(ScalarData sd : sdlist) {
				s+=(sd.sketch==null ? 0 : 1);
				c+=(sd.clade==null ? 0 : 1);
			}
			if(printTime && (ScalarData.makeClade || ScalarData.makeSketch)) {
				if(c>0){t.stopAndStart("Sent "+c+" clades.");}
				else if(s>0){t.stopAndStart("Sent "+s+" sketches.");}
			}
		}else{
			String fname=(cris.fname.replace(",null", ""));
			t.stopAndStart("Processed "+cris.readsIn()+" sequences from "+fname);
			boolean success=sendAndLabel(sdlist);
			if(!success) {new Exception().printStackTrace();}
		}
		
		return sdlist;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn process threads */
	private static ArrayList<ScalarData> spawnThreads(final ConcurrentReadInputStream cris, DisplayParams params, 
		int window, int interval, int minlen, int taxID, boolean breakOnContig){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		final AtomicInteger active=new AtomicInteger(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, params, window, interval, minlen, taxID, breakOnContig, i, active));
		}
		
		//Start the threads and wait for them to finish
		SDAccumulator sda=new SDAccumulator();
		boolean success=ThreadWaiter.startAndWait(alpt, sda);
		errorState&=!success;
		
		//Do anything necessary after processing
		return sda.list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private static class SDAccumulator implements Accumulator<ProcessThread> {
		
		@Override
		public final void accumulate(ProcessThread pt){
			synchronized(pt) {
				list.addAll(pt.sdlist);
				synchronized(ScalarIntervals.class) {
					readsProcessed+=pt.readsProcessedT;
					basesProcessed+=pt.basesProcessedT;
					errorState|=(!pt.success);
				}
			}
		}
		
		@Override
		public final boolean success(){
			return !errorState;
		}

		@Override
		public ReadWriteLock rwlock(){ // TODO Auto-generated method stub
			return null;
		}
		
		ArrayList<ScalarData> list=new ArrayList<ScalarData>();
		
	}
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
	static class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ConcurrentReadInputStream cris_, final DisplayParams params_, 
			int window_, int interval_, int minlen_, int taxID_, boolean breakOnContig_, final int tid_,
			AtomicInteger active_){
			cris=cris_;
			tid=tid_;

			params=params_;
			window=window_;
			interval=interval_;
			minlen=minlen_;
			taxID=taxID_;
			breakOnContig=breakOnContig_;
			active=active_;
		}
		
		//Called by start()
		@Override
		public void run(){
			
			if(ScalarData.makeSketch) {
				SketchTool tool=new SketchTool(2000, params);
				smm=new SketchMakerMini(tool, SketchObject.ONE_SKETCH, params);
			}
			dimers=new KmerTracker(2, window);
			
			//Process the reads
			success=processInner();
			
			//Do anything necessary after processing
		}
		
		/** Iterate through the reads */
		boolean processInner(){
			
			//Grab the first ListNum of reads
			ListNum<Read> ln=cris.nextList();

			//As long as there is a nonempty read list...
			while(ln!=null && ln.size()>0){
//				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");} //Disabled due to non-static access
				
				processList(ln);
				
				//Notify the input stream that the list was used
				cris.returnList(ln);
//				if(verbose){outstream.println("Returned a list.");} //Disabled due to non-static access
				
				//Fetch a new list
				ln=cris.nextList();
			}

			//Notify the input stream that the final list was used
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
			
			if(active.decrementAndGet()<1){
				String fname=(cris.fname.replace(",null", ""));
				t.stopAndStart("Processed "+cris.readsIn()+" sequences from "+fname);
			}
			
			if(sdlist.isEmpty()) {return true;}
			if(sendInThread) {success=sendAndLabel(sdlist);}
			return success;
		}
		
		void processList(ListNum<Read> ln){

			//Grab the actual read list from the ListNum
			final ArrayList<Read> reads=ln.list;
			
			//Loop through each read in the list
			for(int idx=0; idx<reads.size(); idx++){
				final Read r1=reads.get(idx);
				final Read r2=r1.mate;
				
				//Validate reads in worker threads
				if(!r1.validated()){r1.validate(true);}
				if(r2!=null && !r2.validated()){r2.validate(true);}

				//Track the initial length for statistics
				final int initialLength1=r1.length();
				final int initialLength2=r1.mateLength();

				//Increment counters
				readsProcessedT+=r1.pairCount();
				basesProcessedT+=initialLength1+initialLength2;
				processRead(r1);
				processRead(r2);
			}
		}
		
		/**
		 * Process a read or a read pair.
		 * @param r Read
		 * @return True if the read was processed.
		 */
		boolean processRead(final Read r){
			if(r==null || r.length()<minlen) {return false;}
			ScalarData sd=new ScalarData(printName, -1);
			sd.add(r, dimers, smm, interval, minlen, taxID, breakOnContig);
			if(ScalarData.makeSketch) {
				if(r.pairLength()>=ScalarData.minSketchSize) {
					sd.sketch=smm.toSketch(1);
				}else {smm.heap().clear(false);}
			}
			if(sd.clade!=null) {
				assert(!Clade.DELETE_COUNTS);
				sd.clade.finish();
			}
			sdlist.add(sd);
			return true;
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		private final DisplayParams params;
		private final int window;
		private final int interval;
		private final int minlen;
		private final int taxID;
		private final boolean breakOnContig;
		private SketchMakerMini smm;
		private KmerTracker dimers;
		
		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream */
		ArrayList<ScalarData> sdlist=new ArrayList<ScalarData>();
		/** Thread ID */
		final int tid;
		final AtomicInteger active;
		final Timer t=new Timer();
	}
	
	
	/*--------------------------------------------------------------*/

	/** Input file path */
	private ArrayList<String> in=new ArrayList<String>();
	/** Output file path */
	private String out=null;

//	/** Input file format */
//	private final FileFormat ffin;
	/** Output file format */
	private final FileFormat ffout;
	/** Whether to print column headers */
	private boolean header=false;
	/** Whether to print row headers */
	private boolean rowheader=false;
	private boolean raw=false;


	/** Window size for sliding window analysis (0 for global analysis) */
	private int window=50000;
	private int interval=10000;
	private int minlen=500;
	private boolean breakOnContig=true;
	static boolean printName=false;
	private static boolean printPos=false;
	/** Whether to print timing information */
	private static boolean printTime=true;
	static boolean mt=true;
	static boolean sendInThread=true;
	
	static long readsProcessed=0, basesProcessed=0;
	
	/*--------------------------------------------------------------*/

	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;
	/** Whether an error occurred during processing */
	private static boolean errorState=false;

	/*--------------------------------------------------------------*/

	/** Output stream for messages */
	private java.io.PrintStream outstream=System.err;
	/** Whether to print verbose progress messages */
	public static boolean verbose=false;

}
