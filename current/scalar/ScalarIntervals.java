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
import parse.Parse;
import parse.Parser;
import parse.PreParser;
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
 * Emits GC/HH/CAGA, coverage and taxonomy rows sampled from sequence input.
 * Window size counts valid dinucleotides; emission spacing counts valid bases.
 * Mean and population-standard-deviation summaries go separately to stderr.
 * Shared options and counters make concurrent independent invocations unsupported;
 * configure the process-wide state before starting workers.
 *
 * @author Brian Bushnell
 * @contributor Neptune
 * @date Oct 2, 2025
 */
public class ScalarIntervals{

	/**
	 * Main entry point for interval-table generation.
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		ScalarIntervals x=new ScalarIntervals(args);

		//Run the object
		try{x.process(t);}
		finally{Shared.closeStream(x.outstream);}
	}

	/**
	 * Constructs an interval generator and parses command-line arguments.
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
			String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && (b.isEmpty() || b.equalsIgnoreCase("null"))){b=null;}

			if(a.equals("header") || a.equals("colheader") || a.equals("columnheader")){
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
			}else if(a.equals("cov") || a.equals("coverage") || a.equals("covfile")){
				covFile=b;
			}else if(a.equals("depth") || a.equals("depthfile") || a.equals("sam") || a.equals("bam")){
				depthFile=b;
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
			}else if(new File(arg).exists() && FileFormat.isSequence(arg)){
				in.add(arg);
			}else if(parser.out1==stdout && i>0 && Tools.looksLikeOutputStream(arg) && !FileFormat.isSequence(arg)){
				parser.out1=arg;
			}else{
				throw new IllegalArgumentException("Unknown parameter "+args[i]);
			}
		}

		{//Process parser fields
			if(parser.in2!=null){throw new IllegalArgumentException("ScalarIntervals does not support in2; use interleaved input or separate input files.");}
			Parser.processQuality();

			maxReads=parser.maxReads;
			if(parser.in1!=null){in.add(parser.in1);}
			out=parser.out1;
		}


		in=Tools.getFileOrFiles(in, true, false, false, false);
		ScalarData.checkFiles(in, out, parser.overwrite, parser.append);
		ffout=FileFormat.testOutput(out, FileFormat.TXT, null, true, parser.overwrite, parser.append, false);
	}

	/**
	 * Processes input reads and calculates compositional metrics.
	 * Writes interval rows and per-file summary statistics.
	 * @param t Timer for performance tracking
	 */
	void process(Timer t){

		readsProcessed=0;
		basesProcessed=0;
		errorState=false;
		if(verbose){outstream.println("callingToIntervals");}

		if(verbose){outstream.println("Finished reading data; printing to "+out);}

		//Column header now goes into the OUTPUT file via data.print(...) below (gated by header && i==0),
		//matching the data rows exactly. The previous System.err header here mismatched the Mean/STDev rows.
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout);
		try{
			for(int i=0; i<in.size(); i++){
				FileFormat ffin=FileFormat.testInput(in.get(i), FileFormat.FASTA, null, true, true);
				ScalarData data=toIntervals(ffin, window, interval, minlen, breakOnContig, maxReads);
				//scalar/ScalarIntervals#001: header and printName are distinct, ordered arguments.
				data.print(bsw, header && i==0, printName, printPos, interval);
				System.err.print(data.mean(true, ffin.name()));
				System.err.print(data.stdev(true, ffin.name()));
			}
		}finally{
			if(bsw!=null){errorState|=bsw.poisonAndWait();}
		}

		t.stop();
		if(printTime){
			outstream.println();
			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		}
		if(errorState){throw new RuntimeException("An error occurred generating scalar intervals");}
	}

	/** Opens a sequence file and returns its sampled rows; owns the input stream.
	 * maxReads is the input reader's fragment limit, not a base or row limit. */
	public static ScalarData toIntervals(String fname, int window, int interval, int minlen, boolean breakOnContig, long maxReads){
		if(verbose){System.err.println("callingToIntervals(string)");}
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTA, null, true, true);
		return toIntervals(ff, window, interval, minlen, breakOnContig, maxReads);
	}

	/** Loads configured coverage once, then owns a reader for this input descriptor. */
	public static ScalarData toIntervals(FileFormat ff, int window, int interval, int minlen, boolean breakOnContig, long maxReads){
		if(verbose){System.err.println("callingToIntervals(ff)");}
		// Load depth map BEFORE processing reads (static, shared by all ScalarData instances)
		if((covFile!=null || depthFile!=null) && ScalarData.depthMap==null){
			ScalarData temp=new ScalarData(false, -1);
			temp.loadDepth(covFile, depthFile);
		}
		int tid=(ScalarData.parseTID ? bin.BinObject.parseTaxID(ff.name()) : -1);
		final ConcurrentReadInputStream cris;
		cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ff, null);
		cris.start();
		try{return toIntervals(cris, window, interval, minlen, tid, breakOnContig, maxReads);}
		finally{
			//scalar/ScalarIntervals#002: a failed close must reach the API caller too.
			final boolean failed=ReadWrite.closeStreams(cris);
			errorState|=failed;
			if(failed){throw new RuntimeException("Error reading scalar input: "+ff.name());}
		}
	}

	/** Consumes a borrowed, already-started reader; its caller must close it.
	 * Uses per-read labeling tables when clade/sketch creation is enabled. */
	public static ScalarData toIntervals(ConcurrentReadInputStream cris, int window, int interval, int minlen,
		int tid, boolean breakOnContig, long maxReads){
		if(ScalarData.makeSketch || ScalarData.makeClade){
			ArrayList<ScalarData> list=toIntervals_multi(
				cris, window, interval, minlen, tid, breakOnContig, maxReads);
			return collapse(list, mt);
		}
		if(verbose){System.err.println("callingToIntervals(cris)");}
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
					if(r2!=null){data.add(r2, dimers, null, interval, minlen, tid, breakOnContig);}
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
		if(verbose){System.err.println("finished ToIntervals(cris)");}
		return data;
	}

	/** Combines rows into the first table and nulls consumed list slots.
	 * Optionally sorts by ScalarData order; null/empty input yields an empty table. */
	public static ScalarData collapse(ArrayList<ScalarData> list, boolean sort){
		if(list==null || list.isEmpty()){return new ScalarData(printName, -1);}
		if(sort){Shared.sort(list);}
		ScalarData data=list.get(0);
		for(int i=1; i<list.size(); i++){
			ScalarData sd=list.set(i, null);
			data.add(sd);
		}
		return data;
	}

	/** Produces separate sequence tables using the configured serial/worker mode.
	 * Borrows the started reader; optional remote labeling uses process-wide settings. */
	public static ArrayList<ScalarData> toIntervals_multi(ConcurrentReadInputStream cris,
		int window, int interval, int minlen, int tid, boolean breakOnContig, long maxReads){

		if(mt){return toIntervalsMT(cris, window, interval, minlen, tid, breakOnContig, maxReads);}
		if(verbose){System.err.println("callingToIntervals(cris)");}
		Timer t=new Timer();

		if(ScalarData.makeClade){
			Clade.DELETE_COUNTS=false;
			AdjustEntropy.load();
		}
		SketchMakerMini smm=null;
		if(ScalarData.makeSketch){
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
					for(int mate=0; mate<2; mate++){
						final Read r=(mate==0 ? r1 : r1.mate);
						if(r==null || r.length()<minlen){continue;}
						ScalarData sd=new ScalarData(printName, -1);
						sd.add(r, dimers, smm, interval, minlen, tid, breakOnContig);
						if(ScalarData.makeSketch){
							if(r.length()>=ScalarData.minSketchSize){sd.sketch=smm.toSketch(1);}
							else{smm.heap().clear(false);}
						}
						if(sd.clade!=null){
							assert(!Clade.DELETE_COUNTS) : "ScalarIntervals must retain counts until clade labeling";
							sd.clade.finish();
						}
						sdlist.add(sd);
					}
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
		if(verbose){System.err.println("finished ToIntervals(cris)");}

		if(!sendAndLabel(sdlist)){throw new RuntimeException("Scalar taxonomy labeling failed");}

		return sdlist;
	}

	/** Sends nonempty classifier/sketch collections and applies successful assignments.
	 * Returns false when either labeling service reports failure; no requests for empty lists. */
	static boolean sendAndLabel(ArrayList<ScalarData> sdlist){
		if(sdlist==null || sdlist.isEmpty()){return true;}
		Timer t=new Timer();
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		ArrayList<Clade> clades=new ArrayList<Clade>();
		for(ScalarData sd : sdlist){
			if(sd.sketch!=null){sketches.add(sd.sketch);}
			if(sd.clade!=null){clades.add(sd.clade);}
		}
		boolean cladeSuccess=true, sketchSuccess=true;
		if(!clades.isEmpty()){
			cladeSuccess=SendClade.sendAndLabel(clades);
			if(printTime && (!sendInThread || !mt)){t.stopAndStart("Sent "+clades.size()+" clades.");}
			for(ScalarData sd : sdlist){
				if(sd.clade!=null && sd.clade.taxID>0){
					sd.taxIDs.fill(sd.clade.taxID);
				}
			}
		}
		if(!sketches.isEmpty()){
			sketchSuccess=SendSketch.sendAndLabel(sketches, null);
			if(printTime && (!sendInThread || !mt)){t.stopAndStart("Sent "+sketches.size()+" sketches.");}
			for(ScalarData sd : sdlist){
				if(sd.sketch!=null && sd.sketch.taxID>0){
					sd.taxIDs.fill(sd.sketch.taxID);
				}
			}
		}
		return cladeSuccess && sketchSuccess;
	}

	/** Collects worker-owned tables, waits for completion and performs deferred labeling. */
	public static ArrayList<ScalarData> toIntervalsMT(ConcurrentReadInputStream cris,
		int window, int interval, int minlen,
		int tid, boolean breakOnContig, long maxReads){
		if(verbose){System.err.println("callingToIntervals(cris)");}
		Timer t=new Timer();

		if(ScalarData.makeClade){
			Clade.DELETE_COUNTS=false;
			AdjustEntropy.load();
		}
		DisplayParams params=null;
		if(ScalarData.makeSketch){
			SketchObject.processSSU=false;
			SketchObject.postParse();
			params=new DisplayParams();
			params.format=DisplayParams.FORMAT_JSON;
			params.taxLevel=TaxTree.GENUS;
			params.maxRecords=2;
		}

		ArrayList<ScalarData> sdlist=spawnThreads(cris, params,  window,
			interval, minlen, tid, breakOnContig);

		if(sendInThread && sdlist!=null){
			int s=0, c=0;
			for(ScalarData sd : sdlist){
				s+=(sd.sketch==null ? 0 : 1);
				c+=(sd.clade==null ? 0 : 1);
			}
			if(printTime && (ScalarData.makeClade || ScalarData.makeSketch)){
				if(c>0){t.stopAndStart("Sent "+c+" clades.");}
				else if(s>0){t.stopAndStart("Sent "+s+" sketches.");}
			}
		}else{
			String fname=(cris.fname.replace(",null", ""));
			t.stopAndStart("Processed "+cris.readsIn()+" sequences from "+fname);
			boolean success=sendAndLabel(sdlist);
			if(!success){throw new RuntimeException("Scalar taxonomy labeling failed");}
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
		errorState|=!success;
		if(!success){throw new RuntimeException("Scalar interval worker failed; no partial table returned");}

		//Do anything necessary after processing
		return sda.list;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Merges completed worker tables/counters; ThreadWaiter serializes accumulator calls. */
	private static class SDAccumulator implements Accumulator<ProcessThread>{

		@Override
		public final void accumulate(ProcessThread pt){
			synchronized(pt){
				list.addAll(pt.sdlist);
				synchronized(ScalarIntervals.class){
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

		/** No external lock is required by this accumulator's single collector. */
		@Override
		public ReadWriteLock rwlock(){
			return null;
		}

		/** Tables transferred from completed workers, without copying their row columns. */
		ArrayList<ScalarData> list=new ArrayList<ScalarData>();

	}

	/** Owns its dimer/sketch state and table list; borrows the shared input reader. */
	static class ProcessThread extends Thread{

		/** Borrows input/configuration; creates computational state in run(). */
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

			if(ScalarData.makeSketch){
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

			if(sdlist.isEmpty()){return true;}
			return !sendInThread || sendAndLabel(sdlist);
		}

		/** Processes both mates, retaining their individual tables and aggregate counters. */
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
			if(r==null || r.length()<minlen){return false;}
			ScalarData sd=new ScalarData(printName, -1);
			sd.add(r, dimers, smm, interval, minlen, taxID, breakOnContig);
			if(ScalarData.makeSketch){
				if(r.length()>=ScalarData.minSketchSize){
					sd.sketch=smm.toSketch(1);
				}else{smm.heap().clear(false);}
			}
			if(sd.clade!=null){
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

		/** Shared immutable sketch configuration after worker startup. */
		private final DisplayParams params;
		/** Number of valid dimers retained in the rolling window. */
		private final int window;
		/** Valid-base spacing between emitted rows. */
		private final int interval;
		/** Minimum retained read/tail size in the existing sampling policy. */
		private final int minlen;
		/** Input-file taxonomic ID, or negative to permit per-header parsing. */
		private final int taxID;
		/** Clear all dimer state between retained sequences. */
		private final boolean breakOnContig;
		/** Worker-owned sketch accumulator, only allocated when enabled. */
		private SketchMakerMini smm;
		/** Worker-owned rolling dimer accumulator. */
		private KmerTracker dimers;

		/** Shared input stream */
		private final ConcurrentReadInputStream cris;
		/** Worker-owned rows, transferred to the collector after completion. */
		ArrayList<ScalarData> sdlist=new ArrayList<ScalarData>();
		/** Thread ID */
		final int tid;
		/** Shared progress counter; the last successful reader prints elapsed input time. */
		final AtomicInteger active;
		/** Timer beginning at worker construction. */
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
	private boolean header=true;
	/** Reserved compatibility fields; neither currently changes emitted interval rows. */
	private boolean rowheader=false;
	/** Reserved raw-output setting; parsed but not used by interval generation. */
	private boolean raw=false;


	/** Valid-dimer capacity for rolling analysis; zero accumulates without a window. */
	private int window=50000;
	/** Valid-base spacing between samples when window is positive. */
	private int interval=10000;
	/** Minimum retained contig and trailing valid-base count. */
	private int minlen=500;
	/** Reset dimer state at each retained contig. */
	private boolean breakOnContig=true;
	/** Depth/coverage support */
	static String covFile=null;  // Coverage file (pileup or covmaker format)
	/** Optional alignment input used when no coverage table is supplied. */
	static String depthFile=null;
	/** Retain and print input names. */
	static boolean printName=false;
	/** Print stored Start and End, including unlimited per-contig measurements. */
	private static boolean printPos=false;
	/** Whether to print timing information */
	private static boolean printTime=true;
	/** Enable worker processing for the optional per-read labeling path. */
	static boolean mt=true;
	/** Perform labeling in each worker rather than after collection. */
	static boolean sendInThread=true;

	/** Aggregate input sequences and bases, counting both mates. */
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
