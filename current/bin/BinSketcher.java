package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import json.JsonObject;
import json.JsonParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import sketch.DisplayParams;
import sketch.SendSketch;
import sketch.Sketch;
import sketch.SketchMakerMini;
import sketch.SketchObject;
import sketch.SketchTool;
import stream.Read;
import tax.TaxTree;
import template.Accumulator;
import template.ThreadWaiter;

/**
 * Handles sketches and taxonomic assignments for contigs and clusters.
 * Creates MinHash sketches for genomic sequences and queries them against
 * taxonomic databases to determine taxonomic identity. Supports both individual
 * and bulk processing modes with multi-threading.
 *
 * @author Brian Bushnell
 * @date December 11, 2024
 */
public class BinSketcher extends BinObject implements Accumulator<BinSketcher.ProcessThread> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public BinSketcher(int threads_, int minSize_){
		
		threads=Tools.min(threads_, Shared.threads());
		minSize=minSize_;
		
		if(sketchClusters || sketchContigs || sketchOutput){
			SketchObject.AUTOSIZE_LINEAR_DENSITY=sketchDensity;
			SketchObject.AUTOSIZE_LINEAR=true;
			SketchObject.AUTOSIZE=false;
			SketchObject.SET_AUTOSIZE=true;
			SketchObject.minSketchSize=3;
			
//			SketchObject.AUTOSIZE=false;
//			SketchObject.defaultParams.minKeyOccuranceCount=2;
			SketchObject.defaultParams.parse("trackcounts", "trackcounts", null);
//			SketchObject.defaultParams.minProb=0;
			SketchObject.postParse();
			SketchObject.defaultParams.maxRecords=2;
			SketchObject.defaultParams.taxLevel=TaxTree.GENUS;
			tool=new SketchTool(sketchSize, SketchObject.defaultParams);
		}else{
			tool=null;
		}
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
//	void sketchBins(ArrayList<Bin> input, boolean force) {
//		sketch(input, force);
//	}
	
	public void sketch(ArrayList<? extends Sketchable> input, boolean force) {
		ArrayList<Sketchable> updateList=new ArrayList<Sketchable>();
		float mult=(force ? 1 : 2);
		for(Sketchable s : input) {
			if(s.size()>=minSize) {
				synchronized(s) {
					if(s.size()>mult*s.sketchedSize()) {
						s.clearTax();
						updateList.add(s);
					}
				}
			}
		}
		if(updateList.isEmpty()) {return;}
		spawnThreads(updateList);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Spawns and manages ProcessThreads for sketching operations.
	 * Determines optimal thread count based on data size and system resources,
	 * then coordinates parallel processing of sketch generation.
	 * @param list List of Sketchable objects to process
	 */
	private void spawnThreads(ArrayList<? extends Sketchable> list){
		Timer t=new Timer(outstream, true);
		outstream.print("Sketching "+list.size()+" elements: \t");
		
		//Do anything necessary prior to processing
		long bases=0;
		for(Sketchable s : list) {bases+=s.size();}
		
		//Determine how many threads may be used
		final int pthreads=(int)Tools.max(1, Tools.min(threads, Shared.threads(), list.size()/4, bases/40000));
		
		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(pthreads);
		for(int i=0; i<pthreads; i++){
			alpt.add(new ProcessThread(list, i, pthreads));
		}
		assert(alpt.size()==pthreads);
		
		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;
		
		//Do anything necessary after processing
		t.stopAndPrint();
	}
	
	/**
	 * Accumulates results from a completed ProcessThread.
	 * Updates error state based on thread completion status.
	 * @param pt The ProcessThread whose results should be accumulated
	 */
	@Override
	public final void accumulate(ProcessThread pt){
//		linesProcessed+=pt.linesProcessedT;
//		bytesProcessed+=pt.bytesProcessedT;
//		linesOut+=pt.linesOutT;
//		bytesOut+=pt.bytesOutT;
		errorState|=(!pt.success);
		errorState|=(pt.errorStateT);
	}
	
	/** Returns whether all processing completed successfully.
	 * @return true if no errors occurred during processing, false otherwise */
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Worker thread for processing sketches of Sketchable objects.
	 * Handles either individual or bulk sketch generation and taxonomic
	 * assignment via remote database queries.
	 */
	class ProcessThread extends Thread {
		
		//Constructor
		ProcessThread(final ArrayList<? extends Sketchable> contigs_, final int tid_, final int threads_){
			contigs=contigs_;
			tid=tid_;
			threads=threads_;
			params=new DisplayParams();
			params.format=DisplayParams.FORMAT_JSON;
			params.taxLevel=TaxTree.GENUS;
			smm=new SketchMakerMini(tool, SketchObject.ONE_SKETCH, params);
		}
		
		//Called by start()
		/**
		 * Main thread execution method.
		 * Chooses between bulk or individual processing modes based on
		 * data size and configuration, then processes assigned sketches.
		 */
		@Override
		public void run(){
			//Do anything necessary prior to processing
			
			//Process the contigs
			if(sketchInBulk && (1+contigs.size()/threads)>2) {
				processInner_bulk();
			}else {
				processInner_oneByOne();
			}
			
			//Do anything necessary after processing
			
			//Indicate successful exit status
			success=true;
		}
		
		/**
		 * Processes sketches individually, one at a time.
		 * Creates sketch for each assigned object, sends to taxonomic database,
		 * and updates object with taxonomic assignment results.
		 */
		void processInner_oneByOne(){
//			Timer t=new Timer();
			for(int i=tid; i<contigs.size(); i+=threads) {
				Sketchable c=contigs.get(i);
				synchronized(c) {
					Sketch sketch=c.toSketch(smm, dummy);
					if(send) {
						String results=SendSketch.sendSketch(sketch, "refseq", params, 0);
						if(results==null) {continue;}

						JsonObject all=jp.parseJsonObject(results);
						c.setFrom(all);
					}
					assert(c.sketchedSize()==c.size());
				}
			}
//			t.stop("Thread "+tid+" time: ");
		}
		
		/**
		 * Processes sketches in bulk sections for efficiency.
		 * Divides assigned work into sections and processes each section
		 * as a batch to reduce network overhead.
		 */
		void processInner_bulk(){
			final int incr=sectionSize*threads;
			for(int i=tid; i<contigs.size(); i+=incr) {processSection(i, i+incr);}
		}
		
		/**
		 * Processes a section of sketches as a batch.
		 * Creates all sketches in the range, sends them together to the
		 * taxonomic database, then updates objects with results.
		 *
		 * @param from Starting index for this section (inclusive)
		 * @param to Ending index for this section (exclusive)
		 */
		void processSection(final int from, int to){
			ArrayList<Sketch> sketches=new ArrayList<Sketch>(1+contigs.size()/threads);
			for(int i=from; i<contigs.size() && i<to; i+=threads) {
				Sketchable c=contigs.get(i);
				synchronized(c) {
					Sketch sketch=c.toSketch(smm, dummy);
//					assert(sketch!=null) : "Handle null sketches.";//Handled!
					sketches.add(sketch);//Can be null
				}
			}
//			t.stopAndStart("Thread "+tid+" sketch time: ");
			if(!send) {return;}
			ArrayList<JsonObject> results=SendSketch.sendSketches(sketches, "refseq", params);
			assert(results.size() == sketches.size()) : results.size()+", "+sketches.size();
			for(int i=from, j=0; i<contigs.size() && i<to; i+=threads, j++) {
				Sketchable c=contigs.get(i);
				synchronized(c) {
					JsonObject jo=results.get(j);
					c.setFrom(jo);
					assert(jo==null || c.sketchedSize()==c.size());
				}
			}
		}

//		/** Number of reads processed by this thread */
//		protected long linesProcessedT=0;
//		/** Number of bases processed by this thread */
//		protected long bytesProcessedT=0;
//		
//		/** Number of reads retained by this thread */
//		protected long linesOutT=0;
//		/** Number of bases retained by this thread */
//		protected long bytesOutT=0;
		

		final Read dummy=new Read(null, null, null, 0);
		final JsonParser jp=new JsonParser();
		final DisplayParams params;
		
		/** True only if this thread has completed successfully */
		protected boolean errorStateT=false;
		
		boolean success=false;
		
		/** Input list of Sketchable objects to process */
		private final ArrayList<? extends Sketchable> contigs;
		/** Thread ID for this worker */
		final int tid;
		/** Total number of threads in the processing pool */
		final int threads;
		

		final SketchMakerMini smm;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	long linesProcessed=0;
	long linesOut=0;
	long bytesProcessed=0;
	long bytesOut=0;

	
	private final SketchTool tool;
//	private final SketchMakerMini smm;
	private final int threads;
	final int minSize;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	static int sectionSize=100;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns the read-write lock for thread-safe access.
	 * @return ReadWriteLock for synchronization */
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output stream for status messages */
	private PrintStream outstream=System.err;
	/** Flag to enable verbose output messages */
	public static boolean verbose=false;
	/** True if an error was encountered during processing */
	public boolean errorState=false;
	
	public static boolean send=true;
	
}
