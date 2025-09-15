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
 * 
 * @author Brian Bushnell
 * @date December 11, 2024
 *
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
			SketchObject.minSketchSize=5;
			
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
	
	void sketch(ArrayList<? extends Sketchable> input, boolean force) {
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
	
	/** Spawn process threads */
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
	
	@Override
	public final void accumulate(ProcessThread pt){
//		linesProcessed+=pt.linesProcessedT;
//		bytesProcessed+=pt.bytesProcessedT;
//		linesOut+=pt.linesOutT;
//		bytesOut+=pt.bytesOutT;
		errorState|=(!pt.success);
		errorState|=(pt.errorStateT);
	}
	
	@Override
	public final boolean success(){return !errorState;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This class is static to prevent accidental writing to shared variables.
	 * It is safe to remove the static modifier. */
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
		
		/** Iterate through the lines */
		void processInner_oneByOne(){
//			Timer t=new Timer();
			for(int i=tid; i<contigs.size(); i+=threads) {
				Sketchable c=contigs.get(i);
				synchronized(c) {
					assert(c.id()==i);
					Sketch sketch=c.toSketch(smm, dummy);

					String results=SendSketch.sendSketch(sketch, "refseq", params, 0);
					if(results==null) {continue;}

					JsonObject all=jp.parseJsonObject(results);
					c.setFrom(all);
					assert(c.sketchedSize()==c.size());
				}
			}
//			t.stop("Thread "+tid+" time: ");
		}
		
		/** Iterate through the lines */
		void processInner_bulk(){
			final int incr=sectionSize*threads;
			for(int i=tid; i<contigs.size(); i+=incr) {processSection(i, i+incr);}
		}
		
		/** Iterate through the lines */
		void processSection(final int from, int to){
			ArrayList<Sketch> sketches=new ArrayList<Sketch>(1+contigs.size()/threads);
			for(int i=from; i<contigs.size() && i<to; i+=threads) {
				Sketchable c=contigs.get(i);
				synchronized(c) {
					assert(c.id()==i);
					Sketch sketch=c.toSketch(smm, dummy);
					assert(sketch!=null) : "Handle null sketches.";
					sketches.add(sketch);//Note:  Could potentially be null?
				}
			}
//			t.stopAndStart("Thread "+tid+" sketch time: ");
			ArrayList<JsonObject> results=SendSketch.sendSketches(sketches, "refseq", params);
			assert(results.size() == sketches.size()) : results.size()+", "+sketches.size();
			for(int i=from, j=0; i<contigs.size() && i<to; i+=threads, j++) {
				Sketchable c=contigs.get(i);
				synchronized(c) {
					assert(c.id()==i);
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
		
		protected boolean errorStateT=false;
		
		/** True only if this thread has completed successfully */
		boolean success=false;
		
		/** Input */
		private final ArrayList<? extends Sketchable> contigs;
		/** Thread ID */
		final int tid;
		/** Thread ID */
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
	
	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	
}
