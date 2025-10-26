package stream;

import java.io.IOException;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;
import template.ThreadWaiter;

/**
 * Writes SAM text files with parallel conversion and ordered output.
 * 
 * Workers convert Read/SamLine objects to SAM text in parallel.
 * JobQueue ensures output blocks are written in order.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date October 25, 2025
 */
public class SamLineWriter extends SamWriter {

	public static void main(String[] args){
		Timer t=new Timer(), t2=new Timer();
		String in=args[0];
		String out=args[1];

		//Create input streamer
		FileFormat ffin=FileFormat.testInput(in, FileFormat.SAM, null, true, false);
		SamStreamer ss=SamStreamer.makeStreamer(ffin, DEFAULT_THREADS, true, true, -1, false);
		t.stopAndStart("Made streamer");
		
		//Start streamer
		ss.start();
		t.stopAndStart("Started streamer");

		//Create output writer with header from streamer
		FileFormat ffout=FileFormat.testOutput(out, FileFormat.SAM, null, true, false, false, true);
		SamWriter writer=SamWriter.makeWriter(ffout, DEFAULT_THREADS, ss.header, true);
		t.stopAndStart("Made writer");

		//Start writer
		writer.start();
		t.stopAndStart("Started writer");

		//Copy data
		for(ListNum<SamLine> list=ss.nextLines(); list!=null; list=ss.nextLines()){
//			System.err.println("Added list "+list.id);
			writer.addLines(list);
		}
		t.stopAndStart("Finished streamer");

		//Finish
		writer.poisonAndWait();
		t.stopAndStart("Closed writer");

		t2.stop();
		System.err.println();
		System.err.println(Tools.timeReadsBasesProcessed(t2, writer.readsWritten, writer.basesWritten, 8));
		System.err.println(Tools.readsBasesOut(t2.elapsed, writer.readsWritten, writer.basesWritten, 8));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Constructor. */
	SamLineWriter(FileFormat ffout_, int threads_,
		ArrayList<byte[]> header_, boolean useSharedHeader_){
		super(ffout_, threads_, header_, useSharedHeader_);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void start(){spawnThreads();}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn worker and writer threads. */
	void spawnThreads(){
		final int totalThreads=threads+1; //Workers plus writer
		
		alpt=new ArrayList<ProcessThread>(totalThreads);
		for(int i=0; i<totalThreads; i++){alpt.add(new ProcessThread(i));}
		if(verbose){System.err.println("Spawned threads.");}
		
		for(ProcessThread pt : alpt){pt.start();}
		if(verbose){System.err.println("Started threads.");}
	}
	
	private synchronized void writeHeader(){
		if(verbose) {System.err.println("Entered writeHeader");}
		if(headerWritten || supressHeader){return;}
		ArrayList<byte[]> headerLines;
		if(useSharedHeader){
			if(verbose) {System.err.println("Waiting on sharedHeader");}
			headerLines=SamReadInputStream.getSharedHeader(true);
		}else if(header!=null){
			headerLines=header;
		}else {
			// Generate header from Data.scaffoldNames
			headerLines=SamHeader.makeHeaderList(supressHeaderSequences, 
				ReadStreamWriter.MINCHROM, ReadStreamWriter.MAXCHROM);
		}
		if(verbose) {System.err.println("Got shared header");}
		if(headerLines==null) {
			System.err.println("Warning: Header was null, creating empty header");
			headerLines=new ArrayList<byte[]>();
		}
		ByteBuilder bb=new ByteBuilder();
		try{
			for(byte[] line : headerLines) {
				bb.append(line).nl();
				if(bb.length()>=16384) {
					outstream.write(bb.toBytes());
					bb.clear();
				}
			}
			if(bb.length()>=1) {
				outstream.write(bb.toBytes());
				bb.clear();
			}
		}catch(IOException e){
			throw new RuntimeException(e);
		}
		if(verbose) {System.err.println("Wrote header");}
		headerWritten=true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Processing thread - converts reads/lines to SAM text or writes output. */
	private class ProcessThread extends Thread {
		
		/** Constructor. */
		ProcessThread(final int tid_){
			tid=tid_;
		}
		
		/** Called by start(). */
		@Override
		public void run(){
			if(tid==0){
				if(verbose) {System.err.println("Writer tid "+tid+" started.");}
				writeOutput(); //Writer thread
				if(verbose) {System.err.println("Writer tid "+tid+" finished.");}
			}else{
				if(verbose) {System.err.println("Worker tid "+tid+" started.");}
				processJobs(); //Worker thread
				if(verbose) {System.err.println("Worker tid "+tid+" finished.");}
			}
			success=true;
		}
		
		/** Writer thread - outputs ordered blocks to disk. */
		void writeOutput(){
			if(verbose){System.err.println("tid "+tid+" started writeOutput.");}
			
			//Write header first
			writeHeader();
			
			//Write ordered data blocks
			SamWriterOutputJob job=outq.take();
			while(job!=null && !job.last()){
				if(verbose){System.err.println("tid "+tid+" writing block "+job.id);}
				
				try{
					outstream.write(job.bytes);
				}catch(Exception e){
					throw new RuntimeException("Error writing output", e);
				}
				
				if(job.last()){
					if(verbose){System.err.println("tid "+tid+" saw last job.");}
					break;
				}

				if(verbose){System.err.println("tid "+tid+" searching for job "+job.id+1);}
				job=outq.take();
				if(verbose){System.err.println("tid "+tid+" got job "+job.id);}
			}
			if(verbose){System.err.println("tid "+tid+" exited loops on "+job.id);}
			
			//Wait for other threads
			ThreadWaiter.waitForThreadsToFinish(alpt);
			//Accumulate per-thread statistics
			synchronized(SamLineWriter.this) {
				for(ProcessThread pt : alpt){
					synchronized(pt) {
						readsWritten+=pt.readsWrittenT;
						basesWritten+=pt.basesWrittenT;
						errorState=errorState || !pt.success;
					}
				}
			}
			if(verbose){System.err.println("tid "+tid+" exited loops on "+job.id);}
			
			//Close output stream and signal completion
			ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess());
			setFinished(true);
			
			if(verbose){System.err.println("tid "+tid+" finished writeOutput.");}
		}
		
		/** Worker thread - converts reads/lines to formatted bytes. */
		void processJobs(){
			if(verbose){System.err.println("tid "+tid+" started processJobs.");}

			final ByteBuilder bb=new ByteBuilder();
			SamWriterInputJob job=takeJob();
			for(; !job.poison(); job=takeJob()){
				if(verbose){System.err.println("tid "+tid+" processing job "+job.id());}
				
				//Handle LAST job
				if(job.last()){
					if(verbose){System.err.println("tid "+tid+" processing last job "+job.id());}
					//Create LAST output job with empty bytes
					SamWriterOutputJob outJob=new SamWriterOutputJob(job.id(), null, ListNum.LAST);
					outq.add(outJob);
					
					//Inject POISON for other workers
					job=new SamWriterInputJob(null, null, ListNum.POISON, job.id()+1);
					break;
				}
				
				//Convert to SamLines if needed
				ArrayList<SamLine> lines;
				if(job.lines!=null){
					lines=job.lines.list;
				}else{
					lines=new ArrayList<SamLine>(job.reads.list.size());
					lines=toSamLines(job.reads.list);
				}
				if(verbose){System.err.println("tid "+tid+" got "+lines.size()+" lines.");}
				
				//Format SamLines to bytes and count
				for(SamLine sl : lines){
					sl.toBytes(bb);
					bb.nl();
					readsWrittenT++;
					basesWrittenT+=sl.length();
				}
				if(verbose){System.err.println("tid "+tid+" made "+bb.length()+" bytes.");}
				
				//Create output job
				SamWriterOutputJob outJob=new SamWriterOutputJob(job.id(), bb.toBytes(), ListNum.NORMAL);
				outq.add(outJob);
				bb.clear();

				if(verbose){System.err.println("tid "+tid+" fetching job.");}
			}
			if(verbose){System.err.println("tid "+tid+" exited loop on job "+job.id());}
			
			//Re-inject poison for other workers
			assert(job.poison());
			putJob(job);
			if(verbose){System.err.println("tid "+tid+" finished processJobs.");}
		}
		
		/** Number of reads processed by this thread. */
		protected long readsWrittenT=0;
		/** Number of bases processed by this thread. */
		protected long basesWrittenT=0;
		/** True only if this thread completed successfully. */
		boolean success=false;
		/** Thread ID. */
		final int tid;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Thread list for accumulation. */
	private ArrayList<ProcessThread> alpt;
	
}