package stream;

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
 * OrderedQueueSystem ensures output blocks are written in order.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date October 25, 2025
 */
public class SamLineWriter extends SamWriter {

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

		for(ProcessThread pt : alpt){pt.start();}
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
			synchronized(this) {
				if(tid==0){
					writeOutput(); //Writer thread
				}else{
					processJobs(); //Worker thread
				}
				success=true;
			}
		}

		/** Writer thread - outputs ordered blocks to disk. */
		void writeOutput(){
			//Write header first
			writeHeader();

			//Write ordered data blocks
			SamWriterOutputJob job=oqs.getOutput();
			while(job!=null && !job.last()){
				try{
					outstream.write(job.bytes);
				}catch(Exception e){
					throw new RuntimeException("Error writing output", e);
				}

				job=oqs.getOutput();
			}

			//Wait for other threads and accumulate statistics
			ThreadWaiter.waitForThreadsToFinish(alpt);
			synchronized(SamLineWriter.this) {
				for(ProcessThread pt : alpt){
					if(pt!=this) {//This is not successful yet!
						synchronized(pt) {
							readsWritten+=pt.readsWrittenT;
							basesWritten+=pt.basesWrittenT;
							setErrorState(!pt.success);
						}
					}
				}
			}
			
			if(verbose) {System.err.println("Consumer finished accumulating.");}
			//Close output stream and signal completion
			ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess());
			if(verbose) {System.err.println("Consumer finished writing.");}
			oqs.setFinished();
			if(verbose) {System.err.println("Consumer set oqs finished.");}
		}

		/** Worker thread - converts reads/lines to formatted bytes. */
		void processJobs(){
			final ByteBuilder bb=new ByteBuilder();

			SamWriterInputJob job=oqs.getInput();
			while(!job.poison()){
				//Convert to SamLines if needed
				ArrayList<SamLine> lines;
				if(job.lines!=null){
					lines=job.lines.list;
				}else{
					lines=toSamLines(job.reads.list);
				}

				//Format SamLines to bytes and count
				for(SamLine sl : lines){
					sl.toBytes(bb);
					bb.nl();
					readsWrittenT++;
					basesWrittenT+=sl.length();
				}

				//Create output job
				SamWriterOutputJob outJob=new SamWriterOutputJob(job.id(), bb.toBytes(), ListNum.NORMAL);
				oqs.addOutput(outJob);
				bb.clear();

				job=oqs.getInput();
			}

			//Re-inject poison for other workers
			oqs.addInput(job);
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