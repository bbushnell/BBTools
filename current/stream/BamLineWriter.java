package stream;

import java.io.IOException;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import stream.bam.SamToBamConverter;
import structures.ByteBuilder;
import structures.ListNum;
import template.ThreadWaiter;

/**
 * Writes BAM binary files with parallel conversion and ordered output.
 * 
 * Workers convert Read/SamLine objects to BAM binary in parallel.
 * OrderedQueueSystem ensures output blocks are written in order.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date October 25, 2025
 */
public class BamLineWriter extends SamWriter {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructor. */
	BamLineWriter(FileFormat ffout_, int threads_,
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

	@Override
	protected synchronized void writeHeader(){
		if(headerWritten){return;}

		ArrayList<byte[]> headerLines=getHeader();

		// Extract reference names for converter
		ArrayList<String> refNames=new ArrayList<String>();
		for(byte[] line : headerLines) {
			if(line.length > 3 && line[0] == '@' && line[1] == 'S' && line[2] == 'Q') {
				String lineStr = new String(line);
				String[] fields = lineStr.split("\\t");
				for(int i = 1; i < fields.length; i++) {
					if(fields[i].startsWith("SN:")) {
						refNames.add(fields[i].substring(3));
						break;
					}
				}
			}
		}
		// Write BAM header using helper
		stream.bam.BamWriterHelper writer = new stream.bam.BamWriterHelper(outstream);
		try{
			writer.writeHeaderFromLines(headerLines, supressHeader, supressHeaderSequences);
		}catch(IOException e){
			throw new RuntimeException(e);
		}

		// Create converter for workers
		sharedConverter=new SamToBamConverter(refNames.toArray(new String[0]));
		headerWritten=true;
		this.notifyAll();
	}
	
	private SamToBamConverter getConverter() {
		synchronized(this) {
			while(sharedConverter==null) {
				try{this.wait();}
				catch(InterruptedException e){e.printStackTrace();}
			}
		}
		return sharedConverter;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Processing thread - converts reads/lines to BAM binary or writes output. */
	private class ProcessThread extends Thread {

		/** Constructor. */
		ProcessThread(final int tid_){
			tid=tid_;
			if(verbose) {System.err.println("tid "+tid+" created.");}
		}

		/** Called by start(). */
		@Override
		public void run(){
			synchronized(this) {
				if(tid==0){
					writeOutput(); //Writer thread
					if(verbose) {System.err.println("Consumer "+tid+" finished.");}
				}else{
					processJobs(); //Worker thread
					if(verbose) {System.err.println("Worker "+tid+" finished.");}
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
			synchronized(BamLineWriter.this) {
				for(ProcessThread pt : alpt){
					if(pt!=this) {
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
			final SamToBamConverter converter=getConverter();
			while(!job.poison()){
				//Convert to SamLines if needed
				ArrayList<SamLine> lines;
				if(job.lines!=null){
					lines=job.lines.list;
				}else{
					lines=toSamLines(job.reads.list);
				}

				//Format SamLines to BAM bytes and count
				for(SamLine sl : lines){
					try{
//						byte[] bamRecord = converter.convertAlignment(sl);
//						// Write block_size followed by the record
//						bb.appendUint32(bamRecord.length);
//						bb.append(bamRecord);
						converter.appendAlignment(sl, bb);
						readsWrittenT++;
						basesWrittenT+=sl.length();
					}catch(Exception e){
						throw new RuntimeException("Error converting to BAM", e);
					}
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
	/** Converter for SamLine to BAM binary. */
	private volatile SamToBamConverter sharedConverter;

}