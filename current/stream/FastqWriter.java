package stream;

import java.io.OutputStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;
import template.ThreadWaiter;

/**
 * Writes FASTQ files with parallel formatting and ordered output.
 * 
 * Workers convert Read objects to FASTQ text in parallel.
 * OrderedQueueSystem ensures output blocks are written in order.
 * 
 * @author Isla
 * @date October 30, 2025
 */
public class FastqWriter implements Writer {
	
	public static void main(String[] args) {
		Timer t=new Timer();
		String in=args[0];
		String out=(args.length<2 || args[1].equalsIgnoreCase("null") ? null : args[1]);
		int threads=DEFAULT_THREADS;
		if(args.length>2) {threads=Integer.parseInt(args[2]);}
		if(args.length>3) {Shared.SIMD=true;}
		if(args.length>4) {
			ReadWrite.ALLOW_NATIVE_BGZF=ReadWrite.PREFER_NATIVE_BGZF_IN=
				ReadWrite.PREFER_NATIVE_BGZF_OUT=Parse.parseBoolean(args[4]);
		}

//		ByteFile.FORCE_MODE_BF4=false;
		FileFormat ffin=FileFormat.testInput(in, FileFormat.FASTQ, null, true, true);
		FileFormat ffout=FileFormat.testOutput(out, FileFormat.FASTQ, null, true, true, false, true);
		
		SamLine.SET_FROM_OK=ffin.samOrBam();
		ReadStreamByteWriter.USE_ATTACHED_SAMLINE=(ffout!=null && ffout.samOrBam() && ffin.samOrBam());

		final boolean inputReads=(!ffin.samOrBam());
		final boolean outputReads=(ffout!=null && !ffout.samOrBam());
		
		Streamer st=StreamerFactory.makeStreamer(ffin, 0, true, -1, true, outputReads);
		Writer fw=WriterFactory.makeWriter(ffout, true, true, threads, null, ffin.samOrBam());
//		assert(false) : "\n"+ffin+"\n"+ffout+"\n"+st.getClass()+"\n"+fw.getClass();
		process(st, fw, t, inputReads || outputReads);
	}
	
	private static void process(Streamer st, Writer fw, Timer t, boolean readMode) {
		st.start();
		if(fw!=null) {fw.start();}
		long reads=0, bases=0;
		if(readMode) {
			for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()) {
				for(Read r : ln) {
					reads+=r.pairCount();
					bases+=r.pairLength();
				}
				if(fw!=null) {fw.addReads(ln);}
			}
		}else {
			for(ListNum<SamLine> ln=st.nextLines(); ln!=null; ln=st.nextLines()) {
				for(SamLine sl : ln) {
					reads++;
					bases+=sl.length();
				}
				if(fw!=null) {fw.addLines(ln);}
			}
		}
		if(fw!=null) {
			fw.poisonAndWait();
			assert(reads==fw.readsWritten());
			assert(bases==fw.basesWritten());
			reads=fw.readsWritten();
			bases=fw.basesWritten();
		}
		t.stop();
		System.err.println(Tools.timeReadsBasesProcessed(t, reads, bases, 8));
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Constructor. */
	public FastqWriter(String out_, int threads_, boolean writeR1_, boolean writeR2_, boolean overwrite){
		this(FileFormat.testOutput(out_, FileFormat.FASTQ, null, true, overwrite, false, true), 
			threads_, writeR1_, writeR2_);
	}
	
	/** Constructor. */
	public FastqWriter(FileFormat ffout_, int threads_, boolean writeR1_, boolean writeR2_){
		ffout=ffout_;
		fname=ffout_.name();
		threads=Tools.mid(1, threads_<1 ? DEFAULT_THREADS : threads_, Shared.threads());
		writeR1=writeR1_;
		writeR2=writeR2_;
		format=(ffout.format()==UNKNOWN ? FASTQ : ffout.format());
		assert(format==FASTQ || format==FASTA || format==HEADER) : ffout;
		
		assert(writeR1 || writeR2) : "Must write at least one mate";
		
		// Create OQS
		FastqWriterInputJob inputProto=new FastqWriterInputJob(null, 0, ListNum.PROTO);
		FastqWriterOutputJob outputProto=new FastqWriterOutputJob(0, null, ListNum.PROTO);
		oqs=new OrderedQueueSystem<FastqWriterInputJob, FastqWriterOutputJob>(
			threads, true, inputProto, outputProto);
		
		// Open output stream
		outstream=ReadWrite.getOutputStream(fname, false, true, false);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void start(){
		spawnThreads();
	}
	
	@Override
	public void close(){
		// TODO: Unimplemented
	}
	
	@Override
	public long readsWritten(){
		return readsWritten;
	}
	
	@Override
	public long basesWritten(){
		return basesWritten;
	}
	
	public final void add(ArrayList<Read> list, long id) {addReads(new ListNum<Read>(list, id));}
	
	@Override
	public void addReads(ListNum<Read> reads){
		FastqWriterInputJob job=new FastqWriterInputJob(reads, reads.id(), ListNum.NORMAL);
		oqs.addInput(job);
	}
	
	@Override
	public void addLines(ListNum<SamLine> lines){//Should be fairly fast
		ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		for(SamLine sl : lines) {
			reads.add(new Read(sl.seq, sl.qual, sl.qname, -1, false));
		}
		addReads(new ListNum<Read>(reads, lines.id));
	}
	
	@Override
	public void poison(){
		oqs.poison();
	}
	
	@Override
	public boolean waitForFinish(){
		oqs.waitForFinish();
		return errorState;
	}
	
	@Override
	public boolean poisonAndWait(){
		poison();
		return waitForFinish();
	}
	
	@Override
	public boolean errorState(){return errorState;}
	
	@Override
	public boolean finishedSuccessfully() {return !errorState && oqs.finished();}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Spawn worker and writer threads. */
	void spawnThreads(){
		final int totalThreads=threads+1; // Workers plus writer
		
		alpt=new ArrayList<ProcessThread>(totalThreads);
		for(int i=0; i<totalThreads; i++){
			alpt.add(new ProcessThread(i));
		}
		
		for(ProcessThread pt : alpt){
			pt.start();
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Input job for OQS */
	private static class FastqWriterInputJob implements HasID {//TODO: This class should just be a ln
		FastqWriterInputJob(ListNum<Read> reads_, long id_, int type_){
			reads=reads_;
			id=id_;
			type=type_;
		}
		
		@Override public long id(){return id;}
		@Override public boolean poison(){return type==ListNum.POISON;}
		@Override public boolean last(){return type==ListNum.LAST;}
		@Override public FastqWriterInputJob makePoison(long id){
			return new FastqWriterInputJob(null, id, ListNum.POISON);
		}
		@Override public FastqWriterInputJob makeLast(long id){
			return new FastqWriterInputJob(null, id, ListNum.LAST);
		}
		
		final ListNum<Read> reads;
		final long id;
		final int type;
	}
	
	/** Output job for OQS */
	private static class FastqWriterOutputJob implements HasID {
		FastqWriterOutputJob(long id_, byte[] bytes_, int type_){
			id=id_;
			bytes=bytes_;
			type=type_;
		}
		
		@Override public long id(){return id;}
		@Override public boolean poison(){return type==ListNum.POISON;}
		@Override public boolean last(){return type==ListNum.LAST;}
		@Override public FastqWriterOutputJob makePoison(long id){
			return new FastqWriterOutputJob(id, null, ListNum.POISON);
		}
		@Override public FastqWriterOutputJob makeLast(long id){
			return new FastqWriterOutputJob(id, null, ListNum.LAST);
		}
		
		final long id;
		final byte[] bytes;
		final int type;
	}
	
	/** Processing thread - converts reads to FASTQ text or writes output. */
	private class ProcessThread extends Thread {
		
		/** Constructor. */
		ProcessThread(final int tid_){
			tid=tid_;
			setName("FastqWriter-"+(tid==0 ? "Output" : "Worker-"+tid));
		}
		
		/** Called by start(). */
		@Override
		public void run(){
			synchronized(this){
				if(tid==0){
					writeOutput(); // Writer thread
				}else{
					processJobs(); // Worker thread
				}
				success=true;
			}
		}
		
		/** Writer thread - outputs ordered blocks to disk. */
		void writeOutput(){
			// Write ordered data blocks
			FastqWriterOutputJob job=oqs.getOutput();
			while(job!=null && !job.last()){
				try{
					outstream.write(job.bytes);
				}catch(Exception e){
					throw new RuntimeException("Error writing output", e);
				}
				
				job=oqs.getOutput();
			}
			
			// Wait for other threads and accumulate statistics
			ThreadWaiter.waitForThreadsToFinish(alpt);
			synchronized(FastqWriter.this){
				for(ProcessThread pt : alpt){
					if(pt!=this){ // This thread not successful yet!
						synchronized(pt){
							readsWritten+=pt.readsWrittenT;
							basesWritten+=pt.basesWrittenT;
							errorState|=!pt.success;
						}
					}
				}
			}
			
			// Close output stream and signal completion
			ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess());
			oqs.setFinished();
		}
		
		/** Worker thread - converts reads to formatted bytes. */
		void processJobs(){
			final ByteBuilder bb=new ByteBuilder();
			
			FastqWriterInputJob job=oqs.getInput();
			while(!job.poison()){
				ArrayList<Read> reads=job.reads.list;
				
				// Format reads
				if(format==FASTQ) {
					writeFastq(reads, bb);
				}else if(format==FASTA) {
					writeFasta(reads, bb);
				}else if(format==HEADER) {
					writeHeader(reads, bb);
				}else {
					throw new RuntimeException("Bad format: "+format);
				}
				
				// Create output job
				FastqWriterOutputJob outJob=new FastqWriterOutputJob(job.id(), bb.toBytes(), ListNum.NORMAL);
				oqs.addOutput(outJob);
				bb.clear();
				
				job=oqs.getInput();
			}
			
			// Re-inject poison for other workers
			oqs.addInput(job);
		}
		
		private void writeFastq(ArrayList<Read> reads, ByteBuilder bb) {
			for(Read r : reads){
				final Read r1=(r.pairnum()==0 ? r : null);
				final Read r2=(r.pairnum()==1 ? r : r.mate);
				if(writeR1 && r1!=null){
					r1.toFastq(bb);
					bb.nl();
					readsWrittenT++;
					basesWrittenT+=r1.length();
				}
				if(writeR2 && r2!=null){
					r2.toFastq(bb);
					bb.nl();
					readsWrittenT++;
					basesWrittenT+=r2.length();
				}
			}
		}
		
		private void writeFasta(ArrayList<Read> reads, ByteBuilder bb) {
			for(Read r : reads){
				final Read r1=(r.pairnum()==0 ? r : null);
				final Read r2=(r.pairnum()==1 ? r : r.mate);
				if(writeR1 && r1!=null){
					r1.toFasta(bb);
					bb.nl();
					readsWrittenT++;
					basesWrittenT+=r1.length();
				}
				if(writeR2 && r2!=null){
					r2.toFasta(bb);
					bb.nl();
					readsWrittenT++;
					basesWrittenT+=r2.length();
				}
			}
		}
		
		private void writeHeader(ArrayList<Read> reads, ByteBuilder bb) {
			for(Read r : reads){
				final Read r1=(r.pairnum()==0 ? r : null);
				final Read r2=(r.pairnum()==1 ? r : r.mate);
				if(writeR1 && r1!=null){
					bb.appendln(r1.id);
					readsWrittenT++;
				}
				if(writeR2 && r2!=null){
					bb.appendln(r2.id);
					readsWrittenT++;
				}
			}
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
	
	/** Output file path */
	public final String fname;
	/** Output file format */
	final FileFormat ffout;
	/** Output file format as an int */
	public final int format;
	/** Output stream */
	OutputStream outstream;
	/** OQS for coordinating workers and writer */
	final OrderedQueueSystem<FastqWriterInputJob, FastqWriterOutputJob> oqs;
	/** Number of worker threads */
	final int threads;
	/** Write R1 reads (pairnum==0) */
	final boolean writeR1;
	/** Write R2 reads (pairnum==1 or mate) */
	final boolean writeR2;
	/** Thread list for accumulation */
	private ArrayList<ProcessThread> alpt;
	/** Number of reads written */
	protected long readsWritten=0;
	/** Number of bases written */
	protected long basesWritten=0;
	/** True if an error was encountered */
	public boolean errorState=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private static final int FASTQ=FileFormat.FASTQ;
	private static final int FASTA=FileFormat.FASTA;
	private static final int HEADER=FileFormat.HEADER;
	private static final int UNKNOWN=FileFormat.UNKNOWN;
	
	public static int DEFAULT_THREADS=3;
	public static final boolean verbose=false;
	
}