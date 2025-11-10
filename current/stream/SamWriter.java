package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
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
 * @author Brian Bushnell, Isla
 * @date October 25, 2025
 */
public class SamWriter implements Writer {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructor. */
	public SamWriter(FileFormat ffout_, int threads_,
		ArrayList<byte[]> header_, boolean useSharedHeader_){
		ffout=ffout_;
		fname=ffout.name();
		threads=Tools.mid(1, threads_<1 ? DEFAULT_THREADS : threads_, Shared.threads());
		header=header_;
		useSharedHeader=useSharedHeader_;
		supressHeader=(ReadStreamWriter.NO_HEADER || (ffout.append() && ffout.exists()));
		supressHeaderSequences=(ReadStreamWriter.NO_HEADER_SEQUENCES || supressHeader);

		//Create prototype jobs for OrderedQueueSystem
		SamWriterInputJob inputProto=new SamWriterInputJob(null, null, ListNum.PROTO, -1);
		SamWriterOutputJob outputProto=new SamWriterOutputJob(-1, null, ListNum.PROTO);

		oqs=new OrderedQueueSystem<SamWriterInputJob, SamWriterOutputJob>(
			threads, ffout.ordered(), inputProto, outputProto);
		
		if(ffout.bam()){
			outstream=ReadWrite.getBamOutputStream(fname, ffout.append());
		}else {
			outstream=ReadWrite.getOutputStream(fname, ffout.append(), true, ffout.allowSubprocess());
		}
		if(verbose) {System.err.println("outstream="+outstream.getClass());}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public void start(){spawnThreads();}
	
	public final void add(ArrayList<Read> list, long id) {addReads(new ListNum<Read>(list, id));}

	/** Add reads for writing (will be converted to SamLines). */
	public final void addReads(ListNum<Read> reads){
		if(reads==null){return;}
		SamWriterInputJob job=new SamWriterInputJob(reads, null, ListNum.NORMAL, reads.id);
		oqs.addInput(job);
	}

	/** Add already-formatted SamLines for writing. */
	public final void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		SamWriterInputJob job=new SamWriterInputJob(null, lines, ListNum.NORMAL, lines.id);
		oqs.addInput(job);
	}

	/** Signal end of input. */
	public final void poison(){
		oqs.poison();
	}

	/** Wait for all writes to complete. */
	public final boolean waitForFinish(){
		oqs.waitForFinish();
		return errorState();
	}

	/** Convenience method - poison and wait. */
	public final boolean poisonAndWait(){
		poison();
		return waitForFinish();
	}

	@Override
	public void close(){
		//TODO: Interrupt and terminate
	}

	@Override
	public long readsWritten(){return readsWritten;}

	@Override
	public long basesWritten(){return basesWritten;}

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
	
	public static ArrayList<SamLine> toSamLines(ArrayList<Read> reads) {
		ArrayList<SamLine> samLines=new ArrayList<SamLine>();

		for(final Read r1 : reads){
			Read r2=(r1==null ? null : r1.mate);

			SamLine sl1=(r1==null ? null : (ReadStreamWriter.USE_ATTACHED_SAMLINE 
				&& r1.samline!=null ? r1.samline : new SamLine(r1, 0)));
			SamLine sl2=(r2==null ? null : (ReadStreamWriter.USE_ATTACHED_SAMLINE 
				&& r2.samline!=null ? r2.samline : new SamLine(r2, 1)));

			if(!SamLine.KEEP_NAMES && sl1!=null && sl2!=null && ((sl2.qname==null) || 
				!sl2.qname.equals(sl1.qname))){
				sl2.qname=sl1.qname;
			}
			assert(sl1!=null) : r1;
			addSamLine(r1, sl1, samLines);
			addSamLine(r2, sl2, samLines);
		}
		return samLines;
	}

	private static void addSamLine(Read r, SamLine primary, ArrayList<SamLine> samLines) {
		if(r==null || primary==null) {return;}
		
		assert(!ReadStreamWriter.ASSERT_CIGAR || !r.mapped() || primary.cigar!=null) : r;
		samLines.add(primary);

		// Handle secondary alignments
		ArrayList<SiteScore> list=r.sites;
		if(ReadStreamWriter.OUTPUT_SAM_SECONDARY_ALIGNMENTS && list!=null && list.size()>1){
			final Read clone=r.clone();
			for(int i=1; i<list.size(); i++){
				SiteScore ss=list.get(i);
				clone.match=null;
				clone.setFromSite(ss);
				clone.setSecondary(true);
				SamLine secondary=new SamLine(clone, r.pairnum());
				assert(!secondary.primary());
				assert(!ReadStreamWriter.USE_ATTACHED_SAMLINE || secondary.cigar!=null) : r;
				samLines.add(secondary);
			}
		}
	}
	
	ArrayList<byte[]> getHeader(){
		if(verbose) {System.err.println("Fetching header: "+useSharedHeader+","+(header!=null));}
		ArrayList<byte[]> headerLines;
		if(useSharedHeader){
			headerLines=SamReadInputStream.getSharedHeader(true);
		}else if(header!=null){
			headerLines=header;
		}else {
			headerLines=SamHeader.makeHeaderList(supressHeaderSequences, 
				ReadStreamWriter.MINCHROM, ReadStreamWriter.MAXCHROM);
		}
		if(headerLines==null) {
			System.err.println("Warning: Header was null, creating empty header");
			headerLines=new ArrayList<byte[]>();
		}
		if(verbose) {System.err.println("Fetched header: "+(headerLines==null ? "null" : headerLines.size()));}
		return headerLines;
	}

	protected synchronized void writeHeader(){
		if(headerWritten || supressHeader){return;}
		ArrayList<byte[]> headerLines=getHeader();
		
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
		headerWritten=true;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Input job for workers - contains either reads or lines. */
	static class SamWriterInputJob implements HasID {

		public SamWriterInputJob(ListNum<Read> reads_, ListNum<SamLine> lines_, int type_, long id_){
			reads=reads_;
			lines=lines_;
			type=type_;
			id=id_;
			assert(type!=ListNum.NORMAL || ((reads==null) != (lines==null)));
			assert(type==ListNum.NORMAL || type==ListNum.LAST || type==ListNum.POISON || type==ListNum.PROTO);
			assert(reads==null || id==reads.id);
			assert(lines==null || id==lines.id);
		}

		@Override
		public long id(){return id;}

		@Override
		public boolean poison(){return type==ListNum.POISON;}

		@Override
		public boolean last(){return type==ListNum.LAST;}

		@Override
		public SamWriterInputJob makePoison(long id_) {
			return new SamWriterInputJob(null, null, ListNum.POISON, id_);
		}

		@Override
		public SamWriterInputJob makeLast(long id_){
			return new SamWriterInputJob(null, null, ListNum.LAST, id_);
		}

		public final ListNum<Read> reads;
		public final ListNum<SamLine> lines;
		public final int type;
		public final long id;
	}

	/** Output job for writer - ordered formatted bytes. */
	static class SamWriterOutputJob implements HasID {

		public SamWriterOutputJob(long id_, byte[] bytes_, int type_){
			id=id_;
			bytes=bytes_;
			type=type_;
			assert((type==ListNum.NORMAL) == (bytes!=null));
			assert(type==ListNum.NORMAL || type==ListNum.LAST || type==ListNum.PROTO || type==ListNum.POISON);
		}

		@Override
		public long id(){return id;}

		@Override
		public boolean poison(){return false;}

		@Override
		public boolean last(){return type==ListNum.LAST;}

		@Override
		public SamWriterOutputJob makePoison(long id_) {
			return new SamWriterOutputJob(id_, null, ListNum.POISON);
		}

		@Override
		public SamWriterOutputJob makeLast(long id_){
			return new SamWriterOutputJob(id_, null, ListNum.LAST);
		}

		public final long id;
		public final byte[] bytes;
		public final int type;
	}

	/** Processing thread - converts reads/lines to SAM text or writes output. */
	private class ProcessThread extends Thread {

		/** Constructor. */
		ProcessThread(final int tid_){
			tid=tid_;
			setName("SamLineWriter-"+(tid==0 ? "Output" : "Worker-"+tid));
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
			synchronized(SamWriter.this) {
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
	/*----------------     Getters and Setters      ----------------*/
	/*--------------------------------------------------------------*/
	
	synchronized void setErrorState(boolean b){
		errorState|=b;
	}
	
	@Override
	public synchronized boolean errorState() {return errorState;}
	
	@Override
	public boolean finishedSuccessfully() {return !errorState && oqs.finished();}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Output file name. */
	final String fname;
	/** Output file format. */
	final FileFormat ffout;
	/** Number of worker threads. */
	final int threads;
	/** True if header should be pulled from shared input. */
	final boolean useSharedHeader;
	/** True if header should be skipped. */
	final boolean supressHeader;
	/** True if header sequences should be written. */
	final boolean supressHeaderSequences;
	/** True after header is written */
	boolean headerWritten=false;
	/** Header lines to write. */
	final ArrayList<byte[]> header;
	/** Ordered queue system for coordination. */
	final OrderedQueueSystem<SamWriterInputJob, SamWriterOutputJob> oqs;
	/** Output stream. */
	final OutputStream outstream;
	
	/** Thread list for accumulation. */
	private ArrayList<ProcessThread> alpt;
	
	/** Total reads written. */
	public long readsWritten=0;
	/** Total bases written. */
	public long basesWritten=0;
	/** Were any errors encountered */
	private boolean errorState=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	public static int DEFAULT_THREADS=6;

	public static final boolean verbose=false;
	
	/** Print status messages to this output stream */
	protected PrintStream outstream2=System.err;

}