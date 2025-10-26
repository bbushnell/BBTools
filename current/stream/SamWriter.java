package stream;

import java.io.OutputStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ListNum;

/**
 * Multithreaded SAM/BAM writer using parallel conversion and ordered output.
 * 
 * Abstract base class for writing SAM/BAM files with parallel Read-to-SamLine
 * conversion and ordered output via JobQueue. Workers convert reads to formatted
 * bytes in parallel while a single writer thread outputs ordered blocks.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date October 25, 2025
 */
public abstract class SamWriter {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Factory method to create appropriate writer type. */
	public static SamWriter makeWriter(String out, int threads, 
			ArrayList<byte[]> header, boolean useSharedHeader, boolean ordered){
		FileFormat ffout=FileFormat.testOutput(out, FileFormat.SAM, null, true, false, false, ordered);
		return new SamLineWriter(ffout, threads, header, useSharedHeader);
	}
	
	/** Factory method to create appropriate writer type. */
	public static SamWriter makeWriter(FileFormat ffout, int threads,
			ArrayList<byte[]> header, boolean useSharedHeader){
		if(ffout.bam() && ReadWrite.nativeBamOut()) {
			throw new RuntimeException("TODO");
		}else {
			return new SamLineWriter(ffout, threads, header, useSharedHeader);
		}
	}
	
	/** Constructor. */
	public SamWriter(FileFormat ffout_, int threads_, 
			ArrayList<byte[]> header_, boolean useSharedHeader_){
		ffout=ffout_;
		fname=ffout.name();
		ordered=ffout.ordered();
		threads=threads_;
		header=header_;
		useSharedHeader=useSharedHeader_;
		supressHeader=(ReadStreamWriter.NO_HEADER || (ffout.append() && ffout.exists()));
		supressHeaderSequences=(ReadStreamWriter.NO_HEADER_SEQUENCES || supressHeader);
		
		int queueSize=1+2*threads;
		inq=new ArrayBlockingQueue<SamWriterInputJob>(queueSize);
		outq=new JobQueue<SamWriterOutputJob>(queueSize, ordered, true, 0);
		
		try{
			outstream=ReadWrite.getOutputStream(fname, ffout.append(), true, ffout.allowSubprocess());
		}catch(Exception e){
			throw new RuntimeException("Error opening output stream for "+fname, e);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Start worker threads. */
	public abstract void start();
	
	/** Add reads for writing (will be converted to SamLines). */
	public final void addReads(ListNum<Read> reads){
		if(reads==null){return;}
		assert(!reads.poison() && !reads.last()) : "Use poison() to terminate";
		SamWriterInputJob job=new SamWriterInputJob(reads, null, ListNum.NORMAL, reads.id);
		putJob(job);
	}
	
	/** Add already-formatted SamLines for writing. */
	public final void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		assert(!lines.poison() && !lines.last()) : "Use poison() to terminate";
		SamWriterInputJob job=new SamWriterInputJob(null, lines, ListNum.NORMAL, lines.id);
		putJob(job);
	}
	
	/** Signal end of input - creates LAST job and returns immediately. */
	public final synchronized void poison(){
		if(verbose) {System.err.println("Producer called poison().");}
		if(lastSeen){return;} //Prevent multiple calls
		SamWriterInputJob job=new SamWriterInputJob(null, null, ListNum.LAST, maxSeenId+1);
		if(verbose) {System.err.println("Producer adding LAST.");}
		putJob(job);
		if(verbose) {System.err.println("Producer added LAST.");}
		lastSeen=true;
	}
	
	/** Wait for all writes to complete and all threads to terminate. */
	public final synchronized void waitForFinish(){
		if(verbose) {System.err.println("Producer waiting for finish.");}
		while(!finished){
			try{this.wait();}
			catch(InterruptedException e){e.printStackTrace();}
		}
		if(verbose) {System.err.println("Producer wait for finish ended.");}
	}
	
	/** Convenience method - poison and wait. */
	public final void poisonAndWait(){
		poison();
		waitForFinish();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Put job in input queue. */
	final void putJob(SamWriterInputJob job){
		if(verbose) {System.err.println("addJob "+job.id);}
		synchronized(this) {
			assert(!lastSeen || job.poison());
			assert(!job.last() || job.id>maxSeenId);
			assert(!job.poison() || job.id>=maxSeenId);
			maxSeenId=Math.max(job.id, maxSeenId);
		}
		while(job!=null){
			try{
				inq.put(job);
				job=null;
			}catch(InterruptedException e){
				e.printStackTrace();
			}
		}
		if(verbose) {System.err.println("inq now size "+inq.size());}
	}
	
	/** Take job from input queue. */
	final SamWriterInputJob takeJob(){
		if(verbose) {System.err.println("takeJob "+inq.size());}
		SamWriterInputJob job=null;
		while(job==null){
			try{
				job=inq.take();
				assert(job!=null);
			}catch(InterruptedException e){
				e.printStackTrace();
			}
		}
		if(verbose) {System.err.println("took job "+job.id()+": "+inq.size());}
		return job;
	}
	
	/** Set finished flag (synchronized). */
	final synchronized void setFinished(boolean b){
		finished=b;
		this.notifyAll();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Helper Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
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
			assert(type!=ListNum.NORMAL || ((reads==null) != (lines==null))) : 
				"Exactly one must be non-null for non-LAST jobs: "+type+", "+(reads==null);
			assert(type==ListNum.NORMAL || type==ListNum.LAST || type==ListNum.POISON);
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
			assert(type==ListNum.NORMAL || type==ListNum.LAST);
		}
		
		@Override
		public long id(){return id;}
		
		@Override
		public boolean poison(){return false;} //Never poison in output queue
		
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
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Output file name. */
	final String fname;
	/** Output file format. */
	final FileFormat ffout;
	/** True if output should be ordered. */
	final boolean ordered;
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
	/** Input queue for worker threads. */
	final ArrayBlockingQueue<SamWriterInputJob> inq;
	/** Output queue for ordered writing. */
	final JobQueue<SamWriterOutputJob> outq;
	/** Output stream. */
	final OutputStream outstream;
	/** True when writing is complete. */
	private volatile boolean finished=false;
	/** True when poison has been called. */
	private volatile boolean lastSeen=false;
	/** Total reads written. */
	public long readsWritten=0;
	/** Total bases written. */
	public long basesWritten=0;
	/** Were any errors encountered */
	public boolean errorState=false;
	private long maxSeenId=-1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static int DEFAULT_THREADS=4; //TODO: Test, likely input limited at >3.
	
	public static final boolean verbose=false;
	
}