package stream;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Multithreaded SAM/BAM writer using parallel conversion and ordered output.
 * 
 * Abstract base class for writing SAM/BAM files with parallel Read-to-SamLine
 * conversion and ordered output via OrderedQueueSystem. Workers convert reads 
 * to formatted bytes in parallel while a single writer thread outputs ordered blocks.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date October 25, 2025
 */
public abstract class SamWriter {

	public static void main(String[] args){
		Timer t=new Timer(), t2=new Timer();
		String in=args[0];
		String out=args[1];
		ReadWrite.PREFER_NATIVE_BAM_OUT=ReadWrite.ALLOW_NATIVE_BAM_OUT=args.length>2;
		ReadWrite.PREFER_NATIVE_BGZF_OUT=ReadWrite.ALLOW_NATIVE_BGZF=(args.length>3);

		//Create input streamer
		FileFormat ffin=FileFormat.testInput(in, FileFormat.SAM, null, true, false);
		SamStreamer ss=SamStreamer.makeStreamer(ffin, DEFAULT_THREADS, true, true, -1, false);
		t.stopAndStart("Made streamer");

		//Start streamer
		ss.start();
		t.stopAndStart("Started streamer");

		//Create output writer with header from streamer
		FileFormat ffout=FileFormat.testOutput(out, FileFormat.SAM, null, true, false, false, true);
		SamWriter writer=SamWriter.makeWriter(ffout, DEFAULT_THREADS, null, true);
		t.stopAndStart("Made writer - "+writer.getClass());

		//Start writer
		writer.start();
		t.stopAndStart("Started writer");

		//Copy data
		for(ListNum<SamLine> list=ss.nextLines(); list!=null; list=ss.nextLines()){
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
//		if(verbose) {
//			System.err.println("Main finished.");
//			try {
//				// Wait for 1000 milliseconds (1 second)
//				Thread.sleep(1000);
//			} catch (InterruptedException e) {
//				// Handle the interruption if the thread is interrupted while sleeping
//				e.printStackTrace();
//			}
//			Set<Thread> threads = Thread.getAllStackTraces().keySet();
//			for(Thread th : threads){
//				System.err.println("Thread: " + th.getName() + 
//					" daemon=" + th.isDaemon() + " state=" + th.getState());
//			}
//		}
	}

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
			return new BamLineWriter(ffout, threads, header, useSharedHeader);
		}else {
			return new SamLineWriter(ffout, threads, header, useSharedHeader);
		}
	}

	/** Constructor. */
	public SamWriter(FileFormat ffout_, int threads_, 
		ArrayList<byte[]> header_, boolean useSharedHeader_){
		ffout=ffout_;
		fname=ffout.name();
		threads=threads_;
		header=header_;
		useSharedHeader=useSharedHeader_;
		supressHeader=(ReadStreamWriter.NO_HEADER || (ffout.append() && ffout.exists()));
		supressHeaderSequences=(ReadStreamWriter.NO_HEADER_SEQUENCES || supressHeader);
		final boolean nativeBam=(ffout.bam() && ReadWrite.nativeBamOut());
		assert(nativeBam==(getClass()==BamLineWriter.class));

		int queueSize=1+2*threads;

		//Create prototype jobs for OrderedQueueSystem
		SamWriterInputJob inputProto=new SamWriterInputJob(null, null, ListNum.PROTO, -1);
		SamWriterOutputJob outputProto=new SamWriterOutputJob(-1, null, ListNum.PROTO);

		oqs=new OrderedQueueSystem<SamWriterInputJob, SamWriterOutputJob>(
			queueSize, ffout.ordered(), inputProto, outputProto);
		
		if(nativeBam) {
			outstream=ReadWrite.getBgzipStream(fname, false);
//			try{//To see raw bam
//				outstream=new FileOutputStream(fname);
//			}catch(FileNotFoundException e){
//				throw new RuntimeException(e);
//			}
		}else if(ffout.bam()){
			outstream=ReadWrite.getBamOutputStream(fname, ffout.append());
		}else {
			outstream=ReadWrite.getOutputStream(fname, ffout.append(), true, ffout.allowSubprocess());
		}
		if(verbose) {System.err.println("outstream="+outstream.getClass());}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Start worker threads. */
	public abstract void start();

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
	
	ArrayList<byte[]> getHeader(){
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
			assert(type!=ListNum.NORMAL || ((reads==null) != (lines==null))) : 
				"Exactly one must be non-null for NORMAL jobs";
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
			assert(type==ListNum.NORMAL || type==ListNum.LAST || type==ListNum.PROTO);
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

	/*--------------------------------------------------------------*/
	/*----------------     Getters and Setters      ----------------*/
	/*--------------------------------------------------------------*/
	
	synchronized void setErrorState(boolean b){
		errorState|=b;//Can never be unset
	}
	public synchronized boolean errorState() {return errorState;}

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

	/** Total reads written. */
	public long readsWritten=0;
	/** Total bases written. */
	public long basesWritten=0;
	/** Were any errors encountered */
	private boolean errorState=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	public static int DEFAULT_THREADS=8;

	public static final boolean verbose=false;

}