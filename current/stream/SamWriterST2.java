package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Dual-mode SAM writer supporting single-threaded host-driven output (SamWriterST behavior)
 * and optional multi-threaded output with ordering via JobQueue.
 * @author Brian Bushnell
 * @contributor Gemini
 * @date November 18, 2025
 */
public class SamWriterST2 implements Writer {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructor using default configuration. */
	public SamWriterST2(FileFormat ffout_, ArrayList<byte[]> header_, boolean useSharedHeader_){
		this(ffout_, header_, useSharedHeader_, false, 3);
	}
	
	/**
	 * Full constructor with configuration flags.
	 * * @param ffout_ Output file format
	 * @param header_ Optional header lines
	 * @param useSharedHeader_ If true, pulls header from SamReadInputStream
	 * @param threaded_ If true, a separate thread handles writing (Producer-Consumer)
	 * @param queueCapacity_ Size of the JobQueue buffer
	 */
	public SamWriterST2(FileFormat ffout_, ArrayList<byte[]> header_, boolean useSharedHeader_, 
			boolean threaded_, int queueCapacity_){
		
		ffout=ffout_;
		fname=ffout.name();
		header=header_;
		useSharedHeader=useSharedHeader_;
		supressHeader=(ReadStreamWriter.NO_HEADER || (ffout.append() && ffout.exists()));
		supressHeaderSequences=(ReadStreamWriter.NO_HEADER_SEQUENCES || supressHeader);
		
		// Config
		ordered=ffout.ordered(); // Derive ordering requirement from file format
		threaded=threaded_;
		
		// Only create queue if we need ordering or threading
		if(ordered || threaded){
			// JobQueue handles ordering, capacity bounds, and backpressure.
			// The job type is ListNum<SamLine> since SamWriter deals in lines.
			queue=new JobQueue<ListNum<SamLine>>(queueCapacity_, ordered, true, 0);
			queue.name="*SamWriterST2";
		} else{
			// Pure lightweight ST/unordered mode (needs external synchronization if MT)
			queue=null;
		}

		if(ffout.bam()){
			outstream=ReadWrite.getBamOutputStream(fname, ffout.append());
		}else {
			outstream=ReadWrite.getOutputStream(fname, ffout.append(), true, ffout.allowSubprocess());
		}
		if(verbose) {outstream2.println("Made SamWriterST2 (Ordered: "+ordered+", Threaded: "+threaded+")");}
		
//		System.err.println(ffout.ordered()+", "+ordered+", "+threaded+", "+(queue!=null));
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public void start(){
		synchronized(this) {
			if(started) {return;}
			writeHeader();
			started=true;
			if(threaded && queue != null){
				writerThread=new Thread(new WriterRunnable());
				writerThread.start();
			}
		}
	}

	@Override
	public synchronized void poison(){
		// Ensure only one thread executes the shutdown sequence
		if(poisoned){return;}
		poisoned=true;
		
		if(queue != null){
			// 1. Calculate the ID for the poison pill. It must be > maxSeen to guarantee order.
			long poisonID=queue.maxSeen() + 1;
			
			// 2. Create the poison pill (ListNum acts as the HasID object)
			// SamWriter uses ListNum<SamLine>
			ListNum<SamLine> poison=new ListNum<SamLine>(null, poisonID, ListNum.POISON);
			
			// 3. Inject the poison pill into the queue using the explicit JobQueue API.
			queue.poison(poison, false); 
			
			// 4. CRITICAL: If running in Host-Driven (unthreaded) mode, the producer thread 
			// calling poison() MUST perform the final drain using the blocking take().
			if(!threaded){
				// Drain loop terminates when queue.take() returns null (the poison pill)
				while(true){
					ListNum<SamLine> job=queue.take(); 
					
					if(job==null) break;
					
					writeLines(job);
				}
			}
			// If threaded, the worker thread handles the take() loop and exits gracefully.
		}
	}

	@Override
	public synchronized boolean poisonAndWait(){
		poison();
		return waitForFinish();
	}
	
	@Override
	public final synchronized boolean waitForFinish(){
		if(closed) {return errorState;}
		// If threaded, wait for the worker to process the poison pill and exit.
		if(threaded && writerThread != null){
			try{
				writerThread.join();
			} catch (InterruptedException e){
				Thread.currentThread().interrupt();
			}
		}
		
		boolean b=ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess());
		closed=true;
		return errorState|=b;
	}
	
	public final void add(ArrayList<Read> list, long id) {addReads(new ListNum<Read>(list, id));}

	/** Add reads for writing (will be converted to SamLines). */
	@Override
	public final void addReads(ListNum<Read> reads){
		if(reads==null){return;}
		ArrayList<SamLine> lines=toSamLines(reads.list);
		addLines(new ListNum<SamLine>(lines, reads.id));
	}

	/** Add already-formatted SamLines for writing. */
	@Override
	public final void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		if(!started){start();}
		
		if(queue==null){
			// LIGHTWEIGHT MODE (Unordered/Unthreaded): Direct synchronized write.
			writeLines(lines);
		} else{
			// 1. Add to queue (blocks if full/backpressure enabled)
			queue.add(lines);
			
			if(!threaded){
				// HOST-DRIVEN MODE (Ordered/Unthreaded): 
				// Producer acts as the draining consumer immediately after adding.
				// Uses the non-blocking poll() to drain any contiguous jobs that are ready.
				assert(ordered);
				synchronized(queue) {//Essential to maintain ordering
					for(ListNum<SamLine> job=queue.poll(); job!=null; job=queue.poll()){
						if(verbose) {System.err.println("Writing job "+job.id+", expected="+expected);}
						writeLines(job);
					}
				}
			}
			// If threaded, the background thread handles the draining.
		}
	}
	
	/** Write lines to the stream (internal logic). */
	private void writeLines(ListNum<SamLine> lines){
		assert(!ordered || lines.id()==expected++) : lines.id+", "+expected;
		ByteBuilder bb=new ByteBuilder();
		for(SamLine sl : lines){
			if(sl==null) {continue;}
			sl.toBytes(bb);
			bb.nl();
			readsWritten++;
			basesWritten+=sl.length();

			if(bb.length()>=BUFFER_SIZE){
				write(bb);
			}
		}
		if(bb.length()>0){
			write(bb);
		}
	}
	
	/** The actual synchronized I/O call. */
	private void write(ByteBuilder bb) {
		if(bb.length()<0){return;}
		byte[] array=bb.toBytes();
		try{
			// CRITICAL: Synchronize on the shared stream object to prevent file corruption
			synchronized(outstream) {outstream.write(array);}
			bb.clear();
		}catch(IOException e){
			throw new RuntimeException(e);
		}
	}
	
	@Override
	public long readsWritten(){return readsWritten;}

	@Override
	public long basesWritten(){return basesWritten;}

	@Override
	public boolean errorState() {return errorState;}
	
	@Override
	public boolean finishedSuccessfully() {return !errorState && closed;}

	/*--------------------------------------------------------------*/
	/*----------------        Writer Thread         ----------------*/
	/*--------------------------------------------------------------*/
	
	private class WriterRunnable implements Runnable{
		@Override
		public void run(){
			// Loop uses the blocking take() and terminates when it retrieves the poison pill (null)
			for(ListNum<SamLine> job=queue.take(); job!=null; job=queue.take()){
				writeLines(job);
			}
		}
	}


	/*--------------------------------------------------------------*/
	/*----------------         Helper Methods       ----------------*/
	/*--------------------------------------------------------------*/

	// (toSamLines and addSamLine methods remain unchanged)
	public static ArrayList<SamLine> toSamLines(ArrayList<Read> reads) {
		ArrayList<SamLine> samLines=new ArrayList<SamLine>();
		for(final Read r1 : reads){
			if(r1==null) {continue;}
			Read r2=(r1==null ? null : r1.mate);
			SamLine sl1=(r1==null ? null : (ReadStreamWriter.USE_ATTACHED_SAMLINE && r1.samline!=null ? r1.samline : new SamLine(r1, 0)));
			SamLine sl2=(r2==null ? null : (ReadStreamWriter.USE_ATTACHED_SAMLINE && r2.samline!=null ? r2.samline : new SamLine(r2, 1)));
			if(!SamLine.KEEP_NAMES && sl1!=null && sl2!=null && ((sl2.qname==null) || !sl2.qname.equals(sl1.qname))){
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
	
	// (getHeader and writeHeader methods remain unchanged)
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
			outstream2.println("Warning: Header was null, creating empty header");
			headerLines=new ArrayList<byte[]>();
		}
		return headerLines;
	}

	protected void writeHeader(){
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
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	final String fname;
	final FileFormat ffout;
	final boolean useSharedHeader;
	final boolean supressHeader;
	final boolean supressHeaderSequences;
	boolean headerWritten=false;
	final ArrayList<byte[]> header;
	final OutputStream outstream;
	
	final boolean ordered;
	final boolean threaded;
	final JobQueue<ListNum<SamLine>> queue;
	Thread writerThread;
	
	private long readsWritten=0;
	private long basesWritten=0;
	private boolean errorState=false;
	private boolean started=false;
	private boolean closed=false;
	private boolean poisoned=false;

	//Expected next list number, for ordered mode assertions
	private long expected=0;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private static final int BUFFER_SIZE=65536;
	public static final boolean verbose=false;
	protected PrintStream outstream2=System.err;

}