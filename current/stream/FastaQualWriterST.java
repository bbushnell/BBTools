package stream;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Writer for Fasta + Quality files using an internal thread and JobQueue.
 * Ensures ordered output even with multi-threaded producers.
 * Writes two files in lockstep: .fa (bases) and .qual (quality scores).
 * * @author Collei, Brian Bushnell
 * @date November 22, 2025
 */
public class FastaQualWriterST implements Writer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Constructor. */
	public FastaQualWriterST(FileFormat ffFa, String qf, 
			boolean writeR1_, boolean writeR2_){
		ffoutFa=ffFa;
		ffoutQual=FileFormat.testOutput(qf, FileFormat.QUAL, null, true, true, false, false);
		fnameFa=ffFa.name();
		fnameQual=qf;
		
		writeR1=writeR1_;
		writeR2=writeR2_;
		
		assert(writeR1 || writeR2) : "Must write at least one mate";
		
		// Config
		// Assuming we always want ordered output for legacy formats to match input
		ordered=true; 
		queueCapacity=4;
		
		// JobQueue handles ordering
		queue=new JobQueue<ListNum<Read>>(queueCapacity, ordered, true, 0);
		queue.name="*FastaQualWriterST";
		
		// Open output streams
		outstreamFa=ReadWrite.getOutputStream(fnameFa, false, true, ffFa.allowSubprocess());
		outstreamQual=ReadWrite.getOutputStream(fnameQual, false, true, ffoutQual.allowSubprocess()); // Assuming generic handling for .qual
		
		if(verbose){outstream.println("Made FastaQualWriterST for "+fnameFa);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void start(){
		synchronized(this) {
			if(started) {return;}
			started=true;
			writerThread=new Thread(new WriterRunnable());
			writerThread.start();
			if(verbose){outstream.println("Started "+getClass().getName());}
		}
	}
	
	@Override
	public final void add(ArrayList<Read> list, long id) {addReads(new ListNum<Read>(list, id));}
	
	@Override
	public void addReads(ListNum<Read> reads){
		if(reads==null){return;}
		if(!started){start();}
		
		// Add to queue (blocks if full/backpressure enabled)
		queue.add(reads);
	}
	
	@Override
	public void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		for(SamLine sl : lines) {
			reads.add(new Read(sl.seq, sl.qual, sl.qname, -1, false));
		}
		addReads(new ListNum<Read>(reads, lines.id));
	}
	
	@Override
	public synchronized void poison(){
		if(poisoned){return;}
		poisoned=true;
		
		// Calculate poison ID
		long poisonID=queue.maxSeen() + 1;
		ListNum<Read> poison=new ListNum<Read>(null, poisonID, ListNum.POISON);
		queue.poison(poison, false);
		
		if(verbose){outstream.println("Poisoned "+getClass().getName());}
	}
	
	@Override
	public synchronized boolean waitForFinish(){
		if(closed) {return errorState;}
		
		// Ensure poison is sent if not already
		if(!poisoned){poison();}
		
		// Wait for worker thread to finish draining
		if(writerThread != null){
			try{
				writerThread.join();
			} catch (InterruptedException e){
				// Ignore
			}
		}
		
		boolean b=ReadWrite.finishWriting(null, outstreamFa, fnameFa, ffoutFa.allowSubprocess());
		boolean b2=ReadWrite.finishWriting(null, outstreamQual, fnameQual, ffoutQual.allowSubprocess());
		
		outstreamFa=null;
		outstreamQual=null;
		closed=true;
		
		return errorState |= (b || b2);
	}
	
	@Override
	public synchronized boolean poisonAndWait(){
		poison();
		return waitForFinish();
	}
	
	private boolean setError(boolean b) {
		if(b && !errorState) {
			new RuntimeException("Triggered error state. this="+toString()).printStackTrace(outstream);
		}
		errorState|=b;
		return errorState;
	}
	
	@Override
	public boolean errorState(){return errorState;}
	
	@Override
	public boolean finishedSuccessfully() {return !errorState && closed;}
	
	@Override
	public final String fname() {return "("+fnameFa+","+fnameQual+")";}
	
	@Override
	public final String toString() {
		return "FastaQualWriterST"+fname()+" closed="+closed+", poisoned="+poisoned;
	}
	
	@Override
	public long readsWritten(){return readsWritten;}
	
	@Override
	public long basesWritten(){return basesWritten;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Writer Thread         ----------------*/
	/*--------------------------------------------------------------*/
	
	private class WriterRunnable implements Runnable{
		@Override
		public void run(){
			// JobQueue handles ordering. We just take() the next valid job.
			for(ListNum<Read> job=queue.take(); job!=null && !job.poison(); job=queue.take()){
				writeReads(job.list);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Logic          ----------------*/
	/*--------------------------------------------------------------*/
	
	private void writeReads(ArrayList<Read> reads){
		ByteBuilder bbFa=new ByteBuilder();
		ByteBuilder bbQual=new ByteBuilder();
		
		for(Read r : reads){
			if(r==null) {continue;}
			final Read r1=(r.pairnum()==0 ? r : null);
			final Read r2=(r.pairnum()==1 ? r : r.mate);
			
			if(writeR1 && r1!=null){
				writeRead(r1, bbFa, bbQual);
			}
			if(writeR2 && r2!=null){
				writeRead(r2, bbFa, bbQual);
			}
		}

		// Write both buffers to their respective streams
		// Single thread, so no race condition between Fa and Qual writes here
		write(bbFa, outstreamFa);
		write(bbQual, outstreamQual);
	}
	
	private void writeRead(Read r, ByteBuilder bbFa, ByteBuilder bbQual) {
		// Write Fasta
		r.toFasta(bbFa);
		bbFa.nl();
		
		// Write Qual header
		bbQual.append('>');
		bbQual.append(r.id);
		bbQual.nl();
		
		// Write Qual scores
		byte[] quals=r.quality;
		if(quals!=null){
			for(byte b : quals){
				int q = b; 
				bbQual.append(q).append(' ');
			}
			if(quals.length>0) {bbQual.length--;} //Trim trailing space
		}else{
			int fake=Shared.FAKE_QUAL;
			int len=r.length();
			for(int i=0; i<len; i++){
				bbQual.append(fake).append(' ');
			}
			if(len>0) {bbQual.length--;}
		}
		bbQual.nl();
		
		readsWritten++;
		basesWritten+=r.length();
	}
	
	private void write(ByteBuilder bb, OutputStream os) {
		if(bb.length()<0){return;}
		byte[] array=bb.toBytes();
		try{
			os.write(array); //No synchronized(this) needed; only writerThread calls this
			bb.clear();
		}catch(IOException e){
			throw new RuntimeException(e);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private final String fnameFa;
	private final String fnameQual;
	private final FileFormat ffoutFa;
	private final FileFormat ffoutQual;
	
	private OutputStream outstreamFa;
	private OutputStream outstreamQual;
	
	private final boolean writeR1;
	private final boolean writeR2;
	
	private final boolean ordered;
	private final int queueCapacity;
	
	private final JobQueue<ListNum<Read>> queue;
	private Thread writerThread;
	
	private long readsWritten=0;
	private long basesWritten=0;
	
	public boolean errorState=false;
	private boolean started=false;
	private boolean poisoned=false;
	private boolean closed=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final boolean verbose=false;
	
	/** Print status messages to this output stream */
	private final PrintStream outstream=System.err;
	
}