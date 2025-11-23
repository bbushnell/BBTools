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
 * Single-threaded Writer for Fasta + Quality files.
 * Writes two files in lockstep: .fa (bases) and .qual (quality scores).
 * NOT ORDERED.
 * @author Collei, Brian Bushnell
 * @date November 21, 2025
 */
public class FastaQualWriterZT implements Writer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Constructor. */
	public FastaQualWriterZT(FileFormat ffFa, String qf, 
			boolean writeR1_, boolean writeR2_){
		ffoutFa=ffFa;
		fnameFa=ffFa.name();
		fnameQual=qf;
		
		writeR1=writeR1_;
		writeR2=writeR2_;
		
		assert(writeR1 || writeR2) : "Must write at least one mate";
		
		// Open output streams
		outstreamFa=ReadWrite.getOutputStream(fnameFa, false, true, ffFa.allowSubprocess());
		outstreamQual=ReadWrite.getOutputStream(fnameQual, false, true, ffFa.allowSubprocess());
		
		if(verbose){outstream.println("Made FastaQualWriterZT for "+fnameFa);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void start(){
		started=true;
	}
	
	@Override
	public long readsWritten(){
		return readsWritten;
	}
	
	@Override
	public long basesWritten(){
		return basesWritten;
	}
	
	@Override
	public final void add(ArrayList<Read> list, long id) {addReads(new ListNum<Read>(list, id));}
	
	@Override
	public void addReads(ListNum<Read> reads){
		if(reads==null){return;}
		writeReads(reads.list);
	}
	
	@Override
	public void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		for(SamLine sl : lines) {
			reads.add(new Read(sl.seq, sl.qual, sl.qname, -1, false));
		}
		writeReads(reads);
	}
	
	private synchronized void writeReads(ArrayList<Read> reads){
		if(!started){start();}
		
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
				// Assuming Read stores 0-based Phred scores. 
				// .qual files normally expect integers (e.g. "40 40 30").
				int q = b;
				bbQual.append(q).append(' ');
			}
			if(quals.length>0) {bbQual.length--;} //Trim trailing space
		}else{
			// Fake qualities if missing
			// Using a default value (e.g. 30)
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
			synchronized(this) {os.write(array);}
			bb.clear();
		}catch(IOException e){
			throw new RuntimeException(e);
		}
	}
	
	@Override
	public synchronized void poison(){
		if(poisoned) {return;}
		poisoned=true;
		if(verbose) {System.err.println("Set "+getClass().getName()+" poisoned.");}
	}
	
	@Override
	public synchronized boolean waitForFinish(){
		if(closed) {return errorState;}
		assert(poisoned);
		setError(ReadWrite.finishWriting(null, outstreamFa, fnameFa, ffoutFa.allowSubprocess()));
		setError(ReadWrite.finishWriting(null, outstreamQual, fnameQual, ffoutFa.allowSubprocess()));
		outstreamFa=null;
		outstreamQual=null;
		closed=true;
		if(verbose) {System.err.println("Set "+getClass().getName()+" closed.");}
		return errorState;
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
	public synchronized boolean errorState(){return errorState;}
	
	@Override
	public synchronized boolean finishedSuccessfully() {return !errorState && poisoned && closed;}
	
	@Override
	public final String fname() {return "("+fnameFa+","+fnameQual+")";}
	
	@Override
	public final String toString() {
		return "FastaQualWriterZT"+fname()+" closed="+closed+", poisoned="+poisoned;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private final String fnameFa;
	private final String fnameQual;
	private final FileFormat ffoutFa;
	
	private OutputStream outstreamFa;
	private OutputStream outstreamQual;
	
	private final boolean writeR1;
	private final boolean writeR2;
	
	private long readsWritten=0;
	private long basesWritten=0;
	
	private boolean errorState=false;
	private boolean started=false;
	private boolean poisoned=false;
	private boolean closed=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private static final boolean verbose=false;
	
	/** Print status messages to this output stream */
	private final PrintStream outstream=System.err;
	
}