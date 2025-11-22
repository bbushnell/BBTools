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
 * @author Collei, Brian Bushnell
 * @date November 21, 2025
 */
public class FastaQualWriter implements Writer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Constructor. */
	public FastaQualWriter(FileFormat ffFa, String qf, 
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
		
		if(verbose){outstream.println("Made FastaQualWriter for "+fnameFa);}
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
	
	private void writeReads(ArrayList<Read> reads){
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
		poisoned=true;
	}
	
	@Override
	public synchronized boolean waitForFinish(){
		if(closed) {return errorState;}
		boolean b=ReadWrite.finishWriting(null, outstreamFa, fnameFa, ffoutFa.allowSubprocess());
		boolean b2=ReadWrite.finishWriting(null, outstreamQual, fnameQual, ffoutFa.allowSubprocess());
		closed=true;
		return errorState |= (b || b2);
	}
	
	@Override
	public synchronized boolean poisonAndWait(){
		if(poisoned)
		poison();
		return waitForFinish();
	}
	
	@Override
	public boolean errorState(){return errorState;}
	
	@Override
	public boolean finishedSuccessfully() {return !errorState && poisoned;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final String fnameFa;
	public final String fnameQual;
	final FileFormat ffoutFa;
	
	OutputStream outstreamFa;
	OutputStream outstreamQual;
	
	final boolean writeR1;
	final boolean writeR2;
	
	protected long readsWritten=0;
	protected long basesWritten=0;
	
	public boolean errorState=false;
	private boolean started=false;
	private boolean poisoned=false;
	private boolean closed=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final boolean verbose=false;
	
	/** Print status messages to this output stream */
	protected PrintStream outstream=System.err;
	
}