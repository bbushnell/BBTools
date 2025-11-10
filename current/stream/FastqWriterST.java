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
 * Single-threaded FASTQ writer with simple buffering.
 * Simpler alternative to FastqWriter for cases where threading overhead isn't worth it.
 * 
 * @author Isla
 * @date November 10, 2025
 */
public class FastqWriterST implements Writer {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Constructor. */
	public FastqWriterST(String out_, boolean writeR1_, boolean writeR2_, boolean overwrite){
		this(FileFormat.testOutput(out_, FileFormat.FASTQ, null, true, overwrite, false, true), 
			writeR1_, writeR2_);
	}
	
	/** Constructor. */
	public FastqWriterST(FileFormat ffout_, boolean writeR1_, boolean writeR2_){
		ffout=ffout_;
		fname=ffout_.name();
		writeR1=writeR1_;
		writeR2=writeR2_;
		format=(ffout.format()==UNKNOWN ? FASTQ : ffout.format());
		assert(format==FASTQ || format==FASTA || format==HEADER) : ffout;
		
		assert(writeR1 || writeR2) : "Must write at least one mate";
		
		// Open output stream
		outstream=ReadWrite.getOutputStream(fname, false, true, false);
		if(verbose){outstream2.println("Made FastqWriterST");}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void start(){
		started=true;
	}
	
	@Override
	public void close(){
		if(!finished){poison();}
		ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess());
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
		
		ByteBuilder bb=new ByteBuilder();
		try{
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
			
			if(bb.length()>0){
				outstream.write(bb.toBytes());
			}
		}catch(IOException e){
			errorState=true;
			throw new RuntimeException("Error writing FASTQ", e);
		}
	}
	
	@Override
	public void poison(){
		finished=true;
	}
	
	@Override
	public boolean waitForFinish(){
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
	public boolean finishedSuccessfully() {return !errorState && finished;}
	
	/*--------------------------------------------------------------*/
	/*----------------         Helper Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	private void writeFastq(ArrayList<Read> reads, ByteBuilder bb) {
		for(Read r : reads){
			final Read r1=(r.pairnum()==0 ? r : null);
			final Read r2=(r.pairnum()==1 ? r : r.mate);
			if(writeR1 && r1!=null){
				r1.toFastq(bb);
				bb.nl();
				readsWritten++;
				basesWritten+=r1.length();
			}
			if(writeR2 && r2!=null){
				r2.toFastq(bb);
				bb.nl();
				readsWritten++;
				basesWritten+=r2.length();
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
				readsWritten++;
				basesWritten+=r1.length();
			}
			if(writeR2 && r2!=null){
				r2.toFasta(bb);
				bb.nl();
				readsWritten++;
				basesWritten+=r2.length();
			}
		}
	}
	
	private void writeHeader(ArrayList<Read> reads, ByteBuilder bb) {
		for(Read r : reads){
			final Read r1=(r.pairnum()==0 ? r : null);
			final Read r2=(r.pairnum()==1 ? r : r.mate);
			if(writeR1 && r1!=null){
				bb.appendln(r1.id);
				readsWritten++;
			}
			if(writeR2 && r2!=null){
				bb.appendln(r2.id);
				readsWritten++;
			}
		}
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
	/** Write R1 reads (pairnum==0) */
	final boolean writeR1;
	/** Write R2 reads (pairnum==1 or mate) */
	final boolean writeR2;
	/** Number of reads written */
	protected long readsWritten=0;
	/** Number of bases written */
	protected long basesWritten=0;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** True after start() called */
	private boolean started=false;
	/** True after poison() called */
	private boolean finished=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private static final int FASTQ=FileFormat.FASTQ;
	private static final int FASTA=FileFormat.FASTA;
	private static final int HEADER=FileFormat.HEADER;
	private static final int UNKNOWN=FileFormat.UNKNOWN;
	
	public static final boolean verbose=false;
	
	/** Print status messages to this output stream */
	protected PrintStream outstream2=System.err;
	
}