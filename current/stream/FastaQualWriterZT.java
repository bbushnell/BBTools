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
 * Synchronous Writer for separate FASTA bases and numeric QUAL scores.
 * Formats selected records on the caller's thread, then writes the batch's
 * FASTA buffer followed by its QUAL buffer. The two writes are not atomic.
 * Intended for direct construction in a single-producer context; WriterFactory
 * does not select this class. Batch IDs are ignored: writeReads serializes calls
 * in lock-acquisition order, without ID-based reordering. That serialization
 * alone is not a contract for concurrent use of every lifecycle method.
 * @author Collei, Brian Bushnell
 * @date November 21, 2025
 */
public class FastaQualWriterZT implements Writer{
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Opens both output streams immediately, before start() is called.
	 * Both names use ffFa's append and subprocess settings and request buffering.
	 * At least one mate-selection flag must be true.
	 * @param ffFa Nonnull format supplying the FASTA name and shared output options
	 * @param qf Nonnull QUAL output name
	 * @param writeR1_ Select pairnum-0 input reads
	 * @param writeR2_ Select pairnum-1 inputs, or mates of pairnum-0 inputs
	 */
	public FastaQualWriterZT(FileFormat ffFa, String qf, 
			boolean writeR1_, boolean writeR2_){
		ffoutFa=ffFa;
		fnameFa=ffFa.name();
		fnameQual=qf;
		
		writeR1=writeR1_;
		writeR2=writeR2_;
		
		assert(writeR1 || writeR2) : "Must write at least one mate";
		
		// Open output streams
		//ff.append() must be honored: app=t previously truncated (hardcoded false; replicated via stream.sh 2026-09-05)
		outstreamFa=ReadWrite.getOutputStream(fnameFa, ffFa.append(), true, ffFa.allowSubprocess());
		outstreamQual=ReadWrite.getOutputStream(fnameQual, ffFa.append(), true, ffFa.allowSubprocess());
		
		if(verbose){outstream.println("Made FastaQualWriterZT for "+fnameFa);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Marks this writer started; construction already opened the outputs. */
	@Override
	public void start(){started=true;}
	
	/** Returns selected records formatted so far, counted before batch output writes. */
	@Override
	public long readsWritten(){return readsWritten;}
	
	/** Returns bases in selected records formatted so far, without waiting or locking. */
	@Override
	public long basesWritten(){return basesWritten;}
	
	/**
	 * Wraps the existing list in one ListNum and writes it synchronously.
	 * @param list Nonnull read list; null entries are skipped
	 * @param id Batch ID retained by the wrapper but ignored by this writer
	 */
	@Override
	public final void add(ArrayList<Read> list, long id){addReads(new ListNum<Read>(list, id));}
	
	/**
	 * Writes selected reads synchronously; a null batch is ignored.
	 * @param reads Batch with a nonnull list; its ID is ignored
	 */
	@Override
	public void addReads(ListNum<Read> reads){
		if(reads==null){return;}
		writeReads(reads.list);
	}
	
	/**
	 * Converts SAM lines to Reads and writes them synchronously; null batch is ignored.
	 * Borrows sequence/quality arrays and qname without copying and preserves the
	 * SAM pair number for mate selection. Other SAM metadata is discarded; these
	 * converted Reads are not linked to mates.
	 * @param lines Batch of nonnull SAM lines; its ID is ignored
	 */
	@Override
	public void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		for(SamLine sl : lines){
			// Fixed #001: the constructor defaults to pairnum 0, which sent both SAM mates
			// to R1-only output and none to R2-only output. Preserve the selector bit.
			final Read r=new Read(sl.seq, sl.qual, sl.qname, -1, false);
			r.setPairnum(sl.pairnum());
			reads.add(r);
		}
		writeReads(reads);
	}
	
	/**
	 * Starts implicitly, formats selected nonnull records, then writes both buffers.
	 * R1 is a pairnum-0 input; R2 is a pairnum-1 input or the input's mate.
	 * Selection does not deduplicate repeated records or mates already in the list.
	 * @param reads Nonnull list; null entries are skipped
	 */
	private synchronized void writeReads(ArrayList<Read> reads){
		if(!started){start();}
		
		ByteBuilder bbFa=new ByteBuilder();
		ByteBuilder bbQual=new ByteBuilder();
		
		for(Read r : reads){
			if(r==null){continue;}
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
	
	/**
	 * Appends one FASTA/QUAL record and increments formatted-record counters.
	 * FASTA uses Read.toFasta with current Shared.FASTA_WRAP; both headers use r.id
	 * or numericID if the name is null. QUAL writes one space-separated numeric score line.
	 * Absent qualities use current Shared.FAKE_QUAL once per base. Existing quality
	 * arrays are emitted as supplied, without a length check here.
	 * @param r Selected read; a null name uses its numericID in both headers
	 * @param bbFa Destination FASTA batch buffer
	 * @param bbQual Destination QUAL batch buffer
	 */
	private void writeRead(Read r, ByteBuilder bbFa, ByteBuilder bbQual){
		// Write Fasta
		r.toFasta(bbFa);
		bbFa.nl();
		
		// Write Qual header
		// Fixed #002: match Read.toFasta's numeric-ID fallback for an unnamed read.
		//A direct null String append instead emitted the literal null in QUAL.
		bbQual.append('>');
		if(r.id==null){bbQual.append(r.numericID);}else{bbQual.append(r.id);}
		bbQual.nl();
		
		// Write Qual scores
		byte[] quals=r.quality;
		if(quals!=null){
			for(byte b : quals){
				// Read quality bytes are numeric scores, not ASCII-encoded FASTQ characters.
				// QUAL stores these values as integers (e.g. "40 40 30").
				int q=b;
				bbQual.append(q).append(' ');
			}
			if(quals.length>0){bbQual.length--;}//Trim trailing space
		}else{
			// Fake qualities if missing
			// Use the current shared fallback value.
			int fake=Shared.FAKE_QUAL;
			int len=r.length();
			for(int i=0; i<len; i++){
				bbQual.append(fake).append(' ');
			}
			if(len>0){bbQual.length--;}
		}
		bbQual.nl();
		
		readsWritten++;
		basesWritten+=r.length();
	}
	
	/**
	 * Writes a copy of a nonempty buffer to the stream, then clears the buffer.
	 * @param bb Formatted batch buffer; an empty buffer is ignored
	 * @param os Destination output stream
	 * @throws RuntimeException Wraps an IOException on the calling thread
	 */
	private void write(ByteBuilder bb, OutputStream os){
		if(bb.length()==0){return;}//was <0 (never fired); skip empty buffers
		byte[] array=bb.toBytes();
		try{
			synchronized(this){os.write(array);}
			bb.clear();
		}catch(IOException e){
			// Caller-driven: IOException propagates here instead of being lost on an internal worker.
			// There is no separate formatting producer waiting for this method to signal completion.
			throw new RuntimeException(e);
		}
	}
	
	/** Marks end of input once; callers must stop adding before finalization. */
	@Override
	public synchronized void poison(){
		if(poisoned){return;}
		poisoned=true;
		if(verbose){System.err.println("Set "+getClass().getName()+" poisoned.");}
	}
	
	/**
	 * Finalizes both output streams after poison, then releases their references.
	 * There is no internal formatting worker to join. Repeated calls after closure
	 * return the stored error flag without finalizing again.
	 * @return Stored error state including errors reported by finishWriting
	 */
	@Override
	public synchronized boolean waitForFinish(){
		if(closed){return errorState;}
		assert(poisoned);
		setError(ReadWrite.finishWriting(null, outstreamFa, fnameFa, ffoutFa.allowSubprocess()));
		setError(ReadWrite.finishWriting(null, outstreamQual, fnameQual, ffoutFa.allowSubprocess()));
		outstreamFa=null;
		outstreamQual=null;
		closed=true;
		if(verbose){System.err.println("Set "+getClass().getName()+" closed.");}
		return errorState;
	}
	
	/** Marks end of input and finalizes both streams; returns the stored error flag. */
	@Override
	public synchronized boolean poisonAndWait(){
		poison();
		return waitForFinish();
	}

	/**
	 * Marks an external error and finalizes via poisonAndWait.
	 * Formatting/writes happen on the caller's thread, so there is no internal
	 * background backlog to abandon. Finalization still invokes underlying output
	 * APIs; this synchronous design does not guarantee nonblocking I/O.
	 */
	@Override
	public synchronized void finishError(){
		setError(true);
		poisonAndWait();
	}

	/**
	 * Latches a reported error and prints a diagnostic on its first transition.
	 * @param b Error indication to merge; false cannot clear a prior error
	 * @return Updated stored error state
	 */
	private boolean setError(boolean b){
		if(b && !errorState){
			new RuntimeException("Triggered error state. this="+toString()).printStackTrace(outstream);
		}
		errorState|=b;
		return errorState;
	}
	
	/** Returns the stored error flag; write exceptions propagate separately. */
	@Override
	public synchronized boolean errorState(){return errorState;}
	
	/** Returns true when poisoned and closed with no stored error. */
	@Override
	public synchronized boolean finishedSuccessfully(){return !errorState && poisoned && closed;}
	
	/** Returns the configured output names as (fastaName,qualName). */
	@Override
	public final String fname(){return "("+fnameFa+","+fnameQual+")";}
	
	/** Returns this writer type, output names and current closure/poison flags. */
	@Override
	public final String toString(){
		return "FastaQualWriterZT"+fname()+" closed="+closed+", poisoned="+poisoned;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Configured FASTA output name. */
	private final String fnameFa;
	/** Configured QUAL output name. */
	private final String fnameQual;
	/** Primary format supplying append and subprocess options for both files. */
	private final FileFormat ffoutFa;
	
	/** FASTA output opened by construction; cleared after finalization. */
	private OutputStream outstreamFa;
	/** QUAL output opened by construction; cleared after finalization. */
	private OutputStream outstreamQual;
	
	/** Whether to select pairnum-0 inputs. */
	private final boolean writeR1;
	/** Whether to select pairnum-1 inputs or mates. */
	private final boolean writeR2;
	
	/** Selected records formatted, including those in buffers not yet written. */
	private long readsWritten=0;
	/** Bases formatted in selected records. */
	private long basesWritten=0;
	
	/** Latched error reported by finalization or finishError. */
	private boolean errorState=false;
	/** Whether start was called explicitly or by writeReads. */
	private boolean started=false;
	/** Whether end of input was marked. */
	private boolean poisoned=false;
	/** Whether normal finalization reached completion. */
	private boolean closed=false;

	/*--------------------------------------------------------------*/
	/*----------------         Diagnostics          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Enables optional lifecycle diagnostics. */
	private static final boolean verbose=false;
	
	/** Print status messages to this output stream */
	private final PrintStream outstream=System.err;
	
}
