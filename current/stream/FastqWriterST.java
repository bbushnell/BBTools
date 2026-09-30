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
 * Writes FASTQ, FASTA or read names on the submitting thread with per-batch buffering.
 * WriterFactory selects other implementations; this class has no formatting queue or
 * batch-ID reordering. Serialize submissions and lifecycle operations at the caller.
 * addReads/addLines and counters are not synchronized. The lock around individual byte
 * writes does not make concurrent submissions safe. Do not concurrently mutate supplied
 * records or shared serialization settings. Output opens during construction.
 *
 * @author Isla
 * @date November 10, 2025
 */
public class FastqWriterST implements Writer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Tests an output path and opens it with append disabled.
	 * @param out_ Output path; unspecified format defaults to FASTQ
	 * @param writeR1_ Include pairnum0 list entries
	 * @param writeR2_ Include pairnum1 entries or linked mates of pairnum0 entries
	 * @param overwrite Allow replacement of an existing output file
	 */
	public FastqWriterST(String out_, boolean writeR1_, boolean writeR2_, boolean overwrite){
		this(FileFormat.testOutput(out_, FileFormat.FASTQ, null, true, overwrite, false, true),
			writeR1_, writeR2_);
	}

	/** Opens buffered output and retains the descriptor's append setting.
	 * Output opening explicitly disables subprocess allowance. Unknown format becomes
	 * FASTQ; supported formats are FASTQ, FASTA and HEADER. At least one mate is required.
	 * @param ffout_ Nonnull output descriptor, opened immediately
	 * @param writeR1_ Include pairnum0 list entries
	 * @param writeR2_ Include pairnum1 entries or linked mates of pairnum0 entries
	 */
	public FastqWriterST(FileFormat ffout_, boolean writeR1_, boolean writeR2_){
		ffout=ffout_;
		fname=ffout_.name();
		writeR1=writeR1_;
		writeR2=writeR2_;
		format=(ffout.format()==UNKNOWN ? FASTQ : ffout.format());
		assert(format==FASTQ || format==FASTA || format==HEADER) : ffout;

		assert(writeR1 || writeR2) : "Must write at least one mate";

		// Open output stream
		//ffout.append() must be honored: app=t previously truncated (hardcoded false; replicated via stream.sh 2026-09-05)
		outstream=ReadWrite.getOutputStream(fname, ffout.append(), true, false);
		if(verbose){outstream2.println("Made FastqWriterST");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Marks submissions started; starts no formatting thread and emits no header. */
	@Override
	public void start(){started=true;}

	/** Number of selected records formatted, counted before their output succeeds. */
	@Override
	public long readsWritten(){return readsWritten;}

	/** Selected FASTQ/FASTA bases formatted; HEADER mode does not update this counter. */
	@Override
	public long basesWritten(){return basesWritten;}

	/** Wraps a nonnull Read list and submits it immediately.
	 * @param list Read entries; null entries are skipped
	 * @param id Accepted for the Writer API but ignored for ordering
	 */
	@Override
	public final void add(ArrayList<Read> list, long id){addReads(new ListNum<Read>(list, id));}

	/** Formats selected reads in list order; a null wrapper is ignored.
	 * Caller owns serialization of submissions and lifecycle calls.
	 * @param reads Optional wrapper containing a nonnull list; its ID is ignored
	 */
	@Override
	public void addReads(ListNum<Read> reads){//TODO: NOT suitable for MT producers
		if(reads==null){return;}
		writeReads(reads.list);
	}

	/** Converts SAM sequence, quality and name fields to unpaired Read wrappers.
	 * Arrays are shared, not copied, and SAM mate flags are currently lost. Consequently
	 * writeR1 selects all converted entries and writeR2 alone selects none.
	 * @param lines Optional wrapper containing a nonnull list of nonnull SAM records;
	 * its ID is ignored
	 */
	@Override
	public void addLines(ListNum<SamLine> lines){
		if(lines==null){return;}
		ArrayList<Read> reads=new ArrayList<Read>(lines.size());
		//TODO: Probable bug [FastqWriterST#001] - this Read constructor leaves pairnum0
		//and no mate, so writeR2-only drops every SAM entry. Retain SAM mate identity
		//before selection; no direct Java callers found during the 2026-09-30 review.
		for(SamLine sl : lines){
			reads.add(new Read(sl.seq, sl.qual, sl.qname, -1, false));
		}
		writeReads(reads);
	}

	/** Starts lazily, formats a whole batch in memory and submits its bytes once.
	 * Uses current delegated Read/FASTQ serialization settings; does not flush the backend.
	 * @param reads Nonnull list; null entries are skipped
	 */
	private void writeReads(ArrayList<Read> reads){
		if(!started){start();}

		ByteBuilder bb=new ByteBuilder();
		// Format reads
		if(format==FASTQ){
			writeFastq(reads, bb);
		}else if(format==FASTA){
			writeFasta(reads, bb);
		}else if(format==HEADER){
			writeHeader(reads, bb);
		}else{
			throw new RuntimeException("Bad format: "+format);
		}

		write(bb);
		bb=null;
	}

	/** Writes a copy of buffered bytes under this object's output lock and clears on success.
	 * IOException propagates as RuntimeException without updating the cached error flag.
	 * @param bb Buffer to copy and clear, including an empty buffer
	 */
	private void write(ByteBuilder bb){
		if(bb.length()<0){return;}
		byte[] array=bb.toBytes();
		try{
			synchronized(this){outstream.write(array);}
			bb.clear();
		}catch(IOException e){
			throw new RuntimeException(e);
		}
	}

	/** Marks end-of-input for the success flag; does not flush, close or reject later calls. */
	@Override
	public synchronized void poison(){
		poisoned=true;
	}

	/** Finalizes output once and caches the returned status without setting poisoned.
	 * Delegates stream/backend finalization to ReadWrite using the descriptor's
	 * subprocess allowance; output opening separately disabled that allowance.
	 * @return Cached error state, including the output finalizer's result
	 */
	@Override
	public synchronized boolean waitForFinish(){
		if(closed){return errorState;}
		boolean b=ReadWrite.finishWriting(null, outstream, fname, ffout.allowSubprocess());
		closed=true;
		return errorState|=b;
	}

	/** Finalizes output; the current conditional does not mark a fresh writer poisoned.
	 * Call poison separately when the success flag is required until the noted defect is fixed.
	 * @return Cached error state after finalization
	 */
	@Override
	public synchronized boolean poisonAndWait(){
		//TODO: Probable bug [FastqWriterST#002] - "if(poisoned) poison();" only re-poisons
		//an already poisoned writer. Found 2026-09-03 while adding finishError and left
		//outside that change's scope; originally called harmless because poison is idempotent.
		//2026-09-30 source review: a fresh call leaves poisoned=false, so
		//finishedSuccessfully() remains false even when finalization reports no error.
		if(poisoned){poison();}
		return waitForFinish();
	}

	/** Marks an external error, then uses the same finalization path.
	 * There is no own formatting backlog to abandon, but backend finalization may block.
	 */
	@Override
	public synchronized void finishError(){
		errorState=true;
		poisonAndWait();
	}

	/** Cached error flag; direct formatting/write exceptions do not set it automatically. */
	@Override
	public boolean errorState(){return errorState;}

	/** Tests only cached error and poison flags; does not check whether output is closed. */
	@Override
	public boolean finishedSuccessfully(){return !errorState && poisoned;}

	/** Returns the retained output path. */
	@Override
	public final String fname(){return fname;}

	/*--------------------------------------------------------------*/
	/*----------------         Helper Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Appends selected FASTQ records and counts their reads/bases before output.
	 * Selects pairnum0 entries for R1 and pairnum1 entries or linked mates for R2.
	 * @param reads Nonnull list; null entries are skipped
	 * @param bb Destination batch buffer
	 */
	private void writeFastq(ArrayList<Read> reads, ByteBuilder bb){
		for(Read r : reads){
			if(r==null){continue;}
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

	/** Appends selected FASTA records and counts their reads/bases before output.
	 * Selects pairnum0 entries for R1 and pairnum1 entries or linked mates for R2.
	 * @param reads Nonnull list; null entries are skipped
	 * @param bb Destination batch buffer
	 */
	private void writeFasta(ArrayList<Read> reads, ByteBuilder bb){
		for(Read r : reads){
			if(r==null){continue;}
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

	/** Appends selected read names and counts records without updating bases.
	 * Selects pairnum0 entries for R1 and pairnum1 entries or linked mates for R2.
	 * @param reads Nonnull list; null entries are skipped
	 * @param bb Destination batch buffer
	 */
	private void writeHeader(ArrayList<Read> reads, ByteBuilder bb){
		for(Read r : reads){
			if(r==null){continue;}
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

	/** Retained output path. */
	public final String fname;
	/** Retained output descriptor, including append and finalization settings. */
	final FileFormat ffout;
	/** Effective format; UNKNOWN is converted to FASTQ. */
	public final int format;
	/** Backend opened during construction. */
	OutputStream outstream;
	/** Include pairnum0 list entries. */
	final boolean writeR1;
	/** Include pairnum1 entries or linked mates of pairnum0 entries. */
	final boolean writeR2;
	/** Selected records formatted, before successful output. */
	protected long readsWritten=0;
	/** Selected FASTQ/FASTA bases formatted; HEADER does not increment this. */
	protected long basesWritten=0;
	/** Cached finalization/external error flag; not all thrown exceptions set it. */
	public boolean errorState=false;
	/** True after explicit or lazy start. */
	private boolean started=false;
	/** True after poison; the current poisonAndWait conditional may leave it false. */
	private boolean poisoned=false;
	/** True after the output finalizer returns, including a reported error. */
	private boolean closed=false;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** FASTQ output format. */
	private static final int FASTQ=FileFormat.FASTQ;
	/** FASTA output format. */
	private static final int FASTA=FileFormat.FASTA;
	/** Read-name-only output format. */
	private static final int HEADER=FileFormat.HEADER;
	/** Unspecified descriptor format. */
	private static final int UNKNOWN=FileFormat.UNKNOWN;

	/** Enables constructor diagnostics when compiled true. */
	public static final boolean verbose=false;

	/** Destination for optional diagnostic messages. */
	protected PrintStream outstream2=System.err;

}
