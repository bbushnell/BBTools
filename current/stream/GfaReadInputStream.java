package stream;

import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
import shared.Shared;

/**
 * Blockwise reader for sequence fields from GFA segment lines.
 * Construction opens a ByteFile and captures the shared batch length and amino flag.
 * Parsing runs on the caller; the selected backend may use worker threads. Configure
 * shared input settings before use and coordinate iteration, close and restart externally.
 * The synchronized methods alone do not make concurrent use safe.
 *
 * Nonempty lines beginning with byte S supply tab-delimited name and sequence fields;
 * other lines are skipped. Names use US-ASCII decoding and sequence bytes are copied.
 * This reader does not validate the full GFA grammar or trim descriptions. Read
 * construction applies validation and optional header sanitation according to current
 * Read settings and the captured input mode; sanitation is not description trimming.
 *
 * Batches are bounded by read count, not sequence bytes. hasMore may prefetch without
 * consuming; nextList transfers the whole buffered list and returns null at exhaustion.
 * Short fills close the backend before returning buffered reads. An exactly full batch
 * defers the next EOF probe until another refill. Callers should still close explicitly
 * on early termination or failure. Cached error status and thrown exceptions are separate.
 *
 * @author Brian Bushnell
 * @date November 21, 2025
 * @contributor Shinobu (correctness repairs and documentation)
 */
public class GfaReadInputStream extends ReadInputStream{

	/**
	 * Prints the first read, if present, then closes the input in finally.
	 * Output may precede a later close-error rejection; a thrown close can replace
	 * an active processing exception.
	 * @param args Input GFA filename as the first argument
	 * @throws RuntimeException If normal processing completes and input closure reports an error
	 */
	public static void main(final String[] args){
		final GfaReadInputStream fris=new GfaReadInputStream(args[0], true);
		boolean error=false;
		try{
			//[stream/GfaReadInputStream#003] Empty input has no first read; always close after an opened reader is used.
			final ArrayList<Read> list=fris.nextList();
			if(list!=null){System.out.println(list.get(0).toText(false));}
		}finally{
			error=fris.close();
		}
		if(error){throw new RuntimeException("Error closing GFA input.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Resolves the input format and immediately opens the backend.
	 * @param fname Input name; GFA is the default format
	 * @param allowSubprocess_ Whether the format descriptor permits subprocess decompression
	 */
	public GfaReadInputStream(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.GFA, null, allowSubprocess_, false));
	}

	/**
	 * Opens the selected ByteFile and captures the amino flag and stdio indicator.
	 * Warns when the supplied descriptor does not report GFA; it does not convert
	 * another format. This reader always reports unpaired input.
	 * @param ff Input format, name and backend options
	 */
	public GfaReadInputStream(FileFormat ff){
		if(verbose){System.err.println("GfaReadInputStream("+ff+")");}
		flag=(Shared.AMINO_IN ? Read.AAMASK : 0);
		stdin=ff.stdio();
		if(!ff.gfa()){
			System.err.println("Warning: Did not find expected gfa file extension for filename "+ff.name());
		}
		tf=ByteFile.makeByteFile(ff);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Tests for buffered reads, refilling only when the buffer is exhausted and
	 * the backend is open. Prefetching updates generated, not consumed.
	 * An already buffered list remains available after close.
	 * @return True if the current buffer contains an unconsumed entry
	 */
	@Override
	public boolean hasMore(){
		if(buffer==null || next>=buffer.size()){
			//[stream/GfaReadInputStream#004] Closed zero-record input is ordinary EOF, including repeated queries.
			if(tf.isOpen()){
				fillBuffer();
			}
		}
		return (buffer!=null && next<buffer.size());
	}

	/**
	 * Transfers the entire buffered list and clears this reader's reference to it.
	 * Refills when needed, even if the backend has already been closed; an empty
	 * parsed list becomes null. Returned entries are added to consumed.
	 * @return Next nonempty batch, or null when no entries are available
	 * @throws RuntimeException If the vestigial single-read index is nonzero
	 */
	@Override
	public synchronized ArrayList<Read> nextList(){
		if(next!=0){throw new RuntimeException("'next' should not be used when doing blockwise access.");}
		if(buffer==null || next>=buffer.size()){fillBuffer();}
		ArrayList<Read> list=buffer;
		buffer=null;
		if(list!=null && list.size()==0){list=null;}
		consumed+=(list==null ? 0 : list.size());
		return list;
	}

	/**
	 * Replaces an exhausted buffer with up to the captured batch length, advancing IDs and counts.
	 * A short batch closes the backend and folds its returned error flag. A full
	 * batch does not fetch another line merely to test EOF. Parsing or close
	 * exceptions propagate and can interrupt counter updates.
	 */
	private synchronized void fillBuffer(){
		assert(buffer==null || next>=buffer.size());
		buffer=null;
		next=0;
		buffer=toReadList(tf, BUF_LEN, nextReadID, flag);
		int bsize=(buffer==null ? 0 : buffer.size());
		nextReadID+=bsize;
		//[stream/GfaReadInputStream#006] ByteFile.close returns true on error; retain that status even before public close.
		if(bsize<BUF_LEN && tf.close()){errorState=true;}
		generated+=bsize;
		//Defensive: toReadList returns a possibly empty list, never null; nextList maps an empty list to null.
		if(buffer==null){
			if(!errorState){
				errorState=true;
				System.err.println("Null buffer in GfaReadInputStream.");
			}
		}
	}

	/**
	 * Collects reads from nonempty lines whose first byte is S, skipping other lines.
	 * Zero-based tab fields 1 and 2 provide the US-ASCII name and copied sequence; qualities
	 * are passed as null. The initial capacity of at most 400 is not a batch cap.
	 * No complete GFA record/graph validation or description trimming is performed.
	 * @param bf Backend supplying lines on this caller
	 * @param maxReadsToReturn Maximum entries in the returned list
	 * @param numericID ID assigned to the first converted entry
	 * @param flag Read flags captured by the constructor
	 * @return New list, possibly empty, never null on normal return
	 */
	private ArrayList<Read> toReadList(final ByteFile bf, final int maxReadsToReturn,
			long numericID, final int flag){
		ArrayList<Read> list=new ArrayList<Read>(Data.min(400, maxReadsToReturn));
		//Only lines starting with S become Reads; other record types are skipped. IDs advance per constructed entry.
		//[stream/GfaReadInputStream#005] Check capacity before fetching so the next batch's first line is not discarded.
		for(byte[] line=null; list.size()<maxReadsToReturn && (line=bf.nextLine())!=null;){
			if(line.length>0 && line[0]=='S'){
				lp.set(line);
				String id=lp.parseString(1);
				byte[] bases=lp.parseByteArray(2);
				//[stream/GfaReadInputStream#002] Historical nucleotide test rejected '*' via Read validation (reformat exit1).
				//That backstop depends on constructor validation, input mode and Read settings; it is not a universal GFA check.
				Read r=new Read(bases, null, id, numericID++, flag);
				list.add(r);
			}
		}
		return list;
	}

	/**
	 * Closes the backend and folds its returned error flag into the cached status.
	 * Does not discard buffered reads or reset counters. Thrown exceptions propagate.
	 * @return Cached error status after closure
	 */
	@Override
	public boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		errorState|=tf.close();
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		return errorState;
	}

	/**
	 * Clears counters, buffered data and numeric IDs, then delegates backend reset.
	 * Retains the cached error flag and captured configuration. This is not an atomic
	 * rollback if backend reset fails; coordinate it with iteration and close.
	 */
	@Override
	public synchronized void restart(){
		generated=0;
		consumed=0;
		next=0;
		nextReadID=0;
		buffer=null;
		tf.reset();
	}

	/** @return False; this reader does not construct read pairs */
	@Override
	public boolean paired(){return false;}

	/** @return Input name reported by the backend */
	@Override
	public String fname(){return tf.name();}

	/**
	 * Returns the cached flag, including errors reported by automatic or explicit close.
	 * Thrown parsing failures are separate and need not set this flag.
	 * @return True if a cached error has been recorded
	 */
	@Override
	public boolean errorState(){return errorState;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Prefetched batch; ownership transfers to the caller through nextList. */
	private ArrayList<Read> buffer=null;
	/** Vestigial single-read index; blockwise operations leave it at zero. */
	private int next=0;

	/** Backend opened at construction; it may have its own worker threads. */
	private final ByteFile tf;
	/** Captured amino-acid mask, or zero for nucleotide input. */
	private final int flag;

	/** Read-count batch limit captured from Shared at construction. */
	private final int BUF_LEN=Shared.bufferLen();
	/** Captured but unused byte/base limit; batches currently have no such cap. */
	private final long MAX_DATA=Shared.bufferData();//TODO - lot of work for unlikely case of super-long gfa reads.  Must be disabled for paired-ends.

	/** Reused tab-field parser, accessed during caller-side parsing. */
	private final LineParser1 lp=new LineParser1('\t');

	/** Entries counted after completed refills since construction or restart, not batches. */
	public long generated=0;
	/** Entries handed out by nextList since construction or restart, not batches. */
	public long consumed=0;
	/** Starting numeric ID for the next parse batch; reset by restart. */
	private long nextReadID=0;

	/** Whether the constructor's format descriptor marks standard I/O. */
	public final boolean stdin;

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Enables diagnostic messages; configure before use. */
	public static boolean verbose=false;

}
