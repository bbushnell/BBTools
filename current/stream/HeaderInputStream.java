package stream;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.FileFormat;
import shared.Shared;

/**
 * Reads plain header lines as unpaired, name-only records with sequential numeric IDs.
 * Each line returned by ByteFile is decoded as US-ASCII and passed to Read without
 * reader-side trimming or prefix removal; no sequence or quality array is supplied.
 * Read constructor validation may normalize names under its usual header settings.
 * Uses ByteFile1 directly, independently of the ByteFile factory's backend selection.
 * Callers must coordinate access: synchronized batch methods do not synchronize
 * hasMore, close, external counter reads or shared configuration changes.
 *
 * @author Brian Bushnell
 * @date June 1, 2016
 */
public class HeaderInputStream extends ReadInputStream{

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Prints the first record from the first batch and closes the reader on normal completion.
	 * Requires an input containing a record; the returned close status is not examined.
	 * @param args Command-line arguments; args[0] is the input filename
	 */
	public static void main(String[] args){
		
		HeaderInputStream his=new HeaderInputStream(args[0], true);
		
		Read r=his.nextList().get(0);
		System.out.println(r.toText(false));
		his.close();
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Opens a file of plain header lines with optional subprocess support.
	 * Descriptor detection uses FASTQ as its fallback format; this reader still treats
	 * each returned line as a separate header rather than parsing FASTQ records.
	 * @param fname Input filename
	 * @param allowSubprocess_ Whether subprocess decompression is allowed
	 */
	public HeaderInputStream(String fname, boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.FASTQ, null, allowSubprocess_, false));
	}

	/**
	 * Opens ByteFile1 for the descriptor and records whether it denotes a standard stream.
	 * Batch size is captured from Shared during instance initialization.
	 * @param ff Nonnull descriptor of the input source
	 */
	public HeaderInputStream(FileFormat ff){
		//[stream/HeaderInputStream#001] Name the actual header reader in construction diagnostics.
		if(verbose){System.err.println("HeaderInputStream("+ff+")");}
		
		stdin=ff.stdio();
		
		tf=new ByteFile1(ff);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Checks buffered records and fills a new batch if the source is still open.
	 * Can advance input and generated without handing records to the caller.
	 * With assertions enabled, no buffered records and a closed source require prior nonempty input.
	 * @return Whether a buffered record is available
	 */
	@Override
	public boolean hasMore(){
		if(buffer==null || next>=buffer.size()){
			if(tf.isOpen()){
				fillBuffer();
			}else{
				assert(generated>0) : "Was the file empty?";
			}
		}
		return (buffer!=null && next<buffer.size());
	}
	
	/**
	 * Hands off the next batch, filling it as needed and incrementing consumed.
	 * Returned lists and records are not reused by the reader. This is a blockwise API;
	 * the retained cursor guard requires zero, and no public single-record next method exists.
	 * Batch size limits records rather than total header bytes.
	 * @return Nonempty list of name-only reads, or null when no records remain
	 * @throws RuntimeException If the internal block cursor is nonzero
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
	 * Fills one record-limited batch, advances numeric IDs and closes ByteFile1 on a short batch.
	 * The automatic close result is not merged into this reader's local error flag here.
	 */
	private synchronized void fillBuffer(){
		
		assert(buffer==null || next>=buffer.size());
		
		buffer=null;
		next=0;
		
		buffer=toReadList(tf, BUF_LEN, nextReadID);
		int bsize=(buffer==null ? 0 : buffer.size());
		nextReadID+=bsize;
		//A short buffer (fewer than BUF_LEN reads) means toReadList hit EOF (it only stops early on a null line), so close the file eagerly rather than waiting for a later empty read.
		if(bsize<BUF_LEN){tf.close();}

		generated+=bsize;
		//Currently unreachable: toReadList always returns a (possibly empty) list, never null. Defensive guard kept in case that contract changes.
		if(buffer==null){
			if(!errorState){
				errorState=true;
				//[stream/HeaderInputStream#001] Keep the defensive diagnostic consistent with the reader class.
				System.err.println("Null buffer in HeaderInputStream.");
			}
		}
	}
	
	/**
	 * Closes ByteFile1 and merges its retained status into this reader's local flag.
	 * Does not discard an existing buffered batch. ByteFile1 also returns its retained
	 * status when already closed, including after automatic short-batch close.
	 * @return Accumulated local error status after this close call
	 */
	@Override
	public boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		errorState|=tf.close();
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		return errorState;
	}

	/**
	 * Resets counters, numeric IDs and buffer, then delegates close/reopen to ByteFile1.
	 * Retains this reader's local error flag and captured batch size. Rereading depends
	 * on the underlying source; this method does not guarantee replay of a standard stream.
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
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Reads one name-only record per returned line without closing or resetting the source.
	 * Decodes names with US-ASCII and performs no reader-side trimming, prefix removal
	 * or blank-line filtering. Read constructor header processing may still apply.
	 * Stops at the record limit before consuming the next line, preserving batch boundaries.
	 * @param tf Nonnull byte-line source
	 * @param maxReadsToReturn Positive maximum number of records to return
	 * @param numericID Starting numeric ID, incremented once per record
	 * @return Newly allocated, possibly empty list; each Read has null bases and qualities
	 */
	public static ArrayList<Read> toReadList(ByteFile tf, int maxReadsToReturn, long numericID){
		byte[] line=null;
		ArrayList<Read> list=new ArrayList<Read>(Data.min(8192, maxReadsToReturn));
		int added=0;
		
		for(line=tf.nextLine(); line!=null && added<maxReadsToReturn; line=tf.nextLine()){
			
			Read r=new Read(null, null, new String(line, StandardCharsets.US_ASCII), numericID);

			{
				list.add(r);
				added++;
				numericID++;
			}

			//NOT redundant with the loop guard: without this break, the for-loop's increment would call tf.nextLine() one more time and that line would then fail the guard and be discarded - losing one read at every buffer boundary. The break exits before that extra read is consumed.
			if(added>=maxReadsToReturn){break;}
		}
		assert(list.size()<=maxReadsToReturn);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Metadata           ----------------*/
	/*--------------------------------------------------------------*/

	/** Gets the filename of the input source.
	 * @return Filename from the underlying ByteFile */
	@Override
	public String fname(){return tf.name();}

	/** @return false; this reader does not pair its records */
	@Override
	public boolean paired(){return false;}
	
	/**
	 * Returns only this reader's local flag, without querying ByteFile1 or global parser state.
	 * @return The current local error status
	 */
	@Override
	public boolean errorState(){return errorState;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Generated batch retained until handoff; cleared on restart. */
	private ArrayList<Read> buffer=null;
	/** Legacy block cursor; current batch-only access keeps it at zero. */
	private int next=0;
	
	/** ByteFile1 opened at construction and retained across restart. */
	private final ByteFile tf;

	/** Maximum records per generated batch, captured at instance initialization. */
	private final int BUF_LEN=Shared.bufferLen();
	
	/** Individual records placed in generated batches during the current pass. */
	public long generated=0;
	/** Individual records handed to the caller during the current pass. */
	public long consumed=0;
	/** Numeric ID assigned to the first record of the next generated batch. */
	private long nextReadID=0;
	
	/** Whether the input descriptor denotes a standard stream. */
	public final boolean stdin;
	/** Enables construction and close diagnostics on standard error. */
	public static boolean verbose=false;

}
