package stream;

import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.FileFormat;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;

//[stream/OnelineReadInputStream#002] Contracts distinguish captured settings, fragment counts and cached error status.
/**
 * Reads Oneline records as a name and sequence separated by the last tab on each line.
 * Names use the platform default charset; copied bases and null qualities are passed
 * through the current Read constructor settings. Blank or tabless lines are not skipped.
 * Pairing is captured from FASTQ.FORCE_INTERLEAVED when the backend is opened.
 * Each returned list entry is an unpaired read or the first read of a reciprocal pair;
 * entry limits, counters and numeric IDs therefore use fragment units.
 * Parsing runs on the caller; callers must coordinate access, lifecycle and global settings.
 * An unmatched interleaved tail terminates the JVM without a cleanup or flushing guarantee.
 * Other parsing failures may throw without setting the cached error flag; this is not
 * an exhaustive format validator.
 * @author Brian Bushnell, Shinobu
 * @date September 29, 2026
 */
public class OnelineReadInputStream extends ReadInputStream{

	/**
	 * Consumes the first batch and prints only its first read, if present, then closes in finally.
	 * Output may precede a close-error rejection; a thrown close can replace an active exception.
	 * Immediate JVM termination during parsing does not guarantee execution of finally.
	 * @param args Input Oneline filename as the first argument
	 * @throws RuntimeException If normal processing completes and input closure reports an error
	 */
	public static void main(final String[] args){
		final OnelineReadInputStream fris=new OnelineReadInputStream(args[0], true);
		boolean error=false;
		try{
			//[stream/OnelineReadInputStream#003] Empty input has no first read; close after an opened reader is used.
			final ArrayList<Read> list=fris.nextList();
			if(list!=null){System.out.println(list.get(0).toText(false));}
		}finally{
			error=fris.close();
		}
		if(error){throw new RuntimeException("Error closing Oneline input.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Resolves an input descriptor with Oneline as the fallback format and opens its backend.
	 * @param fname Input filename or supported standard-input name
	 * @param allowSubprocess_ Whether descriptor resolution permits subprocess-backed input
	 */
	public OnelineReadInputStream(final String fname, final boolean allowSubprocess_){
		this(FileFormat.testInput(fname, FileFormat.ONELINE, null, allowSubprocess_, false));
	}

	/**
	 * Opens the descriptor's ByteFile and captures forced interleaving and current buffer limits.
	 * A non-Oneline descriptor prints a warning but does not change the parser or infer pairing.
	 * @param ff Nonnull input descriptor
	 */
	public OnelineReadInputStream(final FileFormat ff){
		if(verbose){System.err.println("OnelineReadInputStream("+ff+")");}//[stream/OnelineReadInputStream#001] Corrected the former FastqReadInputStream diagnostic name.
		stdin=ff.stdio();
		if(!ff.oneline()){
			System.err.println("Warning: Did not find expected oneline file extension for filename "+ff.name());
		}
		interleaved=FASTQ.FORCE_INTERLEAVED;
		tf=ByteFile.makeByteFile(ff);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Tests for buffered entries, filling a batch if needed while the backend is open.
	 * Prefetch can advance generated before consumed; buffered entries survive close.
	 * @return true if at least one buffered entry remains
	 */
	@Override
	public boolean hasMore(){
		if(buffer==null || next>=buffer.size()){
			//[stream/OnelineReadInputStream#004] Closed zero-record input is ordinary EOF, including repeated queries.
			if(tf.isOpen()){
				fillBuffer();
			}
		}
		return (buffer!=null && next<buffer.size());
	}

	/**
	 * Hands off the buffered list, parsing a batch if none is available, and counts its entries.
	 * The list and its reads are not copied. An already buffered batch can be returned after close.
	 * @return Nonempty list of unpaired reads or first mates, or null for an empty batch
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

	/** Replaces an exhausted buffer, closes on consumed EOF and counts normally produced entries. */
	private synchronized void fillBuffer(){
		assert(buffer==null || next>=buffer.size());
		buffer=null;
		next=0;
		buffer=toReadList();
		final int bsize=(buffer==null ? 0 : buffer.size());
		//[stream/OnelineReadInputStream#005] A data-limited batch can be short; close only after the parser actually reached EOF.
		//[stream/OnelineReadInputStream#006] ByteFile.close returns true on error; retain that status before public close.
		if(batchReachedEOF && tf.close()){errorState=true;}
		generated+=bsize;
		//Defensive: toReadList returns a nonnull list on normal return; nextList maps an empty list to null.
		if(buffer==null){
			if(!errorState){
				errorState=true;
				System.err.println("Null buffer in OnelineReadInputStream.");//#001: corrected the former FastqReadInputStream name.
			}
		}
	}

	/**
	 * Parses complete units until an entry/byte threshold or consumed EOF, assigning fragment IDs.
	 * MAX_DATA counts both mates' sequence bytes and is checked only after a complete unit.
	 * IDs advance during parsing; a later failure can prevent this batch from reaching generated.
	 * A pending first mate at EOF terminates the JVM with a diagnostic; prior output is not rolled back.
	 * @return Newly allocated, possibly empty list; never null on normal return
	 */
	private ArrayList<Read> toReadList(){
		final ArrayList<Read> list=new ArrayList<Read>(BUF_LEN);
		Read r1=null, r2=null;
		long sum=0;
		byte[] line=tf.nextLine();
		for(; line!=null; line=tf.nextLine()){
			//The last tab permits tabs in names; a tabless line fails at the String constructor.
			final int index=Tools.lastIndexOf(line, (byte)'\t');
			final String id=new String(line, 0, index);
			final byte[] bases=KillSwitch.copyOfRange(line, index+1, line.length);
			sum+=bases.length;
			final Read r=new Read(bases, null, id, nextReadID);
			if(r1==null){
				r1=r;
			}else{
				r2=r;
				//[stream/OnelineReadInputStream#007] Read and CRIS require distinct mate numbers: R1=0, R2=1.
				r2.setPairnum(1);
				r1.mate=r2;
				r2.mate=r1;
			}
			//Emit each unpaired read or completed pair. Threshold breaks clear r1 first,
			//so a pending r1 after the loop can only be an unmatched mate at consumed EOF.
			if(interleaved==(r2!=null)){
				list.add(r1);
				r1=r2=null;
				nextReadID++;
				if(list.size()>=BUF_LEN || sum>=MAX_DATA){break;}
			}
		}
		batchReachedEOF=(line==null);
		//[stream/OnelineReadInputStream#008] Use the existing SCARF fatal-input convention instead of silently dropping an orphan.
		if(r1!=null){
			errorState=true;
			KillSwitch.kill("Incomplete interleaved Oneline pair at end of file '"+tf.name()+
				"': read '"+r1.id+"' has no mate (numeric ID "+r1.numericID+").");
		}
		return list;
	}

	/**
	 * Closes the backend and folds its returned error into cached status, retaining buffers/counters.
	 * A thrown backend exception propagates without a new exceptional-completion status latch.
	 * @return true for a previously cached error or a returned close error
	 */
	@Override
	public boolean close(){
		if(verbose){System.err.println("Closing "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		errorState|=tf.close();
		if(verbose){System.err.println("Closed "+this.getClass().getName()+" for "+tf.name()+"; errorState="+errorState);}
		return errorState;
	}

	/**
	 * Clears counters, IDs and the buffer before asking ByteFile to reset to the beginning.
	 * Cached errors and captured configuration remain. The retained batch EOF field is
	 * overwritten by the next normal parse before use; a backend reset exception propagates.
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

	/** @return Forced interleaving captured at construction; not a detection result */
	@Override
	public boolean paired(){return interleaved;}

	/**
	 * Returns only cached local status; does not poll ByteFile or a separate parser's status.
	 * Inline parsing failures can throw without setting this flag.
	 * @return Cached error status retained across close and restart
	 */
	@Override
	public boolean errorState(){return errorState;}

	/** @return Input name reported by the underlying ByteFile */
	@Override
	public String fname(){return tf.name();}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Prefetched entries awaiting handoff; retained by close, cleared by restart. */
	private ArrayList<Read> buffer=null;
	/** Vestigial cursor; this block-only implementation never increments it. */
	private int next=0;
	/** Backend opened at construction. */
	private final ByteFile tf;
	/** Forced pairing captured at construction, independent of the descriptor's pairing hints. */
	private final boolean interleaved;
	/** Recomputed on every normal parser return before fillBuffer decides whether to close. */
	private boolean batchReachedEOF;
	/** Maximum entries per batch captured from Shared; each pair occupies one entry. */
	private final int BUF_LEN=Shared.bufferLen();
	/** Soft sequence-byte budget captured from Shared, tested after complete unpaired reads or pairs. */
	private final long MAX_DATA=Shared.bufferData();
	/** Entries produced since construction/restart, including prefetch not yet handed to the caller. */
	public long generated=0;
	/** Entries handed to callers since construction/restart, not individual mate count. */
	public long consumed=0;
	/** Fragment ID assigned during parsing; a failure can advance it without advancing generated. */
	private long nextReadID=0;
	/** Whether the original descriptor denotes standard input. */
	public final boolean stdin;
	/** Enables construction/closure diagnostics; configure before use. */
	public static boolean verbose=false;
}
