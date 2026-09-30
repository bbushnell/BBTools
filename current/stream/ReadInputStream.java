package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.KillSwitch;
import structures.ListNum;

/**
 * Base API for batch sequence readers, with static helpers that materialize input.
 * Concrete readers define source handling, restart support and synchronization;
 * this base class does not make instance methods safe for concurrent callers.
 * Static loading helpers use ConcurrentReadInputStream and retain Read references,
 * including attached mates, without copying the Read objects.
 * @author Brian Bushnell
 */
public abstract class ReadInputStream{

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Resolves an input descriptor and delegates to the descriptor-based loader.
	 * A null filename returns immediately without opening an input.
	 * @param fname Input filename, or null
	 * @param defaultFormat Fallback format supplied to FileFormat.testInput
	 * @param maxReads Limit forwarded unchanged; reader-defined units, negative for unlimited
	 * @return Collected entries, possibly empty; null only when fname is null
	 * @see #toReads(FileFormat, long)
	 */
	public static final ArrayList<Read> toReads(final String fname, final int defaultFormat, final long maxReads){
		if(fname==null){return null;}
		final FileFormat ff=FileFormat.testInput(fname, defaultFormat, null, false, true);
		return toReads(ff, maxReads);
	}

	/**
	 * Materializes the descriptor-based loader's entries in a new array.
	 * Read objects and their attached mates are not copied.
	 * @param ff Nonnull primary input descriptor
	 * @param maxReads Limit forwarded unchanged; reader-defined units, negative for unlimited
	 * @return New array, possibly empty; never null on normal return
	 * @see #toReads(FileFormat, long)
	 */
	public static final Read[] toReadArray(final FileFormat ff, final long maxReads){
		final ArrayList<Read> list=toReads(ff, maxReads);
		//#001-fix [stream/ReadInputStream#001]: toReads(ff) never returns null, so this returns a (possibly empty) array, never null. Load-bearing: prok/ProkObject.stripOrganelle iterates array.length with no null guard. Javadoc corrected ("or null if none"->"empty array if none").
		return list==null ? null : list.toArray(new Read[0]);
	}
	
	/**
	 * Starts a ConcurrentReadInputStream with no second input or header-retention request,
	 * collects entries until a null or empty batch, returns borrowed batches and closes it.
	 * Entries are retained in stream order; attached mates are not separately appended.
	 * A reported final input error terminates the process instead of returning potentially
	 * incomplete data, including when assertions are disabled. Process termination does
	 * not guarantee other resource cleanup. Other exceptions propagate without a
	 * finally-based close guarantee; this helper does not validate every input property.
	 * @param ff Nonnull primary input descriptor
	 * @param maxReads Limit forwarded unchanged; reader-defined units, negative for unlimited
	 * @return Newly allocated list, possibly empty; never null on normal return
	 */
	public static final ArrayList<Read> toReads(final FileFormat ff, final long maxReads){
		final ArrayList<Read> list=new ArrayList<Read>();

		/* Start an input stream */
		final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, ff, null);
		cris.start();
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);

		/* Iterate through read lists from the input stream */
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			list.addAll(reads);

			/* Dispose of the old list and fetch a new one */
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		/* Cleanup */
		//Each ListNum is returned exactly once: the loop returns every non-terminal list, this returns the terminal poison/empty buffer (returnList tolerates null). list is never null (may be empty).
		cris.returnList(ln);
		//STR-005: FastqReadInputStream.close folds in the byte reader's latched I/O error here.
		//Match the checked CRIS helper: reported input failures must not become successful reference loads.
		if(ReadWrite.closeStream(cris)){
			KillSwitch.kill("Error: an input error was reported reading "+ff.name()+
				" during ReadInputStream.toReads(); aborting rather than returning potentially incomplete read data.");
		}
		return list;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Retrieves the next batch of read entries from this input.
	 * Buffering, blocking and any synchronization are implementation-specific.
	 * @return Next batch, or null when no more entries remain
	 */
	public abstract ArrayList<Read> nextList();

	/** Reports whether more entries may be available; implementations may read ahead. */
	public abstract boolean hasMore();

	/**
	 * Requests a restart from the beginning of the source, when supported.
	 * Buffer, counter and error-state reset behavior belongs to the concrete reader.
	 * @throws RuntimeException if the implementation does not support restarting
	 */
	public abstract void restart();

	/**
	 * Closes the input and releases its resources.
	 * @return Reported error status after close, which may include earlier read errors
	 */
	public abstract boolean close();

	/**
	 * Indicates whether this reader supplies paired-end data.
	 * @return True when the concrete reader reports paired input
	 */
	public abstract boolean paired();

	/**
	 * Copies array entries, including null entries, into a new mutable list.
	 * @param array Source array; the Read objects themselves are not copied
	 * @return New list in array order, or null for a null or zero-length array
	 */
	protected static final ArrayList<Read> toList(final Read[] array){
		if(array==null || array.length==0){return null;}
		final ArrayList<Read> list=new ArrayList<Read>(array.length);
		for(int i=0; i<array.length; i++){list.add(array[i]);}
		return list;
	}

	/** Returns this instance's stored error flag; subclasses may also consult their inputs. */
	public boolean errorState(){return errorState;}

	/**
	 * Returns the concrete reader's source identifier, such as a filename.
	 * @return Source identifier, or null if unavailable
	 */
	public abstract String fname();

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Local error flag; concrete readers decide when to set or reset it. */
	protected boolean errorState=false;

}
