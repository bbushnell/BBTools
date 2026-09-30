package stream;

import java.util.ArrayList;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import structures.ListNum;

/**
 * Adapts a ConcurrentReadInputStream to single-entry iteration with batch-local backtracking.
 * Returns entries directly without copying reads or separately iterating their mates.
 * Consumed lists are returned to the input stream; exhaustion closes that stream and
 * folds its error status into errorState. Use this mutable iterator from one thread.
 * <p>Null entries are passed through unchanged. Callers using null as an EOF signal
 * must supply batches without null entries. For early termination, the caller must
 * close the underlying stream; automatic closure occurs only when exhaustion is reached.
 *
 * @author Brian Bushnell
 * @date Jul 18, 2014
 */
public class CrisWrapper{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Opens, starts and primes an input stream without separate quality files.
	 * @param maxReads Read/pair limit forwarded to the input factory
	 * @param keepSamHeader Whether to share this input's SAM/BAM header; false leaves prior shared state alone
	 * @param ff1 Nonnull primary input descriptor
	 * @param ff2 Optional second input descriptor
	 */
	public CrisWrapper(final long maxReads, final boolean keepSamHeader, final FileFormat ff1, final FileFormat ff2){
		this(maxReads, keepSamHeader, ff1, ff2, (String)null, (String)null);
	}

	/**
	 * Opens, starts and primes an input stream, forwarding optional quality paths.
	 * @param maxReads Read/pair limit forwarded to the input factory
	 * @param keepSamHeader Whether to share this input's SAM/BAM header; false leaves prior shared state alone
	 * @param ff1 Nonnull primary input descriptor
	 * @param ff2 Optional second input descriptor
	 * @param qf1 Optional primary quality-file path
	 * @param qf2 Optional second quality-file path
	 */
	public CrisWrapper(final long maxReads, final boolean keepSamHeader, final FileFormat ff1,
		final FileFormat ff2, final String qf1, final String qf2){
		//STR-010: forward the caller's preference; inferring it from the format ignored false for SAM/BAM.
		this(ConcurrentReadInputStream.getReadInputStream(maxReads, keepSamHeader, ff1, ff2, qf1, qf2), true);
	}

	/**
	 * Wraps a nonnull input and immediately fetches its first batch, possibly blocking.
	 * @param cris_ Input stream owned by this iteration until exhaustion
	 * @param start Start the input first; otherwise it must already be ready to supply batches
	 */
	public CrisWrapper(final ConcurrentReadInputStream cris_, final boolean start){
		initialize(cris_, start);
	}

	/**
	 * Replaces the input and primes iteration; does not release any previous input/list.
	 * The caller must finish the old iteration before reusing this wrapper. An initially
	 * empty or null batch is returned when present, then the new input is closed.
	 * The index resets to zero, but accumulated errorState is retained.
	 * @param cris_ Nonnull replacement input
	 * @param start Whether to start the replacement before fetching its first batch
	 */
	public void initialize(final ConcurrentReadInputStream cris_, final boolean start){
		cris=cris_;
		if(start){cris.start();}
		ln=cris.nextList();
		reads=(ln==null ? null : ln.list);
		if(reads==null || reads.size()==0){
			reads=null;
			//#001-fix [stream/CrisWrapper#001]: guard ln!=null. cris.nextList() can return null (shutdown), giving reads=null here; ln.id would then NPE. When ln==null there is no list to return - just close.
			if(ln!=null){cris.returnList(ln.id, true);}
			errorState|=ReadWrite.closeStream(cris);
		}
		index=0;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Returns the next list entry, fetching another batch when needed.
	 * Exhaustion returns null and closes the input; null data entries also return null.
	 * Mates remain attached to their returned Read rather than being visited separately.
	 * @return Next entry, or null for exhaustion or a null entry in the supplied batch
	 */
	public Read next(){
		Read r=null;
		if(reads==null || index>=reads.size()){
			if(reads==null){return null;}
			index=0;
			if(reads.size()==0){
				reads=null;
				cris.returnList(ln.id, true);
				errorState|=ReadWrite.closeStream(cris);
				return null;
			}
			cris.returnList(ln.id, false);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
			if(reads==null){
				//#001-fix [stream/CrisWrapper#001]: ln can be null here (nextList returned null on shutdown); guard before ln.id to avoid NPE. No list to return when ln==null - just close.
				if(ln!=null){cris.returnList(ln.id, true);}
				errorState|=ReadWrite.closeStream(cris);
				return null;
			}
		}
		if(index<reads.size()){
			r=reads.get(index);
			index++;
		}else{//A refill can supply an empty terminal list; one re-entry returns that list and closes the input.
			return next();
		}
		return r;
	}

	/**
	 * Moves back one entry within the current batch; requires index greater than zero.
	 * Repeated calls may revisit earlier entries in that batch, but cannot cross a batch
	 * boundary or restore an exhausted input.
	 */
	public void goBack(){assert(index>0); index--;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Currently borrowed batch, returned when advancing to another batch or EOF. */
	private ListNum<Read> ln;
	/** Current batch contents, or null after exhaustion. */
	private ArrayList<Read> reads;
	/** Index of the next entry in the current batch. */
	private int index;
	/** Wrapped input; do not replace it while iteration is in progress. */
	public ConcurrentReadInputStream cris;
	/** Accumulated errors observed on automatic close; initialize does not reset this flag. */
	public boolean errorState=false;

}
