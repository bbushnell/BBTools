package stream;

import java.util.ArrayList;
import java.util.Comparator;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import sort.ReadComparatorClump;
import sort.ReadComparatorTopological5Bit;
import structures.ListNum;

/**
 * Batch cursor around a ConcurrentReadInputStream, with the current first read as a comparison key.
 * Construction fetches an initial batch. Subsequent fetches replace that batch and
 * return the previous batch's ID to the stream; this wrapper does not sort a batch
 * or consume individual reads within it. Two exact comparator classes select optional
 * read preparation that can modify Read metadata, including numericID.
 * Callers must coordinate access and follow the wrapped stream's
 * list ownership and lifecycle rules.
 * @author Brian Bushnell
 */
public class CrisContainer implements Comparable<CrisContainer>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Opens a read stream with a FASTQ fallback descriptor, starts it and fetches its first batch.
	 * Read preprocessing is selected by exact comparator class, not subclass membership.
	 * @param fname Input filename to read from
	 * @param comparator_ Comparator used for head comparisons; null permits fetch/peek without comparison
	 * @param allowSubprocess Whether to allow subprocess execution for file access
	 */
	public CrisContainer(String fname, Comparator<Read> comparator_, boolean allowSubprocess){
		comparator=comparator_;
		//FIXED [stream/CrisContainer#001]: null-comparator NPE. Commit 1bf3ec1a refactored these checks from
		//reference-equality (comparator==ReadComparatorTopological5Bit.comparator, null-safe) to comparator.getClass(),
		//which NPEs when comparator is null. Shuffle2 passes null (it shuffles, never sorts), so its temp-file merge
		//died here before writing any reads. Null-guarding restores the pre-refactor tolerance; getClass() is kept
		//(not reverted to ==) because the refactor split each comparator into dual ascending/descending singletons,
		//so reference-equality would miss the descending instance. -GitHub #12 / G11 2026-07-04
		genKmer=(comparator!=null && comparator.getClass()==ReadComparatorTopological5Bit.class);
		clump=(comparator!=null && comparator.getClass()==ReadComparatorClump.class);
		FileFormat ff=FileFormat.testInput(fname, FileFormat.FASTQ, null, allowSubprocess, true);
		cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ff, null, null, null);
		cris.start();
		fetch();
	}
	
	/**
	 * Wraps an existing stream and immediately fetches its first batch without calling start().
	 * Selects preprocessing by exact comparator class, as in the filename constructor.
	 * @param cris_ Nonnull stream already prepared for nextList calls under its lifecycle contract
	 * @param comparator_ Comparator used for head comparisons; null permits fetch/peek without comparison
	 */
	public CrisContainer(ConcurrentReadInputStream cris_, Comparator<Read> comparator_){
		comparator=comparator_;
		//FIXED [stream/CrisContainer#001] (twin): same null-comparator NPE as the (String,...) ctor above. Currently
		//no caller passes null to this 2-arg form, so it's latent — but the null-guard keeps the twins consistent.
		genKmer=(comparator!=null && comparator.getClass()==ReadComparatorTopological5Bit.class);
		clump=(comparator!=null && comparator.getClass()==ReadComparatorClump.class);
		cris=cris_;
		fetch();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Advances to the next batch and returns the previously held list reference.
	 * After construction, call only while hasMore() is true. The previous batch's ID
	 * is returned to the stream during this call, before its list reference is returned
	 * to the caller; do not assume that reference has a new or independent lifetime.
	 * The newly fetched batch receives the configured preprocessing and supplies peek().
	 * @return Previously held list, or null when no list was held, including construction
	 */
	public ArrayList<Read> fetch(){
		final ArrayList<Read> old=list;
		fetchInner();
		return old;
	}
	
	/**
	 * Fetches and normalizes the next batch, prepares its reads and updates the head.
	 * Returns the retained previous ID with a terminal flag derived from the new list,
	 * then records the newly received ID when a wrapper was supplied.
	 */
	private void fetchInner(){
		ListNum<Read> ln=cris.nextList();
		list=(ln==null ? null : ln.list);
		if(list==null || list.size()<1){list=null;}
		else if(genKmer){
			for(Read r : list){ReadComparatorTopological5Bit.genKmer(r);}
		}else if(clump){
			for(Read r : list){ReadComparatorClump.set(r);}
		}
		read=(list==null ? null : list.get(0));
		//One-batch-later return: hand back the previous ID while fetching its successor; a terminal new list sets poison=true.
		//Do not fetch after hasMore()==false: lastNum may still name an ID already returned. Terminal-wrapper rules belong to cris.
		//TODO: Probable bug [stream/CrisContainer#003]: Shuffle2.mergeAndDump removes a container only when fetch()
		//returns null, so it can fetch after hasMore()==false. Reconcile that caller with this stopping rule;
		//the downstream lifecycle effect is unverified and no behavior change is made here.
		//[stream/CrisContainer#001 FIXED 2026-06-18] The null-list guard above preserves EOF handling; the old code dereferenced a possibly-null list.
		if(lastNum>=0){cris.returnList(lastNum, list==null);}
		if(ln!=null){lastNum=ln.id;}
		assert((read==null)==(list==null || list.size()==0));
	}
	
	/**
	 * Delegates stream closure and status reporting to ReadWrite.closeStream.
	 * Does not explicitly clear this wrapper's current head/list fields.
	 * @return true if the wrapped stream reports an error after closing; false otherwise
	 */
	public boolean close(){return ReadWrite.closeStream(cris);}
	
	/**
	 * Returns the first read of the currently fetched batch without advancing.
	 * @return Current comparison head, or null when the latest fetch found no data
	 */
	public Read peek(){return read;}
	
	/**
	 * Compares current heads using this container's comparator, not the other container's.
	 * @param other Nonnull container with a current read; this container also needs a read and comparator
	 * @return Comparator result: negative, zero or positive according to its ordering
	 */
	@Override
	public int compareTo(CrisContainer other){
		assert(read!=null);
		assert(other.read!=null);
		return comparator.compare(read, other.read);
	}
	
	/**
	 * Compares this container's current read to a specific read.
	 * Requires a configured comparator; acceptable read values follow that comparator's contract.
	 * @param other Read to compare against the current head
	 * @return Negative, zero, or positive value indicating relative order
	 */
	public int compareTo(Read other){
		return comparator.compare(read, other);
	}
	
	/**
	 * Checks the current head only; does not fetch or query whether the stream is open.
	 * @return true if the latest fetch supplied a current head, false otherwise
	 */
	public boolean hasMore(){return read!=null;}
	
	/**
	 * Exposes the wrapped stream without transferring or synchronizing ownership.
	 * @return The same stream retained by this container
	 */
	public ConcurrentReadInputStream cris(){return cris;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Wrapped stream, started here only by the filename constructor. */
	final ConcurrentReadInputStream cris;
	/** First read of the latest nonempty batch, or null after a terminal fetch. */
	private Read read;
	/** ID from the most recently received nonnull ListNum wrapper, used when returning the previous batch. */
	private long lastNum=-1;
	/** Latest nonempty batch list, or null after normalization of a terminal result. */
	private ArrayList<Read> list;
	/** Comparator used for head comparisons; may be null when comparisons are unused. */
	private final Comparator<Read> comparator;
	/** genKmer: generate k-mers for reads (topological 5-bit comparator); clump: apply clumping preprocessing.
	 * [stream/CrisContainer#002 DOC FIXED] were two stacked, mis-ordered javadocs; combined to match the genKmer,clump declaration order. */
	private final boolean genKmer, clump;
	
}
