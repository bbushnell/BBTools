package kmer;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import shared.Tools;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * Linked list implementation of AbstractKmerTable for storing k-mer counts.
 * Each node stores a single k-mer and its count, with overflow handled through
 * chaining to the next node. KmerTable uses these nodes as bucket chains and
 * implements its counting operations separately from the direct-node methods.
 * Structural and count updates require caller coordination. Recursive output
 * and size traversal stop at a nonpositive count, including its following chain.
 *
 * @author Brian Bushnell
 * @date Oct 22, 2013
 */
public class KmerLink extends AbstractKmerTable{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a detached link with count zero and owner -1. */
	public KmerLink(long pivot_){
		pivot=pivot_;
	}

	/** Creates a detached link retaining the supplied scalar count without sign validation. */
	public KmerLink(long pivot_, int value_){
		pivot=pivot_;
		value=value_;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Adds a positive increment and reports structural creation, independent of the count.
	 * An existing zero-count key remains present. Caller must exclude concurrent mutation.
	 * @return 1 for an initialized empty pivot or appended link, otherwise 0 */
	@Override
	public final int incrementAndReturnNumCreated(final long kmer, final int incr){
		//KMR009: find the mutation point once; the resulting count cannot identify creation.
		KmerLink node=this;
		while(node.pivot>=0 && node.pivot!=kmer && node.next!=null){node=node.next;}
		final int created=(node.pivot<0 || node.pivot!=kmer ? 1 : 0);
		node.increment(kmer, incr);
		return created;
	}

	/** Adds a positive increment to a matching link or appends one.
	 * Nonnegative existing counts saturate at Integer.MAX_VALUE on overflow.
	 * @return Resulting count, including the supplied increment for a new key */
	@Override
	public int increment(final long kmer, final int incr){
		if(pivot<0){pivot=kmer; return (value=incr);} //Allows initializing empty nodes to -1
		if(kmer==pivot){
			//KMR019: match KmerTable's saturation when positive addition wraps a nonnegative count.
			if(value<Integer.MAX_VALUE){
				value+=incr;
				if(value<0){value=Integer.MAX_VALUE;}
			}
			return value;
		}
		//KMR009 (kmer/KmerLink#001): a newly appended link contains incr, not necessarily one.
		if(next==null){next=new KmerLink(kmer, incr); return incr;}
		return next.increment(kmer, incr);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Sets count for specified k-mer, creating new node if necessary.
	 * @param kmer The k-mer to set
	 * @param value_ New count value
	 * @return 1 for a new link or initialized empty pivot, otherwise 0
	 */
	@Override
	public int set(long kmer, int value_){
		if(pivot<0){pivot=kmer; value=value_; return 1;} //Allows initializing empty nodes to -1
		if(kmer==pivot){value=value_; return 0;}
		if(next==null){next=new KmerLink(kmer, value_); return 1;}
		return next.set(kmer, value_);
	}

	/**
	 * Sets count for k-mer only if not already present.
	 * @param kmer The k-mer to set
	 * @param value_ Count value to set if k-mer is new
	 * @return 1 for a new link or initialized empty pivot, otherwise 0
	 */
	@Override
	public int setIfNotPresent(long kmer, int value_){
		if(pivot<0){pivot=kmer; value=value_; return 1;} //Allows initializing empty nodes to -1
		if(kmer==pivot){return 0;}
		if(next==null){next=new KmerLink(kmer, value_); return 1;}
		return next.setIfNotPresent(kmer, value_);
	}

	/**
	 * Retrieves the node containing specified k-mer.
	 * @param kmer The k-mer to find
	 * @return KmerLink node containing the k-mer, or null if not found
	 */
	@Override
	KmerLink get(long kmer){
		if(kmer==pivot){return this;}
		return next==null ? null : next.get(kmer);
	}

	/** Appends an existing link when its pivot is absent; retains its following chain too.
	 * The caller must supply a detached, noncyclic chain with suitable distinct pivots. */
	boolean insert(KmerLink n){
		assert(pivot!=-1);
		if(pivot==n.pivot){return false;}
		if(next==null){next=n; return true;}
		return next.insert(n);
	}

	/**
	 * Tests whether specified k-mer is present in the chain.
	 * @param kmer The k-mer to test
	 * @return true if k-mer is present
	 */
	@Override
	public boolean contains(long kmer){
		KmerLink node=get(kmer);
		return node!=null;
	}

	/** Appends live node references tail-to-head; the historical name is retained. */
	void traversePrefix(ArrayList<KmerLink> list){
		if(next!=null){next.traversePrefix(list);}
		list.add(this);
	}

	/** Appends live node references head-to-tail, without count filtering. */
	void traverseInfix(ArrayList<KmerLink> list){
		list.add(this);
		if(next!=null){next.traverseInfix(list);}
	}

	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/

	/** Indicates whether this structure supports resizing.
	 * @return false - linked lists cannot be resized */
	@Override
	boolean canResize(){return false;}

	/** Advertises rebalancing although the retained rebalance entry points currently throw. */
	@Override
	public boolean canRebalance(){return true;}

	/**
	 * Gets array length (unsupported for linked structure).
	 * @return Never returns - throws RuntimeException
	 * @deprecated Linked lists don't have array length
	 */
	@Deprecated
	@Override
	public int arrayLength(){
		throw new RuntimeException("Unsupported.");
	}

	/** Linked storage has no array-resize implementation; always throws. */
	@Deprecated
	@Override
	void resize(){
		throw new RuntimeException("Unsupported.");
	}

	/** Unsupported entry point; the work-list overload is also unimplemented. */
	@Deprecated
	@Override
	public void rebalance(){
		throw new RuntimeException("Please call rebalance(ArrayList<KmerNode>) instead, with an empty list.");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/

	/** Initializes ownership flags for all nodes in chain.
	 * Sets owner to -1 for this node and recursively for all linked nodes. */
	@Override
	public final void initializeOwnership(){
		owner=-1;
		if(next!=null){next.initializeOwnership();}
	}

	/** Clears ownership flags for all nodes in chain.
	 * Delegates to initializeOwnership(). */
	@Override
	public final void clearOwnership(){initializeOwnership();}

	/**
	 * Raises an existing link's owner ID under that link's monitor.
	 * Lookup, the preliminary owner read and the final read are unsynchronized;
	 * callers must separately coordinate structural changes and ownership lifecycle.
	 * @param kmer Existing k-mer whose owner is to be updated
	 * @param newOwner New owner identifier
	 * @return Final owner ID after update attempt
	 */
	@Override
	public final int setOwner(final long kmer, final int newOwner){
		KmerLink n=get(kmer);
		assert(n!=null);
		if(n.owner<=newOwner){
			synchronized(n){
				if(n.owner<newOwner){
					n.owner=newOwner;
				}
			}
		}
		return n.owner;
	}

	/**
	 * Clears owner for k-mer if current owner matches specified owner.
	 * Comparison and replacement share the link's monitor; lookup does not.
	 * @param kmer Existing k-mer whose owner is to be cleared
	 * @param owner Expected current owner ID
	 * @return true if owner was cleared, false if owner didn't match
	 */
	@Override
	public final boolean clearOwner(final long kmer, final int owner){
		KmerLink n=get(kmer);
		assert(n!=null);
		synchronized(n){
			if(n.owner==owner){
				n.owner=-1;
				return true;
			}
		}
		return false;
	}

	/**
	 * Reads an existing link's owner without synchronization.
	 * @param kmer Existing k-mer whose owner is requested
	 * @return Stored owner ID, including -1 for an unowned link
	 */
	@Override
	public final int getOwner(final long kmer){
		KmerLink n=get(kmer);
		assert(n!=null);
		return n.owner;
	}

	/*--------------------------------------------------------------*/

	/**
	 * Sets multiple values for k-mer (unimplemented).
	 * @param kmer The k-mer
	 * @param vals Array of values
	 * @param vlen Length of values
	 * @return Never returns - throws RuntimeException
	 */
	@Override
	public int set(long kmer, int[] vals, int vlen){
		throw new RuntimeException("Unimplemented.");
	}

	/**
	 * Gets count value for specified k-mer.
	 * @param kmer The k-mer to get count for
	 * @return Count value, or -1 if k-mer not found
	 */
	@Override
	public final int getValue(long kmer){
		KmerLink n=get(kmer);
		return n==null ? -1 : n.value;
	}

	/**
	 * Copies the count for a present key into element zero of the supplied array.
	 * @param kmer The k-mer to get count for
	 * @param singleton Caller-owned array with at least one element, required when the key exists
	 * @return singleton array with count, or null if k-mer not found
	 */
	@Override
	public final int[] getValues(long kmer, int[] singleton){
		KmerLink n=get(kmer);
		if(n==null){return null;}
		singleton[0]=n.value;
		return singleton;
	}

	/**
	 * Counts the positive-count prefix, stopping at the first nonpositive node.
	 * @return Number of consecutive positive-count nodes from this link
	 */
	@Override
	public final long size(){
		if(value<1){return 0;}
		long size=1;
		if(next!=null){size+=next.size();}
		return size;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Writes formatted k-mer/count records as bytes, following next links.
	 * Both count bounds are inclusive; a nonpositive node suppresses its following chain.
	 * @param bsw Output stream writer
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @param remaining Optional shared record allowance, decremented before each emitted record
	 * @return Always true; does not report asynchronous writer completion or errors
	 */
	@Override
	public final boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		if(value<1){return true;}
		if(value>=mincount && value<=maxcount){
			if(remaining!=null && remaining.decrementAndGet()<0){return true;}
			bsw.printlnKmer(pivot, value, k);
		}
		if(next!=null){next.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);}
		return true;
	}

	/**
	 * Appends formatted records to an exclusive caller-owned buffer.
	 * Copies a buffer of at least 16000 bytes into a writer job and clears it;
	 * only submission is synchronized. Caller submits the final partial buffer
	 * and prevents concurrent mutation of the chain and buffer.
	 * @param bsw Started writer whose lifecycle remains with the caller
	 * @param bb Exclusive accumulation buffer
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold
	 * @param maxcount Inclusive maximum count threshold
	 * @param remaining Optional shared record allowance, decremented before each emitted record
	 * @return Always true; does not report asynchronous writer completion or errors
	 */
	@Override
	public final boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining){
		if(value<1){return true;}
		//Comprehension [kmer/KmerLink#002 RESOLVED]: dumps now filter on BOTH bounds (value>=mincount && value<=maxcount) in
		//all 3 variants (bytes/_MT/text), matching HashArray/HashForest and the KmerNode1D#001 fix (Brian: "match HashArray").
		//This _MT variant is reachable via KmerTable.dumpKmersAsBytes_MT. Was mincount-only; masked because all callers pass
		//maxcount=Integer.MAX_VALUE, so behavior-neutral today and correct for a real maxcount.
		if(value>=mincount && value<=maxcount){
			if(remaining!=null && remaining.decrementAndGet()<0){return true;}
			toBytes(pivot, value, k, bb);
			bb.nl();
			if(bb.length()>=16000){
				ByteBuilder bb2=new ByteBuilder(bb);
				synchronized(bsw){bsw.addJob(bb2);}
				bb.clear();
			}
		}
		if(next!=null){next.dumpKmersAsBytes_MT(bsw, bb, k, mincount, maxcount, remaining);}
		return true;
	}

	/**
	 * Writes k-mers and counts as text to stream.
	 * @param tsw Text output stream writer
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Inclusive maximum count threshold
	 * @return Always true; does not report asynchronous writer completion or errors
	 */
	@Override
	public final boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
		tsw.print(dumpKmersAsText(new StringBuilder(32), k, mincount, maxcount));
		return true;
	}

	/** Appends the positive-count prefix within both bounds, allocating a buffer only if needed. */
	private final StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount){
		if(value<1){return sb;}
		if(sb==null){sb=new StringBuilder(32);}
		if(value>=mincount && value<=maxcount){sb.append(AbstractKmerTable.toText(pivot, value, k)).append('\n');}
		if(next!=null){next.dumpKmersAsText(sb, k, mincount, maxcount);}
		return sb;
	}

	/**
	 * Fills count histogram array with k-mer count frequencies.
	 * Adds one entry per positive-count link until a nonpositive link ends traversal.
	 * @param ca Histogram array to fill
	 * @param max Highest bin, also receiving counts above max
	 */
	@Override
	public final void fillHistogram(long[] ca, int max){
		if(value<1){return;}
		ca[Tools.min(value, max)]++;
		if(next!=null){next.fillHistogram(ca, max);}
	}

	/**
	 * Fills SuperLongList with individual k-mer count values.
	 * Appends counts until a nonpositive link ends traversal; does not clear the list.
	 * @param sll SuperLongList to add count values to
	 */
	@Override
	public final void fillHistogram(SuperLongList sll){
		if(value<1){return;}
		sll.add(value);
		if(next!=null){next.fillHistogram(sll);}
	}

	/**
	 * Adds each positive-count key's G/C tally once to its count bin.
	 * A nonpositive link ends traversal; existing destination contents are retained.
	 * @param gcCounts Array to accumulate GC counts
	 * @param max Highest count bin, also receiving contributions from counts above max
	 */
	@Override
	public void countGC(long[] gcCounts, int max){
		if(value<1){return;}
		gcCounts[Tools.min(value, max)]+=gc(pivot);
		if(next!=null){next.countGC(gcCounts, max);}
	}


	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Unimplemented work-list rebalancing entry point; always throws. */
	KmerLink rebalance(ArrayList<KmerLink> list){
		throw new RuntimeException("Unsupported.");
	}

	/** Retained unimplemented range helper; always throws. */
	private static KmerLink rebalance(ArrayList<KmerLink> list, int a, int b){
		throw new RuntimeException("Unsupported.");
	}

	/**
	 * Regenerates structure removing zero-value links (unimplemented).
	 * @param limit Regeneration limit
	 * @return Never returns - throws RuntimeException
	 * @deprecated Not yet implemented
	 */
	@Deprecated
	@Override
	public long regenerate(final int limit){
		throw new RuntimeException("TODO - remove zero-value links.");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Encoded key; a negative pivot is used by the direct-node initialization path. */
	long pivot;
	/** Mutable integer count; nonpositive values stop several chain traversals. */
	int value;
	/** Owner marker; -1 is unowned. Plain reads are not synchronized. */
	int owner=-1;
	/** Live next link, or null at the end of the chain. */
	KmerLink next;
}
