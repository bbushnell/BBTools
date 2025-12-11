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
 * chaining to the next node. Optimized for sparse k-mer distributions where
 * most positions have zero or few k-mers.
 *
 * @author Brian Bushnell
 * @date Oct 22, 2013
 */
public class KmerLink extends AbstractKmerTable {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public KmerLink(long pivot_){
		pivot=pivot_;
	}
	
	public KmerLink(long pivot_, int value_){
		pivot=pivot_;
		value=value_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final int incrementAndReturnNumCreated(final long kmer, final int incr) {
		int x=increment(kmer, incr);
		return x==incr ? 1 : 0;
	}
	
	@Override
	public int increment(final long kmer, final int incr){
		if(pivot<0){pivot=kmer; return (value=incr);} //Allows initializing empty nodes to -1
		if(kmer==pivot){
			if(value<Integer.MAX_VALUE){value+=incr;}
			return value;
		}
		if(next==null){next=new KmerLink(kmer, incr); return 1;}
		return next.increment(kmer, incr);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Sets count for specified k-mer, creating new node if necessary.
	 * @param kmer The k-mer to set
	 * @param value_ New count value
	 * @return Number of nodes added (1 if new, 0 if existing)
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
	 * @return Number of nodes added (1 if new, 0 if already present)
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
	
	void traversePrefix(ArrayList<KmerLink> list){
		if(next!=null){next.traversePrefix(list);}
		list.add(this);
	}
	
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
	boolean canResize() {
		return false;
	}
	
	@Override
	public boolean canRebalance() {
		return true;
	}
	
	/**
	 * Gets array length (unsupported for linked structure).
	 * @return Never returns - throws RuntimeException
	 * @deprecated Linked lists don't have array length
	 */
	@Deprecated
	@Override
	public int arrayLength() {
		throw new RuntimeException("Unsupported.");
	}
	
	@Deprecated
	@Override
	void resize() {
		throw new RuntimeException("Unsupported.");
	}
	
	@Deprecated
	@Override
	public void rebalance() {
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
	 * Sets owner for specified k-mer with thread-safe synchronization.
	 * Only updates if new owner ID is greater than current owner.
	 * @param kmer The k-mer to set owner for
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
	 * Uses synchronization to ensure atomic ownership changes.
	 * @param kmer The k-mer to clear owner for
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
	 * Gets current owner ID for specified k-mer.
	 * @param kmer The k-mer to get owner for
	 * @return Owner ID, or -1 if no owner
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
	public int set(long kmer, int[] vals, int vlen) {
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
	 * Gets count value for k-mer as single-element array.
	 * @param kmer The k-mer to get count for
	 * @param singleton Pre-allocated single-element array to populate
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
	 * Counts total number of valid k-mers in chain.
	 * Recursively counts nodes with value >= 1.
	 * @return Total count of stored k-mers
	 */
	@Override
	public final long size() {
		if(value<1){return 0;}
		long size=1;
		if(next!=null){size+=next.size();}
		return size;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Writes k-mers and counts as binary data to stream.
	 * Only outputs k-mers with counts >= mincount.
	 * @param bsw Output stream writer
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold (unused)
	 * @param remaining Counter for remaining k-mers to output
	 * @return Always true
	 */
	@Override
	public final boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		if(value<1){return true;}
		if(value>=mincount){
			if(remaining!=null && remaining.decrementAndGet()<0){return true;}
			bsw.printlnKmer(pivot, value, k);
		}
		if(next!=null){next.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);}
		return true;
	}
	
	/**
	 * Multi-threaded version of binary k-mer output using ByteBuilder buffer.
	 * Accumulates output in buffer and flushes when buffer reaches 16KB.
	 * @param bsw Thread-safe output stream writer
	 * @param bb Buffer for accumulating output
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold
	 * @param maxcount Maximum count threshold (unused)
	 * @param remaining Counter for remaining k-mers to output
	 * @return Always true
	 */
	@Override
	public final boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining){
		if(value<1){return true;}
		if(value>=mincount){
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
	 * @param maxcount Maximum count threshold (unused)
	 * @return Always true
	 */
	@Override
	public final boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount) {
		tsw.print(dumpKmersAsText(new StringBuilder(32), k, mincount, maxcount));
		return true;
	}
	
	private final StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount){
		if(value<1){return sb;}
		if(sb==null){sb=new StringBuilder(32);}
		if(value>=mincount){sb.append(AbstractKmerTable.toText(pivot, value, k)).append('\n');}
		if(next!=null){next.dumpKmersAsText(sb, k, mincount, maxcount);}
		return sb;
	}
	
	/**
	 * Fills count histogram array with k-mer count frequencies.
	 * Increments histogram bin corresponding to k-mer count value.
	 * @param ca Histogram array to fill
	 * @param max Maximum count value to include in histogram
	 */
	@Override
	public final void fillHistogram(long[] ca, int max){
		if(value<1){return;}
		ca[Tools.min(value, max)]++;
		if(next!=null){next.fillHistogram(ca, max);}
	}
	
	/**
	 * Fills SuperLongList with individual k-mer count values.
	 * Adds each valid k-mer's count to the growing list.
	 * @param sll SuperLongList to add count values to
	 */
	@Override
	public final void fillHistogram(SuperLongList sll){
		if(value<1){return;}
		sll.add(value);
		if(next!=null){next.fillHistogram(sll);}
	}
	
	/**
	 * Counts GC content weighted by k-mer frequencies.
	 * Adds GC count of each k-mer multiplied by its frequency to histogram.
	 * @param gcCounts Array to accumulate GC counts
	 * @param max Maximum count value to include
	 */
	@Override
	public void countGC(long[] gcCounts, int max){
		if(value<1){return;}
		gcCounts[Tools.min(value, max)]+=gc(pivot);
		if(next!=null){next.countGC(gcCounts, max);}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	KmerLink rebalance(ArrayList<KmerLink> list){
		throw new RuntimeException("Unsupported.");
	}
	
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
	
	long pivot;
	int value;
	int owner=-1;
	KmerLink next;
}
