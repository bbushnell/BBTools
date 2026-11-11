package kmer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import shared.Primes;
import shared.Tools;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * Hash table implementation using binary search trees for k-mer storage and counting.
 * Stores k-mers in a prime-sized array where each slot contains a binary search tree.
 * Supports automatic resizing, rebalancing, and both 1D and 2D value storage modes.
 * Also serves as the overflow store for HashArray variants. Count nodes and
 * multi-value nodes have different update semantics; the latter retain value lists.
 * Ordinary lookup and mutation do not acquire the forest lock. Callers must
 * coordinate access, and must not resize or rebalance while traversing live nodes.
 * Keys are nonnegative encoded values; bucket selection uses key modulo capacity.
 *
 * @author Brian Bushnell
 * @date Oct 23, 2013
 */
public final class HashForest extends AbstractKmerTable implements Iterable<KmerNode>{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Allocates an empty forest of single-value nodes with the requested growth policy. */
	public HashForest(int initialSize, boolean autoResize_){
		this(initialSize, autoResize_, false);
	}

	/** Allocates buckets rounded upward to a prime, or one bucket for a request below two.
	 * Node storage is selected by twod_; individual nodes are created on insertion. */
	public HashForest(int initialSize, boolean autoResize_, boolean twod_){
		if(initialSize>1){
			initialSize=(int)Tools.min(maxPrime, Primes.primeAtLeast(initialSize));
		}else{
			initialSize=1;
		}
		prime=initialSize;
		sizeLimit=(long)(initialSize*resizeMult);
		array=allocKmerNodeArray(prime);
		autoResize=autoResize_;
		TWOD=twod_;
	}

	/** Creates a detached scalar-initialized node in the configured storage family. */
	private KmerNode makeNode(long kmer, int val){
		return (TWOD ? new KmerNode2D(kmer, val) : new KmerNode1D(kmer, val));
	}

	/** Creates a multi-value node retaining the supplied array; requires TWOD. */
	private KmerNode makeNode(long kmer, int[] vals, int vlen){
		assert(TWOD);
		return new KmerNode2D(kmer, vals, vlen);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Inserts a scalar-initialized node or delegates increment to its concrete node type.
	 * @return The resulting primary scalar value */
	@Override
	public int increment(final long kmer, final int incr){
		final int cell=(int)(kmer%prime);
		KmerNode n=array[cell], prev=null;
		while(n!=null && n.pivot!=kmer){
			prev=n;
			n=(kmer<n.pivot ? n.left : n.right);
		}
		if(n==null){
			n=makeNode(kmer, incr);
			size++;
			if(prev==null){
				array[cell]=n;
			}else{
				if(kmer<prev.pivot){
					prev.left=n;
				}else{
					prev.right=n;
				}
			}
			if(autoResize && size>sizeLimit){resize();}
		}else{
			n.increment(kmer, incr);
		}
		return n.value();
	}

	/** Tracks actual node creation directly, avoiding the direct-node creation heuristic.
	 * @return 1 when a new node is linked, otherwise 0 */
	@Override
	public int incrementAndReturnNumCreated(final long kmer, final int incr){
		final int cell=(int)(kmer%prime);
		KmerNode n=array[cell], prev=null;
		while(n!=null && n.pivot!=kmer){
			prev=n;
			n=(kmer<n.pivot ? n.left : n.right);
		}
		if(n==null){
			n=makeNode(kmer, incr);
			size++;
			if(prev==null){
				array[cell]=n;
			}else{
				if(kmer<prev.pivot){
					prev.left=n;
				}else{
					prev.right=n;
				}
			}
			if(autoResize && size>sizeLimit){resize();}
			return 1;
		}else{
			n.increment(kmer, incr);
			return 0;
		}
	}

//	public final int set_Test(final long kmer, final int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
//			int[] old=getValues(kmer, null);
//			assert(old==null || contains(kmer, old));
//			x=set0(kmer, v);
//			assert(old==null || contains(kmer, old));
//			assert(contains(kmer, v));
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=set0(kmer, v);
//			assert(contains(kmer, v)) : "old="+old+", v="+v+", kmer="+kmer+", get(kmer)="+getValue(kmer);
//			assert(v==old || !contains(kmer, old));
//		}
//		return x;
//	}
//
//	public final int setIfNotPresent_Test(long kmer, int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
////			int[] vals=getValues(kmer, null);
////			assert(vals==null || contains(kmer, vals));
////			x=setIfNotPresent(kmer, v);
////			assert(contains(kmer, vals));
////			assert(contains(kmer, v));
//			x=0;
//			assert(false);
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=setIfNotPresent0(kmer, v);
//			assert((old<1 && contains(kmer, v)) || (old>0 && contains(kmer, old))) : kmer+", "+old+", "+v;
//		}
//		return x;
//	}
//
//	public final int set_Test(final long kmer, final int v[]){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
//			int[] old=getValues(kmer, null);
//			assert(old==null || contains(kmer, old));
//			x=set0(kmer, v);
//			assert(old==null || contains(kmer, old));
//			assert(contains(kmer, v));
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=set0(kmer, v);
//			assert(contains(kmer, v)) : "old="+old+", v="+v+", kmer="+kmer+", get(kmer)="+getValue(kmer);
//			assert(v[0]==old || !contains(kmer, old));
//		}
//		return x;
//	}


	/** Stores a scalar through the configured node family, updating the cached node count.
	 * Count nodes replace counts; multi-value nodes append according to their own policy. */
	@Override
	public int set(long kmer, int value){
		int x=1, cell=(int)(kmer%prime);
		final KmerNode n=array[cell];
		if(n==null){
			array[cell]=makeNode(kmer, value);
		}else{
			x=n.set(kmer, value);
		}
		size+=x;
		if(autoResize && size>sizeLimit){resize();}
		return x;
	}

	/** Inserts or merges multi-value storage; new nodes retain the array and logical length.
	 * Existing nodes follow KmerNode2D's merge semantics, including its recorded limitations. */
	@Override
	public int set(long kmer, int[] vals, int vlen){
		int x=1, cell=(int)(kmer%prime);
		final KmerNode n=array[cell];
		if(n==null){
			array[cell]=makeNode(kmer, vals, vlen);
		}else{
			x=n.set(kmer, vals, vlen);
		}
		size+=x;
		if(autoResize && size>sizeLimit){resize();}
		return x;
	}

	/**
	 * Sets k-mer value only if k-mer is not already present.
	 * @param kmer The k-mer to set
	 * @param value The value to assign
	 * @return 1 if new node created, 0 if k-mer already existed
	 */
	@Override
	public int setIfNotPresent(long kmer, int value){
		int x=1, cell=(int)(kmer%prime);
		final KmerNode n=array[cell];
		if(n==null){
			array[cell]=makeNode(kmer, value);
		}else{
			x=n.setIfNotPresent(kmer, value);
		}
		size+=x;
		if(autoResize && size>sizeLimit){resize();}
		return x;
	}

	/**
	 * Retrieves the concrete node's primary scalar value associated with a k-mer.
	 * @param kmer The k-mer to look up
	 * @return The associated value, or -1 if not found
	 */
	@Override
	public final int getValue(long kmer){
		int cell=(int)(kmer%prime);
		KmerNode n=array[cell];
		return n==null ? -1 : n.getValue(kmer);
	}

	/**
	 * Retrieves values through the concrete node's storage convention.
	 * @param kmer The k-mer to look up
	 * @param singleton One-element destination in single-value mode; ignored in multi-value mode
	 * @return Supplied destination or live multi-value storage; null for absent/unallocated values
	 */
	@Override
	public int[] getValues(long kmer, int[] singleton){
		int cell=(int)(kmer%prime);
		KmerNode n=array[cell];
		return n==null ? null : n.getValues(kmer, singleton);
	}

	/**
	 * Tests whether the table contains a specific k-mer.
	 * @param kmer The k-mer to test
	 * @return true if k-mer is present, false otherwise
	 */
	@Override
	public boolean contains(long kmer){
		return get(kmer)!=null;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/

	/** Resets every node's owner marker without locking; caller must exclude concurrent access. */
	@Override
	public final void initializeOwnership(){
		for(KmerNode n : array){
			if(n!=null){n.initializeOwnership();}
		}
	}

	/** Clears ownership tracking for all nodes in the table */
	@Override
	public final void clearOwnership(){initializeOwnership();}

	/**
	 * Delegates owner update to an existing node; forest lookup does not acquire a lock.
	 * @param kmer Existing k-mer whose owner is to be updated
	 * @param newOwner The new owner identifier
	 * @return Owner after the update; a higher existing owner takes precedence
	 */
	@Override
	public final int setOwner(final long kmer, final int newOwner){
		final int cell=(int)(kmer%prime);
		KmerNode n=array[cell];
		assert(n!=null);
		return n.setOwner(kmer, newOwner);
	}

	/**
	 * Clears ownership for a k-mer if owned by specified owner.
	 * @param kmer Existing k-mer whose owner is to be cleared
	 * @param owner The expected current owner
	 * @return true if ownership was cleared, false otherwise
	 */
	@Override
	public final boolean clearOwner(final long kmer, final int owner){
		final int cell=(int)(kmer%prime);
		KmerNode n=array[cell];
		assert(n!=null);
		return n.clearOwner(kmer, owner);
	}

	/**
	 * Reads an existing node's owner without synchronizing the lookup or field access.
	 * @param kmer Existing k-mer to query
	 * @return Owner identifier, including -1 for an unowned node
	 */
	@Override
	public final int getOwner(final long kmer){
		final int cell=(int)(kmer%prime);
		KmerNode n=array[cell];
		assert(n!=null);
		return n.getOwner(kmer);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Node Access Methods     ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Retrieves the node containing a specific k-mer.
	 * Borrows the live node; callers must exclude structural mutation while using it.
	 * @param kmer The k-mer to find
	 * @return The node containing the k-mer, or null if not found
	 */
	@Override
	public final KmerNode get(long kmer){
		int cell=(int)(kmer%prime);
		KmerNode n=array[cell];
		while(n!=null && n.pivot!=kmer){
			n=(kmer<n.pivot ? n.left : n.right);
		}
		return n;
	}

	/** Borrows the root of a bucket by array index, or null when that bucket is empty. */
	public final KmerNode getNode(int cell){
		KmerNode n=array[cell];
		return n;
	}

	/** Clears a node's child links and reinserts that node without adjusting size.
	 * Rehash callers must first retain every descendant separately and preserve unique pivots. */
	boolean insert(KmerNode n){
		n.left=null;
		n.right=null;
		int cell=(int)(n.pivot%prime);
		if(array[cell]==null){
			array[cell]=n;
			return true;
		}
		return array[cell].insert(n);
	}

	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns whether this table implementation supports resizing */
	@Override
	boolean canResize(){return true;}

	/** Returns whether this table implementation supports rebalancing */
	@Override
	public boolean canRebalance(){return true;}

	/** Returns the cached number of linked nodes, without filtering their scalar values. */
	@Override
	public long size(){return size;}

	/** Returns the length of the underlying hash array */
	@Override
	public int arrayLength(){return array.length;}

	/** Chooses a larger prime bucket count within load heuristics, then relinks existing nodes.
	 * May leave capacity unchanged. Synchronizes resize calls only; ordinary access still
	 * requires the caller's exclusion, and borrowed bucket arrays become stale after growth. */
	@Override
	synchronized void resize(){
//		assert(false);
//		System.err.println("Resizing from "+prime+"; load="+(size*1f/prime));
		sizeLimit=Tools.max((long)(size*1.4), (long)(maxLoadFactor*prime));

		final long maxAllowedByLoadFactor=(long)(size*minLoadMult);
		final long minAllowedByLoadFactor=(long)(size*maxLoadMult);
		assert(maxAllowedByLoadFactor>=minAllowedByLoadFactor);
		if(maxAllowedByLoadFactor<prime){return;}

		long x=10+(long)(prime*resizeMult);
		x=Tools.max(x, minAllowedByLoadFactor);
		x=Tools.min(x, maxAllowedByLoadFactor);

		int prime2=(int)Tools.min(maxPrime, Primes.primeAtLeast(x));

		if(prime2<=prime){return;}

		prime=prime2;
//		System.err.println("Resized to "+prime+"; load="+(size*1f/prime));
		KmerNode[] old=array;
		array=allocKmerNodeArray(prime2);
		ArrayList<KmerNode> list=new ArrayList<KmerNode>(1000);
		for(int i=0; i<old.length; i++){
			if(old[i]!=null){
				old[i].traverseInfix(list);
				for(KmerNode n : list){insert(n);}
				list.clear();
			}
		}
		sizeLimit=Tools.max((long)(size*1.4), (long)(maxLoadFactor*prime));
	}

	/** Rebalances each bucket using a reused work list and installs each returned root. */
	@Override
	public void rebalance(){
		ArrayList<KmerNode> list=new ArrayList<KmerNode>(1000);
		for(int i=0; i<array.length; i++){
			if(array[i]!=null){array[i]=array[i].rebalance(list);}
		}
	}

	/** Drops all bucket roots and resets size, retaining the allocated bucket array. */
	public void clear(){
		size=0;
		Arrays.fill(array, null);
	}

	/**
	 * Regenerate operation not implemented for HashForest.
	 * @param limit Regeneration limit (unused)
	 * @throws RuntimeException Always thrown as operation not supported
	 */
	@Override
	long regenerate(final int limit){
		throw new RuntimeException("Not implemented.");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Writes k-mers and counts as text to the specified writer.
	 * Delegates traversal and filtering to each concrete node family.
	 * @param tsw Writer for text output
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @return Always true; writer completion/errors are not checked here
	 */
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
//		tsw.print("HashForest:\n");
		for(int i=0; i<array.length; i++){
			KmerNode node=array[i];
			if(node!=null){
				node.dumpKmersAsText(tsw, k, mincount, maxcount);
			}
		}
		return true;
	}

	/**
	 * Delegates formatted byte output to every nonnull bucket root.
	 * @param bsw Borrowed writer for formatted output
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @param remaining Optional shared record allowance forwarded to node traversals
	 * @return Always true; writer completion/errors are not checked here
	 */
	@Override
	public boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		//The outer minimum-count gate and extra remaining decrement were removed in0e2547dc: a positive root
		//below mincount can still have qualifying descendants. Leave filtering to the concrete node traversal;
		//its separate nonpositive-value subtree pruning remains unchanged. Do not re-add an outer count gate.
		for(int i=0; i<array.length; i++){
			KmerNode node=array[i];
			if(node!=null){
				node.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);
			}
		}
		return true;
	}

	/**
	 * Delegates formatted byte output using an exclusive caller-owned buffer.
	 * Node methods submit full copied buffers; caller submits the final partial buffer.
	 * This method starts no workers and acquires no forest lock.
	 * @param bsw Borrowed, started writer
	 * @param bb Exclusive buffer for building output
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @param remaining Optional shared record allowance forwarded to node traversals
	 * @return Always true; writer completion/errors are not checked here
	 */
	@Override
	public boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining){
		for(int i=0; i<array.length; i++){
			KmerNode node=array[i];
			if(node!=null){
				node.dumpKmersAsBytes_MT(bsw, bb, k, mincount, maxcount, remaining);
			}
		}
		return true;
	}

	/**
	 * Adds primary-value histogram counts using node-specific subtree pruning.
	 * @param ca Existing accumulation array, not cleared
	 * @param max Highest count bin, also receiving higher counts
	 */
	@Override
	public void fillHistogram(long[] ca, int max){
		for(int i=0; i<array.length; i++){
			KmerNode node=array[i];
			if(node!=null){
				node.fillHistogram(ca, max);
			}
		}
	}

	/** Appends primary values using node-specific subtree pruning.
	 * @param sll Existing accumulation list, not cleared */
	@Override
	public final void fillHistogram(SuperLongList sll){
		for(int i=0; i<array.length; i++){
			KmerNode node=array[i];
			if(node!=null){
				node.fillHistogram(sll);
			}
		}
	}

	/** Adds per-key G/C tallies to primary-value bins through each node's traversal. */
	@Override
	public void countGC(long[] gcCounts, int max){
		for(int i=0; i<array.length; i++){
			KmerNode node=array[i];
			if(node!=null){
				node.countGC(gcCounts, max);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Iteration           ----------------*/
	/*--------------------------------------------------------------*/

	/** Allocates a list of live node references and returns that list's iterator.
	 * Iterator removal affects the temporary list, not this forest. */
	@Override
	public Iterator<KmerNode> iterator(){
		return toList().iterator();
	}

	/** Collects every node in bucket order and root-first subtree order, without value filtering.
	 * The list is independent, but its nodes remain shared and mutable. */
	public ArrayList<KmerNode> toList(){
		assert(size<Integer.MAX_VALUE);
		ArrayList<KmerNode> list=new ArrayList<KmerNode>((int)size);
		for(int i=0; i<array.length; i++){
			if(array[i]!=null){array[i].traverseInfix(list);}
		}
		assert(list.size()==size);
		return list;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrows the bucket array; caller must not replace roots independently or retain it across resizing. */
	public KmerNode[] array(){return array;}

	/** Current bucket roots; replaced when capacity grows. */
	KmerNode[] array;
	/** Current modulus/bucket-array length, normally prime (one is also permitted). */
	int prime;
	/** Cached count of linked entries, maintained by public insertion paths. */
	long size=0;
	/** Entry-count threshold that triggers an automatic resize attempt. */
	long sizeLimit;
	/** Whether public insertion paths attempt growth above sizeLimit. */
	final boolean autoResize;
	/** Selects multi-value nodes instead of scalar count nodes. */
	final boolean TWOD;
	/** Coordination lock exposed to callers; ordinary methods do not acquire it automatically. */
	private final Lock lock=new ReentrantLock();

	/** Returns the caller-managed coordination lock. */
	@Override
	final Lock getLock(){return lock;}

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Largest prime representable by the bucket-index type. */
	final static int maxPrime=(int)Primes.primeAtMost(Integer.MAX_VALUE);
	/** Preferred capacity growth multiplier before load-factor bounds are applied. */
	final static float resizeMult=2.5f;
	/** Target minimum entries per bucket after growth. */
	final static float minLoadFactor=0.75f;
	/** Target maximum entries per bucket after growth. */
	final static float maxLoadFactor=2.5f;
	/** Reciprocal minimum load, used to bound capacity above. */
	final static float minLoadMult=1/minLoadFactor;
	/** Reciprocal maximum load, used to bound capacity below. */
	final static float maxLoadMult=1/maxLoadFactor;



}
