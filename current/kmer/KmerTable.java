package kmer;

import java.util.ArrayList;
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
 * Hash table implementation for storing k-mer counts using separate chaining.
 * Uses a prime-sized array with linked lists to handle collisions.
 * Ordinary methods do not acquire the provided table lock; callers coordinate
 * concurrent access and exclude resizing while using live links. Keys must be
 * nonnegative encoded values. Several legacy traversal operations remain unsupported.
 *
 * @author Brian Bushnell
 * @date Oct 23, 2013
 */
public final class KmerTable extends AbstractKmerTable{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Allocates empty buckets rounded upward to a prime, or one bucket for requests below two. */
	public KmerTable(int initialSize, boolean autoResize_){
		if(initialSize>1){
			initialSize=(int)Tools.min(maxPrime, Primes.primeAtLeast(initialSize));
		}else{
			initialSize=1;
		}
		prime=initialSize;
		sizeLimit=(long)(initialSize*resizeMult);
		array=new KmerLink[prime];
		autoResize=autoResize_;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Updates count fields directly, clamping negative sums to Integer.MAX_VALUE.
	 * Does not use the separate direct-link increment implementation.
	 * @return Current count after insertion or update */
	@Override
	public int increment(final long kmer, final int incr){
		final int cell=(int)(kmer%prime);
		KmerLink n=array[cell], prev=null;
		while(n!=null && n.pivot!=kmer){
			prev=n;
			n=n.next;
		}
		if(n==null){
			n=new KmerLink(kmer, incr);
			size++;
			if(prev==null){
				array[cell]=n;
			}else{
				prev.next=n;
			}
			if(autoResize && size>sizeLimit){resize();}
		}else{
			n.value+=incr;
			if(n.value<0){n.value=Integer.MAX_VALUE;}
		}
		return n.value;
	}

	/**
	 * Increments a k-mer count and returns the number of new entries created.
	 * Used for tracking how many new k-mers are added during batch operations.
	 *
	 * @param kmer The k-mer to increment
	 * @param incr Amount to increment by (must be positive)
	 * @return 1 if a new entry was created, 0 if entry already existed
	 */
	@Override
	public int incrementAndReturnNumCreated(final long kmer, final int incr){
		final int cell=(int)(kmer%prime);
		KmerLink n=array[cell], prev=null;
		while(n!=null && n.pivot!=kmer){
			prev=n;
			n=n.next;
		}
		if(n==null){
			n=new KmerLink(kmer, incr);
			size++;
			if(prev==null){
				array[cell]=n;
			}else{
				prev.next=n;
			}
			if(autoResize && size>sizeLimit){resize();}
			return 1;
		}else{
			n.value+=incr;
			if(n.value<0){n.value=Integer.MAX_VALUE;}
			return 0;
		}
	}

	/**
	 * Sets the count for a k-mer to a specific value.
	 * Creates a new entry if the k-mer is not present in the table.
	 *
	 * @param kmer The k-mer to set
	 * @param value The count value to set
	 * @return Number of new entries created (1 or 0)
	 */
	@Override
	public int set(long kmer, int value){
		int x=1, cell=(int)(kmer%prime);
		final KmerLink n=array[cell];
		if(n==null){
			array[cell]=new KmerLink(kmer, value);
		}else{
			x=n.set(kmer, value);
		}
		size+=x;
		if(autoResize && size>sizeLimit){resize();}
		return x;
	}

	/**
	 * Sets multiple values for a k-mer.
	 * Currently unimplemented in this table type.
	 *
	 * @param kmer The k-mer to set
	 * @param vals Array of values to set
	 * @param vlen Number of valid values in the array
	 * @return Never returns normally
	 * @throws RuntimeException Always thrown as method is unimplemented
	 */
	@Override
	public int set(long kmer, int[] vals, int vlen){
		throw new RuntimeException("Unimplemented.");
	}

	/**
	 * Sets the count for a k-mer only if it is not already present.
	 * Does nothing if the k-mer already exists in the table.
	 *
	 * @param kmer The k-mer to set
	 * @param value The count value to set
	 * @return Number of new entries created (1 or 0)
	 */
	@Override
	public int setIfNotPresent(long kmer, int value){
		int x=1, cell=(int)(kmer%prime);
		final KmerLink n=array[cell];
		if(n==null){
			array[cell]=new KmerLink(kmer, value);
		}else{
			x=n.setIfNotPresent(kmer, value);
		}
		size+=x;
		if(autoResize && size>sizeLimit){resize();}
		return x;
	}

	/**
	 * Retrieves the count value for a k-mer.
	 * @param kmer The k-mer to look up
	 * @return The stored count (including zero), or NOT_PRESENT if the key is absent
	 */
	@Override
	public int getValue(long kmer){
		//KMR006 [kmer/KmerTable#001]: match the abstract absence sentinel without changing stored zero.
		int cell=(int)(kmer%prime);
		KmerLink n=array[cell];
		while(n!=null && n.pivot!=kmer){n=n.next;}
		return n==null ? NOT_PRESENT : n.value;
	}

	/**
	 * Copies getValue's scalar result into the supplied array, including NOT_PRESENT for absence.
	 *
	 * @param kmer The k-mer to look up
	 * @param singleton Nonnull caller-owned array with at least one element
	 * @return The singleton array with the value filled in
	 */
	@Override
	public int[] getValues(long kmer, int[] singleton){
		//KMR005 [kmer/KmerTable#002]: validate the destination, not the nonempty bucket array.
		assert(singleton!=null && singleton.length>0) : "getValues copies one scalar to singleton[0]; length="+
			(singleton==null ? "null" : singleton.length);
		singleton[0]=getValue(kmer);
		return singleton;
	}

	/**
	 * Tests whether a k-mer is present in the table.
	 * @param kmer The k-mer to check
	 * @return true if the k-mer exists in the table, false otherwise
	 */
	@Override
	public boolean contains(long kmer){
		KmerLink node=get(kmer);
		return node!=null;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/

	/** Initializes ownership tracking data structures for all k-mers in the table.
	 * Resets all owner markers to -1 without locking; callers exclude concurrent access. */
	@Override
	public final void initializeOwnership(){
		for(KmerLink n : array){
			if(n!=null){n.initializeOwnership();}
		}
	}

	/** Clears ownership information for all k-mers in the table.
	 * Resets all ownership values to the unowned state. */
	@Override
	public final void clearOwnership(){
		for(KmerLink n : array){
			if(n!=null){n.clearOwnership();}
		}
	}

	/**
	 * Delegates owner update to an existing link; bucket lookup does not acquire a lock.
	 *
	 * @param kmer Existing k-mer whose owner is to be updated
	 * @param newOwner Thread ID of the new owner
	 * @return The actual owner ID after the operation
	 */
	@Override
	public final int setOwner(final long kmer, final int newOwner){
		final int cell=(int)(kmer%prime);
		KmerLink n=array[cell];
		assert(n!=null);
		return n.setOwner(kmer, newOwner);
	}

	/**
	 * Clears an existing link's owner if it matches the supplied ID.
	 * Comparison/replacement use the link monitor; lookup itself is unlocked.
	 *
	 * @param kmer Existing k-mer whose owner is to be cleared
	 * @param owner Thread ID of the supposed owner
	 * @return true if ownership was cleared, false if caller was not the owner
	 */
	@Override
	public final boolean clearOwner(final long kmer, final int owner){
		final int cell=(int)(kmer%prime);
		KmerLink n=array[cell];
		assert(n!=null);
		return n.clearOwner(kmer, owner);
	}

	/**
	 * Reads an existing link's owner without synchronizing lookup or field access.
	 * @param kmer Existing k-mer whose owner is requested
	 * @return Thread ID of the owner, or -1 if unowned
	 */
	@Override
	public final int getOwner(final long kmer){
		final int cell=(int)(kmer%prime);
		KmerLink n=array[cell];
		assert(n!=null);
		return n.getOwner(kmer);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Retrieves the KmerLink node containing a specific k-mer.
	 * Used internally for direct node access.
	 * @param kmer The k-mer to find
	 * @return The KmerLink node, or null if not found
	 */
	@Override
	KmerLink get(long kmer){
		int cell=(int)(kmer%prime);
		KmerLink n=array[cell];
		while(n!=null && n.pivot!=kmer){n=n.next;}
		return n;
	}

	/** Detaches a link from its successor and rehashes it without adjusting size.
	 * Caller must first retain all chain nodes separately and preserve unique pivots. */
	boolean insert(KmerLink n){
		n.next=null;
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

	/** Returns true as this table type supports dynamic resizing */
	@Override
	boolean canResize(){return true;}

	/** Linked buckets have no supported rebalance operation. */
	@Override
	public boolean canRebalance(){return false;}

	/** Returns the cached linked-entry count without filtering nonpositive values. */
	@Override
	public long size(){return size;}

	/** Returns the length of the internal hash array */
	@Override
	public int arrayLength(){return array.length;}

	/**
	 * Selects a larger prime within load heuristics and relinks existing entries.
	 * May leave capacity unchanged. The monitor serializes resize calls only;
	 * callers must separately exclude ordinary access during relinking.
	 */
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
		KmerLink[] old=array;
		array=new KmerLink[prime2];
		ArrayList<KmerLink> list=new ArrayList<KmerLink>(1000);
		for(int i=0; i<old.length; i++){
			if(old[i]!=null){
				old[i].traverseInfix(list);
				for(KmerLink n : list){insert(n);}
				list.clear();
			}
		}
		sizeLimit=Tools.max((long)(size*1.4), (long)(maxLoadFactor*prime));
	}

	/** Retained unsupported traversal: empty tables return, but the first nonempty chain throws. */
	@Override
	public void rebalance(){
		ArrayList<KmerLink> list=new ArrayList<KmerLink>(1000);
		for(int i=0; i<array.length; i++){
			if(array[i]!=null){array[i]=array[i].rebalance(list);}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Dumps k-mers as text format to a stream writer.
	 *
	 * @param tsw The text stream writer to write to
	 * @param k Length of k-mers
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @return Never returns normally
	 * @throws RuntimeException Method is deprecated and unimplemented
	 */
	@Deprecated
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
		throw new RuntimeException("TODO");
	}

	/**
	 * Delegates byte-formatted output to bucket chains using a caller-owned buffer.
	 * Starts no workers here. Visits every nonnull bucket head so a positive head
	 * below mincount cannot suppress a qualifying tail (KMR007). Delegated links
	 * retain their own count bounds and nonpositive-node pruning rules.
	 *
	 * @param bsw The byte stream writer to write to
	 * @param bb Exclusive buffer; caller submits the final partial buffer
	 * @param k Length of k-mers
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @param remaining Atomic counter for limiting output
	 * @return Always true; asynchronous writer completion/errors are not checked here
	 */
	@Override
	public boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining){
		for(int i=0; i<array.length; i++){
			KmerLink node=array[i];
			//KMR007: each link applies count bounds; a below-threshold head does not exclude its tail.
			if(node!=null){
				node.dumpKmersAsBytes_MT(bsw, bb, k, mincount, maxcount, remaining);
			}
		}
		return true;
	}

	/**
	 * Dumps k-mers as bytes to a stream writer.
	 *
	 * @param bsw The byte stream writer to write to
	 * @param k Length of k-mers
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @param remaining Atomic counter for limiting output
	 * @return Never returns normally
	 * @throws RuntimeException Method is deprecated and unimplemented
	 */
	@Deprecated
	@Override
	public boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		throw new RuntimeException("TODO");
	}

	/**
	 * Fills a histogram array with k-mer count distribution.
	 * @param ca Array to fill with histogram counts
	 * @param max Maximum count value to include
	 * @throws RuntimeException Method is deprecated and unimplemented
	 */
	@Deprecated
	@Override
	public void fillHistogram(long[] ca, int max){
		throw new RuntimeException("TODO");
	}

	/** Unsupported list histogram operation; always throws. */
	@Deprecated
	@Override
	public void fillHistogram(SuperLongList sll){
		throw new RuntimeException("TODO");
	}

	/**
	 * Counts GC content distribution across k-mers in the table.
	 * @param gcCounts Array to store GC count results
	 * @param max Maximum count value to process
	 * @throws RuntimeException Method is deprecated and unimplemented
	 */
	@Deprecated
	@Override
	public void countGC(long[] gcCounts, int max){
		throw new RuntimeException("TODO");
	}

	/**
	 * Unsupported removal of entries at or below a count threshold.
	 * @param limit Count threshold (unused)
	 * @return Never returns normally
	 * @throws RuntimeException Always; removal of links is not implemented
	 */
	@Deprecated
	@Override
	public long regenerate(final int limit){
		throw new RuntimeException("TODO - remove zero-value links.");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Bucket heads; replaced on capacity growth. */
	KmerLink[] array;
	/** Current bucket-array length and modulus; one is permitted for a tiny table. */
	int prime;
	/** Cached count of linked entries, maintained by public insertion paths. */
	long size=0;
	/** Entry-count threshold for an automatic resize attempt. */
	long sizeLimit;
	/** Whether public insertion paths attempt growth above sizeLimit. */
	final boolean autoResize;
	/** Caller-managed coordination lock; ordinary methods do not acquire it. */
	private final Lock lock=new ReentrantLock();

	/** Returns the caller-managed coordination lock. */
	@Override
	final Lock getLock(){return lock;}

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Largest prime representable by the bucket-index type. */
	final static int maxPrime=(int)Primes.primeAtMost(Integer.MAX_VALUE);
	/** Preferred growth multiplier before load bounds are applied. */
	final static float resizeMult=2f;
	/** Target minimum entries per bucket after growth. */
	final static float minLoadFactor=0.5f;
	/** Target maximum entries per bucket after growth. */
	final static float maxLoadFactor=0.98f;
	/** Reciprocal minimum load, used to bound capacity above. */
	final static float minLoadMult=1/minLoadFactor;
	/** Reciprocal maximum load, used to bound capacity below. */
	final static float maxLoadMult=1/maxLoadFactor;

}
