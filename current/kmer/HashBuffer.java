package kmer;

import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * Buffered kmer table that distributes kmers across multiple underlying tables.
 * Uses reusable buffers to batch updates to shared in-memory tables.
 * Kmers are distributed to different tables based on hash values for load balancing.
 * Each worker needs its own HashBuffer. Flushes acquire the destination table's
 * lock; reads and ownership operations delegate without flushing or locking, so
 * callers must coordinate them with writes. Pending entries are invisible to reads.
 * Count mode accepts incrementAndReturnNumCreated with increment one; value mode
 * accepts set. Do not mix these modes. Call flush explicitly after normal loading;
 * no automatic final flush or transactional retry is provided.
 *
 * @author Brian Bushnell
 * @date Nov 22, 2013
 */
public class HashBuffer extends AbstractKmerTable{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrows nonempty backend tables and creates private buffers with two-bit DNA routing geometry.
	 * @param tables_ Backend tables implementing the requested count/value operations
	 * @param buflen_ Positive flush threshold in keys
	 * @param k_ Key length used for two-bit routing masks and diagnostics
	 * @param initValues Selects value-setting rather than count-increment buffers
	 * @param setIfNotPresent_ Uses conditional insertion when flushing value buffers */
	public HashBuffer(AbstractKmerTable[] tables_, int buflen_, int k_, boolean initValues, boolean setIfNotPresent_){
		this(tables_, buflen_, k_, initValues, setIfNotPresent_, false);
	}

	/** Creates private buffers using the key alphabet's routing geometry.
	 * Raw keys are retained unchanged; callers supply the same encoding and routing
	 * settings used by their query-side table set.
	 * @param tables_ Borrowed nonempty backend shards
	 * @param buflen_ Positive flush threshold in keys
	 * @param k_ Bases or amino-acid residues per key
	 * @param initValues Selects value-setting rather than count-increment buffers
	 * @param setIfNotPresent_ Uses conditional insertion when flushing value buffers
	 * @param amino Whether keys use five bits per residue rather than two bits per DNA base */
	public HashBuffer(AbstractKmerTable[] tables_, int buflen_, int k_, boolean initValues,
			boolean setIfNotPresent_, boolean amino){
		assert(!amino || (k_>0 && k_<=12)) : "Five-bit amino keys must fit in a nonnegative long for shard routing; k="+k_;
		tables=tables_;
		buflen=buflen_;
		halflen=(int)Math.ceil(buflen*0.5);
		ways=tables.length;
		buffers=new KmerBuffer[ways];
		setIfNotPresent=setIfNotPresent_;
		useValues=initValues;
		final int bitsPerBase=amino ? 5 : 2;
		final int symbolMask=amino ? 31 : 3;
		coreMask=(AbstractKmerTableSet.MASK_CORE ? ~(((-1L)<<(bitsPerBase*(k_-1)))|symbolMask) : -1L);
		middleMask=(AbstractKmerTableSet.MASK_MIDDLE ? makeMiddleMask(k_, amino) : -1L);
		cmMask=coreMask&middleMask;
		for(int i=0; i<ways; i++){
			buffers[i]=new KmerBuffer(buflen, k_, useValues);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Routes a nonnegative packed key using the masks captured at construction. */
	public final int kmerToWay(final long kmer){
		final int way=(int)((kmer&cmMask)%ways);
		return way;
	}

	/**
	 * Queues one count increment in count-only mode and may flush that shard's batch.
	 * Half-full buffers attempt a nonblocking flush at 16-entry intervals; a full
	 * buffer waits for its table lock. The return concerns the whole flushed batch,
	 * including earlier calls, rather than just this key.
	 *
	 * @param kmer The kmer to increment
	 * @param incr Increment value (must be 1)
	 * @return Sum of backend creation returns from this flush, or zero if entries remain buffered
	 */
	@Override
	public int incrementAndReturnNumCreated(final long kmer, final int incr){
		assert(incr==1); //I could just add the kmer multiple times if not true, with addMulti
		final int way=kmerToWay(kmer);
		KmerBuffer buffer=buffers[way];
//		final int size=buffer.addMulti(kmer, incr);
		final int size=buffer.add(kmer);
		//Comprehension (clever, verified): throttled dump — once the buffer is at least half full, attempt a dump every
		//16 adds ((size&SIZEMASK)==0, SIZEMASK=15) with a NON-blocking tryLock (dumpBuffer force=false), but force a
		//blocking dump only when completely full (size>=buflen). This keeps lock contention low while guaranteeing the
		//buffer never exceeds buflen (the force branch drains it at the boundary). incr is asserted ==1 here; bulk incr>1
		//is coalesced in dumpBuffer_inner into incrementAndReturnNumCreated(prev,sum), never into increment().
		if(size>=halflen && (size>=buflen || (size&SIZEMASK)==0)){
			return dumpBuffer(way, size>=buflen);
		}
		return 0;
	}

	/**
	 * Flushes all buffers to their underlying tables.
	 * Forces complete dump of all buffered kmers regardless of buffer fill level.
	 * @return Sum of backend update returns from all batches flushed by this call
	 */
	@Override
	public final long flush(){
		long added=0;
		for(int i=0; i<ways; i++){added+=dumpBuffer(i, true);}
		return added;
	}

	/**
	 * Queues a key/value pair in value mode and may flush that shard's entire batch.
	 * Requires initValues=true at construction. The constructor's setIfNotPresent
	 * flag selects conditional insertion during flush; otherwise backend set is used.
	 *
	 * @param kmer The kmer to set
	 * @param value The value to associate with the kmer
	 * @return Sum of backend update returns from the flushed batch, or zero if still buffered
	 */
	@Override
	public int set(long kmer, int value){
		final int way=kmerToWay(kmer);
		KmerBuffer buffer=buffers[way];
		final int size=buffer.add(kmer, value);
		if(size>=halflen && (size>=buflen || (size&SIZEMASK)==0)){
			return dumpBuffer(way, size>=buflen);
		}
		return 0;
	}

	/**
	 * Sets kmer to multiple values.
	 * Not implemented - this class lacks multi-value buffers.
	 *
	 * @param kmer The kmer to set
	 * @param vals Array of values
	 * @param vlen Length of values to use
	 * @return Never returns - throws RuntimeException
	 * @throws RuntimeException Always thrown as method is unimplemented
	 */
	@Override
	public int set(long kmer, int[] vals, int vlen){
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}

	/**
	 * Unsupported direct conditional-set entry point. For buffered conditional
	 * insertion, construct in value mode with setIfNotPresent=true and call set.
	 *
	 * @param kmer The kmer to conditionally set
	 * @param value The value to set if kmer not present
	 * @return Never returns - throws RuntimeException
	 * @throws RuntimeException Always thrown as method is unimplemented
	 */
	@Override
	public int setIfNotPresent(long kmer, int value){
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}

	/**
	 * Gets the value associated with a kmer.
	 * Delegates to the appropriate underlying table based on kmer hash.
	 * @param kmer The kmer to look up
	 * @return Backend lookup result; pending buffered writes are not included
	 */
	@Override
	public int getValue(long kmer){
		final int way=kmerToWay(kmer);
		return tables[way].getValue(kmer);
	}

	/**
	 * Gets all values associated with a kmer.
	 * Delegates to the appropriate underlying table based on kmer hash.
	 *
	 * @param kmer The kmer to look up
	 * @param singleton Reusable single-element array for efficiency
	 * @return Backend result, potentially borrowed storage or null; pending writes are excluded
	 */
	@Override
	public int[] getValues(long kmer, int[] singleton){
		final int way=kmerToWay(kmer);
		return tables[way].getValues(kmer, singleton);
	}

	/**
	 * Tests whether a kmer is present in the table.
	 * Delegates to the appropriate underlying table based on kmer hash.
	 * @param kmer The kmer to test for presence
	 * @return Whether the backend contains the key; pending writes are excluded
	 */
	@Override
	public boolean contains(long kmer){
		final int way=kmerToWay(kmer);
		return tables[way].contains(kmer);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/

	/** Initializes backend ownership without flushing or acquiring table locks.
	 * Callers finish loading and satisfy each backend's ownership preconditions first. */
	@Override
	public final void initializeOwnership(){
		for(AbstractKmerTable t : tables){t.initializeOwnership();}
	}

	/** Clears backend ownership without flushing or locking; backend post-clear access rules apply. */
	@Override
	public final void clearOwnership(){
		for(AbstractKmerTable t : tables){t.clearOwnership();}
	}

	/** Delegates an ownership claim for an already stored key; returns the backend result. */
	@Override
	public final int setOwner(final long kmer, final int newOwner){
		final int way=kmerToWay(kmer);
		return tables[way].setOwner(kmer, newOwner);
	}

	/**
	 * Clears ownership of a kmer if owned by specified owner.
	 * Delegates to the appropriate underlying table based on kmer hash.
	 *
	 * @param kmer The kmer to clear ownership for
	 * @param owner The expected current owner
	 * @return true if ownership was cleared, false if not owned by specified owner
	 */
	@Override
	public final boolean clearOwner(final long kmer, final int owner){
		final int way=kmerToWay(kmer);
		return tables[way].clearOwner(kmer, owner);
	}

	/**
	 * Gets the current owner of a kmer.
	 * Delegates to the appropriate underlying table based on kmer hash.
	 * @param kmer The kmer to get owner for
	 * @return Backend ownership result; the caller must satisfy its initialization/key preconditions
	 */
	@Override
	public final int getOwner(final long kmer){
		final int way=kmerToWay(kmer);
		return tables[way].getOwner(kmer);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Gets the raw object associated with a kmer.
	 * Delegates without flushing or locking; backends may reject object-valued lookup.
	 * @param kmer The kmer to look up
	 * @return The raw object associated with the kmer
	 */
	@Override
	Object get(long kmer){
		final int way=kmerToWay(kmer);
		return tables[way].get(kmer);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Tries or waits for one table's lock, then flushes pending entries and always unlocks.
	 * Even an empty buffer may acquire the lock. A failed tryLock leaves its batch intact. */
	private int dumpBuffer(final int way, boolean force){
		final KmerBuffer buffer=buffers[way];
		final AbstractKmerTable table=tables[way];
		final int lim=buffer.size();
		if(lim<0){return 0;}

		if(force){table.lock();}
		else if(!table.tryLock()){return 0;}

		try{
			if(SORT_BUFFERS && buffer.values==null){//Can go before or after lock; neither helps much
				buffer.kmers.sortSerial();
			}
			return dumpBuffer_inner(way);
		}finally{
			//KMR102: preserve the original failure without retaining the shared backend lock.
			table.unlock();
		}
	}

	/** Applies entries while holding the backend lock; clears only after normal completion.
	 * Partial failures are not transactional and buffered entries must not be blindly retried. */
	private int dumpBuffer_inner(final int way){
		if(verbose){System.err.println("Dumping buffer for way "+way+" of "+ways);}
		final KmerBuffer buffer=buffers[way];
		final int lim=buffer.size();
		if(lim<1){return 0;}
		final long[] kmers=buffer.kmers.array;
		final int[] values=(buffer.values==null ? null : buffer.values.array);
		if(lim<1){return 0;}
		int added=0;
		final AbstractKmerTable table=tables[way];
//		synchronized(table){
			if(values==null){
//				Arrays.sort(kmers, 0, lim); //Makes it slower
				if(SORT_BUFFERS){
					long prev=-1;
					int sum=0;
					for(int i=0; i<lim; i++){
						final long kmer=kmers[i];
						if(kmer==prev){
							sum++;
						}else{
							if(sum>0){added+=table.incrementAndReturnNumCreated(prev, sum);}
							prev=kmer;
							sum=1;
						}
					}
					if(sum>0){added+=table.incrementAndReturnNumCreated(prev, sum);}
				}else{
					for(int i=0; i<lim; i++){
						final long kmer=kmers[i];
						added+=table.incrementAndReturnNumCreated(kmer, 1);
					}
				}
			}else{
				if(setIfNotPresent){
					for(int i=0; i<lim; i++){
						final long kmer=kmers[i];
						final int value=values[i];
						added+=table.setIfNotPresent(kmer, value);
					}
				}else{
					for(int i=0; i<lim; i++){
						final long kmer=kmers[i];
						final int value=values[i];
						added+=table.set(kmer, value);
//						System.err.println("B: "+kmer+", "+Arrays.toString(((HashArrayHybrid)table).getValues(kmer, new int[1])));
					}
				}
			}
//		}
		buffer.clear();
		uniqueAdded+=added;
		return added;
	}

	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Indicates whether this table supports resizing.
	 * HashBuffer does not support resizing operations.
	 * @return Always false
	 */
	@Override
	final boolean canResize(){return false;}

	/**
	 * Indicates whether this table supports rebalancing.
	 * HashBuffer does not support rebalancing operations.
	 * @return Always false
	 */
	@Override
	public final boolean canRebalance(){return false;}

	/**
	 * Gets the total number of kmers in the table.
	 * Not implemented for HashBuffer as it's a buffering layer.
	 *
	 * @return Never returns - throws RuntimeException
	 * @throws RuntimeException Always thrown as method is unimplemented
	 * @deprecated Not supported for buffering tables
	 */
	@Deprecated
	@Override
	public long size(){
		throw new RuntimeException("Unimplemented.");
	}

	/**
	 * Gets the length of internal arrays.
	 * Not implemented for HashBuffer as it's a buffering layer.
	 *
	 * @return Never returns - throws RuntimeException
	 * @throws RuntimeException Always thrown as method is unimplemented
	 * @deprecated Not supported for buffering tables
	 */
	@Deprecated
	@Override
	public int arrayLength(){
		throw new RuntimeException("Unimplemented.");
	}

	/**
	 * Resizes the internal data structures.
	 * Not implemented for HashBuffer as it doesn't support resizing.
	 * @throws RuntimeException Always thrown as method is unimplemented
	 * @deprecated Not supported for buffering tables
	 */
	@Deprecated
	@Override
	void resize(){
		throw new RuntimeException("Unimplemented.");
	}

	/**
	 * Rebalances the table for better performance.
	 * Not implemented for HashBuffer as it doesn't support rebalancing.
	 * @throws RuntimeException Always thrown as method is unimplemented
	 * @deprecated Not supported for buffering tables
	 */
	@Deprecated
	@Override
	public void rebalance(){
		throw new RuntimeException("Unimplemented.");
	}

	/**
	 * Delegates regeneration without flushing or locking and sums backend results.
	 * The caller excludes conflicting access and respects each backend's support and ownership rules.
	 * @param limit Backend count threshold for removal, not a limit on entries processed
	 * @return Sum of the backends' removal counts; pending buffers are not flushed
	 */
	@Override
	public long regenerate(final int limit){
		long sum=0;
		for(AbstractKmerTable table : tables){
			sum+=table.regenerate(limit);
		}
		return sum;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Dumps kmers as text format to the specified writer.
	 * Delegates to all underlying tables without flushing pending keys or locking.
	 * The caller owns writer completion and error checks; backend return values are ignored.
	 *
	 * @param tsw Text stream writer for output
	 * @param k Kmer length for output formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @return true after all delegations return; does not certify successful output completion
	 */
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
		for(AbstractKmerTable table : tables){
			table.dumpKmersAsText(tsw, k, mincount, maxcount);
		}
		return true;
	}

	/**
	 * Formats kmers as text bytes for the specified writer.
	 * Delegates without flushing pending keys or locking. Count bounds and record
	 * allowance handling belong to each backend. Caller owns writer completion/error checks.
	 *
	 * @param bsw Byte stream writer for output
	 * @param k Kmer length for output formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @param remaining Shared record allowance, or null for unlimited output
	 * @return true after all delegations return; backend return values are ignored
	 */
	@Override
	public boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		for(AbstractKmerTable table : tables){
			table.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);
		}
		return true;
	}

	/**
	 * Unsupported caller-buffered text-byte dumping entry point.
	 *
	 * @param bsw Byte stream writer for output
	 * @param bb Byte builder for formatting
	 * @param k Kmer length for output formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @param remaining Atomic counter for remaining kmers to process
	 * @return Never returns - throws RuntimeException
	 * @throws RuntimeException Always thrown as method is unsupported
	 * @deprecated Not supported for buffering tables
	 */
	@Override
	@Deprecated
	public boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining){
		throw new RuntimeException("Unsupported.");
	}

	/**
	 * Fills histogram array with kmer count distributions.
	 * This wrapper does not implement histogram accumulation.
	 *
	 * @param ca Count array to fill with histogram data
	 * @param max Maximum count value to include
	 * @throws RuntimeException Always thrown as method is unsupported
	 * @deprecated Not supported for buffering tables
	 */
	@Override
	@Deprecated
	public void fillHistogram(long[] ca, int max){
		throw new RuntimeException("Unsupported.");
	}

	/**
	 * Fills histogram list with kmer count distributions.
	 * This wrapper does not implement histogram accumulation.
	 *
	 * @param sll SuperLongList to fill with histogram data
	 * @throws RuntimeException Always thrown as method is unsupported
	 * @deprecated Not supported for buffering tables
	 */
	@Override
	@Deprecated
	public void fillHistogram(SuperLongList sll){
		throw new RuntimeException("Unsupported.");
	}

	/**
	 * Counts GC content distribution across kmers.
	 * Delegates GC counting without flushing pending entries or locking tables.
	 * @param gcCounts Array to accumulate GC count statistics
	 * @param max Backend histogram cap; higher values may be folded into this bin
	 */
	@Override
	public void countGC(long[] gcCounts, int max){
		for(AbstractKmerTable table : tables){
			table.countGC(gcCounts, max);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Increments kmer count by specified amount.
	 * Not supported - use incrementAndReturnNumCreated instead.
	 *
	 * @param kmer The kmer to increment
	 * @param incr The increment amount
	 * @return Never returns - throws RuntimeException
	 * @throws RuntimeException Always thrown as method is unsupported
	 */
	@Override
	public int increment(final long kmer, final int incr){
		throw new RuntimeException("Unsupported");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed shared destination shards; callers keep their references stable. */
	private final AbstractKmerTable[] tables;
	/** Key count that forces a blocking shard flush. */
	private final int buflen;
	/** Ceiling of half the flush threshold, enabling periodic nonblocking flush attempts. */
	private final int halflen;
	/** Number of destination shards and private buffers. */
	private final int ways;
	/** Selects paired value storage during construction. */
	private final boolean useValues;
	/** Thread-confined reusable buffers, one per destination shard. */
	private final KmerBuffer[] buffers;
	/** End-base routing mask captured from MASK_CORE and the constructor width. */
	private final long coreMask;
	/** Middle-symbol routing mask for the selected alphabet, captured at construction. */
	private final long middleMask;
	/** Combined routing mask; does not modify the stored key. */
	private final long cmMask;
	/** Accumulated backend update returns from successfully completed batches. */
	public long uniqueAdded=0;

	/** Tests every sixteenth buffered entry for an opportunistic flush. */
	private static final int SIZEMASK=15;
	/** Uses backend conditional insertion when flushing value-mode entries. */
	private final boolean setIfNotPresent;

	/** Sorts count-only batches and coalesces repeated keys; value batches are never sorted here. */
	public static boolean SORT_BUFFERS=false;

}
