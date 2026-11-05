package ukmer;

import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * Thread-confined buffer that batches copied k-mer words into shared backend tables.
 * Each worker needs its own wrapper and reconstruction Kmer. Flushes lock one
 * backend at a time; lookups, ownership operations and aggregate operations do
 * not flush or lock. Pending keys are invisible to those operations.
 * Public insertion supports count-only buffers (initValues=false); the value
 * flush branch has no supported public method for supplying key/value pairs.
 * Call flush after normal loading. Partial failures are not transactional.
 *
 * @author Brian Bushnell
 * @date Nov 22, 2013
 */
public class HashBufferU extends AbstractKmerTableU{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates ordinary sequence-key buffers; see the five-argument constructor for units and modes. */
	public HashBufferU(AbstractKmerTableU[] tables_, int buflen_, int kbig_, boolean initValues){
		this(tables_, buflen_, kbig_, initValues, false);
	}

	/** Borrows backend tables and allocates private buffers and a reusable reconstruction key.
	 * Packing and hash settings must remain consistent with incoming keys and backends.
	 * @param tables_ Nonempty array of compatible backend tables
	 * @param buflen_ Positive threshold in logical keys; buflen_+1 and word-count products must fit in an int
	 * @param kbig_ Requested total length, normalized by Kmer's current packing geometry
	 * @param initValues Allocates value storage; must be false for supported public count insertion
	 * @param hollow Reconstructs hollow keys; requires two-hash HashArrayH1D backends
	 * @throws IllegalArgumentException If hollow mode has an incompatible backend */
	public HashBufferU(AbstractKmerTableU[] tables_, int buflen_, int kbig_, boolean initValues, boolean hollow){
		tables=tables_;
		buflen=buflen_;
		kmer=new Kmer(kbig_, hollow);
		if(hollow){
			for(AbstractKmerTableU table : tables){
				if(!(table instanceof HashArrayH1D) || ((HashArrayH1D)table).fingerprintOnly()){
					throw new IllegalArgumentException("Hollow buffers require two-hash tables.");
				}
			}
		}
		mult=kmer.mult;
		buflen2=buflen*mult;
		halflen2=((buflen+1)/2)*mult;
		ways=tables.length;
		buffers=new KmerBufferU[ways];
		for(int i=0; i<ways; i++){
			buffers[i]=new KmerBufferU(buflen, kmer.kbig, initValues);
		}
//		tempKey=new long[mult];
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

//	@Override
//	public int incrementAndReturnNumCreated(Kmer kmer){
//		assert(kmer.mult==mult) : kmer.mult+"!="+mult+", kbig="+kmer.kbig+", k="+kmer.k;
//		final int way=getWay(kmer);
//		KmerBufferU buffer=buffers[way];
//		final int size=buffer.add(kmer);
//		if(size==halflen2 || size>=buflen2){
//			return dumpBuffer(way, size>=buflen2);
//		}
//		return 0;
//	}

	/**
	 * Copies one key into its shard's count-only buffer and may flush that whole batch.
	 * The key must match the construction geometry. At or above half capacity,
	 * multiples of 16 buffered long words trigger a try-lock; a full buffer waits
	 * for the lock. Thresholds count words, not calls, when mult is greater than one.
	 *
	 * @param kmer The k-mer to increment
	 * @return Sum of backend creation returns for the whole flushed batch, or zero if still buffered
	 */
	@Override
	public int incrementAndReturnNumCreated(Kmer kmer){
		assert(kmer.mult==mult) : kmer.mult+"!="+mult+", kbig="+kmer.kbig+", k="+kmer.k;
		assert(kmer.hollow==this.kmer.hollow) : "Buffered keys and their reconstruction must use identical hash geometry.";
		final int way=getWay(kmer);
		KmerBufferU buffer=buffers[way];
		final int size=buffer.add(kmer);
		if(size>=halflen2 && (size>=buflen2 || (size&SIZEMASK)==0)){
//		if(size==halflen2 || size>=buflen2){
			return dumpBuffer(way, size>=buflen2);
		}
		return 0;
	}

	/**
	 * Flushes all buffers to their respective backend tables.
	 * Forces all buffered k-mers to be written to the underlying hash tables.
	 * @return Sum of backend update returns for batches flushed by this call
	 */
	@Override
	public final long flush(){
		long added=0;
		for(int i=0; i<ways; i++){added+=dumpBuffer(i, true);}
		return added;
	}

	/**
	 * Sets a k-mer value (unsupported operation).
	 * No public method populates the optional value buffers.
	 *
	 * @param kmer The k-mer to set
	 * @param value The value to assign
	 * @return Never returns, always throws exception
	 * @throws RuntimeException Always thrown as operation is unsupported
	 */
	@Override
	public int set(Kmer kmer, int value){
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}

	/**
	 * Sets a k-mer with multiple values (unsupported operation).
	 * No public method populates the optional value buffers.
	 *
	 * @param kmer The k-mer to set
	 * @param vals Array of values to assign
	 * @return Never returns, always throws exception
	 * @throws RuntimeException Always thrown as operation is unsupported
	 */
	@Override
	public int set(Kmer kmer, int[] vals){
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}

	/**
	 * Sets k-mer value only if not already present (unsupported operation).
	 * No public method populates the optional value buffers.
	 *
	 * @param kmer The k-mer to conditionally set
	 * @param value The value to assign if k-mer not present
	 * @return Never returns, always throws exception
	 * @throws RuntimeException Always thrown as operation is unsupported
	 */
	@Override
	public int setIfNotPresent(Kmer kmer, int value){
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}

	/**
	 * Retrieves the value associated with a k-mer from the appropriate backend table.
	 * Routes k-mer to correct table using hash-based way selection.
	 * @param kmer The k-mer to look up
	 * @return Backend lookup result, excluding pending buffered writes
	 */
	@Override
	public int getValue(Kmer kmer){
		final int way=getWay(kmer);
		return tables[way].getValue(kmer);
	}

	/**
	 * Retrieves value using raw k-mer key and XOR hash.
	 * Routes to appropriate backend table using XOR modulo way count.
	 *
	 * @param key Key representation expected by the backend
	 * @param xor Nonnegative routing hash consistent with the key and stored table geometry
	 * @return Backend lookup result, excluding pending buffered writes
	 */
	@Override
	public int getValue(long[] key, long xor){
		final int way=(int)(xor%ways);
		return tables[way].getValue(key, xor);
	}

	/**
	 * Retrieves multiple values associated with a k-mer.
	 * Delegates to appropriate backend table for multi-value retrieval.
	 *
	 * @param kmer The k-mer to look up
	 * @param singleton Pre-allocated array for single-value results
	 * @return Backend result, possibly borrowed storage or null; pending keys are excluded
	 */
	@Override
	public int[] getValues(Kmer kmer, int[] singleton){
		final int way=getWay(kmer);
		return tables[way].getValues(kmer, singleton);
	}

	/**
	 * Checks if a k-mer is present in any of the backend tables.
	 * Routes k-mer to appropriate table for containment check.
	 * @param kmer The k-mer to check
	 * @return Whether the selected backend contains the key; pending keys are excluded
	 */
	@Override
	public boolean contains(Kmer kmer){
		final int way=getWay(kmer);
		return tables[way].contains(kmer);
	}

	/** Selects the backend by the key's nonnegative hash, updating its cached hash if needed. */
	public final int getWay(Kmer kmer){return (int)(kmer.xor()%ways);}

	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/

	/** Initializes backend ownership without flushing or locking; callers finish loading first. */
	@Override
	public final void initializeOwnership(){
		for(AbstractKmerTableU t : tables){t.initializeOwnership();}
	}

	/** Clears backend ownership without flushing or locking; backend post-clear access rules apply. */
	@Override
	public final void clearOwnership(){
		for(AbstractKmerTableU t : tables){t.clearOwnership();}
	}

	/**
	 * Sets the owner of a k-mer in the appropriate backend table.
	 * @param kmer The k-mer to assign ownership
	 * @param newOwner The new owner ID
	 * @return Backend ownership result, normally the resulting owner rather than the previous owner
	 */
	@Override
	public final int setOwner(final Kmer kmer, final int newOwner){
		final int way=getWay(kmer);
		return tables[way].setOwner(kmer, newOwner);
	}

	/**
	 * Clears ownership of a k-mer if owned by specified owner.
	 * @param kmer The k-mer to clear ownership
	 * @param owner The expected current owner
	 * @return true if ownership was cleared, false if not owned by specified owner
	 */
	@Override
	public final boolean clearOwner(final Kmer kmer, final int owner){
		final int way=getWay(kmer);
		return tables[way].clearOwner(kmer, owner);
	}

	/**
	 * Gets the current owner of a k-mer.
	 * @param kmer The k-mer to check ownership
	 * @return Backend ownership result; its initialization and key-presence preconditions apply
	 */
	@Override
	public final int getOwner(final Kmer kmer){
		final int way=getWay(kmer);
		return tables[way].getOwner(kmer);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Retrieves the object associated with a k-mer from appropriate backend table.
	 * Does not flush or lock; a backend may reject object-valued lookup.
	 * @param kmer The k-mer to look up
	 * @return Object associated with the k-mer
	 */
	@Override
	Object get(Kmer kmer){
		final int way=getWay(kmer);
		return tables[way].get(kmer);
	}

	/**
	 * Retrieves object using raw k-mer array (unsupported operation).
	 * @param kmer Raw k-mer array
	 * @return Never returns, always throws exception
	 * @throws RuntimeException Always thrown as operation is unsupported
	 */
	@Override
	Object get(long[] kmer){
		throw new RuntimeException();
	}

	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Acquires or tries one backend lock, flushes, and releases it even on failure.
	 * An unsuccessful try-lock preserves the batch; an empty batch may still acquire a lock. */
	private int dumpBuffer(final int way, boolean force){
		final KmerBufferU buffer=buffers[way];
		final AbstractKmerTableU table=tables[way];
		final int lim=buffer.size();
		if(lim<0){return 0;}
		if(force){table.lock();}
		else if(!table.tryLock()){return 0;}
		try{return dumpBuffer_inner(way);}
		finally{table.unlock();}
	}

	/** Reconstructs buffered keys under the backend lock and clears only on normal completion.
	 * Earlier updates survive a later exception; blindly retrying the retained batch can double-count. */
	private int dumpBuffer_inner(final int way){
		if(verbose){System.err.println("Dumping buffer for way "+way+" of "+ways);}
		final KmerBufferU buffer=buffers[way];
		final int lim=buffer.size();
		if(lim<1){return 0;}
		final long[] kmers=buffer.kmers.array;
		final int[] values=(buffer.values==null ? null : buffer.values.array);
		int added=0;
		final AbstractKmerTableU table=tables[way];
		final long array1[]=kmer.array1();
//		synchronized(table){
			if(values==null){
//				System.err.println("way="+way);
				for(int j=0; j<lim;){
					for(int x=0; x<mult; x++, j++){
						if(verbose){System.err.println("x="+x+", j="+j);}
						array1[x]=kmers[j];
					}
					kmer.fillArray2();
					if(verbose){System.err.println("Incrementing "+kmer+"; xor="+kmer.xor());}
//					assert(kmer.mod(ways)==way) : kmer+", "+way+", "+ways+", "+kmer.mod(ways)+", "+kmer.xor()+"\n"+
//						Arrays.toString(kmer.array1())+"\n"+Arrays.toString(kmer.array2())+"\n"+Arrays.toString(kmer.key());
//					assert(kmer.verify(false));
					int x=table.incrementAndReturnNumCreated(kmer);
					added+=x;
				}
			}else{
				for(int i=0, j=0; j<lim; i++){
					for(int x=0; x<mult; x++, j++){
						array1[x]=kmers[j];
					}
					kmer.fillArray2();
					added+=table.setIfNotPresent(kmer, values[i]);
				}
			}
//		}
		buffer.clear();
		return added;
	}

	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns false as this implementation does not support resizing */
	@Override
	final boolean canResize(){return false;}

	/** Returns false as this implementation does not support rebalancing */
	@Override
	public final boolean canRebalance(){return false;}

	/**
	 * Gets total size (deprecated/unsupported operation).
	 * @return Never returns, always throws exception
	 * @throws RuntimeException Always thrown as operation is unimplemented
	 * @deprecated This operation is not supported
	 */
	@Deprecated
	@Override
	public long size(){
		throw new RuntimeException("Unimplemented.");
	}

	/**
	 * Gets array length (deprecated/unsupported operation).
	 * @return Never returns, always throws exception
	 * @throws RuntimeException Always thrown as operation is unimplemented
	 * @deprecated This operation is not supported
	 */
	@Deprecated
	@Override
	public int arrayLength(){
		throw new RuntimeException("Unimplemented.");
	}

	/**
	 * Resizes the hash table (deprecated/unsupported operation).
	 * @throws RuntimeException Always thrown as operation is unimplemented
	 * @deprecated This operation is not supported
	 */
	@Deprecated
	@Override
	void resize(){
		throw new RuntimeException("Unimplemented.");
	}

	/**
	 * Rebalances the hash table (deprecated/unsupported operation).
	 * @throws RuntimeException Always thrown as operation is unimplemented
	 * @deprecated This operation is not supported
	 */
	@Deprecated
	@Override
	public void rebalance(){
		throw new RuntimeException("Unimplemented.");
	}

	/**
	 * Delegates regeneration without flushing or locking and sums backend results.
	 * Callers exclude conflicting access and respect backend support and ownership rules.
	 * @param limit Backend count threshold for removal, not a limit on entries processed
	 * @return Sum of backend removal counts, excluding pending buffered entries
	 */
	@Override
	public long regenerate(final int limit){
		long sum=0;
		for(AbstractKmerTableU table : tables){
			sum+=table.regenerate(limit);
		}
		return sum;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Dumps k-mers from all backend tables as text format.
	 * Delegates without flushing pending keys or locking; backend return values are ignored.
	 * The caller owns writer completion and error checks.
	 *
	 * @param tsw Text stream writer for output
	 * @param k K-mer length for output formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @return true after all delegations return; does not certify output completion
	 */
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
		for(AbstractKmerTableU table : tables){
			table.dumpKmersAsText(tsw, k, mincount, maxcount);
		}
		return true;
	}

	/**
	 * Formats backend k-mers as text bytes without flushing pending keys or locking.
	 * Count filtering and record allowances follow each backend's contract.
	 * The caller owns writer completion and error checks.
	 *
	 * @param bsw Byte stream writer for text output
	 * @param k K-mer length for output formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @param remaining Shared record allowance, or null for unlimited output
	 * @return true after all delegations return; backend return values are ignored
	 */
	@Override
	public boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		for(AbstractKmerTableU table : tables){
			table.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);
		}
		return true;
	}

	/**
	 * Caller-buffered text-byte dumping (deprecated/unsupported operation).
	 *
	 * @param bsw Byte stream writer
	 * @param bb Byte builder for output formatting
	 * @param k K-mer length
	 * @param mincount Minimum count threshold
	 * @param maxcount Maximum count threshold
	 * @param remaining Atomic counter for remaining items
	 * @return Never returns, always throws exception
	 * @throws RuntimeException Always thrown as operation is unsupported
	 * @deprecated This operation is not supported
	 */
	@Override
	@Deprecated
	public boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining){
		throw new RuntimeException("Unsupported.");
	}

	/**
	 * Fills histogram array (deprecated/unsupported operation).
	 *
	 * @param ca Count array for histogram
	 * @param max Maximum count value
	 * @throws RuntimeException Always thrown as operation is unsupported
	 * @deprecated This operation is not supported
	 */
	@Override
	@Deprecated
	public void fillHistogram(long[] ca, int max){
		throw new RuntimeException("Unsupported.");
	}

	/**
	 * Fills histogram using SuperLongList (deprecated/unsupported operation).
	 * @param sll SuperLongList for histogram data
	 * @throws RuntimeException Always thrown as operation is unsupported
	 * @deprecated This operation is not supported
	 */
	@Override
	@Deprecated
	public void fillHistogram(SuperLongList sll){
		throw new RuntimeException("Unsupported.");
	}

	/**
	 * Counts GC content distribution across all backend tables.
	 * Delegates without flushing pending entries or locking tables.
	 * @param gcCounts Array to store GC count distribution
	 * @param max Backend histogram cap; higher counts may be folded into this bin
	 */
	@Override
	public void countGC(long[] gcCounts, int max){
		for(AbstractKmerTableU table : tables){
			table.countGC(gcCounts, max);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Increments k-mer count (unsupported operation).
	 * Use incrementAndReturnNumCreated instead.
	 *
	 * @param kmer The k-mer to increment
	 * @return Never returns, always throws exception
	 * @throws RuntimeException Always thrown as operation is unsupported
	 */
	@Override
	public int increment(Kmer kmer){
		throw new RuntimeException("Unsupported");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed shared destination shards; callers keep their references stable. */
	private final AbstractKmerTableU[] tables;
	/** Blocking flush threshold in logical keys. */
	private final int buflen;
	/** Blocking flush threshold in stored long words. */
	private final int buflen2;
	/** Ceiling of half the key threshold, expressed in stored long words. */
	private final int halflen2;
	/** Stored long words per logical key. */
	private final int mult;
	/** Number of backend shards and private buffers. */
	private final int ways;
	/** Owned reusable word buffers, one per shard. */
	private final KmerBufferU[] buffers;
	/** Thread-confined scratch key used to reconstruct words before backend updates. */
	private final Kmer kmer;

	/** Enables opportunistic flush attempts at multiples of 16 occupied long words. */
	private static final int SIZEMASK=15;

}
