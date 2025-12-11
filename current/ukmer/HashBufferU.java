package ukmer;

import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * Multi-way buffered k-mer hash table for efficient k-mer tracking and management.
 * Distributes k-mers across multiple backend hash tables using hash-based routing.
 * Provides smart buffer flushing with force and try-lock mechanisms for optimal performance.
 *
 * @author Brian Bushnell
 * @date Nov 22, 2013
 */
public class HashBufferU extends AbstractKmerTableU {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public HashBufferU(AbstractKmerTableU[] tables_, int buflen_, int kbig_, boolean initValues){
		tables=tables_;
		buflen=buflen_;
		kmer=new Kmer(kbig_);
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
//	public int incrementAndReturnNumCreated(Kmer kmer) {
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
	 * Increments the count for a k-mer and returns the number of new k-mers created.
	 * Routes k-mer to appropriate buffer based on hash value and triggers buffer dumps
	 * when thresholds are reached.
	 *
	 * @param kmer The k-mer to increment
	 * @return Number of new k-mers created during this operation
	 */
	@Override
	public int incrementAndReturnNumCreated(Kmer kmer) {
		assert(kmer.mult==mult) : kmer.mult+"!="+mult+", kbig="+kmer.kbig+", k="+kmer.k;
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
	 * @return Total number of k-mers added to tables during flush
	 */
	@Override
	public final long flush(){
		long added=0;
		for(int i=0; i<ways; i++){added+=dumpBuffer(i, true);}
		return added;
	}
	
	/**
	 * Sets a k-mer value (unsupported operation).
	 * This class lacks value buffers and cannot set arbitrary k-mer values.
	 *
	 * @param kmer The k-mer to set
	 * @param value The value to assign
	 * @return Never returns, always throws exception
	 * @throws RuntimeException Always thrown as operation is unsupported
	 */
	@Override
	public int set(Kmer kmer, int value) {
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}
	
	/**
	 * Sets a k-mer with multiple values (unsupported operation).
	 * This class lacks value buffers and cannot set arbitrary k-mer values.
	 *
	 * @param kmer The k-mer to set
	 * @param vals Array of values to assign
	 * @return Never returns, always throws exception
	 * @throws RuntimeException Always thrown as operation is unsupported
	 */
	@Override
	public int set(Kmer kmer, int[] vals) {
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}
	
	/**
	 * Sets k-mer value only if not already present (unsupported operation).
	 * This class lacks value buffers and cannot set arbitrary k-mer values.
	 *
	 * @param kmer The k-mer to conditionally set
	 * @param value The value to assign if k-mer not present
	 * @return Never returns, always throws exception
	 * @throws RuntimeException Always thrown as operation is unsupported
	 */
	@Override
	public int setIfNotPresent(Kmer kmer, int value) {
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}
	
	/**
	 * Retrieves the value associated with a k-mer from the appropriate backend table.
	 * Routes k-mer to correct table using hash-based way selection.
	 * @param kmer The k-mer to look up
	 * @return The value associated with the k-mer, or default if not found
	 */
	@Override
	public int getValue(Kmer kmer) {
		final int way=getWay(kmer);
		return tables[way].getValue(kmer);
	}
	
	/**
	 * Retrieves value using raw k-mer key and XOR hash.
	 * Routes to appropriate backend table using XOR modulo way count.
	 *
	 * @param key Raw k-mer key array
	 * @param xor XOR hash value for routing
	 * @return The value associated with the k-mer key
	 */
	@Override
	public int getValue(long[] key, long xor) {
		final int way=(int)(xor%ways);
		return tables[way].getValue(key, xor);
	}
	
	/**
	 * Retrieves multiple values associated with a k-mer.
	 * Delegates to appropriate backend table for multi-value retrieval.
	 *
	 * @param kmer The k-mer to look up
	 * @param singleton Pre-allocated array for single-value results
	 * @return Array of values associated with the k-mer
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
	 * @return true if k-mer is present, false otherwise
	 */
	@Override
	public boolean contains(Kmer kmer) {
		final int way=getWay(kmer);
		return tables[way].contains(kmer);
	}
	
	public final int getWay(Kmer kmer){return (int)(kmer.xor()%ways);}
	
	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Initializes ownership tracking for all backend tables */
	@Override
	public final void initializeOwnership(){
		for(AbstractKmerTableU t : tables){t.initializeOwnership();}
	}
	
	/** Clears ownership tracking for all backend tables */
	@Override
	public final void clearOwnership(){
		for(AbstractKmerTableU t : tables){t.clearOwnership();}
	}
	
	/**
	 * Sets the owner of a k-mer in the appropriate backend table.
	 * @param kmer The k-mer to assign ownership
	 * @param newOwner The new owner ID
	 * @return Previous owner ID
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
	 * @return Owner ID, or default value if unowned
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
	 * @param kmer The k-mer to look up
	 * @return Object associated with the k-mer
	 */
	@Override
	Object get(Kmer kmer) {
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
	Object get(long[] kmer) {
		throw new RuntimeException();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private int dumpBuffer(final int way, boolean force){
		final KmerBufferU buffer=buffers[way];
		final AbstractKmerTableU table=tables[way];
		final int lim=buffer.size();
		if(lim<0){return 0;}
		if(force){table.lock();}
		else if(!table.tryLock()){return 0;}
		final int x=dumpBuffer_inner(way);
		table.unlock();
		return x;
	}
	
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
	final boolean canResize() {return false;}
	
	/** Returns false as this implementation does not support rebalancing */
	@Override
	public final boolean canRebalance() {return false;}
	
	/**
	 * Gets total size (deprecated/unsupported operation).
	 * @return Never returns, always throws exception
	 * @throws RuntimeException Always thrown as operation is unimplemented
	 * @deprecated This operation is not supported
	 */
	@Deprecated
	@Override
	public long size() {
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
	public int arrayLength() {
		throw new RuntimeException("Unimplemented.");
	}
	
	/**
	 * Resizes the hash table (deprecated/unsupported operation).
	 * @throws RuntimeException Always thrown as operation is unimplemented
	 * @deprecated This operation is not supported
	 */
	@Deprecated
	@Override
	void resize() {
		throw new RuntimeException("Unimplemented.");
	}
	
	/**
	 * Rebalances the hash table (deprecated/unsupported operation).
	 * @throws RuntimeException Always thrown as operation is unimplemented
	 * @deprecated This operation is not supported
	 */
	@Deprecated
	@Override
	public void rebalance() {
		throw new RuntimeException("Unimplemented.");
	}
	
	/**
	 * Regenerates all backend tables with the specified limit.
	 * Delegates regeneration to each backend table and sums results.
	 * @param limit Maximum number of elements to regenerate
	 * @return Total number of elements regenerated across all tables
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
	 * Delegates to each backend table for text output generation.
	 *
	 * @param tsw Text stream writer for output
	 * @param k K-mer length for output formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @return true when dumping completes successfully
	 */
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
		for(AbstractKmerTableU table : tables){
			table.dumpKmersAsText(tsw, k, mincount, maxcount);
		}
		return true;
	}
	
	/**
	 * Dumps k-mers from all backend tables as binary format.
	 * Delegates to each backend table for binary output generation.
	 *
	 * @param bsw Byte stream writer for binary output
	 * @param k K-mer length for output formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @param remaining Atomic counter for remaining items to process
	 * @return true when dumping completes successfully
	 */
	@Override
	public boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		for(AbstractKmerTableU table : tables){
			table.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);
		}
		return true;
	}
	
	/**
	 * Multi-threaded k-mer dumping as bytes (deprecated/unsupported operation).
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
	 * Delegates GC counting to each backend table and aggregates results.
	 * @param gcCounts Array to store GC count distribution
	 * @param max Maximum count value to consider
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
	public int increment(Kmer kmer) {
		throw new RuntimeException("Unsupported");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private final AbstractKmerTableU[] tables;
	private final int buflen;
	private final int buflen2;
	private final int halflen2;
	private final int mult;
	private final int ways;
	private final KmerBufferU[] buffers;
	private final Kmer kmer;
	
	private static final int SIZEMASK=15;

}
