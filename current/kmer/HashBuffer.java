package kmer;

import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * Buffered kmer table that distributes kmers across multiple underlying tables.
 * Uses in-memory buffers to batch kmer operations before flushing to persistent storage.
 * Kmers are distributed to different tables based on hash values for load balancing.
 *
 * @author Brian Bushnell
 * @date Nov 22, 2013
 */
public class HashBuffer extends AbstractKmerTable {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public HashBuffer(AbstractKmerTable[] tables_, int buflen_, int k_, boolean initValues, boolean setIfNotPresent_){
		tables=tables_;
		buflen=buflen_;
		halflen=(int)Math.ceil(buflen*0.5);
		ways=tables.length;
		buffers=new KmerBuffer[ways];
		setIfNotPresent=setIfNotPresent_;
		useValues=initValues;
		coreMask=(AbstractKmerTableSet.MASK_CORE ? ~(((-1L)<<(2*(k_-1)))|3) : -1L);
		middleMask=(AbstractKmerTableSet.MASK_MIDDLE ? makeMiddleMask(k_, false) : -1L); //Note - this does not support amino acids.
		cmMask=coreMask&middleMask;
		for(int i=0; i<ways; i++){
			buffers[i]=new KmerBuffer(buflen, k_, useValues);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public final int kmerToWay(final long kmer){
		final int way=(int)((kmer&cmMask)%ways);
		return way;
	}
	
	/**
	 * Adds kmer to buffer and returns number of new kmers created.
	 * Buffers the kmer and flushes to underlying table if buffer fills up.
	 *
	 * @param kmer The kmer to increment
	 * @param incr Increment value (must be 1)
	 * @return Number of new kmers created after potential buffer flush
	 */
	@Override
	public int incrementAndReturnNumCreated(final long kmer, final int incr) {
		assert(incr==1); //I could just add the kmer multiple times if not true, with addMulti
		final int way=kmerToWay(kmer);
		KmerBuffer buffer=buffers[way];
//		final int size=buffer.addMulti(kmer, incr);
		final int size=buffer.add(kmer);
		if(size>=halflen && (size>=buflen || (size&SIZEMASK)==0)){
			return dumpBuffer(way, size>=buflen);
		}
		return 0;
	}
	
	/**
	 * Flushes all buffers to their underlying tables.
	 * Forces complete dump of all buffered kmers regardless of buffer fill level.
	 * @return Total number of new kmers added across all tables
	 */
	@Override
	public final long flush(){
		long added=0;
		for(int i=0; i<ways; i++){added+=dumpBuffer(i, true);}
		return added;
	}
	
	/**
	 * Sets kmer to specific value via buffering.
	 * Buffers the kmer-value pair and flushes if buffer fills up.
	 *
	 * @param kmer The kmer to set
	 * @param value The value to associate with the kmer
	 * @return Number of new kmers created after potential buffer flush
	 */
	@Override
	public int set(long kmer, int value) {
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
	public int set(long kmer, int[] vals, int vlen) {
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}
	
	/**
	 * Sets kmer value only if not already present.
	 * Not implemented - this class lacks value buffers for conditional setting.
	 *
	 * @param kmer The kmer to conditionally set
	 * @param value The value to set if kmer not present
	 * @return Never returns - throws RuntimeException
	 * @throws RuntimeException Always thrown as method is unimplemented
	 */
	@Override
	public int setIfNotPresent(long kmer, int value) {
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}
	
	/**
	 * Gets the value associated with a kmer.
	 * Delegates to the appropriate underlying table based on kmer hash.
	 * @param kmer The kmer to look up
	 * @return The value associated with the kmer, or default if not found
	 */
	@Override
	public int getValue(long kmer) {
		final int way=kmerToWay(kmer);
		return tables[way].getValue(kmer);
	}
	
	/**
	 * Gets all values associated with a kmer.
	 * Delegates to the appropriate underlying table based on kmer hash.
	 *
	 * @param kmer The kmer to look up
	 * @param singleton Reusable single-element array for efficiency
	 * @return Array of values associated with the kmer
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
	 * @return true if kmer is present, false otherwise
	 */
	@Override
	public boolean contains(long kmer) {
		final int way=kmerToWay(kmer);
		return tables[way].contains(kmer);
	}
	

	
	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final void initializeOwnership(){
		for(AbstractKmerTable t : tables){t.initializeOwnership();}
	}
	
	/** Clears ownership tracking for all underlying tables.
	 * Delegates ownership clearing to each underlying table. */
	@Override
	public final void clearOwnership(){
		for(AbstractKmerTable t : tables){t.clearOwnership();}
	}
	
	/** Initializes ownership tracking for all underlying tables.
	 * Delegates ownership initialization to each underlying table. */
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
	 * @return The owner identifier, or default if unowned
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
	 * Package-private method that delegates to underlying table.
	 * @param kmer The kmer to look up
	 * @return The raw object associated with the kmer
	 */
	@Override
	Object get(long kmer) {
		final int way=kmerToWay(kmer);
		return tables[way].get(kmer);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private int dumpBuffer(final int way, boolean force){
		final KmerBuffer buffer=buffers[way];
		final AbstractKmerTable table=tables[way];
		final int lim=buffer.size();
		if(lim<0){return 0;}
		
		if(force){table.lock();}
		else if(!table.tryLock()){return 0;}
		
		if(SORT_BUFFERS && buffer.values==null){//Can go before or after lock; neither helps much
			buffer.kmers.sortSerial();
		}
		
		final int x=dumpBuffer_inner(way);
		table.unlock();
		return x;
	}
	
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
	final boolean canResize() {return false;}
	
	/**
	 * Indicates whether this table supports rebalancing.
	 * HashBuffer does not support rebalancing operations.
	 * @return Always false
	 */
	@Override
	public final boolean canRebalance() {return false;}
	
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
	public long size() {
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
	public int arrayLength() {
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
	void resize() {
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
	public void rebalance() {
		throw new RuntimeException("Unimplemented.");
	}
	
	/**
	 * Regenerates kmers in all underlying tables up to specified limit.
	 * Delegates regeneration to each underlying table and sums results.
	 * @param limit Maximum number of kmers to regenerate per table
	 * @return Total number of kmers regenerated across all tables
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
	 * Delegates to all underlying tables to write their kmers.
	 *
	 * @param tsw Text stream writer for output
	 * @param k Kmer length for output formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @return Always true indicating successful completion
	 */
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
		for(AbstractKmerTable table : tables){
			table.dumpKmersAsText(tsw, k, mincount, maxcount);
		}
		return true;
	}
	
	/**
	 * Dumps kmers as binary format to the specified writer.
	 * Delegates to all underlying tables to write their kmers.
	 *
	 * @param bsw Byte stream writer for output
	 * @param k Kmer length for output formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @param remaining Atomic counter for remaining kmers to process
	 * @return Always true indicating successful completion
	 */
	@Override
	public boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		for(AbstractKmerTable table : tables){
			table.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);
		}
		return true;
	}
	
	/**
	 * Multi-threaded version of binary kmer dumping.
	 * Not supported for HashBuffer due to buffering complexity.
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
	 * Not supported for HashBuffer due to buffering complexity.
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
	 * Not supported for HashBuffer due to buffering complexity.
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
	 * Delegates GC counting to all underlying tables.
	 * @param gcCounts Array to accumulate GC count statistics
	 * @param max Maximum count value to include in statistics
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
	public int increment(final long kmer, final int incr) {
		throw new RuntimeException("Unsupported");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private final AbstractKmerTable[] tables;
	private final int buflen;
	private final int halflen;
	private final int ways;
	private final boolean useValues;
	private final KmerBuffer[] buffers;
	private final long coreMask;
	private final long middleMask;
	private final long cmMask;
	public long uniqueAdded=0;
	
	private static final int SIZEMASK=15;
	private final boolean setIfNotPresent;
	
	public static boolean SORT_BUFFERS=false;

}
