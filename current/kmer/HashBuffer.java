package kmer;

import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * @author Brian Bushnell
 * @date Nov 22, 2013
 *
 */
public class HashBuffer extends AbstractKmerTable {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs a HashBuffer with specified underlying tables and buffer configuration.
	 * Initializes individual buffers for each table and calculates masking parameters.
	 *
	 * @param tables_ Array of underlying kmer tables to distribute data to
	 * @param buflen_ Size of each internal buffer before flushing
	 * @param k_ Kmer length for mask calculations
	 * @param initValues Whether to initialize value buffers
	 * @param setIfNotPresent_ Whether to use setIfNotPresent mode
	 */
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
	
	@Override
	public final long flush(){
		long added=0;
		for(int i=0; i<ways; i++){added+=dumpBuffer(i, true);}
		return added;
	}
	
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
	
	@Override
	public int set(long kmer, int[] vals, int vlen) {
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}
	
	@Override
	public int setIfNotPresent(long kmer, int value) {
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}
	
	@Override
	public int getValue(long kmer) {
		final int way=kmerToWay(kmer);
		return tables[way].getValue(kmer);
	}
	
	@Override
	public int[] getValues(long kmer, int[] singleton){
		final int way=kmerToWay(kmer);
		return tables[way].getValues(kmer, singleton);
	}
	
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
	
	@Override
	public final void clearOwnership(){
		for(AbstractKmerTable t : tables){t.clearOwnership();}
	}
	
	@Override
	public final int setOwner(final long kmer, final int newOwner){
		final int way=kmerToWay(kmer);
		return tables[way].setOwner(kmer, newOwner);
	}
	
	@Override
	public final boolean clearOwner(final long kmer, final int owner){
		final int way=kmerToWay(kmer);
		return tables[way].clearOwner(kmer, owner);
	}
	
	@Override
	public final int getOwner(final long kmer){
		final int way=kmerToWay(kmer);
		return tables[way].getOwner(kmer);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	@Override
	final boolean canResize() {return false;}
	
	@Override
	public final boolean canRebalance() {return false;}
	
	@Deprecated
	@Override
	public long size() {
		throw new RuntimeException("Unimplemented.");
	}
	
	@Deprecated
	@Override
	public int arrayLength() {
		throw new RuntimeException("Unimplemented.");
	}
	
	@Deprecated
	@Override
	void resize() {
		throw new RuntimeException("Unimplemented.");
	}
	
	@Deprecated
	@Override
	public void rebalance() {
		throw new RuntimeException("Unimplemented.");
	}
	
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
	
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
		for(AbstractKmerTable table : tables){
			table.dumpKmersAsText(tsw, k, mincount, maxcount);
		}
		return true;
	}
	
	@Override
	public boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		for(AbstractKmerTable table : tables){
			table.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);
		}
		return true;
	}
	
	@Override
	@Deprecated
	public boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining){
		throw new RuntimeException("Unsupported.");
	}
	
	@Override
	@Deprecated
	public void fillHistogram(long[] ca, int max){
		throw new RuntimeException("Unsupported.");
	}
	
	@Override
	@Deprecated
	public void fillHistogram(SuperLongList sll){
		throw new RuntimeException("Unsupported.");
	}
	
	@Override
	public void countGC(long[] gcCounts, int max){
		for(AbstractKmerTable table : tables){
			table.countGC(gcCounts, max);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public int increment(final long kmer, final int incr) {
		throw new RuntimeException("Unsupported");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Array of underlying kmer tables that store the actual data */
	private final AbstractKmerTable[] tables;
	/** Maximum size of each buffer before forced flushing */
	private final int buflen;
	/** Half of buffer length used as soft threshold for conditional flushing */
	private final int halflen;
	/** Number of underlying tables and corresponding buffers */
	private final int ways;
	/** Whether to use value buffers for storing kmer-associated values */
	private final boolean useValues;
	/** Array of individual buffers corresponding to each underlying table */
	private final KmerBuffer[] buffers;
	/** Bit mask for core kmer hashing that excludes terminal bases */
	private final long coreMask;
	/** Bit mask for middle region hashing to improve distribution */
	private final long middleMask;
	/** Combined core and middle mask for final kmer hashing */
	private final long cmMask;
	/** Total count of unique kmers added across all buffer operations */
	public long uniqueAdded=0;
	
	/** Bit mask used for periodic buffer size checking optimization */
	private static final int SIZEMASK=15;
	/** Whether to use conditional setting mode when flushing buffers */
	private final boolean setIfNotPresent;
	
	/**
	 * Whether to sort buffer contents before flushing for duplicate consolidation
	 */
	public static boolean SORT_BUFFERS=false;

}
