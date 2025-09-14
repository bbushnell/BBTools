package ukmer;

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
	
	@Override
	public final long flush(){
		long added=0;
		for(int i=0; i<ways; i++){added+=dumpBuffer(i, true);}
		return added;
	}
	
	@Override
	public int set(Kmer kmer, int value) {
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}
	
	@Override
	public int set(Kmer kmer, int[] vals) {
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}
	
	@Override
	public int setIfNotPresent(Kmer kmer, int value) {
		throw new RuntimeException("Unimplemented method; this class lacks value buffers");
	}
	
	@Override
	public int getValue(Kmer kmer) {
		final int way=getWay(kmer);
		return tables[way].getValue(kmer);
	}
	
	@Override
	public int getValue(long[] key, long xor) {
		final int way=(int)(xor%ways);
		return tables[way].getValue(key, xor);
	}
	
	@Override
	public int[] getValues(Kmer kmer, int[] singleton){
		final int way=getWay(kmer);
		return tables[way].getValues(kmer, singleton);
	}
	
	@Override
	public boolean contains(Kmer kmer) {
		final int way=getWay(kmer);
		return tables[way].contains(kmer);
	}
	
	public final int getWay(Kmer kmer){return (int)(kmer.xor()%ways);}
	
	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final void initializeOwnership(){
		for(AbstractKmerTableU t : tables){t.initializeOwnership();}
	}
	
	@Override
	public final void clearOwnership(){
		for(AbstractKmerTableU t : tables){t.clearOwnership();}
	}
	
	@Override
	public final int setOwner(final Kmer kmer, final int newOwner){
		final int way=getWay(kmer);
		return tables[way].setOwner(kmer, newOwner);
	}
	
	@Override
	public final boolean clearOwner(final Kmer kmer, final int owner){
		final int way=getWay(kmer);
		return tables[way].clearOwner(kmer, owner);
	}
	
	@Override
	public final int getOwner(final Kmer kmer){
		final int way=getWay(kmer);
		return tables[way].getOwner(kmer);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	Object get(Kmer kmer) {
		final int way=getWay(kmer);
		return tables[way].get(kmer);
	}
	
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
		for(AbstractKmerTableU table : tables){
			sum+=table.regenerate(limit);
		}
		return sum;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
		for(AbstractKmerTableU table : tables){
			table.dumpKmersAsText(tsw, k, mincount, maxcount);
		}
		return true;
	}
	
	@Override
	public boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		for(AbstractKmerTableU table : tables){
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
		for(AbstractKmerTableU table : tables){
			table.countGC(gcCounts, max);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
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
