package ukmer;

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
 * @author Brian Bushnell
 * @date Oct 23, 2013
 *
 */
public final class HashForestU extends AbstractKmerTableU implements Iterable<KmerNodeU> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
//	public HashForestU(int initialSize, boolean autoResize_){
//		this(initialSize, autoResize_, false);
//	}
	
	/**
	 * Constructs a hash forest with specified parameters.
	 * Initializes hash table with prime-sized array and configures k-mer parameters.
	 *
	 * @param initialSize Initial size of hash table (will be adjusted to nearest prime)
	 * @param k_ K-mer length
	 * @param autoResize_ Whether to automatically resize when load factor exceeds limit
	 * @param twod_ Whether to use 2D value storage mode
	 */
	public HashForestU(int initialSize, int k_, boolean autoResize_, boolean twod_){
		if(initialSize>1){
			initialSize=(int)Tools.min(maxPrime, Primes.primeAtLeast(initialSize));
		}else{
			initialSize=1;
		}
		prime=initialSize;
		sizeLimit=(long) (initialSize*resizeMult);
		array=allocKmerNodeArray(prime);
		k=k_;
		coreMask=Kmer.toCoreMask(k);
		autoResize=autoResize_;
		TWOD=twod_;
	}
	
	/**
	 * Creates a new node with single integer value.
	 * @param kmer K-mer to store in the node
	 * @param val Integer value to associate with the k-mer
	 * @return New node containing the k-mer and value
	 */
	private KmerNodeU makeNode(Kmer kmer, int val){return makeNode(kmer.key(), val);}
	/**
	 * Creates a new node with array of integer values.
	 * @param kmer K-mer to store in the node
	 * @param vals Array of integer values to associate with the k-mer
	 * @return New node containing the k-mer and values
	 */
	private KmerNodeU makeNode(Kmer kmer, int[] vals){return makeNode(kmer.key(), vals);}
	
	/**
	 * Creates a new node with single integer value from raw k-mer key.
	 * Chooses 1D or 2D node type based on TWOD configuration.
	 *
	 * @param kmer Raw k-mer key as long array
	 * @param val Integer value to associate with the k-mer
	 * @return New node containing the k-mer and value
	 */
	private KmerNodeU makeNode(long[] kmer, int val){
		return (TWOD ? new KmerNodeU2D(kmer, val) : new KmerNodeU1D(kmer, val));
	}
	
	/**
	 * Creates a new 2D node with array of integer values from raw k-mer key.
	 * @param kmer Raw k-mer key as long array
	 * @param vals Array of integer values to associate with the k-mer
	 * @return New 2D node containing the k-mer and values
	 */
	private KmerNodeU makeNode(long[] kmer, int[] vals){
		assert(TWOD);
		return new KmerNodeU2D(kmer, vals);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public KmerNodeU findParent(Kmer kmer, final int cell){return findParent(kmer.key(), cell);}
	
	public KmerNodeU findParent(final long[] kmer, final int cell){
		KmerNodeU n=array[cell], prev=null;
		int cmp=(n==null ? 0 : compare(kmer, n.pivot()));
		while(cmp!=0){
			prev=n;
			n=(cmp<0 ? n.left : n.right);
			cmp=(n==null ? 0 : compare(kmer, n.pivot()));
		}
		return prev;
	}
	
	@Override
	public int increment(Kmer kmer){
		final int cell=kmer.mod(prime);
		KmerNodeU n=array[cell], prev=null;
		final long[] key=kmer.key();
		int cmp=(n==null ? 0 : compare(key, n.pivot()));
		while(cmp!=0){
			prev=n;
			n=(cmp<0 ? n.left : n.right);
			cmp=(n==null ? 0 : compare(key, n.pivot()));
		}
		if(n==null){
			n=makeNode(kmer, 1);
			size++;
			if(prev==null){
				array[cell]=n;
			}else{
				if(compare(key, prev.pivot)<0){
					prev.left=n;
				}else{
					prev.right=n;
				}
			}
			if(autoResize && size>sizeLimit){resize();}
		}else{
			n.increment(kmer);
		}
		return n.value();
	}
	
	@Override
	public int incrementAndReturnNumCreated(Kmer kmer){
//		assert(kmer.verify(false));
////		Kmer old=kmer.clone(); //123
////		System.err.println("cell should be "+kmer.mod(prime)+"; prime="+prime);
//		int a=getValue(kmer);
//		int x=incrementAndReturnNumCreated0(kmer);
////		System.err.println("cell should be "+kmer.mod(prime)+"; prime="+prime);
//		int b=getValue(kmer);
////		System.err.println("cell should be "+kmer.mod(prime)+"; prime="+prime);
////		assert(old.equals(kmer));
//		assert(Tools.max(a, 0)+1==b) : a+", "+b+", "+x+", "+kmer+", "+kmer.arraysToString();
//		return x;
//	}
//
//	public int incrementAndReturnNumCreated0(Kmer kmer){//123
		final int cell=kmer.mod(prime);
		if(verbose){System.err.println("Placed in cell "+cell+":  "+Arrays.toString(kmer.key()));}
//		assert(cell==kmer.xor()%prime);
		KmerNodeU n=array[cell], prev=null;
		final long[] key=kmer.key();
		int cmp=(n==null ? 0 : compare(key, n.pivot()));
		while(cmp!=0){
			prev=n;
			n=(cmp<0 ? n.left : n.right);
			cmp=(n==null ? 0 : compare(key, n.pivot()));
		}
		if(n==null){
			n=makeNode(kmer, 1);
			size++;
			if(prev==null){
				array[cell]=n;
			}else{
				if(compare(key, prev.pivot)<0){
					prev.left=n;
				}else{
					prev.right=n;
				}
			}
			if(autoResize && size>sizeLimit){resize();}
			return 1;
		}else{
			n.increment(kmer);
			return 0;
		}
	}
	
//	public final int set_Test(final long[] kmer, final int v){
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
//	public final int setIfNotPresent_Test(Kmer kmer, int v){
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
//	public final int set_Test(final long[] kmer, final int v[]){
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
	
	
	@Override
	public int set(Kmer kmer, int value){
		int x=1, cell=kmer.mod(prime);
		final KmerNodeU n=array[cell];
		if(n==null){
			array[cell]=makeNode(kmer, value);
		}else{
			x=n.set(kmer, value);
		}
		size+=x;
		if(autoResize && size>sizeLimit){resize();}
		return x;
	}
	
	@Override
	public int set(Kmer kmer, int[] vals) {
		int x=1, cell=kmer.mod(prime);
		final KmerNodeU n=array[cell];
		if(n==null){
			array[cell]=makeNode(kmer, vals);
		}else{
			x=n.set(kmer, vals);
		}
		size+=x;
		if(autoResize && size>sizeLimit){resize();}
		return x;
	}
	
	@Override
	public int setIfNotPresent(Kmer kmer, int value){
		int x=1, cell=kmer.mod(prime);
		final KmerNodeU n=array[cell];
		if(n==null){
			array[cell]=makeNode(kmer, value);
		}else{
			x=n.setIfNotPresent(kmer, value);
		}
		size+=x;
		if(autoResize && size>sizeLimit){resize();}
		return x;
	}
	
	@Override
	public final int getValue(Kmer kmer){
		return getValue(kmer.key(), kmer.xor());
	}
	
//	int getValue(KmerNodeU n){
//		return getValue(n.pivot, n.xor());
//	}
	
	@Override
	public int getValue(long[] key, long xor) {
		int cell=(int)(xor%prime);
		if(verbose){System.err.println("Looking in cell "+cell+": "+array[cell]);}
		KmerNodeU n=array[cell];
		return n==null ? -1 : n.getValue(key);
	}
	
	@Override
	Object get(long[] key) {
		throw new RuntimeException("Unimplemented.");
	}
	
	@Override
	public int[] getValues(Kmer kmer, int[] singleton){
		int cell=kmer.mod(prime);
		KmerNodeU n=array[cell];
		return n==null ? null : n.getValues(kmer, singleton);
	}
	
	@Override
	public boolean contains(Kmer kmer){
		return get(kmer)!=null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final void initializeOwnership(){
		for(KmerNodeU n : array){
			if(n!=null){n.initializeOwnership();}
		}
	}
	
	@Override
	public final void clearOwnership(){initializeOwnership();}
	
	@Override
	public final int setOwner(final Kmer kmer, final int newOwner){
		final int cell=kmer.mod(prime);
		KmerNodeU n=array[cell];
		assert(n!=null);
		return n.setOwner(kmer, newOwner);
	}
	
	@Override
	public final boolean clearOwner(final Kmer kmer, final int owner){
		final int cell=kmer.mod(prime);
		KmerNodeU n=array[cell];
		assert(n!=null);
		return n.clearOwner(kmer, owner);
	}
	
	@Override
	public final int getOwner(final Kmer kmer){
		final int cell=kmer.mod(prime);
		KmerNodeU n=array[cell];
		assert(n!=null);
		return n.getOwner(kmer);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	final KmerNodeU get(Kmer kmer){
		int cell=kmer.mod(prime);
		KmerNodeU n=array[cell];
		final long[] key=kmer.key();
		int cmp=(n==null ? 0 : compare(key, n.pivot()));
		while(cmp!=0){
			n=(cmp<0 ? n.left : n.right);
			cmp=(n==null ? 0 : compare(key, n.pivot()));
		}
		return n;
	}
	
	/**
	 * Gets the root node for a specific hash table cell.
	 * @param cell Hash table cell index
	 * @return Root node for the cell, or null if cell is empty
	 */
	public final KmerNodeU getNode(int cell){
		KmerNodeU n=array[cell];
		return n;
	}
	
	/**
	 * Inserts a node into the appropriate hash table cell.
	 * Clears existing left/right pointers and finds correct position in tree.
	 * @param n Node to insert
	 * @return true if insertion was successful
	 */
	boolean insert(KmerNodeU n){
		n.left=null;
		n.right=null;
		int cell=(int)(Kmer.xor(n.pivot(), coreMask)%prime);
		if(array[cell]==null){
			array[cell]=n;
			return true;
		}
		return array[cell].insert(n);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	boolean canResize() {return true;}
	
	@Override
	public boolean canRebalance() {return true;}
	
	@Override
	public long size() {return size;}
	
	@Override
	public int arrayLength() {return array.length;}
	
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
		KmerNodeU[] old=array;
		array=allocKmerNodeArray(prime2);
		ArrayList<KmerNodeU> list=new ArrayList<KmerNodeU>(1000);
		for(int i=0; i<old.length; i++){
			if(old[i]!=null){
				old[i].traverseInfix(list);
				for(KmerNodeU n : list){
					insert(n);
//					assert(getValue(n)==n.value());//123 slow
				}
				list.clear();
			}
		}
		sizeLimit=Tools.max((long)(size*1.4), (long)(maxLoadFactor*prime));
	}
	
	@Override
	public void rebalance(){
		ArrayList<KmerNodeU> list=new ArrayList<KmerNodeU>(1000);
		for(int i=0; i<array.length; i++){
			if(array[i]!=null){array[i]=array[i].rebalance(list);}
		}
	}
	
	/** Removes all k-mers from the table and resets size to zero */
	public void clear() {
		size=0;
		Arrays.fill(array, null);
	}
	
	@Override
	long regenerate(final int limit) {
		throw new RuntimeException("Not implemented.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
//		tsw.print("HashForest:\n");
		for(int i=0; i<array.length; i++){
			KmerNodeU node=array[i];
			if(node!=null && node.value()>=mincount){
//				StringBuilder sb=new StringBuilder();
//				tsw.print(node.dumpKmersAsText(sb, k, mincount, maxcount));
				node.dumpKmersAsText(tsw, k, mincount, maxcount);
			}
		}
		return true;
	}
	
	@Override
	public boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
//		tsw.print("HashForest:\n");
		for(int i=0; i<array.length; i++){
			KmerNodeU node=array[i];
			if(node!=null && node.value()>=mincount){
//				StringBuilder sb=new StringBuilder();
//				tsw.print(node.dumpKmersAsText(sb, k, mincount, maxcount));
				if(remaining!=null && remaining.decrementAndGet()<0){return true;}
				node.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);
			}
		}
		return true;
	}
	
	@Override
	public boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining){
		for(int i=0; i<array.length; i++){
			KmerNodeU node=array[i];
			if(node!=null && node.value()>=mincount){
				if(remaining!=null && remaining.decrementAndGet()<0){return true;}
				node.dumpKmersAsBytes_MT(bsw, bb, k, mincount, maxcount, remaining);
			}
		}
		return true;
	}
	
	@Override
	public void fillHistogram(long[] ca, int max){
		for(int i=0; i<array.length; i++){
			KmerNodeU node=array[i];
			if(node!=null){
				node.fillHistogram(ca, max);
			}
		}
	}
	
	@Override
	public void fillHistogram(SuperLongList sll){
		for(int i=0; i<array.length; i++){
			KmerNodeU node=array[i];
			if(node!=null){
				node.fillHistogram(sll);
			}
		}
	}
	
	@Override
	public void countGC(long[] gcCounts, int max){
		for(int i=0; i<array.length; i++){
			KmerNodeU node=array[i];
			if(node!=null){
				node.countGC(gcCounts, max);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Iteration           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public Iterator<KmerNodeU> iterator() {
		return toList().iterator();
	}
	
	/**
	 * Converts all nodes in the table to an ArrayList.
	 * Performs infix traversal of all binary trees.
	 * @return ArrayList containing all nodes in the table
	 */
	public ArrayList<KmerNodeU> toList(){
		assert(size<Integer.MAX_VALUE);
		ArrayList<KmerNodeU> list=new ArrayList<KmerNodeU>((int)size);
		for(int i=0; i<array.length; i++){
			if(array[i]!=null){array[i].traverseInfix(list);}
		}
		assert(list.size()==size);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns the underlying hash table array for direct access */
	public KmerNodeU[] array() {return array;}
	
	/** Hash table array where each cell contains root of binary search tree */
	KmerNodeU[] array;
	/** Current size of hash table array (always prime number) */
	int prime;
	/** Total number of k-mers currently stored in the table */
	long size=0;
	/** Size threshold that triggers automatic resizing when exceeded */
	long sizeLimit;
	/** K-mer length used for this table */
	final int k;
	/** Bit mask for extracting core k-mer bits during hashing */
	final long coreMask;
	/** Whether table automatically resizes when size limit is exceeded */
	final boolean autoResize;
	/** Whether table operates in 2D mode with multiple values per k-mer */
	final boolean TWOD;
	/** Reentrant lock for thread-safe operations */
	private final Lock lock=new ReentrantLock();
	
	@Override
	final Lock getLock(){return lock;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Largest prime number that fits in an integer */
	final static int maxPrime=(int)Primes.primeAtMost(Integer.MAX_VALUE);
	/** Minimum multiplier for table size when resizing */
	final static float resizeMult=2.5f; //Resize by a minimum of this much
	/** Target load factor after resizing to avoid immediate re-resize */
	final static float minLoadFactor=0.75f; //Resize by enough to get the load above this factor
	/** Maximum load factor before triggering resize operation */
	final static float maxLoadFactor=2.5f; //Resize by enough to get the load under this factor
	/** Inverse of minimum load factor for resize calculations */
	final static float minLoadMult=1/minLoadFactor;
	/** Inverse of maximum load factor for resize calculations */
	final static float maxLoadMult=1/maxLoadFactor;
	

	
}
