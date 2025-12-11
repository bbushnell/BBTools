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
 * Hash table implementation using binary search trees for k-mer storage.
 * Each hash bucket stores a binary search tree of KmerNodeU entries keyed
 * by k-mer, with optional 2D value storage and automatic resizing.
 * @author Brian Bushnell
 * @date Oct 23, 2013
 */
public final class HashForestU extends AbstractKmerTableU implements Iterable<KmerNodeU> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
//	public HashForestU(int initialSize, boolean autoResize_){
//		this(initialSize, autoResize_, false);
//	}
	
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
	
	private KmerNodeU makeNode(Kmer kmer, int val){return makeNode(kmer.key(), val);}
	private KmerNodeU makeNode(Kmer kmer, int[] vals){return makeNode(kmer.key(), vals);}
	
	private KmerNodeU makeNode(long[] kmer, int val){
		return (TWOD ? new KmerNodeU2D(kmer, val) : new KmerNodeU1D(kmer, val));
	}
	
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
	
	/**
	 * Increments the count for a k-mer by 1, creating a new node if needed and
	 * triggering resize when the load factor exceeds the configured limit.
	 * @param kmer K-mer to increment
	 * @return New count value after increment
	 */
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
	
	/**
	 * Increments the count for a k-mer and reports whether a new node was
	 * created for it.
	 * @param kmer K-mer to increment
	 * @return 1 if a new node was created, 0 if an existing node was incremented
	 */
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
	
	
	/**
	 * Sets the value for a k-mer, creating a node if necessary, and updates the
	 * table size and resize threshold as needed.
	 * @param kmer K-mer to set value for
	 * @param value Value to associate with the k-mer
	 * @return 1 if a new node was created, 0 if an existing node was updated
	 */
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
	
	/**
	 * Sets an array of values for a k-mer, creating a node if necessary, and
	 * updates the table size and resize threshold as needed.
	 * @param kmer K-mer to set values for
	 * @param vals Array of values to associate with the k-mer
	 * @return 1 if a new node was created, 0 if an existing node was updated
	 */
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
	
	/**
	 * Sets the value for a k-mer only if it is not already present, leaving
	 * existing values unchanged.
	 * @param kmer K-mer to conditionally set value for
	 * @param value Value to associate with the k-mer
	 * @return 1 if a new node was created, 0 if the k-mer already existed
	 */
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
	
	/**
	 * Gets the integer value associated with a k-mer by looking up its key and
	 * precomputed xor hash.
	 * @param kmer K-mer to look up
	 * @return Associated value, or -1 if the k-mer is not found
	 */
	@Override
	public final int getValue(Kmer kmer){
		return getValue(kmer.key(), kmer.xor());
	}
	
//	int getValue(KmerNodeU n){
//		return getValue(n.pivot, n.xor());
//	}
	
	/**
	 * Gets the integer value associated with a raw k-mer key and xor hash,
	 * searching the appropriate hash bucket tree.
	 * @param key Raw k-mer key as long array
	 * @param xor XOR hash value for the k-mer
	 * @return Associated value, or -1 if the k-mer is not found
	 */
	@Override
	public int getValue(long[] key, long xor) {
		int cell=(int)(xor%prime);
		if(verbose){System.err.println("Looking in cell "+cell+": "+array[cell]);}
		KmerNodeU n=array[cell];
		return n==null ? -1 : n.getValue(key);
	}
	
	/**
	 * Unsupported object-valued lookup; this implementation throws a RuntimeException.
	 * @param key Raw k-mer key as long array
	 * @return Never returns normally
	 * @throws RuntimeException Always thrown as the method is not implemented
	 */
	@Override
	Object get(long[] key) {
		throw new RuntimeException("Unimplemented.");
	}
	
	/**
	 * Gets the array of values associated with a k-mer in 2D mode, optionally
	 * reusing a singleton buffer for single-value nodes.
	 * @param kmer K-mer to look up
	 * @param singleton Reusable array for single values to avoid allocation
	 * @return Array of values, or null if the k-mer is not found
	 */
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
	
	/** Initializes ownership tracking for all nodes in the table by resetting
	 * their owner metadata. */
	@Override
	public final void initializeOwnership(){
		for(KmerNodeU n : array){
			if(n!=null){n.initializeOwnership();}
		}
	}
	
	/**
	 * Clears ownership tracking for all nodes by reinitializing ownership state.
	 */
	@Override
	public final void clearOwnership(){initializeOwnership();}
	
	/**
	 * Sets the owner ID for a k-mer node in the ownership fields.
	 * @param kmer K-mer whose ownership to set
	 * @param newOwner ID of the new owner
	 * @return Previous owner ID
	 */
	@Override
	public final int setOwner(final Kmer kmer, final int newOwner){
		final int cell=kmer.mod(prime);
		KmerNodeU n=array[cell];
		assert(n!=null);
		return n.setOwner(kmer, newOwner);
	}
	
	/**
	 * Clears ownership of a k-mer node if it is currently owned by the specified
	 * owner.
	 * @param kmer K-mer whose ownership to clear
	 * @param owner Expected current owner ID
	 * @return true if ownership was successfully cleared
	 */
	@Override
	public final boolean clearOwner(final Kmer kmer, final int owner){
		final int cell=kmer.mod(prime);
		KmerNodeU n=array[cell];
		assert(n!=null);
		return n.clearOwner(kmer, owner);
	}
	
	/**
	 * Gets the current owner ID for a k-mer node.
	 * @param kmer K-mer to check ownership for
	 * @return Current owner ID
	 */
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
	
	/**
	 * Gets the node containing the specified k-mer by traversing the binary tree
	 * in the appropriate hash bucket.
	 * @param kmer K-mer to search for
	 * @return Node containing the k-mer, or null if not found
	 */
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
	
	public final KmerNodeU getNode(int cell){
		KmerNodeU n=array[cell];
		return n;
	}
	
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
	
	/** Returns whether the table supports resizing operations. */
	@Override
	boolean canResize() {return true;}
	
	/** Returns whether the table supports tree rebalancing operations. */
	@Override
	public boolean canRebalance() {return true;}
	
	/** Returns the total number of k-mers stored in the table. */
	@Override
	public long size() {return size;}
	
	/** Returns the length of the underlying hash table array. */
	@Override
	public int arrayLength() {return array.length;}
	
	/** Resizes the hash table to accommodate more entries by choosing a new prime
	 * size and rehashing all existing nodes when the load factor exceeds limits. */
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
	
	/** Rebalances all binary trees in the hash table, reducing tree depth to
	 * improve lookup performance. */
	@Override
	public void rebalance(){
		ArrayList<KmerNodeU> list=new ArrayList<KmerNodeU>(1000);
		for(int i=0; i<array.length; i++){
			if(array[i]!=null){array[i]=array[i].rebalance(list);}
		}
	}
	
	public void clear() {
		size=0;
		Arrays.fill(array, null);
	}
	
	/**
	 * Regenerates table contents up to the specified limit; currently not
	 * implemented and throws a RuntimeException.
	 * @param limit Maximum number of entries to regenerate
	 * @return Never returns normally
	 * @throws RuntimeException Always thrown as the method is not implemented
	 */
	@Override
	long regenerate(final int limit) {
		throw new RuntimeException("Not implemented.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Dumps k-mers as text to a stream writer, emitting only those with counts
	 * within the specified range.
	 * @param tsw Text stream writer for output
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @return true when dumping is complete
	 */
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
	
	/**
	 * Dumps k-mers as bytes to a stream writer, decrementing an optional counter
	 * to support early termination on large datasets.
	 * @param bsw Byte stream writer for output
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @param remaining Counter for limiting output size
	 * @return true when dumping is complete or the limit is reached
	 */
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
	
	/**
	 * Multi-threaded version of k-mer byte dumping that uses a thread-local
	 * ByteBuilder for efficient string construction.
	 * @param bsw Byte stream writer for output
	 * @param bb Thread-local byte builder for string construction
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for output
	 * @param maxcount Maximum count threshold for output
	 * @param remaining Counter for limiting output size
	 * @return true when dumping is complete or the limit is reached
	 */
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
	
	/**
	 * Fills a count histogram array with the frequency distribution of k-mers,
	 * accumulating counts up to the specified maximum.
	 * @param ca Count array to populate (index = count, value = frequency)
	 * @param max Maximum count value to include in histogram
	 */
	@Override
	public void fillHistogram(long[] ca, int max){
		for(int i=0; i<array.length; i++){
			KmerNodeU node=array[i];
			if(node!=null){
				node.fillHistogram(ca, max);
			}
		}
	}
	
	/** Fills a histogram using an expandable SuperLongList for large count ranges.
	 * @param sll SuperLongList to populate with count frequency data */
	@Override
	public void fillHistogram(SuperLongList sll){
		for(int i=0; i<array.length; i++){
			KmerNodeU node=array[i];
			if(node!=null){
				node.fillHistogram(sll);
			}
		}
	}
	
	/**
	 * Counts GC content distribution across all stored k-mers and increments the
	 * provided GC histogram.
	 * @param gcCounts Array to store GC content frequencies
	 * @param max Maximum GC count to track
	 */
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
	
	/**
	 * Returns an iterator over all nodes in the table by materializing them into
	 * an ArrayList via infix traversal.
	 * @return Iterator over all k-mer nodes
	 */
	@Override
	public Iterator<KmerNodeU> iterator() {
		return toList().iterator();
	}
	
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
	
	public KmerNodeU[] array() {return array;}
	
	KmerNodeU[] array;
	int prime;
	long size=0;
	long sizeLimit;
	final int k;
	final long coreMask;
	final boolean autoResize;
	final boolean TWOD;
	private final Lock lock=new ReentrantLock();
	
	/** Returns the reentrant lock used for thread synchronization. */
	@Override
	final Lock getLock(){return lock;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	final static int maxPrime=(int)Primes.primeAtMost(Integer.MAX_VALUE);
	final static float resizeMult=2.5f; //Resize by a minimum of this much
	final static float minLoadFactor=0.75f; //Resize by enough to get the load above this factor
	final static float maxLoadFactor=2.5f; //Resize by enough to get the load under this factor
	final static float minLoadMult=1/minLoadFactor;
	final static float maxLoadMult=1/maxLoadFactor;
	

	
}
