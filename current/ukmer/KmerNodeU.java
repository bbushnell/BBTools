package ukmer;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.TextStreamWriter;
import shared.Tools;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * Abstract binary tree node for k-mer storage with tree-based organization.
 * Provides binary search tree functionality for k-mer data structures with support
 * for incremental counting, value storage, and tree balancing operations.
 *
 * @author Brian Bushnell
 * @date Oct 22, 2013
 */
public abstract class KmerNodeU extends AbstractKmerTableU {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	protected KmerNodeU(long[] pivot_){
		pivot=pivot_.clone();
	}
	
	public abstract KmerNodeU makeNode(long[] pivot_, int value_);
	public abstract KmerNodeU makeNode(long[] pivot_, int[] values_);
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Increments the count for the specified k-mer.
	 * Delegates to the long[] array version for actual processing.
	 * @param kmer The k-mer to increment
	 * @return The new count after incrementing
	 */
	@Override
	public final int increment(Kmer kmer){return increment(kmer.key());}
	
	public final int increment(long[] kmer){
		final int cmp=compare(kmer, pivot);
		if(cmp<0){
			if(left==null){left=makeNode(kmer, 1); return 1;}
			return left.increment(kmer);
		}else if(cmp>0){
			if(right==null){right=makeNode(kmer, 1); return 1;}
			return right.increment(kmer);
		}else{
			if(value()<Integer.MAX_VALUE){set(value()+1);}
			return value();
		}
	}
	
	/**
	 * Increments a k-mer count and returns whether a new node was created.
	 * Delegates to the long[] array version for actual processing.
	 * @param kmer The k-mer to increment
	 * @return 1 if a new node was created, 0 if existing node was incremented
	 */
	@Override
	public final int incrementAndReturnNumCreated(Kmer kmer){return incrementAndReturnNumCreated(kmer.key());}
	
	public final int incrementAndReturnNumCreated(long[] kmer) {
		int x=increment(kmer);
		return x==1 ? 1 : 0;
	}
	
	/**
	 * Sets the value for a k-mer in the binary search tree.
	 * Creates new nodes as needed for k-mers not already in the tree.
	 * Includes verbose debugging output when enabled.
	 *
	 * @param kmer The k-mer array to set
	 * @param value The value to assign to the k-mer
	 * @return Number of nodes added (1 if new node created, 0 if existing updated)
	 */
	public final int set(long[] kmer, int value){
		if(verbose){System.err.println("Set0: kmer="+Arrays.toString(kmer)+", v="+value+", old="+Arrays.toString(values(new int[1])));}
		if(verbose){System.err.println("A");}
		final int cmp=compare(kmer, pivot);
		if(cmp<0){
			if(verbose){System.err.println("B");}
			if(left==null){left=makeNode(kmer, value); return 1;}
			if(verbose){System.err.println("C");}
			return left.set(kmer, value);
		}else if(cmp>0){
			if(verbose){System.err.println("D");}
			if(right==null){right=makeNode(kmer, value); return 1;}
			if(verbose){System.err.println("E");}
			return right.set(kmer, value);
		}else{
			if(verbose){System.err.println("F");}
			set(value);
		}
		if(verbose){System.err.println("G");}
		return 0;
	}
	
	
	/**
	 * Sets the value for a k-mer only if it's not already present in the tree.
	 * Creates new nodes for k-mers not found, but leaves existing k-mers unchanged.
	 *
	 * @param kmer The k-mer array to set
	 * @param value The value to assign if k-mer is not present
	 * @return Number of nodes added (1 if new node created, 0 if k-mer already exists)
	 */
	public final int setIfNotPresent(long[] kmer, int value){
		if(verbose){System.err.println("setIfNotPresent0: kmer="+kmer+", v="+value+", old="+Arrays.toString(values(new int[0])));}
		final int cmp=compare(kmer, pivot);
		if(cmp<0){
			if(left==null){left=makeNode(kmer, value); return 1;}
			return left.setIfNotPresent(kmer, value);
		}else if(cmp>0){
			if(right==null){right=makeNode(kmer, value); return 1;}
			return right.setIfNotPresent(kmer, value);
		}
		return 0;
	}
	
	public final int getValue(long[] kmer){
		KmerNodeU n=get(kmer);
		return n==null ? -1 : n.value();
	}
	
	public final int[] getValues(long[] kmer, int[] singleton){
		KmerNodeU n=get(kmer);
		return n==null ? null : n.values(singleton);
	}
	
	public final boolean contains(long[] kmer){
		KmerNodeU node=get(kmer);
		return node!=null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	public KmerNodeU left(){return left;}
	public KmerNodeU right(){return right;}
	public long[] pivot(){return pivot;}
	public int owner(){return owner;}
	
	public Kmer fillKmer(Kmer x){
		assert(pivot!=null);
		long[] key=x.array1();
		for(int i=0; i<pivot.length; i++){
			key[i]=pivot[i];
		}
		x.fillArray2();
		return x;
	}
	
	public int count(){return value();}
	protected abstract int value();
	protected abstract int[] values(int[] singleton);
	/**
	 * Sets the primary value for this node.
	 * @param value_ The new value to store
	 * @return The new value after setting
	 */
	public abstract int set(int value_);
	protected abstract int set(int[] values_);
	
	/**
	 * Retrieves the node containing the specified k-mer using iterative search.
	 * More efficient than recursive approach for deep trees.
	 * @param kmer The k-mer array to search for
	 * @return The node containing the k-mer, or null if not found
	 */
	@Override
	final KmerNodeU get(final long[] kmer){
//		if(kmer<pivot){
//			return left==null ? null : left.get(kmer);
//		}else if(kmer>pivot){
//			return right==null ? null : right.get(kmer);
//		}else{
//			return this;
//		}
		KmerNodeU n=this;
		int cmp=compare(kmer, n.pivot);
		while(cmp!=0){
			n=(cmp<0 ? n.left : n.right);
			cmp=(n==null ? 0 : compare(kmer, n.pivot));
		}
		return n;
	}
	
	final KmerNodeU getNodeOrParent(long[] kmer){
		final int cmp=compare(kmer, pivot);
		if(cmp==0){return this;}
		if(cmp<0){return left==null ? this : left.getNodeOrParent(kmer);}
		return right==null ? this : right.getNodeOrParent(kmer);
	}
	
	final boolean insert(KmerNodeU n){
		assert(pivot!=null);
		final int cmp=compare(n.pivot, pivot);
		if(cmp<0){
			if(left==null){left=n; return true;}
			return left.insert(n);
		}else if(cmp>0){
			if(right==null){right=n; return true;}
			return right.insert(n);
		}else{
			return false;
		}
	}
	
	final void traversePrefix(ArrayList<KmerNodeU> list){
		if(left!=null){left.traversePrefix(list);}
		list.add(this);
		if(right!=null){right.traversePrefix(list);}
	}
	
	final void traverseInfix(ArrayList<KmerNodeU> list){
		list.add(this);
		if(left!=null){left.traverseInfix(list);}
		if(right!=null){right.traverseInfix(list);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Calculates the number of nodes in the subtree rooted at this node.
	 * Only counts nodes with value >= 1.
	 * @return Total number of nodes in this subtree
	 */
	@Override
	public final long size() {
		if(value()<1){return 0;}
		long size=1;
		if(left!=null){size+=left.size();}
		if(right!=null){size+=right.size();}
		return size;
	}
	
	final KmerNodeU rebalance(ArrayList<KmerNodeU> list){
		assert(list.isEmpty());
		traversePrefix(list);
		KmerNodeU n=this;
		if(list.size()>2){
			n=rebalance(list, 0, list.size()-1);
		}
		list.clear();
		return n;
	}
	
	private static final KmerNodeU rebalance(ArrayList<KmerNodeU> list, int a, int b){
		final int size=b-a+1;
		final int middle=a+size/2;
		final KmerNodeU n=list.get(middle);
		if(size<4){
			if(size==1){
				n.left=n.right=null;
			}else if(size==2){
				KmerNodeU n1=list.get(a);
				n.left=n1;
				n.right=null;
				n1.left=n1.right=null;
			}else{
				assert(size==3);
				KmerNodeU n1=list.get(a), n2=list.get(b);
				n.left=n1;
				n.right=n2;
				n1.left=n1.right=null;
				n2.left=n2.right=null;
			}
		}else{
			n.left=rebalance(list, a, middle-1);
			n.right=rebalance(list, middle+1, b);
		}
		return n;
	}
	
	/**
	 * Regeneration operation not supported for tree-based storage.
	 * @param limit Unused parameter
	 * @throws RuntimeException Always thrown as operation is not supported
	 */
	@Override
	public long regenerate(final int limit){
		throw new RuntimeException("Not supported.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Dumps k-mers as text to the specified writer.
	 * Delegates to abstract method for format-specific implementation.
	 *
	 * @param tsw Text stream writer for output
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for inclusion
	 * @param maxcount Maximum count threshold for inclusion
	 * @return true if dump completed successfully
	 */
	@Override
	public final boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount) {
		tsw.print(dumpKmersAsText(new StringBuilder(32), k, mincount, maxcount));
		return true;
	}
	
	protected abstract StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount);
	
	protected abstract ByteBuilder dumpKmersAsText(ByteBuilder bb, int k, int mincount, int maxcount);
	
	/**
	 * Fills a histogram array with k-mer count frequencies.
	 * Recursively processes all nodes in the subtree.
	 * @param ca Count array where index represents count value
	 * @param max Maximum count value to include in histogram
	 */
	@Override
	public final void fillHistogram(long[] ca, int max){
		final int value=value();
		if(value<1){return;}
		ca[Tools.min(value, max)]++;
		if(left!=null){left.fillHistogram(ca, max);}
		if(right!=null){right.fillHistogram(ca, max);}
	}
	
	/**
	 * Fills a SuperLongList with k-mer count values for histogram generation.
	 * Recursively processes all nodes in the subtree.
	 * @param sll SuperLongList to store count values
	 */
	@Override
	public final void fillHistogram(SuperLongList sll){
		final int value=value();
		if(value<1){return;}
		sll.add(value);
		if(left!=null){left.fillHistogram(sll);}
		if(right!=null){right.fillHistogram(sll);}
	}
	
	/**
	 * Counts GC content in k-mers, accumulating results by count value.
	 * Processes all k-mer components in the pivot array.
	 * @param gcCounts Array to accumulate GC counts, indexed by k-mer count
	 * @param max Maximum count value to process
	 */
	@Override
	public final void countGC(long[] gcCounts, int max){
		final int value=value();
		if(value<1){return;}
		int index=Tools.min(value, max);
		for(long x : pivot){
			gcCounts[index]+=gc(x);
		}
		if(left!=null){left.countGC(gcCounts, max);}
		if(right!=null){right.countGC(gcCounts, max);}
	}
	
	/** Returns string representation of the pivot k-mer array */
	@Override
	public String toString(){return Arrays.toString(pivot);}

	abstract boolean TWOD();
	abstract int numValues();
	
	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Initializes ownership to unowned (-1) for this node and all children */
	@Override
	public final void initializeOwnership(){
		owner=-1;
		if(left!=null){left.initializeOwnership();}
		if(right!=null){right.initializeOwnership();}
	}
	
	/**
	 * Clears ownership for this node and all children, same as initializeOwnership
	 */
	@Override
	public final void clearOwnership(){initializeOwnership();}
	
	
	/**
	 * Sets the owner for a k-mer node using thread-safe operations.
	 * Only updates owner if new owner ID is higher than current.
	 *
	 * @param kmer The k-mer to set ownership for
	 * @param newOwner The new owner ID
	 * @return The actual owner ID after the operation
	 */
	public final int setOwner(final long[] kmer, final int newOwner){
		KmerNodeU n=get(kmer);
		assert(n!=null);
		if(n.owner<=newOwner){
			synchronized(n){
				if(n.owner<newOwner){
					n.owner=newOwner;
				}
			}
		}
		return n.owner;
	}
	
	
	public final boolean clearOwner(final long[] kmer, final int owner){
		KmerNodeU n=get(kmer);
		assert(n!=null);
		synchronized(n){
			if(n.owner==owner){
				n.owner=-1;
				return true;
			}
		}
		return false;
	}
	
	
	public final int getOwner(final long[] kmer){
		KmerNodeU n=get(kmer);
		assert(n!=null);
		return n.owner;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Recall Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	abstract int set(long[] kmer, int[] vals);
	
	@Override
	public int set(Kmer kmer, int value) {
		return set(kmer.key(), value);
	}
	
	/**
	 * Sets multiple values for a k-mer using Kmer wrapper object.
	 * @param kmer The Kmer object containing the k-mer to set
	 * @param vals Array of values to assign
	 * @return Implementation-specific return value
	 */
	@Override
	public int set(Kmer kmer, int[] vals) {
		return set(kmer.key(), vals);
	}
	
	/**
	 * Sets value for k-mer only if not already present, using Kmer wrapper.
	 * @param kmer The Kmer object containing the k-mer to set
	 * @param value The value to assign if k-mer is not present
	 * @return Number of nodes added
	 */
	@Override
	public int setIfNotPresent(Kmer kmer, int value) {
		return setIfNotPresent(kmer.key(), value);
	}
	
	/**
	 * Gets the value for a k-mer using Kmer wrapper object.
	 * @param kmer The Kmer object containing the k-mer to look up
	 * @return The value associated with the k-mer, or -1 if not found
	 */
	@Override
	public int getValue(Kmer kmer) {
		return getValue(kmer.key());
	}
	
	/**
	 * Gets values array for a k-mer using Kmer wrapper object.
	 * @param kmer The Kmer object containing the k-mer to look up
	 * @param singleton Reusable array for single values
	 * @return Array of values, or null if not found
	 */
	@Override
	public int[] getValues(Kmer kmer, int[] singleton) {
		return getValues(kmer.key(), singleton);
	}
	
	/**
	 * Checks if a k-mer is present using Kmer wrapper object.
	 * @param kmer The Kmer object containing the k-mer to search for
	 * @return true if k-mer is found, false otherwise
	 */
	@Override
	public boolean contains(Kmer kmer) {
		return contains(kmer.key());
	}
	
	/**
	 * Gets value for k-mer with XOR parameter (ignored in this implementation).
	 * @param key The k-mer array to look up
	 * @param xor XOR value (unused in tree implementation)
	 * @return The value associated with the k-mer, or -1 if not found
	 */
	@Override
	public int getValue(long[] key, long xor) {
		return getValue(key);
	}
	
	/**
	 * Sets owner for k-mer using Kmer wrapper object.
	 * @param kmer The Kmer object containing the k-mer
	 * @param newOwner The new owner ID
	 * @return The actual owner ID after the operation
	 */
	@Override
	public int setOwner(Kmer kmer, int newOwner) {
		return setOwner(kmer.key(), newOwner);
	}
	
	/**
	 * Clears owner for k-mer using Kmer wrapper object.
	 * @param kmer The Kmer object containing the k-mer
	 * @param owner The expected current owner ID
	 * @return true if ownership was cleared, false otherwise
	 */
	@Override
	public boolean clearOwner(Kmer kmer, int owner) {
		return clearOwner(kmer.key(), owner);
	}
	
	/**
	 * Gets owner ID for k-mer using Kmer wrapper object.
	 * @param kmer The Kmer object containing the k-mer
	 * @return The owner ID, or -1 if unowned
	 */
	@Override
	public int getOwner(Kmer kmer) {
		return getOwner(kmer.key());
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	final long[] pivot;
	int owner=-1;
	KmerNodeU left, right;
	
}
