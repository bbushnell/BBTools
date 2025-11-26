package ukmer;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.TextStreamWriter;
import shared.Tools;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * @author Brian Bushnell
 * @date Oct 22, 2013
 *
 */
public abstract class KmerNodeU extends AbstractKmerTableU {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Constructs a KmerNodeU with the specified pivot k-mer.
	 * Clones the pivot array to prevent external modification.
	 * @param pivot_ The k-mer that serves as the pivot for this node
	 */
	protected KmerNodeU(long[] pivot_){
		pivot=pivot_.clone();
	}
	
	/**
	 * Factory method to create a new node with specified pivot and value.
	 * Implementation varies by concrete subclass to create appropriate node type.
	 *
	 * @param pivot_ The k-mer pivot for the new node
	 * @param value_ The initial value for the new node
	 * @return A new KmerNodeU instance
	 */
	public abstract KmerNodeU makeNode(long[] pivot_, int value_);
	/**
	 * Factory method to create a new node with specified pivot and values array.
	 * Used for multi-dimensional k-mer storage implementations.
	 *
	 * @param pivot_ The k-mer pivot for the new node
	 * @param values_ Array of values for the new node
	 * @return A new KmerNodeU instance
	 */
	public abstract KmerNodeU makeNode(long[] pivot_, int[] values_);
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final int increment(Kmer kmer){return increment(kmer.key());}
	
	/**
	 * Increments the count for the specified k-mer in the binary search tree.
	 * Creates new nodes as needed for k-mers not already in the tree.
	 * Prevents overflow by capping count at Integer.MAX_VALUE.
	 *
	 * @param kmer The k-mer array to increment
	 * @return The new count after incrementing
	 */
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
	
	@Override
	public final int incrementAndReturnNumCreated(Kmer kmer){return incrementAndReturnNumCreated(kmer.key());}
	
	/**
	 * Increments a k-mer count and returns whether a new node was created.
	 * @param kmer The k-mer array to increment
	 * @return 1 if a new node was created (count was 1), 0 otherwise
	 */
	public final int incrementAndReturnNumCreated(long[] kmer) {
		int x=increment(kmer);
		return x==1 ? 1 : 0;
	}
	
	/** Returns number of nodes added */
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
	
	
	/** Returns number of nodes added */
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
	
	/**
	 * Retrieves the value associated with a k-mer.
	 * @param kmer The k-mer array to look up
	 * @return The value associated with the k-mer, or -1 if not found
	 */
	public final int getValue(long[] kmer){
		KmerNodeU n=get(kmer);
		return n==null ? -1 : n.value();
	}
	
	/**
	 * Retrieves the values array associated with a k-mer.
	 * Uses singleton array for single-value nodes to avoid allocation.
	 *
	 * @param kmer The k-mer array to look up
	 * @param singleton Reusable array for single values
	 * @return Array of values associated with the k-mer, or null if not found
	 */
	public final int[] getValues(long[] kmer, int[] singleton){
		KmerNodeU n=get(kmer);
		return n==null ? null : n.values(singleton);
	}
	
	/**
	 * Checks whether a k-mer is present in the tree.
	 * @param kmer The k-mer array to search for
	 * @return true if the k-mer is found, false otherwise
	 */
	public final boolean contains(long[] kmer){
		KmerNodeU node=get(kmer);
		return node!=null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Gets the left child node */
	public KmerNodeU left(){return left;}
	/** Gets the right child node */
	public KmerNodeU right(){return right;}
	/** Gets the pivot k-mer for this node */
	public long[] pivot(){return pivot;}
	/** Gets the owner identifier for this node */
	public int owner(){return owner;}
	
	/**
	 * Fills a Kmer object with this node's pivot k-mer data.
	 * Copies pivot array contents into the Kmer's internal arrays.
	 * @param x The Kmer object to fill with pivot data
	 * @return The filled Kmer object
	 */
	public Kmer fillKmer(Kmer x){
		assert(pivot!=null);
		long[] key=x.array1();
		for(int i=0; i<pivot.length; i++){
			key[i]=pivot[i];
		}
		x.fillArray2();
		return x;
	}
	
	/** Gets the count value for this node, same as value() */
	public int count(){return value();}
	/** Gets the primary value stored in this node */
	protected abstract int value();
	/**
	 * Gets all values stored in this node as an array.
	 * @param singleton Reusable array for single values to avoid allocation
	 * @return Array containing all values for this node
	 */
	protected abstract int[] values(int[] singleton);
	/** Returns new value */
	public abstract int set(int value_);
	/**
	 * Sets multiple values for this node.
	 * @param values_ Array of values to store in this node
	 * @return Implementation-specific return value
	 */
	protected abstract int set(int[] values_);
	
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
	
	/**
	 * Gets the node containing a k-mer, or its parent if the k-mer is not found.
	 * Used for insertion operations to find the insertion point.
	 * @param kmer The k-mer array to search for
	 * @return The node containing the k-mer, or the parent where it should be inserted
	 */
	final KmerNodeU getNodeOrParent(long[] kmer){
		final int cmp=compare(kmer, pivot);
		if(cmp==0){return this;}
		if(cmp<0){return left==null ? this : left.getNodeOrParent(kmer);}
		return right==null ? this : right.getNodeOrParent(kmer);
	}
	
	/**
	 * Inserts a node into the binary search tree.
	 * Maintains tree ordering based on k-mer comparison.
	 * @param n The node to insert
	 * @return true if insertion succeeded, false if k-mer already exists
	 */
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
	
	/**
	 * Performs in-order traversal of the tree, adding nodes to the list.
	 * Results in nodes sorted by k-mer value.
	 * @param list ArrayList to store traversed nodes in sorted order
	 */
	final void traversePrefix(ArrayList<KmerNodeU> list){
		if(left!=null){left.traversePrefix(list);}
		list.add(this);
		if(right!=null){right.traversePrefix(list);}
	}
	
	/**
	 * Performs pre-order traversal of the tree, adding nodes to the list.
	 * Visits current node before its children.
	 * @param list ArrayList to store traversed nodes in pre-order
	 */
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
	
	@Override
	public final long size() {
		if(value()<1){return 0;}
		long size=1;
		if(left!=null){size+=left.size();}
		if(right!=null){size+=right.size();}
		return size;
	}
	
	/**
	 * Rebalances the binary search tree to improve performance.
	 * Uses the provided list as workspace for tree reconstruction.
	 * @param list Temporary ArrayList for tree reconstruction (must be empty)
	 * @return Root of the rebalanced tree
	 */
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
	
	/**
	 * Recursively rebalances a range of nodes from a sorted list.
	 * Creates a balanced binary tree from the sorted node list.
	 *
	 * @param list Sorted list of nodes to rebalance
	 * @param a Start index of the range to rebalance
	 * @param b End index of the range to rebalance
	 * @return Root node of the rebalanced subtree
	 */
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
	
	@Override
	public long regenerate(final int limit){
		throw new RuntimeException("Not supported.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount) {
		tsw.print(dumpKmersAsText(new StringBuilder(32), k, mincount, maxcount));
		return true;
	}
	
	/**
	 * Appends k-mer text representation to StringBuilder.
	 * Implementation varies by concrete subclass for different formats.
	 *
	 * @param sb StringBuilder to append text to
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for inclusion
	 * @param maxcount Maximum count threshold for inclusion
	 * @return The StringBuilder with appended k-mer text
	 */
	protected abstract StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount);
	
	/**
	 * Appends k-mer text representation to ByteBuilder for efficient output.
	 * Implementation varies by concrete subclass for different formats.
	 *
	 * @param bb ByteBuilder to append text to
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for inclusion
	 * @param maxcount Maximum count threshold for inclusion
	 * @return The ByteBuilder with appended k-mer text
	 */
	protected abstract ByteBuilder dumpKmersAsText(ByteBuilder bb, int k, int mincount, int maxcount);
	
	@Override
	public final void fillHistogram(long[] ca, int max){
		final int value=value();
		if(value<1){return;}
		ca[Tools.min(value, max)]++;
		if(left!=null){left.fillHistogram(ca, max);}
		if(right!=null){right.fillHistogram(ca, max);}
	}
	
	@Override
	public final void fillHistogram(SuperLongList sll){
		final int value=value();
		if(value<1){return;}
		sll.add(value);
		if(left!=null){left.fillHistogram(sll);}
		if(right!=null){right.fillHistogram(sll);}
	}
	
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
	
	@Override
	public String toString(){return Arrays.toString(pivot);}

	/** Returns whether this node supports two-dimensional storage */
	abstract boolean TWOD();
	/** Returns the number of values stored in this node */
	abstract int numValues();
	
	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final void initializeOwnership(){
		owner=-1;
		if(left!=null){left.initializeOwnership();}
		if(right!=null){right.initializeOwnership();}
	}
	
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
	
	
	/**
	 * Clears ownership for a k-mer node if currently owned by specified owner.
	 * Thread-safe operation that only clears if owner matches.
	 *
	 * @param kmer The k-mer to clear ownership for
	 * @param owner The expected current owner ID
	 * @return true if ownership was cleared, false if owner didn't match
	 */
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
	
	
	/**
	 * Gets the current owner ID for a k-mer node.
	 * @param kmer The k-mer to get ownership for
	 * @return The owner ID, or -1 if unowned
	 */
	public final int getOwner(final long[] kmer){
		KmerNodeU n=get(kmer);
		assert(n!=null);
		return n.owner;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Recall Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Sets multiple values for a k-mer in the tree.
	 * Implementation varies by concrete subclass.
	 *
	 * @param kmer The k-mer to set values for
	 * @param vals Array of values to set
	 * @return Implementation-specific return value
	 */
	abstract int set(long[] kmer, int[] vals);
	
	@Override
	public int set(Kmer kmer, int value) {
		return set(kmer.key(), value);
	}
	
	@Override
	public int set(Kmer kmer, int[] vals) {
		return set(kmer.key(), vals);
	}
	
	@Override
	public int setIfNotPresent(Kmer kmer, int value) {
		return setIfNotPresent(kmer.key(), value);
	}
	
	@Override
	public int getValue(Kmer kmer) {
		return getValue(kmer.key());
	}
	
	@Override
	public int[] getValues(Kmer kmer, int[] singleton) {
		return getValues(kmer.key(), singleton);
	}
	
	@Override
	public boolean contains(Kmer kmer) {
		return contains(kmer.key());
	}
	
	@Override
	public int getValue(long[] key, long xor) {
		return getValue(key);
	}
	
	@Override
	public int setOwner(Kmer kmer, int newOwner) {
		return setOwner(kmer.key(), newOwner);
	}
	
	@Override
	public boolean clearOwner(Kmer kmer, int owner) {
		return clearOwner(kmer.key(), owner);
	}
	
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
	
	/** The k-mer that serves as the pivot for binary search tree ordering */
	final long[] pivot;
	/** Owner identifier for thread-safe concurrent access, -1 indicates unowned */
	int owner=-1;
	KmerNodeU left, right;
	
}
