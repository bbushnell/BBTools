package ukmer;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.TextStreamWriter;
import shared.Tools;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * Binary-search-tree node used by HashForestU victim tables.
 * Each node retains a cloned long-kmer pivot, exposes that live storage through
 * public pivot(), links same-family children by the current
 * AbstractKmerTableU.compare order and delegates value storage to 1D/2D
 * subclasses. UKM106 tracks the signed-word comparator's order defect. Mutating
 * callers provide any required table-level synchronization.
 *
 * @author Brian Bushnell
 * @date Oct 22, 2013
 */
public abstract class KmerNodeU extends AbstractKmerTableU{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Copies the supplied key words so future caller reuse cannot retarget this node. */
	protected KmerNodeU(long[] pivot_){
		pivot=pivot_.clone();
	}

	/** Creates a same-family child node for one value. */
	public abstract KmerNodeU makeNode(long[] pivot_, int value_);
	/** Creates a same-family child node for one value-list payload. */
	public abstract KmerNodeU makeNode(long[] pivot_, int[] values_);

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Increments the count for a wrapped key.
	 * @param kmer The k-mer to increment
	 * @return count reported by the matching or newly created node
	 */
	@Override
	public final int increment(Kmer kmer){return increment(kmer.key());}

	/**
	 * Increments a raw key in this subtree, creating a child node initialized to one.
	 * Matching nodes attempt {@code set(value()+1)} only below saturation; 2D nodes
	 * may deduplicate that value rather than append another copy.
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

	/**
	 * Increments a wrapped key and reports creation using this node's scalar return.
	 * @param kmer The k-mer to increment
	 * @return 1 when the increment return is one, otherwise 0
	 */
	@Override
	public final int incrementAndReturnNumCreated(Kmer kmer){return incrementAndReturnNumCreated(kmer.key());}

	/**
	 * Increments a raw key and maps a returned count of one to creation.
	 * @param kmer raw retained key words
	 * @return 1 when {@link #increment(long[])} returns one, otherwise 0
	 */
	public final int incrementAndReturnNumCreated(long[] kmer){
		int x=increment(kmer);
		return x==1 ? 1 : 0;
	}

	/**
	 * Sets a scalar value for a raw key, creating a same-family child when absent.
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
	 * Sets a scalar value for a raw key only when that key is absent.
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

	/** Returns the primary value for a raw retained key, or -1 when absent. */
	public final int getValue(long[] kmer){
		KmerNodeU n=get(kmer);
		return n==null ? -1 : n.value();
	}

	/** Returns the full value list for a raw retained key, or null when absent. */
	public final int[] getValues(long[] kmer, int[] singleton){
		KmerNodeU n=get(kmer);
		return n==null ? null : n.values(singleton);
	}

	/** Returns true when this subtree contains a raw retained key. */
	public final boolean contains(long[] kmer){
		KmerNodeU node=get(kmer);
		return node!=null;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the left child, or null. */
	public KmerNodeU left(){return left;}
	/** Returns the right child, or null. */
	public KmerNodeU right(){return right;}
	/** Returns this node's live retained key words; callers must not mutate them. */
	public long[] pivot(){return pivot;}
	/** Returns this node's current owner id. */
	public int owner(){return owner;}

	/** Copies this node's retained key into a reusable Kmer and rebuilds its reverse strand. */
	public Kmer fillKmer(Kmer x){
		assert(pivot!=null);
		long[] key=x.array1();
		for(int i=0; i<pivot.length; i++){
			key[i]=pivot[i];
		}
		x.fillArray2();
		return x;
	}

	/** Returns this node's primary value. */
	public int count(){return value();}
	/** Returns this node's primary value. */
	protected abstract int value();
	/** Returns live values for list nodes, or fills singleton storage for 1D nodes. */
	protected abstract int[] values(int[] singleton);
	/**
	 * Stores or merges one scalar value in this node.
	 * @param value_ The value to store
	 * @return implementation-specific scalar return
	 */
	public abstract int set(int value_);
	/**
	 * Stores or merges a value list in this node.
	 * @param values_ value list, stopping at the first negative entry when merged
	 * @return implementation-specific insertion/update status
	 */
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

	/** Returns a matching node or the terminal parent to which a missing key would attach. */
	final KmerNodeU getNodeOrParent(long[] kmer){
		final int cmp=compare(kmer, pivot);
		if(cmp==0){return this;}
		if(cmp<0){return left==null ? this : left.getNodeOrParent(kmer);}
		return right==null ? this : right.getNodeOrParent(kmer);
	}

	/** Inserts a detached node into this subtree, returning false for an equal pivot. */
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

	//TODO: Possible bug [ukmer/KmerNodeU#001] - misleading names:
	//traversePrefix does L-self-R (INORDER), traverseInfix below does self-L-R
	//(PREORDER) - swapped vs standard BST terms. Usage currently correct
	//(rebalance needs sorted->uses this inorder one; resize/toList order-agnostic),
	//so LOW/latent, but a trap for future callers. UKM106 separately tracks the
	//signed-word compare order these traversals preserve.
	/** Adds nodes in left/self/right order despite the historical prefix name. */
	final void traversePrefix(ArrayList<KmerNodeU> list){
		if(left!=null){left.traversePrefix(list);}
		list.add(this);
		if(right!=null){right.traversePrefix(list);}
	}

	/** Adds nodes in self/left/right order despite the historical infix name. */
	final void traverseInfix(ArrayList<KmerNodeU> list){
		list.add(this);
		if(left!=null){left.traverseInfix(list);}
		if(right!=null){right.traverseInfix(list);}
	}

	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Counts this node and descendants only when this node's first value is positive;
	 * a first value below one prunes the entire current subtree.
	 * @return positive-first-value nodes reachable without hitting a pruned root
	 */
	@Override
	public final long size(){
		if(value()<1){return 0;}
		long size=1;
		if(left!=null){size+=left.size();}
		if(right!=null){size+=right.size();}
		return size;
	}

	/** Returns this subtree rebalanced from current left/self/right traversal order. */
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

	/** Recursively turns a traversal-ordered node range into a balanced subtree. */
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
	 * Appends implementation-specific text records to the writer.
	 *
	 * @param tsw Text stream writer for output
	 * @param k total kbig decoding length
	 * @param mincount Minimum count threshold for inclusion
	 * @param maxcount Maximum count threshold for inclusion
	 * @return true after local append; writer completion is caller-owned
	 */
	@Override
	public final boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
		tsw.print(dumpKmersAsText(new StringBuilder(32), k, mincount, maxcount));
		return true;
	}

	/** Appends qualifying subtree entries as text, using {@code k} as total length. */
	protected abstract StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount);

	/** Appends qualifying subtree entries as text bytes, using {@code k} as total length. */
	protected abstract ByteBuilder dumpKmersAsText(ByteBuilder bb, int k, int mincount, int maxcount);

	/**
	 * Adds this subtree's positive first values to a dense histogram; a first value
	 * below one prunes the entire current subtree.
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
	 * Adds this subtree's positive first values to a sparse histogram; a first value
	 * below one prunes the entire current subtree.
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
	 * Adds this subtree's pivot GC totals by positive first value; a first value
	 * below one prunes the entire current subtree.
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

	/** Returns a debug representation of this node's retained key words. */
	@Override
	public String toString(){return Arrays.toString(pivot);}

	/** Returns true for list-valued 2D nodes. */
	abstract boolean TWOD();
	/** Returns this node's logical value count. */
	abstract int numValues();

	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/

	/** Initializes ownership to unowned (-1) for this node and all children. */
	@Override
	public final void initializeOwnership(){
		owner=-1;
		if(left!=null){left.initializeOwnership();}
		if(right!=null){right.initializeOwnership();}
	}

	/** Clears ownership for this node and all children, same as initializeOwnership. */
	@Override
	public final void clearOwnership(){initializeOwnership();}


	/**
	 * Claims a raw key for a higher-numbered owner.
	 * @param kmer The k-mer to set ownership for
	 * @param newOwner candidate owner ID
	 * @return resulting owner ID
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
	 * Clears a raw key's owner only when it equals the supplied owner.
	 * @param kmer The k-mer key
	 * @param owner expected current owner ID
	 * @return true if ownership was cleared
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
	 * Returns the owner for a raw key.
	 * @param kmer The k-mer key
	 * @return owner ID, or -1 when unowned
	 */
	public final int getOwner(final long[] kmer){
		KmerNodeU n=get(kmer);
		assert(n!=null);
		return n.owner;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Recall Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Stores multiple values for a raw retained key. */
	abstract int set(long[] kmer, int[] vals);

	/** Delegates wrapped-key scalar set calls to the retained raw-key path. */
	@Override
	public int set(Kmer kmer, int value){
		return set(kmer.key(), value);
	}

	/**
	 * Sets multiple values for a k-mer using Kmer wrapper object.
	 * @param kmer The Kmer object containing the k-mer to set
	 * @param vals Array of values to assign
	 * @return Implementation-specific return value
	 */
	@Override
	public int set(Kmer kmer, int[] vals){
		return set(kmer.key(), vals);
	}

	/**
	 * Sets value for k-mer only if not already present, using Kmer wrapper.
	 * @param kmer The Kmer object containing the k-mer to set
	 * @param value The value to assign if k-mer is not present
	 * @return Number of nodes added
	 */
	@Override
	public int setIfNotPresent(Kmer kmer, int value){
		return setIfNotPresent(kmer.key(), value);
	}

	/**
	 * Gets the value for a k-mer using Kmer wrapper object.
	 * @param kmer The Kmer object containing the k-mer to look up
	 * @return The value associated with the k-mer, or -1 if not found
	 */
	@Override
	public int getValue(Kmer kmer){
		return getValue(kmer.key());
	}

	/**
	 * Gets values array for a k-mer using Kmer wrapper object.
	 * @param kmer The Kmer object containing the k-mer to look up
	 * @param singleton Reusable array for single values
	 * @return Array of values, or null if not found
	 */
	@Override
	public int[] getValues(Kmer kmer, int[] singleton){
		return getValues(kmer.key(), singleton);
	}

	/**
	 * Checks if a k-mer is present using Kmer wrapper object.
	 * @param kmer The Kmer object containing the k-mer to search for
	 * @return true if k-mer is found, false otherwise
	 */
	@Override
	public boolean contains(Kmer kmer){
		return contains(kmer.key());
	}

	/**
	 * Gets value for k-mer with XOR parameter (ignored in this implementation).
	 * @param key The k-mer array to look up
	 * @param xor XOR value (unused in tree implementation)
	 * @return The value associated with the k-mer, or -1 if not found
	 */
	@Override
	public int getValue(long[] key, long xor){
		return getValue(key);
	}

	/**
	 * Sets owner for k-mer using Kmer wrapper object.
	 * @param kmer The Kmer object containing the k-mer
	 * @param newOwner The new owner ID
	 * @return The actual owner ID after the operation
	 */
	@Override
	public int setOwner(Kmer kmer, int newOwner){
		return setOwner(kmer.key(), newOwner);
	}

	/**
	 * Clears owner for k-mer using Kmer wrapper object.
	 * @param kmer The Kmer object containing the k-mer
	 * @param owner The expected current owner ID
	 * @return true if ownership was cleared, false otherwise
	 */
	@Override
	public boolean clearOwner(Kmer kmer, int owner){
		return clearOwner(kmer.key(), owner);
	}

	/**
	 * Gets owner ID for k-mer using Kmer wrapper object.
	 * @param kmer The Kmer object containing the k-mer
	 * @return The owner ID, or -1 if unowned
	 */
	@Override
	public int getOwner(Kmer kmer){
		return getOwner(kmer.key());
	}


	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Live retained key words for this tree node; constructed as a caller-key clone. */
	final long[] pivot;
	/** Ownership metadata used by ownership sweep workers; -1 means unowned. */
	int owner=-1;
	/** Binary-search-tree children. */
	KmerNodeU left, right;

}
