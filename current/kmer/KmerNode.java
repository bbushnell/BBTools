package kmer;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.TextStreamWriter;
import shared.Tools;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * Abstract binary search tree node for storing k-mer counts and associated data.
 * Each node stores a k-mer pivot and subclass-defined value storage, with ownership
 * tracking and in-place rebalancing of the child links.
 * Child links and value storage are mutable and exposed through live references.
 * Tree updates require exclusive access or the enclosing table's lock. Monitor-guarded
 * owner writes do not synchronize the plain reads, structural updates or value mutations.
 * @author Brian Bushnell
 * @date Oct 22, 2013
 */
public abstract class KmerNode extends AbstractKmerTable{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Stores the pivot; child links remain null and ownership starts unclaimed. */
	protected KmerNode(long pivot_){
		pivot=pivot_;
	}

	/** Creates a detached node of the concrete value-storage family. */
	public abstract KmerNode makeNode(long pivot_, int value_);
	/** Creates a multi-value node where supported; the concrete type may retain the array. */
	public abstract KmerNode makeNode(long pivot_, int[] values_, int vlen_);

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Count-oriented insertion/update through the concrete node's scalar accessors.
	 * In single-value nodes, positive-increment overflow saturates at Integer.MAX_VALUE.
	 * @return Current scalar value, or the supplied increment when a child is created */
	@Override
	public final int increment(final long kmer, final int incr){
		if(pivot<0){pivot=kmer; return set(incr);} //Allows initializing empty nodes to -1
		if(kmer<pivot){
			if(left==null){left=makeNode(kmer, incr); return incr;}
			return left.increment(kmer, incr);
		}else if(kmer>pivot){
			if(right==null){right=makeNode(kmer, incr); return incr;}
			return right.increment(kmer, incr);
		}else{
			//KMR002: match primary count cells by saturating positive-increment overflow in victim nodes.
			final int count=value();
			if(count<Integer.MAX_VALUE){
				final int sum=count+incr;
				set(sum<0 ? Integer.MAX_VALUE : sum);
			}
			return value();
		}
	}

	/** Applies increment through the concrete node's existing scalar semantics.
	 * Creation depends on key presence, not the resulting count; an existing zero-count
	 * key is still present. Requires exclusive access, like other tree mutations.
	 * @return 1 for an initialized empty pivot or new child, otherwise 0 */
	@Override
	public final int incrementAndReturnNumCreated(final long kmer, final int incr){
		//KMR008: locate the structural insertion/update point before increment changes its count or pivot.
		final KmerNode node=getNodeOrParent(kmer);
		final int created=(node.pivot<0 || node.pivot!=kmer ? 1 : 0);
		node.increment(kmer, incr);
		return created;
	}

//	public final int set_Test(final long kmer, final int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD()){
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
//	public final int setIfNotPresent_Test(long kmer, int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD()){
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


	/** Locates or creates an entry, delegating replacement versus append semantics to its node type.
	 * @return 1 for a new child or initialized empty pivot, otherwise 0 */
	@Override
	public final int set(long kmer, int value){
		if(verbose){System.err.println("Set0: kmer="+kmer+", v="+value+", old="+Arrays.toString(values(new int[1])));}
		if(pivot<0){pivot=kmer; set(value); return 1;} //Allows initializing empty nodes to -1
		if(verbose){System.err.println("A");}
		if(kmer<pivot){
			if(verbose){System.err.println("B");}
			if(left==null){left=makeNode(kmer, value); return 1;}
			if(verbose){System.err.println("C");}
			return left.set(kmer, value);
		}else if(kmer>pivot){
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
	 * Initializes a scalar value for a k-mer only if its pivot is not already present.
	 * Creates a new node with the specified value if the k-mer is absent.
	 * @param kmer The k-mer to set
	 * @param value Initial scalar payload if the pivot is absent
	 * @return 1 for a new child or initialized empty pivot, otherwise 0
	 */
	@Override
	public final int setIfNotPresent(long kmer, int value){
		if(verbose){System.err.println("setIfNotPresent0: kmer="+kmer+", v="+value+", old="+Arrays.toString(values(new int[0])));}
		if(pivot<0){pivot=kmer; set(value); return 1;} //Allows initializing empty nodes to -1
		if(kmer<pivot){
			if(left==null){left=makeNode(kmer, value); return 1;}
			return left.setIfNotPresent(kmer, value);
		}else if(kmer>pivot){
			if(right==null){right=makeNode(kmer, value); return 1;}
			return right.setIfNotPresent(kmer, value);
		}
		return 0;
	}

	/**
	 * Retrieves the concrete node's primary scalar value for a k-mer.
	 * @param kmer The k-mer to look up
	 * @return Primary value (the count in single-value nodes), or -1 if absent
	 */
	@Override
	public final int getValue(long kmer){
		KmerNode n=get(kmer);
		return n==null ? -1 : n.value();
	}

	/**
	 * Retrieves values through the concrete node's borrowing/copying convention.
	 * @param kmer The k-mer to look up
	 * @param singleton One-element destination for single-value nodes; ignored by multi-value nodes
	 * @return Supplied single-value destination or live multi-value storage; null for absent/unallocated values
	 */
	@Override
	public final int[] getValues(long kmer, int[] singleton){
		KmerNode n=get(kmer);
		return n==null ? null : n.values(singleton);
	}

	/**
	 * Tests whether a k-mer is present in the tree.
	 * @param kmer The k-mer to search for
	 * @return true if the k-mer is present, false otherwise
	 */
	@Override
	public final boolean contains(long kmer){
		KmerNode node=get(kmer);
		return node!=null;
	}

	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the live left child, or null. */
	public KmerNode left(){return left;}
	/** Returns the live right child, or null. */
	public KmerNode right(){return right;}
	/** Returns the stored encoded pivot, including an empty-pivot marker. */
	public long pivot(){return pivot;}
	/** Reads the current owner marker without synchronization. */
	public int owner(){return owner;}

	/** Returns the concrete node's primary scalar value. */
	public int count(){return value();}
	/** Supplies the primary scalar value according to the node's storage convention. */
	protected abstract int value();
	/** Copies a scalar into singleton or lends backing multi-value storage. */
	protected abstract int[] values(int[] singleton);
	/**
	 * Replaces a scalar count or attempts to append a multi-value entry.
	 * @param value_ Scalar payload to store
	 * @return Supplied value, including when an append finds a duplicate
	 */
	public abstract int set(int value_);
	/** Adopts/merges multi-value storage where supported; single-value nodes reject it. */
	protected abstract int set(int[] values_, int vlen);

	/**
	 * Finds the node containing the specified k-mer.
	 * Walks left/right child links iteratively according to pivot comparisons.
	 * @param kmer The k-mer to search for
	 * @return The node containing the k-mer, or null if not found
	 */
	@Override
	final KmerNode get(long kmer){
//		if(kmer<pivot){
//			return left==null ? null : left.get(kmer);
//		}else if(kmer>pivot){
//			return right==null ? null : right.get(kmer);
//		}else{
//			return this;
//		}
		KmerNode n=this;
		while(n!=null && n.pivot!=kmer){
			n=(kmer<n.pivot ? n.left : n.right);
		}
		return n;
	}

	/** Returns the matching/empty node, or the parent whose appropriate child is absent. */
	final KmerNode getNodeOrParent(long kmer){
		if(pivot==kmer || pivot<0){return this;}
		if(kmer<pivot){return left==null ? this : left.getNodeOrParent(kmer);}
		return right==null ? this : right.getNodeOrParent(kmer);
	}

	/** Links an existing node at a vacant position, retaining any descendants.
	 * Caller must preserve distinct ordered pivots and avoid shared/cyclic subtrees. */
	final boolean insert(KmerNode n){
		assert(pivot!=-1);
		if(n.pivot<pivot){
			if(left==null){left=n; return true;}
			return left.insert(n);
		}else if(n.pivot>pivot){
			if(right==null){right=n; return true;}
			return right.insert(n);
		}else{
			return false;
		}
	}

	/** Appends nodes in ascending-pivot order (left, root, right), despite the legacy name. */
	final void traversePrefix(ArrayList<KmerNode> list){
		if(left!=null){left.traversePrefix(list);}
		list.add(this);
		if(right!=null){right.traversePrefix(list);}
	}

	/** Appends nodes in root-first order (root, left, right), despite the legacy name. */
	final void traverseInfix(ArrayList<KmerNode> list){
		list.add(this);
		if(left!=null){left.traverseInfix(list);}
		if(right!=null){right.traverseInfix(list);}
	}

	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts nodes with positive primary values, pruning each nonpositive node's subtree.
	 * @return Number of nodes counted under that pruning rule */
	@Override
	public final long size(){
		if(value()<1){return 0;}
		long size=1;
		if(left!=null){size+=left.size();}
		if(right!=null){size+=right.size();}
		return size;
	}

	/** Collects nodes in pivot order and rewires trees with more than two nodes.
	 * Clears the supplied empty work list; caller must install the returned root. */
	final KmerNode rebalance(ArrayList<KmerNode> list){
		assert(list.isEmpty());
		traversePrefix(list);
		KmerNode n=this;
		if(list.size()>2){
			n=rebalance(list, 0, list.size()-1);
		}
		list.clear();
		return n;
	}

	/** Reuses nodes from the inclusive, nonempty sorted range to build a balanced subtree. */
	private static final KmerNode rebalance(ArrayList<KmerNode> list, int a, int b){
		final int size=b-a+1;
		final int middle=a+size/2;
		final KmerNode n=list.get(middle);
		if(size<4){
			if(size==1){
				n.left=n.right=null;
			}else if(size==2){
				KmerNode n1=list.get(a);
				n.left=n1;
				n.right=null;
				n1.left=n1.right=null;
			}else{
				assert(size==3);
				KmerNode n1=list.get(a), n2=list.get(b);
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

	/** Node-level regeneration is unsupported; an enclosing table must manage removal. */
	@Override
	public long regenerate(final int limit){
		throw new RuntimeException("Not supported.");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/

	/** Formats through the concrete node traversal and submits text to a borrowed writer.
	 * @return Always true, without waiting for writer completion or reporting its errors */
	@Override
	public final boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
		tsw.print(dumpKmersAsText(new StringBuilder(32), k, mincount, maxcount));
		return true;
	}

	/** Appends through subclass traversal/filtering rules and returns the accumulation buffer. */
	protected abstract StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount);

	/** Byte-buffer counterpart of the subclass-specific text traversal. */
	protected abstract ByteBuilder dumpKmersAsText(ByteBuilder bb, int k, int mincount, int maxcount);

	/** Adds one node per primary-value bin, clamping at max and pruning nonpositive subtrees. */
	@Override
	public final void fillHistogram(long[] ca, int max){
		final int value=value();
		if(value<1){return;}
		ca[Tools.min(value, max)]++;
		if(left!=null){left.fillHistogram(ca, max);}
		if(right!=null){right.fillHistogram(ca, max);}
	}

	/** Appends positive primary values; nonpositive nodes suppress their descendants. */
	@Override
	public final void fillHistogram(SuperLongList sll){
		final int value=value();
		if(value<1){return;}
		sll.add(value);
		if(left!=null){left.fillHistogram(sll);}
		if(right!=null){right.fillHistogram(sll);}
	}

	/** Adds each key's G/C tally once to its primary-value bin; prunes nonpositive subtrees. */
	@Override
	public void countGC(long[] gcCounts, int max){
		final int value=value();
		if(value<1){return;}
		gcCounts[Tools.min(value, max)]+=gc(pivot);
		if(left!=null){left.countGC(gcCounts, max);}
		if(right!=null){right.countGC(gcCounts, max);}
	}

	/** Reports whether the concrete node uses multi-value storage. */
	abstract boolean TWOD();
	/** Returns the concrete node's logical entry count, not its backing-array capacity. */
	abstract int numValues();

	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/

	/** Resets every linked node's owner to -1 without filtering values or taking locks. */
	@Override
	public final void initializeOwnership(){
		owner=-1;
		if(left!=null){left.initializeOwnership();}
		if(right!=null){right.initializeOwnership();}
	}

	/** Resets all owner markers using initializeOwnership. */
	@Override
	public final void clearOwnership(){initializeOwnership();}

	/**
	 * Raises an existing node's owner under its monitor when the requested ID is larger.
	 * Lookup and preliminary/final owner reads are unsynchronized; callers coordinate
	 * tree changes and ownership initialization/clearing separately.
	 * @param kmer Existing k-mer whose ownership is to be updated
	 * @param newOwner The new owner identifier
	 * @return The final owner value after the operation
	 */
	@Override
	public final int setOwner(final long kmer, final int newOwner){
		KmerNode n=get(kmer);
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

	/** Clears a matching existing node's owner under its monitor; lookup itself is unlocked. */
	@Override
	public final boolean clearOwner(final long kmer, final int owner){
		KmerNode n=get(kmer);
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
	 * Reads an existing node's owner without synchronization.
	 * @param kmer Existing k-mer whose owner is requested
	 * @return The owner identifier, or -1 if unowned
	 */
	@Override
	public final int getOwner(final long kmer){
		KmerNode n=get(kmer);
		assert(n!=null);
		return n.owner;
	}

	/** Estimates 48 bytes per linked node; excludes subclass value arrays and is not a heap measurement. */
	public long calcMem(){//48 is a guess
		return 48+(left==null ? 0 : left.calcMem())+(right==null ? 0 : right.calcMem());
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Encoded key, or a negative pivot awaiting direct-node initialization. */
	long pivot;
	/** Plain owner marker; -1 is unowned. */
	int owner=-1;
	/** Live ordered child links, rewired by insertion and rebalancing. */
	KmerNode left, right;

}
