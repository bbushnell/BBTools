package kmer;

import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/**
 * One-dimensional k-mer node for storing a single integer value with a k-mer pivot.
 * Part of a binary search tree structure used for k-mer frequency tracking and analysis.
 * Optimized for scenarios where each k-mer requires only one value (typically count).
 * The inherited tree operations mutate this node and its descendants. Histogram
 * and output traversal stops at nodes with nonpositive counts; callers that mark
 * entries for deletion must restore table invariants before consuming those traversals.
 *
 * @author Brian Bushnell
 * @date Oct 22, 2013
 */
public class KmerNode1D extends KmerNode{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a detached node with count zero; inherited links initially remain null. */
	public KmerNode1D(long pivot_){
		super(pivot_);
	}

	/** Creates a detached node with the supplied integer value, without validating its sign. */
	public KmerNode1D(long pivot_, int value_){
		super(pivot_);
		value=value_;
	}

	/** Creates a detached single-value child for inherited tree insertion. */
	@Override
	public final KmerNode makeNode(long pivot_, int value_){
		return new KmerNode1D(pivot_, value_);
	}

	/** Rejects multi-value child construction for this single-value node type. */
	@Override
	public final KmerNode makeNode(long pivot_, int[] values_, int vlen){
		throw new RuntimeException("Unimplemented");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Rejects multi-value insertion, including insertion into descendant nodes. */
	@Override
	public final int set(long kmer, int[] vals, int vlen){
		throw new RuntimeException("Unimplemented.");
	}

	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns this node's stored integer, including any nonpositive marker. */
	@Override
	public int value(){return value;}

	/**
	 * Fills the provided array with this node's single value.
	 * Copies the stored value to the first position of the singleton array.
	 * @param singleton Array of length 1 to receive the value
	 * @return The provided singleton array with value at position 0
	 */
	@Override
	protected int[] values(int[] singleton){
		assert(singleton.length==1);
		singleton[0]=value;
		return singleton;
	}

	/**
	 * Sets the value stored in this node.
	 * @param value_ New value to store
	 * @return The newly set value
	 */
	@Override
	public int set(int value_){return value=value_;}

	/** Rejects replacement of the scalar value with a value array. */
	@Override
	protected int set(int[] values_, int vlen){
		throw new RuntimeException("Unimplemented");
	}

	/** Returns the number of valid values stored in this node.
	 * @return 1 if value is positive, 0 if value is less than 1 */
	@Override
	int numValues(){return value<1 ? 0 : 1;}

	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/

	/** Indicates whether this node type supports resizing operations.
	 * @return Always false for KmerNode1D */
	@Override
	boolean canResize(){return false;}

	/** Indicates whether this node type supports rebalancing operations.
	 * @return Always true for KmerNode1D */
	@Override
	public boolean canRebalance(){return true;}

	/** Nodes have no backing array; this inherited table operation is unsupported. */
	@Deprecated
	@Override
	public int arrayLength(){
		throw new RuntimeException("Unsupported.");
	}

	/** Nodes grow through links rather than array resizing; always throws. */
	@Deprecated
	@Override
	void resize(){
		throw new RuntimeException("Unsupported.");
	}

	/** Requires the inherited rebalance overload that accepts an empty work list. */
	@Deprecated
	@Override
	public void rebalance(){
		throw new RuntimeException("Please call rebalance(ArrayList<KmerNode>) instead, with an empty list.");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Dumps k-mers from this node and its subtree as bytes to a stream writer.
	 * Visits the current node before left and right subtrees, within inclusive count bounds.
	 * A nonpositive node suppresses its subtree; a depleted remaining counter also stops
	 * traversal below the matching node. The writer remains owned by the caller.
	 *
	 * @param bsw Stream writer for output
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for inclusion
	 * @param maxcount Maximum count threshold for inclusion
	 * @param remaining Optional shared output allowance, decremented before each emitted record
	 * @return Always true; does not report asynchronous writer completion or errors
	 */
	@Override
	public final boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		if(value<1){return true;}
		//Comprehension [kmer/KmerNode1D#001 RESOLVED]: dumps now filter on BOTH bounds (value>=mincount && value<=maxcount),
		//matching HashArray dumps (which require readCellValue<=maxcount). Previously only mincount was applied here and in
		//_MT / both dumpKmersAsText variants; that was masked because all callers pass maxcount=Integer.MAX_VALUE, but a
		//real maxcount would have leaked high-count kmers from HashForest-stored entries (standalone or HashArray victims).
		if(value>=mincount && value<=maxcount){
			if(remaining!=null && remaining.decrementAndGet()<0){return true;}
			bsw.printlnKmer(pivot, value, k);
		}
		if(left!=null){left.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);}
		if(right!=null){right.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);}
		return true;
	}

	/**
	 * Traverses current, left and right nodes using a caller-owned buffer.
	 * Copies a buffer of at least 16000 bytes into a writer job and clears it;
	 * only job submission is synchronized. The caller submits any final partial buffer
	 * and must prevent concurrent mutation of the tree or this buffer.
	 *
	 * @param bsw Stream writer for synchronized output
	 * @param bb Buffer for accumulating k-mer data
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for inclusion
	 * @param maxcount Maximum count threshold for inclusion
	 * @param remaining Optional shared output allowance, decremented before each emitted record
	 * @return Always true; does not report asynchronous writer completion or errors
	 */
	@Override
	public final boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining){
		if(value<1){return true;}
		if(value>=mincount && value<=maxcount){
			if(remaining!=null && remaining.decrementAndGet()<0){return true;}
			toBytes(pivot, value, k, bb);
			bb.nl();
			if(bb.length()>=16000){
				ByteBuilder bb2=new ByteBuilder(bb);
				synchronized(bsw){bsw.addJob(bb2);}
				bb.clear();
			}
		}
		if(left!=null){left.dumpKmersAsBytes_MT(bsw, bb, k, mincount, maxcount, remaining);}
		if(right!=null){right.dumpKmersAsBytes_MT(bsw, bb, k, mincount, maxcount, remaining);}
		return true;
	}

	/**
	 * Recursively dumps k-mers from this node and subtree as text.
	 * Visits current, left and right nodes, appending entries within inclusive count bounds.
	 * Nonpositive nodes suppress their subtrees.
	 *
	 * @param sb Accumulation buffer, or null to allocate one at the first positive node
	 * @param k K-mer length for text formatting
	 * @param mincount Minimum count threshold for inclusion
	 * @param maxcount Maximum count threshold for inclusion
	 * @return The accumulation buffer, possibly null when this node is nonpositive
	 */
	@Override
	protected final StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount){
		if(value<1){return sb;}
		if(sb==null){sb=new StringBuilder(32);}
		if(value>=mincount && value<=maxcount){sb.append(AbstractKmerTable.toText(pivot, value, k)).append('\n');}
		if(left!=null){left.dumpKmersAsText(sb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(sb, k, mincount, maxcount);}
		return sb;
	}

	/**
	 * Recursively dumps k-mers from this node and subtree using ByteBuilder.
	 * Visits current, left and right nodes, appending entries within inclusive count bounds.
	 * Nonpositive nodes suppress their subtrees.
	 *
	 * @param bb Accumulation buffer, or null to allocate one at the first positive node
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for inclusion
	 * @param maxcount Maximum count threshold for inclusion
	 * @return The accumulation buffer, possibly null when this node is nonpositive
	 */
	@Override
	protected final ByteBuilder dumpKmersAsText(ByteBuilder bb, int k, int mincount, int maxcount){
		if(value<1){return bb;}
		if(bb==null){bb=new ByteBuilder(32);}
		if(value>=mincount && value<=maxcount){bb.append(AbstractKmerTable.toBytes(pivot, value, k)).append('\n');}
		if(left!=null){left.dumpKmersAsText(bb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(bb, k, mincount, maxcount);}
		return bb;
	}

	/** Indicates whether this node supports two-dimensional value storage.
	 * @return Always false for one-dimensional node implementation */
	@Override
	final boolean TWOD(){return false;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Scalar payload; nonpositive values are skipped by output and histogram traversal. */
	int value;

}
