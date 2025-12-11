package kmer;

import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/**
 * One-dimensional k-mer node for storing a single integer value with a k-mer pivot.
 * Part of a binary search tree structure used for k-mer frequency tracking and analysis.
 * Optimized for scenarios where each k-mer requires only one value (typically count).
 *
 * @author Brian Bushnell
 * @date Oct 22, 2013
 */
public class KmerNode1D extends KmerNode {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public KmerNode1D(long pivot_){
		super(pivot_);
	}
	
	public KmerNode1D(long pivot_, int value_){
		super(pivot_);
		value=value_;
	}
	
	@Override
	public final KmerNode makeNode(long pivot_, int value_){
		return new KmerNode1D(pivot_, value_);
	}
	
	@Override
	public final KmerNode makeNode(long pivot_, int[] values_, int vlen){
		throw new RuntimeException("Unimplemented");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final int set(long kmer, int[] vals, int vlen) {
		throw new RuntimeException("Unimplemented.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
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
	
	@Override
	protected int set(int[] values_, int vlen){
		throw new RuntimeException("Unimplemented");
	}
	
	/** Returns the number of valid values stored in this node.
	 * @return 1 if value is positive, 0 if value is less than 1 */
	@Override
	int numValues(){return value<1 ? 0 : 1;}
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Indicates whether this node type supports resizing operations.
	 * @return Always false for KmerNode1D */
	@Override
	boolean canResize() {
		return false;
	}
	
	/** Indicates whether this node type supports rebalancing operations.
	 * @return Always true for KmerNode1D */
	@Override
	public boolean canRebalance() {
		return true;
	}
	
	@Deprecated
	@Override
	public int arrayLength() {
		throw new RuntimeException("Unsupported.");
	}
	
	@Deprecated
	@Override
	void resize() {
		throw new RuntimeException("Unsupported.");
	}
	
	@Deprecated
	@Override
	public void rebalance() {
		throw new RuntimeException("Please call rebalance(ArrayList<KmerNode>) instead, with an empty list.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Dumps k-mers from this node and its subtree as bytes to a stream writer.
	 * Performs in-order traversal, outputting k-mers that meet count criteria.
	 *
	 * @param bsw Stream writer for output
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for inclusion
	 * @param maxcount Maximum count threshold for inclusion
	 * @param remaining Atomic counter for limiting output quantity
	 * @return Always true
	 */
	@Override
	public final boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		if(value<1){return true;}
		if(value>=mincount){
			if(remaining!=null && remaining.decrementAndGet()<0){return true;}
			bsw.printlnKmer(pivot, value, k);
		}
		if(left!=null){left.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);}
		if(right!=null){right.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);}
		return true;
	}
	
	/**
	 * Multi-threaded version of k-mer dumping with buffered output.
	 * Uses ByteBuilder for efficient batching before writing to stream.
	 * Automatically flushes buffer when it reaches 16KB for memory efficiency.
	 *
	 * @param bsw Stream writer for synchronized output
	 * @param bb Buffer for accumulating k-mer data
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for inclusion
	 * @param maxcount Maximum count threshold for inclusion
	 * @param remaining Atomic counter for limiting output quantity
	 * @return Always true
	 */
	@Override
	public final boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining){
		if(value<1){return true;}
		if(value>=mincount){
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
	 * Performs in-order traversal, appending k-mers meeting count criteria.
	 *
	 * @param sb StringBuilder to accumulate text output
	 * @param k K-mer length for text formatting
	 * @param mincount Minimum count threshold for inclusion
	 * @param maxcount Maximum count threshold for inclusion
	 * @return The provided StringBuilder with appended k-mer text
	 */
	@Override
	protected final StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount){
		if(value<1){return sb;}
		if(sb==null){sb=new StringBuilder(32);}
		if(value>=mincount){sb.append(AbstractKmerTable.toText(pivot, value, k)).append('\n');}
		if(left!=null){left.dumpKmersAsText(sb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(sb, k, mincount, maxcount);}
		return sb;
	}
	
	/**
	 * Recursively dumps k-mers from this node and subtree using ByteBuilder.
	 * Performs in-order traversal, appending k-mers meeting count criteria.
	 *
	 * @param bb ByteBuilder to accumulate output
	 * @param k K-mer length for formatting
	 * @param mincount Minimum count threshold for inclusion
	 * @param maxcount Maximum count threshold for inclusion
	 * @return The provided ByteBuilder with appended k-mer data
	 */
	@Override
	protected final ByteBuilder dumpKmersAsText(ByteBuilder bb, int k, int mincount, int maxcount){
		if(value<1){return bb;}
		if(bb==null){bb=new ByteBuilder(32);}
		if(value>=mincount){bb.append(AbstractKmerTable.toBytes(pivot, value, k)).append('\n');}
		if(left!=null){left.dumpKmersAsText(bb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(bb, k, mincount, maxcount);}
		return bb;
	}
	
	/** Indicates whether this node supports two-dimensional value storage.
	 * @return Always false for one-dimensional node implementation */
	@Override
	final boolean TWOD(){return false;}
	
	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	int value;
	
}
