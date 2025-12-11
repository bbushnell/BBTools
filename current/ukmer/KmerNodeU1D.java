package ukmer;

import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/**
 * Single-dimensional k-mer node for optimized linear storage structures.
 * Provides k-mer node implementation optimized for single-dimensional storage arrays
 * with minimal memory overhead and cache-efficient access patterns.
 * This specialized node type stores only a single integer value per k-mer.
 *
 * @author Brian Bushnell
 * @date Oct 22, 2013
 */
public class KmerNodeU1D extends KmerNodeU {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public KmerNodeU1D(long[] pivot_){
		super(pivot_);
	}
	
	public KmerNodeU1D(long[] pivot_, int value_){
		super(pivot_);
		value=value_;
	}
	
	@Override
	public final KmerNodeU makeNode(long[] pivot_, int value_){
		return new KmerNodeU1D(pivot_, value_);
	}
	
	@Override
	public final KmerNodeU makeNode(long[] pivot_, int[] values_){
		throw new RuntimeException("Unimplemented");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Sets multiple values for a k-mer.
	 * This method is not implemented for single-dimensional nodes.
	 *
	 * @param kmer The k-mer to set values for
	 * @param vals Array of values to set
	 * @return Never returns - throws RuntimeException
	 * @throws RuntimeException Always thrown as this operation is not supported
	 */
	@Override
	public final int set(long[] kmer, int[] vals) {
		throw new RuntimeException("Unimplemented.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Gets the integer value stored in this node */
	@Override
	public int value(){return value;}
	
	/**
	 * Fills a singleton array with this node's value.
	 * @param singleton Single-element array to fill with the node's value
	 * @return The same array passed as parameter, now containing the node's value
	 */
	@Override
	protected int[] values(int[] singleton){
		assert(singleton.length==1);
		singleton[0]=value;
		return singleton;
	}
	
	/**
	 * Sets the value for this node.
	 * @param value_ The new value to set
	 * @return The newly set value
	 */
	@Override
	public int set(int value_){return value=value_;}
	
	/**
	 * Sets multiple values for this node.
	 * This method is not implemented for single-dimensional nodes.
	 *
	 * @param values_ Array of values to set
	 * @return Never returns - throws RuntimeException
	 * @throws RuntimeException Always thrown as this operation is not supported
	 */
	@Override
	protected int set(int[] values_){
		throw new RuntimeException("Unimplemented");
	}
	
	/** Gets the number of values stored (0 if value < 1, otherwise 1) */
	@Override
	int numValues(){return value<1 ? 0 : 1;}
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Indicates whether this node type supports resizing (always false) */
	@Override
	boolean canResize() {
		return false;
	}
	
	/** Indicates whether this node type supports rebalancing (always true) */
	@Override
	public boolean canRebalance() {
		return true;
	}
	
	/**
	 * Gets the array length for this node type.
	 * This method is deprecated and not supported.
	 *
	 * @return Never returns - throws RuntimeException
	 * @throws RuntimeException Always thrown as this operation is not supported
	 * @deprecated This method is not supported for this node type
	 */
	@Deprecated
	@Override
	public int arrayLength() {
		throw new RuntimeException("Unsupported.");
	}
	
	/**
	 * Resizes the node's internal storage.
	 * This method is deprecated and not supported.
	 * @throws RuntimeException Always thrown as this operation is not supported
	 * @deprecated This method is not supported for this node type
	 */
	@Deprecated
	@Override
	void resize() {
		throw new RuntimeException("Unsupported.");
	}
	
	/**
	 * Rebalances the node structure.
	 * This method is deprecated - use the ArrayList version instead.
	 * @throws RuntimeException Always thrown with instruction to use ArrayList version
	 * @deprecated Use rebalance(ArrayList<KmerNode>) instead
	 */
	@Deprecated
	@Override
	public void rebalance() {
		throw new RuntimeException("Please call rebalance(ArrayList<KmerNode>) instead, with an empty list.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Dumps k-mers as bytes to a ByteStreamWriter with filtering by count thresholds.
	 * Recursively processes left and right child nodes in tree traversal order.
	 *
	 * @param bsw The ByteStreamWriter to write k-mer data to
	 * @param k The k-mer length for output formatting
	 * @param mincount Minimum count threshold for k-mer inclusion
	 * @param maxcount Maximum count threshold for k-mer inclusion
	 * @param remaining Atomic counter for limiting total k-mers dumped
	 * @return Always returns true
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
	 * Multi-threaded version of k-mer dumping using ByteBuilder for buffering.
	 * Uses synchronized block for thread-safe ByteStreamWriter access and automatic
	 * buffer flushing when buffer exceeds 16000 bytes.
	 *
	 * @param bsw The ByteStreamWriter for final output
	 * @param bb The ByteBuilder buffer for accumulating k-mer data
	 * @param k The k-mer length for output formatting
	 * @param mincount Minimum count threshold for k-mer inclusion
	 * @param maxcount Maximum count threshold for k-mer inclusion
	 * @param remaining Atomic counter for limiting total k-mers dumped
	 * @return Always returns true
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
	 * Dumps k-mers as text to a StringBuilder with count filtering.
	 * Creates StringBuilder if null is passed and recursively processes child nodes.
	 *
	 * @param sb The StringBuilder to append k-mer text to (created if null)
	 * @param k The k-mer length for text conversion
	 * @param mincount Minimum count threshold for k-mer inclusion
	 * @param maxcount Maximum count threshold for k-mer inclusion
	 * @return The StringBuilder containing k-mer text data
	 */
	@Override
	protected final StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount){
		if(value<1){return sb;}
		if(sb==null){sb=new StringBuilder(32);}
		if(value>=mincount){sb.append(AbstractKmerTableU.toText(pivot, value, k)).append('\n');}
		if(left!=null){left.dumpKmersAsText(sb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(sb, k, mincount, maxcount);}
		return sb;
	}
	
	/**
	 * Dumps k-mers as text to a ByteBuilder with count filtering.
	 * Creates ByteBuilder if null is passed and recursively processes child nodes.
	 *
	 * @param bb The ByteBuilder to append k-mer data to (created if null)
	 * @param k The k-mer length for text conversion
	 * @param mincount Minimum count threshold for k-mer inclusion
	 * @param maxcount Maximum count threshold for k-mer inclusion
	 * @return The ByteBuilder containing k-mer text data
	 */
	@Override
	protected final ByteBuilder dumpKmersAsText(ByteBuilder bb, int k, int mincount, int maxcount){
		if(value<1){return bb;}
		if(bb==null){bb=new ByteBuilder(32);}
		if(value>=mincount){bb.append(AbstractKmerTableU.toBytes(pivot, value, k)).append('\n');}
		if(left!=null){left.dumpKmersAsText(bb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(bb, k, mincount, maxcount);}
		return bb;
	}
	
	/** Indicates this is not a two-dimensional node (always returns false) */
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
