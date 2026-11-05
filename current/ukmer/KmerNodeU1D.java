package ukmer;

import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/**
 * HashForestU count node that stores one scalar value for a long-kmer pivot.
 * Child creation keeps the 1D node family, list insertion is unsupported while
 * singleton retrieval works, and dump methods emit self/left/right from each
 * positive-count root when count thresholds pass.
 *
 * @author Brian Bushnell
 * @date Oct 22, 2013
 */
public class KmerNodeU1D extends KmerNodeU{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a count node with no positive count yet. */
	public KmerNodeU1D(long[] pivot_){
		super(pivot_);
	}

	/** Creates a count node with an initial value. */
	public KmerNodeU1D(long[] pivot_, int value_){
		super(pivot_);
		value=value_;
	}

	/** Creates a same-family child node for one value. */
	@Override
	public final KmerNodeU makeNode(long[] pivot_, int value_){
		return new KmerNodeU1D(pivot_, value_);
	}

	/** Rejects multi-value child creation for this one-value node type. */
	@Override
	public final KmerNodeU makeNode(long[] pivot_, int[] values_){
		throw new RuntimeException("Unimplemented");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Rejects value-list insertion; 1D forest nodes store one scalar count.
	 * @throws RuntimeException always, because this node family is scalar-only
	 */
	@Override
	public final int set(long[] kmer, int[] vals){
		throw new RuntimeException("Unimplemented.");
	}

	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the stored scalar count. */
	@Override
	public int value(){return value;}

	/**
	 * Returns the supplied singleton scratch after storing this node's count in slot zero.
	 * @param singleton one-element scratch array
	 * @return the supplied singleton array
	 */
	@Override
	protected int[] values(int[] singleton){
		assert(singleton.length==1);
		singleton[0]=value;
		return singleton;
	}

	/**
	 * Stores one scalar count.
	 * @param value_ new count
	 * @return the stored count
	 */
	@Override
	public int set(int value_){return value=value_;}

	/**
	 * Rejects value-list replacement; 1D nodes cannot hold association lists.
	 * @throws RuntimeException always, because this node family is scalar-only
	 */
	@Override
	protected int set(int[] values_){
		throw new RuntimeException("Unimplemented");
	}

	/** Gets the number of values stored (0 if {@code value<1}, otherwise 1). */
	@Override
	int numValues(){return value<1 ? 0 : 1;}

	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/

	/** This tree node cannot resize a backing array. */
	@Override
	boolean canResize(){
		return false;
	}

	/** This tree node can be rebalanced by materializing and rebuilding links. */
	@Override
	public boolean canRebalance(){
		return true;
	}

	/**
	 * Array length is unsupported for tree nodes.
	 * @throws RuntimeException always
	 * @deprecated this node has no backing hash array
	 */
	@Deprecated
	@Override
	public int arrayLength(){
		throw new RuntimeException("Unsupported.");
	}

	/**
	 * Resize is unsupported for tree nodes.
	 * @throws RuntimeException always
	 * @deprecated this node has no backing hash array
	 */
	@Deprecated
	@Override
	void resize(){
		throw new RuntimeException("Unsupported.");
	}

	/**
	 * Parameterless rebalance is unsupported; callers must pass reusable scratch storage.
	 * @throws RuntimeException always
	 * @deprecated use the package-private ArrayList overload instead
	 */
	@Deprecated
	@Override
	public void rebalance(){
		throw new RuntimeException("Please call rebalance(ArrayList<KmerNode>) instead, with an empty list.");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Emits self/left/right text-byte scalar records using total kbig length {@code k};
	 * a nonpositive value prunes the whole current subtree, {@code true} may include
	 * quota early-stop and writer completion is caller-owned.
	 */
	@Override
	public final boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		if(value<1){return true;}
		if(value>=mincount && value<=maxcount){
			if(remaining!=null && remaining.decrementAndGet()<0){return true;}
			bsw.printlnKmer(pivot, value, k);
		}
		if(left!=null){left.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);}
		if(right!=null){right.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);}
		return true;
	}

	/**
	 * Appends self/left/right text-byte scalar records to the caller-thread-exclusive
	 * buffer, submitting jobs when it reaches or exceeds 16000 bytes; a nonpositive
	 * value prunes the whole current subtree, {@code true} may include quota
	 * early-stop, and final-tail
	 * submission plus writer completion are caller-owned.
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

	/** Appends self/left/right text records, pruning this subtree at nonpositive roots. */
	@Override
	protected final StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount){
		if(value<1){return sb;}
		if(sb==null){sb=new StringBuilder(32);}
		if(value>=mincount && value<=maxcount){sb.append(AbstractKmerTableU.toText(pivot, value, k)).append('\n');}
		if(left!=null){left.dumpKmersAsText(sb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(sb, k, mincount, maxcount);}
		return sb;
	}

	/** Appends self/left/right text bytes, pruning this subtree at nonpositive roots. */
	@Override
	protected final ByteBuilder dumpKmersAsText(ByteBuilder bb, int k, int mincount, int maxcount){
		if(value<1){return bb;}
		if(bb==null){bb=new ByteBuilder(32);}
		if(value>=mincount && value<=maxcount){bb.append(AbstractKmerTableU.toBytes(pivot, value, k)).append('\n');}
		if(left!=null){left.dumpKmersAsText(bb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(bb, k, mincount, maxcount);}
		return bb;
	}

	/** Indicates this is not a two-dimensional node. */
	@Override
	final boolean TWOD(){return false;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Single stored count; values below one mark a logically empty node. */
	int value;

}
