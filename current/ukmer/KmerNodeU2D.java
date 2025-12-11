package ukmer;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Binary tree node for k-mer storage supporting multiple values per k-mer.
 * Extends KmerNodeU to provide two-dimensional value storage where each k-mer
 * can be associated with an array of integer values rather than a single value.
 * Supports dynamic array resizing and duplicate value prevention.
 *
 * @author Brian Bushnell
 * @date Nov 7, 2014
 */
public class KmerNodeU2D extends KmerNodeU {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public KmerNodeU2D(long[] pivot_){
		super(pivot_);
	}
	
	public KmerNodeU2D(long[] pivot_, int value_){
		super(pivot_);
		assert(value_>=0 || value_==-1);
		values=new int[] {value_, -1};
	}
	
	public KmerNodeU2D(long[] pivot_, int[] vals_){
		super(pivot_);
		values=vals_;
	}
	
	@Override
	public final KmerNodeU makeNode(long[] pivot_, int value_){
		return new KmerNodeU2D(pivot_, value_);
	}
	
	@Override
	public final KmerNodeU makeNode(long[] pivot_, int[] values_){
		return new KmerNodeU2D(pivot_, values_);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
//	public final int set_Test(final long[] kmer, final int v[]){
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
//			assert(v[0]==old || !contains(kmer, old));
//		}
//		return x;
//	}
	
	@Override
	public int set(long[] kmer, int vals[]){
		final int cmp=compare(kmer, pivot);
		if(cmp<0){
			if(left==null){left=new KmerNodeU2D(kmer, vals); return 1;}
			return left.set(kmer, vals);
		}else if(cmp>0){
			if(right==null){right=new KmerNodeU2D(kmer, vals); return 1;}
			return right.set(kmer, vals);
		}else{
			insertValue(vals);
		}
		return 0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Gets the first value from the values array.
	 * @return First value in array, or 0 if values array is null */
	@Override
	protected int value(){return values==null ? 0 : values[0];}
	
	/**
	 * Gets the complete values array for this node.
	 * @param singleton Unused parameter for compatibility with parent class
	 * @return The values array associated with this k-mer
	 */
	@Override
	protected int[] values(int[] singleton){
		return values;
	}
	
	/**
	 * Inserts a single value into this node's values array.
	 * @param value_ Value to insert
	 * @return The inserted value
	 */
	@Override
	public int set(int value_){
		insertValue(value_);
		return value_;
	}
	
	/**
	 * Sets the values array for this node.
	 * @param values_ Array of values to set
	 * @return 1 if values was null before setting, 0 otherwise
	 */
	@Override
	protected int set(int[] values_){
		int ret=(values==null ? 1 : 0);
		insertValue(values_);
		return ret;
	}
	
	/**
	 * Counts the number of valid values in the values array.
	 * Values are valid until the first -1 terminator is encountered.
	 * @return Number of valid values (non-negative values before first -1)
	 */
	@Override
	int numValues(){
		if(values==null){return 0;}
		for(int i=0; i<values.length; i++){
			if(values[i]==-1){return i;}
		}
		return values.length;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Inserts a value into the values array, preventing duplicates.
	 * If array is full, doubles its size. Uses -1 as array terminator.
	 * @param v Value to insert
	 * @return 1 if value was added, 0 if duplicate found
	 */
	private int insertValue(int v){
		for(int i=0; i<values.length; i++){
			if(values[i]==v){return 0;}
			if(values[i]==-1){values[i]=v;return 1;}
		}
		final int oldSize=values.length;
		final int newSize=(int)Tools.min(Shared.MAX_ARRAY_LEN, oldSize*2L);
		assert(newSize>values.length) : "Overflow.";
		values=Arrays.copyOf(values, newSize);
		values[oldSize]=v;
		Arrays.fill(values, oldSize+1, newSize, -1);
		return 1;
	}
	
	/**
	 * Inserts multiple values from an array into this node's values array.
	 * If values array is null, assigns the input array directly.
	 * Otherwise, inserts each value individually until -1 terminator.
	 *
	 * @param vals Array of values to insert
	 * @return 1 if values array was null and assigned, 0 otherwise
	 */
	private int insertValue(int[] vals){
		if(values==null){
			values=vals;
			return 1;
		}
		for(int v : vals){
			if(v<0){break;}
			insertValue(v);
		}
		return 0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Indicates whether this node type supports resizing operations.
	 * @return false, as KmerNodeU2D does not support resizing */
	@Override
	boolean canResize() {
		return false;
	}
	
	/** Indicates whether this node type supports rebalancing operations.
	 * @return true, as KmerNodeU2D supports tree rebalancing */
	@Override
	public boolean canRebalance() {
		return true;
	}

	/**
	 * Gets the array length for this node type.
	 * This method is deprecated and not supported.
	 * @throws RuntimeException Always thrown as operation is unsupported
	 * @deprecated This operation is not supported for KmerNodeU2D
	 */
	@Deprecated
	@Override
	public int arrayLength() {
		throw new RuntimeException("Unsupported.");
	}

	/**
	 * Resizes the internal data structures of this node.
	 * This method is deprecated and not supported.
	 * @throws RuntimeException Always thrown as operation is unsupported
	 * @deprecated This operation is not supported for KmerNodeU2D
	 */
	@Deprecated
	@Override
	void resize() {
		throw new RuntimeException("Unsupported.");
	}

	/**
	 * Rebalances the tree structure for this node.
	 * This method is deprecated; use rebalance(ArrayList) instead.
	 * @throws RuntimeException Always thrown with instruction to use alternative method
	 * @deprecated Use rebalance(ArrayList<KmerNode>) with empty list instead
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
	 * Dumps k-mers and their values as bytes to output stream in single-threaded mode.
	 * Performs in-order traversal of the tree structure.
	 *
	 * @param bsw Output stream writer for k-mer data
	 * @param k Length of k-mers being dumped
	 * @param mincount Minimum count threshold (unused in this implementation)
	 * @param maxcount Maximum count threshold (unused in this implementation)
	 * @param remaining Counter for remaining k-mers to process (may be null)
	 * @return true when dumping completes successfully
	 */
	@Override
	public final boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		if(values==null){return true;}
		if(remaining!=null && remaining.decrementAndGet()<0){return true;}
		bsw.printlnKmer(pivot, values, k);
		if(left!=null){left.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);}
		if(right!=null){right.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);}
		return true;
	}
	
	/**
	 * Dumps k-mers and their values as bytes to output stream in multi-threaded mode.
	 * Uses ByteBuilder for buffering and submits jobs to writer when buffer reaches 16KB.
	 *
	 * @param bsw Thread-safe output stream writer
	 * @param bb Byte buffer for accumulating output before submission
	 * @param k Length of k-mers being dumped
	 * @param mincount Minimum count threshold (unused in this implementation)
	 * @param maxcount Maximum count threshold (unused in this implementation)
	 * @param remaining Counter for remaining k-mers to process (may be null)
	 * @return true when dumping completes successfully
	 */
	@Override
	public final boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining){
		if(values==null){return true;}
		if(remaining!=null && remaining.decrementAndGet()<0){return true;}
		toBytes(pivot, values, k, bb);
		bb.nl();
		if(bb.length()>=16000){
			ByteBuilder bb2=new ByteBuilder(bb);
			synchronized(bsw){bsw.addJob(bb2);}
			bb.clear();
		}
		if(left!=null){left.dumpKmersAsBytes_MT(bsw, bb, k, mincount, maxcount, remaining);}
		if(right!=null){right.dumpKmersAsBytes_MT(bsw, bb, k, mincount, maxcount, remaining);}
		return true;
	}
	
	/**
	 * Dumps k-mers and their values as text using StringBuilder.
	 * Performs in-order traversal and appends formatted k-mer data.
	 *
	 * @param sb StringBuilder to append output to (created if null)
	 * @param k Length of k-mers being dumped
	 * @param mincount Minimum count threshold (unused in this implementation)
	 * @param maxcount Maximum count threshold (unused in this implementation)
	 * @return StringBuilder containing the formatted k-mer output
	 */
	@Override
	protected final StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount){
		if(values==null){return sb;}
		if(sb==null){sb=new StringBuilder(32);}
		sb.append(AbstractKmerTableU.toText(pivot, values, k)).append('\n');
		if(left!=null){left.dumpKmersAsText(sb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(sb, k, mincount, maxcount);}
		return sb;
	}
	
	/**
	 * Dumps k-mers and their values as text using ByteBuilder.
	 * Performs in-order traversal and appends formatted k-mer data.
	 *
	 * @param bb ByteBuilder to append output to (created if null)
	 * @param k Length of k-mers being dumped
	 * @param mincount Minimum count threshold (unused in this implementation)
	 * @param maxcount Maximum count threshold (unused in this implementation)
	 * @return ByteBuilder containing the formatted k-mer output
	 */
	@Override
	protected final ByteBuilder dumpKmersAsText(ByteBuilder bb, int k, int mincount, int maxcount){
		if(values==null){return bb;}
		if(bb==null){bb=new ByteBuilder(32);}
		bb.append(AbstractKmerTableU.toBytes(pivot, values, k)).append('\n');
		if(left!=null){left.dumpKmersAsText(bb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(bb, k, mincount, maxcount);}
		return bb;
	}
	
	/** Identifies this as a two-dimensional k-mer node type.
	 * @return true, indicating this node supports multiple values per k-mer */
	@Override
	final boolean TWOD(){return true;}
	
	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	int[] values;
	
}
