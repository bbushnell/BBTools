package kmer;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * A k-mer tree node that supports multiple values per k-mer.
 * Extends KmerNode to allow storage of arrays of integers rather than single values.
 * Uses dynamic array resizing and duplicate detection for efficient multi-value storage.
 *
 * @author Brian Bushnell
 * @date Nov 7, 2014
 */
public class KmerNode2D extends KmerNode {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public KmerNode2D(long pivot_){
		super(pivot_);
	}
	
	public KmerNode2D(long pivot_, int value_){
		super(pivot_);
		assert(value_>=0 || value_==-1);
		values=new int[] {value_, -1};
		numValues=1;
	}
	
	public KmerNode2D(long pivot_, int[] vals_, int vlen){
		super(pivot_);
		values=vals_;
		numValues=vlen;
		assert(values!=null || vlen==0);
		assert(values==null || (vlen<=values.length && vlen>=0));
//		assert(countValues(values)==vlen) : countValues(values)+", "+vlen; //TODO: Slow assertion //123
	}
	
	/**
	 * Factory method to create new KmerNode2D with single value.
	 * @param pivot_ The k-mer value for the new node
	 * @param value_ Initial value for the new node
	 * @return New KmerNode2D instance
	 */
	@Override
	public final KmerNode makeNode(long pivot_, int value_){
		return new KmerNode2D(pivot_, value_);
	}
	
	/**
	 * Factory method to create new KmerNode2D with multiple values.
	 * @param pivot_ The k-mer value for the new node
	 * @param values_ Array of values for the new node
	 * @param vlen Number of valid values in the array
	 * @return New KmerNode2D instance
	 */
	@Override
	public final KmerNode makeNode(long pivot_, int[] values_, int vlen){
		return new KmerNode2D(pivot_, values_, vlen);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
//	public final int set_Test(final long kmer, final int v[]){
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
	
	/**
	 * Sets multiple values for a k-mer in the tree, creating nodes as needed.
	 * Recursively traverses the tree to find the correct position.
	 * @param kmer The k-mer to store values for
	 * @param vals Array of values to store
	 * @param vlen Number of valid values in the array
	 * @return Number of new nodes added to the tree
	 */
	@Override
	public int set(long kmer, int vals[], int vlen){
		if(pivot<0){pivot=kmer; insertValue(vals, vlen); return 1;} //Allows initializing empty nodes to -1
		if(kmer<pivot){
			if(left==null){left=new KmerNode2D(kmer, vals, vlen); return 1;}
			return left.set(kmer, vals, vlen);
		}else if(kmer>pivot){
			if(right==null){right=new KmerNode2D(kmer, vals, vlen); return 1;}
			return right.set(kmer, vals, vlen);
		}else{
			insertValue(vals, vlen);
		}
		return 0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns the first value stored in this node.
	 * @return First value in the values array, or 0 if no values stored */
	@Override
	protected int value(){return values==null ? 0 : values[0];}
	
	/**
	 * Returns the values array for this node.
	 * @param singleton Unused parameter for compatibility
	 * @return The values array stored in this node
	 */
	@Override
	protected int[] values(int[] singleton){
		return values;
	}
	
	/**
	 * Adds a single value to this node.
	 * @param value_ The value to add
	 * @return The value that was added
	 */
	@Override
	public int set(int value_){
		insertValue(value_);
		return value_;
	}
	
	/**
	 * Sets multiple values for this node.
	 * @param values_ Array of values to store
	 * @param vlen Number of valid values in the array
	 * @return 1 if values were null before, 0 otherwise
	 */
	@Override
	protected int set(int[] values_, int vlen){
		int ret=(values==null ? 1 : 0);
		insertValue(values_, vlen);
		return ret;
	}
	
	/** Returns the number of values stored in this node.
	 * @return Count of valid values in the values array */
	@Override
	int numValues(){
//		assert(countValues(values)==numValues) : countValues(values)+", "+numValues; //TODO: Slow assertion //123
		return numValues;
//		asdf
//		if(values==null){return 0;}
//		for(int i=0; i<values.length; i++){
//			if(values[i]==-1){return i;}
//		}
//		return values.length;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Inserts a single value into this node's value list.
	 * @param v The value to insert
	 * @return Number of values added (0 if duplicate, 1 if new)
	 */
	private int insertValue(int v){
		return insertIntoList(v);
	}
	
	/**
	 * Inserts multiple values into this node's value list.
	 * Skips duplicate values and stops at negative values.
	 * @param vals Array of values to insert
	 * @param vlen Number of valid values to process
	 * @return 1 if values array was null before, 0 otherwise
	 */
	private int insertValue(int[] vals, int vlen){
//		assert(countValues(vals)==vlen) : countValues(vals)+", "+vlen; //TODO: Slow assertion //123
		assert(vals!=null || vlen==0);
		assert(vals==null || (vlen<=vals.length && vlen>=0));
		if(values==null){
			values=vals;
			numValues=vlen;
			return 1;
		}
		for(int v : vals){
			if(v<0){break;}
			insertIntoList(v);
		}
		return 0;
	}
	
	private final int countValues(int[] vals){
		if(vals==null) {return 0;}
		int count=0;
		for(int v : vals){
			if(v>=0){
				count++;
			}else{
				break;
			}
		}
		return count;
	}
	
	private final int insertIntoList(final int v){
//		assert(countValues(values)==numValues) : countValues(values)+", "+numValues; //TODO: Slow assertion //123
		assert(v>=0);
		
		if(values==null){
			values=new int[] {v, -1};
			numValues=1;
			return 1;
		}
		
		for(int i=numValues-1, lim=Tools.max(0, numValues-slowAddLimit); i>=lim; i--){//This is the slow bit
			if(values[i]==v){return 0;}
			if(values[i]<0){
				values[i]=v;
				numValues++;
				return 1;
			}
		}
		//At this point the size is big, or the element was not found
		
		if(numValues>=values.length){//resize
			assert(numValues==values.length);
			final int oldSize=values.length;
			final int newSize=(int)Tools.min(Shared.MAX_ARRAY_LEN, oldSize*2L);
			assert(newSize>values.length) : "Overflow.";
			values=KillSwitch.copyOf(values, newSize);
			Arrays.fill(values, oldSize, newSize, -1);
		}
		
		//quick add
		assert(values[numValues]<0);
		values[numValues]=v;
		numValues++;

//		assert(countValues(values)==numValues) : countValues(values)+", "+numValues; //TODO: Slow assertion //123
		return 1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Indicates whether this node type supports resizing operations.
	 * @return false - KmerNode2D does not support resizing */
	@Override
	boolean canResize() {
		return false;
	}
	
	/** Indicates whether this node type supports rebalancing operations.
	 * @return true - KmerNode2D supports rebalancing */
	@Override
	public boolean canRebalance() {
		return true;
	}

	/**
	 * Deprecated method that throws RuntimeException.
	 * @return Never returns - always throws exception
	 * @throws RuntimeException Always thrown as method is unsupported
	 */
	@Deprecated
	@Override
	public int arrayLength() {
		throw new RuntimeException("Unsupported.");
	}

	/** Deprecated method that throws RuntimeException.
	 * @throws RuntimeException Always thrown as method is unsupported */
	@Deprecated
	@Override
	void resize() {
		throw new RuntimeException("Unsupported.");
	}

	/** Deprecated method that throws RuntimeException.
	 * @throws RuntimeException Always thrown - use rebalance(ArrayList) instead */
	@Deprecated
	@Override
	public void rebalance() {
		throw new RuntimeException("Please call rebalance(ArrayList<KmerNode>) instead, with an empty list.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Outputs k-mer and values as bytes to a stream writer, recursively traversing subtrees.
	 * @param bsw Stream writer for output
	 * @param k Length of k-mer for formatting
	 * @param mincount Minimum count threshold (unused in this implementation)
	 * @param maxcount Maximum count threshold (unused in this implementation)
	 * @param remaining Counter for remaining items to process
	 * @return true if operation completed successfully
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
	 * Multi-threaded version of k-mer dumping using ByteBuilder for buffering.
	 * Flushes buffer to writer when it exceeds 16KB to maintain memory efficiency.
	 * @param bsw Thread-safe stream writer for output
	 * @param bb ByteBuilder for accumulating output before writing
	 * @param k Length of k-mer for formatting
	 * @param mincount Minimum count threshold (unused)
	 * @param maxcount Maximum count threshold (unused)
	 * @param remaining Counter for remaining items to process
	 * @return true if operation completed successfully
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
	 * Outputs k-mer and values as text to StringBuilder, recursively traversing subtrees.
	 * @param sb StringBuilder to append output to (created if null)
	 * @param k Length of k-mer for formatting
	 * @param mincount Minimum count threshold (unused)
	 * @param maxcount Maximum count threshold (unused)
	 * @return StringBuilder containing formatted output
	 */
	@Override
	protected final StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount){
		if(values==null){return sb;}
		if(sb==null){sb=new StringBuilder(32);}
		sb.append(AbstractKmerTable.toText(pivot, values, k)).append('\n');
		if(left!=null){left.dumpKmersAsText(sb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(sb, k, mincount, maxcount);}
		return sb;
	}
	
	/**
	 * Outputs k-mer and values as text to ByteBuilder, recursively traversing subtrees.
	 * @param bb ByteBuilder to append output to (created if null)
	 * @param k Length of k-mer for formatting
	 * @param mincount Minimum count threshold (unused)
	 * @param maxcount Maximum count threshold (unused)
	 * @return ByteBuilder containing formatted output
	 */
	@Override
	protected final ByteBuilder dumpKmersAsText(ByteBuilder bb, int k, int mincount, int maxcount){
		if(values==null){return bb;}
		if(bb==null){bb=new ByteBuilder(32);}
		bb.append(AbstractKmerTable.toBytes(pivot, values, k)).append('\n');
		if(left!=null){left.dumpKmersAsText(bb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(bb, k, mincount, maxcount);}
		return bb;
	}
	
	/** Indicates this node type supports two-dimensional (multi-value) storage.
	 * @return true - this is a 2D node type */
	@Override
	final boolean TWOD(){return true;}
	
	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	int[] values;
	private int numValues;
	private static final int slowAddLimit=4;
	
}
