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
 * Uses dynamic array resizing and bounded duplicate detection: single-value
 * insertion checks only the last four logical entries, not the entire list.
 * Array-taking construction retains the supplied array and logical length.
 * Callers must coordinate mutation and must not modify borrowed storage independently.
 *
 * @author Brian Bushnell
 * @date Nov 7, 2014
 */
public class KmerNode2D extends KmerNode{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a detached node with no value array and logical length zero. */
	public KmerNode2D(long pivot_){
		super(pivot_);
	}

	/** Creates a detached node with value storage and a trailing -1 sentinel.
	 * A nonnegative value occupies one logical slot; -1 leaves the logical prefix empty. */
	public KmerNode2D(long pivot_, int value_){
		super(pivot_);
		assert(value_>=0 || value_==-1);
		values=new int[] {value_, -1};
		//KMR018: -1 is padding, not a logical value; later insertion must start at slot zero.
		numValues=(value_>=0 ? 1 : 0);
	}

	/** Retains the value array; vlen must be its nonnegative logical entry count. */
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
	 * @param values_ Value storage retained directly by the new node
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
	 * Merges values into the matching node or creates a node that retains the array.
	 * Empty-pivot initialization also reports a newly created entry. Existing value
	 * storage uses the bounded duplicate check and negative terminator described below.
	 * @param kmer The k-mer to store values for
	 * @param vals Values to merge, or storage to retain when a node has no array
	 * @param vlen Logical length when retaining the array; upper bound when merging its nonnegative prefix
	 * @return 1 for a new node or initialized empty pivot, otherwise 0
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

	/** Returns the first array element, independent of the cached logical length.
	 * @return Zero for null storage, otherwise element zero (which must exist) */
	@Override
	protected int value(){return values==null ? 0 : values[0];}

	/**
	 * Borrows the backing array, including any spare capacity and sentinel entries.
	 * @param singleton Unused parameter for compatibility
	 * @return Live backing storage, or null when unallocated; do not mutate independently
	 */
	@Override
	protected int[] values(int[] singleton){
		return values;
	}

	/**
	 * Attempts to add one nonnegative value using bounded duplicate detection.
	 * @param value_ The value to add
	 * @return The supplied value, including when it was already present in the checked suffix
	 */
	@Override
	public int set(int value_){
		insertValue(value_);
		return value_;
	}

	/**
	 * Adopts an array when unallocated, otherwise merges up to vlen entries or the first negative value.
	 * @param values_ Array of values to store
	 * @param vlen Logical length when adopting storage, or maximum number of entries to merge
	 * @return 1 if values were null before, 0 otherwise
	 */
	@Override
	protected int set(int[] values_, int vlen){
		int ret=(values==null ? 1 : 0);
		insertValue(values_, vlen);
		return ret;
	}

	/** Returns the cached logical length without rescanning storage.
	 * @return Stored length; scalar -1 initialization contributes no logical value */
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
	 * @return 0 for a duplicate in the checked suffix, otherwise 1
	 */
	private int insertValue(int v){
		return insertIntoList(v);
	}

	/**
	 * Retains the supplied array when unallocated; otherwise merges its nonnegative logical prefix.
	 * The merge applies bounded duplicate detection and stops at vlen or the first negative entry.
	 * @param vals Array to retain or merge; null is allowed only with zero vlen
	 * @param vlen Logical length when retaining the array, or maximum number of entries to merge
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
		//KMR016: backing storage beyond the supplied logical prefix must not become stored values.
		for(int i=0; i<vlen; i++){
			final int v=vals[i];
			if(v<0){break;}
			insertIntoList(v);
		}
		return 0;
	}

	/** Counts a nonnegative prefix for the retained, currently disabled diagnostic assertions. */
	private final int countValues(int[] vals){
		if(vals==null){return 0;}
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

	/** Adds one value after checking at most slowAddLimit trailing logical entries.
	 * Grows full storage and fills newly allocated capacity with -1 sentinels. */
	private final int insertIntoList(final int v){
//		assert(countValues(values)==numValues) : countValues(values)+", "+numValues; //TODO: Slow assertion //123
		assert(v>=0);

		if(values==null){
			values=new int[] {v, -1};
			numValues=1;
			return 1;
		}

		//KMR017 policy question (previously "KmerNode2D dedup, likely by-design"): only the last four logical
		//entries are checked, avoiding a full scan but permitting older values to repeat. HashArray2D scans
		//its full list instead. Caller tolerance of this difference remains unverified; preserve it for now.
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

	/** Rejects table-level resizing; internal value-list growth remains supported.
	 * @return Always false */
	@Override
	boolean canResize(){return false;}

	/** Indicates whether this node type supports rebalancing operations.
	 * @return true - KmerNode2D supports rebalancing */
	@Override
	public boolean canRebalance(){return true;}

	/**
	 * Deprecated method that throws RuntimeException.
	 * @return Never returns - always throws exception
	 * @throws RuntimeException Always thrown as method is unsupported
	 */
	@Deprecated
	@Override
	public int arrayLength(){
		throw new RuntimeException("Unsupported.");
	}

	/** Deprecated method that throws RuntimeException.
	 * @throws RuntimeException Always thrown as method is unsupported */
	@Deprecated
	@Override
	void resize(){
		throw new RuntimeException("Unsupported.");
	}

	/** Deprecated method that throws RuntimeException.
	 * @throws RuntimeException Always thrown - use rebalance(ArrayList) instead */
	@Deprecated
	@Override
	public void rebalance(){
		throw new RuntimeException("Please call rebalance(ArrayList<KmerNode>) instead, with an empty list.");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Outputs k-mer and values as bytes to a stream writer, recursively traversing subtrees.
	 * Visits current, left and right nodes. Null storage or exhausted allowance suppresses
	 * the current subtree. Formatting consumes the value array, not numValues.
	 * @param bsw Stream writer for output
	 * @param k Length of k-mer for formatting
	 * @param mincount Minimum count threshold (unused in this implementation)
	 * @param maxcount Maximum count threshold (unused in this implementation)
	 * @param remaining Optional shared allowance, decremented before each node emission attempt;
	 * a negative result suppresses the current subtree
	 * @return Always true; does not report asynchronous writer completion or errors
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
	 * Visits current, left and right nodes with caller-owned mutable output storage.
	 * At 16000 bytes, copies the buffer into a writer job and clears it; only job
	 * submission is synchronized. Caller submits the final partial buffer and keeps
	 * the tree and buffer stable against concurrent mutation.
	 * @param bsw Already-started writer whose lifecycle remains with the caller
	 * @param bb Exclusive accumulation buffer
	 * @param k Length of k-mer for formatting
	 * @param mincount Minimum count threshold (unused)
	 * @param maxcount Maximum count threshold (unused)
	 * @param remaining Optional shared allowance, decremented before each node emission attempt;
	 * a negative result suppresses the current subtree
	 * @return Always true; does not report asynchronous writer completion or errors
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
	 * Visits current, left and right nodes; null value storage suppresses its subtree.
	 * @param sb Accumulation buffer, allocated if null and this node has value storage
	 * @param k Length of k-mer for formatting
	 * @param mincount Minimum count threshold (unused)
	 * @param maxcount Maximum count threshold (unused)
	 * @return Accumulation buffer, possibly null when this node has no value storage
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
	 * Visits current, left and right nodes; null value storage suppresses its subtree.
	 * @param bb Accumulation buffer, allocated if null and this node has value storage
	 * @param k Length of k-mer for formatting
	 * @param mincount Minimum count threshold (unused)
	 * @param maxcount Maximum count threshold (unused)
	 * @return Accumulation buffer, possibly null when this node has no value storage
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
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Live retained or grown storage; spare capacity normally contains -1 sentinels. */
	int[] values;
	/** Cached logical length, independent of the backing array's allocated length. */
	private int numValues;
	/** Maximum suffix length searched for duplicates on each insertion. */
	private static final int slowAddLimit=4;

}
