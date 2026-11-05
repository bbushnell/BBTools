package ukmer;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * HashForestU association node that stores a borrowed live value list per
 * long-kmer pivot. Child creation keeps the 2D node family, scalar sets require
 * initialized nonempty storage and skip already-present values, and dump methods
 * emit self/left/right lists without count-threshold filtering.
 *
 * @author Brian Bushnell
 * @date Nov 7, 2014
 */
public class KmerNodeU2D extends KmerNodeU{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a multi-value node with no value list yet. */
	public KmerNodeU2D(long[] pivot_){
		super(pivot_);
	}

	/** Creates a multi-value node containing one value and a -1 terminator. */
	public KmerNodeU2D(long[] pivot_, int value_){
		super(pivot_);
		assert(value_>=0 || value_==-1);
		values=new int[] {value_, -1};
	}

	/** Creates a multi-value node that borrows the supplied live value list. */
	public KmerNodeU2D(long[] pivot_, int[] vals_){
		super(pivot_);
		values=vals_;
	}

	/** Creates a same-family child node for a single initial value. */
	@Override
	public final KmerNodeU makeNode(long[] pivot_, int value_){
		return new KmerNodeU2D(pivot_, value_);
	}

	/** Creates a same-family child node that borrows an existing value list. */
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

	/** Inserts a borrowed value list by current tree order, merging until negative on equality. */
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

	/**
	 * Returns the first stored value.
	 * @return first value in the list, or 0 when this node has no value list
	 */
	@Override
	protected int value(){return values==null ? 0 : values[0];}

	/**
	 * Returns the borrowed live value list backing this node; callers must not
	 * mutate adopted or returned arrays independently.
	 * @param singleton unused scratch array for interface compatibility
	 * @return the direct value array associated with this k-mer
	 */
	@Override
	protected int[] values(int[] singleton){
		return values;
	}

	/**
	 * Attempts to add one distinct value to initialized nonempty storage.
	 * @param value_ value to insert
	 * @return the supplied value, even when it was already present
	 */
	@Override
	public int set(int value_){
		insertValue(value_);
		return value_;
	}

	/**
	 * Directly adopts this node's first value list without deduplication, or merges
	 * values until the first negative entry.
	 * @param values_ array of values to insert
	 * @return 1 if this node did not yet have a value list, 0 otherwise
	 */
	@Override
	protected int set(int[] values_){
		int ret=(values==null ? 1 : 0);
		insertValue(values_);
		return ret;
	}

	/** Counts entries before -1, returns the full length when unterminated, or 0 for null. */
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
	 * Inserts one value if absent, requiring initialized nonempty storage and growing
	 * full lists that have no -1 terminator.
	 * @param v value to insert
	 * @return 1 if value was added, 0 if it was already present
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
	 * Directly adopts the first value list for this node without deduplication, or
	 * merges values until the first negative entry.
	 * @param vals array of values to insert
	 * @return 1 if the input list was adopted, 0 if it was merged
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
	 * Emits self/left/right text-byte value-list records using total kbig length
	 * {@code k}; a null value list prunes the whole current subtree, {@code true}
	 * may include quota early-stop and writer completion is caller-owned.
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
	 * Appends self/left/right text-byte value-list records to the caller-thread-
	 * exclusive buffer, submitting jobs when it reaches or exceeds 16000 bytes; a
	 * null value list prunes
	 * the whole current subtree, {@code true} may include quota early-stop, and
	 * final-tail submission plus writer completion are caller-owned.
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

	/** Appends self/left/right text records, pruning this subtree at null list roots. */
	@Override
	protected final StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount){
		if(values==null){return sb;}
		if(sb==null){sb=new StringBuilder(32);}
		sb.append(AbstractKmerTableU.toText(pivot, values, k)).append('\n');
		if(left!=null){left.dumpKmersAsText(sb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(sb, k, mincount, maxcount);}
		return sb;
	}

	/** Appends self/left/right text bytes, pruning this subtree at null list roots. */
	@Override
	protected final ByteBuilder dumpKmersAsText(ByteBuilder bb, int k, int mincount, int maxcount){
		if(values==null){return bb;}
		if(bb==null){bb=new ByteBuilder(32);}
		bb.append(AbstractKmerTableU.toBytes(pivot, values, k)).append('\n');
		if(left!=null){left.dumpKmersAsText(bb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(bb, k, mincount, maxcount);}
		return bb;
	}

	/** Identifies this as a two-dimensional node storing lists of values. */
	@Override
	final boolean TWOD(){return true;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrowed live values; may be full without a terminator. Do not mutate independently. */
	int[] values;

}
