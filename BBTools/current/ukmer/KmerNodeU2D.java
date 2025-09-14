package ukmer;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Allows multiple values per kmer.
 * @author Brian Bushnell
 * @date Nov 7, 2014
 *
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
	
	/** Returns number of nodes added */
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
	
	@Override
	protected int value(){return values==null ? 0 : values[0];}
	
	@Override
	protected int[] values(int[] singleton){
		return values;
	}
	
	@Override
	public int set(int value_){
		insertValue(value_);
		return value_;
	}
	
	@Override
	protected int set(int[] values_){
		int ret=(values==null ? 1 : 0);
		insertValue(values_);
		return ret;
	}
	
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
	
	/** Returns number of values added */
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
	
	/** Returns number of values added */
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
	
	@Override
	boolean canResize() {
		return false;
	}
	
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
	
	@Override
	public final boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		if(values==null){return true;}
		if(remaining!=null && remaining.decrementAndGet()<0){return true;}
		bsw.printlnKmer(pivot, values, k);
		if(left!=null){left.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);}
		if(right!=null){right.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);}
		return true;
	}
	
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
	
	@Override
	protected final StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount){
		if(values==null){return sb;}
		if(sb==null){sb=new StringBuilder(32);}
		sb.append(AbstractKmerTableU.toText(pivot, values, k)).append('\n');
		if(left!=null){left.dumpKmersAsText(sb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(sb, k, mincount, maxcount);}
		return sb;
	}
	
	@Override
	protected final ByteBuilder dumpKmersAsText(ByteBuilder bb, int k, int mincount, int maxcount){
		if(values==null){return bb;}
		if(bb==null){bb=new ByteBuilder(32);}
		bb.append(AbstractKmerTableU.toBytes(pivot, values, k)).append('\n');
		if(left!=null){left.dumpKmersAsText(bb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(bb, k, mincount, maxcount);}
		return bb;
	}
	
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
