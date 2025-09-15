package kmer;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Allows multiple values per kmer.
 * @author Brian Bushnell
 * @date Nov 7, 2014
 *
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
	
	@Override
	public final KmerNode makeNode(long pivot_, int value_){
		return new KmerNode2D(pivot_, value_);
	}
	
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
	
	/** Returns number of nodes added */
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
	protected int set(int[] values_, int vlen){
		int ret=(values==null ? 1 : 0);
		insertValue(values_, vlen);
		return ret;
	}
	
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
	
	/** Returns number of values added */
	private int insertValue(int v){
		return insertIntoList(v);
	}
	
	/** Returns number of values added */
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
		sb.append(AbstractKmerTable.toText(pivot, values, k)).append('\n');
		if(left!=null){left.dumpKmersAsText(sb, k, mincount, maxcount);}
		if(right!=null){right.dumpKmersAsText(sb, k, mincount, maxcount);}
		return sb;
	}
	
	@Override
	protected final ByteBuilder dumpKmersAsText(ByteBuilder bb, int k, int mincount, int maxcount){
		if(values==null){return bb;}
		if(bb==null){bb=new ByteBuilder(32);}
		bb.append(AbstractKmerTable.toBytes(pivot, values, k)).append('\n');
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
	private int numValues;
	private static final int slowAddLimit=4;
	
}
