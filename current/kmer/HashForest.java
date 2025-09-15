package kmer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import shared.Primes;
import shared.Tools;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * @author Brian Bushnell
 * @date Oct 23, 2013
 *
 */
public final class HashForest extends AbstractKmerTable implements Iterable<KmerNode> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public HashForest(int initialSize, boolean autoResize_){
		this(initialSize, autoResize_, false);
	}
	
	public HashForest(int initialSize, boolean autoResize_, boolean twod_){
		if(initialSize>1){
			initialSize=(int)Tools.min(maxPrime, Primes.primeAtLeast(initialSize));
		}else{
			initialSize=1;
		}
		prime=initialSize;
		sizeLimit=(long) (initialSize*resizeMult);
		array=allocKmerNodeArray(prime);
		autoResize=autoResize_;
		TWOD=twod_;
	}
	
	private KmerNode makeNode(long kmer, int val){
		return (TWOD ? new KmerNode2D(kmer, val) : new KmerNode1D(kmer, val));
	}
	
	private KmerNode makeNode(long kmer, int[] vals, int vlen){
		assert(TWOD);
		return new KmerNode2D(kmer, vals, vlen);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public int increment(final long kmer, final int incr){
		final int cell=(int)(kmer%prime);
		KmerNode n=array[cell], prev=null;
		while(n!=null && n.pivot!=kmer){
			prev=n;
			n=(kmer<n.pivot ? n.left : n.right);
		}
		if(n==null){
			n=makeNode(kmer, incr);
			size++;
			if(prev==null){
				array[cell]=n;
			}else{
				if(kmer<prev.pivot){
					prev.left=n;
				}else{
					prev.right=n;
				}
			}
			if(autoResize && size>sizeLimit){resize();}
		}else{
			n.increment(kmer, incr);
		}
		return n.value();
	}
	
	@Override
	public int incrementAndReturnNumCreated(final long kmer, final int incr){
		final int cell=(int)(kmer%prime);
		KmerNode n=array[cell], prev=null;
		while(n!=null && n.pivot!=kmer){
			prev=n;
			n=(kmer<n.pivot ? n.left : n.right);
		}
		if(n==null){
			n=makeNode(kmer, incr);
			size++;
			if(prev==null){
				array[cell]=n;
			}else{
				if(kmer<prev.pivot){
					prev.left=n;
				}else{
					prev.right=n;
				}
			}
			if(autoResize && size>sizeLimit){resize();}
			return 1;
		}else{
			n.increment(kmer, incr);
			return 0;
		}
	}
	
//	public final int set_Test(final long kmer, final int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
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
//			assert(v==old || !contains(kmer, old));
//		}
//		return x;
//	}
//
//	public final int setIfNotPresent_Test(long kmer, int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
////			int[] vals=getValues(kmer, null);
////			assert(vals==null || contains(kmer, vals));
////			x=setIfNotPresent(kmer, v);
////			assert(contains(kmer, vals));
////			assert(contains(kmer, v));
//			x=0;
//			assert(false);
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=setIfNotPresent0(kmer, v);
//			assert((old<1 && contains(kmer, v)) || (old>0 && contains(kmer, old))) : kmer+", "+old+", "+v;
//		}
//		return x;
//	}
//
//	public final int set_Test(final long kmer, final int v[]){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
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
	public int set(long kmer, int value){
		int x=1, cell=(int)(kmer%prime);
		final KmerNode n=array[cell];
		if(n==null){
			array[cell]=makeNode(kmer, value);
		}else{
			x=n.set(kmer, value);
		}
		size+=x;
		if(autoResize && size>sizeLimit){resize();}
		return x;
	}
	
	@Override
	public int set(long kmer, int[] vals, int vlen) {
		int x=1, cell=(int)(kmer%prime);
		final KmerNode n=array[cell];
		if(n==null){
			array[cell]=makeNode(kmer, vals, vlen);
		}else{
			x=n.set(kmer, vals, vlen);
		}
		size+=x;
		if(autoResize && size>sizeLimit){resize();}
		return x;
	}
	
	@Override
	public int setIfNotPresent(long kmer, int value){
		int x=1, cell=(int)(kmer%prime);
		final KmerNode n=array[cell];
		if(n==null){
			array[cell]=makeNode(kmer, value);
		}else{
			x=n.setIfNotPresent(kmer, value);
		}
		size+=x;
		if(autoResize && size>sizeLimit){resize();}
		return x;
	}
	
	@Override
	public final int getValue(long kmer){
		int cell=(int)(kmer%prime);
		KmerNode n=array[cell];
		return n==null ? -1 : n.getValue(kmer);
	}
	
	@Override
	public int[] getValues(long kmer, int[] singleton){
		int cell=(int)(kmer%prime);
		KmerNode n=array[cell];
		return n==null ? null : n.getValues(kmer, singleton);
	}
	
	@Override
	public boolean contains(long kmer){
		return get(kmer)!=null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final void initializeOwnership(){
		for(KmerNode n : array){
			if(n!=null){n.initializeOwnership();}
		}
	}
	
	@Override
	public final void clearOwnership(){initializeOwnership();}
	
	@Override
	public final int setOwner(final long kmer, final int newOwner){
		final int cell=(int)(kmer%prime);
		KmerNode n=array[cell];
		assert(n!=null);
		return n.setOwner(kmer, newOwner);
	}
	
	@Override
	public final boolean clearOwner(final long kmer, final int owner){
		final int cell=(int)(kmer%prime);
		KmerNode n=array[cell];
		assert(n!=null);
		return n.clearOwner(kmer, owner);
	}
	
	@Override
	public final int getOwner(final long kmer){
		final int cell=(int)(kmer%prime);
		KmerNode n=array[cell];
		assert(n!=null);
		return n.getOwner(kmer);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	final KmerNode get(long kmer){
		int cell=(int)(kmer%prime);
		KmerNode n=array[cell];
		while(n!=null && n.pivot!=kmer){
			n=(kmer<n.pivot ? n.left : n.right);
		}
		return n;
	}
	
	public final KmerNode getNode(int cell){
		KmerNode n=array[cell];
		return n;
	}
	
	boolean insert(KmerNode n){
		n.left=null;
		n.right=null;
		int cell=(int)(n.pivot%prime);
		if(array[cell]==null){
			array[cell]=n;
			return true;
		}
		return array[cell].insert(n);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	boolean canResize() {return true;}
	
	@Override
	public boolean canRebalance() {return true;}
	
	@Override
	public long size() {return size;}
	
	@Override
	public int arrayLength() {return array.length;}
	
	@Override
	synchronized void resize(){
//		assert(false);
//		System.err.println("Resizing from "+prime+"; load="+(size*1f/prime));
		sizeLimit=Tools.max((long)(size*1.4), (long)(maxLoadFactor*prime));

		final long maxAllowedByLoadFactor=(long)(size*minLoadMult);
		final long minAllowedByLoadFactor=(long)(size*maxLoadMult);
		assert(maxAllowedByLoadFactor>=minAllowedByLoadFactor);
		if(maxAllowedByLoadFactor<prime){return;}
		
		long x=10+(long)(prime*resizeMult);
		x=Tools.max(x, minAllowedByLoadFactor);
		x=Tools.min(x, maxAllowedByLoadFactor);
		
		int prime2=(int)Tools.min(maxPrime, Primes.primeAtLeast(x));
		
		if(prime2<=prime){return;}
		
		prime=prime2;
//		System.err.println("Resized to "+prime+"; load="+(size*1f/prime));
		KmerNode[] old=array;
		array=allocKmerNodeArray(prime2);
		ArrayList<KmerNode> list=new ArrayList<KmerNode>(1000);
		for(int i=0; i<old.length; i++){
			if(old[i]!=null){
				old[i].traverseInfix(list);
				for(KmerNode n : list){insert(n);}
				list.clear();
			}
		}
		sizeLimit=Tools.max((long)(size*1.4), (long)(maxLoadFactor*prime));
	}
	
	@Override
	public void rebalance(){
		ArrayList<KmerNode> list=new ArrayList<KmerNode>(1000);
		for(int i=0; i<array.length; i++){
			if(array[i]!=null){array[i]=array[i].rebalance(list);}
		}
	}
	
	public void clear() {
		size=0;
		Arrays.fill(array, null);
	}
	
	@Override
	long regenerate(final int limit) {
		throw new RuntimeException("Not implemented.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
//		tsw.print("HashForest:\n");
		for(int i=0; i<array.length; i++){
			KmerNode node=array[i];
			if(node!=null && node.value()>=mincount){
//				StringBuilder sb=new StringBuilder();
//				tsw.print(node.dumpKmersAsText(sb, k, mincount, maxcount));
				node.dumpKmersAsText(tsw, k, mincount, maxcount);
			}
		}
		return true;
	}
	
	@Override
	public boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		for(int i=0; i<array.length; i++){
			KmerNode node=array[i];
			if(node!=null && node.value()>=mincount){
				if(remaining!=null && remaining.decrementAndGet()<0){return true;}
				node.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);
			}
		}
		return true;
	}
	
	@Override
	public boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining){
		for(int i=0; i<array.length; i++){
			KmerNode node=array[i];
			if(node!=null && node.value()>=mincount){
				if(remaining!=null && remaining.decrementAndGet()<0){return true;}
				node.dumpKmersAsBytes_MT(bsw, bb, k, mincount, maxcount, remaining);
			}
		}
		return true;
	}
	
	@Override
	public void fillHistogram(long[] ca, int max){
		for(int i=0; i<array.length; i++){
			KmerNode node=array[i];
			if(node!=null){
				node.fillHistogram(ca, max);
			}
		}
	}
	
	@Override
	public final void fillHistogram(SuperLongList sll){
		for(int i=0; i<array.length; i++){
			KmerNode node=array[i];
			if(node!=null){
				node.fillHistogram(sll);
			}
		}
	}
	
	@Override
	public void countGC(long[] gcCounts, int max){
		for(int i=0; i<array.length; i++){
			KmerNode node=array[i];
			if(node!=null){
				node.countGC(gcCounts, max);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Iteration           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public Iterator<KmerNode> iterator() {
		return toList().iterator();
	}
	
	public ArrayList<KmerNode> toList(){
		assert(size<Integer.MAX_VALUE);
		ArrayList<KmerNode> list=new ArrayList<KmerNode>((int)size);
		for(int i=0; i<array.length; i++){
			if(array[i]!=null){array[i].traverseInfix(list);}
		}
		assert(list.size()==size);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public KmerNode[] array() {return array;}
	
	KmerNode[] array;
	int prime;
	long size=0;
	long sizeLimit;
	final boolean autoResize;
	final boolean TWOD;
	private final Lock lock=new ReentrantLock();
	
	@Override
	final Lock getLock(){return lock;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	final static int maxPrime=(int)Primes.primeAtMost(Integer.MAX_VALUE);
	final static float resizeMult=2.5f; //Resize by a minimum of this much
	final static float minLoadFactor=0.75f; //Resize by enough to get the load above this factor
	final static float maxLoadFactor=2.5f; //Resize by enough to get the load under this factor
	final static float minLoadMult=1/minLoadFactor;
	final static float maxLoadMult=1/maxLoadFactor;
	

	
}
