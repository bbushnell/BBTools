package kmer;

import java.util.ArrayList;
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
public final class KmerTable extends AbstractKmerTable {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public KmerTable(int initialSize, boolean autoResize_){
		if(initialSize>1){
			initialSize=(int)Tools.min(maxPrime, Primes.primeAtLeast(initialSize));
		}else{
			initialSize=1;
		}
		prime=initialSize;
		sizeLimit=(long) (initialSize*resizeMult);
		array=new KmerLink[prime];
		autoResize=autoResize_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public int increment(final long kmer, final int incr){
		final int cell=(int)(kmer%prime);
		KmerLink n=array[cell], prev=null;
		while(n!=null && n.pivot!=kmer){
			prev=n;
			n=n.next;
		}
		if(n==null){
			n=new KmerLink(kmer, incr);
			size++;
			if(prev==null){
				array[cell]=n;
			}else{
				prev.next=n;
			}
			if(autoResize && size>sizeLimit){resize();}
		}else{
			n.value+=incr;
			if(n.value<0){n.value=Integer.MAX_VALUE;}
		}
		return n.value;
	}
	
	@Override
	public int incrementAndReturnNumCreated(final long kmer, final int incr){
		final int cell=(int)(kmer%prime);
		KmerLink n=array[cell], prev=null;
		while(n!=null && n.pivot!=kmer){
			prev=n;
			n=n.next;
		}
		if(n==null){
			n=new KmerLink(kmer, incr);
			size++;
			if(prev==null){
				array[cell]=n;
			}else{
				prev.next=n;
			}
			if(autoResize && size>sizeLimit){resize();}
			return 1;
		}else{
			n.value+=incr;
			if(n.value<0){n.value=Integer.MAX_VALUE;}
			return 0;
		}
	}
	
	@Override
	public int set(long kmer, int value){
		int x=1, cell=(int)(kmer%prime);
		final KmerLink n=array[cell];
		if(n==null){
			array[cell]=new KmerLink(kmer, value);
		}else{
			x=n.set(kmer, value);
		}
		size+=x;
		if(autoResize && size>sizeLimit){resize();}
		return x;
	}
	
	@Override
	public int set(long kmer, int[] vals, int vlen) {
		throw new RuntimeException("Unimplemented.");
	}
	
	@Override
	public int setIfNotPresent(long kmer, int value){
		int x=1, cell=(int)(kmer%prime);
		final KmerLink n=array[cell];
		if(n==null){
			array[cell]=new KmerLink(kmer, value);
		}else{
			x=n.setIfNotPresent(kmer, value);
		}
		size+=x;
		if(autoResize && size>sizeLimit){resize();}
		return x;
	}
	
	@Override
	public int getValue(long kmer){
		int cell=(int)(kmer%prime);
		KmerLink n=array[cell];
		while(n!=null && n.pivot!=kmer){n=n.next;}
		return n==null ? 0 : n.value;
	}
	
	@Override
	public int[] getValues(long kmer, int[] singleton){
		assert(array.length==0);
		singleton[0]=getValue(kmer);
		return singleton;
	}
	
	@Override
	public boolean contains(long kmer){
		KmerLink node=get(kmer);
		return node!=null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final void initializeOwnership(){
		for(KmerLink n : array){
			if(n!=null){n.initializeOwnership();}
		}
	}
	
	@Override
	public final void clearOwnership(){
		for(KmerLink n : array){
			if(n!=null){n.clearOwnership();}
		}
	}
	
	@Override
	public final int setOwner(final long kmer, final int newOwner){
		final int cell=(int)(kmer%prime);
		KmerLink n=array[cell];
		assert(n!=null);
		return n.setOwner(kmer, newOwner);
	}
	
	@Override
	public final boolean clearOwner(final long kmer, final int owner){
		final int cell=(int)(kmer%prime);
		KmerLink n=array[cell];
		assert(n!=null);
		return n.clearOwner(kmer, owner);
	}
	
	@Override
	public final int getOwner(final long kmer){
		final int cell=(int)(kmer%prime);
		KmerLink n=array[cell];
		assert(n!=null);
		return n.getOwner(kmer);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	KmerLink get(long kmer){
		int cell=(int)(kmer%prime);
		KmerLink n=array[cell];
		while(n!=null && n.pivot!=kmer){n=n.next;}
		return n;
	}
	
	boolean insert(KmerLink n){
		n.next=null;
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
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	boolean canResize() {return true;}
	
	@Override
	public boolean canRebalance() {return false;}
	
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
		KmerLink[] old=array;
		array=new KmerLink[prime2];
		ArrayList<KmerLink> list=new ArrayList<KmerLink>(1000);
		for(int i=0; i<old.length; i++){
			if(old[i]!=null){
				old[i].traverseInfix(list);
				for(KmerLink n : list){insert(n);}
				list.clear();
			}
		}
		sizeLimit=Tools.max((long)(size*1.4), (long)(maxLoadFactor*prime));
	}
	
	@Override
	public void rebalance(){
		ArrayList<KmerLink> list=new ArrayList<KmerLink>(1000);
		for(int i=0; i<array.length; i++){
			if(array[i]!=null){array[i]=array[i].rebalance(list);}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Deprecated
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
		throw new RuntimeException("TODO");
	}
	
	@Override
	public boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining){
		for(int i=0; i<array.length; i++){
			KmerLink node=array[i];
			if(node!=null && node.value>=mincount){
				node.dumpKmersAsBytes_MT(bsw, bb, k, mincount, maxcount, remaining);
			}
		}
		return true;
	}
	
	@Deprecated
	@Override
	public boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		throw new RuntimeException("TODO");
	}
	
	@Deprecated
	@Override
	public void fillHistogram(long[] ca, int max){
		throw new RuntimeException("TODO");
	}
	
	@Deprecated
	@Override
	public void fillHistogram(SuperLongList sll){
		throw new RuntimeException("TODO");
	}
	
	@Deprecated
	@Override
	public void countGC(long[] gcCounts, int max){
		throw new RuntimeException("TODO");
	}
	
	@Deprecated
	@Override
	public long regenerate(final int limit){
		throw new RuntimeException("TODO - remove zero-value links.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	KmerLink[] array;
	int prime;
	long size=0;
	long sizeLimit;
	final boolean autoResize;
	private final Lock lock=new ReentrantLock();
	
	@Override
	final Lock getLock(){return lock;}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	final static int maxPrime=(int)Primes.primeAtMost(Integer.MAX_VALUE);
	final static float resizeMult=2f; //Resize by a minimum of this much
	final static float minLoadFactor=0.5f; //Resize by enough to get the load above this factor
	final static float maxLoadFactor=0.98f; //Resize by enough to get the load under this factor
	final static float minLoadMult=1/minLoadFactor;
	final static float maxLoadMult=1/maxLoadFactor;
	
}
