package kmer;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;

import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import shared.Tools;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * @author Brian Bushnell
 * @date Oct 22, 2013
 *
 */
public class KmerLink extends AbstractKmerTable {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public KmerLink(long pivot_){
		pivot=pivot_;
	}
	
	public KmerLink(long pivot_, int value_){
		pivot=pivot_;
		value=value_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final int incrementAndReturnNumCreated(final long kmer, final int incr) {
		int x=increment(kmer, incr);
		return x==incr ? 1 : 0;
	}
	
	@Override
	public int increment(final long kmer, final int incr){
		if(pivot<0){pivot=kmer; return (value=incr);} //Allows initializing empty nodes to -1
		if(kmer==pivot){
			if(value<Integer.MAX_VALUE){value+=incr;}
			return value;
		}
		if(next==null){next=new KmerLink(kmer, incr); return 1;}
		return next.increment(kmer, incr);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Returns number of nodes added */
	@Override
	public int set(long kmer, int value_){
		if(pivot<0){pivot=kmer; value=value_; return 1;} //Allows initializing empty nodes to -1
		if(kmer==pivot){value=value_; return 0;}
		if(next==null){next=new KmerLink(kmer, value_); return 1;}
		return next.set(kmer, value_);
	}
	
	/** Returns number of nodes added */
	@Override
	public int setIfNotPresent(long kmer, int value_){
		if(pivot<0){pivot=kmer; value=value_; return 1;} //Allows initializing empty nodes to -1
		if(kmer==pivot){return 0;}
		if(next==null){next=new KmerLink(kmer, value_); return 1;}
		return next.setIfNotPresent(kmer, value_);
	}
	
	@Override
	KmerLink get(long kmer){
		if(kmer==pivot){return this;}
		return next==null ? null : next.get(kmer);
	}
	
	boolean insert(KmerLink n){
		assert(pivot!=-1);
		if(pivot==n.pivot){return false;}
		if(next==null){next=n; return true;}
		return next.insert(n);
	}
	
	@Override
	public boolean contains(long kmer){
		KmerLink node=get(kmer);
		return node!=null;
	}
	
	void traversePrefix(ArrayList<KmerLink> list){
		if(next!=null){next.traversePrefix(list);}
		list.add(this);
	}
	
	void traverseInfix(ArrayList<KmerLink> list){
		list.add(this);
		if(next!=null){next.traverseInfix(list);}
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
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final void initializeOwnership(){
		owner=-1;
		if(next!=null){next.initializeOwnership();}
	}
	
	@Override
	public final void clearOwnership(){initializeOwnership();}
	
	@Override
	public final int setOwner(final long kmer, final int newOwner){
		KmerLink n=get(kmer);
		assert(n!=null);
		if(n.owner<=newOwner){
			synchronized(n){
				if(n.owner<newOwner){
					n.owner=newOwner;
				}
			}
		}
		return n.owner;
	}
	
	@Override
	public final boolean clearOwner(final long kmer, final int owner){
		KmerLink n=get(kmer);
		assert(n!=null);
		synchronized(n){
			if(n.owner==owner){
				n.owner=-1;
				return true;
			}
		}
		return false;
	}
	
	@Override
	public final int getOwner(final long kmer){
		KmerLink n=get(kmer);
		assert(n!=null);
		return n.owner;
	}
	
	/*--------------------------------------------------------------*/
	
	@Override
	public int set(long kmer, int[] vals, int vlen) {
		throw new RuntimeException("Unimplemented.");
	}
	
	@Override
	public final int getValue(long kmer){
		KmerLink n=get(kmer);
		return n==null ? -1 : n.value;
	}
	
	@Override
	public final int[] getValues(long kmer, int[] singleton){
		KmerLink n=get(kmer);
		if(n==null){return null;}
		singleton[0]=n.value;
		return singleton;
	}
	
	@Override
	public final long size() {
		if(value<1){return 0;}
		long size=1;
		if(next!=null){size+=next.size();}
		return size;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		if(value<1){return true;}
		if(value>=mincount){
			if(remaining!=null && remaining.decrementAndGet()<0){return true;}
			bsw.printlnKmer(pivot, value, k);
		}
		if(next!=null){next.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);}
		return true;
	}
	
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
		if(next!=null){next.dumpKmersAsBytes_MT(bsw, bb, k, mincount, maxcount, remaining);}
		return true;
	}
	
	@Override
	public final boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount) {
		tsw.print(dumpKmersAsText(new StringBuilder(32), k, mincount, maxcount));
		return true;
	}
	
	private final StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount){
		if(value<1){return sb;}
		if(sb==null){sb=new StringBuilder(32);}
		if(value>=mincount){sb.append(AbstractKmerTable.toText(pivot, value, k)).append('\n');}
		if(next!=null){next.dumpKmersAsText(sb, k, mincount, maxcount);}
		return sb;
	}
	
	@Override
	public final void fillHistogram(long[] ca, int max){
		if(value<1){return;}
		ca[Tools.min(value, max)]++;
		if(next!=null){next.fillHistogram(ca, max);}
	}
	
	@Override
	public final void fillHistogram(SuperLongList sll){
		if(value<1){return;}
		sll.add(value);
		if(next!=null){next.fillHistogram(sll);}
	}
	
	@Override
	public void countGC(long[] gcCounts, int max){
		if(value<1){return;}
		gcCounts[Tools.min(value, max)]+=gc(pivot);
		if(next!=null){next.countGC(gcCounts, max);}
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	KmerLink rebalance(ArrayList<KmerLink> list){
		throw new RuntimeException("Unsupported.");
	}
	
	private static KmerLink rebalance(ArrayList<KmerLink> list, int a, int b){
		throw new RuntimeException("Unsupported.");
	}
	
	@Deprecated
	@Override
	public long regenerate(final int limit){
		throw new RuntimeException("TODO - remove zero-value links.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	long pivot;
	int value;
	int owner=-1;
	KmerLink next;
}
