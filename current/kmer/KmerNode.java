package kmer;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.TextStreamWriter;
import shared.Tools;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * @author Brian Bushnell
 * @date Oct 22, 2013
 *
 */
public abstract class KmerNode extends AbstractKmerTable {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	protected KmerNode(long pivot_){
		pivot=pivot_;
	}
	
	public abstract KmerNode makeNode(long pivot_, int value_);
	public abstract KmerNode makeNode(long pivot_, int[] values_, int vlen_);
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final int increment(final long kmer, final int incr){
		if(pivot<0){pivot=kmer; return set(incr);} //Allows initializing empty nodes to -1
		if(kmer<pivot){
			if(left==null){left=makeNode(kmer, incr); return incr;}
			return left.increment(kmer, incr);
		}else if(kmer>pivot){
			if(right==null){right=makeNode(kmer, incr); return incr;}
			return right.increment(kmer, incr);
		}else{
			if(value()<Integer.MAX_VALUE){set(value()+incr);}
			return value();
		}
	}
	
	@Override
	public final int incrementAndReturnNumCreated(final long kmer, final int incr) {
		int x=increment(kmer, incr);
		return x==1 ? 1 : 0;
	}
	
//	public final int set_Test(final long kmer, final int v){
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
//			assert(v==old || !contains(kmer, old));
//		}
//		return x;
//	}
//
//	public final int setIfNotPresent_Test(long kmer, int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD()){
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
	
	
	/** Returns number of nodes added */
	@Override
	public final int set(long kmer, int value){
		if(verbose){System.err.println("Set0: kmer="+kmer+", v="+value+", old="+Arrays.toString(values(new int[1])));}
		if(pivot<0){pivot=kmer; set(value); return 1;} //Allows initializing empty nodes to -1
		if(verbose){System.err.println("A");}
		if(kmer<pivot){
			if(verbose){System.err.println("B");}
			if(left==null){left=makeNode(kmer, value); return 1;}
			if(verbose){System.err.println("C");}
			return left.set(kmer, value);
		}else if(kmer>pivot){
			if(verbose){System.err.println("D");}
			if(right==null){right=makeNode(kmer, value); return 1;}
			if(verbose){System.err.println("E");}
			return right.set(kmer, value);
		}else{
			if(verbose){System.err.println("F");}
			set(value);
		}
		if(verbose){System.err.println("G");}
		return 0;
	}
	
	
	/** Returns number of nodes added */
	@Override
	public final int setIfNotPresent(long kmer, int value){
		if(verbose){System.err.println("setIfNotPresent0: kmer="+kmer+", v="+value+", old="+Arrays.toString(values(new int[0])));}
		if(pivot<0){pivot=kmer; set(value); return 1;} //Allows initializing empty nodes to -1
		if(kmer<pivot){
			if(left==null){left=makeNode(kmer, value); return 1;}
			return left.setIfNotPresent(kmer, value);
		}else if(kmer>pivot){
			if(right==null){right=makeNode(kmer, value); return 1;}
			return right.setIfNotPresent(kmer, value);
		}
		return 0;
	}
	
	@Override
	public final int getValue(long kmer){
		KmerNode n=get(kmer);
		return n==null ? -1 : n.value();
	}
	
	@Override
	public final int[] getValues(long kmer, int[] singleton){
		KmerNode n=get(kmer);
		return n==null ? null : n.values(singleton);
	}
	
	@Override
	public final boolean contains(long kmer){
		KmerNode node=get(kmer);
		return node!=null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	public KmerNode left(){return left;}
	public KmerNode right(){return right;}
	public long pivot(){return pivot;}
	public int owner(){return owner;}
	
	public int count(){return value();}
	protected abstract int value();
	protected abstract int[] values(int[] singleton);
	/** Returns new value */
	public abstract int set(int value_);
	protected abstract int set(int[] values_, int vlen);
	
	@Override
	final KmerNode get(long kmer){
//		if(kmer<pivot){
//			return left==null ? null : left.get(kmer);
//		}else if(kmer>pivot){
//			return right==null ? null : right.get(kmer);
//		}else{
//			return this;
//		}
		KmerNode n=this;
		while(n!=null && n.pivot!=kmer){
			n=(kmer<n.pivot ? n.left : n.right);
		}
		return n;
	}
	
	final KmerNode getNodeOrParent(long kmer){
		if(pivot==kmer || pivot<0){return this;}
		if(kmer<pivot){return left==null ? this : left.getNodeOrParent(kmer);}
		return right==null ? this : right.getNodeOrParent(kmer);
	}
	
	final boolean insert(KmerNode n){
		assert(pivot!=-1);
		if(n.pivot<pivot){
			if(left==null){left=n; return true;}
			return left.insert(n);
		}else if(n.pivot>pivot){
			if(right==null){right=n; return true;}
			return right.insert(n);
		}else{
			return false;
		}
	}
	
	final void traversePrefix(ArrayList<KmerNode> list){
		if(left!=null){left.traversePrefix(list);}
		list.add(this);
		if(right!=null){right.traversePrefix(list);}
	}
	
	final void traverseInfix(ArrayList<KmerNode> list){
		list.add(this);
		if(left!=null){left.traverseInfix(list);}
		if(right!=null){right.traverseInfix(list);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Private Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final long size() {
		if(value()<1){return 0;}
		long size=1;
		if(left!=null){size+=left.size();}
		if(right!=null){size+=right.size();}
		return size;
	}
	
	final KmerNode rebalance(ArrayList<KmerNode> list){
		assert(list.isEmpty());
		traversePrefix(list);
		KmerNode n=this;
		if(list.size()>2){
			n=rebalance(list, 0, list.size()-1);
		}
		list.clear();
		return n;
	}
	
	private static final KmerNode rebalance(ArrayList<KmerNode> list, int a, int b){
		final int size=b-a+1;
		final int middle=a+size/2;
		final KmerNode n=list.get(middle);
		if(size<4){
			if(size==1){
				n.left=n.right=null;
			}else if(size==2){
				KmerNode n1=list.get(a);
				n.left=n1;
				n.right=null;
				n1.left=n1.right=null;
			}else{
				assert(size==3);
				KmerNode n1=list.get(a), n2=list.get(b);
				n.left=n1;
				n.right=n2;
				n1.left=n1.right=null;
				n2.left=n2.right=null;
			}
		}else{
			n.left=rebalance(list, a, middle-1);
			n.right=rebalance(list, middle+1, b);
		}
		return n;
	}
	
	@Override
	public long regenerate(final int limit){
		throw new RuntimeException("Not supported.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount) {
		tsw.print(dumpKmersAsText(new StringBuilder(32), k, mincount, maxcount));
		return true;
	}
	
	protected abstract StringBuilder dumpKmersAsText(StringBuilder sb, int k, int mincount, int maxcount);
	
	protected abstract ByteBuilder dumpKmersAsText(ByteBuilder bb, int k, int mincount, int maxcount);
	
	@Override
	public final void fillHistogram(long[] ca, int max){
		final int value=value();
		if(value<1){return;}
		ca[Tools.min(value, max)]++;
		if(left!=null){left.fillHistogram(ca, max);}
		if(right!=null){right.fillHistogram(ca, max);}
	}
	
	@Override
	public final void fillHistogram(SuperLongList sll){
		final int value=value();
		if(value<1){return;}
		sll.add(value);
		if(left!=null){left.fillHistogram(sll);}
		if(right!=null){right.fillHistogram(sll);}
	}
	
	@Override
	public void countGC(long[] gcCounts, int max){
		final int value=value();
		if(value<1){return;}
		gcCounts[Tools.min(value, max)]+=gc(pivot);
		if(left!=null){left.countGC(gcCounts, max);}
		if(right!=null){right.countGC(gcCounts, max);}
	}

	abstract boolean TWOD();
	abstract int numValues();
	
	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final void initializeOwnership(){
		owner=-1;
		if(left!=null){left.initializeOwnership();}
		if(right!=null){right.initializeOwnership();}
	}
	
	@Override
	public final void clearOwnership(){initializeOwnership();}
	
	@Override
	public final int setOwner(final long kmer, final int newOwner){
		KmerNode n=get(kmer);
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
		KmerNode n=get(kmer);
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
		KmerNode n=get(kmer);
		assert(n!=null);
		return n.owner;
	}
	
	public long calcMem() {//48 is a guess
		return 48+(left==null ? 0 : left.calcMem())+(right==null ? 0 : right.calcMem());
	}
	
	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	long pivot;
	int owner=-1;
	KmerNode left, right;
	
}
