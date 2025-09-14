package ukmer;

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
public abstract class KmerNodeU extends AbstractKmerTableU {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	protected KmerNodeU(long[] pivot_){
		pivot=pivot_.clone();
	}
	
	public abstract KmerNodeU makeNode(long[] pivot_, int value_);
	public abstract KmerNodeU makeNode(long[] pivot_, int[] values_);
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final int increment(Kmer kmer){return increment(kmer.key());}
	
	public final int increment(long[] kmer){
		final int cmp=compare(kmer, pivot);
		if(cmp<0){
			if(left==null){left=makeNode(kmer, 1); return 1;}
			return left.increment(kmer);
		}else if(cmp>0){
			if(right==null){right=makeNode(kmer, 1); return 1;}
			return right.increment(kmer);
		}else{
			if(value()<Integer.MAX_VALUE){set(value()+1);}
			return value();
		}
	}
	
	@Override
	public final int incrementAndReturnNumCreated(Kmer kmer){return incrementAndReturnNumCreated(kmer.key());}
	
	public final int incrementAndReturnNumCreated(long[] kmer) {
		int x=increment(kmer);
		return x==1 ? 1 : 0;
	}
	
	/** Returns number of nodes added */
	public final int set(long[] kmer, int value){
		if(verbose){System.err.println("Set0: kmer="+Arrays.toString(kmer)+", v="+value+", old="+Arrays.toString(values(new int[1])));}
		if(verbose){System.err.println("A");}
		final int cmp=compare(kmer, pivot);
		if(cmp<0){
			if(verbose){System.err.println("B");}
			if(left==null){left=makeNode(kmer, value); return 1;}
			if(verbose){System.err.println("C");}
			return left.set(kmer, value);
		}else if(cmp>0){
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
	public final int setIfNotPresent(long[] kmer, int value){
		if(verbose){System.err.println("setIfNotPresent0: kmer="+kmer+", v="+value+", old="+Arrays.toString(values(new int[0])));}
		final int cmp=compare(kmer, pivot);
		if(cmp<0){
			if(left==null){left=makeNode(kmer, value); return 1;}
			return left.setIfNotPresent(kmer, value);
		}else if(cmp>0){
			if(right==null){right=makeNode(kmer, value); return 1;}
			return right.setIfNotPresent(kmer, value);
		}
		return 0;
	}
	
	public final int getValue(long[] kmer){
		KmerNodeU n=get(kmer);
		return n==null ? -1 : n.value();
	}
	
	public final int[] getValues(long[] kmer, int[] singleton){
		KmerNodeU n=get(kmer);
		return n==null ? null : n.values(singleton);
	}
	
	public final boolean contains(long[] kmer){
		KmerNodeU node=get(kmer);
		return node!=null;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	public KmerNodeU left(){return left;}
	public KmerNodeU right(){return right;}
	public long[] pivot(){return pivot;}
	public int owner(){return owner;}
	
	public Kmer fillKmer(Kmer x){
		assert(pivot!=null);
		long[] key=x.array1();
		for(int i=0; i<pivot.length; i++){
			key[i]=pivot[i];
		}
		x.fillArray2();
		return x;
	}
	
	public int count(){return value();}
	protected abstract int value();
	protected abstract int[] values(int[] singleton);
	/** Returns new value */
	public abstract int set(int value_);
	protected abstract int set(int[] values_);
	
	@Override
	final KmerNodeU get(final long[] kmer){
//		if(kmer<pivot){
//			return left==null ? null : left.get(kmer);
//		}else if(kmer>pivot){
//			return right==null ? null : right.get(kmer);
//		}else{
//			return this;
//		}
		KmerNodeU n=this;
		int cmp=compare(kmer, n.pivot);
		while(cmp!=0){
			n=(cmp<0 ? n.left : n.right);
			cmp=(n==null ? 0 : compare(kmer, n.pivot));
		}
		return n;
	}
	
	final KmerNodeU getNodeOrParent(long[] kmer){
		final int cmp=compare(kmer, pivot);
		if(cmp==0){return this;}
		if(cmp<0){return left==null ? this : left.getNodeOrParent(kmer);}
		return right==null ? this : right.getNodeOrParent(kmer);
	}
	
	final boolean insert(KmerNodeU n){
		assert(pivot!=null);
		final int cmp=compare(n.pivot, pivot);
		if(cmp<0){
			if(left==null){left=n; return true;}
			return left.insert(n);
		}else if(cmp>0){
			if(right==null){right=n; return true;}
			return right.insert(n);
		}else{
			return false;
		}
	}
	
	final void traversePrefix(ArrayList<KmerNodeU> list){
		if(left!=null){left.traversePrefix(list);}
		list.add(this);
		if(right!=null){right.traversePrefix(list);}
	}
	
	final void traverseInfix(ArrayList<KmerNodeU> list){
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
	
	final KmerNodeU rebalance(ArrayList<KmerNodeU> list){
		assert(list.isEmpty());
		traversePrefix(list);
		KmerNodeU n=this;
		if(list.size()>2){
			n=rebalance(list, 0, list.size()-1);
		}
		list.clear();
		return n;
	}
	
	private static final KmerNodeU rebalance(ArrayList<KmerNodeU> list, int a, int b){
		final int size=b-a+1;
		final int middle=a+size/2;
		final KmerNodeU n=list.get(middle);
		if(size<4){
			if(size==1){
				n.left=n.right=null;
			}else if(size==2){
				KmerNodeU n1=list.get(a);
				n.left=n1;
				n.right=null;
				n1.left=n1.right=null;
			}else{
				assert(size==3);
				KmerNodeU n1=list.get(a), n2=list.get(b);
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
	public final void countGC(long[] gcCounts, int max){
		final int value=value();
		if(value<1){return;}
		int index=Tools.min(value, max);
		for(long x : pivot){
			gcCounts[index]+=gc(x);
		}
		if(left!=null){left.countGC(gcCounts, max);}
		if(right!=null){right.countGC(gcCounts, max);}
	}
	
	@Override
	public String toString(){return Arrays.toString(pivot);}

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
	
	
	public final int setOwner(final long[] kmer, final int newOwner){
		KmerNodeU n=get(kmer);
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
	
	
	public final boolean clearOwner(final long[] kmer, final int owner){
		KmerNodeU n=get(kmer);
		assert(n!=null);
		synchronized(n){
			if(n.owner==owner){
				n.owner=-1;
				return true;
			}
		}
		return false;
	}
	
	
	public final int getOwner(final long[] kmer){
		KmerNodeU n=get(kmer);
		assert(n!=null);
		return n.owner;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Recall Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	abstract int set(long[] kmer, int[] vals);
	
	@Override
	public int set(Kmer kmer, int value) {
		return set(kmer.key(), value);
	}
	
	@Override
	public int set(Kmer kmer, int[] vals) {
		return set(kmer.key(), vals);
	}
	
	@Override
	public int setIfNotPresent(Kmer kmer, int value) {
		return setIfNotPresent(kmer.key(), value);
	}
	
	@Override
	public int getValue(Kmer kmer) {
		return getValue(kmer.key());
	}
	
	@Override
	public int[] getValues(Kmer kmer, int[] singleton) {
		return getValues(kmer.key(), singleton);
	}
	
	@Override
	public boolean contains(Kmer kmer) {
		return contains(kmer.key());
	}
	
	@Override
	public int getValue(long[] key, long xor) {
		return getValue(key);
	}
	
	@Override
	public int setOwner(Kmer kmer, int newOwner) {
		return setOwner(kmer.key(), newOwner);
	}
	
	@Override
	public boolean clearOwner(Kmer kmer, int owner) {
		return clearOwner(kmer.key(), owner);
	}
	
	@Override
	public int getOwner(Kmer kmer) {
		return getOwner(kmer.key());
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------       Invalid Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	final long[] pivot;
	int owner=-1;
	KmerNodeU left, right;
	
}
