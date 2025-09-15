package kmer;

import java.util.ArrayList;
import java.util.Arrays;

import shared.KillSwitch;
import shared.Primes;
import shared.Shared;
import shared.Tools;
import structures.IntList2;

/**
 * Stores kmers in a long[] and counts in an int[], with a victim cache.
 * @author Brian Bushnell
 * @date Oct 25, 2013
 *
 */
public final class HashArrayHybrid extends HashArray {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public HashArrayHybrid(int[] schedule_, long coreMask_){
		super(schedule_, coreMask_, true);
		values=allocInt1D(prime+extra);
		setList=new IntList2();
		setList.add(null, 0);
		setList.add(null, 0);
	}
	
//	public HashArrayHybrid(int initialSize, int maxSize, long mask, boolean autoResize_){
//		super(initialSize, maxSize, mask, autoResize_, true);
//		values=allocInt1D(prime+extra);
//		setList=new IntList2();
//		setList.add(null);
//		setList.add(null);
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final int increment(final long kmer, final int incr){
		int cell=kmerToCell(kmer);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			assert(n>-2);
			if(n==kmer){
				assert(values[cell]>=0);
				values[cell]+=incr;
				if(values[cell]<0){values[cell]=Integer.MAX_VALUE;}
				return values[cell];
			}else if(n==NOT_PRESENT){
				array[cell]=kmer;
				size++;
				values[cell]=incr;
				if(autoResize && size+victims.size>sizeLimit){resize();}
				return incr;
			}
		}
		int x=victims.increment(kmer, incr);
		if(autoResize && size+victims.size>sizeLimit){resize();}
		return x;
	}
	
	@Override
	public final int incrementAndReturnNumCreated(final long kmer, final int incr){
		int cell=kmerToCell(kmer);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			assert(n>-2);
			if(n==kmer){
				assert(values[cell]>=0);
				values[cell]+=incr;
				if(values[cell]<0){values[cell]=Integer.MAX_VALUE;}
				return 0;
			}else if(n==NOT_PRESENT){
				array[cell]=kmer;
				size++;
				values[cell]=incr;
				if(autoResize && size+victims.size>sizeLimit){resize();}
				return 1;
			}
		}
		int x=victims.incrementAndReturnNumCreated(kmer, incr);
		if(autoResize && size+victims.size>sizeLimit){resize();}
		return x;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	protected final int readCellValue(int cell) {
		final int x=values[cell];
		if(x>-2){return x;}
		return setList.get(0-x)[0];
	}
	
	@Override
	protected final int[] readCellValues(int cell, int[] singleton) {
		final int x=values[cell];
		if(x>-2){
			singleton[0]=values[cell];
			return singleton;
		}
		return setList.get(0-x);
	}
	
	@Override
	protected final void insertValue(long kmer, int[] vals, int cell, int vlen) {
		if(verbose){System.err.println("insertValue("+kmer+", "+Arrays.toString(vals)+", "+cell+"); old="+values[cell]);}
		assert(array[cell]==kmer);
		if(vals.length==1){
			if(verbose){System.err.println("A: length=1");}
			insertValue(kmer, vals[0], cell);
			return;
		}
		final int old=values[cell];
		if(old==vals[0] && vals[1]==NOT_PRESENT){//Maybe this is an attempt to accelerate the most common case?  It's only for when there is no array.
			if(verbose){System.err.println("B: old==vals[0] && vals[1]==-1");}
			return; //Nothing to do
		}else if(old<-1){//An array already exists
			if(verbose){System.err.println("C: old<-1");}
			for(int i : vals){
				if(i==NOT_PRESENT){break;}
				insertIntoList(i, -old);
			}
		}else{//Make a new array
			final int[] temp;
			if(old>0){//Move the old value to a new array.  Note that this will probably never be used.
				if(verbose){System.err.println("D: old>0");}
				temp=allocInt1D(vals.length+1);
				temp[0]=old;
				for(int i=0; i<vals.length; i++){temp[i+1]=vals[i];}
				vlen++;
			}else{
				if(verbose){System.err.println("E: old>0");}
				temp=vals;
			}
			values[cell]=-setList.size;
			setList.add(temp, vlen);
		}
	}
	
	@Override
	protected final void insertValue(long kmer, int v, int cell) {
		assert(array[cell]==kmer);
		assert(v>0);
		final int cc=values[cell];
		if(cc==v){
			return;
		}else if(cc<-1){
			insertIntoList(v, -cc);
		}else if(cc>0){
			values[cell]=-setList.size;
			setList.add(new int[] {cc, v, -1, -1}, 2);
		}else{
			values[cell]=v;
		}
	}
	
	@Deprecated
	private final int insertIntoListOld(final int v, final int loc){
		
		if(loc>=setList.size){
			assert(loc==setList.size);
			setList.add(null);
		}
			
		int[] set=setList.get(loc);
		if(set==null){
			set=new int[] {-1, -1};
			setList.set(loc, set);
		}
		
		for(int i=0; i<set.length; i++){//This is the slow bit
			if(set[i]==v){return 0;}
			if(set[i]<0){set[i]=v;return 1;}
		}
		final int oldSize=set.length;
		final int newSize=(int)Tools.min(Shared.MAX_ARRAY_LEN, oldSize*2L);
		assert(newSize>set.length) : "Overflow.";
		set=KillSwitch.copyOf(set, newSize);
		set[oldSize]=v;
		Arrays.fill(set, oldSize+1, newSize, -1);
		setList.set(loc, set);
		return 1;
	}
	
	private final int insertIntoList(final int v, final int loc){
		return setList.insertIntoList(v, loc);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final boolean canRebalance() {return false;}
	
	@Override
	protected synchronized void resize(){
		
		if(verbose){
			System.err.println("Resizing from "+prime+"; load="+(size*1f/prime));
		}
		
//		assert(TESTMODE);
//		if(TESTMODE){
//			for(int i=0; i<ll.size; i++){
//				assert(contains(ll.get(i), il.get(i)));
//				assert(!contains(ll.get(i), Integer.MAX_VALUE));
//			}
//		}
		
//		System.err.println("Resizing from "+prime+"; load="+(size*1f/prime));
		if(prime>=maxPrime){
//			sizeLimit=0xFFFFFFFFFFFFL;
//			return;
			KillSwitch.memKill(new OutOfMemoryError());
		}
		
		final long oldSize=size, oldVSize=victims.size;
		if(schedule!=null){
			final long oldPrime=prime;
			prime=nextScheduleSize();
			if(prime<=oldPrime){KillSwitch.memKill(new OutOfMemoryError());}
			sizeLimit=(long)((atMaxSize() ? maxLoadFactorFinal : maxLoadFactor)*prime);
		}else{//Old method
			final long totalSize=oldSize+oldVSize;

			final long maxAllowedByLoadFactor=(long)(totalSize*minLoadMult);
			final long minAllowedByLoadFactor=(long)(totalSize*maxLoadMult);

			//		sizeLimit=Tools.min((long)(maxLoadFactor*prime), maxPrime);

			assert(maxAllowedByLoadFactor>=minAllowedByLoadFactor);
			if(maxAllowedByLoadFactor<prime){
				sizeLimit=(long)(maxLoadFactor*prime);
				return;
			}

			long x=10+(long)(prime*resizeMult);
			x=Tools.max(x, minAllowedByLoadFactor);
			x=Tools.min(x, maxAllowedByLoadFactor);

			int prime2=(int)Tools.min(maxPrime, Primes.primeAtLeast(x));

			if(prime2<=prime){
				sizeLimit=(long)(maxLoadFactor*prime);
				assert(prime2==prime) : "Resizing to smaller array? "+totalSize+", "+prime+", "+x;
				return;
			}

			prime=prime2;
			sizeLimit=(long)(maxLoadFactor*prime);
		}
		
//		System.err.println("Resized to "+prime+"; load="+(size*1f/prime));
		long[] oldKmers=array;
		int[] oldValues=values;
		IntList2 oldList=setList;
		setList=new IntList2();
		setList.add(null, 0);
		setList.add(null, 0);
		KmerNode[] oldVictims=victims.array;
		array=allocLong1D(prime+extra);
		Arrays.fill(array, -1);
		values=allocInt1D(prime+extra);
		ArrayList<KmerNode> nodeList=new ArrayList<KmerNode>((int)(victims.size)); //Can fail if more than Integer.MAX_VALUE
		for(int i=0; i<oldVictims.length; i++){
			if(oldVictims[i]!=null){oldVictims[i].traverseInfix(nodeList);}
		}
		Arrays.fill(oldVictims, null);
		victims.size=0;
		size=0;
		
//		long added=0;
		for(int i=0; i<oldKmers.length; i++){
			final long kmer=oldKmers[i];
			if(kmer!=-1){
//				final int[] old=getValues(kmer, new int[1]);
//				final long oldsize=(size+victims.size);
//				assert(old==null);
				final int v=oldValues[i];
//				added++;
				
//				System.err.println("Found "+kmer+"->"+v);
				
				assert(v<-1 || v>0);
				if(v>=0){
//					assert(!contains(kmer));
//					long olds=size+victims.size; //123
					set(kmer, v);
//					assert(contains(kmer));
//					assert(size+victims.size==olds+1);
				}else{
//					if(verbose){
//						System.err.println("i="+i+", v="+v+", old="+Arrays.toString(oldList.get(-v))+", current="+Arrays.toString(setList.get(-v))+
//								", get()="+Arrays.toString(getValues(kmer, new int[1])));
//					}
//					assert(!contains(kmer));
//					long olds=size+victims.size; //123
					set(kmer, oldList.get(-v), oldList.getLen(-v));
//					assert(contains(kmer));
//					assert(size+victims.size==olds+1);
				}
			}
		}
//		assert(added==oldSize);
		
		final int[] singleton=new int[1];
//		added=0;
		for(KmerNode n : nodeList){
			if(n.pivot>-1){
//				added++;
//				final int[] old=getValues(n.pivot, new int[1]);
//				assert(old==null);
				if(n.numValues()>1){
//					assert(!contains(n.pivot()));
//					long olds=size+victims.size; //123
					set(n.pivot, n.values(singleton), n.numValues());
//					assert(size+victims.size==olds+1);
//					assert(contains(n.pivot()));
				}else{
//					assert(!contains(n.pivot()));
//					long olds=size+victims.size; //123
					set(n.pivot, n.value());
//					assert(size+victims.size==olds+1);
//					assert(contains(n.pivot()));
				}
//				assert(old==null || contains(n.pivot, old));
//				assert(contains(n.pivot, n.value()));
			}
		}
//		assert(added==oldVSize);
		
		assert(oldSize+oldVSize==size+victims.size) : oldSize+" + "+oldVSize+" = "+(oldSize+oldVSize)+" -> "+size+" + "+victims.size+" = "+(size+victims.size);
		
		if(verbose){System.err.println("Resized to "+prime+". "+oldSize+" + "+oldVSize+" = "+(oldSize+oldVSize)+" -> "+size+" + "+victims.size+" = "+(size+victims.size));}
		
//		assert(TESTMODE);
//		if(TESTMODE){
//			for(int i=0; i<ll.size; i++){
//				long kmer=ll.get(i);
//				int v=il.get(i);
//				assert(contains(kmer, v)) : i+", "+ll.size+", "+kmer+", "+v+", "+Arrays.toString(getValues(kmer, new int[1]));
//				assert(!contains(kmer, Integer.MAX_VALUE));
//			}
//		}
	}
	
	@Deprecated
	@Override
	public void rebalance(){
		throw new RuntimeException("Unimplemented.");
	}
	
	@Deprecated
	@Override
	public long regenerate(final int limit){
		throw new RuntimeException("Not supported.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int[] values;
	private IntList2 setList;
	
}
