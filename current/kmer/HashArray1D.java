package kmer;

import java.util.ArrayList;
import java.util.Arrays;

import shared.KillSwitch;
import shared.Primes;
import shared.Tools;
import structures.SuperLongList;

/**
 * Stores kmers in a long[] and counts in an int[], with a victim cache.
 * @author Brian Bushnell
 * @date Oct 25, 2013
 *
 */
public final class HashArray1D extends HashArray {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public HashArray1D(int[] schedule_, long coreMask_){
		super(schedule_, coreMask_, false);
		values=allocInt1D(prime+extra);
	}
	
	public HashArray1D(int initialSize, long coreMask, boolean autoResize_){
		super(initialSize, coreMask, autoResize_, false);
		values=allocInt1D(prime+extra);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final int increment(final long kmer, final int incr){
		int cell=kmerToCell(kmer);
		
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				values[cell]+=incr;
				if(values[cell]<0){values[cell]=Integer.MAX_VALUE;}
				return values[cell];
			}else if(n==NOT_PRESENT){
				array[cell]=kmer;
				size++;
				values[cell]=incr;
				if(autoResize && size+victims.size>sizeLimit){resize();}
				return 1;
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
			if(n==kmer){
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
		return victims.incrementAndReturnNumCreated(kmer, incr);
	}
	
	@Override
	public void fillHistogram(SuperLongList sll){
		for(int i=0; i<values.length; i++){
			int count=values[i];
			if(count>0){sll.add(count);}
		}
		if(victims!=null){
			victims.fillHistogram(sll);
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public final int readCellValue(int cell) {
		return values[cell];
	}
	
	@Override
	protected final int[] readCellValues(int cell, int[] singleton) {
		singleton[0]=values[cell];
		return singleton;
	}
	
	@Override
	protected final void insertValue(long kmer, int v, int cell) {
		assert(array[cell]==kmer);
		values[cell]=v;
	}
	
	@Override
	protected final void insertValue(long kmer, int[] vals, int cell, int vlen) {
		assert(array[cell]==kmer);
		assert(vals.length==1);
		values[cell]=vals[0];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final boolean canRebalance() {return false;}
	
//	@Override
//	protected synchronized void resize_old(){
////		assert(false);
////		System.err.println("Resizing from "+prime+"; load="+(size*1f/prime));
//		if(prime>=maxPrime){
//			sizeLimit=0xFFFFFFFFFFFFL;
//			return;
//		}
//		
//		final long oldSize=size, oldVSize=victims.size;
//		final long totalSize=oldSize+oldVSize;
//		
//		final long maxAllowedByLoadFactor=(long)(totalSize*minLoadMult);
//		final long minAllowedByLoadFactor=(long)(totalSize*maxLoadMult);
//
////		sizeLimit=Tools.min((long)(maxLoadFactor*prime), maxPrime);
//		
//		assert(maxAllowedByLoadFactor>=minAllowedByLoadFactor);
//		if(maxAllowedByLoadFactor<prime){
//			sizeLimit=(long)(maxLoadFactor*prime);
//			return;
//		}
//		
//		long x=10+(long)(prime*resizeMult);
//		x=Tools.max(x, minAllowedByLoadFactor);
//		x=Tools.min(x, maxAllowedByLoadFactor);
//		
//		int prime2=(int)Tools.min(maxPrime, Primes.primeAtLeast(x));
//		
//		if(prime2<=prime){
//			sizeLimit=(long)(maxLoadFactor*prime);
//			assert(prime2==prime) : "Resizing to smaller array? "+totalSize+", "+prime+", "+x;
//			return;
//		}
//		
//		prime=prime2;
////		System.err.println("Resized to "+prime+"; load="+(size*1f/prime));
//		long[] oldk=array;
//		int[] oldc=values;
//		KmerNode[] oldv=victims.array;
//		array=allocLong1D(prime2+extra);
//		Arrays.fill(array, NOT_PRESENT);
//		values=allocInt1D(prime2+extra);
//		ArrayList<KmerNode> list=victims.toList();
//		Arrays.fill(oldv, null);
//		victims.size=0;
//		size=0;
//		sizeLimit=Long.MAX_VALUE;
//		
//		if(TWO_PASS_RESIZE){
//			for(int i=0; i<oldk.length; i++){
//				if(oldk[i]>NOT_PRESENT && oldc[i]>1){set(oldk[i], oldc[i]);}
//			}
//			for(KmerNode n : list){
//				if(n.pivot>NOT_PRESENT && n.value()>1){set(n.pivot, n.value());}
//			}
//			for(int i=0; i<oldk.length; i++){
//				if(oldk[i]>NOT_PRESENT && oldc[i]<=1){set(oldk[i], oldc[i]);}
//			}
//			for(KmerNode n : list){
//				if(n.pivot>NOT_PRESENT && n.value()<=1){set(n.pivot, n.value());}
//			}
//		}else{
//			for(int i=0; i<oldk.length; i++){
//				if(oldk[i]>NOT_PRESENT){set(oldk[i], oldc[i]);}
//			}
//			for(KmerNode n : list){
//				if(n.pivot>NOT_PRESENT){set(n.pivot, n.value());}
//			}
//		}
//		
//		assert(oldSize+oldVSize==size+victims.size) : oldSize+", "+oldVSize+" -> "+size+", "+victims.size;
//		
//		sizeLimit=(long)(maxLoadFactor*prime);
//	}
	
	@Override
	protected synchronized void resize(){
//		assert(false);
//		System.err.println("Resizing from "+prime+"; load="+(size*1f/prime));
		if(prime>=maxPrime){
//			sizeLimit=0xFFFFFFFFFFFFL;
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
		long[] oldk=array;
		int[] oldc=values;
		KmerNode[] oldv=victims.array;
		array=allocLong1D(prime+extra);
		Arrays.fill(array, NOT_PRESENT);
//		System.err.print(prime+" ");//123
		values=allocInt1D(prime+extra);
		ArrayList<KmerNode> list=victims.toList();
		Arrays.fill(oldv, null);
		victims.size=0;
		size=0;
		
		if(TWO_PASS_RESIZE){
			for(int i=0; i<oldk.length; i++){
				if(oldk[i]>NOT_PRESENT && oldc[i]>1){set(oldk[i], oldc[i]);}
			}
			for(KmerNode n : list){
				if(n.pivot>NOT_PRESENT && n.value()>1){set(n.pivot, n.value());}
			}
			for(int i=0; i<oldk.length; i++){
				if(oldk[i]>NOT_PRESENT && oldc[i]<=1){set(oldk[i], oldc[i]);}
			}
			for(KmerNode n : list){
				if(n.pivot>NOT_PRESENT && n.value()<=1){set(n.pivot, n.value());}
			}
		}else{
			for(int i=0; i<oldk.length; i++){
				if(oldk[i]>NOT_PRESENT){set(oldk[i], oldc[i]);}
			}
			for(KmerNode n : list){
				if(n.pivot>NOT_PRESENT){set(n.pivot, n.value());}
			}
		}
		
		assert(oldSize+oldVSize==size+victims.size) : oldSize+", "+oldVSize+" -> "+size+", "+victims.size;
	}
	
	@Deprecated
	@Override
	public void rebalance(){
		throw new RuntimeException("Unimplemented.");
	}
	
	@Override
	public long regenerate(final int limit){
		long sum=0;
		assert(owners==null) : "Clear ownership before regeneration.";
		for(int pos=0; pos<values.length; pos++){
			final long key=array[pos];
			if(key>=0){
				final int value=values[pos];
				values[pos]=NOT_PRESENT;
				array[pos]=NOT_PRESENT;
				size--;
				if(value>limit){
					set(key, value);
				}else{
					sum++;
				}
			}
		}
		
		ArrayList<KmerNode> nodes=victims.toList();
		victims.clear();
		for(KmerNode node : nodes){
			int value=node.value();
			if(value<=limit){
				sum++;
			}else{
				set(node.pivot, node.value());
			}
		}
		
		return sum;
	}
	
	@Override
	public String toString(){
		return Arrays.toString(array);
	}
	
	public Walker walk(){
		return new Walker1D();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int[] values;
	
	public int[] values(){return values;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Walker            ----------------*/
	/*--------------------------------------------------------------*/
	
	public class Walker1D extends Walker {
		
		Walker1D(){
			ha=HashArray1D.this;
			victims=ha.victims().toList();
		}
		
		/** 
		 * Fills this object with the next key and value.
		 * @return True if successful.
		 */
		public boolean next(){
			while(i<values.length && values[i]==NOT_XPRESENT){i++;}
			if(i<values.length){
				kmer=array[i];
				value=values[i];
				assert(value!=NOT_XPRESENT);
				i++;
				return true;
			}
			if(i2<victims.size()){
				KmerNode kn=victims.get(i2);
				kmer=kn.pivot;
				value=kn.value();
				assert(value!=NOT_XPRESENT);
				i2++;
				return true;
			}
			kmer=-1;
			value=NOT_XPRESENT;
			return false;
		}
		
		public long kmer(){return kmer;}
		public int value(){return value;}
		
		/** Hash map over which this is walking */
		private HashArray1D ha;
		/** Victim list of the hash map */
		private ArrayList<KmerNode> victims;
		
		private long kmer;
		private int value;

		/** Potential next kmer cell; may point to an empty cell */
		private int i=0;
		/** Next victim in list */
		private int i2=0;
	}
	
	//TODO: Remove after fixing array initialization
	private static final int NOT_XPRESENT=0;

	public long calcMem() {
		long mem=0;
		mem+=(array.length*8);
		mem+=(values.length*4);
		mem+=(owners==null ? 0 : owners.length()*4);
		for(KmerNode kn : victims.array){
			mem+=8;
			if(kn!=null){mem+=kn.calcMem();}
		}
		// TODO Auto-generated method stub
		return mem;
	}

	
}
