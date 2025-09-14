package ukmer;

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
public final class HashArrayUHybrid extends HashArrayU {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public HashArrayUHybrid(int[] schedule_, int k_, int kbig_){
		super(schedule_, k_, kbig_, true);
		values=allocInt1D(prime+extra);
		setList=new IntList2();
		setList.add(null);
		setList.add(null);
	}
	
//	public HashArrayUHybrid(int initialSize, int k_, int kbig_, boolean autoResize_){
//		super(initialSize, k_, kbig_, autoResize_, true);
//		values=allocInt1D(prime+extra);
//		setList=new IntList2();
//		setList.add(null);
//		setList.add(null);
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final int increment(final Kmer kmer){
		final int cell=findKmerOrEmpty(kmer);
		
		if(cell==HASH_COLLISION){
			int x=victims.increment(kmer);
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return x;
		}else if(arrays[0][cell]==NOT_PRESENT){
			setKmer(kmer.key(), cell);
			size++;
			values[cell]=1;
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return 1;
		}else{
			values[cell]++;
			if(values[cell]<0){values[cell]=Integer.MAX_VALUE;}
			return values[cell];
		}
	}
	
	@Override
	public final int incrementAndReturnNumCreated(final Kmer kmer){
		final int cell=findKmerOrEmpty(kmer);
		
		if(cell==HASH_COLLISION){
			int x=victims.incrementAndReturnNumCreated(kmer);
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return x;
		}else if(arrays[0][cell]==NOT_PRESENT){
			setKmer(kmer.key(), cell);
			size++;
			values[cell]=1;
			if(autoResize && size+victims.size>sizeLimit){resize();}
			return 1;
		}else{
			values[cell]++;
			if(values[cell]<0){values[cell]=Integer.MAX_VALUE;}
			return 0;
		}
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
	protected final void insertValue(long[] kmer, int[] vals, int cell) {
		if(verbose){System.err.println("insertValue("+kmer+", "+Arrays.toString(vals)+", "+cell+"); old="+values[cell]);}
		assert(matches(kmer, cell));
		if(vals.length==1){
			if(verbose){System.err.println("A: length=1");}
			insertValue(kmer, vals[0], cell);
			return;
		}
		final int old=values[cell];
		if(old==vals[0] && vals[1]==NOT_PRESENT){
			if(verbose){System.err.println("B: old==vals[0] && vals[1]==-1");}
			return; //Nothing to do
		}else if(old<-1){//An array already exists
			if(verbose){System.err.println("C: old<-1");}
			for(int i : vals){
				if(i==-1){break;}
				insertIntoList(i, -old);
			}
		}else{//Add the list
			final int[] temp;
			if(old>0){//Move the old value to a new array.  Note that this will probably never be used.
				if(verbose){System.err.println("D: old>0");}
				temp=allocInt1D(vals.length+1);
				temp[0]=old;
				for(int i=0; i<vals.length; i++){temp[i+1]=vals[i];}
			}else{
				if(verbose){System.err.println("E: old>0");}
				temp=vals;
			}
			values[cell]=-setList.size;
			setList.add(temp);
		}
	}
	
	@Override
	protected final void insertValue(long[] kmer, int v, int cell) {
		assert(matches(kmer, cell));
		assert(v>0);
		final int cc=values[cell];
		if(cc==v){
			return;
		}else if(cc<-1){
			insertIntoList(v, -cc);
		}else if(cc>0){
			values[cell]=-setList.size;
			setList.add(new int[] {cc, v, -1, -1});
		}else{
			values[cell]=v;
		}
	}
	
	private final int insertIntoList(final int v, final int loc){
		
		if(loc>=setList.size){
			assert(loc==setList.size);
			setList.add(null);
		}
			
		int[] set=setList.get(loc);
		if(set==null){
			set=new int[] {-1, -1};
			setList.set(loc, set);
		}
		
		for(int i=0; i<set.length; i++){
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
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final boolean canRebalance() {return false;}
	
//	@Override
//	protected synchronized void resize(){
//		
//		if(verbose){
//			System.err.println("Resizing from "+prime+"; load="+(size*1f/prime));
//		}
//		
////		assert(TESTMODE);
////		if(TESTMODE){
////			for(int i=0; i<ll.size; i++){
////				assert(contains(ll.get(i), il.get(i)));
////				assert(!contains(ll.get(i), Integer.MAX_VALUE));
////			}
////		}
//		
////		System.err.println("Resizing from "+prime+"; load="+(size*1f/prime));
//		if(prime>=maxPrime){
//			sizeLimit=0xFFFFFFFFFFFFL;
//			return;
//		}
//		
//		final long oldSize=size, oldVSize=victims.size;
//		final long totalSize=this.size+victims.size;
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
//		long[][] oldk=arrays;
//		int[] oldc=values;
//		arrays=KillSwitch.allocLong2D(mult, prime2+extra);
//		for(int i=0; i<mult; i++){
//			Arrays.fill(arrays[i], NOT_PRESENT);
//		}
//		IntList2 oldList=setList;
//		setList=new IntList2();
//		setList.add(null);
//		setList.add(null);//TODO: May have to add 3 of them to avoid HASH_COLLISION at -2
//		values=allocInt1D(prime2+extra);
//		ArrayList<KmerNodeU> list=victims.toList();
//		victims.clear();
//		size=0;
//		sizeLimit=Long.MAX_VALUE;
//		
//		final int[] singleton=new int[] {NOT_PRESENT};
//		final Kmer kmer=new Kmer(kbig);
//		{
//			for(int i=0; i<oldk.length; i++){
//				if(oldk[0][i]>NOT_PRESENT){
//					final int v=oldc[i];
//					fillKmer(i, kmer, oldk);
//					if(v>=0){
//						set(kmer, v);
//					}else{
//						set(kmer, oldList.get(-v));
//					}
//				}
//			}
//		}
//		
//		for(KmerNodeU n : list){
//			if(n.pivot[0]>NOT_PRESENT){
//				kmer.setFrom(n.pivot());
//				if(n.numValues()>1){
//					set(kmer, n.values(singleton));
//				}else{
//					set(kmer, n.value());
//				}
//			}else{assert(false);}
//		}
//		
//		assert(oldSize+oldVSize==size+victims.size) : oldSize+" + "+oldVSize+" = "+(oldSize+oldVSize)+" -> "+size+" + "+victims.size+" = "+(size+victims.size);
//		
//		if(verbose){System.err.println("Resized to "+prime+". "+oldSize+" + "+oldVSize+" = "+(oldSize+oldVSize)+" -> "+size+" + "+victims.size+" = "+(size+victims.size));}
//		
//		sizeLimit=(long)(maxLoadFactor*prime);
//		
////		assert(TESTMODE);
////		if(TESTMODE){
////			for(int i=0; i<ll.size; i++){
////				long[] kmer=ll.get(i);
////				int v=il.get(i);
////				assert(contains(kmer, v)) : i+", "+ll.size+", "+kmer+", "+v+", "+Arrays.toString(getValues(kmer, new int[1]));
////				assert(!contains(kmer, Integer.MAX_VALUE));
////			}
////		}
//	}

	
	@Override
	protected synchronized void resize(){
		if(verbose){System.err.println("Resizing from "+prime+"; load="+(size*1f/prime));}
		if(prime>=maxPrime){
//			sizeLimit=0xFFFFFFFFFFFFL;
//			return;
			KillSwitch.memKill(new OutOfMemoryError());
		}

		final int oldPrime=prime;
		final long oldSize=size, oldVSize=victims.size;
		final long totalSize=oldSize+oldVSize;
		
		if(schedule!=null){
			prime=nextScheduleSize();
			if(prime<=oldPrime){KillSwitch.memKill(new OutOfMemoryError());}
			sizeLimit=(long)((atMaxSize() ? maxLoadFactorFinal : maxLoadFactor)*prime);
		}else{//Old method
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
		long[][] oldk=arrays;
		int[] oldc=values;
		arrays=KillSwitch.allocLong2D(mult, prime+extra);
		for(int i=0; i<mult; i++){
			Arrays.fill(arrays[i], NOT_PRESENT);
		}
		IntList2 oldList=setList;
		setList=new IntList2();
		setList.add(null);
		setList.add(null);//TODO: May have to add 3 of them to avoid HASH_COLLISION at -2
		values=allocInt1D(prime+extra);
		ArrayList<KmerNodeU> list=victims.toList();
		victims.clear();
		size=0;
		
		final int[] singleton=new int[] {NOT_PRESENT};
		final Kmer kmer=new Kmer(kbig);
		{
			for(int i=0; i<oldk.length; i++){
				if(oldk[0][i]>NOT_PRESENT){
					final int v=oldc[i];
					fillKmer(i, kmer, oldk);
					if(v>=0){
						set(kmer, v);
					}else{
						set(kmer, oldList.get(-v));
					}
				}
			}
		}
		
		for(KmerNodeU n : list){
			if(n.pivot[0]>NOT_PRESENT){
				kmer.setFrom(n.pivot());
				if(n.numValues()>1){
					set(kmer, n.values(singleton));
				}else{
					set(kmer, n.value());
				}
			}else{assert(false);}
		}
		
		assert(oldSize+oldVSize==size+victims.size) : oldSize+" + "+oldVSize+" = "+(oldSize+oldVSize)+" -> "+size+" + "+victims.size+" = "+(size+victims.size);
		
		if(verbose){System.err.println("Resized to "+prime+". "+oldSize+" + "+oldVSize+" = "+(oldSize+oldVSize)+" -> "+size+" + "+victims.size+" = "+(size+victims.size));}
		
		sizeLimit=(long)(maxLoadFactor*prime);
		
//		assert(TESTMODE);
//		if(TESTMODE){
//			for(int i=0; i<ll.size; i++){
//				long[] kmer=ll.get(i);
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
