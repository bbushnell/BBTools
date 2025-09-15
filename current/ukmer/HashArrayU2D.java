package ukmer;

import java.util.ArrayList;
import java.util.Arrays;

import shared.KillSwitch;
import shared.Primes;
import shared.Shared;
import shared.Tools;

/**
 * Stores kmers in a long[] and values in an int[][], with a victim cache.
 * @author Brian Bushnell
 * @date Nov 7, 2014
 *
 */
public final class HashArrayU2D extends HashArrayU {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public HashArrayU2D(int[] schedule_, int k_, int kbig_){
		super(schedule_, k_, kbig_, true);
		values=allocInt2D(prime+extra);
	}
	
//	public HashArrayU2D(int initialSize, int k_, int kbig_, boolean autoResize_){
//		super(initialSize, k_, kbig_, autoResize_, true);
//		values=allocInt2D(prime+extra);
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Deprecated
	@Override
	public int increment(final Kmer kmer){
		throw new RuntimeException("Unsupported.");
	}
	
	@Deprecated
	@Override
	public int incrementAndReturnNumCreated(final Kmer kmer){
		throw new RuntimeException("Unsupported.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	protected final int readCellValue(int cell) {
		int[] set=values[cell];
		return set==null ? 0 : set[0];
	}
	
	@Override
	protected final int[] readCellValues(int cell, int[] singleton) {
		return values[cell];
	}
	
	/** Returns number of values added */
	@Override
	protected final void insertValue(final long[] kmer, final int v, final int cell){
		assert(matches(kmer, cell));
		if(values[cell]==null){
			values[cell]=new int[] {v, NOT_PRESENT};
			return;
		}
		int[] set=values[cell];
		assert(set!=null);
		
		for(int i=0; i<set.length; i++){
			if(set[i]==v){return;}
			else if(set[i]<0){set[i]=v;return;}
		}
		final int oldSize=set.length;
		final int newSize=(int)Tools.min(Shared.MAX_ARRAY_LEN, oldSize*2L);
		assert(newSize>set.length) : "Overflow.";
		set=KillSwitch.copyOf(set, newSize);
		set[oldSize]=v;
		Arrays.fill(set, oldSize+1, newSize, NOT_PRESENT);
		values[cell]=set;
	}
	
	/** Returns number of values added */
	@Override
	protected final void insertValue(final long[] kmer, final int[] vals, final int cell){
		assert(matches(kmer, cell));
		if(values[cell]==null){
			values[cell]=vals;
		}else{
			for(int v : vals){
				if(v<0){break;}
				insertValue(kmer, v, cell);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public final boolean canRebalance() {return false;}
	
//	@Override
//	protected synchronized void resize(){
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
////		System.err.println("Resizing from "+prime+" to "+prime2+"; size="+size);
//		
//		prime=prime2;
////		System.err.println("Resized to "+prime+"; load="+(size*1f/prime));
//		long[][] oldk=arrays;
//		int[][] oldc=values;
//		arrays=allocLong2D(mult, prime2+extra);
//		for(int i=0; i<mult; i++){
//			Arrays.fill(arrays[i], NOT_PRESENT);
//		}
//		values=allocInt2D(prime2+extra);
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
//					set(fillKmer(i, kmer, oldk), oldc[i]);
//				}
//			}
//		}
//		
//		for(KmerNodeU n : list){
//			if(n.pivot[0]>NOT_PRESENT){
//				kmer.setFrom(n.pivot());
//				set(kmer, n.values(singleton));
//			}
//			else{assert(false);}
//		}
//		
//		assert(oldSize+oldVSize==size+victims.size) : oldSize+", "+oldVSize+" -> "+size+", "+victims.size;
//		
//		sizeLimit=(long)(maxLoadFactor*prime);
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
		int[][] oldc=values;
		arrays=allocLong2D(mult, prime+extra);
		for(int i=0; i<mult; i++){
			Arrays.fill(arrays[i], NOT_PRESENT);
		}
		values=allocInt2D(prime+extra);
		ArrayList<KmerNodeU> list=victims.toList();
		victims.clear();
		size=0;
		
		final int[] singleton=new int[] {NOT_PRESENT};
		final Kmer kmer=new Kmer(kbig);
		{
			for(int i=0; i<oldk.length; i++){
				if(oldk[0][i]>NOT_PRESENT){
					set(fillKmer(i, kmer, oldk), oldc[i]);
				}
			}
		}
		
		for(KmerNodeU n : list){
			if(n.pivot[0]>NOT_PRESENT){
				kmer.setFrom(n.pivot());
				set(kmer, n.values(singleton));
			}
			else{assert(false);}
		}
		
		assert(oldSize+oldVSize==size+victims.size) : oldSize+", "+oldVSize+" -> "+size+", "+victims.size;
	}
	
	@Deprecated
	@Override
	public void rebalance(){
		throw new RuntimeException("Unimplemented.");
	}
	
	@Deprecated
	@Override
	public long regenerate(final int limit){
		assert(false) : "This is not tested or intended for use.";
		long sum=0;
		assert(owners==null) : "Clear ownership before regeneration.";
		final Kmer kmer=new Kmer(kbig);
		for(int pos=0; pos<values.length; pos++){
			Kmer key=fillKmer(pos, kmer);
			if(key!=null){
				final int[] value=values[pos];
				values[pos]=null;
				arrays[0][pos]=NOT_PRESENT;
				size--;
				if(value!=null){
					assert(value[0]>0);
					set(key, value);
				}else{
					sum++;
				}
			}
		}
		
		ArrayList<KmerNodeU> nodes=victims.toList();
		victims.clear();
		for(KmerNodeU node : nodes){
			int value=node.value();
			if(value<1){
				sum++;
			}else{
				kmer.setFrom(node.pivot());
				set(kmer, node.values(null));//TODO: Probably unsafe or unwise.  Should test for singletons, etc.
			}
		}
		
		return sum;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int[][] values;
	

	
}
