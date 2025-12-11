package kmer;

import java.util.ArrayList;
import java.util.Arrays;

import shared.KillSwitch;
import shared.Primes;
import shared.Shared;
import shared.Tools;

/**
 * Stores k-mers in a long[] and values in a 2D int[][] array, supporting multiple
 * integer values per k-mer with dynamic array growth. Uses linear probing for
 * collision resolution and a victim cache for overflow handling.
 *
 * @author Brian Bushnell
 * @date Nov 7, 2014
 */
public final class HashArray2D extends HashArray {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public HashArray2D(int[] schedule_, long coreMask_){
		super(schedule_, coreMask_, true);
		values=allocInt2D(prime+extra);
	}
	
//	public HashArray2D(int initialSize, int maxSize, long mask, boolean autoResize_){
//		super(initialSize, maxSize, mask, autoResize_, true);
//		values=allocInt2D(prime+extra);
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Deprecated increment operation; not supported for 2D value storage.
	 * Always throws RuntimeException to prevent misuse.
	 *
	 * @param kmer K-mer key
	 * @param incr Increment amount (ignored)
	 * @return never returns
	 * @throws RuntimeException always thrown
	 * @deprecated Use insertValue instead of incrementing
	 */
	@Deprecated
	@Override
	public int increment(final long kmer, final int incr){
		throw new RuntimeException("Unsupported.");
	}
	
	/**
	 * Deprecated increment operation; not supported for 2D value storage.
	 * Always throws RuntimeException to prevent misuse.
	 *
	 * @param kmer K-mer key
	 * @param incr Increment amount (ignored)
	 * @return never returns
	 * @throws RuntimeException always thrown
	 * @deprecated Use insertValue instead of incrementing
	 */
	@Deprecated
	@Override
	public int incrementAndReturnNumCreated(final long kmer, final int incr){
		throw new RuntimeException("Unsupported.");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Reads the first stored value for the given cell, or 0 if no values exist.
	 * @param cell Array index
	 * @return First integer value or 0 if empty
	 */
	@Override
	protected final int readCellValue(int cell) {
		int[] set=values[cell];
		return set==null ? 0 : set[0];
	}
	
	/**
	 * Returns the full value array stored in the cell (may be null).
	 * @param cell Array index
	 * @param singleton Ignored placeholder to match superclass signature
	 * @return Value array or null
	 */
	@Override
	protected final int[] readCellValues(int cell, int[] singleton) {
		return values[cell];
	}
	
	/**
	 * Inserts a single integer value into the specified cell, expanding storage as needed while avoiding duplicates.
	 * @param kmer The k-mer key (must match array[cell])
	 * @param v The integer value to insert
	 * @param cell Array index where the k-mer is stored
	 */
	@Override
	protected final void insertValue(final long kmer, final int v, final int cell){
		assert(array[cell]==kmer);
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
	
	/**
	 * Inserts multiple integer values into the specified cell for the given k-mer.
	 * If cell is empty, directly assigns the values array. Otherwise, inserts
	 * each valid value (non-negative) individually to prevent duplicates.
	 *
	 * @param kmer The k-mer key (must match array[cell])
	 * @param vals Array of integer values to insert
	 * @param cell Array index where the k-mer is stored
	 * @param vlen Length of valid values in the array
	 */
	@Override
	protected final void insertValue(final long kmer, final int[] vals, final int cell, final int vlen){
		assert(array[cell]==kmer);
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
	
	/** Returns false as 2D arrays do not support rebalancing operations */
	@Override
	public final boolean canRebalance() {return false;}
	
	/**
	 * Resizes the hash table when load factor exceeds limits or memory pressure occurs.
	 * Allocates new arrays with larger prime-sized capacity, then rehashes all k-mers
	 * from both main array and victim cache. Uses schedule-based or multiplicative
	 * growth depending on configuration. Kills process if maximum size exceeded.
	 */
	@Override
	protected synchronized void resize(){
//		assert(false);
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
			//		System.err.println("Resizing from "+prime+" to "+prime2+"; size="+size);

			prime=prime2;
			sizeLimit=(long)(maxLoadFactor*prime);
		}

//		System.err.println("Resized to "+prime+"; load="+(size*1f/prime));
		long[] oldk=array;
		int[][] oldc=values;
		KmerNode[] oldv=victims.array;
		array=allocLong1D(prime+extra);
		Arrays.fill(array, NOT_PRESENT);
		values=allocInt2D(prime+extra);
		ArrayList<KmerNode> list=new ArrayList<KmerNode>((int)(victims.size)); //Can fail if more than Integer.MAX_VALUE
		for(int i=0; i<oldv.length; i++){
			if(oldv[i]!=null){oldv[i].traverseInfix(list);}
		}
		Arrays.fill(oldv, null);
		victims.size=0;
		size=0;
		
		final int[] singleton=new int[] {NOT_PRESENT};
		
		for(int i=0; i<oldk.length; i++){
			if(oldk[i]>NOT_PRESENT){
//				assert(!contains(oldk[i]));
				set(oldk[i], oldc[i], -1);
//				assert(contains(oldk[i]));
//				assert(Tools.equals(getValues(oldk[i], singleton), oldc[i]));
			}
		}
		
		for(KmerNode n : list){
			if(n.pivot>NOT_PRESENT){
//				assert(!contains(n.pivot));
				set(n.pivot, n.values(singleton), n.numValues());
//				assert(contains(n.pivot));
//				assert(Tools.equals(getValues(n.pivot, singleton), n.values(singleton)));
			}
		}
		
		assert(oldSize+oldVSize==size+victims.size) : oldSize+", "+oldVSize+" -> "+size+", "+victims.size;
	}
	
	/**
	 * Deprecated rebalance method; unimplemented for 2D value storage and always throws.
	 * @throws RuntimeException Always thrown
	 * @deprecated Not supported for HashArray2D
	 */
	@Deprecated
	@Override
	public void rebalance(){
		throw new RuntimeException("Unimplemented.");
	}
	
	/**
	 * Deprecated rebalance operation - not implemented for 2D arrays.
	 * 2D hash arrays do not support rebalancing due to complexity of value management.
	 * @throws RuntimeException Always thrown as operation is unimplemented
	 * @deprecated 2D arrays do not support rebalancing
	 */
	@Deprecated
	@Override
	public long regenerate(final int limit){
		assert(false) : "This is not tested or intended for use.";
		long sum=0;
		assert(owners==null) : "Clear ownership before regeneration.";
		for(int pos=0; pos<values.length; pos++){
			final long key=array[pos];
			if(key>=0){
				final int[] value=values[pos];
				values[pos]=null;
				array[pos]=NOT_PRESENT;
				size--;
				if(value!=null){
					assert(value[0]>0);
					set(key, value, -1);
				}else{
					sum++;
				}
			}
		}
		
		ArrayList<KmerNode> nodes=victims.toList();
		victims.clear();
		for(KmerNode node : nodes){
			set(node.pivot, node.values(null), node.numValues());//TODO: Probably unsafe or unwise.  Should test for singletons, etc.
		}
		
		return sum;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int[][] values;
	

	
}
