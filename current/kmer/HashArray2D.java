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
 * Main cells retain packed value arrays and scan their full used prefixes for
 * duplicates. Overflow nodes have a distinct bounded duplicate scan; this class
 * does not establish a global uniqueness guarantee across both storage paths.
 * Counting increments are unsupported. Shared mutations require caller coordination.
 *
 * @author Brian Bushnell
 * @date Nov 7, 2014
 */
public final class HashArray2D extends HashArray{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Allocates primary value slots using the borrowed resize schedule and a 2D victim forest. */
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
	 * @deprecated Use set to store values; counting is unsupported
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
	 * @deprecated Use set to store values; counting is unsupported
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
	protected final int readCellValue(int cell){
		int[] set=values[cell];
		return set==null ? 0 : set[0];
	}

	/**
	 * Returns the live value array stored in the cell, including unused negative padding.
	 * May return null; callers must not mutate retained storage.
	 * @param cell Array index
	 * @param singleton Ignored placeholder to match superclass signature
	 * @return Value array or null
	 */
	@Override
	protected final int[] readCellValues(int cell, int[] singleton){
		return values[cell];
	}

	/**
	 * Inserts a nonnegative value into the specified cell, expanding storage as needed.
	 * Scans the full used prefix for an existing value before insertion.
	 * @param kmer The k-mer key (must match array[cell])
	 * @param v Nonnegative value; negative values are reserved for unused slots
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

		//Comprehension (verified invariant): a cell's value set is always PACKED — real values (scaffold IDs, >=0) in a
		//prefix, then NOT_PRESENT(-1) padding, no gaps. Every insert either fills the first negative slot or appends at
		//oldSize and fills the rest with NOT_PRESENT, so the prefix stays contiguous. Hence scanning for v (dedup) and
		//stopping at the first negative slot is complete: v cannot appear after a negative, so a negative means "v absent".
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
	 * If cell is empty, retains the values array without copying or deduplicating it.
	 * The caller must supply a packed nonnegative prefix and not mutate retained storage.
	 * Otherwise, scans until the first
	 * negative value and inserts each preceding value individually.
	 *
	 * @param kmer The k-mer key (must match array[cell])
	 * @param vals Array of integer values to insert
	 * @param cell Array index where the k-mer is stored
	 * @param vlen Logical value count supplied by the caller; this implementation does not inspect it
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
	public final boolean canRebalance(){return false;}

	/**
	 * Grows according to the schedule or legacy load heuristics, rehashing primary and victim keys.
	 * Primary reinsertion retains value arrays directly; only overflow insertion scans
	 * their used prefixes (KMR001). Unscheduled growth can return without resizing;
	 * exhausted scheduled growth invokes memory-failure handling. This is not triggered
	 * by a memory-pressure monitor. Callers exclude ordinary access and ownership use;
	 * the resize monitor does not protect other methods or remap owner markers.
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
			final long key=oldk[i];
			if(key>NOT_PRESENT){
//				assert(!contains(oldk[i]));
				final int[] vals=oldc[i];
				final int cell=findKmerOrEmpty(key);
				if(cell==HASH_COLLISION){
					//KMR001: only a new victim node needs the used prefix length; primary cells retain the array directly.
					int vlen=0;
					if(vals!=null){while(vlen<vals.length && vals[vlen]>=0){vlen++;}}
					victims.set(key, vals, vlen);
				}else{
					assert(array[cell]==NOT_PRESENT) : "Old primary keys are unique during resize; an occupied target would double-count size: key="+
						key+", cell="+cell+", existing="+array[cell];
					array[cell]=key;
					values[cell]=vals;
					size++;
				}
				//Retain set()'s growth check, including schedules that require another resize during reinsertion.
				if(autoResize && size+victims.size>sizeLimit){resize();}
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
	 * Unsupported experimental regeneration, rejected immediately when assertions are enabled.
	 * With assertions disabled, reinserts nonnull primary value arrays and all victims;
	 * it ignores limit and counts only removed null-valued primary entries.
	 * This retained code is not an implemented count-threshold filter.
	 * @param limit Ignored by this experimental implementation
	 * @return Removed null-valued primary entries if the unsupported path completes
	 * @throws AssertionError Always with assertions enabled
	 * @deprecated Not tested or intended for use
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

	/** Primary value arrays aligned with key cells; references can be borrowed from callers. */
	private int[][] values;



}
