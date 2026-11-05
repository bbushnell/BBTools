package kmer;

import java.util.ArrayList;
import java.util.Arrays;

import shared.KillSwitch;
import shared.Primes;
import shared.Tools;
import structures.SuperLongList;

/**
 * A 1-dimensional hash table for storing k-mers with single integer counts.
 * Uses a flat array structure with linear probing and a victim cache for overflow.
 * Optimized for memory efficiency by storing counts in a separate int[] array
 * rather than multidimensional arrays used by other hash table implementations.
 * Main-array counts saturate at Integer.MAX_VALUE on positive overflow. The
 * victim-node path now also saturates positive additions (KMR002).
 * Shared mutation requires external locking; exposed arrays are borrowed storage.
 *
 * @author Brian Bushnell
 * @date Oct 25, 2013
 */
public final class HashArray1D extends HashArray{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates aligned key/count arrays using the borrowed resize schedule. */
	public HashArray1D(int[] schedule_, long coreMask_){
		super(schedule_, coreMask_, false);
		values=allocInt1D(prime+extra);
	}

	/** Creates aligned key/count arrays without a schedule; initialSize is used directly. */
	public HashArray1D(int initialSize, long coreMask, boolean autoResize_){
		super(initialSize, coreMask, autoResize_, false);
		values=allocInt1D(prime+extra);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Adds a positive increment, returning the resulting count. Existing nonnegative
	 * primary counts saturate on overflow; new entries retain incr directly. */
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
				//Comprehension [kmer/HashArray1D#001 RESOLVED]: return incr, not 1 — increment() is contracted to
				//"Returns count", and a NEW entry's count IS incr. Now consistent with every sibling (HashForest.increment
				//-> n.value()==incr; KmerTable.increment -> n.value==incr; HashArrayHybrid.increment:74 -> return incr) and
				//with this class's own victims path below (return x). Was previously `return 1`, correct only for incr=1
				//(every current caller passes incr=1, so this is behavior-neutral today and hardens the incr>1 path).
				return incr;
			}
		}
		int x=victims.increment(kmer, incr);
		if(autoResize && size+victims.size>sizeLimit){resize();}
		return x;
	}

	/** Adds a positive increment, returning one for a new key and zero for an existing key. */
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
		int x=victims.incrementAndReturnNumCreated(kmer, incr);
		if(autoResize && size+victims.size>sizeLimit){resize();}
		return x;
	}

	/**
	 * Appends positive primary counts, then delegates victim-node counts to the list.
	 * Primary scanning uses count cells rather than testing key occupancy.
	 * @param sll List to receive all count values for histogram generation
	 */
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

	/**
	 * Reads the count value from a specific array cell.
	 * @param cell Array index to read from
	 * @return Count value at the specified cell
	 */
	@Override
	public final int readCellValue(int cell){
		return values[cell];
	}

	/**
	 * Reads count value from a cell into a provided array.
	 * For 1D arrays, only the first element is populated.
	 *
	 * @param cell Array index to read from
	 * @param singleton Single-element array to receive the value
	 * @return The singleton array with value at index 0
	 */
	@Override
	protected final int[] readCellValues(int cell, int[] singleton){
		singleton[0]=values[cell];
		return singleton;
	}

	/**
	 * Inserts a single integer value at the specified cell.
	 * Assumes the k-mer is already stored in the array at the given cell.
	 *
	 * @param kmer K-mer that should already be stored at this cell
	 * @param v Value to store
	 * @param cell Array index where the value should be stored
	 */
	@Override
	protected final void insertValue(long kmer, int v, int cell){
		assert(array[cell]==kmer);
		values[cell]=v;
	}

	/**
	 * Inserts value from an array into the specified cell.
	 * Requires a one-element array with assertions enabled; vlen is not consulted.
	 *
	 * @param kmer K-mer that should already be stored at this cell
	 * @param vals Array containing values to insert (only index 0 used)
	 * @param cell Array index where the value should be stored
	 * @param vlen Length of valid data in vals array (ignored, assumes 1)
	 */
	@Override
	protected final void insertValue(long kmer, int[] vals, int cell, int vlen){
		assert(array[cell]==kmer);
		assert(vals.length==1);
		values[cell]=vals[0];
	}

	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Indicates whether this hash table supports rebalancing operations.
	 * 1D hash arrays do not support rebalancing.
	 * @return Always false for HashArray1D
	 */
	@Override
	public final boolean canRebalance(){return false;}

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

	/**
	 * Resizes the hash table when load factor exceeds threshold.
	 * On growth, replaces arrays and rehashes all existing keys/counts, including victims.
	 * Unscheduled growth may return without resizing; exhausted scheduled growth
	 * invokes memory-failure handling. Ownership markers are not remapped.
	 * Uses either scheduled sizes or load-factor based calculation for new size.
	 * Supports two-pass resizing where high-count entries are moved first.
	 * Its monitor excludes other resize calls, not ordinary table access.
	 */
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

	/** Rebalancing operation is not supported for 1D hash arrays.
	 * @throws RuntimeException Always thrown as operation is unimplemented */
	@Deprecated
	@Override
	public void rebalance(){
		throw new RuntimeException("Unimplemented.");
	}

	/**
	 * Removes entries with counts at or below the specified limit.
	 * Entries with higher counts are preserved. Returns count of removed entries.
	 * Requires callers to clear ownership first and exclude concurrent access.
	 * Surviving keys are reinserted to close probe gaps; this method does not clear owners.
	 *
	 * @param limit Maximum count value to remove (inclusive)
	 * @return Number of entries that were removed
	 */
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

	/** Returns string representation of the k-mer array contents.
	 * @return String showing array values */
	@Override
	public String toString(){
		return Arrays.toString(array);
	}

	/** Creates a primary-then-victim walker that skips nonpositive primary counts.
	 * Keep the table stable during traversal; victim nodes are captured by reference. */
	public Walker walk(){
		return new Walker1D();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Count cells aligned with the key array, including its probe tail. */
	private int[] values;

	/** Returns live count storage, replaced during resizing; caller must preserve key alignment. */
	public int[] values(){return values;}

	/*--------------------------------------------------------------*/
	/*----------------            Walker            ----------------*/
	/*--------------------------------------------------------------*/

	/** Traverses live primary storage and a captured list of live victim nodes. */
	public class Walker1D extends Walker{

		/** Captures victim node references; primary arrays remain those of the enclosing table. */
		Walker1D(){
			ha=HashArray1D.this;
			victims=ha.victims().toList();
		}

		/**
		 * Advances to the next k-mer entry and fills kmer/value fields.
		 * Iterates through main array first, then victim cache entries.
		 * Skips nonpositive primary counts: fresh empty cells contain zero and
		 * regenerate() leaves -1 in vacated cells. Victim traversal is unchanged.
		 * No snapshot or synchronization protects table mutation.
		 * @return true for an emitted entry, false at exhaustion
		 */
		public boolean next(){
			//KMR004 (kmer/HashArray1D#002): both zero-initialized and regenerated empty cells must be skipped.
			while(i<values.length && values[i]<=NOT_XPRESENT){i++;}
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

		/** Returns the last emitted key, or -1 after exhaustion; undefined before next(). */
		public long kmer(){return kmer;}
		/** Returns the last emitted value, or zero after exhaustion; undefined before next(). */
		public int value(){return value;}

		/** Reference to the hash array being walked */
		private HashArray1D ha;
		/** Captured overflow node references; node values are read when visited. */
		private ArrayList<KmerNode> victims;

		/** Last emitted key, replaced with -1 at exhaustion. */
		private long kmer;
		/** Last emitted count, replaced with NOT_XPRESENT at exhaustion. */
		private int value;

		/** Next primary-array slot to examine. */
		private int i=0;
		/** Next captured victim node to visit. */
		private int i2=0;
	}

	/** Zero-initialized count cells and the exhausted walker's value. Regeneration uses -1. */
	private static final int NOT_XPRESENT=0;

	/** Estimates array payload and victim-node bytes, omitting JVM object/array headers.
	 * Uses long arithmetic for array payloads; this remains an estimate, not a measured heap requirement. */
	public long calcMem(){
		long mem=0;
		//KMR020: widen before multiplication so large array payloads cannot wrap as int products.
		mem+=(array.length*8L);
		mem+=(values.length*4L);
		mem+=(owners==null ? 0 : owners.length()*4L);
		for(KmerNode kn : victims.array){
			mem+=8;
			if(kn!=null){mem+=kn.calcMem();}
		}
		// TODO Auto-generated method stub
		return mem;
	}


}
