package bloom;

import java.lang.Thread.State;
import java.util.ArrayList;
import java.util.Arrays;
import shared.Random;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import shared.KillSwitch;
import shared.Primes;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;


/**
 * Atomic prime-sized Bloom array with multi-lane hashing and optional prefilter gating.
 * Occupancy describes physical main-table cells, not counts returned from the prefilter.
 * No-argument occupancy methods cache their first scan; updates do not invalidate it.
 * The first shutdown refreshes that cache for a completed table. Callers must finish
 * writers before finalization and safely publish the table; shutdown neither stops
 * writers nor prevents later updates. Later mutation can leave cached statistics stale.
 * Positive-depth occupancy overloads rescan, but are not atomic snapshots during writes.
 * The legacy decrement operation has an unresolved deletion/prefilter contract.
 *
 * @author Brian Bushnell
 * @date Aug 17, 2012
 */
public final class KCountArray7MTA extends KCountArray{

	/** Serialization identifier. */
	private static final long serialVersionUID=568264681638739631L;
	
	/** Demonstrates basic increment and read operations for selected keys. */
	public static void main(String[] args){
		long cells=Long.parseLong(args[0]);
		int bits=Integer.parseInt(args[1]);
		int hashes=Integer.parseInt(args[2]);
		
		verbose=false;
		
		KCountArray7MTA kca=new KCountArray7MTA(cells, bits, hashes, null, 0);
		
		System.out.println(kca.read(0));
		kca.increment(0);
		System.out.println(kca.read(0));
		kca.increment(0);
		System.out.println(kca.read(0));
		System.out.println();
		
		System.out.println(kca.read(1));
		kca.increment(1);
		System.out.println(kca.read(1));
		kca.increment(1);
		System.out.println(kca.read(1));
		System.out.println();
		
		System.out.println(kca.read(100));
		kca.increment(100);
		System.out.println(kca.read(100));
		kca.increment(100);
		System.out.println(kca.read(100));
		kca.increment(100);
		System.out.println(kca.read(100));
		System.out.println();
		

		System.out.println(kca.read(150));
		kca.increment(150);
		System.out.println(kca.read(150));
		System.out.println();
		
	}
	
	/**
	 * Constructs a prime-partitioned atomic array with optional prefilter gating.
	 *
	 * @param cells_ target number of cells, rounded down per partition to a prime
	 * @param bits_ bits per packed count cell
	 * @param hashes_ number of Bloom lanes
	 * @param prefilter_ optional first-stage count array
	 * @param prefilterLimit_ count below which the prefilter result is returned directly
	 */
	public KCountArray7MTA(long cells_, int bits_, int hashes_, KCountArray prefilter_, int prefilterLimit_){
		super(getPrimeCells(cells_, bits_), bits_, getDesiredArrays(cells_, bits_));
//		verbose=false;
//		System.out.println(cells);
		cellsPerArray=cells/numArrays;
		wordsPerArray=(int)((cellsPerArray%cellsPerWord)==0 ? (cellsPerArray/cellsPerWord) : (cellsPerArray/cellsPerWord+1));
		cellMod=cellsPerArray;
		hashes=hashes_;
		prefilter=prefilter_;
		prefilterLimit=(prefilter==null ? 0 : Tools.min(prefilter.maxValue, prefilterLimit_));
//		System.out.println("cells="+cells+", words="+words+", wordsPerArray="+wordsPerArray+", numArrays="+numArrays+", hashes="+hashes);
		
		matrix=allocMatrix(numArrays, wordsPerArray);
		
		useLocks=(cellBits>1 && hashes>1) && (LOCKED_INCREMENT || !SET_LOCKED_INCREMENT);
		if(useLocks){
			locks=new Lock[NUM_LOCKS];
			for(int i=0; i<NUM_LOCKS; i++){locks[i]=new ReentrantLock();}
		}else{locks=null;}
				
//		matrix=new AtomicIntegerArray[numArrays];
//		for(int i=0; i<matrix.length; i++){
//			matrix[i]=new AtomicIntegerArray(wordsPerArray);
//		}
		
		assert(hashes>0 && hashes<=hashMasks.length);
	}
	
	/** Returns the partition count needed to keep each AtomicIntegerArray addressable. */
	private static int getDesiredArrays(long desiredCells, int bits){
		
		long words=Tools.max((desiredCells*bits+31)/32, minArrays);
		int arrays=minArrays;
		while(words/arrays>=Integer.MAX_VALUE){
			arrays*=2;
		}
//		assert(false) : arrays;
		return arrays;
//		return Tools.max(arrays, Data.LOGICAL_PROCESSORS*4);
	}
	
	/** Returns the largest prime-per-partition cell count not exceeding the target. */
	private static long getPrimeCells(long desiredCells, int bits){
		
		int arrays=getDesiredArrays(desiredCells, bits);
		
		long x=(desiredCells+arrays-1)/arrays;
		long x2=Primes.primeAtMost(x);
		return x2*arrays;
	}
	
	@Override
	public final int read(final long rawKey){
		if(verbose){System.err.println("Reading raw key "+rawKey);}
		if(prefilter!=null){
			int pre=prefilter.read(rawKey);
			if(pre<prefilterLimit){return pre;}
		}
		long key2=hash(rawKey, 0);
		int min=readHashed(key2);
		for(int i=1; i<hashes && min>0; i++){
			if(verbose){System.err.println("Reading. i="+i+", key2="+key2);}
			key2=Long.rotateRight(key2, hashBits);
			key2=hash(key2, i);
			if(verbose){System.err.println("Rot/hash. i="+i+", key2="+key2);}
			min=min(min, readHashed(key2));
		}
		return min;
	}

	@Override
	public final int read(final long rawKey, final long rawKey2){
		if(verbose){System.err.println("Reading raw keys "+rawKey+", "+rawKey2);}
		if(prefilter!=null){
			int pre=prefilter.read(rawKey, rawKey2);
			if(pre<prefilterLimit){return pre;}
		}
		long key=hash(rawKey, 0);
		int min=readHashed(key);
		for(int i=1; i<hashes && min>0; i++){
			key=nextHash(key, rawKey2, i);
			min=min(min, readHashed(key));
		}
		return min;
	}
	
	@Override
	public final int read(final long[] rawKeys){
		if(verbose){System.err.println("Reading raw key "+Arrays.toString(rawKeys));}
		if(prefilter!=null){
			int pre=prefilter.read(rawKeys);
			if(pre<prefilterLimit){return pre;}
		}
		long key2=hash(rawKeys[0], (int)(1+(rawKeys[0])%5));
		int min=maxValue;
		for(int i=0; i<hashes; i++){
			for(int keynum=0; keynum<rawKeys.length; keynum++){
				if(verbose){System.err.println("Reading. i="+i+", key2="+key2);}
				key2=hash(key2^rawKeys[keynum], i);
				if(verbose){System.err.println("Rot/hash. i="+i+", key2="+key2);}
			}
			min=min(min, readHashed(key2));
			key2=Long.rotateRight(key2, hashBits);
		}
		return min;
	}
	
	/** Returns the four left-neighbor count sum, clamped at Integer.MAX_VALUE. */
	@Override
	public final int readLeft(final long key, final int k, boolean makeCanonical){
		assert(k<=32);
		final long key2=key>>>2;
		final int shift=2*(k-1);
		final long akey=key2|(0L<<shift);
		final long ckey=key2|(1L<<shift);
		final long gkey=key2|(2L<<shift);
		final long tkey=key2|(3L<<shift);
		final int a=read(makeCanonical ? makeCanonical2(akey, k) : akey);
		final int c=read(makeCanonical ? makeCanonical2(ckey, k) : ckey);
		final int g=read(makeCanonical ? makeCanonical2(gkey, k) : gkey);
		final int t=read(makeCanonical ? makeCanonical2(tkey, k) : tkey);
		final long sum=(long)a+c+g+t;
		return sum>Integer.MAX_VALUE ? Integer.MAX_VALUE : (int)sum;
	}
	
	/** Returns the four right-neighbor count sum, clamped at Integer.MAX_VALUE. */
	@Override
	public final int readRight(final long key, final int k, boolean makeCanonical){
		assert(k<=32);
		final long mask=(k>=32 ? -1L : ~((-1L)<<(2*k)));
		final long key2=(key<<2)&mask;
		final long akey=key2|0L;
		final long ckey=key2|1L;
		final long gkey=key2|2L;
		final long tkey=key2|3L;
		final int a=read(makeCanonical ? makeCanonical2(akey, k) : akey);
		final int c=read(makeCanonical ? makeCanonical2(ckey, k) : ckey);
		final int g=read(makeCanonical ? makeCanonical2(gkey, k) : gkey);
		final int t=read(makeCanonical ? makeCanonical2(tkey, k) : tkey);
		final long sum=(long)a+c+g+t;
		return sum>Integer.MAX_VALUE ? Integer.MAX_VALUE : (int)sum;
	}
	
	@Override
	public final int[] readAllLeft(final long key, final int k, boolean makeCanonical, int[] rvec){
		assert(k<=32);
		if(rvec==null){rvec=KillSwitch.allocInt1D(4);}
		final long key2=key>>>2;
		final int shift=2*(k-1);
		final long akey=key2|(0L<<shift);
		final long ckey=key2|(1L<<shift);
		final long gkey=key2|(2L<<shift);
		final long tkey=key2|(3L<<shift);
		rvec[0]=read(makeCanonical ? makeCanonical2(akey, k) : akey);
		rvec[1]=read(makeCanonical ? makeCanonical2(ckey, k) : ckey);
		rvec[2]=read(makeCanonical ? makeCanonical2(gkey, k) : gkey);
		rvec[3]=read(makeCanonical ? makeCanonical2(tkey, k) : tkey);
		return rvec;
	}
	
	@Override
	public final int[] readAllRight(final long key, final int k, boolean makeCanonical, int[] rvec){
		assert(k<=32);
		final long mask=(k>=32 ? -1L : ~((-1L)<<(2*k)));
		final long key2=(key<<2)&mask;
		final long akey=key2|0L;
		final long ckey=key2|1L;
		final long gkey=key2|2L;
		final long tkey=key2|3L;
		rvec[0]=read(makeCanonical ? makeCanonical2(akey, k) : akey);
		rvec[1]=read(makeCanonical ? makeCanonical2(ckey, k) : ckey);
		rvec[2]=read(makeCanonical ? makeCanonical2(gkey, k) : gkey);
		rvec[3]=read(makeCanonical ? makeCanonical2(tkey, k) : tkey);
		return rvec;
	}
	
	/** Returns the packed count for an already lane-hashed key. */
	private final int readHashed(long key){
		if(verbose){System.err.print("Reading hashed key "+key);}
//		System.out.println("key="+key);
		int arrayNum=(int)(key&arrayMask);
		key=(key>>>arrayBits)%(cellMod);
//		key=(key>>>(arrayBits+1))%(cellMod);
//		System.out.println("array="+arrayNum);
//		System.out.println("key2="+key);
		AtomicIntegerArray array=matrix[arrayNum];
		int index=(int)(key>>>indexShift);
//		assert(false) : indexShift;
//		System.out.println("index="+index);
		int word=array.get(index);
//		System.out.println("word="+Integer.toHexString(word));
		assert(word>>>(cellBits*key) == word>>>(cellBits*(key&cellMask)));
//		int cellShift=(int)(cellBits*(key&cellMask));
		int cellShift=(int)(cellBits*key);
		if(verbose){System.err.println(", array="+arrayNum+", index="+index+", cellShift="+(cellShift%32)+", value="+((int)((word>>>cellShift)&valueMask)));}
//		System.out.println("cellShift="+cellShift);
		return (int)((word>>>cellShift)&valueMask);
	}
	
	@Override
	public final void write(final long key, int value){
		throw new RuntimeException("Not allowed for this class.");
	}
	
	@Override
	public final void increment(long[] keys){
//		assert(false) : "This method is not really needed.";
		for(int i=0; i<keys.length; i++){
			increment(keys[i]);
//			keys[i]=hash(keys[i], 0);
//			incrementPartiallyHashed(hash(keys[i], 0));
		}
	}
	
	@Override
	public final void increment(final long rawKey, final int amt){
//		if(verbose){System.err.println("\n*** Incrementing raw key "+rawKey+" ***");}
//		
//		if(prefilter!=null){
//			int x=prefilter.read(rawKey);
//			if(x<prefilterLimit){return;/* x;*/}
//		}
//		
//		if(useLocks){incrementLocked(rawKey, amt);}
//		else{incrementUnlocked(rawKey, amt);}
		incrementAndReturnUnincremented(rawKey, amt);
	}

	@Override
	public final void increment(final long rawKey, final long rawKey2, final int amt){
		incrementAndReturnUnincremented(rawKey, rawKey2, amt);
	}
	
	/*
	private void incrementUnlocked(long rawKey, int amt){
		long key2=rawKey;
//		int min=maxValue;
		for(int i=0; i<hashes; i++){
			key2=hash(key2, i);
			if(verbose){System.err.println("key2="+key2+", value="+readHashed(key2));}
//			min=min(min, incrementHashedLocal(key2, amt));
			incrementHashedLocal(key2, amt);
			key2=Long.rotateRight(key2, hashBits);
		}
//		return min;
	}
	
	private void incrementLocked(long rawKey, int amt){
		assert(amt>0) : "Don't increment by zero, and negatives should use decrement amt="+amt;
		final Lock lock=getLock(rawKey);
		lock.lock();
		final int min=read(rawKey);
		final long newMinL=Tools.min(min+amt, maxValue);
		final int newMin=(int)newMinL;
		assert(newMin==newMinL && newMin>0) : newMin+", "+newMinL+", "+min+", "+amt;
		if(newMin<=min){
			assert(newMin==min) : min;
			lock.unlock();
			return;
		}
		long key2=rawKey;
		for(int i=0; i<hashes; i++){
			key2=hash(key2, i);
			if(verbose){System.err.println("key2="+key2+", value="+readHashed(key2));}
			//incrementHashedLocal_ifAtMost(key2, min);
			incrementHashedLocal_toAtLeast(key2, newMin);
			key2=Long.rotateRight(key2, hashBits);
		}
		lock.unlock();
//		return max(min, newMin);
	}
	*/
	
	/**
	 * Legacy primitive-lane subtraction; its public deletion contract remains unresolved.
	 * TODO: Resolve deletion semantics and the historical prefilter assertion before relying on this operation.
	 * When enabled, the existing assertion requires a nonnull prefilter, although this method neither
	 * consults nor changes it. That assertion is unresolved, not a validated precondition.
	 * For positive amounts the body subtracts from main-table lanes with per-word CAS,
	 * clamping each cell at zero. It does not take the increment stripe lock or refresh
	 * cached occupancy. Reads below the prefilter threshold may still return unchanged
	 * prefilter counts. Conservative increments, saturation and collisions prevent a
	 * general inverse-of-insertion guarantee; concurrent multi-lane deletion is not atomic.
	 * @param rawKey Primitive key whose main-table lanes are visited
	 * @param amt Subtraction amount; valid public deletion semantics/range remain unresolved
	 */
	@Override
	public final void decrement(final long rawKey, int amt){
		if(verbose){System.err.println("\n*** Decrementing raw key "+rawKey+" ***");}
		
		assert(prefilter!=null); //Why?
		
		long key2=rawKey;
//		int min=maxValue;
		for(int i=0; i<hashes; i++){
			key2=hash(key2, i);
			if(verbose){System.err.println("key2="+key2+", value="+readHashed(key2));}

			decrementHashedLocal(key2, amt);
//			min=min(min, decrementHashedLocal(key2, amt));
			key2=Long.rotateRight(key2, hashBits);
		}
//		return min;
	}
	
	@Override
	public int incrementAndReturnUnincremented(final long rawKey, final int incr){

		if(verbose){System.err.println("\n*** Incrementing raw key "+rawKey+" ***");}
		
		if(prefilter!=null){
			int x=prefilter.read(rawKey);
			if(x<prefilterLimit){return x;}
		}
		
		if(useLocks){return incrementAndReturnUnincrementedLocked(rawKey, incr);}
		else{return incrementAndReturnUnincrementedUnlocked(rawKey, incr);}
	}

	@Override
	public int incrementAndReturnUnincremented(final long rawKey, final long rawKey2, final int incr){
		if(verbose){System.err.println("\n*** Incrementing raw keys "+rawKey+", "+rawKey2+" ***");}
		if(prefilter!=null){
			int x=prefilter.read(rawKey, rawKey2);
			if(x<prefilterLimit){return x;}
		}
		if(useLocks){return incrementAndReturnUnincrementedLocked(rawKey, rawKey2, incr);}
		return incrementAndReturnUnincrementedUnlocked(rawKey, rawKey2, incr);
	}
	
	/** Raises every primitive-key lane to the locked saturated target and returns the previous minimum. */
	private int incrementAndReturnUnincrementedLocked(final long rawKey, final int incr){
		assert(incr>0) : incr;
		final Lock lock=getLock(rawKey);
		lock.lock();
		final int min=read(rawKey);
		if(min>=maxValue){
			lock.unlock();
			return min;
		}
		final int newMin=(int)Tools.min(min+(long)incr, maxValue);
		long key2=rawKey;
		int value=maxValue;
		for(int i=0; i<hashes; i++){
			key2=hash(key2, i);
			if(verbose){System.err.println("key2="+key2+", value="+readHashed(key2));}
			int x=incrementHashedLocalAndReturnUnincremented_toAtLeast(key2, newMin);
			value=min(value, x);
			key2=Long.rotateRight(key2, hashBits);
		}
		lock.unlock();
		return value;
	}

	/** Raises every dual-hash lane to the locked saturated target and returns the previous minimum. */
	private int incrementAndReturnUnincrementedLocked(final long rawKey, final long rawKey2, final int incr){
		assert(incr>0) : incr;
		final Lock lock=getLock(rawKey);
		lock.lock();
		final int min=read(rawKey, rawKey2);
		if(min>=maxValue){
			lock.unlock();
			return min;
		}
		final int newMin=(int)Tools.min(min+(long)incr, maxValue);
		long key=hash(rawKey, 0);
		int value=incrementHashedLocalAndReturnUnincremented_toAtLeast(key, newMin);
		for(int i=1; i<hashes; i++){
			key=nextHash(key, rawKey2, i);
			value=min(value, incrementHashedLocalAndReturnUnincremented_toAtLeast(key, newMin));
		}
		lock.unlock();
		return value;
	}
	
	/** Increments every primitive-key lane with CAS and returns the previous minimum. */
	private int incrementAndReturnUnincrementedUnlocked(final long rawKey, final int incr){
		long key2=rawKey;
		int value=maxValue;
		for(int i=0; i<hashes; i++){
			key2=hash(key2, i);
			if(verbose){System.err.println("key2="+key2+", value="+readHashed(key2));}
			int x=incrementHashedLocalAndReturnUnincremented(key2, incr);
			value=min(value, x);
			key2=Long.rotateRight(key2, hashBits);
		}
		return value;
	}

	/** Increments every dual-hash lane with CAS and returns the previous minimum. */
	private int incrementAndReturnUnincrementedUnlocked(final long rawKey, final long rawKey2, final int incr){
		long key=hash(rawKey, 0);
		int value=incrementHashedLocalAndReturnUnincremented(key, incr);
		for(int i=1; i<hashes; i++){
			key=nextHash(key, rawKey2, i);
			value=min(value, incrementHashedLocalAndReturnUnincremented(key, incr));
		}
		return value;
	}

	/** Preserves the original lane recurrence when rawKey2 is zero. */
	private long nextHash(final long previous, final long rawKey2, final int row){
		final long seed=Long.rotateRight(previous, hashBits)^Long.rotateLeft(rawKey2, 25*row);
		return hash(seed, row);
	}

	/** Replays row 0 through {@code row} with raw primary and secondary inputs for
	 * diagnostics. Valid rows satisfy 0&lt;={@code row}&lt;{@code hashes}; the returned
	 * value is the mixed long before cell mapping and no table state is changed. */
	public final long hashForLane(final long rawKey, final long rawKey2, final int row){
		assert(row>=0 && row<hashes) : row;
		long key=hash(rawKey, 0);
		for(int i=1; i<=row; i++){key=nextHash(key, rawKey2, i);}
		return key;
	}
	
	@Override
	public int incrementAndReturnUnincremented(final long[] rawKeys, final int incr){

		if(verbose){System.err.println("\n*** Incrementing raw keys "+Arrays.toString(rawKeys)+" ***");}
		
		if(prefilter!=null){
			int x=prefilter.read(rawKeys);
			if(x<prefilterLimit){return x;}
		}
		
		long key2=hash(rawKeys[0], (int)(1+(rawKeys[0])%5));
		int value=maxValue;
		for(int i=0; i<hashes; i++){
			for(int keynum=0; keynum<rawKeys.length; keynum++){
				key2=hash(key2^rawKeys[keynum], i);
				if(verbose){System.err.println("key2="+key2+", value="+readHashed(key2));}
			}
			int x=incrementHashedLocalAndReturnUnincremented(key2, incr);
			value=min(value, x);
			key2=Long.rotateRight(key2, hashBits);
		}
//		assert(value+1==read(rawKeys) || value==maxValue) : value+", "+read(rawKeys);
		return value;
	}
	
	@Override
	public long[] transformToFrequency(){
		return transformToFrequency(matrix);
	}
	
	@Override
	public ByteBuilder toContentsString(){
		ByteBuilder sb=new ByteBuilder();
		sb.append('[');
		String comma="";
		for(AtomicIntegerArray array : matrix){
			for(int i=0; i<array.length(); i++){
				int word=array.get(i);
				for(int j=0; j<cellsPerWord; j++){
					int x=word&valueMask;
					sb.append(comma);
					sb.append(x);
					word>>>=cellBits;
					comma=", ";
				}
			}
		}
		sb.append(']');
		return sb;
	}
	
	/** Returns cached positive main-cell occupancy divided by logical cell capacity.
	 * Shares the lifecycle and post-mutation limitations of {@link #cellsUsed()}. */
	@Override
	public double usedFraction(){return cellsUsed()/(double)cells;}
	
	/** Rescans main cells at or above a positive depth and divides by logical capacity.
	 * Does not consult or refresh the no-argument cache; not a snapshot during writes. */
	@Override
	public double usedFraction(int mindepth){return cellsUsed(mindepth)/(double)cells;}
	
	/**
	 * Rescans main-table cells at or above the requested positive depth on every call.
	 * Ignores prefilter-visible counts and does not consult or update the no-argument
	 * occupancy cache. Stop writers for a stable count; atomic word reads do not make
	 * the multi-array scan a coherent snapshot during mutation.
	 * @param mindepth Positive inclusive minimum stored count
	 * @return Observed main-table cells meeting the threshold
	 */
	@Override
	public long cellsUsed(int mindepth){
		return cellsUsedMT(mindepth);
	}
	
	/** Rescans cells at or above a positive threshold with one worker per array.
	 * Joins scan workers before returning; callers separately coordinate writers.
	 * Does not update cached occupancy or inspect the prefilter. */
	public long cellsUsedMT(int mindepth){
//		assert(false) : matrix.length;
		ArrayList<CountUsedThread> list=new ArrayList<CountUsedThread>(matrix.length);
		for(AtomicIntegerArray aia : matrix){
			CountUsedThread ctt=new CountUsedThread(aia, mindepth);
			ctt.start();
			list.add(ctt);
		}
		long x=0;
		for(CountUsedThread ctt : list){
			while(ctt.getState()!=State.TERMINATED){
				try{
					ctt.join();
				}catch(InterruptedException e){
					e.printStackTrace();
				}
			}
			x+=ctt.count;
		}
		return x;
	}
	
	/** Counts cells at or above a positive threshold in one AtomicIntegerArray row. */
	private class CountUsedThread extends Thread{
		/** Binds this worker to one row and a positive count threshold. */
		public CountUsedThread(AtomicIntegerArray a_, int mindepth_){
			array=a_;
			mindepth=mindepth_;
		}
		@Override
		public void run(){
			long temp=0;
			if(array!=null){
//				System.out.println("C");
//				assert(false) : Integer.toBinaryString(valueMask);
				if(cellBits==32){
					for(int i=0, max=array.length(); i<max; i++){
						int word=array.get(i);
						if(word!=0){
							int x=word&valueMask;
							if(x>=mindepth){temp++;}
						}
					}
				}else{
					for(int i=0, max=array.length(); i<max; i++){
						//					System.out.println("D: "+Integer.toHexString(word));
						int word=array.get(i);
						while(word!=0){
							int x=word&valueMask;
							//						System.out.println("E: "+x+", "+mindepth);
							if(x>=mindepth){temp++;}
							word=word>>>cellBits;
						}
					}
				}
			}
			count=temp;
		}
		/** Row to scan for cell counting. */
		private final AtomicIntegerArray array;
		/** Positive count threshold for inclusion. */
		private final int mindepth;
		/** Number of cells that met the threshold. */
		public long count;
	}
	
	
	/** Applies one raw diagnostic mix step without rotation or lane replay.
	 * Row 0 performs an extra mix; valid rows follow {@code hashMasks.length}, not
	 * the live {@code hashes} count. The returned value is not a cell index. */
	public final long diagnosticHash(final long key, final int row){return hash(key, row);}

	@Override
	final long hash(long key, int row){
		int cell=(int)((Long.MAX_VALUE&key)%(hashArrayLength-1));
//		int cell=(int)(hashCellMask&(key));
		
		if(row==0){//Doublehash only first time
			key=key^hashMasks[(row+4)&7][cell];
//			key=key^hashMasks[4][cell];
			cell=(int)(hashCellMask&(key>>5));
//			cell=(int)(hashCellMask&(key>>hashBits));
//			cell=(int)((Long.MAX_VALUE&key)%(hashArrayLength-1));
		}

		assert(row>=0 && row<hashMasks.length) : row+", "+hashMasks.length;
		assert(cell>=0 && cell<hashMasks[0].length) : cell+", "+hashMasks[0].length;
		return key^hashMasks[row][cell];
	}
	
	
//	@Override
//	final long hash(long key, int row){
//		long code=key;
//		for(int i=0; i<6; i++){
//			int row2=((row+i)&7);
//			int cell=(int)(key&hashCellMask);
//			code^=hashMasks[row2][cell];
//			key>>=6;
//		}
//		return code;
//	}
	
	/** Returns randomized hash masks for lane-specific key mixing. */
	private static long[][] makeMasks(int rows, int cols){
		
		long seed;
		synchronized(KCountArray7MTA.class){
			seed=counter^SEEDMASK;
			counter+=7;
		}
		
		Timer t=new Timer();
		long[][] r=new long[rows][cols];
		Random randy=Shared.threadLocalRandom(seed);
		for(int i=0; i<r.length; i++){
			fillMasks(r[i], randy);
		}
		t.stop();
		if(t.elapsed>200000000L){System.out.println("Mask-creation time: "+t);}
		return r;
	}
	
	/** Fills a hash-mask row with 16-bit half-balanced masks and unique low/high slots. */
	private static void fillMasks(long[] r, Random randy){
//		for(int i=0; i<r.length; i++){
//			long x=0;
//			while(Long.bitCount(x&0xFFFFFFFF)!=16){
//				x=randy.nextLong();
//			}
//			r[i]=(x&Long.MAX_VALUE);
//		}
		
		final int hlen=(1<<hashBits);
		assert(r.length==hlen);
		int[] count1=new int[hlen];
		int[] count2=new int[hlen];
		final long mask=hlen-1;

		for(int i=0; i<r.length; i++){
			long x=0;
			int y=0;
			int z=0;
			while(Long.bitCount(x&0xFFFFFFFFL)!=16){
				x=randy.nextLong();
				while(Long.bitCount(x&0xFFFFFFFFL)<16){
					x|=(1L<<randy.nextInt(32));
				}
				while(Long.bitCount(x&0xFFFFFFFFL)>16){
					x&=(~(1L<<randy.nextInt(32)));
				}
				while(Long.bitCount(x&0xFFFFFFFF00000000L)<16){
					x|=(1L<<(randy.nextInt(32)+32));
				}
				while(Long.bitCount(x&0xFFFFFFFF00000000L)>16){
					x&=(~(1L<<(randy.nextInt(32)+32)));
				}
				
//				System.out.print(".");
//				y=(((int)(x&mask))^i);
				y=(((int)(x&mask)));
				z=(int)((x>>hashBits)&mask);
				if(count1[y]>0 || count2[z]>0){
					x=0;
				}
			}
//			System.out.println(Long.toBinaryString(x));
			r[i]=(x&Long.MAX_VALUE);
			count1[y]++;
			count2[z]++;
		}
		
	}
	
	
	/** No-op: construction initializes storage; this does not reset finalization or occupancy. */
	@Override
	public void initialize(){}
	
	/**
	 * Refreshes cached positive-cell occupancy on the first call, then marks it finalized.
	 * Later calls do not rescan, even after mutation. The caller must stop/join writers
	 * first; this method does not stop them or make the table immutable. Subsequent
	 * updates can change readable counts while leaving cached occupancy stale; use a
	 * positive-depth scan for an uncached observation, with writers stopped for stability.
	 */
	@Override
	public void shutdown(){
		if(finished){return;}
		synchronized(this){
			if(finished){return;}
			
			cellsUsed=-1;
//			for(int i=0; i<numArrays; i++){
//				cellsUsed+=cellsUsedPersonal.get(i);
//			}
			cellsUsed();
			
			assert(!finished);
			finished=true;
		}
	}
	
	/** Selects the striped lock for a raw primitive key. */
	private final Lock getLock(long rawKey){
		final Lock lock=locks[(int)((rawKey&Long.MAX_VALUE)%(NUM_LOCKS))];
		return lock;
	}
	
	/** Increments one lane's packed cell with CAS and returns the new saturated value. */
	private int incrementHashedLocal(long key, int amt){
		final int num=(int)(key&arrayMask);
		final AtomicIntegerArray array=matrix[num];
		key=(key>>>arrayBits)%(cellMod);
//		key=(key>>>(arrayBits+1))%(cellMod);
		int index=(int)(key>>>indexShift);
		int cellShift=(int)(cellBits*key);
		int value, word, word2;
		do{
			assert(index>=0) : key+", "+cellMod+", "+cellBits+", "+valueMask+", "+arrayMask+", "+index+", "+num;
			word=array.get(index);
			value=((word>>>cellShift)&valueMask);
			value=(int)min(value+(long)amt, maxValue);
			word2=(value<<cellShift)|(word&~((valueMask)<<cellShift));
		}while(word!=word2 && !array.compareAndSet(index, word, word2));
//		if(value==1){cellsUsedPersonal.incrementAndGet(num);}
		return value;
	}
	
//	private int incrementHashedLocal_ifAtMost(long key, final int thresh){
//		final int num=(int)(key&arrayMask);
//		final AtomicIntegerArray array=matrix[num];
//		key=(key>>>arrayBits)%(cellMod);
////		key=(key>>>(arrayBits+1))%(cellMod);
//		int index=(int)(key>>>indexShift);
//		int cellShift=(int)(cellBits*key);
//		int value, word, word2;
//		do{
//			assert(index>=0) : key+", "+cellMod+", "+cellBits+", "+valueMask+", "+arrayMask+", "+index+", "+num;
//			word=array.get(index);
//			value=((word>>>cellShift)&valueMask);
//			if(value>thresh){return value;}//Too high; don't increment
//			value=min(value+1, maxValue);
//			word2=(value<<cellShift)|(word&~((valueMask)<<cellShift));
//		}while(word!=word2 && !array.compareAndSet(index, word, word2));
//		return value;
//	}
	
	/** Raises one lane's packed cell by one when it is below the locked target. */
	private int incrementHashedLocal_toAtLeast(long key, final int newMin){
		assert(newMin<=maxValue);
		final int num=(int)(key&arrayMask);
		final AtomicIntegerArray array=matrix[num];
		key=(key>>>arrayBits)%(cellMod);
//		key=(key>>>(arrayBits+1))%(cellMod);
		int index=(int)(key>>>indexShift);
		int cellShift=(int)(cellBits*key);
		int value, word, word2;
		do{
			assert(index>=0) : key+", "+cellMod+", "+cellBits+", "+valueMask+", "+arrayMask+", "+index+", "+num;
			word=array.get(index);
			value=((word>>>cellShift)&valueMask);
			if(value>=newMin){return value;}//Too high; don't increment
			value=(int)min(value+1L, maxValue);
			word2=(value<<cellShift)|(word&~((valueMask)<<cellShift));
		}while(word!=word2 && !array.compareAndSet(index, word, word2));
		return value;
	}
	
	/** Increments one lane's packed cell with CAS and returns its previous value. */
	private int incrementHashedLocalAndReturnUnincremented(long key, int incr){
		assert(incr>=0);
		final int num=(int)(key&arrayMask);
		final AtomicIntegerArray array=matrix[num];
		key=(key>>>arrayBits)%(cellMod);
//		key=(key>>>(arrayBits+1))%(cellMod);
		int index=(int)(key>>>indexShift);
		int cellShift=(int)(cellBits*key);
		int value, word, word2;
		do{
			word=array.get(index);
			value=((word>>>cellShift)&valueMask);
			int value2=(int)min(value+(long)incr, maxValue);//was min(value+incr,...) (int add): overflows at cellBits=32 when value near 2^31-1, wrapping the maxValue cap negative. Matches the guarded twin incrementHashedLocal. Latent (needs 32-bit cells + ~2^31 count of one kmer).
			word2=(value2<<cellShift)|(word&~((valueMask)<<cellShift));
		}while(word!=word2 && !array.compareAndSet(index, word, word2));
//		if(value==1){cellsUsedPersonal.incrementAndGet(num);}
		return value;
	}
	
//	private int incrementHashedLocalAndReturnUnincremented_ifAtMost(long key, int incr, final int thresh){
//		assert(incr>=0);
//		final int num=(int)(key&arrayMask);
//		final AtomicIntegerArray array=matrix[num];
//		key=(key>>>arrayBits)%(cellMod);
////		key=(key>>>(arrayBits+1))%(cellMod);
//		int index=(int)(key>>>indexShift);
//		int cellShift=(int)(cellBits*key);
//		int value, word, word2;
//		do{
//			word=array.get(index);
//			value=((word>>>cellShift)&valueMask);
//			if(value>thresh){return value;}//Too high; don't increment
//			int value2=min(value+incr, maxValue);
//			word2=(value2<<cellShift)|(word&~((valueMask)<<cellShift));
//		}while(word!=word2 && !array.compareAndSet(index, word, word2));
////		if(value==1){cellsUsedPersonal.incrementAndGet(num);}
//		return value;
//	}
	
	/** Raises one lane's packed cell to a target and returns its previous value. */
	private int incrementHashedLocalAndReturnUnincremented_toAtLeast(long key, int newMin){
		assert(newMin>0 && newMin<=maxValue) : newMin;
		final int num=(int)(key&arrayMask);
		final AtomicIntegerArray array=matrix[num];
		key=(key>>>arrayBits)%(cellMod);
//		key=(key>>>(arrayBits+1))%(cellMod);
		int index=(int)(key>>>indexShift);
		int cellShift=(int)(cellBits*key);
		int value, word, word2;
		final int value2s=newMin<<cellShift;
		do{
			word=array.get(index);
			value=((word>>>cellShift)&valueMask);
			if(value>=newMin){return value;}//Too high; don't increment
			word2=value2s|(word&~((valueMask)<<cellShift));
		}while(word!=word2 && !array.compareAndSet(index, word, word2));
//		if(value==1){cellsUsedPersonal.incrementAndGet(num);}
		return value;
	}
	
	/** Decrements one lane's packed cell with CAS and clamps it at zero. */
	private int decrementHashedLocal(long key, int amt){
		final int num=(int)(key&arrayMask);
		final AtomicIntegerArray array=matrix[num];
		key=(key>>>arrayBits)%(cellMod);
//		key=(key>>>(arrayBits+1))%(cellMod);
		int index=(int)(key>>>indexShift);
		int cellShift=(int)(cellBits*key);
		int value, word, word2;
		do{
			word=array.get(index);
			value=((word>>>cellShift)&valueMask);
			value=max(value-amt, 0);
			word2=(value<<cellShift)|(word&~((valueMask)<<cellShift));
		}while(word!=word2 && !array.compareAndSet(index, word, word2));
//		if(value==1){cellsUsedPersonal.incrementAndGet(num);}
		return value;
	}
	
	/**
	 * Returns the cached number of positive main-table cells, scanning lazily if unset.
	 * An early call can cache a partial construction state. Updates never invalidate
	 * this cache; the first shutdown explicitly refreshes it, but later shutdown calls
	 * do not. usedFraction() and inherited no-depth summaries/estimates share it.
	 * For a fresh scan use {@link #cellsUsed(int)} with depth1; neither scan provides
	 * a coherent snapshot while writers mutate the table. Prefilter cells are excluded.
	 */
	public long cellsUsed(){
		if(cellsUsed<0){
			synchronized(this){
				if(cellsUsed<0){
					cellsUsed=cellsUsed(1);
				}
			}
		}
		return cellsUsed;
	}
	
	@Override
	public KCountArray prefilter(){
		return prefilter;
	}
	
	@Override
	public void purgeFilter(){
		prefilter=null;
	}
	
	/** Sets the hash-mask seed, or chooses a random seed when {@code seed} is negative. */
	public static synchronized void setSeed(long seed){
		if(seed>=0){SEEDMASK=seed;}
		else{
			Random randy=shared.Shared.random();
			SEEDMASK=randy.nextLong();
		}
	}
	
	/** Whether the first shutdown finalized cached occupancy; not a mutation guard. */
	private boolean finished=false;
	
	/** Cached positive main-cell count, or -1 before scanning; updates do not invalidate it. */
	private long cellsUsed=-1;
	/** Atomic packed-count storage rows. */
	private final AtomicIntegerArray[] matrix;
	/** Number of Bloom lanes. */
	private final int hashes;
	/** Number of 32-bit words in each row. */
	private final int wordsPerArray;
	/** Number of packed cells addressed by each row. */
	private final long cellsPerArray;
	/** Row-local cell-address modulus. */
	private final long cellMod;
	/** Precomputed lane-mixing masks. */
	private final long[][] hashMasks=makeMasks(8, hashArrayLength);
	/** Prefilter result below which the main table is bypassed. */
	private final int prefilterLimit;
	
	/** Bits used for hash-mask slot indexing. */
	private static final int hashBits=6;
	/** Number of hash masks per row. */
	private static final int hashArrayLength=1<<hashBits;
	/** Mask for hash-mask slot selection. */
	private static final int hashCellMask=hashArrayLength-1;
	
	/** Optional first-stage count array. */
	private KCountArray prefilter;
	
	/** Whether primitive-key increments raise every lane under a striped lock. */
	private final transient boolean useLocks;
	/** Striped locks for multi-lane primitive-key updates. */
	private final transient Lock[] locks;
	/** Number of striped locks. */
	private static final transient int NUM_LOCKS=1999;

	/** Seed counter for mask generation. */
	private static long counter=0;
	/** User-configurable seed xor mask. */
	private static long SEEDMASK=0;
	
}
