package ukmer;

import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * Hash-only table for long kmers.  The resizable form stores the independent
 * pair (xor, xor2); the fixed-capacity form stores only xor2 and uses the
 * query's xor to select the probe chain.  Neither form can reconstruct or
 * enumerate kmers, so it is intended for sequence-driven operations such as
 * bridging, error correction, and filtering rather than contig seeding.
 *
 * Occupancy is encoded exclusively by values: zero is virgin, positive is
 * occupied, and negative is a fixed-table deletion tombstone.  Every 63-bit
 * hash value, including zero, is therefore valid.
 *
 * @author Noelle
 * @date August 30, 2026
 */
public final class HashArrayH1D extends AbstractKmerTableU {

	/** Makes a resizable two-hash table from a normal size schedule. */
	public HashArrayH1D(int[] schedule_){
		this(schedule_, true);
	}

	/**
	 * @param schedule_ Prime capacities; fixed tables use only schedule_[0]
	 * @param storePlacementHash Store xor as well as xor2, enabling resize
	 */
	public HashArrayH1D(int[] schedule_, boolean storePlacementHash){
		if(schedule_==null || schedule_.length<1 || schedule_[0]<1){
			throw new IllegalArgumentException("A positive hash-table schedule is required.");
		}
		if(!storePlacementHash && schedule_.length!=1){
			throw new IllegalArgumentException("Fingerprint-only tables cannot resize; use a one-element schedule.");
		}
		schedule=schedule_;
		prime=schedule[0];
		hash1=(storePlacementHash ? allocLong1D(prime) : null);
		hash2=allocLong1D(prime);
		values=allocInt1D(prime);
		resizeLimit=(long)((schedule.length==1 ? maxLoadFactorFinal : maxLoadFactor)*prime);
	}

	/** True for the 8-byte fingerprint-only, deliberately non-resizable form. */
	public boolean fingerprintOnly(){return hash1==null;}

	@Override
	public int increment(Kmer kmer){
		final int cell=findOrEmpty(kmer.xor(), kmer.xor2());
		if(cell<0){return retryIncrement(kmer);}
		if(values[cell]>0){
			if(values[cell]<Integer.MAX_VALUE){values[cell]++;}
			return values[cell];
		}
		insert(cell, kmer.xor(), kmer.xor2(), 1);
		growIfNeeded();
		return 1;
	}

	@Override
	public int incrementAndReturnNumCreated(Kmer kmer){
		final int cell=findOrEmpty(kmer.xor(), kmer.xor2());
		if(cell<0){return retryIncrementCreated(kmer);}
		if(values[cell]>0){
			if(values[cell]<Integer.MAX_VALUE){values[cell]++;}
			return 0;
		}
		insert(cell, kmer.xor(), kmer.xor2(), 1);
		growIfNeeded();
		return 1;
	}

	@Override
	public int set(Kmer kmer, int value){
		checkValue(value);
		final long x1=kmer.xor(), x2=kmer.xor2();
		int cell=findOrEmpty(x1, x2);
		if(cell<0){
			resizeOrFail();
			cell=findOrEmpty(x1, x2);
			if(cell<0){throw capacityFailure();}
		}
		final boolean created=values[cell]<=0;
		if(created){insert(cell, x1, x2, value);}
		else{values[cell]=value;}
		growIfNeeded();
		return created ? 1 : 0;
	}

	@Override
	public int set(Kmer kmer, int[] vals){
		if(vals==null || vals.length!=1){
			throw new UnsupportedOperationException("HashArrayH1D stores exactly one integer value per kmer.");
		}
		return set(kmer, vals[0]);
	}

	@Override
	public int setIfNotPresent(Kmer kmer, int value){
		checkValue(value);
		final long x1=kmer.xor(), x2=kmer.xor2();
		int cell=findOrEmpty(x1, x2);
		if(cell<0){
			resizeOrFail();
			cell=findOrEmpty(x1, x2);
			if(cell<0){throw capacityFailure();}
		}
		if(values[cell]>0){return 0;}
		insert(cell, x1, x2, value);
		growIfNeeded();
		return 1;
	}

	@Override
	public int getValue(Kmer kmer){return getValue(kmer.xor(), kmer.xor2());}

	@Override
	public int getValue(long[] key, long xor){return getValue(xor, Kmer.xor2(key));}

	public int getValue(long xor, long xor2){
		final int cell=find(xor, xor2);
		return cell<0 ? NOT_PRESENT : values[cell];
	}

	@Override
	public int[] getValues(Kmer kmer, int[] singleton){
		singleton[0]=getValue(kmer);
		return singleton;
	}

	@Override
	public boolean contains(Kmer kmer){return find(kmer.xor(), kmer.xor2())>=0;}

	@Override
	public int setOwner(Kmer kmer, int newOwner){
		final int cell=find(kmer.xor(), kmer.xor2());
		if(cell<0){throw new IllegalArgumentException("Cannot own an absent kmer.");}
		if(owners==null){throw new IllegalStateException("Ownership was not initialized.");}
		int current=owners.get(cell);
		while(current<newOwner && !owners.compareAndSet(cell, current, newOwner)){current=owners.get(cell);}
		return current<newOwner ? newOwner : current;
	}

	@Override
	public boolean clearOwner(Kmer kmer, int owner){
		final int cell=find(kmer.xor(), kmer.xor2());
		if(cell<0){throw new IllegalArgumentException("Cannot clear ownership of an absent kmer.");}
		return owners.compareAndSet(cell, owner, NO_OWNER);
	}

	@Override
	public int getOwner(Kmer kmer){
		final int cell=find(kmer.xor(), kmer.xor2());
		if(cell<0){throw new IllegalArgumentException("Absent kmer has no owner.");}
		return owners.get(cell);
	}

	@Override
	public void initializeOwnership(){
		if(owners!=null){throw new IllegalStateException("Ownership was already initialized.");}
		owners=allocAtomicInt(prime);
		for(int i=0; i<prime; i++){owners.set(i, NO_OWNER);}
	}

	@Override
	public void clearOwnership(){owners=null;}

	@Override
	public void fillHistogram(long[] counts, int max){
		for(int value : values){if(value>0){counts[Math.min(value, max)]++;}}
	}

	@Override
	public void fillHistogram(SuperLongList counts){
		for(int value : values){if(value>0){counts.add(value);}}
	}

	@Override
	long regenerate(int limit){
		if(owners!=null){throw new IllegalStateException("Clear ownership before regeneration.");}
		long removed=0;
		if(fingerprintOnly()){
			for(int i=0; i<prime; i++){
				if(values[i]>0 && values[i]<=limit){values[i]=TOMBSTONE; size--; removed++;}
			}
			return removed;
		}
		final long[] old1=hash1, old2=hash2;
		final int[] oldValues=values;
		hash1=allocLong1D(prime);
		hash2=allocLong1D(prime);
		values=allocInt1D(prime);
		final long oldSize=size;
		size=0;
		for(int i=0; i<oldValues.length; i++){
			if(oldValues[i]>limit){insertPair(old1[i], old2[i], oldValues[i]);}
			else if(oldValues[i]>0){removed++;}
		}
		if(oldSize!=size+removed){throw new IllegalStateException("Hash regeneration lost entries: "+oldSize+" -> "+size+" + "+removed);}
		return removed;
	}

	@Override
	synchronized void resize(){
		if(fingerprintOnly()){throw new UnsupportedOperationException("A fingerprint-only hash table cannot resize because xor was not stored.");}
		if(owners!=null){throw new IllegalStateException("Clear ownership before resizing.");}
		if(schedulePos>=schedule.length-1){throw capacityFailure();}
		final long[] old1=hash1, old2=hash2;
		final int[] oldValues=values;
		final long oldSize=size;
		prime=schedule[++schedulePos];
		hash1=allocLong1D(prime);
		hash2=allocLong1D(prime);
		values=allocInt1D(prime);
		size=0;
		resizeLimit=(long)((schedulePos==schedule.length-1 ? maxLoadFactorFinal : maxLoadFactor)*prime);
		for(int i=0; i<oldValues.length; i++){
			if(oldValues[i]>0){insertPair(old1[i], old2[i], oldValues[i]);}
		}
		if(size!=oldSize){throw new IllegalStateException("Hash resize lost entries: "+oldSize+" -> "+size);}
	}

	@Override
	boolean canResize(){return !fingerprintOnly();}

	@Override
	public boolean canRebalance(){return false;}

	@Override
	public void rebalance(){throw new UnsupportedOperationException("HashArrayH1D does not rebalance.");}

	@Override
	public long size(){return size;}

	@Override
	public int arrayLength(){return prime;}

	@Override
	Object get(long[] key){throw unavailable("lookup without the placement hash");}

	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){throw unavailable("kmer dump");}

	@Override
	public boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){throw unavailable("kmer dump");}

	@Override
	public boolean dumpKmersAsBytes_MT(ByteStreamWriter bsw, ByteBuilder bb, int k, int mincount, int maxcount, AtomicLong remaining){throw unavailable("kmer dump");}

	@Override
	public void countGC(long[] gcCounts, int max){throw unavailable("GC counting");}

	@Override
	Lock getLock(){return lock;}

	public long hash1At(int cell){return fingerprintOnly() ? -1 : hash1[cell];}
	public long hash2At(int cell){return hash2[cell];}
	public int valueAt(int cell){return values[cell];}

	private int retryIncrement(Kmer kmer){
		resizeOrFail();
		return increment(kmer);
	}

	private int retryIncrementCreated(Kmer kmer){
		resizeOrFail();
		return incrementAndReturnNumCreated(kmer);
	}

	private int find(long x1, long x2){
		int cell=(int)(x1%prime);
		for(int probes=0; probes<prime; probes++){
			final int value=values[cell];
			if(value==EMPTY){return NOT_PRESENT;}
			if(value>0 && hash2[cell]==x2 && (fingerprintOnly() || hash1[cell]==x1)){return cell;}
			if(++cell==prime){cell=0;}
		}
		return NOT_PRESENT;
	}

	private int findOrEmpty(long x1, long x2){
		int cell=(int)(x1%prime), tombstone=NOT_PRESENT;
		for(int probes=0; probes<prime; probes++){
			final int value=values[cell];
			if(value==EMPTY){return tombstone>=0 ? tombstone : cell;}
			if(value==TOMBSTONE){if(tombstone<0){tombstone=cell;}}
			else if(hash2[cell]==x2 && (fingerprintOnly() || hash1[cell]==x1)){return cell;}
			if(++cell==prime){cell=0;}
		}
		return tombstone;
	}

	private void insert(int cell, long x1, long x2, int value){
		if(!fingerprintOnly()){hash1[cell]=x1;}
		hash2[cell]=x2;
		values[cell]=value;
		size++;
	}

	private void insertPair(long x1, long x2, int value){
		final int cell=findOrEmpty(x1, x2);
		if(cell<0){throw capacityFailure();}
		if(values[cell]>0){throw new IllegalStateException("Duplicate hash pair during rehash.");}
		insert(cell, x1, x2, value);
	}

	private void growIfNeeded(){
		if(!fingerprintOnly() && size>resizeLimit){resize();}
	}

	private void resizeOrFail(){
		if(fingerprintOnly() || schedulePos>=schedule.length-1){throw capacityFailure();}
		resize();
	}

	private IllegalStateException capacityFailure(){
		return new IllegalStateException("Hash-only kmer table capacity exhausted at "+size+" entries in "+prime+" cells"
			+(fingerprintOnly() ? "; fingerprint-only tables cannot resize, so increase preallocation" : "; no larger scheduled table is available"));
	}

	private static void checkValue(int value){
		if(value<=0){throw new IllegalArgumentException("HashArrayH1D reserves non-positive values for occupancy state: "+value);}
	}

	private static UnsupportedOperationException unavailable(String operation){
		return new UnsupportedOperationException(operation+" requires explicit kmer storage and is unavailable in a hash-only table.");
	}

	private final int[] schedule;
	private int schedulePos=0;
	private int prime;
	private long[] hash1;
	private long[] hash2;
	private int[] values;
	private AtomicIntegerArray owners;
	private long size=0, resizeLimit;
	private final Lock lock=new ReentrantLock();

	private static final int EMPTY=0, TOMBSTONE=-2;
	private static final float maxLoadFactor=0.88f, maxLoadFactorFinal=0.95f;
}
