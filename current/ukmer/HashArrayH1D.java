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
 * pair (xor, xor2) plus a count; the fingerprint-only fixed-capacity form
 * stores xor2 plus a count and uses the query's xor to select the probe chain.
 * Neither form can reconstruct or enumerate kmers, so it is intended for
 * sequence-driven operations such as bridging, error correction, and filtering
 * rather than contig seeding.
 *
 * Occupancy is encoded exclusively by values: zero is virgin, positive is
 * occupied, and negative is a fixed-table deletion tombstone.  Every 63-bit
 * hash value, including zero, is therefore valid.
 *
 * Identity is by hash, not literal k-mer equality: resizable tables compare
 * (xor, xor2), while fingerprint-only tables use the query's xor to choose the
 * probe chain and compare only stored xor2.  A collision in the compared hashes
 * can merge distinct literal keys.  Mutating and lookup methods do not take
 * {@link #getLock()} themselves; callers provide shard-level synchronization
 * when multiple threads can touch the same table.
 *
 * @author Noelle
 * @date August 30, 2026
 */
public final class HashArrayH1D extends AbstractKmerTableU{

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

	/** True when each cell stores only xor2 and the count word, not xor itself. */
	public boolean fingerprintOnly(){return hash1==null;}

	/** Increments the single count for {@code kmer}, saturating at {@link Integer#MAX_VALUE}. */
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

	/** Increments {@code kmer} and reports whether this call created a new entry. */
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

	/**
	 * Stores exactly one positive value for {@code kmer}.
	 * @return 1 when a new hash pair was created, 0 when an existing pair was overwritten
	 */
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

	/**
	 * Stores the only element of {@code vals}; hash-only tables do not support value lists.
	 * @return 1 when a new hash pair was created, 0 when an existing pair was overwritten
	 */
	@Override
	public int set(Kmer kmer, int[] vals){
		if(vals==null || vals.length!=1){
			throw new UnsupportedOperationException("HashArrayH1D stores exactly one integer value per kmer.");
		}
		return set(kmer, vals[0]);
	}

	/**
	 * Stores {@code value} only when {@code kmer} is not already present.
	 * @return 1 when a new hash pair was created, 0 when the pair was already present
	 */
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

	/** Returns the count for a full {@link Kmer}, or {@link #NOT_PRESENT}. */
	@Override
	public int getValue(Kmer kmer){return getValue(kmer.xor(), kmer.xor2());}

	/** Returns the count for an explicit key and nonnegative placement hash, or {@link #NOT_PRESENT}. */
	@Override
	public int getValue(long[] key, long xor){return getValue(xor, Kmer.xor2(key));}

	/** Returns the count for a pair of nonnegative retained hashes, or {@link #NOT_PRESENT}. */
	public int getValue(long xor, long xor2){
		final int cell=find(xor, xor2);
		return cell<0 ? NOT_PRESENT : values[cell];
	}

	/** Writes the single count for {@code kmer} into caller-owned {@code singleton}. */
	@Override
	public int[] getValues(Kmer kmer, int[] singleton){
		singleton[0]=getValue(kmer);
		return singleton;
	}

	/** Returns true when the table contains {@code kmer}. */
	@Override
	public boolean contains(Kmer kmer){return find(kmer.xor(), kmer.xor2())>=0;}

	/** Raises a present k-mer's initialized owner to {@code newOwner} atomically. */
	@Override
	public int setOwner(Kmer kmer, int newOwner){
		final int cell=find(kmer.xor(), kmer.xor2());
		if(cell<0){throw new IllegalArgumentException("Cannot own an absent kmer.");}
		if(owners==null){throw new IllegalStateException("Ownership was not initialized.");}
		int current=owners.get(cell);
		while(current<newOwner && !owners.compareAndSet(cell, current, newOwner)){current=owners.get(cell);}
		return current<newOwner ? newOwner : current;
	}

	/** Clears {@code owner} from a present k-mer with initialized ownership when it is still current. */
	@Override
	public boolean clearOwner(Kmer kmer, int owner){
		final int cell=find(kmer.xor(), kmer.xor2());
		if(cell<0){throw new IllegalArgumentException("Cannot clear ownership of an absent kmer.");}
		return owners.compareAndSet(cell, owner, NO_OWNER);
	}

	/** Returns the current initialized owner for a present k-mer. */
	@Override
	public int getOwner(Kmer kmer){
		final int cell=find(kmer.xor(), kmer.xor2());
		if(cell<0){throw new IllegalArgumentException("Absent kmer has no owner.");}
		return owners.get(cell);
	}

	/** Allocates ownership state for every primary cell. */
	@Override
	public void initializeOwnership(){
		if(owners!=null){throw new IllegalStateException("Ownership was already initialized.");}
		owners=allocAtomicInt(prime);
		for(int i=0; i<prime; i++){owners.set(i, NO_OWNER);}
	}

	/** Releases ownership state after caller-side synchronization is complete. */
	@Override
	public void clearOwnership(){owners=null;}

	/** Replaces entries with the requested final ownership state by tombstones. */
	public long removeByOwner(final int owner){
		if(owners==null){throw new IllegalStateException("Ownership was not initialized.");}
		long removed=0;
		for(int i=0; i<prime; i++){
			if(values[i]>0 && owners.get(i)==owner){
				values[i]=TOMBSTONE;
				size--;
				removed++;
			}
		}
		return removed;
	}

	/** Adds occupied counts, clamping each value to the caller-supplied {@code max} bin. */
	@Override
	public void fillHistogram(long[] counts, int max){
		for(int value : values){if(value>0){counts[Math.min(value, max)]++;}}
	}

	/** Appends all positive counts to {@code counts}. */
	@Override
	public void fillHistogram(SuperLongList counts){
		for(int value : values){if(value>0){counts.add(value);}}
	}

	/** Removes counts at or below {@code limit}; rehashes retained pairs when xor is available. */
	@Override
	public long regenerate(int limit){
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

	/** Advances to the next scheduled prime and reinserts all live hash pairs. */
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

	/** Returns true only when retained placement hashes make resizing possible in principle. */
	@Override
	boolean canResize(){return !fingerprintOnly();}

	/** Returns false because hash-only tables do not retain explicit keys for tree rebalancing. */
	@Override
	public boolean canRebalance(){return false;}

	/** Unsupported: open-addressed hash tables have no tree nodes to rebalance. */
	@Override
	public void rebalance(){throw new UnsupportedOperationException("HashArrayH1D does not rebalance.");}

	/** Returns the number of occupied hash slots. */
	@Override
	public long size(){return size;}

	/** Returns the primary array capacity. */
	@Override
	public int arrayLength(){return prime;}

	/** Unsupported: only the retained hash pair can identify a cell. */
	@Override
	Object get(long[] key){throw unavailable("lookup without the placement hash");}

	/** Unsupported: hash-only cells cannot reconstruct literal kmers for text dumping. */
	@Override
	public boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){throw unavailable("kmer dump");}

	/** Unsupported: hash-only cells cannot reconstruct literal kmers for byte dumping. */
	@Override
	public boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){throw unavailable("kmer dump");}

	/** Unsupported: hash-only cells cannot reconstruct literal kmers for threaded byte dumping. */
	@Override
	public boolean dumpKmersAsBytes_MT(ByteStreamWriter bsw, ByteBuilder bb, int k, int mincount, int maxcount, AtomicLong remaining){throw unavailable("kmer dump");}

	/** Unsupported: hash-only cells cannot count GC without explicit keys. */
	@Override
	public void countGC(long[] gcCounts, int max){throw unavailable("GC counting");}

	@Override
	Lock getLock(){return lock;}

	/** Returns the stored placement hash, or {@code -1} in fingerprint-only mode. */
	public long hash1At(int cell){return fingerprintOnly() ? -1 : hash1[cell];}

	/** Returns the stored secondary hash. */
	public long hash2At(int cell){return hash2[cell];}

	/** Returns the raw occupancy/count word for {@code cell}. */
	public int valueAt(int cell){return values[cell];}

	/** Forces one resize before retrying an insertion that found no empty cell. */
	private int retryIncrement(Kmer kmer){
		resizeOrFail();
		return increment(kmer);
	}

	/** Forces one resize before retrying a creation-counting insertion. */
	private int retryIncrementCreated(Kmer kmer){
		resizeOrFail();
		return incrementAndReturnNumCreated(kmer);
	}

	/** Finds an occupied cell matching the retained hash pair, stopping at a virgin cell. */
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

	/** Finds an existing cell or the first reusable tombstone or virgin cell. */
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

	/** Writes one occupied cell; callers must have proven it was not already live. */
	private void insert(int cell, long x1, long x2, int value){
		if(!fingerprintOnly()){hash1[cell]=x1;}
		hash2[cell]=x2;
		values[cell]=value;
		size++;
	}

	/** Reinserts an already-positive value during same-capacity rebuilds and resizes. */
	private void insertPair(long x1, long x2, int value){
		final int cell=findOrEmpty(x1, x2);
		if(cell<0){throw capacityFailure();}
		if(values[cell]>0){throw new IllegalStateException("Duplicate hash pair during rehash.");}
		insert(cell, x1, x2, value);
	}

	/** Expands only the resizable two-hash form after passing its load threshold. */
	private void growIfNeeded(){
		if(!fingerprintOnly() && size>resizeLimit){resize();}
	}

	/** Resizes when possible, otherwise throws the shared capacity failure. */
	private void resizeOrFail(){
		if(fingerprintOnly() || schedulePos>=schedule.length-1){throw capacityFailure();}
		resize();
	}

	/** Creates the capacity error used when no more scheduled cells exist. */
	private IllegalStateException capacityFailure(){
		return new IllegalStateException("Hash-only kmer table capacity exhausted at "+size+" entries in "+prime+" cells"
			+(fingerprintOnly() ? "; fingerprint-only tables cannot resize, so increase preallocation" : "; no larger scheduled table is available"));
	}

	/** Requires positive counts because zero and negative cells encode table state. */
	private static void checkValue(int value){
		if(value<=0){throw new IllegalArgumentException("HashArrayH1D reserves non-positive values for occupancy state: "+value);}
	}

	/** Creates an exception for APIs that need literal k-mer recovery. */
	private static UnsupportedOperationException unavailable(String operation){
		return new UnsupportedOperationException(operation+" requires explicit kmer storage and is unavailable in a hash-only table.");
	}

	/** Scheduled prime capacities; the fingerprint-only form must have one entry. */
	private final int[] schedule;
	/** Current index into {@link #schedule}. */
	private int schedulePos=0;
	/** Current primary capacity. */
	private int prime;
	/** Optional placement hashes; null in fingerprint-only mode. */
	private long[] hash1;
	/** Secondary hashes retained in both modes. */
	private long[] hash2;
	/** Occupancy and positive count words. */
	private int[] values;
	/** Optional per-cell owner words for ownership passes. */
	private AtomicIntegerArray owners;
	/** Occupied cell count and growth threshold. */
	private long size=0, resizeLimit;
	/** Shard lock returned through the table API. */
	private final Lock lock=new ReentrantLock();

	/** Virgin cell and deleted-cell sentinels in {@link #values}. */
	private static final int EMPTY=0, TOMBSTONE=-2;
	/** Normal and final-table load factors for resizable tables. */
	private static final float maxLoadFactor=0.88f, maxLoadFactorFinal=0.95f;
}
