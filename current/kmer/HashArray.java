package kmer;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import fileIO.ByteStreamWriter;
import fileIO.TextStreamWriter;
import shared.Primes;
import shared.Tools;
import structures.ByteBuilder;
import structures.SuperLongList;

/**
 * Shared key storage for array-backed short-key table variants.
 * A masked key selects a primary cell; up to 60 consecutive cells are probed
 * before delegating to a HashForest overflow store. Subclasses supply scalar,
 * multi-value or hybrid value storage. size() counts primary keys only;
 * victims().size() supplies the remaining keys.
 * Callers coordinate shared mutations through the table lock or exclusive access.
 * Array and victim accessors expose live storage, which resizing may replace.
 * Ownership atomics protect owner values only; initialize them after loading and
 * exclude structural mutation while claiming owners or traversing cells.
 * @author Brian Bushnell
 * @date Nov 7, 2014
 *
 */
public abstract class HashArray extends AbstractKmerTable{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructs a HashArray with a resizing schedule.
	 * Retains the caller's nonempty schedule without copying; callers must not mutate it.
	 * Uses its first size as the modulus and allocates an additional probe tail.
	 *
	 * @param schedule_ Array of sizes for automatic resizing
	 * @param coreMask_ Mask for computing hash cell positions
	 * @param twod_ Whether this table stores 2D values (arrays)
	 */
	HashArray(int[] schedule_, long coreMask_, boolean twod_){
		schedule=schedule_;
		autoResize=schedule.length>1;
		prime=schedule[0];

		sizeLimit=(long)((schedule.length==1 ? maxLoadFactorFinal : maxLoadFactor)*prime);
		array=allocLong1D(prime+extra);
		victims=new HashForest(Tools.max(10, prime/victimRatio), autoResize, twod_);
		Arrays.fill(array, NOT_PRESENT);
		twoD=twod_;
		coreMask=coreMask_;
//		coreMask2=coreMask_|3;
		coreMask2=coreMask_; //Simplifies fast fill
	}

	/**
	 * Constructs a HashArray without a schedule. The positive initial size is used
	 * directly, without prime rounding; subclass resize logic chooses later sizes.
	 *
	 * @param initialSize Initial hash table capacity
	 * @param coreMask_ Mask for computing hash cell positions
	 * @param autoResize_ Whether to automatically resize when load factor exceeded
	 * @param twod_ Whether this table stores 2D values (arrays)
	 */
	HashArray(int initialSize, long coreMask_, boolean autoResize_, boolean twod_){
		schedule=null;
		prime=initialSize;
		sizeLimit=(long)(maxLoadFactor*prime);
		array=allocLong1D(prime+extra);
		victims=new HashForest(Tools.max(10, initialSize/victimRatio), autoResize_, twod_);
		Arrays.fill(array, NOT_PRESENT);
		autoResize=autoResize_;
		twoD=twod_;
		coreMask=coreMask_;
//		coreMask2=coreMask_|3;
		coreMask2=coreMask_; //Simplifies fast fill
	}

	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/*--------------------------------------------------------------*/
	/*----------------         Test Methods         ----------------*/
	/*--------------------------------------------------------------*/

//	public final int set_Test(final long kmer, final int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
//			int[] old=getValues(kmer, new int[1]);
//			assert(old==null || contains(kmer, old));
//			if(verbose){System.err.println("Fetched "+Arrays.toString(old));}
//			x=set0(kmer, v);
//			assert(old==null || contains(kmer, old)) : "old="+Arrays.toString(old)+", v="+v+", kmer="+kmer+
//				", get(kmer)="+(Arrays.toString(getValues(kmer, new int[1])));
//			assert(contains(kmer, v));
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=set0(kmer, v);
//			assert(contains(kmer, v)) : "old="+old+", v="+v+", kmer="+kmer+", get(kmer)="+getValue(kmer);
//			assert(v==old || !contains(kmer, old));
//		}
//		return x;
//	}
//
//	public final int set_Test(final long kmer, final int v[]){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
//			final int[] singleton=new int[1];
//			int[] old=getValues(kmer, singleton);
//			assert(old==null || contains(kmer, old));
//			if(verbose){System.err.println("Before: old="+Arrays.toString(old)+", v="+Arrays.toString(v));}
//			x=set0(kmer, v);
//			if(verbose){System.err.println("After:  old="+Arrays.toString(old)+", v="+Arrays.toString(v)+", get()="+Arrays.toString(getValues(kmer, singleton)));}
//			assert(old==null || contains(kmer, old)) : "old="+Arrays.toString(old)+", v="+Arrays.toString(v)+", kmer="+kmer+
//				", get(kmer)="+(Arrays.toString(getValues(kmer, new int[1])));
//			assert(contains(kmer, v)) : "old="+Arrays.toString(old)+", v="+Arrays.toString(v)+", kmer="+kmer+
//				", get(kmer)="+(Arrays.toString(getValues(kmer, new int[1])));
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=set0(kmer, v);
//			assert(contains(kmer, v)) : "old="+old+", v="+v+", kmer="+kmer+", get(kmer)="+getValue(kmer);
//			assert(v[0]==old || !contains(kmer, old));
//		}
//		return x;
//	}
//
//	public final int setIfNotPresent_Test(long kmer, int v){
//		assert(TESTMODE);
//		final int x;
//		if(TWOD){
////			int[] vals=getValues(kmer, null);
////			assert(vals==null || contains(kmer, vals));
////			x=setIfNotPresent(kmer, v);
////			assert(contains(kmer, vals));
////			assert(contains(kmer, v));
//			x=0;
//			assert(false);
//		}else{
//			int old=getValue(kmer);
//			assert(old==0 || old==-1 || contains(kmer, old));
//			x=setIfNotPresent0(kmer, v);
//			assert((old<1 && contains(kmer, v)) || (old>0 && contains(kmer, old))) : kmer+", "+old+", "+v;
//		}
//		return x;
//	}

	/*--------------------------------------------------------------*/

	/**
	 * Maps a k-mer to its primary hash table cell position.
	 * @param kmer The k-mer to hash
	 * @return Cell index in the hash table array
	 */
	public final int kmerToCell(long kmer){
		final int cell=(int)((kmer&coreMask)%prime);
		return cell;
	}

	/** Inserts or updates values through the subclass contract, spilling a full probe
	 * window to the victim forest. Returns one for a new key and zero for an existing key.
	 * Copying and valid-prefix handling depend on the chosen value representation. */
	@Override
	public final int set(final long kmer, final int[] v, final int vlen){
		int cell=kmerToCell(kmer);

		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				if(verbose){System.err.println("A2: Adding "+kmer+", "+Arrays.toString(v)+", "+cell);}
				insertValue(kmer, v, cell, vlen);
				if(verbose){System.err.println("A2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
				return 0;
			}else if(n==NOT_PRESENT){
				if(verbose){System.err.println("B2: Adding "+kmer+", "+Arrays.toString(v)+", "+cell);}
				array[cell]=kmer;
				insertValue(kmer, v, cell, vlen);
				if(verbose){System.err.println("B2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
				size++;
				if(autoResize && size+victims.size>sizeLimit){resize();}
				return 1;
			}
		}
		if(verbose){System.err.println("C2: Adding "+kmer+", "+v+", "+cell);}
		final int x=victims.set(kmer, v, vlen);
		if(autoResize && size+victims.size>sizeLimit){resize();}
		if(verbose){System.err.println("C2: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
		return x;
	}

	/** Inserts or updates one value through the subclass contract; returns the new-key count.
	 * Scalar tables replace values, while multi-value tables can merge them. */
	@Override
	public final int set(final long kmer, final int v){
		int cell=kmerToCell(kmer);

//		assert(TESTMODE);
//		ll.add(kmer);
//		il.add(v);

		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				if(verbose){System.err.println("A1: Adding "+kmer+", "+v+", "+cell);}
				insertValue(kmer, v, cell);
				if(verbose){System.err.println("A1: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
				return 0;
			}else if(n==NOT_PRESENT){
				if(verbose){System.err.println("B1: Adding "+kmer+", "+v+", "+cell);}
				array[cell]=kmer;
				insertValue(kmer, v, cell);
				if(verbose){System.err.println("B1: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1])));}
				size++;
				if(autoResize && size+victims.size>sizeLimit){resize();}
				return 1;
			}
		}
		if(verbose){System.err.println("C1: Adding "+kmer+", "+v+", "+cell+
				"; victims.get(kmer)="+Arrays.toString(victims.getValues(kmer, new int[1])));}
		final int x=victims.set(kmer, v);
		if(autoResize && size+victims.size>sizeLimit){resize();}
		if(verbose){System.err.println("C1: getValues("+kmer+") = "+Arrays.toString(getValues(kmer, new int[1]))+
				"; victims.get(kmer)="+Arrays.toString(victims.getValues(kmer, new int[1])));}
		return x;
	}

	/** Inserts only an absent key, preserving all values of an existing key.
	 * @return One for a new key, zero otherwise */
	@Override
	public final int setIfNotPresent(long kmer, int value){
		int cell=kmerToCell(kmer);
		for(final int max=cell+extra; cell<max; cell++){
			long n=array[cell];
			if(n==kmer){
				return 0;
			}else if(n==NOT_PRESENT){
				array[cell]=kmer;
				insertValue(kmer, value, cell);
				size++;
				if(autoResize && size+victims.size>sizeLimit){resize();}
				return 1;
			}
		}
//		System.err.println("size="+size+", prime="+prime+", limit="+sizeLimit);
		int x=victims.setIfNotPresent(kmer, value);
		if(autoResize && size+victims.size>sizeLimit){resize();}
		return x;
	}

	/** Returns the representation's scalar view, or NOT_PRESENT for an absent key. */
	@Override
	public final int getValue(long kmer){
		int cell=findKmer(kmer);
		if(cell==NOT_PRESENT){return NOT_PRESENT;}
		if(cell==HASH_COLLISION){return victims.getValue(kmer);}
		return readCellValue(cell);
	}

	/**
	 * Retrieves the single value associated with a k-mer, starting search from specified cell.
	 * @param kmer The k-mer to look up
	 * @param startCell Primary cell computed for this key in the current table; not a prior probe hit
	 * @return The value associated with the k-mer, or NOT_PRESENT if not found
	 */
	public final int getValue(long kmer, int startCell){
		int cell=findKmer(kmer, startCell);
		if(cell==NOT_PRESENT){return NOT_PRESENT;}
		if(cell==HASH_COLLISION){return victims.getValue(kmer);}
		return readCellValue(cell);
	}

	/** Returns value storage or fills the caller's singleton (length at least one).
	 * Absence sets singleton[0] to NOT_PRESENT. Returned multi-value arrays may be
	 * live table storage; callers must not modify them. */
	@Override
	public final int[] getValues(long kmer, int[] singleton){
		int cell=findKmer(kmer);
		if(cell==NOT_PRESENT){
			singleton[0]=NOT_PRESENT;
			return singleton;
		}
		if(cell==HASH_COLLISION){return victims.getValues(kmer, singleton);}
		return readCellValues(cell, singleton);
	}

	/** Tests key presence in the primary array and, when needed, the victim forest. */
	@Override
	public final boolean contains(long kmer){
		int cell=findKmer(kmer);
		if(cell==NOT_PRESENT){return false;}
		if(cell==HASH_COLLISION){return victims.contains(kmer);}
		return true;
	}

	/**
	 * Retrieves the k-mer stored at a specific cell position.
	 * @param cell The cell index
	 * @return The k-mer at that position
	 */
	public final long getKmer(int cell){
		return array[cell];
	}

	/*--------------------------------------------------------------*/
	/*----------------          Ownership           ----------------*/
	/*--------------------------------------------------------------*/

	/** Allocates/reset owner markers to NO_OWNER after loading, including victims.
	 * Requires exclusive access and no existing primary owner array. */
	@Override
	public final void initializeOwnership(){
		assert(owners==null);
		owners=allocAtomicInt(array.length);
		for(int i=0; i<array.length; i++){
			owners.set(i, NO_OWNER);
		}
		victims.initializeOwnership();
	}

	/** Releases primary owner storage and resets victim ownership under exclusive access. */
	@Override
	public final void clearOwnership(){
		owners=null;
		victims.clearOwnership();
	}

	/** Raises an existing key's owner; requires initialized owners and stable structure. */
	@Override
	public final int setOwner(final long kmer, final int newOwner){
		final int cell=findKmer(kmer);
		assert(cell!=NOT_PRESENT);
		if(cell==HASH_COLLISION){return victims.setOwner(kmer, newOwner);}
		return setOwner(kmer, newOwner, cell);
	}

	/**
	 * Sets the thread owner for a k-mer at a specific cell.
	 * Raises the owner with compare-and-set; concurrent lower claims do not reduce it.
	 * Requires initialized ownership storage and a stable key-to-cell mapping.
	 * This does not lock the key/value arrays or exclude resizing.
	 *
	 * @param kmer The k-mer to claim ownership of
	 * @param newOwner Thread ID of the new owner
	 * @param cell Current primary-array cell containing the k-mer
	 * @return An owner observed or installed by this call (may be higher than newOwner)
	 */
	public final int setOwner(final long kmer, final int newOwner, final int cell){
		assert(array[cell]==kmer);
		final int original=owners.get(cell);
		int current=original;
		while(current<newOwner){
			boolean success=owners.compareAndSet(cell, current, newOwner);
			if(!success){current=owners.get(cell);}
			else{current=newOwner;}
		}
		assert(current>=original) : "original="+original+", current="+current+", newOwner="+newOwner+", re-read="+owners.get(cell);
		return current;
	}

	/** Atomically clears a matching owner of an existing key, with stable structure. */
	@Override
	public final boolean clearOwner(final long kmer, final int owner){
		final int cell=findKmer(kmer);
		assert(cell!=NOT_PRESENT);
		if(cell==HASH_COLLISION){return victims.clearOwner(kmer, owner);}
		return clearOwner(kmer, owner, cell);
	}

	/**
	 * Clears ownership of a k-mer at a specific cell if owned by the specified thread.
	 *
	 * @param kmer The k-mer to release
	 * @param owner Thread ID that should currently own the k-mer
	 * @param cell Current primary-array cell containing the k-mer
	 * @return true if ownership was successfully cleared, false otherwise
	 */
	public final boolean clearOwner(final long kmer, final int owner, final int cell){
		assert(array[cell]==kmer);
		boolean success=owners.compareAndSet(cell, owner, NO_OWNER);
		return success;
	}

	/** Reads an existing key's owner; requires initialized owners and stable structure. */
	@Override
	public final int getOwner(final long kmer){
		final int cell=findKmer(kmer);
		assert(cell!=NOT_PRESENT);
		if(cell==HASH_COLLISION){return victims.getOwner(kmer);}
		return getCellOwner(cell);
	}

	/**
	 * Gets the thread owner of a specific cell.
	 * @param cell The cell index to check
	 * @return Thread ID of the current owner
	 */
	public final int getCellOwner(final int cell){
		return owners.get(cell);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Nonpublic Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Updates a hybrid table's scalar count in its two-dimensional victim forest.
	 * Existing nodes must contain exactly one nonnegative count. Replace that value
	 * in place, preserving its logical length; ordinary 2D node setters merge values.
	 * Absent keys use the forest's existing insertion, accounting and resize path.
	 * Requires the caller's mutation exclusion, just like primary-array counting.
	 * @param kmer Key whose primary probe window is full
	 * @param incr Positive amount to add, saturating at Integer.MAX_VALUE
	 * @param returnNumCreated Whether to return a new-key indicator instead of the count
	 * @return Updated count, or one for a new key and zero for an existing key
	 */
	protected final int incrementScalarVictim(final long kmer, final int incr, final boolean returnNumCreated){
		assert(twoD && incr>0) : "Hybrid scalar victim counting requires 2D storage and a positive increment; twoD="+twoD+", incr="+incr;
		final KmerNode node=victims.get(kmer);
		if(node==null){
			return returnNumCreated ? victims.incrementAndReturnNumCreated(kmer, incr) : victims.increment(kmer, incr);
		}
		assert(node.numValues()==1) : "Hybrid counting cannot interpret association lists as counts; kmer="+kmer+", numValues="+node.numValues();
		final int[] vals=node.values(null);
		assert(vals!=null && vals.length>0 && vals[0]>=0) : "Hybrid scalar victims must retain one nonnegative count for lookup and resize; kmer="+kmer;
		//KMR021: Node2D.set would append the sum and leave the old count in slot zero.
		final int sum=vals[0]+incr;
		vals[0]=(sum<0 ? Integer.MAX_VALUE : sum);
		return returnNumCreated ? 0 : vals[0];
	}

	/**
	 * Inserts a single value for a k-mer at a specific cell.
	 * Implementation varies by subclass depending on storage strategy.
	 *
	 * @param kmer The k-mer being stored
	 * @param v The value to insert
	 * @param cell The cell position where the k-mer is stored
	 */
	protected abstract void insertValue(final long kmer, final int v, final int cell);

//	protected abstract void insertValue(final long kmer, final int[] vals, final int cell);

	/** Stores values in an already assigned primary cell. Subclasses define copy/merge
	 * behavior and use of the valid prefix; vlen supports HashArrayHybridFast's IntList3.
	 * Overflow values instead go through HashForest.set. */
	protected abstract void insertValue(final long kmer, final int[] vals, final int cell, final int vlen);

	/**
	 * Reads the single value stored at a specific cell.
	 * @param cell The cell position to read from
	 * @return The value stored at that position
	 */
	protected abstract int readCellValue(int cell);
	/**
	 * Reads all values stored at a specific cell.
	 * @param cell The cell position to read from
	 * @param singleton Array to reuse for single values
	 * @return Subclass-owned value storage or the supplied singleton; callers must not mutate retained arrays
	 */
	protected abstract int[] readCellValues(int cell, int[] singleton);

	/** Unsupported node access: primary entries do not have node objects. */
	@Override
	final Object get(long kmer){
		throw new RuntimeException("Unimplemented.");
	}

	/**
	 * Finds the cell position of a k-mer in the hash table.
	 * The returned cell is borrowed and remains valid only while the table structure is stable.
	 * @param kmer The k-mer to find
	 * @return Matching cell, NOT_PRESENT at the first empty cell, or HASH_COLLISION
	 * when the full probe window is occupied; the latter requires a separate victim lookup
	 */
	public final int findKmer(long kmer){
		return findKmer(kmer, kmerToCell(kmer));
	}

	/**
	 * Finds the cell position of a k-mer starting from a specific cell.
	 * Uses linear probing within the extra cell window.
	 *
	 * @param kmer The k-mer to find
	 * @param startCell Current primary hash cell; the full probe window must fit the array
	 * @return Cell position, NOT_PRESENT if not found, or HASH_COLLISION if probing limit reached
	 */
	final int findKmer(final long kmer, final int startCell){
		//Comprehension (clever, verified): three-way return is the crux of the whole design. A hit returns the cell; a
		//NOT_PRESENT sentinel returns immediately (linear insertion guarantees the kmer cannot be past the first gap);
		//exhausting the `extra`-wide probe window returns HASH_COLLISION — meaning "if present at all, it overflowed to
		//the victim forest." set() spills to victims under exactly the same full-window condition, so the two stay in
		//lockstep: a kmer is in the array iff findable within the window, else in victims, never lost between them.
		int cell=startCell;
		for(final int max=cell+extra; cell<max; cell++){
			final long n=array[cell];
			if(n==kmer){return cell;}
			else if(n==NOT_PRESENT){return NOT_PRESENT;}
		}
		return HASH_COLLISION;
	}

	/**
	 * Finds a k-mer or the first empty cell suitable for storing it.
	 * @param kmer The k-mer to find or place
	 * @return Cell position of k-mer or empty cell, or HASH_COLLISION if no space available
	 */
	final int findKmerOrEmpty(long kmer){
		int cell=kmerToCell(kmer);
		for(final int max=cell+extra; cell<max; cell++){
			final long n=array[cell];
			if(n==kmer || n==NOT_PRESENT){return cell;}
		}
		return HASH_COLLISION;
	}

	/*--------------------------------------------------------------*/
	/*----------------   Resizing and Rebalancing   ----------------*/
	/*--------------------------------------------------------------*/

	/** Reports resize capability, independently of this instance's automatic-growth setting. */
	@Override
	final boolean canResize(){return true;}

	/** Returns the primary-array key count, excluding the victim forest. */
	@Override
	final public long size(){return size;}

	/** Returns allocated key slots including the extra probe tail. */
	@Override
	final public int arrayLength(){return array.length;}

	/** Rebuilds key/value storage according to the subclass policy; callers exclude access. */
	@Override
	protected abstract void resize();

	/*--------------------------------------------------------------*/
	/*----------------         Info Dumping         ----------------*/
	/*--------------------------------------------------------------*/

	/** Emits primary entries then victims without closing or checking the writer.
	 * Scalar primary values use count bounds; two-dimensional primary values do not.
	 * @return Always true; completion and error reporting remain the caller's responsibility */
	@Override
	public final boolean dumpKmersAsText(TextStreamWriter tsw, int k, int mincount, int maxcount){
		if(twoD){
			//Comprehension [kmer/HashArray 2D-dump, NOT a bug]: the 2D branch intentionally applies NEITHER mincount NOR
			//maxcount (all three dump variants do the same) — a 2D entry holds an int[] of values (e.g. multiple scaffold
			//IDs, Seal/BBDuk ref tables), so a scalar count threshold has no well-defined meaning here. Masked in any case:
			//every 2D-table dump caller (Seal:766/1543, BBDuk ref) passes mincount<=1 and maxcount=Integer.MAX_VALUE. The
			//1D branch below is the count-filtered path; keep 2D unfiltered unless a per-value count semantic is defined.
			final int[] singleton=new int[1];
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=NOT_PRESENT){
					tsw.print(toText(kmer, readCellValues(i, singleton), k).append('\n'));
				}
			}
		}else{
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=NOT_PRESENT && (mincount<2 || readCellValue(i)>=mincount) && readCellValue(i)<=maxcount){
					tsw.print(toText(kmer, readCellValue(i), k).append('\n'));
				}
			}
		}
		if(victims!=null){
			victims.dumpKmersAsText(tsw, k, mincount, maxcount);
		}
		return true;
	}

	/** Emits primary entries then victims, reserving each record from remaining when supplied.
	 * A negative reservation stops output. Scalar primary values use count bounds;
	 * two-dimensional primary values do not. Does not finish/check the writer.
	 * @return Always true */
	@Override
	public final boolean dumpKmersAsBytes(ByteStreamWriter bsw, int k, int mincount, int maxcount, AtomicLong remaining){
		if(twoD){
			final int[] singleton=new int[1];
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=NOT_PRESENT){
					if(remaining!=null && remaining.decrementAndGet()<0){return true;}
					bsw.printlnKmer(kmer, readCellValues(i, singleton), k);
				}
			}
		}else{
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=NOT_PRESENT && (mincount<2 || readCellValue(i)>=mincount) && readCellValue(i)<=maxcount){
					if(remaining!=null && remaining.decrementAndGet()<0){return true;}
					bsw.printlnKmer(kmer, readCellValue(i), k);
				}
			}
		}
		if(victims!=null){
			victims.dumpKmersAsBytes(bsw, k, mincount, maxcount, remaining);
		}
		return true;
	}

	/** Appends records to an exclusive caller buffer, submitting copies at 16000 bytes.
	 * The caller must submit the final partial buffer and finish/check the writer.
	 * Uses the same primary filtering/reservation rules as the direct byte dump;
	 * creates no worker threads and delegates overflow traversal to the victim forest.
	 * @return Always true */
	@Override
	public final boolean dumpKmersAsBytes_MT(final ByteStreamWriter bsw, final ByteBuilder bb, final int k, final int mincount, int maxcount, AtomicLong remaining){
		if(twoD){
			final int[] singleton=new int[1];
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=NOT_PRESENT){
					if(remaining!=null && remaining.decrementAndGet()<0){return true;}
					toBytes(kmer, readCellValues(i, singleton), k, bb);
					bb.nl();
					if(bb.length()>=16000){
						ByteBuilder bb2=new ByteBuilder(bb);
						synchronized(bsw){bsw.addJob(bb2);}
						bb.clear();
					}
				}
			}
		}else{
			for(int i=0; i<array.length; i++){
				long kmer=array[i];
				if(kmer!=NOT_PRESENT && (mincount<2 || readCellValue(i)>=mincount) && readCellValue(i)<=maxcount){
					if(remaining!=null && remaining.decrementAndGet()<0){return true;}
					toBytes(kmer, readCellValue(i), k, bb);
					bb.nl();
					if(bb.length()>=16000){
						ByteBuilder bb2=new ByteBuilder(bb);
						synchronized(bsw){bsw.addJob(bb2);}
						bb.clear();
					}
				}
			}
		}
		if(victims!=null){
			victims.dumpKmersAsBytes_MT(bsw, bb, k, mincount, maxcount, remaining);
		}
		return true;
	}

	/** Adds one observation per key to min(scalarValue, max), then visits victims.
	 * Caller supplies nonnegative scalar values and an array covering index max. */
	@Override
	public final void fillHistogram(long[] ca, int max){
		for(int i=0; i<array.length; i++){
			long kmer=array[i];
			if(kmer!=NOT_PRESENT){
				int count=Tools.min(readCellValue(i), max);
				ca[count]++;
			}
		}
		if(victims!=null){
			victims.fillHistogram(ca, max);
		}
	}

	/** Appends each primary scalar value, then delegates victim values to the list. */
	@Override
	public void fillHistogram(SuperLongList sll){
		for(int i=0; i<array.length; i++){
			long kmer=array[i];
			if(kmer!=NOT_PRESENT){
				int count=readCellValue(i);
				sll.add(count);
			}
		}
		if(victims!=null){
			victims.fillHistogram(sll);
		}
	}

	/** Adds each key's GC bases to its capped scalar-value bin, including victims.
	 * Requires nonnegative values and a gcCounts array covering index max. */
	@Override
	public final void countGC(long[] gcCounts, int max){
		for(int i=0; i<array.length; i++){
			long kmer=array[i];
			if(kmer!=NOT_PRESENT){
				int count=readCellValue(i);
				int index=Tools.min(count, max);
				gcCounts[index]+=gc(kmer);
			}
		}
		if(victims!=null){
			victims.countGC(gcCounts, max);
		}
	}

	/** Returns the victim cache used for collision handling.
	 * @return The HashForest victim cache */
	public HashForest victims(){
		return victims;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Per-cell owner IDs, allocated by initializeOwnership and released by clearOwnership. */
	AtomicIntegerArray owners;
	/** Main hash table array storing k-mers */
	long[] array;
	/** Current primary modulus, excluding the allocated probe tail. */
	int prime;
	/** Number of k-mers currently stored in main table */
	long size=0;
	/** Combined primary-plus-victim key threshold for an automatic resize attempt. */
	long sizeLimit;
	/** Victim cache for handling hash collisions */
	final HashForest victims;
	/** Whether automatic resizing is enabled */
	final boolean autoResize;
	/** Whether this table stores 2D values (arrays) */
	public final boolean twoD;
	/** Caller-managed coordination lock; ordinary table methods do not acquire it. */
	private final Lock lock=new ReentrantLock();
	/** Mask for extracting core bits from k-mer for hashing */
	private final long coreMask;//for ways
	/** Retained secondary-mask state; currently equal to coreMask and unused by this class. */
	private final long coreMask2;//for cells

	/** Advances within a nonnull schedule, or repeats its last entry at exhaustion.
	 * @return The next scheduled size for resizing */
	protected int nextScheduleSize(){
		if(schedulePos<schedule.length-1){schedulePos++;}
		return schedule[schedulePos];
	}

	/** Checks whether the current position is the last entry of a nonnull schedule.
	 * @return true if at maximum size, false if further resizing is possible */
	protected boolean atMaxSize(){
		return schedulePos>=schedule.length-1;
	}

	/** Borrowed nonempty resize schedule, or null for the unscheduled constructor. */
	protected final int[] schedule;
	/** Current position in the resize schedule */
	private int schedulePos=0;

	/** Returns live key storage; the reference becomes stale after resizing.
	 * @return The long array storing k-mers */
	public long[] array(){return array;}

	/** Returns live ownership storage, or null when ownership is not initialized.
	 * @return The atomic integer array tracking thread ownership */
	public AtomicIntegerArray owners(){return owners;}
	/** Returns the caller-managed coordination lock. */
	@Override
	final Lock getLock(){return lock;}

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Initial divisor for victim cache size relative to main table */
	final static int victimRatio=16; //Initial divisor for victim cache size; it self-resizes.
	/** Number of extra cells allocated for linear probing */
	final static int extra=60; //Amazingly, increasing this gave increasing returns past 300.  Old default was 21.  Could allow higher maxLoadFactorFinal and smaller victim cache.
	/** Maximum prime number usable for hash table size */
	final static int maxPrime=Primes.primeAtMost(Integer.MAX_VALUE-extra-20);
	/** Minimum resize multiplier when growing the table */
	final static float resizeMult=2f; //Resize by a minimum of this much; not needed for schedule
	/** Minimum load factor after resizing */
	final static float minLoadFactor=0.58f; //Resize by enough to get the load above this factor; not needed for schedule
	/** Load factor threshold that triggers automatic resizing */
	final static float maxLoadFactor=0.88f; //Reaching this load triggers resizing
	/** Maximum load factor for final size without further resizing */
	final static float maxLoadFactorFinal=0.95f; //Reaching this load triggers killing
	/** Inverse of minimum load factor for calculations */
	final static float minLoadMult=1/minLoadFactor;
	/** Inverse of maximum load factor for calculations */
	final static float maxLoadMult=1/maxLoadFactor;

}
