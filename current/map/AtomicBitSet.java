package map;

import java.util.concurrent.atomic.AtomicIntegerArray;

//NOTE [map/AtomicBitSet#003] (LOW/latent, NOT patched - 0 callers): this 1-bit lock-free set does NOT
//override the base merge methods add(AbstractBitSet)/add(RawBitSet)/setToMax(AbstractBitSet). Those
//dispatch only on RawBitSet/MultiBitSet (else throw "Bad class") and cast `this` to RawBitSet - so calling
//them on an AtomicBitSet throws (RuntimeException or ClassCastException). The factory never returns this
//type; whole-set merges inherited from AbstractBitSet are unsupported for this atomic variant.
/** One-bit elements with atomic per-cell updates. Individual reads/ORs can overlap,
 * but clear() and cardinality() are cellwise operations, not atomic snapshots.
 * Construction, capacity changes and publication require external coordination.
 * Whole-set AbstractBitSet merges are unsupported; use the cell-level methods.
 * @author Brian Bushnell
 */
public class AtomicBitSet extends AbstractBitSet{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an empty atomic set; zero capacity leaves storage null.
	 * @param capacity_ Nonnegative logical capacity; rounded cell count must fit an int */
	public AtomicBitSet(long capacity_){
		setCapacity(capacity_, 0);
	}

	/** Creates an empty atomic set with optional spare cells.
	 * @param capacity_ Nonnegative logical bit capacity
	 * @param extra Nonnegative spare cells; rounded active cells plus extra must fit an int */
	public AtomicBitSet(long capacity_, int extra){
		setCapacity(capacity_, extra);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Atomic Updates       ----------------*/
	/*--------------------------------------------------------------*/

	/** Atomically ORs one word with retry on contention; an unchanged word needs no CAS.
	 * Concurrent per-cell updates are supported while storage and capacity remain fixed.
	 * @param cell Storage-cell index in [0, length())
	 * @param mask Bits to set; padding outside logical capacity must remain zero */
	@Override
	public void addToCell(final int cell, final int mask){
		int old=array.get(cell);
		int update=old|mask;
		while(update!=old && !array.compareAndSet(cell, old, update)){
			old=array.get(cell);
			update=old|mask;
		}
	}

	/** Takes elementwise maxima by atomically ORing one word.
	 * @param cell Storage-cell index in [0, length())
	 * @param mask Bits to set; padding outside logical capacity must remain zero */
	@Override
	public void setToMax(final int cell, final int mask){
		addToCell(cell, mask);
	}

	/** Atomically sets one bit; the amount does not matter beyond its positive precondition.
	 * @param x Logical bit index in [0, capacity())
	 * @param amt Positive increment; the count saturates at one */
	@Override
	public void increment(int x, int amt){
		assert(amt>0);
		//[map/AtomicBitSet#001] FIXED: was x<=capacity (off-by-one, shared family anomaly with Raw/MultiBitSet);
		//valid indices are [0,capacity-1]. 1-bit so amt is ignored (set-bit). Unreached today (0 callers).
		assert(x>=0 && x<capacity);
		final int cell=x/32;
		final int bit=x&31;
		final int mask=1<<bit;
		int old=array.get(cell);
		int update=old|mask;
		while(update!=old && !array.compareAndSet(cell, old, update)){
			old=array.get(cell);
			update=old|mask;
		}
	}

	/**
	 * Returns the value of the bit at position x.
	 * @param x Logical bit index in [0, capacity())
	 * @return 1 if the bit is set, 0 if clear
	 */
	@Override
	public int getCount(int x){
		assert(x>=0 && x<capacity); //[map/AtomicBitSet#001] FIXED: was x<=capacity (off-by-one, see increment)
		final int cell=x/32;
		final int bit=x&31;
		final int mask=1<<bit;
		int value=array.get(cell);
		return (value&mask)==mask ? 1 : 0;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Storage Lifecycle      ----------------*/
	/*--------------------------------------------------------------*/

	/** Clears each active cell atomically, leaving spare cells untouched.
	 * The whole clear is not atomic: concurrent updates can survive or be overwritten. */
	@Override
	public void clear(){
		for(int i=0; i<length; i++){
			array.set(i, 0);
		}
	}

	/**
	 * Returns the number of bits set to 1 in this bit set.
	 * Includes last-cell padding if callers have set it through the packed-cell API.
	 * @return Count observed while visiting active cells; not a concurrent snapshot
	 */
	@Override
	public long cardinality(){
		long sum=0;
		for(int i=0; i<length; i++){
			int value=array.get(i);
			sum+=Integer.bitCount(value);
		}
		return sum;
	}

	/**
	 * Sets the capacity of the bit set and reallocates array if needed.
	 * Calculates required array length as (capacity+31)/32 to accommodate all bits.
	 * Only reallocates if new capacity exceeds current maximum capacity. Reused
	 * storage is not cleared; replacement storage starts empty. Requires exclusive
	 * access and safe publication before concurrent readers or writers resume.
	 *
	 * Extra cells do not extend the reuse bound. A shrink does not clear hidden cells.
	 * @param capacity_ Nonnegative logical bit capacity
	 * @param extra Nonnegative spare cells; rounded active cells plus extra must fit an int
	 */
	@Override
	public void setCapacity(long capacity_, int extra){
		capacity=capacity_;
		length=(int)((capacity+31)/32);
		if(maxCapacity<capacity){
			maxLength=length+extra;
			maxCapacity=(long)length*32;//Wide arithmetic preserves reusable storage at large bit capacities.
			array=new AtomicIntegerArray(maxLength);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Accessors          ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns logical bit capacity, excluding spare storage.
	 * @return Current logical bit capacity */
	@Override
	public long capacity(){return capacity;}

	/** Returns active cell count, which can be smaller than the backing array length.
	 * @return Rounded-up logical capacity divided by 32 */
	@Override
	public int length(){return length;}

	/** Returns the number of bits per logical element.
	 * @return One */
	@Override
	public final int bits(){return 1;}

	/** Returns borrowed atomic storage, including spare cells; resizing may replace it.
	 * @return Backing array, or null before the first positive-capacity allocation */
	public AtomicIntegerArray array(){return array;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Logical bits covered by reusable storage, excluding extra cells. */
	private long maxCapacity=0;
	/** Current logical bit capacity. */
	private long capacity=0;
	/** Allocated backing cells including extra cells. */
	private int maxLength=0;
	/** Active cells for the current capacity. */
	private int length=0;
	/** Atomic packed cells; null before the first positive allocation. */
	private AtomicIntegerArray array;

}
