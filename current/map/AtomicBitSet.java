package map;

import java.util.concurrent.atomic.AtomicIntegerArray;

//NOTE [map/AtomicBitSet#003] (LOW/latent, NOT patched - 0 callers): this 1-bit lock-free set does NOT
//override the base merge methods add(AbstractBitSet)/add(RawBitSet)/setToMax(AbstractBitSet). Those
//dispatch only on RawBitSet/MultiBitSet (else throw "Bad class") and cast `this` to RawBitSet - so calling
//them on an AtomicBitSet throws (RuntimeException or ClassCastException). The factory never returns this
//type; it's constructed directly and currently unused. If wired up, either override merge here or document
//that cross-set merge is unsupported for the atomic variant.
public class AtomicBitSet extends AbstractBitSet {

	public AtomicBitSet(long capacity_){
		setCapacity(capacity_, 0);
	}

	public AtomicBitSet(long capacity_, int extra){
		setCapacity(capacity_, extra);
	}

	//CLEVER [verified]: lock-free OR via CAS retry. `update!=old` short-circuits a no-op (mask adds no new
	//bits) so no CAS is issued; on contention it re-reads and retries until the OR lands. 1-bit: OR==union.
	@Override
	public void addToCell(final int cell, final int mask){
		int old=array.get(cell);
		int update=old|mask;
		while(update!=old && !array.compareAndSet(cell, old, update)){
			old=array.get(cell);
			update=old|mask;
		}
	}
	
	@Override
	public void setToMax(final int cell, final int mask){
		addToCell(cell, mask);
	}

	@Override
	public void increment(int x, int amt) {
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
	 * @param x The bit position to check (0-based)
	 * @return 1 if the bit is set, 0 if clear
	 */
	@Override
	public int getCount(int x) {
		assert(x>=0 && x<capacity); //[map/AtomicBitSet#001] FIXED: was x<=capacity (off-by-one, see increment)
		final int cell=x/32;
		final int bit=x&31;
		final int mask=1<<bit;
		int value=array.get(cell);
		return (value&mask)==mask ? 1 : 0;
	}

	/** Atomically clears all bits in the set to zero.
	 * Sets each cell in the underlying array to 0. */
	@Override
	public void clear(){
		for(int i=0; i<length; i++){
			array.set(i, 0);
		}
	}

	/**
	 * Returns the number of bits set to 1 in this bit set.
	 * Counts set bits across all cells using Integer.bitCount().
	 * @return The count of set bits in the entire bit set
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
	 * Only reallocates if new capacity exceeds current maximum capacity.
	 *
	 * @param capacity_ The new bit capacity
	 * @param extra Additional cells to allocate beyond minimum required
	 */
	@Override
	public void setCapacity(long capacity_, int extra){
		capacity=capacity_;
		length=(int)((capacity+31)/32);
		if(maxCapacity<capacity){
			maxLength=length+extra;
			maxCapacity=length*32; //NOTE (latent LOW, =RawBitSet#002): int overflow for length>~67M; should be (long)length*32
			array=new AtomicIntegerArray(maxLength);
		}
	}

	/** Returns the maximum number of bits this set can hold */
	@Override
	public long capacity(){return capacity;}

	/** Returns the length of the underlying integer array */
	@Override
	public int length(){return length;}

	/** Returns the number of bits per element (always 1 for bit sets) */
	@Override
	public final int bits(){return 1;}
	
	public AtomicIntegerArray array(){return array;}
	
	private long maxCapacity=0;
	private long capacity=0;
	private int maxLength=0;
	private int length=0;
	private AtomicIntegerArray array;

}
