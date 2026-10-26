package map;

/** Unsynchronized, one-bit elements packed into an int array.
 * Capacity changes reuse storage when possible and do not clear reused bits.
 * A replacement allocation starts empty; callers must clear before logical reuse.
 * @author Brian Bushnell
 */
public class RawBitSet extends AbstractBitSet{

	/** Creates storage for the requested number of logical bits. */
	RawBitSet(long capacity_){
		setCapacity(capacity_, 0);
	}

	/** Creates storage with extra backing cells for later use. */
	RawBitSet(long capacity_, int extra){
		setCapacity(capacity_, extra);
	}
	
	/** ORs a packed word into an active storage cell. */
	@Override
	public void addToCell(final int cell, final int mask){
		int old=array[cell];
		int update=old|mask;
		array[cell]=update;
	}
	
	/** For one-bit elements, per-element maximum is a packed OR. */
	@Override
	public void setToMax(final int cell, final int mask){
		addToCell(cell, mask);
	}
	
	//1-bit set: `amt` is ignored (a bit can't count past 1) - increment == set-bit. assert(amt>0) just
	//rejects a no-op/negative caller. Saturating by nature.
	@Override
	public void increment(int x, int amt){
		assert(amt>0);
		//TODO: Possible bug [map/RawBitSet#001] FIXED - bound was `x<=capacity`, off by one: valid indices
		//are [0,capacity-1]; x==capacity reaches a phantom bit (last-cell slack -> miscounted by cardinality)
		//or AIOOBE when capacity%32==0. Shared family anomaly with MultiBitSet (same fix). -ea-only sharpening;
		//valid callers (Sketch.countMatches: i<a.length==capacity) already satisfy x<capacity.
		assert(x>=0 && x<capacity);
		final int cell=x/32;
		final int bit=x&31;
		final int mask=1<<bit;
		final int old=array[cell];
		final int update=old|mask;
		array[cell]=update;
	}
	
	/**
	 * Returns 1 if the bit at position x is set, 0 if it is clear.
	 * Efficiently checks a single bit using bitwise operations.
	 * @param x The bit position to check (0-based index)
	 * @return 1 if bit is set, 0 if bit is clear
	 */
	@Override
	public int getCount(int x){
		assert(x>=0 && x<capacity); //[map/RawBitSet#001] FIXED: was x<=capacity (off-by-one, see increment)
		final int cell=x/32;
		final int bit=x&31;
		final int mask=1<<bit;
		final int value=array[cell];
		return (value&mask)==mask ? 1 : 0;
	}
	
	/** Clears all bits in the set by zeroing all integer cells.
	 * Resets the entire bit set to an empty state. */
	@Override
	public void clear(){
		for(int i=0; i<length; i++){
			array[i]=0;
		}
	}
	
	/**
	 * Returns the number of set bits in the entire bit set.
	 * Uses Integer.bitCount() for efficient population count of each 32-bit cell.
	 * @return Total count of bits set to 1
	 */
	@Override
	public long cardinality(){
		long sum=0;
		for(int i=0; i<length; i++){
			int value=array[i];
			sum+=Integer.bitCount(value);
		}
		return sum;
	}
	
	/**
	 * Sets the capacity and allocates backing array if needed.
	 * Only reallocates if the new capacity exceeds the current maximum capacity.
	 * @param capacity_ New maximum number of bits to support
	 * @param extra Additional integer cells to allocate beyond minimum required
	 */
	@Override
	public void setCapacity(long capacity_, int extra){
		capacity=capacity_;
		length=(int)((capacity+31)/32); //cells = ceil(capacity/32)
		if(maxCapacity<capacity){ //only (re)allocate when growing past what we've ever held
			maxLength=length+extra;
			//Keep capacity arithmetic wide: int multiplication overflow previously forced needless reallocation.
			//Extra cells remain outside this conservative reuse bound. [map/RawBitSet#002]
			maxCapacity=(long)length*32;
			array=new int[maxLength];
		}
		//NOTE: a shrink/same-size reuse does NOT clear stale bits - caller must clear() before reuse.
	}

	/** Returns the maximum number of bits this set can hold */
	@Override
	public long capacity(){return capacity;}

	/** Returns the number of 32-bit integer cells currently in use */
	@Override
	public int length(){return length;}

	/** Returns 1, indicating this bit set stores 1 bit per element */
	@Override
	public final int bits(){return 1;}
	
	/** Returns borrowed backing storage, including spare cells; resizing may replace it. */
	public int[] array(){return array;}
	
	/** Logical bits covered by reusable storage, excluding extra cells. */
	private long maxCapacity=0;
	/** Current logical bit capacity. */
	private long capacity=0;
	/** Allocated storage cells, including extra cells. */
	private int maxLength=0;
	/** Active storage cells for the current logical capacity. */
	private int length=0;
	/** Packed words; null until a positive capacity first allocates storage. */
	private int[] array;

}
