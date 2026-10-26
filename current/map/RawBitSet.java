package map;

/** Unsynchronized, one-bit elements packed into an int array.
 * Capacity changes reuse storage when possible and do not clear reused bits.
 * A replacement allocation starts empty; callers must clear before logical reuse.
 * @author Brian Bushnell
 */
public class RawBitSet extends AbstractBitSet{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an empty set; zero capacity leaves the backing array null.
	 * @param capacity_ Nonnegative logical bit capacity; rounded cell count must fit an int */
	RawBitSet(long capacity_){
		setCapacity(capacity_, 0);
	}

	/** Creates an empty set with optional spare storage; extra does not extend the reuse bound.
	 * @param capacity_ Nonnegative logical bit capacity
	 * @param extra Nonnegative extra cells; rounded cell count plus extra must fit an int */
	RawBitSet(long capacity_, int extra){
		setCapacity(capacity_, extra);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Updates        ----------------*/
	/*--------------------------------------------------------------*/

	/** ORs a packed word into an active storage cell, leaving other cells unchanged.
	 * @param cell Storage-cell index in [0, length())
	 * @param mask Bits to set; padding outside logical capacity must remain zero */
	@Override
	public void addToCell(final int cell, final int mask){
		int old=array[cell];
		int update=old|mask;
		array[cell]=update;
	}

	/** Takes elementwise maxima by ORing one packed word.
	 * @param cell Storage-cell index in [0, length())
	 * @param mask Bits to set; padding outside logical capacity must remain zero */
	@Override
	public void setToMax(final int cell, final int mask){
		addToCell(cell, mask);
	}

	/** Sets one bit to one; a positive amount has no effect beyond setting the bit.
	 * @param x Logical bit index in [0, capacity())
	 * @param amt Positive increment; the count saturates at one */
	@Override
	public void increment(int x, int amt){
		assert(amt>0);
		//Fixed boundary [map/RawBitSet#001]: x==capacity would access padding or the next cell.
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
	 * @param x Logical bit index in [0, capacity())
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

	/*--------------------------------------------------------------*/
	/*----------------       Storage Lifecycle      ----------------*/
	/*--------------------------------------------------------------*/

	/** Zeroes active cells and their padding; spare cells beyond length() are untouched. */
	@Override
	public void clear(){
		for(int i=0; i<length; i++){
			array[i]=0;
		}
	}

	/**
	 * Counts set bits in active cells, including any nonzero last-cell padding.
	 * @return Number of set bits; equals logical cardinality when padding is zero
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
	 * Changes logical capacity, retaining old cells when the reuse bound permits.
	 * Otherwise allocates an empty replacement, discarding prior contents. A shrink
	 * does not clear hidden bits: clear before reuse if those bits must not reappear.
	 * Extra cells are allocated only on replacement and do not extend the reuse bound.
	 * @param capacity_ Nonnegative logical bit capacity
	 * @param extra Nonnegative spare cells; rounded active cells plus extra must fit an int
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

	/*--------------------------------------------------------------*/
	/*----------------           Accessors          ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the logical bit capacity, excluding spare storage.
	 * @return Current logical bit capacity */
	@Override
	public long capacity(){return capacity;}

	/** Returns the number of active cells, including the partially used final cell.
	 * @return Rounded-up logical capacity divided by 32 */
	@Override
	public int length(){return length;}

	/** Returns the number of bits per logical element.
	 * @return One */
	@Override
	public final int bits(){return 1;}

	/** Returns borrowed mutable storage, including spare cells; resizing may replace it.
	 * @return Backing array, or null before the first positive-capacity allocation */
	public int[] array(){return array;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

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
