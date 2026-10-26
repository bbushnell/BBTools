package map;

import shared.Tools;

/** Unsynchronized packed counters with two bits per logical element (values 0–3).
 * Positive increments and packed additions saturate each element independently.
 * Capacity changes may replace storage with an empty array; reused storage is
 * not cleared. Callers own clearing and exclusive access during mutations.
 * @author Brian Bushnell
 */
public class MultiBitSet extends AbstractBitSet{

	/** Demonstrates per-element saturation, printing successive states to stderr.
	 * @param args Unused */
	public static void main(String[] args){
		MultiBitSet mbs=new MultiBitSet(20);
		System.err.println(mbs);
		mbs.increment(5, 2);
		System.err.println(mbs);
		mbs.increment(5, 1);
		System.err.println(mbs);
		mbs.increment(1, 1);
		System.err.println(mbs);
		mbs.increment(2, 1);
		System.err.println(mbs);
		mbs.increment(7, 3);
		System.err.println(mbs);
		mbs.increment(8, 4);
		System.err.println(mbs);
		mbs.increment(15, 2);
		System.err.println(mbs);
		mbs.increment(15, 2);
		System.err.println(mbs);
		mbs.increment(5, 2);
		System.err.println(mbs);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates an empty set; zero capacity leaves the backing array null.
	 * @param capacity_ Nonnegative counter capacity; rounded cell count must fit an int */
	MultiBitSet(long capacity_){
		setCapacity(capacity_, 0);
	}

	/** Creates an empty set with optional spare backing cells.
	 * @param capacity_ Nonnegative logical counter capacity
	 * @param extra Nonnegative spare cells; rounded active cells plus extra must fit an int */
	MultiBitSet(long capacity_, int extra){
		setCapacity(capacity_, extra);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Packed Updates        ----------------*/
	/*--------------------------------------------------------------*/

	/** Adds 16 packed counts to one cell, saturating each count independently at three.
	 * @param cell Storage-cell index in [0, length())
	 * @param other Packed two-bit counts; padding outside logical capacity must remain zero */
	@Override
	public void addToCell(final int cell, final int other){
		if(other==0){return;}
		int a=array[cell], b=other, c=0;
		for(int i=0; i<elementsPerCell; i++){
			int av=a&elementMask;
			int bv=b&elementMask;
			int sum=Tools.min(av+bv, elementMask);
			int mask=sum<<(bits*i);
			c|=mask;
			a>>>=bits;
			b>>>=bits;
		}
		array[cell]=c;
	}

	/** Takes the elementwise maximum of 16 packed counts, not the maximum of the words.
	 * @param cell Storage-cell index in [0, length())
	 * @param other Packed two-bit counts; padding outside logical capacity must remain zero */
	@Override
	public void setToMax(final int cell, final int other){
		if(other==0){return;}
		int a=array[cell], b=other, c=0;
		for(int i=0; i<elementsPerCell; i++){
			int av=a&elementMask;
			int bv=b&elementMask;
			int max=Tools.max(av, bv);
			int mask=max<<(bits*i);
			c|=mask;
			a>>>=bits;
			b>>>=bits;
		}
		array[cell]=c;
	}

	/** Increments one counter, saturating at three without affecting neighboring counters.
	 * @param x Logical counter index in [0, capacity())
	 * @param amt Positive amount to add, including amounts greater than three */
	@Override
	public void increment(int x, int amt){
		assert(amt>0);
		//Fixed boundary [map/MultiBitSet#001]: x==capacity would access padding or the next cell.
		assert(x>=0 && x<capacity);
		final int cell=x/elementsPerCell;
		final int element=x&modulus;
		final int shift=bits*element;
		final int mask=elementMask<<shift;
		final int contents=array[cell];
		final int currentValue=(contents&mask)>>>shift;
		//Compare remaining capacity before adding; a large positive amount must not overflow into adjacent elements.
		final int newValue=(amt>=elementMask-currentValue ? elementMask : currentValue+amt);
		if(newValue==currentValue){return;}
		final int update=((contents&~mask)|(newValue<<shift));
		array[cell]=update;
	}

	/** Reads one logical counter without modifying storage.
	 * @param x Logical counter index in [0, capacity())
	 * @return Count in [0, 3] */
	@Override
	public int getCount(int x){
		assert(x>=0 && x<capacity); //[map/MultiBitSet#001] FIXED: was x<=capacity (off-by-one, see increment)
		final int cell=x/elementsPerCell;
		final int element=x&modulus;
		final int shift=bits*element;
		final int contents=array[cell];
		final int currentValue=(contents>>>shift)&elementMask;
		return currentValue;
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
	 * Returns the number of non-zero elements in this bit set.
	 * Includes nonzero last-cell padding; callers using packed-cell access must keep it zero.
	 * @return Number of nonzero packed elements, not the sum of their counts
	 */
	@Override
	public long cardinality(){
		long sum=0;
		for(int i=0; i<length; i++){
			int value=array[i];
			//counts NON-ZERO 2-bit elements (not their values). while(value!=0) early-exits once all remaining
			//high elements are zero; a non-zero high element keeps value!=0 until it shifts down and is counted.
			while(value!=0){
				if((value&elementMask)!=0){sum++;}
				value>>>=bits;
			}
		}
		return sum;
	}

	/**
	 * Sets the capacity and reallocates the underlying array if necessary.
	 * Calculates required cell count and allocates additional extra cells if specified.
	 * Only reallocates if new capacity exceeds the reuse bound. Replacement storage
	 * starts empty and discards prior counts; reused storage retains its contents.
	 * Shrinking does not clear hidden counters. Extra cells do not extend the reuse bound.
	 *
	 * @param capacity_ Nonnegative logical counter capacity
	 * @param extra Nonnegative spare cells; rounded active cells plus extra must fit an int
	 */
	@Override
	public void setCapacity(long capacity_, int extra){
		capacity=capacity_;
		length=(int)((capacity+elementsPerCell-1)/elementsPerCell); //cells = ceil(capacity/elementsPerCell)
		if(maxCapacity<capacity){
			maxLength=length+extra;
			//Wide arithmetic preserves same-capacity reuse above Integer.MAX_VALUE logical elements.
			maxCapacity=(long)length*elementsPerCell;
			array=new int[maxLength];
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Accessors          ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns logical counter capacity, excluding spare storage.
	 * @return Current logical counter capacity */
	@Override
	public long capacity(){return capacity;}

	/** Returns active cell count, which can be smaller than the backing array length.
	 * @return Rounded-up logical capacity divided by 16 */
	@Override
	public int length(){return length;}

	/** Returns the number of bits used per element.
	 * @return Always returns 2 for MultiBitSet */
	@Override
	public final int bits(){return bits;}

	/** Returns borrowed mutable storage, including spare cells; resizing may replace it.
	 * @return Backing array, or null before the first positive-capacity allocation */
	public int[] array(){return array;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Logical elements covered by reusable storage, excluding extra cells. */
	private long maxCapacity=0;
	/** Current logical element capacity. */
	private long capacity=0;
	/** Allocated storage cells, including extra cells. */
	private int maxLength=0;
	/** Active cells for the current capacity. */
	private int length=0;
	/** Packed storage; null before the first positive allocation. */
	private int[] array;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Width of one saturating counter. */
	public static final int bits=2;
	/** Maximum counter value and mask for its low bits. */
	public static final int elementMask=~((-1)<<bits);
	/** Counters packed into one int. */
	public static final int elementsPerCell=32/bits;
	/** Mask selecting the counter index within a cell. */
	public static final int modulus=elementsPerCell-1;

}
