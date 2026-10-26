package map;

import shared.Tools;

/** Unsynchronized packed counters with two bits per logical element (values 0–3).
 * Positive increments and packed additions saturate each element independently.
 * Capacity changes may replace storage with an empty array; reused storage is
 * not cleared. Callers own clearing and exclusive access during mutations.
 * @author Brian Bushnell
 */
public class MultiBitSet extends AbstractBitSet{

	/** Demonstrates per-element saturation; arguments are unused. */
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
	
	/** Allocates the requested number of logical counters. */
	MultiBitSet(long capacity_){
		setCapacity(capacity_, 0);
	}

	/** Allocates counters with extra backing cells. */
	MultiBitSet(long capacity_, int extra){
		setCapacity(capacity_, extra);
	}

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

	@Override
	public void increment(int x, int amt){
		assert(amt>0);
		//TODO: Possible bug [map/MultiBitSet#001] FIXED - bound was `x<=capacity`, off by one (valid indices
		//[0,capacity-1]; x==capacity -> phantom element or AIOOBE when capacity%elementsPerCell==0). Shared
		//family anomaly with RawBitSet (same fix). -ea-only; valid callers already satisfy x<capacity.
		assert(x>=0 && x<capacity);
		final int cell=x/elementsPerCell;
		final int element=x&modulus;
//		System.err.println(cell+", "+element);
		final int shift=bits*element;
//		System.err.println(shift);
		final int mask=elementMask<<shift;
//		System.err.println(mask);
		final int contents=array[cell];
//		System.err.println(contents);
		final int currentValue=(contents&mask)>>>shift;
//		System.err.println(currentValue);
		//Compare remaining capacity before adding; a large positive amount must not overflow into adjacent elements.
		final int newValue=(amt>=elementMask-currentValue ? elementMask : currentValue+amt);
//		System.err.println(newValue);
		if(newValue==currentValue){return;}
		final int update=((contents&~mask)|(newValue<<shift));
//		System.err.println(update);
		array[cell]=update;
	}

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

	/** Sets all elements to zero by clearing all cells in the underlying array. */
	@Override
	public void clear(){
		for(int i=0; i<length; i++){
			array[i]=0;
		}
	}

	/**
	 * Returns the number of non-zero elements in this bit set.
	 * Iterates through all cells and counts elements with non-zero values.
	 * @return Count of positions that contain non-zero values
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
	 * Only reallocates if new capacity exceeds current maximum capacity.
	 *
	 * @param capacity_ New capacity in number of elements
	 * @param extra Additional cells to allocate beyond required
	 */
	@Override
	public void setCapacity(long capacity_, int extra){
		capacity=capacity_;
		length=(int)((capacity+elementsPerCell-1)/elementsPerCell); //cells = ceil(capacity/elementsPerCell)
//		assert(false) : (capacity+elementsPerCell-1)+", "+elementsPerCell+", "+((capacity+elementsPerCell-1)/elementsPerCell);
		if(maxCapacity<capacity){
			maxLength=length+extra;
			//Wide arithmetic preserves same-capacity reuse above Integer.MAX_VALUE logical elements.
			maxCapacity=(long)length*elementsPerCell;
			array=new int[maxLength];
		}
//		assert(false) : capacity+", "+length;
	}

	/** Returns the current capacity of this bit set.
	 * @return Maximum number of elements this bit set can hold */
	@Override
	public long capacity(){return capacity;}

	/** Returns the number of cells in the underlying array.
	 * @return Number of 32-bit integers used to store the bit set */
	@Override
	public int length(){return length;}

	/** Returns the number of bits used per element.
	 * @return Always returns 2 for MultiBitSet */
	@Override
	public final int bits(){return bits;}
	
	/** Returns borrowed packed storage, including spare cells; resizing may replace it. */
	public int[] array(){return array;}
	
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
	
	/** Width of one saturating counter. */
	public static final int bits=2;
	/** Maximum counter value and mask for its low bits. */
	public static final int elementMask=~((-1)<<bits);
	/** Counters packed into one int. */
	public static final int elementsPerCell=32/bits;
	/** Mask selecting the counter index within a cell. */
	public static final int modulus=elementsPerCell-1;

}
