package structures;

import shared.Tools;

/**
 * A space-efficient bit set implementation that stores multi-bit values (0-3) for each position.
 * Uses 2 bits per element, packing 16 elements into each 32-bit integer for memory efficiency.
 * Supports increment operations with saturation at maximum value and various set operations.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class MultiBitSet extends AbstractBitSet {

	/**
	 * Test method demonstrating MultiBitSet functionality.
	 * Creates a bit set and performs various increment operations to show behavior.
	 * @param args Command-line arguments (unused)
	 */
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
	
	/** Creates a MultiBitSet with the specified capacity and no extra space.
	 * @param capacity_ Maximum number of elements this bit set can hold */
	MultiBitSet(long capacity_){
		setCapacity(capacity_, 0);
	}

	/**
	 * Creates a MultiBitSet with the specified capacity plus extra space.
	 * Extra space allows for expansion without reallocation.
	 * @param capacity_ Maximum number of elements this bit set can hold
	 * @param extra Additional cells to allocate beyond required capacity
	 */
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
		assert(x>=0 && x<=capacity);
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
		final int newValue=Tools.min(currentValue+amt, elementMask);
//		System.err.println(newValue);
		if(newValue==currentValue){return;}
		final int update=((contents&~mask)|(newValue<<shift));
//		System.err.println(update);
		array[cell]=update;
	}

	@Override
	public int getCount(int x){
		assert(x>=0 && x<=capacity);
		final int cell=x/elementsPerCell;
		final int element=x&modulus;
		final int shift=bits*element;
		final int contents=array[cell];
		final int currentValue=(contents>>>shift)&elementMask;
		return currentValue;
	}

	@Override
	public void clear(){
		for(int i=0; i<length; i++){
			array[i]=0;
		}
	}

	@Override
	public long cardinality(){
		long sum=0;
		for(int i=0; i<length; i++){
			int value=array[i];
			while(value!=0){
				if((value&elementMask)!=0){sum++;}
				value>>>=bits;
			}
		}
		return sum;
	}

	@Override
	public void setCapacity(long capacity_, int extra){
		capacity=capacity_;
		length=(int)((capacity+elementsPerCell-1)/elementsPerCell);
//		assert(false) : (capacity+elementsPerCell-1)+", "+elementsPerCell+", "+((capacity+elementsPerCell-1)/elementsPerCell);
		if(maxCapacity<capacity){
			maxLength=length+extra;
			maxCapacity=length*elementsPerCell;
			array=new int[maxLength];
		}
//		assert(false) : capacity+", "+length;
	}

	@Override
	public long capacity(){return capacity;}

	@Override
	public int length(){return length;}

	@Override
	public final int bits(){return bits;}
	
	/** Returns direct access to the underlying storage array.
	 * @return The array of integers storing the packed bit set data */
	public int[] array(){return array;}
	
	/** Maximum capacity for which space has been allocated */
	private long maxCapacity=0;
	/** Current capacity in number of elements */
	private long capacity=0;
	/** Maximum length of allocated array */
	private int maxLength=0;
	/** Current length of array in use */
	private int length=0;
	/** Underlying storage array where each int holds 16 packed 2-bit elements */
	private int[] array;
	
	/** Number of bits per element (always 2 for MultiBitSet) */
	public static final int bits=2;
	/** Bit mask for extracting a single element value (0b11 = 3) */
	public static final int elementMask=~((-1)<<bits);
	/** Number of elements stored in each 32-bit cell (16 for 2-bit elements) */
	public static final int elementsPerCell=32/bits;
	/** Modulus for computing element index within a cell (elementsPerCell - 1) */
	public static final int modulus=elementsPerCell-1;

}
