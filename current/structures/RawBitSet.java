package structures;

public class RawBitSet extends AbstractBitSet {

	RawBitSet(long capacity_){
		setCapacity(capacity_, 0);
	}

	RawBitSet(long capacity_, int extra){
		setCapacity(capacity_, extra);
	}
	
	@Override
	public void addToCell(final int cell, final int mask){
		int old=array[cell];
		int update=old|mask;
		array[cell]=update;
	}
	
	@Override
	public void setToMax(final int cell, final int mask){
		addToCell(cell, mask);
	}
	
	@Override
	public void increment(int x, int amt){
		assert(amt>0);
		assert(x>=0 && x<=capacity);
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
		assert(x>=0 && x<=capacity);
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
		length=(int)((capacity+31)/32);
		if(maxCapacity<capacity){
			maxLength=length+extra;
			maxCapacity=length*32;
			array=new int[maxLength];
		}
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
	
	public int[] array(){return array;}
	
	private long maxCapacity=0;
	private long capacity=0;
	private int maxLength=0;
	private int length=0;
	private int[] array;

}
