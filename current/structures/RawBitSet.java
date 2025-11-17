package structures;

/**
 * Low-level bit manipulation structure with minimal memory overhead.
 * Stores one bit per element using packed integer arrays for efficient memory usage.
 * Provides basic set operations for tracking boolean states across large ranges.
 *
 * @author Brian Bushnell
 * @date 2014
 */
public class RawBitSet extends AbstractBitSet {

	/** Creates a RawBitSet with the specified capacity.
	 * @param capacity_ Maximum number of bits that can be stored */
	RawBitSet(long capacity_){
		setCapacity(capacity_, 0);
	}

	/**
	 * Creates a RawBitSet with the specified capacity and extra buffer space.
	 * @param capacity_ Maximum number of bits that can be stored
	 * @param extra Additional integer cells to allocate beyond minimum required
	 */
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
	
	@Override
	public int getCount(int x){
		assert(x>=0 && x<=capacity);
		final int cell=x/32;
		final int bit=x&31;
		final int mask=1<<bit;
		final int value=array[cell];
		return (value&mask)==mask ? 1 : 0;
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
			sum+=Integer.bitCount(value);
		}
		return sum;
	}
	
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

	@Override
	public long capacity(){return capacity;}

	@Override
	public int length(){return length;}

	@Override
	public final int bits(){return 1;}
	
	/** Returns direct access to the underlying integer array storing the bits */
	public int[] array(){return array;}
	
	/** Maximum capacity for which array space has been allocated */
	private long maxCapacity=0;
	/** Current maximum number of bits this set can hold */
	private long capacity=0;
	/** Length of the allocated integer array */
	private int maxLength=0;
	/** Number of integer cells currently needed for the current capacity */
	private int length=0;
	/** Backing array storing bits packed into 32-bit integers */
	private int[] array;

}
