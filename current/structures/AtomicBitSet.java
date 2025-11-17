package structures;

import java.util.concurrent.atomic.AtomicIntegerArray;

/**
 * Thread-safe bit set implementation using atomic operations for concurrent access.
 * Extends AbstractBitSet to provide lock-free bit manipulation operations
 * suitable for high-performance multi-threaded applications.
 * Each bit position represents a boolean value with atomic read/write guarantees.
 *
 * @author Brian Bushnell
 * @date 2013
 */
public class AtomicBitSet extends AbstractBitSet {

	/** Constructs an AtomicBitSet with the specified capacity.
	 * @param capacity_ The number of bits this set can hold */
	public AtomicBitSet(long capacity_){
		setCapacity(capacity_, 0);
	}

	/**
	 * Constructs an AtomicBitSet with specified capacity plus extra space.
	 * The extra parameter allows pre-allocation of additional cells to avoid
	 * reallocation during growth.
	 *
	 * @param capacity_ The number of bits this set can hold
	 * @param extra Additional cells to pre-allocate beyond minimum required
	 */
	public AtomicBitSet(long capacity_, int extra){
		setCapacity(capacity_, extra);
	}

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
		assert(x>=0 && x<=capacity);
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

	@Override
	public int getCount(int x) {
		assert(x>=0 && x<=capacity);
		final int cell=x/32;
		final int bit=x&31;
		final int mask=1<<bit;
		int value=array.get(cell);
		return (value&mask)==mask ? 1 : 0;
	}

	@Override
	public void clear(){
		for(int i=0; i<length; i++){
			array.set(i, 0);
		}
	}

	@Override
	public long cardinality(){
		long sum=0;
		for(int i=0; i<length; i++){
			int value=array.get(i);
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
			array=new AtomicIntegerArray(maxLength);
		}
	}

	@Override
	public long capacity(){return capacity;}

	@Override
	public int length(){return length;}

	@Override
	public final int bits(){return 1;}
	
	/** Returns the underlying AtomicIntegerArray for direct access */
	public AtomicIntegerArray array(){return array;}
	
	/** Maximum capacity allocated in the underlying array */
	private long maxCapacity=0;
	/** Current capacity (number of bits this set can hold) */
	private long capacity=0;
	/** Maximum length allocated in the underlying array */
	private int maxLength=0;
	/** Current length of the underlying array */
	private int length=0;
	/** Underlying atomic integer array storing the bit data */
	private AtomicIntegerArray array;

}
