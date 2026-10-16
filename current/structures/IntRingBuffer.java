package structures;

import java.util.Arrays;

import shared.Timer;

/**
 * Fixed-capacity circular buffer for int values, using a power-of-two physical
 * array and masks for wraparound. Supported logical capacity is {@code 1..1<<30}; the
 * physical backing array has at least two slots and is rounded up to a power of
 * two.
 * Raw accessors do not validate empty reads or offsets; callers must prefill the
 * buffer or restrict reads to written logical positions.
 *@author Brian Bushnell
 *@date October 2, 2025
 */
public final class IntRingBuffer{
	
	public static void main(String[] args){
		int size=Integer.parseInt(args[0]);
		long iters=Long.parseLong(args[1]), sum=0;
		Timer t=new Timer();
		IntRingBuffer ring=new IntRingBuffer(size);
		for(long i=0; i<iters; i++){
			ring.add((int)i);
			sum+=ring.getOldestUnchecked();
		}
		t.stop("Sum="+sum);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a buffer with the requested logical capacity.
	 * The physical array has at least two slots and is rounded up to the next power
	 * of two.
	 * @param size_ Fixed logical capacity in the supported {@code 1..1<<30} range.
	 *        Values outside this range are rejected when assertions are enabled.
	 */
	public IntRingBuffer(int size_){
		assert(size_>=1 && size_<=(1<<30)) :
				"IntRingBuffer size="+size_+" is outside 1..1<<30; power-of-two "+
				"int backing sizes above 1<<30 cannot stay positive and Java shift distances wrap";
		size=size_;
		int bits=1;
		while(1<<bits<size){bits++;}
		array=new int[1<<bits];
		mask=array.length-1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds a value to the buffer, overwriting the oldest value if full.
	 * @param value The value to add.
	 */
	public final void add(int value){
		count++;
		array[pos]=value;
		pos=(pos+1)&mask;
	}

	/**
	 * Gets the raw value at the next insertion position.
	 * This slot may be unwritten unless the buffer has been filled or wrapped.
	 * @return The value at the next insertion position.
	 */
	public final long getCurrent(){
		return array[pos];
	}

	/**
	 * Gets the most recently added value.
	 * The returned slot is raw and may be unwritten when the buffer is empty.
	 * @return The most recent value.
	 */
	public final long getPrev(){
		return array[(pos-1)&mask];
	}

	/**
	 * Gets the oldest logical value.
	 * Returns the first slot until the logical buffer has filled.
	 * @return The oldest logical value.
	 */
	public final long getOldest(){
		return array[(count<size) ? 0 : (pos-size)&mask];
	}

	/**
	 * Gets the slot that will be oldest once the logical buffer is full.
	 * Assumes the buffer has been filled or pre-filled with safe values.
	 * @return The oldest full-buffer value.
	 */
	public final long getOldestUnchecked(){
		return array[(pos-size)&mask];//Faster; be sure to pre-fill with a safe value
	}

	/**
	 * Gets a value at a specified offset from the most recent value.
	 * @param offset Positions back from the most recent value; must be within
	 *        written or pre-filled logical contents, with 0 as most recent.
	 * @return The value at the specified offset.
	 */
	public final long get(int offset){
		return array[(pos-offset-1)&mask];
	}

	/**
	 * Fills the entire buffer with a specified value.
	 * Does not change the insertion position or valid logical size.
	 * @param value The value to fill the buffer with.
	 */
	public final void fill(int value){
		Arrays.fill(array, value);
	}
	
	/**
	 * Returns the number of valid elements in the buffer.
	 * @return The number of elements, limited by buffer capacity.
	 */
	public final int size(){
		return (int)Math.min(count, size);
	}
	
	/** Resets the insertion position and valid size without clearing contents. */
	public void clear(){
		pos=0;
		count=0;
		//fill(0); //Not needed
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int[] array;
	private final int mask;
	private final int size;
	
	private int pos=0;
	private long count=0;//Optional
	
}
