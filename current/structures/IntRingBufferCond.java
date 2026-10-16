package structures;

import java.util.Arrays;

import shared.Timer;

/**
 * Fixed-capacity circular buffer for int values, using conditionals for
 * wraparound. Raw accessors do not validate empty reads or offsets; callers
 * must prefill the buffer or restrict reads to written positions.
 * 
 * @author Brian Bushnell
 * @date May 8, 2025
 */
public final class IntRingBufferCond{
	
	public static void main(String[] args){
		int size=Integer.parseInt(args[0]);
		long iters=Long.parseLong(args[1]), sum=0;
		Timer t=new Timer();
		IntRingBufferCond ring=new IntRingBufferCond(size);
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
	 * Creates a buffer of specified size.
	 * @param size_ Positive fixed capacity of the buffer.
	 */
	public IntRingBufferCond(int size_){
		size=size_;
		limit=size-1;
		array=new int[size];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds a value to the buffer, overwriting the oldest value if full.
	 * The return value uses -1 as a no-eviction sentinel; callers that test for
	 * old>=0 must use nonnegative payloads.
	 * @param value The value to add.
	 * @return Evicted value or -1 if none.
	 */
	public final int add(int value){
		count++;
		int old=(count>size ? array[pos] : -1);//Prefilling with -1 avoids conditional
		array[pos]=value;
		pos=(pos>=limit) ? 0 : pos+1;
//		System.err.println("Count="+count+", pos="+pos+", limit="+limit+", size="+size+", got "+value+", returning "+old);
//		assert(count<200);
		return old;
	}

	/**
	 * Gets the raw value at the next insertion position.
	 * This slot may be unwritten unless the buffer has been filled or wrapped.
	 * @return The value at the next insertion position.
	 */
	public final int getCurrent(){
		return array[pos];
	}

	/**
	 * Gets the most recently added value.
	 * The returned slot is raw and may be unwritten when the buffer is empty.
	 * @return The most recent value.
	 */
	public final int getPrev(){
		return array[(pos-1+size)%size];
	}

	/**
	 * Gets the oldest logical value.
	 * Returns the first slot until the buffer has filled.
	 * @return The oldest logical value.
	 */
	public final int getOldest(){
		return (count<=size) ? array[0] : array[pos];
	}

	/**
	 * Gets the slot that is oldest once the buffer is full.
	 * Assumes the buffer has been filled or pre-filled with safe values.
	 * @return The oldest full-buffer value.
	 */
	public final int getOldestUnchecked(){
		return array[pos]; // Position of oldest when buffer is full
	}

	/**
	 * Gets a value at a specified offset from the most recent value.
	 * @param offset Positions back from the most recent value; must be within
	 *        written or pre-filled contents, with 0 as most recent.
	 * @return The value at the specified offset.
	 */
	public final int get(int offset){
		return array[(pos-1-offset+size)%size];
	}

	/**
	 * Fills the entire buffer with a specified value.
	 * Does not change the insertion position or valid size.
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
	
	/** @return True if the buffer is full. */
	public final boolean isFull(){
		return count>=size;
	}
	
	/** @return True if the buffer has overflowed. */
	public final boolean overflowed(){
		return count>size;
	}
	
	/** Resets the insertion position and valid size without clearing contents. */
	public void clear(){
		pos=0;
		count=0;
		//fill(-1); //Not needed, but can be faster
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int[] array;
	private final int size;
	private final int limit;
	
	private int pos=0;
	private long count=0;//Optional
}
