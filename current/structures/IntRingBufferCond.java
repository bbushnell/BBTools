package structures;

import java.util.Arrays;

import shared.Timer;

/**
 * A circular buffer of fixed size for storing int values.
 * <p>
 * This implementation uses conditionals instead of modulus.
 * 
 * @author Brian Bushnell
 * @date May 8, 2025
 */
public final class IntRingBufferCond {
	
	public static void main(String[] args) {
		int size=Integer.parseInt(args[0]);
		long iters=Long.parseLong(args[1]), sum=0;
		Timer t=new Timer();
		IntRingBufferCond ring=new IntRingBufferCond(size);
		for(long i=0; i<iters; i++) {
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
	 * @param size The fixed capacity of the buffer.
	 */
	public IntRingBufferCond(int size_) {
		size=size_;
		limit=size-1;
		array=new int[size];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds a value to the buffer, overwriting the oldest value if full.
	 * @param value The value to add.
	 * @return Evicted value or -1 if none.
	 */
	public final int add(int value) {
		count++;
		int old=(count>size ? array[pos] : -1);//Prefilling with -1 avoids conditional
		array[pos]=value;
		pos=(pos>=limit) ? 0 : pos+1;
//		System.err.println("Count="+count+", pos="+pos+", limit="+limit+", size="+size+", got "+value+", returning "+old);
//		assert(count<200);
		return old;
	}

	/**
	 * Gets the value at the current insertion position.
	 * @return The value at the current position.
	 */
	public final int getCurrent() {
		return array[pos];
	}

	/**
	 * Gets the most recently added value.
	 * @return The most recent value.
	 */
	public final int getPrev() {
		return array[(pos-1+size)%size];
	}

	/**
	 * Gets the oldest value in the buffer with bounds checking.
	 * Returns the first element if the buffer isn't yet full.
	 * @return The oldest value.
	 */
	public final int getOldest() {
		return (count<=size) ? array[0] : array[pos];
	}

	/**
	 * Gets the oldest value without bounds checking.
	 * Assumes the buffer has been pre-filled with safe values.
	 * @return The oldest value.
	 */
	public final int getOldestUnchecked() {
		return array[pos]; // Position of oldest when buffer is full
	}

	/**
	 * Gets a value at a specified offset from the most recent value.
	 * @param offset The number of positions back from the most recent value (0=most recent).
	 * @return The value at the specified offset.
	 */
	public final int get(int offset) {
	    return array[(pos-1-offset+size)%size];
	}

	/**
	 * Fills the entire buffer with a specified value.
	 * @param value The value to fill the buffer with.
	 */
	public final void fill(int value) {
		Arrays.fill(array, value);
	}
	
	/**
	 * Returns the number of valid elements in the buffer.
	 * @return The number of elements, limited by buffer capacity.
	 */
	public final int size() {
		return (int)Math.min(count, size);
	}
	
	/** @return True if the buffer is full. */
	public final boolean isFull() {
		return count>=size;
	}
	
	/** @return True if the buffer has overflowed. */
	public final boolean overflowed() {
		return count>size;
	}
	
	public void clear() {
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
	private int count=0;//Optional
}