package structures;

import java.util.Arrays;

import shared.Timer;

/**
 * Fixed-size circular buffer for longs using power-of-two physical storage and bit masking.
 * Supports arbitrary logical sizes with fast add/get operations.
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 8, 2025
 */
public final class RingBuffer {
	
	/** Benchmark entry point measuring add/get throughput.
	 * Args: [buffer_size] [iteration_count] */
	public static void main(String[] args) {
		int size=Integer.parseInt(args[0]);
		long iters=Long.parseLong(args[1]), sum=0;
		Timer t=new Timer();
		RingBuffer ring=new RingBuffer(size);
		for(long i=0; i<iters; i++) {
			ring.add(i);
			sum+=ring.getOldestUnchecked();
		}
		t.stop("Sum="+sum);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a buffer of specified logical size; physical array is rounded to next power of two.
	 * @param size_ Logical capacity of the buffer */
	public RingBuffer(int size_) {
		size=size_;
		int bits=1;
		while(1<<bits<size) {bits++;}
		array=new long[1<<bits];
		mask=array.length-1;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds a value to the buffer, overwriting the oldest value if full.
	 * Uses bit masking for fast circular wraparound.
	 * @param value The value to store in the buffer
	 */
	public final void add(long value) {
		count++;
		array[pos]=value;
		pos=(pos+1)&mask;
	}

	/**
	 * Gets the value at the current insertion position.
	 * Returns the value that would be overwritten by the next add().
	 * @return The value at the current write position
	 */
	public final long getCurrent() {
		return array[pos];
	}

	/**
	 * Gets the most recently added value.
	 * Uses bit masking to wrap around buffer boundaries.
	 * @return The most recently stored value
	 */
	public final long getPrev() {
		return array[(pos-1)&mask];
	}

	/**
	 * Gets the oldest value in the buffer with bounds checking.
	 * Returns the first element if the buffer isn't yet full.
	 * Safe version of getOldestUnchecked().
	 * @return The oldest value in the buffer
	 */
	public final long getOldest() {
		return array[(count<size) ? 0 : (pos-size)&mask];
	}

	/**
	 * Gets the oldest value without bounds checking.
	 * Assumes the buffer has been pre-filled with safe values.
	 * Faster than getOldest() but requires careful initialization.
	 * @return The oldest value in the buffer
	 */
	public final long getOldestUnchecked() {
		return array[(pos-size)&mask];//Faster; be sure to pre-fill with a safe value
	}

	/**
	 * Gets a value at a specified offset from the most recent value.
	 * Offset 0 returns the most recent value, 1 returns the previous, etc.
	 * @param offset The number of positions back from the most recent value
	 * @return The value at the specified offset
	 */
	public final long get(int offset) {
		return array[(pos-offset-1)&mask];
	}

	/**
	 * Fills the entire physical buffer with a specified value.
	 * Useful for initializing buffer to safe values for getOldestUnchecked().
	 * @param value The value to fill all buffer positions with
	 */
	public final void fill(long value) {
		Arrays.fill(array, value);
	}
	
	/**
	 * Returns the number of valid elements currently in the buffer.
	 * Limited by the logical buffer capacity specified in constructor.
	 * @return The number of elements, capped at buffer capacity
	 */
	public final int size() {
		return (int)Math.min(count, size);
	}
	
	/** Clears buffer state (resets position; contents left as-is). */
	public void clear() {
		pos=0;
		//fill(0); //Not needed
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final long[] array;
	private final int mask;
	private final int size;
	
	private int pos=0;
	private long count=0;//Optional
	
}
