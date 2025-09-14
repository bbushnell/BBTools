package structures;

import java.util.Arrays;

import shared.Timer;

/**
 * A circular buffer of fixed size for storing long values using modulo arithmetic.
 * Uses masks for speed, but allows arbitrary logical buffer size.
 * The physical buffer size is rounded up to the next power of two.
 *@author Brian Bushnell
 *@contributor Isla (Highly-customized Claude instance)
 *@date May 8, 2025
 */
public final class RingBuffer {
	
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

	/**
	 * Creates a buffer of specified size.
	 * @param size The fixed capacity of the buffer.
	 */
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
	 * @param value The value to add.
	 */
	public final void add(long value) {
		count++;
		array[pos]=value;
		pos=(pos+1)&mask;
	}

	/**
	 * Gets the value at the current insertion position.
	 * @return The value at the current position.
	 */
	public final long getCurrent() {
		return array[pos];
	}

	/**
	 * Gets the most recently added value.
	 * @return The most recent value.
	 */
	public final long getPrev() {
		return array[(pos-1)&mask];
	}

	/**
	 * Gets the oldest value in the buffer with bounds checking.
	 * Returns the first element if the buffer isn't yet full.
	 * @return The oldest value.
	 */
	public final long getOldest() {
		return array[(count<size) ? 0 : (pos-size)&mask];
	}

	/**
	 * Gets the oldest value without bounds checking.
	 * Assumes the buffer has been pre-filled with safe values.
	 * @return The oldest value.
	 */
	public final long getOldestUnchecked() {
		return array[(pos-size)&mask];//Faster; be sure to pre-fill with a safe value
	}

	/**
	 * Gets a value at a specified offset from the most recent value.
	 * @param offset The number of positions back from the most recent value (0=most recent).
	 * @return The value at the specified offset.
	 */
	public final long get(int offset) {
		return array[(pos-offset-1)&mask];
	}

	/**
	 * Fills the entire buffer with a specified value.
	 * @param value The value to fill the buffer with.
	 */
	public final void fill(long value) {
		Arrays.fill(array, value);
	}
	
	/**
	 * Returns the number of valid elements in the buffer.
	 * @return The number of elements, limited by buffer capacity.
	 */
	public final int size() {
		return (int)Math.min(count, mask+1);
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
