package structures;

import java.util.Arrays;

import shared.Timer;

/**
 * Circular buffer for storing long values using modulo arithmetic for any size.
 * Provides fixed-capacity storage with automatic overwriting of oldest values.
 * Uses standard modulo operations instead of bit masking for intuitive code.
 * Optimized for scenarios requiring arbitrary buffer sizes.
 *
 * @author Brian Bushnell
 * @contributor Isla (Highly-customized Claude instance)
 * @date May 8, 2025
 */
public final class RingBufferMod {
	
	/** Benchmark entry point: fills the buffer repeatedly and sums oldest values.
	 * Arguments: args[0]=size, args[1]=iterations. */
	public static void main(String[] args) {
		int size=Integer.parseInt(args[0]);
		long iters=Long.parseLong(args[1]), sum=0;
		Timer t=new Timer();
		RingBufferMod ring=new RingBufferMod(size);
		for(long i=0; i<iters; i++) {
			ring.add(i);
			sum+=ring.getOldestUnchecked();
		}
		t.stop("Sum="+sum);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates a ring buffer with specified fixed capacity.
	 * @param size_ The maximum number of elements the buffer can hold */
	public RingBufferMod(int size_) {
		size=size_;
		array=new long[size];
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Adds a value to the buffer, overwriting oldest value when full.
	 * Advances the insertion position and increments the total count.
	 * @param value The long value to store in the buffer
	 */
	public final void add(long value) {
		count++;
		array[pos]=value;
		pos=(pos+1)%size;
	}

	/**
	 * Gets the value at the current insertion position.
	 * This position contains the oldest value when buffer is full.
	 * @return The value at the current insertion position
	 */
	public final long getCurrent() {
		return array[pos];
	}

	/**
	 * Gets the most recently added value from the buffer.
	 * Calculates position using modulo arithmetic with wraparound handling.
	 * @return The value that was most recently added
	 */
	public final long getPrev() {
		return array[(pos-1+size)%size];
	}

	/**
	 * Gets the oldest value in the buffer with bounds checking.
	 * Returns first element if buffer not yet full, otherwise oldest element.
	 * Safe for use regardless of buffer fill state.
	 * @return The oldest value currently stored in the buffer
	 */
	public final long getOldest() {
		return (count<=size) ? array[0] : array[pos];
	}

	/**
	 * Gets the oldest value without bounds checking for performance.
	 * Assumes buffer has been filled with safe values. Use with caution.
	 * @return The oldest value at the current position
	 */
	public final long getOldestUnchecked() {
		return array[pos]; // Position of oldest when buffer is full
	}

	/**
	 * Gets a value at specified offset from the most recent value.
	 * Uses modulo arithmetic to handle wraparound correctly.
	 * @param offset Number of positions back from most recent (0=most recent)
	 * @return The value at the specified offset position
	 */
	public final long get(int offset) {
	    return array[(pos-1-offset+size)%size];
	}

	/**
	 * Fills the entire buffer array with a specified value.
	 * Useful for initialization or resetting buffer contents.
	 * @param value The value to fill all buffer positions with
	 */
	public final void fill(long value) {
		Arrays.fill(array, value);
	}
	
	/**
	 * Returns the number of valid elements currently in the buffer.
	 * Limited by buffer capacity, never exceeds the fixed size.
	 * @return The count of valid elements, capped at buffer capacity
	 */
	public final int size() {
		return (int)Math.min(count, size);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final long[] array;
	private final int size;
	
	private int pos=0;
	private long count=0;//Optional
}