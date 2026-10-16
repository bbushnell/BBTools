package structures;

import java.util.Arrays;

import shared.Timer;

/**
 * Fixed-capacity circular buffer for long values, using modulo arithmetic for
 * wraparound. Supported logical capacity is any positive size, and the physical
 * array has exactly that many slots.
 * Raw accessors do not validate empty reads or offsets; callers must either
 * tolerate the current backing value or restrict reads to written logical positions.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 8, 2025
 */
public final class RingBufferMod{
	
	/** Benchmark entry point: fills the buffer repeatedly and sums oldest values.
	 * Arguments: args[0]=size, args[1]=iterations. */
	public static void main(String[] args){
		int size=Integer.parseInt(args[0]);
		long iters=Long.parseLong(args[1]), sum=0;
		Timer t=new Timer();
		RingBufferMod ring=new RingBufferMod(size);
		for(long i=0; i<iters; i++){
			ring.add(i);
			sum+=ring.getOldestUnchecked();
		}
		t.stop("Sum="+sum);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Init             ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a ring buffer with the requested fixed capacity.
	 * @param size_ Positive logical and physical capacity.
	 */
	public RingBufferMod(int size_){
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
	public final void add(long value){
		count++;
		array[pos]=value;
		pos=(pos+1)%size;
	}

	/**
	 * Gets the value at the current insertion position.
	 * This position contains the oldest value when buffer is full.
	 * @return The value at the current insertion position
	 */
	public final long getCurrent(){
		return array[pos];
	}

	/**
	 * Gets the most recently added value from the buffer.
	 * Calculates position using modulo arithmetic with wraparound handling.
	 * @return The value that was most recently added
	 */
	public final long getPrev(){
		return array[(pos-1+size)%size];
	}

	/**
	 * Gets the raw oldest-slot value for the current fill state.
	 * Selects array[0] before a full logical cycle rather than the next slot to age
	 * out, but does not reject reads from an empty buffer.
	 * @return The oldest-slot value for the current fill state.
	 */
	public final long getOldest(){
		return (count<=size) ? array[0] : array[pos];
	}

	/**
	 * Gets the oldest value without bounds checking for performance.
	 * Assumes buffer has been filled with safe values. Use with caution.
	 * @return The oldest value at the current position
	 */
	public final long getOldestUnchecked(){
		return array[pos]; // Position of oldest when buffer is full
	}

	/**
	 * Gets a value at specified offset from the most recent value.
	 * Uses modulo arithmetic to handle wraparound correctly.
	 * @param offset Number of positions back from most recent (0=most recent)
	 * @return The value at the specified offset position
	 */
	public final long get(int offset){
		return array[(pos-1-offset+size)%size];
	}

	/**
	 * Fills the entire buffer array with a specified value.
	 * Useful for initialization or resetting buffer contents.
	 * @param value The value to fill all buffer positions with
	 */
	public final void fill(long value){
		Arrays.fill(array, value);
	}
	
	/**
	 * Returns the number of valid elements currently in the buffer.
	 * Limited by buffer capacity, never exceeds the fixed size.
	 * @return The count of valid elements, capped at buffer capacity
	 */
	public final int size(){
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