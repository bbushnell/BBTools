package structures;

import java.util.Arrays;

import shared.Timer;

/**
 * Fixed-capacity circular buffer for long values, using a power-of-two physical
 * array and masks for wraparound. Supported logical capacity is {@code 1..1<<30};
 * the physical backing array has at least two slots and is rounded up to a power
 * of two.
 * Raw accessors do not validate empty reads or offsets; callers must either
 * tolerate the current backing value or restrict reads to written logical positions.
 * @author Brian Bushnell
 * @contributor Isla
 * @date May 8, 2025
 */
public final class RingBuffer{
	
	/** Benchmark entry point measuring add/get throughput.
	 * Args: [buffer_size] [iteration_count] */
	public static void main(String[] args){
		int size=Integer.parseInt(args[0]);
		long iters=Long.parseLong(args[1]), sum=0;
		Timer t=new Timer();
		RingBuffer ring=new RingBuffer(size);
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
	 * Creates a buffer with the requested logical capacity.
	 * The physical array has at least two slots and is rounded up to the next power
	 * of two.
	 * @param size_ Fixed logical capacity in the supported {@code 1..1<<30} range.
	 *        Values outside this range are rejected when assertions are enabled.
	 */
	public RingBuffer(int size_){
		assert(size_>=1 && size_<=(1<<30)) :
				"RingBuffer size="+size_+" is outside 1..1<<30; power-of-two "+
				"int backing sizes above 1<<30 cannot stay positive and Java shift distances wrap";
		size=size_;
		int bits=1;
		while(1<<bits<size){bits++;}
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
	public final void add(long value){
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
	 * Uses bit masking to wrap around buffer boundaries.
	 * @return The most recently stored value
	 */
	public final long getPrev(){
		return array[(pos-1)&mask];
	}

	/**
	 * Gets the raw oldest-slot value for the current fill state.
	 * Selects array[0] before a full logical cycle rather than the next slot to age
	 * out, but does not reject reads from an empty buffer.
	 * @return The oldest-slot value for the current fill state.
	 */
	public final long getOldest(){
		return array[(count<size) ? 0 : (pos-size)&mask];
	}

	/**
	 * Gets the raw value at the slot that will age out on the next add once the
	 * logical buffer has wrapped. Before the first full cycle, or after clear(),
	 * this can read backing-array values not written in the current cycle.
	 * @return The raw value at the oldest unchecked slot.
	 */
	public final long getOldestUnchecked(){
		return array[(pos-size)&mask];//Faster; be sure to pre-fill with a safe value
	}

	/**
	 * Gets a value at a specified offset from the most recent value.
	 * Offset 0 returns the most recent value, 1 returns the previous, etc.
	 * @param offset The number of positions back from the most recent value
	 * @return The value at the specified offset
	 */
	public final long get(int offset){
		return array[(pos-offset-1)&mask];
	}

	/**
	 * Fills the entire physical buffer with a specified value.
	 * Useful for initializing buffer to safe values for getOldestUnchecked().
	 * @param value The value to fill all buffer positions with
	 */
	public final void fill(long value){
		Arrays.fill(array, value);
	}
	
	/**
	 * Returns the number of valid elements currently in the buffer.
	 * Limited by the logical buffer capacity specified in constructor.
	 * @return The number of elements, capped at buffer capacity
	 */
	public final int size(){
		return (int)Math.min(count, size);
	}
	
	/** Clears buffer state (resets position; contents left as-is). */
	public void clear(){
		pos=0;
		count=0;
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
