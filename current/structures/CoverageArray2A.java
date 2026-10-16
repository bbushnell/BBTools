package structures;
import java.util.concurrent.atomic.AtomicIntegerArray;

import shared.KillSwitch;
import shared.Tools;

/**
 * Atomic CoverageArray implementation with two unsigned 16-bit coverage values
 * packed per integer. Individual updates use compare-and-set loops and saturate
 * positive increments above 0xFFFF. Inherited aggregate-statistic caches are not
 * updated atomically with these element updates.
 *
 * @author Brian Bushnell
 * @date Sep 20, 2014
 */
public class CoverageArray2A extends CoverageArray{

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Serial version UID for serialization compatibility. */
	private static final long serialVersionUID=98483952072098494L;

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		//TODO
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public CoverageArray2A(int chrom, int len){
		super(chrom, len);
		int intLen=intIdx(len)+1;
		array=KillSwitch.allocAtomicInt(intLen);
		minIndex=0;
		maxIndex=len-1;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the packed integer index containing a logical coverage index. */
	private static final int intIdx(int idx){
		return (idx+1)/2;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Atomically increments coverage at the specified location by 1.
	 * @param loc Location to increment.
	 */
	@Override
	public void increment(int loc){
		increment(loc, 1);
	}

	/**
	 * Atomically increments coverage at the specified location by the given amount.
	 * Uses compare-and-exchange operations to ensure thread safety.
	 * Coverage values are capped at 0xFFFF to prevent overflow.
	 *
	 * @param loc Location to increment
	 * @param amt Amount to increment by (must be non-negative)
	 */
	@Override
	public void increment(final int loc, final int amt){
		assert(amt>=0) : "This does not currently allow negative increments.";
		if(amt==0){return;}
		final int intIdx=intIdx(loc);
		boolean overflow=((loc&1)==1 ? incrementUpper(intIdx, amt) : incrementLower(intIdx, amt));
		if(overflow && !OVERFLOWED){
			 System.err.println("Note: Coverage capped at "+0xFFFF);
			 OVERFLOWED=true;
		}
	}

	/** Saturating compare-and-set increment of the lower packed 16-bit value. */
	private boolean incrementLower(final int intIdx, final int amt){
		for(;;){
			final int oldVal=array.get(intIdx);
			int lower=oldVal&lowerMask, upper=oldVal&upperMask;
			int charVal=lower;
			// Widen before adding so saturation sees sums above Integer.MAX_VALUE.
			int charVal2=(int)Tools.min(0xFFFF, (long)charVal+amt);
			boolean overflow=(charVal2>=0xFFFF);
			int newVal=(charVal2)|upper;
			if(array.compareAndSet(intIdx, oldVal, newVal)){return overflow;}
		}
	}

	/** Saturating compare-and-set increment of the upper packed 16-bit value. */
	private boolean incrementUpper(int intIdx, int amt){
		for(;;){
			final int oldVal=array.get(intIdx);
			int lower=oldVal&lowerMask, upper=oldVal&upperMask;
			int charVal=upper>>>16;
			// Widen before adding so saturation sees sums above Integer.MAX_VALUE.
			int charVal2=(int)Tools.min(0xFFFF, (long)charVal+amt);
			boolean overflow=(charVal2>=0xFFFF);
			int newVal=(charVal2<<16)|lower;
			if(array.compareAndSet(intIdx, oldVal, newVal)){return overflow;}
		}
	}

	/**
	 * Increments coverage values over a range of positions.
	 * For atomic arrays, synchronization is not needed so this calls incrementRange directly.
	 *
	 * @param min Starting position (inclusive)
	 * @param max Ending position (inclusive)
	 * @param amt Amount to increment by
	 */
	@Override
	public void incrementRangeSynchronized(int min, int max, int amt){
		incrementRange(min, max, amt);//Synchronized is not needed
	}

	/**
	 * Increments a range by visiting each position.
	 * @param min Starting position; negative starts are clamped to 0.
	 * @param max Ending position; positions above maxIndex are clamped.
	 * @param amt Amount to increment by.
	 */
	public void incrementRangeSlow(int min, int max, int amt){
		if(min<0){min=0;}
		if(max>maxIndex){max=maxIndex;}
		for(int loc=min; loc<=max; loc++){
			increment(loc, amt);
		}
	}

	/**
	 * Atomically increments coverage values over a range of positions.
	 * Currently uses per-position increments because the packed fast path below
	 * has known boundary issues.
	 *
	 * @param min Starting position (inclusive)
	 * @param max Ending position (inclusive)
	 * @param amt Amount to increment by
	 */
	@Override
	public void incrementRange(int min, int max, int amt){
		if(amt>0xFFF || true){
			incrementRangeSlow(min, max, amt);
			return;
		}
		//TODO:  This should be 2x as fast, but currently gives slightly wrong results.
		//Off by ~1% so probably a boundary issue
		//Try printing range and aborting
		if(max>maxIndex){max=maxIndex;}
		if((min&1)==1){increment(min, amt);}
		if((max&1)==0){increment(max, amt);}
		int minIdx=intIdx(min+1);
		int maxIdx=intIdx(max-1);
		for(int i=minIdx; i<=maxIdx; i++){
			for(;;){
				final int oldVal=array.get(i);
				int lower=oldVal&lowerMask, upper=(oldVal&upperMask)>>16;
				int lower2=Tools.min(0xFFFF, lower+amt);
				int upper2=Tools.min(0xFFFF, upper+amt);
				int newVal=lower2|(upper2<<16);
				if(array.compareAndSet(i, oldVal, newVal)){break;}
			}
		}
	}

	/**
	 * Atomically sets the coverage value at the specified location.
	 * Value is capped at 0xFFFF to prevent overflow of 16-bit storage.
	 * @param loc Location to set
	 * @param val0 Value to set (will be capped at 0xFFFF)
	 */
	@Override
	public void set(int loc, int val0){
		assert(val0>=0) : "This does not currently allow negative values.";
		final int intIdx=intIdx(loc);
		final int val=Tools.min(val0, 0xFFFF);
		boolean overflow=((loc&1)==1 ? setUpper(intIdx, val) : setLower(intIdx, val));
		if(val0!=val && !OVERFLOWED){
			 System.err.println("Note: Coverage capped at "+0xFFFF);
			 OVERFLOWED=true;
		}
	}

	/** Compare-and-set write of the lower packed 16-bit value. */
	private boolean setLower(final int intIdx, final int amt){
		for(;;){
			final int oldVal=array.get(intIdx);
			int lower=amt, upper=oldVal&upperMask;
			int newVal=lower|upper;
			if(array.compareAndSet(intIdx, oldVal, newVal)){return false;}
		}
	}

	/** Compare-and-set write of the upper packed 16-bit value. */
	private boolean setUpper(int intIdx, int amt){
		for(;;){
			final int oldVal=array.get(intIdx);
			int lower=oldVal&lowerMask, upper=amt<<16;
			int newVal=lower|upper;
			if(array.compareAndSet(intIdx, oldVal, newVal)){return false;}
		}
	}

	/**
	 * Atomically retrieves the coverage value at the specified location.
	 * Extracts the appropriate 16-bit value from the packed integer storage.
	 * Returns 0 when the computed packed integer index is outside the backing
	 * array bounds.
	 *
	 * @param loc Location to retrieve coverage from
	 * @return Coverage value at the specified location (0-65535)
	 */
	@Override
	public int get(int loc){
		final int intIdx=intIdx(loc);
		final int intVal=intIdx<0 || intIdx>=array.length() ? 0 : array.get(intIdx);
		return (loc&1)==1 ? (intVal>>>16) : (intVal&lowerMask);
	}

	/**
	 * Resize operation is not supported for atomic coverage arrays.
	 * @param newlen New length (ignored)
	 * @throws RuntimeException Always thrown as resize is unsupported
	 */
	@Override
	public void resize(int newlen){
		throw new RuntimeException("Resize: Unsupported.");
	}

	/**
	 * Creates a string representation of the coverage array.
	 * Returns coverage values in comma-separated format within square brackets.
	 * @return String representation of the coverage array
	 */
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		for(int i=0; i<length(); i++){
			if(i>0){sb.append(", ");}
			sb.append(get(i));
		}
		sb.append(']');
		return sb.toString();
	}

	/**
	 * Converts the coverage array to a character array.
	 * Each logical coverage value is cast to a char.
	 * @return Character array containing coverage values
	 */
	@Override
	public char[] toArray(){return copyToArray();}

	/** Returns independent writable storage for median sorting.
	 * The copy is assembled element by element, not as one array-wide atomic snapshot.
	 */
	@Override
	protected char[] toMedianArray(){return copyToArray();}

	/** Copies logical values into an independent writable array.
	 * Concurrent updates can make the result differ from any single point in time.
	 */
	private char[] copyToArray(){
		char[] array2=new char[length()];
		for(int i=0; i<array2.length; i++){
			array2[i]=(char)get(i);
		}
		return array2;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public final AtomicIntegerArray array;
//	@Override
//	public int length(){return maxIndex-minIndex+1;}
	/**
	 * Returns the length of the underlying atomic integer array.
	 * Values are packed two per integer, with one extra packed cell allocated by
	 * the constructor.
	 * @return Length of the underlying integer array
	 */
	@Override
	public int arrayLength(){return array.length();}

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private static boolean OVERFLOWED=false;

	private static final int lowerMask=0x0000FFFF;
	private static final int upperMask=0xFFFF0000;

}
