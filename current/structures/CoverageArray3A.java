package structures;
import java.util.concurrent.atomic.AtomicIntegerArray;

import shared.KillSwitch;

/**
 * Thread-safe atomic version of CoverageArray using AtomicIntegerArray.
 * Provides concurrent access to coverage data with overflow protection.
 * Values are capped at Integer.MAX_VALUE to prevent overflow errors.
 *
 * @author Brian Bushnell
 * @date Sep 20, 2014
 */
public class CoverageArray3A extends CoverageArray {
	
	/** Serial version UID for serialization compatibility */
	private static final long serialVersionUID = 98483952072098494L;
	
	public static void main(String[] args){
		//TODO
	}
	
	public CoverageArray3A(int chrom, int len){
		super(chrom, len);
		array=KillSwitch.allocAtomicInt(len);
		minIndex=0;
		maxIndex=len-1;
	}
	
	@Override
	public void increment(int loc){
		increment(loc, 1);
	}
	
	@Override
	public void increment(int loc, int amt) {
		int val=array.addAndGet(loc, (int)amt);
//		assert(val>=0 || amt<0) : "Overflow!";
		if(val<0 && amt>0){
			if(!OVERFLOWED){
				 System.err.println("Note: Coverage capped at "+Integer.MAX_VALUE);
				 OVERFLOWED=true;
			}
			array.set(loc, Integer.MAX_VALUE);
		}
	}

	/**
	 * Increments coverage across a range of positions.
	 * Synchronization is not needed due to atomic operations, so delegates to incrementRange.
	 *
	 * @param min Start position (inclusive)
	 * @param max End position (inclusive)
	 * @param amt Amount to add to each position
	 */
	@Override
	public void incrementRangeSynchronized(int min, int max, int amt) {
		incrementRange(min, max, amt);//Synchronized is not needed
	}
	
	/**
	 * Atomically increments coverage across a range of positions.
	 * Bounds-checks the range and prevents overflow at each position.
	 *
	 * @param min Start position (inclusive, clamped to 0)
	 * @param max End position (inclusive, clamped to maxIndex)
	 * @param amt Amount to add to each position
	 */
	@Override
	public void incrementRange(int min, int max, int amt){
		if(min<0){min=0;}
		if(max>maxIndex){max=maxIndex;}
		boolean over=false;
		for(int loc=min; loc<=max; loc++){
			int val=array.addAndGet(loc, (int)amt);
			if(val<0 && amt>0){
				over=true;
				array.set(loc, Integer.MAX_VALUE);
			}
		}
		if(over && !OVERFLOWED){
			synchronized(CoverageArray3A.class){
				if(!OVERFLOWED){
					System.err.println("Note: Coverage capped at "+Integer.MAX_VALUE);
					OVERFLOWED=true;
				}
			}
		}
	}
	
	/**
	 * Atomically sets the coverage value at the specified location.
	 * @param loc Array position to set
	 * @param val New coverage value
	 */
	@Override
	public void set(int loc, int val){
//		if(loc<0 || loc>=maxIndex){return;}
		array.set(loc, val);
	}
	
	/**
	 * Gets the coverage value at the specified location.
	 * Returns 0 for out-of-bounds positions.
	 * @param loc Array position to query
	 * @return Coverage value at the position, or 0 if out of bounds
	 */
	@Override
	public int get(int loc){
		return loc<0 || loc>=array.length() ? 0 : array.get(loc);
	}
	
	/**
	 * Resize operation is not supported for atomic arrays.
	 * @param newlen New length (ignored)
	 * @throws RuntimeException Always thrown as resize is unsupported
	 */
	@Override
	public void resize(int newlen){
		throw new RuntimeException("Resize: Unsupported.");
	}
	
	/**
	 * Creates a string representation of the coverage array.
	 * Returns values in array notation: [val1, val2, val3, ...]
	 * @return String representation of all coverage values
	 */
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		for(int i=0; i<=length(); i++){
			if(i>0){sb.append(", ");}
			sb.append((int)array.get(i));
		}
		sb.append(']');
		return sb.toString();
	}
	
	/**
	 * Converts the atomic coverage array to a regular int array.
	 * Creates a snapshot of current values.
	 * @return New int array containing current coverage values
	 */
	@Override
	public int[] toArray() {
		int[] array2=new int[length()];
		for(int i=0; i<array2.length; i++) {
			array2[i]=get(i);
		}
		return array2;
	}
	
	
	public final AtomicIntegerArray array;
//	@Override
//	public int length(){return maxIndex-minIndex+1;}
	/** Returns the length of the underlying atomic array.
	 * @return Length of the coverage array */
	@Override
	public int arrayLength(){return array.length();}
	
	private static boolean OVERFLOWED=false;
	
}
