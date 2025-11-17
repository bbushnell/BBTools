package structures;
import java.util.concurrent.atomic.AtomicIntegerArray;

import shared.KillSwitch;
import shared.Tools;

/**
 * Atomic version 
 * @author Brian Bushnell
 * @date Sep 20, 2014
 *
 */
public class CoverageArray2A extends CoverageArray {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 98483952072098494L;
	
	/** Program entry point.
	 * @param args Command-line arguments */
	public static void main(String[] args){
		//TODO
	}
	
	/**
	 * Constructs an atomic coverage array for the specified chromosome and length.
	 * Allocates an AtomicIntegerArray to store packed 16-bit coverage values.
	 * @param chrom Chromosome identifier
	 * @param len Length of the coverage array in bases
	 */
	public CoverageArray2A(int chrom, int len){
		super(chrom, len);
		int intLen=intIdx(len)+1;
		array=KillSwitch.allocAtomicInt(intLen);
		minIndex=0;
		maxIndex=len-1;
	}
	
	/**
	 * Converts a base position to the corresponding integer array index.
	 * Two 16-bit values are packed per integer, so array index is (idx+1)/2.
	 * @param idx Base position in the coverage array
	 * @return Integer array index containing the packed coverage value
	 */
	private static final int intIdx(int idx) {
		return (idx+1)/2;
	}
	
	/**
	 * @param loc
	 */
	@Override
	public void increment(int loc){
		increment(loc, 1);
	}
	
	@Override
	public void increment(final int loc, final int amt) {
		assert(amt>=0) : "This does not currently allow negative increments.";
		final int intIdx=intIdx(loc);
		boolean overflow=((loc&1)==1 ? incrementUpper(intIdx, amt) : incrementLower(intIdx, amt));
		if(overflow && !OVERFLOWED){
			 System.err.println("Note: Coverage capped at "+0xFFFF);
			 OVERFLOWED=true;
		}
	}
	
	private boolean incrementLower(final int intIdx, final int amt) {
		boolean overflow=false;
		for(int oldVal=0, actual=amt; oldVal!=actual; ) {
			oldVal=array.get(intIdx);
			int lower=oldVal&lowerMask, upper=oldVal&upperMask;
			int charVal=lower;
			int charVal2=Tools.min(0xFFFF, charVal+amt);
			overflow=(charVal2>=0xFFFF);
			int newVal=(charVal2)|upper;
			actual=array.compareAndExchange(intIdx, oldVal, newVal);
		}
		return overflow;
	}
	
	/**
	 * Atomically increments the upper 16-bit value at the specified integer index.
	 * Uses compare-and-exchange loop to ensure atomic updates.
	 * Preserves the lower 16-bit value while updating the upper value.
	 *
	 * @param intIdx Integer array index
	 * @param amt Amount to increment by
	 * @return true if overflow occurred (value reached 0xFFFF), false otherwise
	 */
	private boolean incrementUpper(int intIdx, int amt) {
		boolean overflow=false;
		for(int oldVal=0, actual=amt; oldVal!=actual; ) {
			oldVal=array.get(intIdx);
			int lower=oldVal&lowerMask, upper=oldVal&upperMask;
			int charVal=upper>>>16;
			int charVal2=Tools.min(0xFFFF, charVal+amt);
			overflow=(charVal2>=0xFFFF);
			int newVal=(charVal2<<16)|lower;
			actual=array.compareAndExchange(intIdx, oldVal, newVal);
		}
		return overflow;
	}

	@Override
	public void incrementRangeSynchronized(int min, int max, int amt) {
		incrementRange(min, max, amt);//Synchronized is not needed
	}
	
	/**
	 * Increments coverage values over a range using individual increment operations.
	 * Safe but slower method that calls increment() for each position in the range.
	 *
	 * @param min Starting position (inclusive)
	 * @param max Ending position (inclusive)
	 * @param amt Amount to increment by
	 */
	public void incrementRangeSlow(int min, int max, int amt){
		if(min<0){min=0;}
		if(max>maxIndex){max=maxIndex;}
		for(int loc=min; loc<=max; loc++){
			increment(loc, amt);
		}
	}
	
	@Override
	public void incrementRange(int min, int max, int amt){
		if(amt>0xFFF || true) {
			incrementRangeSlow(min, max, amt);
			return;
		}
		//TODO:  This should be 2x as fast, but currently gives slightly wrong results.
		//Off by ~1% so probably a boundary issue
		//Try printing range and aborting
		if(max>maxIndex){max=maxIndex;}
		if((min&1)==1) {increment(min, amt);}
		if((max&1)==0) {increment(max, amt);}
		int minIdx=intIdx(min+1);
		int maxIdx=intIdx(max-1);
		for(int i=minIdx; i<=maxIdx; i++) {
			for(int oldVal=0, actual=amt; oldVal!=actual; ) {
				oldVal=array.get(i);
				int lower=oldVal&lowerMask, upper=(oldVal&upperMask)>>16;
				int lower2=Tools.min(0xFFFF, lower+amt);
				int upper2=Tools.min(0xFFFF, upper+amt);
				int newVal=lower2|(upper2<<16);
				actual=array.compareAndExchange(i, oldVal, newVal);
			}
		}
	}
	
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
	
	private boolean setLower(final int intIdx, final int amt) {
		for(int oldVal=0, actual=amt; oldVal!=actual; ) {
			oldVal=array.get(intIdx);
			int lower=amt, upper=oldVal&upperMask;
			int newVal=lower|upper;
			actual=array.compareAndExchange(intIdx, oldVal, newVal);
		}
		return false;
	}
	
	/**
	 * Atomically sets the upper 16-bit value at the specified integer index.
	 * Uses compare-and-exchange loop to ensure atomic updates.
	 * Preserves the lower 16-bit value while setting the upper value.
	 *
	 * @param intIdx Integer array index
	 * @param amt Value to set in the upper 16 bits
	 * @return Always returns false (no overflow possible with set operations)
	 */
	private boolean setUpper(int intIdx, int amt) {
		for(int oldVal=0, actual=amt; oldVal!=actual; ) {
			oldVal=array.get(intIdx);
			int lower=oldVal&lowerMask, upper=amt<<16;
			int newVal=lower|upper;
			actual=array.compareAndExchange(intIdx, oldVal, newVal);
		}
		return false;
	}
	
	@Override
	public int get(int loc){
		final int intIdx=intIdx(loc);
		final int intVal=intIdx<0 || intIdx>=array.length() ? 0 : array.get(intIdx);
		return (loc&1)==1 ? (intVal>>>16) : (intVal&lowerMask);
	}
	
	@Override
	public void resize(int newlen){
		throw new RuntimeException("Resize: Unsupported.");
	}
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		for(int i=0; i<=length(); i++){
			if(i>0){sb.append(", ");}
			sb.append(get(i));
		}
		sb.append(']');
		return sb.toString();
	}
	
	@Override
	public char[] toArray() {
		char[] array2=new char[length()];
		for(int i=0; i<array2.length; i++) {
			array2[i]=(char)get(i);
		}
		return array2;
	}
	
	/** Atomic integer array storing packed 16-bit coverage values */
	public final AtomicIntegerArray array;
//	@Override
//	public int length(){return maxIndex-minIndex+1;}
	@Override
	public int arrayLength(){return array.length();}
	
	/**
	 * Flag indicating whether any coverage value has overflowed the 16-bit limit
	 */
	private static boolean OVERFLOWED=false;
	
	/** Bit mask (0x0000FFFF) for extracting lower 16-bit coverage values */
	private static final int lowerMask=0x0000FFFF;
	/** Bit mask (0xFFFF0000) for extracting upper 16-bit coverage values */
	private static final int upperMask=0xFFFF0000;
	
}
