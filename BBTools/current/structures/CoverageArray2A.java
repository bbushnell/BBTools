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
	
	public static void main(String[] args){
		//TODO
	}
	
	public CoverageArray2A(int chrom, int len){
		super(chrom, len);
		int intLen=intIdx(len)+1;
		array=KillSwitch.allocAtomicInt(intLen);
		minIndex=0;
		maxIndex=len-1;
	}
	
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
	
	public final AtomicIntegerArray array;
//	@Override
//	public int length(){return maxIndex-minIndex+1;}
	@Override
	public int arrayLength(){return array.length();}
	
	private static boolean OVERFLOWED=false;
	
	private static final int lowerMask=0x0000FFFF;
	private static final int upperMask=0xFFFF0000;
	
}
