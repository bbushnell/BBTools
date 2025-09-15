package structures;
import java.util.concurrent.atomic.AtomicIntegerArray;

import shared.KillSwitch;

/**
 * Atomic version 
 * @author Brian Bushnell
 * @date Sep 20, 2014
 *
 */
public class CoverageArray3A extends CoverageArray {
	
	/**
	 * 
	 */
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
	
	/**
	 * @param loc
	 */
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

	@Override
	public void incrementRangeSynchronized(int min, int max, int amt) {
		incrementRange(min, max, amt);//Synchronized is not needed
	}
	
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
	
	@Override
	public void set(int loc, int val){
//		if(loc<0 || loc>=maxIndex){return;}
		array.set(loc, val);
	}
	
	@Override
	public int get(int loc){
		return loc<0 || loc>=array.length() ? 0 : array.get(loc);
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
			sb.append((int)array.get(i));
		}
		sb.append(']');
		return sb.toString();
	}
	
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
	@Override
	public int arrayLength(){return array.length();}
	
	private static boolean OVERFLOWED=false;
	
}
