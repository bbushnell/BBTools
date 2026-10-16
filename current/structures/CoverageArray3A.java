package structures;
import java.util.concurrent.atomic.AtomicIntegerArray;

import shared.KillSwitch;

/**
 * Fixed-length 32-bit CoverageArray backed by an AtomicIntegerArray.
 * Individual AtomicIntegerArray accesses are atomic; range increments and the
 * add-then-clamp positive-overflow path are compound and can be observed
 * partially done.
 *
 * @author Brian Bushnell
 * @date Sep 20, 2014
 */
public class CoverageArray3A extends CoverageArray{

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

	/**
	 * Creates a fixed-length atomic 32-bit coverage array.
	 * @param chrom Chromosome identifier.
	 * @param len Logical and backing length.
	 */
	public CoverageArray3A(int chrom, int len){
		super(chrom, len);
		array=KillSwitch.allocAtomicInt(len);
		minIndex=0;
		maxIndex=len-1;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Increments one location by 1.
	 * @param loc Location to increment; must be in bounds.
	 */
	@Override
	public void increment(int loc){
		increment(loc, 1);
	}

	/**
	 * Increments one location with an atomic add and a separate overflow clamp.
	 * @param loc Location to increment; must be in bounds.
	 * @param amt Amount to add; positive overflow saturates.
	 */
	@Override
	public void increment(int loc, int amt){
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
	 * Delegates to incrementRange without synchronizing the full range.
	 * @param min First position to increment.
	 * @param max Last position to increment.
	 * @param amt Amount to add at each position.
	 */
	@Override
	public void incrementRangeSynchronized(int min, int max, int amt){
		incrementRange(min, max, amt);//Synchronized is not needed
	}

	/**
	 * Increments an inclusive fixed-length range.
	 * Each array access is atomic, but the whole range and any overflow clamp are
	 * not one atomic transaction.
	 * @param min First position to increment; negative starts are clamped to 0.
	 * @param max Last position to increment; high endpoints are clamped to maxIndex.
	 * @param amt Amount to add at each position; positive overflow saturates.
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
	 * Atomically sets one location.
	 * @param loc Location to set; must be in bounds.
	 * @param val New value stored without saturation.
	 */
	@Override
	public void set(int loc, int val){
//		if(loc<0 || loc>=maxIndex){return;}
		array.set(loc, val);
	}

	/**
	 * Gets one location.
	 * @param loc Location to query.
	 * @return Stored value, or 0 for negative and out-of-bounds positions.
	 */
	@Override
	public int get(int loc){
		return loc<0 || loc>=array.length() ? 0 : array.get(loc);
	}

	/**
	 * Reports that this fixed-length atomic implementation cannot resize.
	 * @param newlen Requested length.
	 * @throws RuntimeException Always thrown.
	 */
	@Override
	public void resize(int newlen){
		throw new RuntimeException("Resize: Unsupported.");
	}

	/**
	 * Formats the full logical length.
	 * @return String representation of all logical coverage values.
	 */
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		for(int i=0; i<length(); i++){
			if(i>0){sb.append(", ");}
			sb.append((int)array.get(i));
		}
		sb.append(']');
		return sb.toString();
	}

	/**
	 * Copies logical values into a new array one element at a time.
	 * Concurrent updates can make the result differ from any single point in time.
	 * @return Per-element copy of logical coverage values.
	 */
	@Override
	public int[] toArray(){return copyToArray();}

	/** Returns independent writable storage for median sorting.
	 * The copy is assembled element by element, not as one array-wide atomic snapshot.
	 */
	@Override
	protected int[] toMedianArray(){return copyToArray();}

	/** Copies logical values into an independent writable array.
	 * Concurrent updates can make the result differ from any single point in time.
	 */
	private int[] copyToArray(){
		int[] array2=new int[length()];
		for(int i=0; i<array2.length; i++){
			array2[i]=get(i);
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
	 * Returns the fixed backing-array length.
	 */
	@Override
	public int arrayLength(){return array.length();}

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private static boolean OVERFLOWED=false;

}
