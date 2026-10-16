package structures;
import dna.Data;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Tools;


/**
 * CoverageArray backed by a growing signed 32-bit integer array.
 * Positive writes saturate above Integer.MAX_VALUE, high nonnegative locations
 * grow the backing array, and negative locations are ignored by set().
 */
public class CoverageArray3 extends CoverageArray{

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private static final long serialVersionUID=-4216985130070239610L;

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		runSpeedTest(args);

//		translateGenomeBuild(args);
	}

	public static void runSpeedTest(String[] args){

		long time1=System.nanoTime();

		CoverageArray3 ca=(CoverageArray3)read(args[1]);
		ca.chromosome=Byte.parseByte(args[0]);
		long time2=System.nanoTime();

//		int dot=args[1].lastIndexOf(".");
//		String outfile=args[1].substring(0,dot)+".ca";

		args[1]=args[1].replace('\\', '/');
		int slash=args[1].lastIndexOf('/');
		String outfile;
		if(slash<1){
			outfile="coverage-chr"+ca.chromosome+"-build"+Data.GENOME_BUILD+".ca";
		}else{
			outfile=args[1].substring(0,slash+1)+"coverage-chr"+ca.chromosome+"-build"+Data.GENOME_BUILD+".ca";
		}

		System.out.println("minIndex="+ca.minIndex+", maxIndex="+ca.maxIndex+", length="+ca.array.length+
				"; time="+Tools.format("%.3f seconds", (time2-time1)/1000000000d));

		long time3=System.nanoTime();
		ReadWrite.write(ca, outfile, false);
		ca=null;
		System.gc();
		ca=(CoverageArray3)read(outfile);
		long time4=System.nanoTime();

		System.out.println("minIndex="+ca.minIndex+", maxIndex="+ca.maxIndex+", length="+ca.array.length+
				"; time="+Tools.format("%.3f seconds", (time4-time3)/1000000000d));


	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a resizable 32-bit coverage array with the requested logical length.
	 * @param chrom Chromosome identifier.
	 * @param len Initial logical and backing length.
	 */
	public CoverageArray3(int chrom, int len){
		super(chrom, len);
		array=KillSwitch.allocInt1D(len);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Increments one location by a long amount, growing the backing array when needed.
	 * @param loc Location to increment; negative values are ignored by set().
	 * @param amt Amount to add before saturation to Integer.MAX_VALUE.
	 */
	public void increment(int loc, long amt){
		set(loc, get(loc)+amt);
	}

	@Override
	public void increment(int loc){
		set(loc, get(loc)+1L);
	}

	@Override
	public void increment(int loc, int amt){
		increment(loc, (long)amt);
	}

	/**
	 * Synchronized wrapper for inclusive range increments on the growing int array.
	 * @param min First position to increment.
	 * @param max Last position to increment.
	 * @param amt Amount to add at each position.
	 */
	@Override
	public synchronized void incrementRangeSynchronized(int min, int max, int amt){
		incrementRange(min, max, amt);
	}

	/**
	 * Increments an inclusive range by an int amount.
	 * @param min First position to increment; negative starts are clamped to 0.
	 * @param max Last position to increment; ranges ending below 0 are ignored,
	 * and high endpoints grow the array.
	 * @param amt Amount to add at each position.
	 */
	@Override
	public void incrementRange(int min, int max, int amt){
		incrementRange(min, max, (long)amt);
	}

	/**
	 * Increments an inclusive range by a long amount.
	 * @param min First position to increment; negative starts are clamped to 0.
	 * @param max Last position to increment; ranges ending below 0 are ignored,
	 * and high endpoints grow the array.
	 * @param amt Amount to add at each position before positive saturation.
	 */
	public void incrementRange(int min, int max, long amt){
		if(max<0){return;}
		if(min<0){min=0;}
		if(max>=array.length){//Increase size
			int newlen=1+(7*max(array.length, max))/4;
			assert(newlen>max);
			resize(newlen);
			assert(array.length==newlen);
		}
		for(int i=min; i<=max; i++){
			long val=array[i]+amt;
			if(val>Integer.MAX_VALUE){
				val=Integer.MAX_VALUE;
				 if(!OVERFLOWED){
					 System.err.println("Note: Coverage capped at "+Integer.MAX_VALUE);
					 OVERFLOWED=true;
				 }
			}
			array[i]=(int)val;
		}
		minIndex=min(min, minIndex);
		maxIndex=max(max, maxIndex);
	}

	/**
	 * Sets one location to an int value.
	 * @param loc Location to set; negative values are ignored.
	 * @param val Value to store before positive saturation.
	 */
	@Override
	public void set(int loc, int val){
		set(loc, (long)val);
	}

	/**
	 * Sets one location to a long value, growing the backing array when needed.
	 * @param loc Location to set; negative values are ignored.
	 * @param val Value to store before positive saturation to Integer.MAX_VALUE.
	 */
	public void set(int loc, long val){

		if(loc>=array.length){//Increase size
			int newlen=1+(7*max(array.length, loc))/4;
			assert(newlen>loc);
			resize(newlen);
			assert(array.length==newlen);
		}else if(loc<0){
//			minIndex=min(0, minIndex);
//			maxIndex=max(0, maxIndex);
			return;
		}

		if(val>Integer.MAX_VALUE && !OVERFLOWED){
			System.err.println("Note: Coverage capped at "+Integer.MAX_VALUE);
			OVERFLOWED=true;
		}
		array[loc]=(val>Integer.MAX_VALUE ? Integer.MAX_VALUE : (int)val);
		minIndex=min(loc, minIndex);
		maxIndex=max(loc, maxIndex);
	}

	/**
	 * Gets one location.
	 * @param loc Location to retrieve.
	 * @return Stored value, or 0 for negative and beyond-backing positions.
	 */
	@Override
	public int get(int loc){
		return loc>=array.length || loc<0 ? 0 : array[loc];
	}

	/**
	 * Replaces the backing array while preserving copied positions.
	 * @param newlen New backing-array capacity.
	 */
	@Override
	public void resize(int newlen){
//		System.err.println("Resized CoverageArray "+chromosome+" to "+newlen);
		int[] temp=KillSwitch.allocInt1D(newlen);
		int lim=min(array.length, newlen);
		assert(lim>maxIndex) : lim+","+maxIndex;
		for(int i=0; i<lim; i++){
			temp[i]=array[i];
		}
		array=temp;
	}

	/**
	 * Formats the full logical length, including untouched trailing zeroes.
	 */
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		for(int i=0; i<length(); i++){
			if(i>0){sb.append(", ");}
			sb.append((int)array[i]);
		}
		sb.append(']');
		return sb.toString();
	}

	/**
	 * Returns the mutable internal storage directly.
	 * @return Live backing array.
	 */
	@Override
	public int[] toArray(){return array;}


	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public int[] array;
//	@Override
//	public int length(){return maxIndex-minIndex+1;}
	/**
	 * Returns the current backing-array capacity, which can exceed length().
	 */
	@Override
	public int arrayLength(){return array.length;}

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private static boolean OVERFLOWED=false;

}
