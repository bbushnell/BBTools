package structures;
import dna.Data;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Tools;

/**
 * Growable 16-bit coverage array backed by a public live char array.
 * Values above {@link Character#MAX_VALUE} are capped when detected, but
 * negative values and integer overflow can still wrap through the char cast.
 * Mutators are unsynchronized except for incrementRangeSynchronized(); callers
 * that share an instance must coordinate other mutation and direct array access.
 * set() and increment methods update touched-coordinate bounds but leave
 * inherited aggregate caches stale unless callers reset them.
 * The inherited logical length can differ from the resized backing array length.
 */
public class CoverageArray2 extends CoverageArray{

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private static final long serialVersionUID=8242586595591123194L;

	/*--------------------------------------------------------------*/
	/*----------------             Main             ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs the legacy serialized-array timing helper. */
	public static void main(String[] args){
		runSpeedTest(args);

//		translateGenomeBuild(args);
	}

	/**
	 * Reads a serialized CoverageArray2 from args[1], sets its chromosome from
	 * args[0], writes a genome-build-named .ca file in the input or current
	 * directory, and rereads it while printing load and reread timings.
	 */
	public static void runSpeedTest(String[] args){

		long time1=System.nanoTime();

		CoverageArray2 ca=(CoverageArray2)read(args[1]);
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
		ca=(CoverageArray2)read(outfile);
		long time4=System.nanoTime();

		System.out.println("minIndex="+ca.minIndex+", maxIndex="+ca.maxIndex+", length="+ca.array.length+
				"; time="+Tools.format("%.3f seconds", (time4-time3)/1000000000d));


	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Allocates 16-bit coverage storage for chrom with inherited logical length len. */
	public CoverageArray2(int chrom, int len){
		super(chrom, len);
		array=KillSwitch.allocChar1D(len);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Adds amt to the stored value at loc via int arithmetic and stores it
	 * through set(int, int).
	 * @param loc Genomic position to increment.
	 * @param amt Amount to add to the current coverage value.
	 */
	@Override
	public void increment(int loc, int amt){
		set(loc, get(loc)+amt);
	}

	/** Adds one to loc through set(int, int). */
	@Override
	public void increment(int loc){
		set(loc, get(loc)+1);
	}

	/** Synchronized wrapper around incrementRange(int, int, int). */
	@Override
	public synchronized void incrementRangeSynchronized(int min, int max, int amt){
		incrementRange(min, max, amt);
	}

	/**
	 * Adds amt to every position in [min,max], resizing when max exceeds the
	 * backing array and capping detected values above Character.MAX_VALUE.
	 * Negative min is clamped to zero; a negative max has no effect.
	 * @param min First genomic position to increment.
	 * @param max Last genomic position to increment, inclusive.
	 * @param amt Amount to add to each current coverage value.
	 */
	@Override
	public void incrementRange(int min, int max, int amt){
		if(max<0){return;}
		if(min<0){min=0;}
		if(max>=array.length){//Increase size
			int newlen=1+(7*max(array.length, max))/4;
			assert(newlen>max);
			resize(newlen);
			assert(array.length==newlen);
		}
		for(int i=min; i<=max; i++){
			int val=array[i]+amt;
			if(val>Character.MAX_VALUE){
				val=Character.MAX_VALUE;
				if(!OVERFLOWED){
					System.err.println("Note: Coverage capped at "+(int)(Character.MAX_VALUE)+"; please use the flag 32bit for higher values.");
					OVERFLOWED=true;
				}
			}
			array[i]=(char)val;
		}
		minIndex=min(min, minIndex);
		maxIndex=max(max, maxIndex);
	}


	/** Sets loc to val, resizing for high positions and ignoring negative positions.
	 * Values above Character.MAX_VALUE are capped; all other ints are cast to char.
	 * @param loc Genomic position to set.
	 * @param val Coverage value to store.
	 */
	@Override
	public void set(int loc, int val){

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

		if(val>Character.MAX_VALUE && !OVERFLOWED){
			System.err.println("Note: Coverage capped at "+(int)(Character.MAX_VALUE)+"; please use the flag 32bit for higher values.");
			OVERFLOWED=true;
		}
		array[loc]=(val>Character.MAX_VALUE ? Character.MAX_VALUE : (char)val);
		minIndex=min(loc, minIndex);
		maxIndex=max(loc, maxIndex);
	}

	/** Gets the coverage value at loc.
	 * @param loc Genomic position to query.
	 * @return Coverage value at loc, or 0 for negative/out-of-bounds loc.
	 */
	@Override
	public int get(int loc){
		return loc>=array.length || loc<0 ? 0 : array[loc];
	}

	/** Replaces the backing array with newlen cells and preserves the copied prefix.
	 * @param newlen New backing array length; asserted to exceed maxIndex.
	 */
	@Override
	public void resize(int newlen){
//		System.err.println("Resized CoverageArray "+chromosome+" to "+newlen);
		char[] temp=KillSwitch.allocChar1D(newlen);
		int lim=min(array.length, newlen);
		assert(lim>maxIndex) : lim+","+maxIndex;
		for(int i=0; i<lim; i++){
			temp[i]=array[i];
		}
		array=temp;
	}

	/** Returns coverage values from 0 to maxIndex as a comma-separated list.
	 * @return String representation of populated coverage positions.
	 */
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		for(int i=0; i<=maxIndex; i++){
			if(i>0){sb.append(", ");}
			sb.append((int)array[i]);
		}
		sb.append(']');
		return sb.toString();
	}

	/** Returns the live backing array; caller modifications alter this CoverageArray2.
	 * @return Live backing array of 16-bit coverage values.
	 */
	@Override
	public char[] toArray(){return array;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Live 16-bit coverage storage exposed by toArray(). */
	public char[] array;
//	@Override
//	public int length(){return maxIndex-minIndex+1;}

	/** Returns the allocated backing length, not the inherited logical length. */
	@Override
	public int arrayLength(){return array.length;}

	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** True after any CoverageArray2 emitted the 16-bit cap warning. */
	private static boolean OVERFLOWED=false;
//	private static final long serialVersionUID = -7493066925636540386L;

}
