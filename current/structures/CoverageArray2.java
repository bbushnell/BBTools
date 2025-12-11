package structures;
import dna.Data;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Tools;


public class CoverageArray2 extends CoverageArray {
	
	private static final long serialVersionUID = 8242586595591123194L;
	
	public static void main(String[] args){
		runSpeedTest(args);
		
//		translateGenomeBuild(args);
	}
	
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
	
	public CoverageArray2(int chrom, int len){
		super(chrom, len);
		array=KillSwitch.allocChar1D(len);
	}
	
	/**
	 * Increments coverage at the specified location by the given amount.
	 * Delegates to set() method with current value plus increment amount.
	 * @param loc Genomic position to increment
	 * @param amt Amount to add to current coverage value
	 */
	@Override
	public void increment(int loc, int amt) {
		set(loc, get(loc)+amt);
	}
	
	@Override
	public void increment(int loc) {
		set(loc, get(loc)+1);
	}

	@Override
	public synchronized void incrementRangeSynchronized(int min, int max, int amt) {
		incrementRange(min, max, amt);
	}

	@Override
	public void incrementRange(int min, int max, int amt) {
		if(min<0){min=0;}
		if(max>=array.length){//Increase size
			int newlen=1+(7*max(array.length, max))/4;
			assert(newlen>max);
			resize(newlen);
			assert(array.length==newlen);
		}else if(max<0){max=-1;}
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
	
	
	/**
	 * Sets coverage value at the specified genomic position.
	 * Automatically resizes array if position exceeds current length.
	 * Caps values at Character.MAX_VALUE with overflow warning on first occurrence.
	 * Updates minIndex and maxIndex bounds tracking.
	 *
	 * @param loc Genomic position to set (negative positions are ignored)
	 * @param val Coverage value to set (capped at Character.MAX_VALUE)
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
	
	/**
	 * Gets coverage value at the specified genomic position.
	 * Returns 0 for out-of-bounds positions.
	 * @param loc Genomic position to query
	 * @return Coverage value at position, or 0 if out of bounds
	 */
	@Override
	public int get(int loc){
		return loc>=array.length || loc<0 ? 0 : array[loc];
	}
	
	/**
	 * Resizes the internal char array to accommodate larger genomic coordinates.
	 * Preserves existing coverage data up to the minimum of old and new lengths.
	 * @param newlen New array length (must be greater than maxIndex)
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
	
	/**
	 * Returns string representation of coverage values from 0 to maxIndex.
	 * Formats as comma-separated list enclosed in brackets.
	 * @return String representation of coverage array contents
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
	
	/** Returns direct reference to the internal char array.
	 * @return Internal char array containing coverage values */
	@Override
	public char[] toArray() {return array;}
	
	public char[] array;
//	@Override
//	public int length(){return maxIndex-minIndex+1;}
	/**
	 * Returns the current allocated length of the internal array.
	 * May be larger than the actual data range.
	 * @return Length of internal char array
	 */
	@Override
	public int arrayLength(){return array.length;}
	
	private static boolean OVERFLOWED=false;
	/**
	 * 
	 */
//	private static final long serialVersionUID = -7493066925636540386L;
	
}
