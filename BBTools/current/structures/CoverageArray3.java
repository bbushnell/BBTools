package structures;
import dna.Data;
import fileIO.ReadWrite;
import shared.KillSwitch;
import shared.Tools;


public class CoverageArray3 extends CoverageArray {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -4216985130070239610L;
	
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
	
	public CoverageArray3(int chrom, int len){
		super(chrom, len);
		array=KillSwitch.allocInt1D(len);
	}
	
	/**
	 * @param loc
	 * @param amt
	 */
	public void increment(int loc, long amt) {
		set(loc, get(loc)+amt);
	}
	
	/**
	 * @param loc
	 */
	@Override
	public void increment(int loc) {
		set(loc, get(loc)+1L);
	}
	
	@Override
	public void increment(int loc, int amt) {
		increment(loc, (long)amt);
	}

	@Override
	public synchronized void incrementRangeSynchronized(int min, int max, int amt) {
		incrementRange(min, max, amt);
	}
	
	@Override
	public void incrementRange(int min, int max, int amt) {
		incrementRange(min, max, (long)amt);
	}

	public void incrementRange(int min, int max, long amt) {
		if(min<0){min=0;}
		if(max>=array.length){//Increase size
			int newlen=1+(7*max(array.length, max))/4;
			assert(newlen>max);
			resize(newlen);
			assert(array.length==newlen);
		}else if(max<0){max=-1;}
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
	
	@Override
	public void set(int loc, int val){
		set(loc, (long)val);
	}
	
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
	
	@Override
	public int get(int loc){
		return loc>=array.length || loc<0 ? 0 : array[loc];
	}
	
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
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		sb.append('[');
		for(int i=0; i<=length(); i++){
			if(i>0){sb.append(", ");}
			sb.append((int)array[i]);
		}
		sb.append(']');
		return sb.toString();
	}
	
	@Override
	public int[] toArray() {return array;}
	
	
	public int[] array;
//	@Override
//	public int length(){return maxIndex-minIndex+1;}
	@Override
	public int arrayLength(){return array.length;}
	
	private static boolean OVERFLOWED=false;
	
}
