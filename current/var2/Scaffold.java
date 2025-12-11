package var2;

import shared.Parse;
import shared.Shared;
import shared.Tools;
import stream.SamLine;
import structures.CoverageArray;
import structures.CoverageArray2;
import structures.CoverageArray3;
import structures.CoverageArray3A;

/**
 * Represents a single scaffold (chromosome/contig) in a reference genome.
 * Handles coverage tracking, sequence storage, and variant-related calculations.
 * Supports optional strand-specific coverage tracking and lazy initialization
 * of coverage arrays for memory efficiency.
 *
 * @author Brian Bushnell
 * @contributor Isla Winglet
 */
public class Scaffold {
	
	/**
	 * Constructs a Scaffold by parsing a SAM header line.
	 * Expects SAM format: @SQ	SN:scaffold_0	LN:1785514	AS:build 9
	 * @param line SAM header line as byte array
	 * @param scafnum Scaffold number to assign
	 */
	public Scaffold(byte[] line, int scafnum){
		assert(Tools.startsWith(line, "@SQ\t")) : new String(line);
		number=scafnum;
		int a=0, b=0;
		
		// Skip @SQ field
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		b++;
		a=b;
		
		// Parse SN: field (scaffold name)
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		assert(Tools.startsWith(line, "SN:", a));
		name=new String(line, a+3, b-a-3);
		b++;
		a=b;
		
		// Parse LN: field (length)
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		assert(Tools.startsWith(line, "LN:", a));
		length=Parse.parseInt(line, a+3, b);
		b++;
		a=b;
	}
	
	public Scaffold(String name_, int scafnum_, int len_){
		name=name_;
		number=scafnum_;
		length=len_;
	}
	
	/**
	 * Adds coverage information from a SAM alignment line.
	 * Extracts alignment coordinates and updates coverage arrays.
	 * @param sl SAM line containing alignment information
	 */
	public void add(SamLine sl){
		int start=sl.pos-1;
		int stop=sl.stop(start, false, false);
		increment(start, stop, sl.strand());
	}
	
	/**
	 * Increments coverage for a specified range, with optional strand tracking.
	 * Uses lazy initialization to create coverage arrays only when needed.
	 * Thread-safe through synchronized initialization block.
	 *
	 * @param from Start position (inclusive)
	 * @param to End position (exclusive)
	 * @param strand Strand information (+ or -)
	 */
	public void increment(int from, int to, int strand){
		// Lazy initialization with double-checked locking pattern
		if(!initialized()){
			synchronized(this){
				if(!initialized()){
					// Choose coverage array implementation based on static settings
					ca=useCA3A ? new CoverageArray3A(number, length) : useCA3 ? new CoverageArray3(number, length) : new CoverageArray2(number, length);
					if(trackStrand){
						caMinus=useCA3A ? new CoverageArray3A(number, length) : useCA3 ? new CoverageArray3(number, length) : new CoverageArray2(number, length);
					}
				}
				initialized=true;
			}
		}
		assert(initialized());
		assert(ca!=null);
		
		// Update coverage arrays
		ca.incrementRangeSynchronized(from, to, 1);
		if(trackStrand && strand==Shared.MINUS){caMinus.incrementRangeSynchronized(from, to, 1);}
	}
	
	/**
	 * Legacy synchronized version of increment method.
	 * Less efficient than current implementation but provided for compatibility.
	 *
	 * @param from Start position (inclusive)
	 * @param to End position (exclusive)
	 * @param strand Strand information (+ or -)
	 */
	public synchronized void incrementOld(int from, int to, int strand){
		if(ca==null){
			ca=useCA3 ? new CoverageArray3(number, length) : new CoverageArray2(number, length);
		}
		ca.incrementRange(from, to);
		if(trackStrand && strand==Shared.MINUS){
			if(caMinus==null){
				caMinus=useCA3 ? new CoverageArray3(number, length) : new CoverageArray2(number, length);
			}
			caMinus.incrementRange(from, to);
		}
	}
	
	public String getSequence(SamLine sl) {
		int start=sl.start(false, false);
		int stop=sl.stop(start, false, false);
		return getSequence(start, stop);
	}
	
	public String getSequence(int start, int stop) {
		assert(bases!=null) : this;
		start=Tools.max(0, start);
		stop=Tools.min(bases.length-1, stop);
		String s=new String(bases, start, stop-start+1);
		return s;
	}
	
	public int calcCoverage(Var v){
		return calcCoverage(v, ca);
	}
	
	/**
	 * Calculates minus-strand coverage at a variant position.
	 * Only available when strand tracking is enabled.
	 * @param v Var object defining the position
	 * @return Average minus-strand coverage across the variant region
	 */
	public int minusCoverage(Var v){
		assert(trackStrand);
		return calcCoverage(v, caMinus);
	}
	
	/**
	 * Calculates coverage for a variant using the specified coverage array.
	 * Handles different variant types with appropriate coverage calculation strategies.
	 *
	 * @param v Var object defining the position and type
	 * @param ca Coverage array to query
	 * @return Average coverage appropriate for the variant type
	 */
	public int calcCoverage(Var v, CoverageArray ca){
		final int a=v.start;
		final int b=v.stop;
		if(ca==null || ca.maxIndex<a){return 0;}
		final int type=v.type();
		final int avg;
		final int rlen=v.reflen();
		long sum=0;
		
		if(type==Var.SUB || type==Var.NOCALL || type==Var.DEL){
			// For substitutions and deletions, average coverage across the reference span
			for(int i=a; i<b; i++){
				sum+=ca.get(i);
			}
			avg=(int)Math.round(sum/(double)rlen);
		}else if(type==Var.INS){
			// For insertions, interpolate between flanking positions
			assert(rlen==0 && a==b);
			if(b>=ca.maxIndex){
				sum=2*ca.get(ca.maxIndex);
				avg=(int)(sum/2);
			}else{
				sum=ca.get(a)+ca.get(b);
				avg=(int)Math.ceil(sum/2);
			}
		}else if(type==Var.LJUNCT){
			// Left junction: take coverage from right side
			avg=ca.get(Tools.min(ca.maxIndex, a+1));
		}else if(type==Var.RJUNCT){
			// Right junction: take coverage from left side
			avg=ca.get(Tools.max(0, a-1));
		}else{
			throw new RuntimeException("Unknown type "+type+"\n"+v);
		}
		return avg;
	}
	
	@Override
	public String toString(){
		return "@SQ\tSN:"+name+"\tLN:"+length+"\tID:"+number;
	}
	
	/** Clears coverage arrays to free memory.
	 * Thread-safe operation. */
	public synchronized void clearCoverage(){
		ca=null;
		caMinus=null;
		initialized=false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	public final String name;
	public final int number;
	public final int length;
	private CoverageArray ca;
	private CoverageArray caMinus;
	public byte[] bases;
	private boolean initialized(){return initialized;};
	private boolean initialized;

	/*--------------------------------------------------------------*/
	/*----------------      Static Methods          ----------------*/
	/*--------------------------------------------------------------*/
	
	public static void setCA3(boolean b){useCA3=b;}
	
	public static void setCA3A(boolean b){useCA3A=b;}
	
	public static void setTrackStrand(boolean b){trackStrand=b;}
	
	public static boolean trackStrand(){return trackStrand;}

	/*--------------------------------------------------------------*/
	/*----------------      Static Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	private static boolean useCA3=false;
	private static boolean useCA3A=true;
	private static boolean trackStrand=false;
}