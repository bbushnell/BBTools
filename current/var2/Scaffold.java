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
	 * 
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
	
	/**
	 * Constructs a Scaffold with explicit parameters.
	 * 
	 * @param name_ Scaffold name
	 * @param scafnum_ Scaffold number
	 * @param len_ Scaffold length in bases
	 */
	public Scaffold(String name_, int scafnum_, int len_){
		name=name_;
		number=scafnum_;
		length=len_;
	}
	
	/**
	 * Adds coverage information from a SAM alignment line.
	 * Extracts alignment coordinates and updates coverage arrays.
	 * 
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
	
	/**
	 * Extracts reference sequence covered by a SAM alignment.
	 * 
	 * @param sl SAM line defining the region
	 * @return Reference sequence as String
	 */
	public String getSequence(SamLine sl) {
		int start=sl.start(false, false);
		int stop=sl.stop(start, false, false);
		return getSequence(start, stop);
	}
	
	/**
	 * Extracts reference sequence for a specified coordinate range.
	 * 
	 * @param start Start position (inclusive)
	 * @param stop Stop position (inclusive)
	 * @return Reference sequence as String
	 */
	public String getSequence(int start, int stop) {
		assert(bases!=null) : this;
		start=Tools.max(0, start);
		stop=Tools.min(bases.length-1, stop);
		String s=new String(bases, start, stop-start+1);
		return s;
	}
	
	/**
	 * Calculates total coverage at a variant position.
	 * 
	 * @param v Var object defining the position
	 * @return Average coverage across the variant region
	 */
	public int calcCoverage(Var v){
		return calcCoverage(v, ca);
	}
	
	/**
	 * Calculates minus-strand coverage at a variant position.
	 * Only available when strand tracking is enabled.
	 * 
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
	
	/**
	 * Returns SAM header representation of this scaffold.
	 * 
	 * @return SAM @SQ header line as String
	 */
	@Override
	public String toString(){
		return "@SQ\tSN:"+name+"\tLN:"+length+"\tID:"+number;
	}
	
	/**
	 * Clears coverage arrays to free memory.
	 * Thread-safe operation.
	 */
	public synchronized void clearCoverage(){
		ca=null;
		caMinus=null;
		initialized=false;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Scaffold name (chromosome/contig identifier) */
	public final String name;
	/** Numeric scaffold identifier */
	public final int number;
	/** Length of scaffold in base pairs */
	public final int length;
	/** Primary coverage array for tracking read depth */
	private CoverageArray ca;
	/** Minus-strand coverage array (only used when strand tracking enabled) */
	private CoverageArray caMinus;
	/** Reference sequence bases (may be null if not loaded) */
	public byte[] bases;
	/** Returns whether coverage arrays have been initialized */
	private boolean initialized(){return initialized;};
	/** Initialization status flag */
	private boolean initialized;

	/*--------------------------------------------------------------*/
	/*----------------      Static Methods          ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Sets whether to use CoverageArray3 implementation.
	 * @param b true to use CoverageArray3, false for CoverageArray2 */
	public static void setCA3(boolean b){useCA3=b;}
	
	/** Sets whether to use CoverageArray3A implementation.
	 * @param b true to use CoverageArray3A */
	public static void setCA3A(boolean b){useCA3A=b;}
	
	/** Enables or disables strand-specific coverage tracking.
	 * @param b true to track strand-specific coverage */
	public static void setTrackStrand(boolean b){trackStrand=b;}
	
	/** Returns whether strand-specific coverage tracking is enabled.
	 * @return true if tracking strand-specific coverage */
	public static boolean trackStrand(){return trackStrand;}

	/*--------------------------------------------------------------*/
	/*----------------      Static Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Whether to use CoverageArray3 implementation */
	private static boolean useCA3=false;
	/** Whether to use CoverageArray3A implementation */
	private static boolean useCA3A=true;
	/** Whether to track strand-specific coverage */
	private static boolean trackStrand=false;
}