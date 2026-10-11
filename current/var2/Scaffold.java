package var2;

import parse.Parse;
import shared.Shared;
import shared.Tools;
import stream.SamLine;
import structures.CoverageArray;
import structures.CoverageArray2;
import structures.CoverageArray3;
import structures.CoverageArray3A;

/**
 * Reference-scaffold metadata, optional sequence, and lazily allocated coverage.
 * Coverage updates and sequence slices use zero-based inclusive endpoints;
 * variant spans use {@link Var}'s zero-based exclusive stop.
 * The public sequence array is supplied by the caller, not copied or allocated here.
 * Configure coverage storage and strand tracking before use, and finish coverage
 * updates before querying or clearing it. Lazy increment allocation publishes its
 * arrays through a volatile flag; reset/read lifecycles still require caller
 * sequencing and are not a general thread-safety guarantee.
 *
 * @author Brian Bushnell
 * @contributor Isla
 * @contributor Sucrose
 */
public class Scaffold{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses the first three tab-delimited fields in the order {@code @SQ},
	 * {@code SN:name}, {@code LN:length}. Additional fields are ignored.
	 * @param line Header bytes with name and length in those fixed positions
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

	/** Creates metadata without allocating coverage or sequence storage.
	 * @param name_ Reference name
	 * @param scafnum_ Reference number used by variants and coverage arrays
	 * @param len_ Reference length in bases */
	public Scaffold(String name_, int scafnum_, int len_){
		name=name_;
		number=scafnum_;
		length=len_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts the entire inclusive span from POS-1 through SamLine.stop,
	 * excluding clip contributions. Does not check mapping status or split the
	 * span into aligned blocks; callers select the records to count.
	 * @param sl Alignment on this scaffold */
	public void add(SamLine sl){
		int start=sl.pos-1;
		int stop=sl.stop(start, false, false);
		increment(start, stop, sl.strand());
	}

	/**
	 * Adds one to total coverage across an inclusive range, and to minus coverage
	 * when strand tracking is enabled and strand is {@link Shared#MINUS}.
	 * Allocates arrays on first use: atomic int storage takes precedence over
	 * ordinary int storage, otherwise char storage is used. Range bounds and
	 * saturation are handled by the selected coverage-array implementation.
	 *
	 * @param from Zero-based start position (inclusive)
	 * @param to Zero-based end position (inclusive)
	 * @param strand Numeric strand code; only Shared.MINUS adds minus coverage
	 */
	public void increment(int from, int to, int strand){
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
	 * Legacy update path, synchronized on this scaffold. Chooses int or char
	 * storage from useCA3, ignoring useCA3A, and allocates minus storage only on
	 * a minus-strand update. Does not set the newer initialization flag.
	 * Do not mix this path with increment between coverage resets: increment
	 * can replace arrays created here while initialized is still false.
	 *
	 * @param from Zero-based start position (inclusive)
	 * @param to Zero-based end position (inclusive)
	 * @param strand Numeric strand code; only Shared.MINUS adds minus coverage
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

	/** Returns a reference slice for the alignment span, excluding clip contributions.
	 * Delegates coordinate handling to SamLine and getSequence(int, int).
	 * @param sl Alignment on this scaffold
	 * @return New string decoded from the reference bytes */
	public String getSequence(SamLine sl){
		int start=sl.start(false, false);
		int stop=sl.stop(start, false, false);
		return getSequence(start, stop);
	}

	/** Returns a new string for an inclusive, zero-based reference interval.
	 * Clamps start upward to zero and stop downward to bases.length-1; the resulting
	 * offset and length must be valid for String's byte-array constructor.
	 * Requires bases to be present; decoding uses the platform default charset.
	 * @param start First base, inclusive
	 * @param stop Last base, inclusive
	 * @return Reference sequence in stored orientation */
	public String getSequence(int start, int stop){
		assert(bases!=null) : this;
		start=Tools.max(0, start);
		stop=Tools.min(bases.length-1, stop);
		String s=new String(bases, start, stop-start+1);
		return s;
	}

	/** Calculates variant coverage from the total-coverage array.
	 * @param v Variant on this scaffold
	 * @return Type-specific coverage, or zero if no applicable storage exists */
	public int calcCoverage(Var v){
		return calcCoverage(v, ca);
	}

	/**
	 * Calculates variant coverage from the minus-strand array.
	 * Requires strand tracking to be enabled.
	 * @param v Variant on this scaffold
	 * @return Type-specific minus coverage, or zero if no applicable storage exists
	 */
	public int minusCoverage(Var v){
		assert(trackStrand);
		return calcCoverage(v, caMinus);
	}

	/**
	 * Calculates coverage from an explicit array without changing either argument.
	 * SUB, NOCALL and DEL use the rounded mean across [start, stop).
	 * INS has equal start and stop and samples that position, not distinct flanks.
	 * LJUNCT samples start+1 capped at maxIndex; RJUNCT samples start-1 floored at zero.
	 * Returns zero for a null array or a start beyond maxIndex. Other out-of-range
	 * positions are delegated to the array's get method.
	 *
	 * @param v Variant on this scaffold, with coordinates appropriate for its type
	 * @param ca Coverage array to query, or null
	 * @return Coverage for the variant type
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
			// Insertions have zero reference span, so a and b denote the same position.
			// Widen before arithmetic so supported int coverage cannot overflow (VAR2-003).
			assert(rlen==0 && a==b);
			if(b>=ca.maxIndex){
				sum=2L*ca.get(ca.maxIndex);
				avg=(int)(sum/2);
			}else{
				sum=(long)ca.get(a)+ca.get(b);
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

	/** Returns an SQ-style header containing name, length, and scaffold number. */
	@Override
	public String toString(){
		return "@SQ\tSN:"+name+"\tLN:"+length+"\tID:"+number;
	}

	/** Releases coverage arrays and resets lazy initialization, retaining metadata
	 * and bases. Call only after coverage readers and writers have finished;
	 * their accesses do not all acquire this method's monitor. */
	public synchronized void clearCoverage(){
		ca=null;
		caMinus=null;
		initialized=false;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/

	/** Reference name. */
	public final String name;
	/** Reference number used to associate variants with this scaffold. */
	public final int number;
	/** Declared reference length in bases; not checked against bases.length here. */
	public final int length;
	/** Total coverage, allocated on demand and discarded by clearCoverage. */
	private CoverageArray ca;
	/** Minus-strand coverage, present only when allocated by a strand-tracking path. */
	private CoverageArray caMinus;
	/** Caller-supplied reference sequence; null until assigned and not copied here. */
	public byte[] bases;
	/** Returns the volatile lazy-increment publication flag. */
	private boolean initialized(){return initialized;}
	/** Publishes increment's lazy coverage arrays; not set by incrementOld. */
	private volatile boolean initialized;

	/*--------------------------------------------------------------*/
	/*----------------      Static Methods          ----------------*/
	/*--------------------------------------------------------------*/

	/** Selects int rather than char storage for future non-atomic allocations. */
	public static void setCA3(boolean b){useCA3=b;}

	/** Selects atomic int storage for future increment allocations, overriding useCA3. */
	public static void setCA3A(boolean b){useCA3A=b;}

	/** Configures minus-strand tracking before allocation; does not retrofit existing arrays. */
	public static void setTrackStrand(boolean b){trackStrand=b;}

	/** Returns the current process-wide strand-tracking setting. */
	public static boolean trackStrand(){return trackStrand;}

	/*--------------------------------------------------------------*/
	/*----------------      Static Fields           ----------------*/
	/*--------------------------------------------------------------*/

	/** Use ordinary int coverage when atomic storage is not selected. */
	private static boolean useCA3=false;
	/** Use atomic int coverage in increment; ignored by incrementOld. */
	private static boolean useCA3A=true;
	/** Allocate/update minus-strand coverage in addition to total coverage. */
	private static boolean trackStrand=false;
}
