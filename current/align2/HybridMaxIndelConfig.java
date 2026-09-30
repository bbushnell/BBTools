package align2;

/** Immutable per-invocation search bounds for selective max-indel retry.
 * Values configure a worker's BBIndex between attempts; they do not mutate
 * process-wide defaults or impose a final CIGAR indel-length filter.
 * Each bound is a base-coordinate distance. Primary and secondary correspond
 * to maxindel and maxindel2 respectively, with larger values for the retry.
 * @author Collei */
public final class HybridMaxIndelConfig {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates ordered narrow/wide bounds without an additional retry MAPQ cutoff. */
	public HybridMaxIndelConfig(int lowPrimary_, int lowSum_, int highPrimary_, int highSum_){
		this(lowPrimary_, lowSum_, highPrimary_, highSum_, 0);
	}

	/** Validates bounds before workers can use this configuration.
	 * @param lowPrimary_ Positive primary bound for the initial search
	 * @param lowSum_ Initial secondary bound, at least lowPrimary_
	 * @param highPrimary_ Retry primary bound, strictly greater than lowPrimary_
	 * @param highSum_ Retry secondary bound, at least highPrimary_ and greater than lowSum_
	 * @param retryMinMapq_ Additional retry acceptance cutoff in [0,255]; zero disables it
	 * @throws IllegalArgumentException If bounds are nonpositive, unordered, or not wider,
	 * or if the MAPQ cutoff is outside its supported range */
	public HybridMaxIndelConfig(int lowPrimary_, int lowSum_, int highPrimary_, int highSum_, int retryMinMapq_){
		if(lowPrimary_<1 || lowSum_<lowPrimary_){
			throw new IllegalArgumentException("Low max-indel bounds must be positive and ordered: "+
				lowPrimary_+"/"+lowSum_);
		}
		//BBIndex widens upper search bounds before addition; camelWalk caps its int horizon.
		//Seed penalties likewise widen before multiplication and then apply their existing cap.
		//These arithmetic safeguards do not validate all mapper behavior at extreme distances.
		if(highPrimary_<=lowPrimary_ || highSum_<=lowSum_ || highSum_<highPrimary_){
			throw new IllegalArgumentException("Retry max-indel bounds must be larger and ordered: low="+
				lowPrimary_+"/"+lowSum_+", high="+highPrimary_+"/"+highSum_);
		}
		if(retryMinMapq_<0 || retryMinMapq_>255){
			throw new IllegalArgumentException("Retry minimum MAPQ must be in [0,255]: "+retryMinMapq_);
		}
		lowPrimary=lowPrimary_; lowSum=lowSum_;
		highPrimary=highPrimary_; highSum=highSum_;
		retryMinMapq=retryMinMapq_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Describes initial and retry search bounds; intentionally omits the separate MAPQ cutoff. */
	@Override
	public String toString(){return lowPrimary+"/"+lowSum+"->"+highPrimary+"/"+highSum;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Initial maxindel distance in bases. */
	public final int lowPrimary;
	/** Initial maxindel2 distance in bases. */
	public final int lowSum;
	/** Retry maxindel distance in bases. */
	public final int highPrimary;
	/** Retry maxindel2 distance in bases. */
	public final int highSum;
	/** Additional acceptance cutoff; geometry and other guards still apply at zero. */
	public final int retryMinMapq;
}
