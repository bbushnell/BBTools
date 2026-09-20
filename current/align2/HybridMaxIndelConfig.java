package align2;

/** Immutable per-invocation bounds for selective max-indel retry.
 * @author Collei
 */
public final class HybridMaxIndelConfig {

	public HybridMaxIndelConfig(int lowPrimary_,int lowSum_,int highPrimary_,int highSum_){
		this(lowPrimary_,lowSum_,highPrimary_,highSum_,0);
	}

	public HybridMaxIndelConfig(int lowPrimary_,int lowSum_,int highPrimary_,int highSum_,int retryMinMapq_){
		if(lowPrimary_<1 || lowSum_<lowPrimary_){
			throw new IllegalArgumentException("Low max-indel bounds must be positive and ordered: "+
				lowPrimary_+"/"+lowSum_);
		}
		if(highPrimary_<=lowPrimary_ || highSum_<=lowSum_ || highSum_<highPrimary_){
			throw new IllegalArgumentException("Retry max-indel bounds must be larger and ordered: low="+
				lowPrimary_+"/"+lowSum_+", high="+highPrimary_+"/"+highSum_);
		}
		if(retryMinMapq_<0 || retryMinMapq_>255){
			throw new IllegalArgumentException("Retry minimum MAPQ must be in [0,255]: "+retryMinMapq_);
		}
		lowPrimary=lowPrimary_;lowSum=lowSum_;
		highPrimary=highPrimary_;highSum=highSum_;
		retryMinMapq=retryMinMapq_;
	}

	@Override
	public String toString(){return lowPrimary+"/"+lowSum+"->"+highPrimary+"/"+highSum;}

	public final int lowPrimary;
	public final int lowSum;
	public final int highPrimary;
	public final int highSum;
	public final int retryMinMapq;
}
