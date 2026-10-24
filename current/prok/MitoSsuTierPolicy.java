package prok;

/** Coverage-library controls selected from matched genomic accuracy/cost runs.
 * Original16S verification still precedes this mitochondrial-only policy.
 * @author Brian, Raiden */
public final class MitoSsuTierPolicy {
	public MitoSsuTierPolicy(final int maxModels_, final int minHits_, final float queryFraction_, final float topFraction_){
		if(maxModels_<1 || minHits_<0 || !Float.isFinite(queryFraction_) || queryFraction_<0 || queryFraction_>1
			|| !Float.isFinite(topFraction_) || topFraction_<0 || topFraction_>1){
			throw new IllegalArgumentException("Mito tier requires positive model cap, nonnegative hit floor and fractions in[0,1]");
		}
		maxModels=maxModels_;minHits=minHits_;queryFraction=queryFraction_;topFraction=topFraction_;
	}
	public static MitoSsuTierPolicy configured(){
		return new MitoSsuTierPolicy(Integer.parseInt(System.getProperty("bbtools.mito.maxModels","3")),
			Integer.parseInt(System.getProperty("bbtools.mito.minHits","80")),
			Float.parseFloat(System.getProperty("bbtools.mito.queryFraction","0")),
			Float.parseFloat(System.getProperty("bbtools.mito.topFraction","0")));
	}
	public int entryFloor(final int queryUnique){
		assert(queryUnique>=0):"The entry fraction uses the actual deduplicated window9-mer count";
		return Math.max(minHits,(int)Math.ceil((double)queryFraction*queryUnique));
	}
	public int subsequentFloor(final int topShared, final int queryUnique){
		assert(topShared>=0):"Relative eligibility is based on this window's top mito score";
		return Math.max(entryFloor(queryUnique),(int)Math.ceil((double)topFraction*topShared));
	}
	public final int maxModels,minHits;
	public final float queryFraction,topFraction;
}
