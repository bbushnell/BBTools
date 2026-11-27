package prok;

/** Member-local native gates over one combined index's already counted query.
 * This changes eligibility only; DDL ranking, the shared cap and stopping remain
 * the caller's responsibility. One instance belongs to one worker.
 * @author Brian, Raiden
 */
public final class NcrnaMemberNativeEligibility{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaMemberNativeEligibility(final NcrnaFamily[] families){
		if(families==null || families.length!=3){
			throw new IllegalArgumentException("Member eligibility needs the three ordered cobalamin families");
		}
		for(int member=0; member<3; member++){
			final NcrnaFamily family=families[member];
			if(family==null || !NAMES[member].equals(family.name) || family.library==null
				|| family.library.length<2 || family.indexK!=9 || family.indexTopN<family.library.length
				|| family.strictIndexCutoff || family.indexScoreMargin!=-1){
				throw new IllegalArgumentException("Member eligibility requires ordered multi-model K9 families whose separate caps retain every eligible model");
			}
			TrnaKmerIndex.validateConfiguration(9, family.adaptFloor, family.adaptTopFrac,
				family.adaptQFrac, family.fixedMinHits);
			ends[member+1]=Math.addExact(ends[member], family.library.length);
			adaptive[member]=family.adaptive;
			floor[member]=family.adaptFloor;
			topFraction[member]=family.adaptTopFrac;
			queryFraction[member]=family.adaptQFrac;
			fixed[member]=family.fixedMinHits;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Ascending global IDs, with each member's own maximum and keep-best fallback.
	 * TrnaKmerIndex.shortlist applies these gates after counting unique query words
	 * with model posting multiplicity. No second sequence scan/index is needed. */
	public int select(final TrnaKmerIndex index, final int[] destination){
		if(index==null || destination==null || destination.length<models()){
			throw new IllegalArgumentException("Member eligibility needs current native counts and scratch for the complete model axis");
		}
		int written=0;
		final int queryWords=Math.max(1, index.lastUniqueQueryKmers());
		for(int member=0; member<3; member++){
			int best=-1, maximum=0;
			for(int model=ends[member]; model<ends[member+1]; model++){
				final int count=index.lastSharedCount(model);
				if(count>maximum){
					maximum=count;
					best=model;
				}
			}
			final int cutoff=adaptive[member] ? (int)Math.ceil(Math.max(floor[member],
				Math.max(topFraction[member]*maximum, queryFraction[member]*queryWords))) : fixed[member];
			final int start=written;
			for(int model=ends[member]; model<ends[member+1]; model++){
				final int count=index.lastSharedCount(model);
				if(count>0 && count>=cutoff){destination[written++]=model;}
			}
			if(written==start && best>=0){destination[written++]=best;}
		}
		assert(written<=models()) : shared.KillSwitch.assertDie("Disjoint member ranges must fit the combined DDL model axis");
		return written;
	}
	public int models(){return ends[3];}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final int[] ends=new int[4], fixed=new int[3];
	private final boolean[] adaptive=new boolean[3];
	private final float[] floor=new float[3], topFraction=new float[3], queryFraction=new float[3];
	private static final String[] NAMES={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
}
