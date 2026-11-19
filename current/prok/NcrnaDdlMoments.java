package prok;

/** Equal-locus population moments from model/event summaries.
 * A reporter computes each event's model mean and second moment once, then
 * sums those event summaries in primitive arrays indexed by original parent.
 * This accumulator reduces each parent's sums with equal event weight and
 * gives every eligible parent exactly one weight, regardless of window count.
 * Raw model/event counts remain separate from the statistical denominator.
 * @author Brian, Raiden
 */
public final class NcrnaDdlMoments{

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	/** Model weights have already been normalized separately within each event.
	 * A paired event has both nonempty signal/background classes. The caller
	 * records exclusions and never passes an empty class here as a zero score. */
	public void addLocus(final double eventMeans, final double eventSecondMoments,
			final int eventCount, final long modelCount){
		require(eventCount>0 && modelCount>=eventCount && finite(eventMeans) && eventMeans>=0
			&& finite(eventSecondMoments) && eventSecondMoments>=0,
			"An eligible locus requires nonempty model/event classes and finite nonnegative bucket-hit moments");
		final double mean=eventMeans/eventCount;
		final double second=eventSecondMoments/eventCount;
		variance(mean, second);
		first+=mean;
		secondMoment+=second;
		loci++;
		events+=eventCount;
		models+=modelCount;
	}

	/** Population mean: every eligible locus contributes total weight one. */
	public double mean(){return loci==0 ? Double.NaN : first/loci;}

	/** Retains within-event model variation; not the spread of event means. */
	public double standardDeviation(){
		return loci==0 ? Double.NaN : Math.sqrt(variance(mean(), secondMoment/loci));
	}

	public long loci(){return loci;}
	public long events(){return events;}
	public long models(){return models;}

	/** Reset a reusable report cell; no per-model or per-event objects. */
	public void clear(){
		first=0;
		secondMoment=0;
		loci=0;
		events=0;
		models=0;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** The only clamp is negative roundoff at relative tolerance1e-12.
	 * A positive variance, however small, is never replaced by a floor. */
	private static double variance(final double mean, final double second){
		final double square=mean*mean;
		final double value=second-square;
		require(finite(square) && finite(value) && value>=-1e-12*Math.max(1, Math.max(second, square)),
			"Weighted second moments must be at least squared means; a material negative variance means the model/event/locus weights disagree");
		return Math.max(0, value);
	}

	private static boolean finite(final double value){return !Double.isNaN(value) && !Double.isInfinite(value);}
	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private double first, secondMoment;
	private long loci, events, models;
}
