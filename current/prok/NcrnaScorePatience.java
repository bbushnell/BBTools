package prok;

/** Default-off shared-walk experiment. Strict score improvements reset
 * stagnation whether or not the alignment passes. Expired patience stops once
 * the ranked shortlist threshold is reached OR a model has passed.
 * A worker resets this state separately for primary and physical-edge walks.
 * @author Ganyu
 */
public final class NcrnaScorePatience{

	public static final int CEILING=384;
	public static final int STAGNATION=1, HARD_CEILING=2, LIST_EXHAUSTED=3;

	public NcrnaScorePatience(final int shortlist_, final int patience_){
		if(shortlist_<0 || shortlist_>CEILING || patience_<1 || patience_==Integer.MAX_VALUE){
			throw new IllegalArgumentException("Score patience needs a ranked shortlist threshold in0..384 and a positive finite stagnation budget");
		}
		shortlist=shortlist_;
		patience=patience_;
		reset();
	}

	public void reset(){
		bestPresent=false;
		bestScore=0;
		sinceImproved=0;
		visited=0;
		eligible=0;
		lastObserved=0;
		stopped=false;
		passed=false;
	}

	/** Called before visiting each ranked entry, including entries rejected by
	 * the unchanged length prefilter. They do not consume the alignment budget.
	 * Checking here also honors an expired budget when the threshold ends in a skip. */
	public boolean stopBefore(final int rank){
		assert(!stopped && rank==visited+1 && rank<=CEILING)
			: "One ordered capped list must visit each rank once before alignment; rank="+rank+", visited="+visited;
		if(bestPresent && sinceImproved>=patience && (visited>=shortlist || passed)){
			stopped=true;
			return true;
		}
		visited=rank;
		return false;
	}

	/** Only actual span-eligible alignments call this, with the same score as
	 * the winner comparator. Negative scores and failing alignments are real. */
	public void observe(final double score, final boolean currentPass){
		assert(!stopped && visited>lastObserved && Double.isFinite(score))
			: "Each visited model contributes at most one finite span-eligible score";
		eligible++;
		lastObserved=visited;
		passed|=currentPass;
		if(!bestPresent || score>bestScore){
			bestPresent=true;
			bestScore=score;
			sinceImproved=0;
		}else{sinceImproved++;}
	}

	/** A walk reaching its complete list ends by exhaustion; rank384 is reported
	 * as ceiling even if stagnation also becomes sufficient on that final rank. */
	public int reason(final int listSize){
		assert(listSize>=0 && listSize<=CEILING && visited<=listSize
			&& (stopped ? visited<listSize : visited==listSize))
			: "The stop reason must describe the actual terminal position of this capped walk";
		return visited==CEILING ? HARD_CEILING : stopped ? STAGNATION : LIST_EXHAUSTED;
	}
	public int visited(){return visited;}
	public int eligible(){return eligible;}
	public boolean ceilingReached(){return visited==CEILING;}
	public static String reasonName(final int reason){
		switch(reason){
			case STAGNATION:return "stagnation";
			case HARD_CEILING:return "ceiling";
			case LIST_EXHAUSTED:return "list_exhausted";
			default:throw new IllegalArgumentException("Unknown score-patience stop reason: "+reason);
		}
	}
	public final int shortlist, patience;
	private int sinceImproved, visited, eligible, lastObserved;
	private double bestScore;
	private boolean bestPresent, stopped, passed;
}
