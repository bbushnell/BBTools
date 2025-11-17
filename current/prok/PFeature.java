package prok;

import java.util.Comparator;

import structures.Feature;

/**
 * Represents a genomic feature such as a gene, with start, stop, and strand.
 * @author Brian Bushnell
 * @date Sep 24, 2018
 *
 */
abstract class PFeature extends ProkObject implements Comparable<PFeature>, Feature {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Constructs a genomic feature with coordinate and strand information.
	 *
	 * @param scafName_ Scaffold name containing this feature
	 * @param start_ 0-based start position of feature
	 * @param stop_ 0-based stop position of feature (inclusive)
	 * @param strand_ Strand orientation (0=forward, 1=reverse)
	 * @param scaflen_ Total length of the containing scaffold
	 */
	public PFeature(String scafName_, int start_, int stop_, int strand_, int scaflen_){
		scafName=scafName_;
		start=start_;
		stop=stop_;
		strand=strand_;
		scaflen=scaflen_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Flips the feature coordinates to the opposite strand.
	 * Recalculates start and stop positions relative to scaffold end
	 * and toggles the flipped state.
	 */
	public final void flip(){
		int a=scaflen-start-1;
		int b=scaflen-stop-1;
		start=b;
		stop=a;
		flipped=flipped^1;
	}
	
	/** Gets the current effective strand after any flipping operations.
	 * @return Current strand (original strand XOR flipped state) */
	public final int currentStrand(){
		return strand^flipped;
	}
	
	/** Calculates the length of this feature in bases.
	 * @return Feature length (stop - start + 1) */
	public final int length(){
		return stop-start+1;
	}
	
	@Override
	public final int compareTo(PFeature f) {
		int x=scafName.compareTo(f.scafName);
		if(x!=0){return x;}
		if(stop!=f.stop){return stop-f.stop;}
		return start-f.start;
	}
	
	/** Gets the flipped state of this feature */
	public final int flipped(){return flipped;}
	
	/**
	 * Gets the score for this feature.
	 * Implementation varies by feature type.
	 * @return Feature-specific score value
	 */
	public abstract float score();
	
	@Override
	public final int start() {return start;}
	
	@Override
	public final int stop() {return stop;}
	
	@Override
	public final int strand() {return strand;}
	
	@Override
	public final String name() {return null;}
	
	@Override
	public final String seqid() {return scafName;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Name of the scaffold containing this feature */
	public final String scafName;
	/** Original strand orientation (0=forward, 1=reverse) */
	public final int strand;
	/** Total length of the containing scaffold in bases */
	public final int scaflen;
	
	/** 0-based position of first base of feature **/
	public int start;
	/** 0-based position of last base of feature **/
	public int stop;
	/** Tracks whether coordinates have been flipped (0=not flipped, 1=flipped) */
	private int flipped=0;
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Static comparator instance for sorting features by score (high to low) */
	@SuppressWarnings("synthetic-access")
	public static final FeatureComparatorScore featureComparatorScore=new FeatureComparatorScore();
	
	//Sorts so that high scores are first.
	/** Comparator for sorting features by score in descending order.
	 * Higher scoring features are ordered first, with positional comparison as tiebreaker. */
	private static class FeatureComparatorScore implements Comparator<PFeature> {

		/** Private constructor to enforce singleton pattern */
		private FeatureComparatorScore(){}
		
		@Override
		public int compare(PFeature a, PFeature b) {
			float sa=a.score(), sb=b.score();
			if(sa<sb){return 1;}
			if(sb<sa){return -1;}
			return a.compareTo(b);
		}
		
	}
	
}
