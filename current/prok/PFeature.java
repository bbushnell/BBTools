package prok;

import java.util.Comparator;

import structures.Feature;

/**
 * Abstract base class for genomic features with coordinate and strand management.
 * Provides coordinate flipping between forward and reverse strands for prokaryotic
 * sequence processing and implements consistent feature ordering by position.
 * @author Brian Bushnell
 * @date Sep 24, 2018
 */
abstract class PFeature extends ProkObject implements Comparable<PFeature>, Feature {
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/

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
	
	/** Reverse-strands coordinates in place: mirror each (p-&gt;scaflen-1-p), swap so the smaller
	 * mirrored value stays in start, toggle flip parity. Assumes 0&lt;=start&lt;=stop&lt;scaflen (caller
	 * contract, not enforced here). */
	public final void flip(){
		int a=scaflen-start-1;
		int b=scaflen-stop-1;
		start=b;//cross-assign: start<=stop => b<=a, so b->start keeps start<=stop
		stop=a;
		flipped=flipped^1;//flipped in {0,1}: init 0, only toggled here
	}
	
	/** Effective strand after flips: original strand XOR flip-parity (an odd number of flips
	 * inverts it). Uses prok's 0/1 convention (0=fwd, 1=rev). NOTE: the Feature interface javadoc
	 * documents +1/-1/0, which disagrees with this -- see the TODO. */
	public final int currentStrand(){
		//TODO: Possible bug [prok/PFeature#001] - assumes strand in {0,1}; Feature.strand() doc says +1/-1/0. Confirm prok's convention via Orf/GeneCaller; likely a wrong interface doc, not a logic bug here.
		return strand^flipped;
	}
	
	public final int length(){
		return stop-start+1;
	}
	
	@Override
	public final int compareTo(PFeature f) {
		int x=scafName.compareTo(f.scafName);
		if(x!=0){return x;}
		//ascending stop then start; subtraction overflow-safe for well-formed prok coords (small genomes), not enforced
		if(stop!=f.stop){return stop-f.stop;}
		return start-f.start;
	}
	
	public final int flipped(){return flipped;}
	
	public abstract float score();
	
	@Override
	public final int start() {return start;}
	
	@Override
	public final int stop() {return stop;}
	
	@Override
	public final int strand() {return strand;}
	
	/** Gets the feature name (always returns null for base implementation) */
	@Override
	public final String name() {return null;}
	
	@Override
	public final String seqid() {return scafName;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final String scafName;
	public final int strand;
	public final int scaflen;
	
	public int start;
	public int stop;
	private int flipped=0;
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	@SuppressWarnings("synthetic-access")
	public static final FeatureComparatorScore featureComparatorScore=new FeatureComparatorScore();
	
	//Sorts so that high scores are first.
	private static class FeatureComparatorScore implements Comparator<PFeature> {

		private FeatureComparatorScore(){}
		
		@Override
		public int compare(PFeature a, PFeature b) {
			float sa=a.score(), sb=b.score();
			if(sa<sb){return 1;}
			if(sb<sa){return -1;}
			//NaN reaches here (both '<' comparisons false): positional tie via compareTo, no throw
			return a.compareTo(b);
		}
		
	}
	
}
