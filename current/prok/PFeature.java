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
	
	public final void flip(){
		int a=scaflen-start-1;
		int b=scaflen-stop-1;
		start=b;
		stop=a;
		flipped=flipped^1;
	}
	
	public final int currentStrand(){
		return strand^flipped;
	}
	
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
	
	public final int flipped(){return flipped;}
	
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
	
	public final String scafName;
	public final int strand;
	public final int scaflen;
	
	/** 0-based position of first base of feature **/
	public int start;
	/** 0-based position of last base of feature **/
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
			return a.compareTo(b);
		}
		
	}
	
}
