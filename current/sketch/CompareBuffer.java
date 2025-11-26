package sketch;

import shared.Tools;
import structures.AbstractBitSet;

/**
 * Tracks and calculates various hit and comparison metrics for k-mer sequence analysis.
 * Manages computational buffer for tracking sequence hits, depth, and similarity metrics
 * across different k-mer lengths with optional BitSet allocation for memory efficiency.
 *
 * @author Brian Bushnell
 * @date September 1, 2025
 */
public class CompareBuffer extends SketchObject{
	
	/** Constructs a CompareBuffer with optional BitSet allocation.
	 * @param makeBS Whether to create an AbstractBitSet for k-mer tracking */
	public CompareBuffer(boolean makeBS){
		if(makeBS){
			cbs=AbstractBitSet.make(0, bitSetBits);
		}else{
			cbs=null;
		}
	}
	
	/**
	 * Sets all comparison metrics in a single operation.
	 * @param hits_ Total k-mer hits between query and reference
	 * @param multiHits_ Number of k-mers with multiple hits
	 * @param unique2_ Number of unique k-mers with depth 2
	 * @param unique3_ Number of unique k-mers with depth 3
	 * @param noHits_ Number of query k-mers with no reference hits
	 * @param contamHits_ Number of contamination hits detected
	 * @param contam2Hits_ Number of contamination hits with depth 2
	 * @param multiContamHits_ Number of multi-occurrence contamination hits
	 * @param queryDivisor_ Total number of valid k-mers in query
	 * @param refDivisor_ Total number of valid k-mers in reference
	 * @param querySize_ Size of query sketch
	 * @param refSize_ Size of reference sketch
	 * @param depthSum_ Sum of hit depths for coverage calculation
	 * @param depthSum2_ Sum of squared hit depths
	 * @param refHitSum_ Sum of reference hits for average calculation
	 * @param k1hits_ Number of hits for first k-mer size
	 * @param k1seenQ_ Number of k-mers seen in query for first k-size
	 * @param k1seenR_ Number of k-mers seen in reference for first k-size
	 */
	void set(final int hits_, final int multiHits_, final int unique2_, final int unique3_, final int noHits_,
			final int contamHits_, final int contam2Hits_, final int multiContamHits_,
			final int queryDivisor_, final int refDivisor_, final int querySize_, final int refSize_, 
			final long depthSum_, final double depthSum2_, final long refHitSum_,
			final int k1hits_, final int k1seenQ_, final int k1seenR_){
		hits=hits_;
		multiHits=multiHits_;
		unique2=unique2_;
		unique3=unique3_;
		noHits=noHits_;
		
		contamHits=contamHits_;
		contam2Hits=contam2Hits_;
		multiContamHits=multiContamHits_;
		
		queryDivisor=queryDivisor_;
		refDivisor=refDivisor_;
		
		querySize=querySize_;
		refSize=refSize_;

		depthSum=depthSum_;
		depthSum2=(float)depthSum2_;
		refHitSum=refHitSum_;

		hits1=k1hits_;
		qSeen1=k1seenQ_;
		rSeen1=k1seenR_;
	}
	
	/** Resets all comparison metrics to zero. */
	void clear(){
		hits=multiHits=0;
		unique2=unique3=noHits=0;
		contamHits=contam2Hits=multiContamHits=0;
		refDivisor=queryDivisor=0;
		refSize=querySize=0;
		depthSum=0;
		depthSum2=0;
		refHitSum=0;
		hits1=qSeen1=rSeen1=0;
	}
	
	/** Calculates average depth per hit.
	 * @return Average depth (depthSum/hits), or 0 if no hits */
	float depth(){
		return depthSum<1 ? 0 : depthSum/Tools.max(1.0f, hits);
	}
	
	/** Calculates average squared depth per hit.
	 * @return Average squared depth (depthSum2/hits), or 0 if no hits */
	float depth2(){
		return depthSum2<=0 ? 0 : depthSum2/Tools.max(1.0f, hits);
	}
	
	/** Calculates average reference hits per query hit.
	 * @return Average reference hits (refHitSum/hits), or 0 if no hits */
	float avgRefHits(){
		return refHitSum<1 ? 0 : refHitSum/Tools.max(1.0f, hits);
	}
	
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return "hits="+hits+", refDivisor="+refDivisor+", queryDivisor="+queryDivisor+", refSize="+refSize+", querySize="+querySize+
				", contamHits="+contamHits+", contam2Hits="+contam2Hits+", multiContamHits="+multiContamHits+", depthSum="+depthSum+", depthSum2="+depthSum2+
				", hits="+hits+", multiHits="+multiHits+", unique2="+unique2+", unique3="+unique3+", noHits="+noHits;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Calculates weighted k-mer identity using minimum divisor.
	 * @return Weighted k-mer identity (hits/minDivisor) */
	public final float wkid(){
		final int div=minDivisor();
		return hits/(float)div;
	}
	/** Calculates k-mer identity using maximum divisor.
	 * @return K-mer identity (hits/maxDivisor) */
	final float kid(){
		final int div=maxDivisor();
		return hits/(float)div;
	}
	/** Calculates legacy Average Nucleotide Identity from weighted k-mer identity.
	 * @return ANI value converted from wkid */
	final float aniOld(){
		float wkid=wkid();
		final float ani=wkidToAni(wkid);
		return ani;
	}
	/**
	 * Calculates Average Nucleotide Identity using dual k-mer analysis when available.
	 * Uses combined ani1/ani2 calculation if k2>0, otherwise falls back to aniOld.
	 * @return ANI value optimized using available k-mer data
	 */
	public final float ani(){
		final float ani;
		if(k2>0 && useToValue2){
			float ani1=ani1();
			float ani2=ani2();
//			ani=0.5f*(ani1+ani2);
			ani=0.5f*(Tools.max(0.9f*ani2, ani1)+Tools.max(0.8f*ani1, ani2));
//			return (ani1*qSeen1+ani2*qSeen2())/queryDivisor;
		}else{
			ani=aniOld();
		}
		
//		System.err.println("ani="+ani+"aniOld="+aniOld()+", ani1="+ani1()+", ani2="+ani2()+", anid="+(float)aniDual()+"\n"
////				+"gf="+(float)gf+", wkid1="+wkid1+", wkid2="+wkid2+"\n"
//						+ "k1f="+k1Fraction()+", hits="+hits+", hits1="+hits1+", hits2="+hits2()+", qSeen1()="+qSeen1()+", rSeen1()="+rSeen1()+"\n"
//								+ "qSeen2()="+qSeen2()+", rSeen2()="+rSeen2()+", minDivisor1()="+minDivisor1()+", minDivisor2()="+minDivisor2()+"\n");
		return ani;
	}

	/** Calculates weighted k-mer identity for first k-mer size.
	 * @return Weighted k-mer identity using minimum divisor for k1 */
	final float wkid1(){
		final int div=minDivisor1();
		return hits1()/(float)div;
	}
	/** Calculates k-mer identity for first k-mer size using maximum divisor.
	 * @return K-mer identity for k1 using maximum divisor */
	final float kid1(){
		final int div=maxDivisor1();
		return hits1()/(float)div;
	}
	/** Calculates Average Nucleotide Identity for first k-mer size.
	 * @return ANI value converted from wkid1 using exact conversion for k */
	final float ani1(){
		float wkid=wkid1();
		final float ani=wkidToAniExact(wkid, k);
		return ani;
	}

	/** Calculates weighted k-mer identity for second k-mer size.
	 * @return Weighted k-mer identity using minimum divisor for k2 */
	final float wkid2(){
		final int div=minDivisor2();
		return hits2()/(float)div;
	}
	/** Calculates k-mer identity for second k-mer size using maximum divisor.
	 * @return K-mer identity for k2 using maximum divisor */
	final float kid2(){
		final int div=maxDivisor2();
		return hits2()/(float)div;
	}
	/** Calculates Average Nucleotide Identity for second k-mer size.
	 * @return ANI value converted from wkid2 using exact conversion for k2 */
	final float ani2(){
		assert(k2>0);
		float wkid=wkid2();
		final float ani=wkidToAniExact(wkid, k2);
		return ani;
	}
	
	/**
	 * Calculates dual k-mer Average Nucleotide Identity using mathematical relationship
	 * between different k-mer sizes. Uses power relationship to estimate ANI from
	 * wkid1/wkid2 ratio.
	 * @return ANI estimate using dual k-mer mathematical relationship
	 */
	final float aniDual(){
		assert(k2>0);
		float wkid1=wkid1();
		float wkid2=wkid2();
		float ratio=(wkid1/wkid2);
		float exp=1f/(k-k2);//TODO - make this initialized
		double ani=Math.pow(ratio, exp);
		double gf=wkid2/Math.pow(ani, k2);
		
//		System.err.println("ani="+ani()+"aniOld="+aniOld()+", ani1="+ani1()+", ani2="+ani2()+", anid="+(float)ani+"\n"
//				+"gf="+(float)gf+", wkid1="+wkid1+", wkid2="+wkid2+"\n"
//						+ "k1f="+k1Fraction()+", hits="+hits+", hits1="+hits1+", hits2="+hits2()+", qSeen1()="+qSeen1()+", rSeen1()="+rSeen1()+"\n"
//								+ "qSeen2()="+qSeen2()+", rSeen2()="+rSeen2()+", minDivisor1()="+minDivisor1()+", minDivisor2()="+minDivisor2()+"\n");
		
		return (float)ani;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Returns total k-mer hits between query and reference. */
	public int hits(){return hits;}
	/** Returns number of k-mers with multiple hits. */
	int multiHits(){return multiHits;}
	/** Returns number of query k-mers with no reference hits. */
	int noHits(){return noHits;}
	/** Returns number of unique k-mers with depth 2. */
	int unique2(){return unique2;}
	/** Returns number of unique k-mers with depth 3. */
	int unique3(){return unique3;}

	/** Returns number of contamination hits detected. */
	int contamHits(){return contamHits;}
	/** Returns number of contamination hits with depth 2. */
	int contam2Hits(){return contam2Hits;}
	/** Returns number of multi-occurrence contamination hits. */
	int multiContamHits(){return multiContamHits;}
	
	/** Returns total number of valid k-mers in query. */
	int queryDivisor(){return queryDivisor;}
	/** Returns total number of valid k-mers in reference. */
	int refDivisor(){return refDivisor;}
	
	/** Returns size of query sketch. */
	int querySize(){return querySize;}
	/** Returns size of reference sketch. */
	int refSize(){return refSize;}

	/** Returns sum of hit depths for coverage calculation. */
	long depthSum(){return depthSum;}
	/** Returns sum of squared hit depths. */
	float depthSum2(){return depthSum2;}
	/** Returns sum of reference hits for average calculation. */
	long refHitSum(){return refHitSum;}
	
	/*--------------------------------------------------------------*/

	/** Returns number of hits for first k-mer size. */
	int hits1(){return hits1;}
	/** Returns number of k-mers seen in query for first k-size. */
	int qSeen1(){return qSeen1;}
	/** Returns number of k-mers seen in reference for first k-size. */
	int rSeen1(){return rSeen1;}
	/** Returns minimum divisor for first k-mer size calculations.
	 * @return Minimum of qSeen1 and rSeen1, at least 1 */
	int minDivisor1(){return Tools.max(1, Tools.min(qSeen1, rSeen1));}
	/** Returns maximum divisor for first k-mer size calculations.
	 * @return Maximum of qSeen1 and rSeen1, at least 1 */
	int maxDivisor1(){return Tools.max(1, qSeen1, rSeen1);}

	/** Returns number of hits for second k-mer size.
	 * @return Difference between total hits and first k-mer hits */
	int hits2(){return hits-hits1;}
	/** Returns number of k-mers seen in query for second k-size.
	 * @return Difference between queryDivisor and qSeen1 */
	int qSeen2(){return queryDivisor-qSeen1;}
	/** Returns number of k-mers seen in reference for second k-size.
	 * @return Difference between refDivisor and rSeen1 */
	int rSeen2(){return refDivisor-rSeen1;}
	/** Returns minimum divisor for second k-mer size calculations.
	 * @return Minimum of qSeen2 and rSeen2, at least 1 */
	int minDivisor2(){return Tools.max(1, Tools.min(qSeen2(), rSeen2()));}
	/** Returns maximum divisor for second k-mer size calculations.
	 * @return Maximum of qSeen2 and rSeen2, at least 1 */
	int maxDivisor2(){return Tools.max(1, qSeen2(), rSeen2());}
	
	/*--------------------------------------------------------------*/

	//For WKID
	/** Returns minimum divisor for weighted k-mer identity calculations.
	 * @return Minimum of queryDivisor and refDivisor, at least 1 */
	int minDivisor(){return Tools.max(1, Tools.min(queryDivisor, refDivisor));}
	//For KID
	/** Returns maximum divisor for k-mer identity calculations.
	 * @return Maximum of queryDivisor and refDivisor, at least 1 */
	int maxDivisor(){return Tools.max(1, queryDivisor, refDivisor);}
	/** Returns minimum sketch size.
	 * @return Minimum of querySize and refSize, at least 1 */
	int minSize(){return Tools.max(1, Tools.min(querySize, refSize));}
	/** Returns maximum sketch size.
	 * @return Maximum of querySize and refSize, at least 1 */
	int maxSize(){return Tools.max(1, querySize, refSize);}

	/** Returns number of unique hits (non-multi-hit k-mers).
	 * @return Total hits minus multi-hits */
	int uniqueHits(){return hits-multiHits;}
	/** Returns number of unique contamination hits.
	 * @return Contamination hits minus multi-contamination hits */
	int uniqueContamHits(){return contamHits-multiContamHits;}
	
	/** Calculates fraction of query k-mers seen in first k-mer size.
	 * @return Ratio of qSeen1 to total queryDivisor */
	float k1Fraction(){
		return qSeen1/Tools.max(queryDivisor, 1f);
	}
	
	/*--------------------------------------------------------------*/
	
	/** Total k-mer hits between query and reference. */
	private int hits;
	/** Number of k-mers with multiple hits. */
	private int multiHits;
	/** Number of query k-mers with no reference hits. */
	private int noHits;
	/** Number of unique k-mers with depth 2. */
	private int unique2;
	/** Number of unique k-mers with depth 3. */
	private int unique3;

	/** Number of contamination hits detected. */
	private int contamHits;
	/** Number of contamination hits with depth 2. */
	private int contam2Hits;
	/** Number of multi-occurrence contamination hits. */
	private int multiContamHits;
	
	/** Total number of valid k-mers in query. */
	private int queryDivisor;
	/** Total number of valid k-mers in reference. */
	private int refDivisor;
	
	/** Size of query sketch. */
	private int querySize;
	/** Size of reference sketch. */
	private int refSize;

	/** Sum of hit depths for coverage calculation. */
	private long depthSum;
	/** Sum of squared hit depths. */
	private float depthSum2;
	/** Sum of reference hits for average calculation. */
	private long refHitSum;

	/** Number of hits for first k-mer size. */
	private int hits1;
	/** Number of k-mers seen in query for first k-size. */
	private int qSeen1;
	/** Number of k-mers seen in reference for first k-size. */
	private int rSeen1;
	
	/*--------------------------------------------------------------*/

	/**
	 * Optional BitSet for k-mer tracking, used only for comparisons not indexing.
	 */
	public final AbstractBitSet cbs; //Only for comparisons, not index
	
}
