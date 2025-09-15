package sketch;

import shared.Tools;
import structures.AbstractBitSet;

public class CompareBuffer extends SketchObject{
	
	public CompareBuffer(boolean makeBS){
		if(makeBS){
			cbs=AbstractBitSet.make(0, bitSetBits);
		}else{
			cbs=null;
		}
	}
	
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
	
	float depth(){
		return depthSum<1 ? 0 : depthSum/Tools.max(1.0f, hits);
	}
	
	float depth2(){
		return depthSum2<=0 ? 0 : depthSum2/Tools.max(1.0f, hits);
	}
	
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
	
	final float wkid(){
		final int div=minDivisor();
		return hits/(float)div;
	}
	final float kid(){
		final int div=maxDivisor();
		return hits/(float)div;
	}
	final float aniOld(){
		float wkid=wkid();
		final float ani=wkidToAni(wkid);
		return ani;
	}
	final float ani(){
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

	final float wkid1(){
		final int div=minDivisor1();
		return hits1()/(float)div;
	}
	final float kid1(){
		final int div=maxDivisor1();
		return hits1()/(float)div;
	}
	final float ani1(){
		float wkid=wkid1();
		final float ani=wkidToAniExact(wkid, k);
		return ani;
	}

	final float wkid2(){
		final int div=minDivisor2();
		return hits2()/(float)div;
	}
	final float kid2(){
		final int div=maxDivisor2();
		return hits2()/(float)div;
	}
	final float ani2(){
		assert(k2>0);
		float wkid=wkid2();
		final float ani=wkidToAniExact(wkid, k2);
		return ani;
	}
	
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
	
	int hits(){return hits;}
	int multiHits(){return multiHits;}
	int noHits(){return noHits;}
	int unique2(){return unique2;}
	int unique3(){return unique3;}

	int contamHits(){return contamHits;}
	int contam2Hits(){return contam2Hits;}
	int multiContamHits(){return multiContamHits;}
	
	int queryDivisor(){return queryDivisor;}
	int refDivisor(){return refDivisor;}
	
	int querySize(){return querySize;}
	int refSize(){return refSize;}

	long depthSum(){return depthSum;}
	float depthSum2(){return depthSum2;}
	long refHitSum(){return refHitSum;}
	
	/*--------------------------------------------------------------*/

	int hits1(){return hits1;}
	int qSeen1(){return qSeen1;}
	int rSeen1(){return rSeen1;}
	int minDivisor1(){return Tools.max(1, Tools.min(qSeen1, rSeen1));}
	int maxDivisor1(){return Tools.max(1, qSeen1, rSeen1);}

	int hits2(){return hits-hits1;}
	int qSeen2(){return queryDivisor-qSeen1;}
	int rSeen2(){return refDivisor-rSeen1;}
	int minDivisor2(){return Tools.max(1, Tools.min(qSeen2(), rSeen2()));}
	int maxDivisor2(){return Tools.max(1, qSeen2(), rSeen2());}
	
	/*--------------------------------------------------------------*/

	//For WKID
	int minDivisor(){return Tools.max(1, Tools.min(queryDivisor, refDivisor));}
	//For KID
	int maxDivisor(){return Tools.max(1, queryDivisor, refDivisor);}
	int minSize(){return Tools.max(1, Tools.min(querySize, refSize));}
	int maxSize(){return Tools.max(1, querySize, refSize);}

	int uniqueHits(){return hits-multiHits;}
	int uniqueContamHits(){return contamHits-multiContamHits;}
	
	float k1Fraction(){
		return qSeen1/Tools.max(queryDivisor, 1f);
	}
	
	/*--------------------------------------------------------------*/
	
	private int hits;
	private int multiHits;
	private int noHits;
	private int unique2;
	private int unique3;

	private int contamHits;
	private int contam2Hits;
	private int multiContamHits;
	
	private int queryDivisor;
	private int refDivisor;
	
	private int querySize;
	private int refSize;

	private long depthSum;
	private float depthSum2;
	private long refHitSum;

	private int hits1;
	private int qSeen1;
	private int rSeen1;
	
	/*--------------------------------------------------------------*/

	public final AbstractBitSet cbs; //Only for comparisons, not index
	
}
