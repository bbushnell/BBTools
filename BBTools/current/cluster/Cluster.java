package cluster;

import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicLongArray;

import stream.Read;

/**
 * @author Brian Bushnell
 * @date Mar 24, 2014
 *
 */
public class Cluster{
	
	public Cluster(int id_, int k1_, int k2_, int arraylen1_, int arraylen2_){
		
		id=id_;
		k1=k1_;
		k2=k2_;
		arraylen1=arraylen1_;
		arraylen2=arraylen2_;
		
		kmerArray1=new AtomicLongArray(arraylen1);
		kmerProbArray1=new float[arraylen1];

		kmerArray2=new AtomicLongArray(arraylen2);
		kmerProbArray2=new float[arraylen2];
	}
	
	/*--------------------------------------------------------------*/

	public void recalculate(){
		gc=(float)(gcCount.doubleValue()/baseCount.doubleValue());

		if(k1>0){
			long kmerCount=0;
			for(int i=0; i<arraylen1; i++){
				kmerCount+=kmerArray1.get(i);
			}
			double extra=(0.05/arraylen1);
			double mult=(0.95/kmerCount);
			for(int i=0; i<arraylen1; i++){
				kmerProbArray1[i]=(float)(kmerArray1.get(i)*mult+extra);
			}
		}
		if(k2>0){
			long kmerCount=0;
			for(int i=0; i<arraylen2; i++){
				kmerCount+=kmerArray2.get(i);
			}
			double extra=(0.05/arraylen2);
			double mult=(0.95/kmerCount);
			for(int i=0; i<arraylen2; i++){
				kmerProbArray2[i]=(float)(kmerArray2.get(i)*mult+extra);
			}
		}
	}
	
	public void resetAtomics(){
		for(int i=0; i<arraylen1; i++){
			kmerArray1.set(i, 0);
		}
		for(int i=0; i<arraylen2; i++){
			kmerArray2.set(i, 0);
		}
		depthsum1.set(0);
		depthsum2.set(0);
		readCount.set(0);
		baseCount.set(0);
		gcCount.set(0);
	}
	
	public void add(Read r){
		if(r==null){return;}
		ReadTag rt=(ReadTag)r.obj;
		assert(rt!=null);
		final byte[] bases=r.bases;
		
		readCount.incrementAndGet();
		baseCount.addAndGet(bases.length);
		gcCount.addAndGet(rt.gcCount);
		
		if(rt.strand==0){
			depthsum1.addAndGet(rt.depth);
		}else{
			depthsum2.addAndGet(rt.depth);
		}
		
		if(k1>0){
			int[] kmers=rt.kmerArray1(k1);
			int kmer=-1, run=0;
			for(int i=0; i<kmers.length; i++){
				int x=kmers[i];
				if(x==kmer){
					run++;
				}else{
					if(run>0){kmerArray1.addAndGet(kmer, run);}
					kmer=x;
					run=1;
				}
			}
			if(run>0){kmerArray1.addAndGet(kmer, run);}
		}

		if(k2>0){
			int[] kmers=rt.kmerArray2(k2);
			for(int kmer=0; kmer<kmers.length; kmer++){
				int x=kmers[kmer];
				if(x>0){kmerArray2.addAndGet(kmer, x);}
			}
		}
	}
	
	public float score(Read r) {
		if(r==null){return 0;}
		return r.mate==null ? scoreSingle(r) : scorePaired(r);
	}
	
	public float scoreSingle(Read r) {
		if(r==null){return 0;}
		ReadTag rt=(ReadTag)r.obj;
		
		assert(false) : "TODO";
		float depthScore=scoreDepthSingle(rt);
		float gcScore=scoreGcSingle(rt);
		float kmerScore=scoreKmer1(rt);
		assert(false);
		float depthWeight=.2f;
		float gcWeight=.2f;
		float kmerWeight=.6f;
		
		return depthWeight*depthScore+gcWeight*gcScore+kmerWeight*kmerScore;
	}
	
	/**
	 * @param rt
	 * @return
	 */
	private float scoreKmer1(ReadTag rt) {
		int[] kmers=rt.kmerArray1(k1);
		
		float score=0;
		if(scoreMode1==SCORE_MODE_AND){
			float f=ClusterTools.andCount(kmers, kmerArray1);
			assert(false);
		}else if(scoreMode1==SCORE_MODE_MULT){
			float f=ClusterTools.innerProduct(kmers, kmerProbArray1);
			assert(false);
		}else{
			throw new RuntimeException(""+scoreMode1);
		}
		
		return score;
	}
	
	/**
	 * @param rt
	 * @return
	 */
	private float scoreKmer2(ReadTag rt) {
		int[] kmers=rt.kmerArray2(k2);
		float[] probs=rt.kmerFreq2(k2);
		
		float score=0;
		if(scoreMode2==SCORE_MODE_AND){
			float f=ClusterTools.andCount(kmers, kmerArray2);
			assert(false);
		}else if(scoreMode2==SCORE_MODE_MULT){
			float f=ClusterTools.innerProduct(kmers, kmerProbArray2);
			assert(false);
		}else if(scoreMode2==SCORE_MODE_DIF){
			float f=ClusterTools.absDif(probs, kmerProbArray2);
			assert(false);
		}else if(scoreMode2==SCORE_MODE_RMS){
			float f=ClusterTools.rmsDif(probs, kmerProbArray2);
			assert(false);
		}else if(scoreMode2==SCORE_MODE_KS){
			float f=ClusterTools.ksFunction(probs, kmerProbArray2);
			assert(false);
		}else{
			throw new RuntimeException(""+scoreMode2);
		}
		
		return score;
	}

	/**
	 * @param rt
	 * @return
	 */
	private float scoreGcSingle(ReadTag rt) {
		assert(false) : "TODO";
		// TODO Auto-generated method stub
		return 0;
	}

	/**
	 * @param rt
	 * @return
	 */
	private float scoreDepthSingle(ReadTag rt) {
		assert(false) : "TODO";
		// TODO Auto-generated method stub
		return 0;
	}
	
	public float scorePaired(Read r) {
		assert(false) : "TODO";
		if(r==null){return 0;}
		ReadTag rt=(ReadTag)r.obj;
		
//		ReadTag rt1=rt.r
		
		return 0;
	}
	
	/*--------------------------------------------------------------*/
	
	public final int id;
	
	/** 'big' kmer */
	public final int k1;
	/** 'small' kmer */
	public final int k2;

	public final int arraylen1;
	public final int arraylen2;
	
	/*--------------------------------------------------------------*/
	
	public float gc;
	public int depth1, depth2;
	
	final AtomicLongArray kmerArray1;
	final float[] kmerProbArray1;
	
	final AtomicLongArray kmerArray2;
	final float[] kmerProbArray2;
	
	final AtomicLong depthsum1=new AtomicLong(0);
	final AtomicLong depthsum2=new AtomicLong(0);
	
	final AtomicLong readCount=new AtomicLong(0);
	final AtomicLong baseCount=new AtomicLong(0);
//	final AtomicLong kmerCount=new AtomicLong(0);
	final AtomicLong gcCount=new AtomicLong(0);
	
	/*--------------------------------------------------------------*/

	public static final int SCORE_MODE_DIF=0;
	public static final int SCORE_MODE_RMS=1;
	public static final int SCORE_MODE_AND=2;
	public static final int SCORE_MODE_MULT=3;
	public static final int SCORE_MODE_KS=4;
	
	public static int scoreMode1=SCORE_MODE_MULT;
	public static int scoreMode2=SCORE_MODE_RMS;
}
