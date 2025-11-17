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
	
	/**
	 * Creates a new cluster with specified parameters and initializes k-mer arrays.
	 *
	 * @param id_ Unique identifier for this cluster
	 * @param k1_ Length of 'big' k-mers for analysis
	 * @param k2_ Length of 'small' k-mers for analysis
	 * @param arraylen1_ Size of the big k-mer array
	 * @param arraylen2_ Size of the small k-mer array
	 */
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

	/**
	 * Recalculates GC content and k-mer probability distributions.
	 * Updates kmerProbArray1 and kmerProbArray2 based on accumulated k-mer counts
	 * with smoothing (95% from counts, 5% uniform distribution).
	 */
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
	
	/** Resets all atomic counters and k-mer arrays to zero.
	 * Used to clear cluster statistics before reprocessing. */
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
	
	/**
	 * Adds a read to this cluster, updating all statistical counters and k-mer arrays.
	 * Extracts k-mer counts, depth information, and GC content from the read's ReadTag.
	 * @param r The read to add to this cluster
	 */
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
	
	/**
	 * Calculates similarity score between a read and this cluster.
	 * Delegates to scoreSingle() or scorePaired() based on read type.
	 * @param r The read to score against this cluster
	 * @return Similarity score (higher values indicate better matches)
	 */
	public float score(Read r) {
		if(r==null){return 0;}
		return r.mate==null ? scoreSingle(r) : scorePaired(r);
	}
	
	/**
	 * Calculates similarity score for a single (unpaired) read.
	 * Combines depth, GC content, and k-mer similarity scores with weights.
	 * Note: Implementation is incomplete (contains TODO assertions).
	 *
	 * @param r The single read to score
	 * @return Weighted similarity score
	 */
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
	
	/**
	 * Calculates similarity score for paired reads.
	 * Note: Implementation is incomplete (marked as TODO).
	 * @param r Paired read to score
	 * @return Paired read similarity score
	 */
	public float scorePaired(Read r) {
		assert(false) : "TODO";
		if(r==null){return 0;}
		ReadTag rt=(ReadTag)r.obj;
		
//		ReadTag rt1=rt.r
		
		return 0;
	}
	
	/*--------------------------------------------------------------*/
	
	/** Unique identifier for this cluster */
	public final int id;
	
	/** 'big' kmer */
	public final int k1;
	/** 'small' kmer */
	public final int k2;

	/** Size of the big k-mer array (kmerArray1) */
	public final int arraylen1;
	/** Size of the small k-mer array (kmerArray2) */
	public final int arraylen2;
	
	/*--------------------------------------------------------------*/
	
	/** GC content ratio of reads in this cluster */
	public float gc;
	/** Coverage depth for strand 2 */
	/** Coverage depth for strand 1 */
	public int depth1, depth2;
	
	/** Thread-safe array storing counts for big k-mers (k1) */
	final AtomicLongArray kmerArray1;
	/** Probability distribution for big k-mers derived from counts */
	final float[] kmerProbArray1;
	
	/** Thread-safe array storing counts for small k-mers (k2) */
	final AtomicLongArray kmerArray2;
	/** Probability distribution for small k-mers derived from counts */
	final float[] kmerProbArray2;
	
	/** Thread-safe accumulator for total depth on strand 1 */
	final AtomicLong depthsum1=new AtomicLong(0);
	/** Thread-safe accumulator for total depth on strand 2 */
	final AtomicLong depthsum2=new AtomicLong(0);
	
	/** Thread-safe counter for total reads added to this cluster */
	final AtomicLong readCount=new AtomicLong(0);
	/** Thread-safe counter for total bases in all reads in this cluster */
	final AtomicLong baseCount=new AtomicLong(0);
//	final AtomicLong kmerCount=new AtomicLong(0);
	/** Thread-safe counter for total GC bases across all reads */
	final AtomicLong gcCount=new AtomicLong(0);
	
	/*--------------------------------------------------------------*/

	/** Scoring mode using absolute difference between probability distributions */
	public static final int SCORE_MODE_DIF=0;
	/** Scoring mode using root mean square difference between distributions */
	public static final int SCORE_MODE_RMS=1;
	/** Scoring mode using logical AND operation on k-mer presence */
	public static final int SCORE_MODE_AND=2;
	/** Scoring mode using multiplication (inner product) of distributions */
	public static final int SCORE_MODE_MULT=3;
	/**
	 * Scoring mode using Kolmogorov-Smirnov function for distribution comparison
	 */
	public static final int SCORE_MODE_KS=4;
	
	/** Active scoring mode for big k-mers (k1), defaults to SCORE_MODE_MULT */
	public static int scoreMode1=SCORE_MODE_MULT;
	/** Active scoring mode for small k-mers (k2), defaults to SCORE_MODE_RMS */
	public static int scoreMode2=SCORE_MODE_RMS;
}
