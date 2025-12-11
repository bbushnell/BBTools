package cluster;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLongArray;

import dna.AminoAcid;
import jgi.Dedupe;
import shared.Tools;

/**
 * Utility methods for k-mer analysis and clustering operations in sequence processing.
 * Provides k-mer conversion, canonical encoding, frequency analysis, and similarity
 * scoring functions for bioinformatics clustering pipelines.
 *
 * @author Brian Bushnell
 * @date Mar 24, 2014
 */
public class ClusterTools {
	
	public static int[] toKmerCounts(byte[] bases, Object object, int k) {
		// TODO Auto-generated method stub
		return null;
	}

	public static int[] toKmers(final byte[] bases, int[] array_, final int k){
		if(bases==null || bases.length<k){return null;}
		final int alen=bases.length-k+1;
		final int[] array=(array_!=null && array_.length==alen ? array_ : new int[alen]);
		
		final int shift=2*k;
		final int shift2=shift-2;
		final int mask=~((-1)<<shift);
		
		int kmer=0;
		int rkmer=0;
		int len=0;
		
		for(int i=0, j=0; i<bases.length; i++){
			byte b=bases[i];
			int x=Dedupe.baseToNumber[b];
			int x2=Dedupe.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
//			if(b=='N'){len=0; rkmer=0;}else{len++;} //This version will transform 'N' into 'A'
			if(verbose){System.err.println("Scanning2 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k), Tools.min(i+1, k)));}
			if(len>=k){
				array[j]=Tools.min(kmer, rkmer);
				j++;
			}
		}
		
		Arrays.sort(array);
		return array;
	}
	
	public static int[] toKmerCounts(final byte[] bases, int[] array_, final int k, final int alen){
		if(bases==null || bases.length<k){return null;}
		final int[] array=(array_!=null && array_.length==alen ? array_ : new int[alen]);
		
		final int shift=2*k;
		final int shift2=shift-2;
		final int mask=~((-1)<<shift);
		
		int kmer=0;
		int rkmer=0;
		int len=0;
		
		for(int i=0, j=0; i<bases.length; i++){
			byte b=bases[i];
			int x=Dedupe.baseToNumber[b];
			int x2=Dedupe.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask;
//			if(b=='N'){len=0; rkmer=0;}else{len++;} //This version will transform 'N' into 'A'
			if(verbose){System.err.println("Scanning2 i="+i+", kmer="+kmer+", rkmer="+rkmer+", bases="+new String(bases, Tools.max(0, i-k), Tools.min(i+1, k)));}
			if(len>=k){
				array[Tools.min(kmer, rkmer)]++;
			}
		}
		
		Arrays.sort(array);
		return array;
	}
	
	public static int maxCanonicalKmer(int k){
		final int bits=2*k;
		final int max=(int)((1L<<bits)-1);
		int high=0;
		for(int kmer=0; kmer<=max; kmer++){
			int canon=Tools.min(kmer, AminoAcid.reverseComplementBinaryFast(kmer, k));
			high=Tools.max(canon, high);
		}
		return high;
	}
	
	/**
	 * Calculates fraction of k-mers present in cluster counts.
	 * Counts how many k-mers from the read have non-zero counts
	 * in the cluster's atomic count array.
	 *
	 * @param kmers Read k-mers array
	 * @param counts Cluster k-mer counts
	 * @return Fraction of k-mers present (0.0 to 1.0)
	 */
	static final float andCount(int[] kmers, AtomicLongArray counts){
		int sum=0;
		for(int i=0; i<kmers.length; i++){
			int kmer=kmers[i];
			long count=counts.get(kmer);
			if(count>0){sum++;}
		}
		return sum/(float)kmers.length;
	}
	
	/**
	 * Computes inner product between read k-mers and cluster probabilities.
	 * Sums the probability values for each k-mer present in the read,
	 * skipping negative k-mer values (invalid k-mers).
	 *
	 * @param kmers Read k-mers array
	 * @param probs Cluster k-mer probability distribution
	 * @return Sum of probabilities for matching k-mers
	 */
	static final float innerProduct(int[] kmers, float[] probs){
		float sum=0;
		for(int kmer : kmers){
			if(kmer>=0){
				sum+=probs[kmer];
			}
		}
		return sum;
	}
	
	/**
	 * Calculates sum of absolute differences between two frequency arrays.
	 * Compares corresponding elements and accumulates absolute differences
	 * to measure dissimilarity between k-mer frequency distributions.
	 *
	 * @param a First k-mer frequency array
	 * @param b Second k-mer frequency array
	 * @return Sum of absolute differences between arrays
	 */
	static final float absDif(float[] a, float[] b){
		assert(a.length==b.length);
		double sum=0;
		for(int i=0; i<a.length; i++){
			sum+=Tools.absdif((double)a[i], (double)b[i]);
		}

		return (float)sum;
	}
	
	/**
	 * Calculates root mean square difference between two frequency arrays.
	 * Computes squared differences, averages them, and returns the square root
	 * to measure the magnitude of differences between distributions.
	 *
	 * @param a First k-mer frequency array
	 * @param b Second k-mer frequency array
	 * @return Root mean square difference between arrays
	 */
	static final float rmsDif(float[] a, float[] b){
		assert(a.length==b.length);
		double sum=0;
		for(int i=0; i<a.length; i++){
//			double d=Tools.absdif((double)a[i], (double)b[i]);
			double d=((double)a[i])-((double)b[i]);
			sum+=d*d;
		}

		return (float)Math.sqrt(sum/a.length);
	}
	
	/**
	 * Computes Kullback-Leibler-like divergence between frequency distributions.
	 * Calculates sum of a[i] * log(a[i]/b[i]) with small epsilon added to prevent
	 * division by zero, measuring asymmetric difference between distributions.
	 *
	 * @param a First k-mer frequency array
	 * @param b Second k-mer frequency array
	 * @return KL-like divergence score
	 */
	static final float ksFunction(float[] a, float[] b){
		assert(a.length==b.length);
		double sum=0;
		for(int i=0; i<a.length; i++){
			double ai=a[i]+0.00001;
			double bi=b[i]+0.00001;
			double d=(double)ai*Math.log(ai/bi);
			sum+=d;
		}
		
		return (float)sum;
	}
	
	public static boolean verbose=false;
	
}
