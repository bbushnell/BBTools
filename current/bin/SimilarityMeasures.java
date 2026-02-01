package bin;

import java.util.Arrays;

import parse.Parse;
import shared.Shared;
import shared.Tools;
import simd.Vector;

/** Mostly written by ChatGPT and modified by me */
public class SimilarityMeasures {
	
    public static void main(String[] args) {
        float[] sample1={0.1f, 0.2f, 0.3f, 0.4f};
        float[] sample2={0.1f, 0.2f, 0.4f, 0.3f};
        int[] sample1i={1, 2, 3, 4};
        int[] sample2i={1, 2, 4, 3};
        int[] sample3i={2, 4, 6, 8};
        int[] sample4i={8, 6, 4, 2};

        // Print the similarity vector
        System.out.println("Difference Vector Float12: "+Arrays.toString(calculateDifferenceVector(sample1, sample2)));
        System.out.println("Difference Vector Int12:   "+Arrays.toString(calculateDifferenceVector(sample1i, sample2i)));
        System.out.println("Difference Vector Int13:   "+Arrays.toString(calculateDifferenceVector(sample1i, sample3i)));
        System.out.println("Difference Vector Int14:   "+Arrays.toString(calculateDifferenceVector(sample1i, sample4i)));
    }
	
    /**
     * Parses command-line arguments to configure similarity measure flags.
     * Sets boolean flags for enabling different similarity measures.
     *
     * @param arg The full argument string (not used)
     * @param a The parameter name
     * @param b The parameter value to parse
     * @return true if the parameter was recognized and parsed
     */
    public static boolean parse(String arg, String a, String b){
    	if(a.equals("null")){
    		//Do nothing
    	}else if(a.equals("cosine") || a.equals("cos")){
    		COSINE=Parse.parseBoolean(b);
    	}else if(a.equals("gccompensated")){
    		GC_COMPENSATED=Parse.parseBoolean(b);
    	}else if(a.equals("euclid") || a.equals("euc")){
    		EUCLID=Parse.parseBoolean(b);
    	}else if(a.equals("absolute") || a.equals("abs")){
    		ABSOLUTE=Parse.parseBoolean(b);
    	}else if(a.equals("jsd")){
    		JSD=Parse.parseBoolean(b);
    	}else if(a.equals("hellinger") || a.equals("hell") || a.equals("hel")){
    		HELLINGER=Parse.parseBoolean(b);
    	}else if(a.equals("ks") || a.equals("kst")){
    		KST=Parse.parseBoolean(b);
    	}else {
    		return false;
    	}
    	
    	return true;
    }

    /**
     * Computes all similarity measures between two float arrays.
     * Returns a vector containing cosine difference, Euclidean distance,
     * absolute difference, Jensen-Shannon divergence, Hellinger distance,
     * and Kolmogorov-Smirnov test statistic.
     *
     * @param a First probability distribution
     * @param b Second probability distribution
     * @return Array of 6 similarity measure values
     */
    public static float[] calculateDifferenceVector(float[] a, float[] b) {
//        float cosineSimilarity=cosineSimilarity(a, b);
        float cosineDifference=cosineDifference(a, b);
        float euclideanDistance=euclideanDistance(a, b);
        float absoluteDifference=absDif(a, b);
        float jensenShannonDivergence=jensenShannonDivergence(a, b);
        float hellingerDistance=hellingerDistance(a, b);
        float ksDifference=ksTest(a, b);

        return new float[] {
            cosineDifference,
            euclideanDistance,
            absoluteDifference,
            jensenShannonDivergence,
            hellingerDistance,
            ksDifference
        };
    }

    //For setting thresholds before neural net is implemented
    /**
     * Calculates weighted average of enabled similarity measures for thresholding.
     * Only includes measures with enabled flags in the average.
     * Normalizes input arrays by their sums before calculation.
     *
     * @param a First frequency histogram
     * @param b Second frequency histogram
     * @return Average similarity score, or 0 if result is invalid
     */
    public static float calculateDifferenceAverage(int[] a, int[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
        float cosineDifference=(COSINE ? cosineDifference(a, b, inva, invb) : 0);
        float euclideanDistance=(EUCLID ? euclideanDistance(a, b, inva, invb) : 0);
        float absoluteDifference=(ABSOLUTE ? absDif(a, b, inva, invb) : 0);
        float jensenShannonDivergence=(JSD ? jensenShannonDivergence(a, b, inva, invb) : 0);
        float hellingerDistance=(HELLINGER? hellingerDistance(a, b, inva, invb) : 0);
        float ksDifference=(KST ? ksTest(a, b, inva, invb) : 0);
        int div=(COSINE ? 1 : 0)+(EUCLID ? 1 : 0)+(ABSOLUTE ? 1 : 0)+(JSD ? 1 : 0)+(HELLINGER ? 1 : 0)+(KST ? 1 : 0);
        float ret=(cosineDifference+euclideanDistance+absoluteDifference+
        		jensenShannonDivergence+hellingerDistance+ksDifference)/div;
        return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

    /**
     * Computes all similarity measures between two integer arrays.
     * Normalizes arrays by their sums and returns cosine difference,
     * Euclidean distance, absolute difference, Jensen-Shannon divergence,
     * Hellinger distance, and Kolmogorov-Smirnov test statistic.
     *
     * @param a First frequency histogram
     * @param b Second frequency histogram
     * @return Array of 6 similarity measure values
     */
    public static float[] calculateDifferenceVector(int[] a, int[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
//        float cosineSimilarity=cosineSimilarity(a, b, inva, invb);
        float cosineDifference=cosineDifference(a, b, inva, invb);
        float euclideanDistance=euclideanDistance(a, b, inva, invb);
        float absoluteDifference=absDif(a, b, inva, invb);
        float jensenShannonDivergence=jensenShannonDivergence(a, b, inva, invb);
        float hellingerDistance=hellingerDistance(a, b, inva, invb);
        float ksDifference=ksTest(a, b, inva, invb);

        return new float[] {
                cosineDifference,
            euclideanDistance,
            absoluteDifference,
            jensenShannonDivergence,
            hellingerDistance,
            ksDifference
        };
    }

    /**
     * Computes cosine difference (1 - cosine similarity) for float arrays.
     * @param a First vector
     * @param b Second vector
     * @return Cosine difference value between 0 and 2
     */
    public static float cosineDifference(float[] a, float[] b) {
    	return 1-cosineSimilarity(a, b);
    }
    
    /**
     * Computes cosine similarity between two float arrays.
     * Calculates dot product divided by product of vector magnitudes.
     * Returns 0 for invalid or negative results.
     *
     * @param a First vector
     * @param b Second vector
     * @return Cosine similarity value between 0 and 1
     */
    public static float cosineSimilarity(float[] a, float[] b) {
        float dotProduct=0f;
        float normVec1=0f;
        float normVec2=0f;

        for (int i=0; i<a.length; i++) {
        	float ai=a[i], bi=b[i];
            dotProduct+=ai*bi;
            normVec1+=ai*ai;
            normVec2+=bi*bi;
        }

        float ret=(float)(dotProduct/(Math.sqrt(normVec1)*Math.sqrt(normVec2)));
        return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

    /**
     * Computes cosine difference for integer arrays with normalization.
     * Normalizes arrays by their sums before calculating similarity.
     *
     * @param a First frequency array
     * @param b Second frequency array
     * @return Cosine difference, or 0 if result is invalid
     */
    public static float cosineDifference(int[] a, int[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=1-cosineSimilarity(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

    /**
     * Computes cosine difference with precomputed normalization factors.
     *
     * @param a First frequency array
     * @param b Second frequency array
     * @param inva Inverse of sum of array a (1/sum(a))
     * @param invb Inverse of sum of array b (1/sum(b))
     * @return Cosine difference value
     */
    public static float cosineDifference(int[] a, int[] b, float inva, float invb) {
    	return 1-cosineSimilarity(a, b, inva, invb);
    }

    /**
     * Computes GC-compensated cosine difference using k-mer GC content mapping.
     * Applies GC bias correction before calculating cosine similarity.
     *
     * @param a First k-mer frequency array
     * @param b Second k-mer frequency array
     * @param k K-mer length for GC mapping
     * @return GC-compensated cosine difference
     */
    public static float cosineDifferenceCompensated(int[] a, int[] b, int k) {
    	return 1-cosineSimilarityCompensated(a, b, k, BinObject.gcmapMatrix[k]);
    }
    
    /**
     * Computes cosine similarity for integer arrays with normalization factors.
     * Uses SIMD acceleration when available or GC compensation if enabled.
     * Applies numerical stability by limiting minimum norm values.
     *
     * @param a First frequency array
     * @param b Second frequency array
     * @param inva Inverse of sum of array a
     * @param invb Inverse of sum of array b
     * @return Cosine similarity between 0 and 1
     */
    public static float cosineSimilarity(int[] a, int[] b, float inva, float invb) {
    	if(GC_COMPENSATED) {return cosineSimilarityCompensated(a, b, 4);}
    	if(Shared.SIMD) {return Vector.cosineSimilarity(a, b, inva, invb);}
        float dotProduct=0f;
        float normVec1=0f;
        float normVec2=0f;

        for (int i=0; i<a.length; i++) {
        	float ai=a[i]*inva, bi=b[i]*invb;
            dotProduct+=ai*bi;
            normVec1+=ai*ai;
            normVec2+=bi*bi;
        }

	    normVec1=Math.max(normVec1, 1e-15f);
	    normVec2=Math.max(normVec2, 1e-15f);
        return (float)(dotProduct/(Math.sqrt(normVec1)*Math.sqrt(normVec2)));
    }

    /**
     * Computes cosine difference for long arrays with normalization.
     * @param a First frequency array
     * @param b Second frequency array
     * @return Cosine difference, or 0 if result is invalid
     */
    public static float cosineDifference(long[] a, long[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=1-cosineSimilarity(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

    /**
     * Computes cosine difference for long arrays with precomputed normalization.
     *
     * @param a First frequency array
     * @param b Second frequency array
     * @param inva Inverse of sum of array a
     * @param invb Inverse of sum of array b
     * @return Cosine difference value
     */
    public static float cosineDifference(long[] a, long[] b, float inva, float invb) {
    	return 1-cosineSimilarity(a, b, inva, invb);
    }
    
    /**
     * Computes cosine similarity for long arrays with normalization factors.
     * Applies numerical stability by limiting minimum norm values.
     *
     * @param a First frequency array
     * @param b Second frequency array
     * @param inva Inverse of sum of array a
     * @param invb Inverse of sum of array b
     * @return Cosine similarity between 0 and 1
     */
    public static float cosineSimilarity(long[] a, long[] b, float inva, float invb) {
        float dotProduct=0f;
        float normVec1=0f;
        float normVec2=0f;

        for (int i=0; i<a.length; i++) {
        	float ai=a[i]*inva, bi=b[i]*invb;
            dotProduct+=ai*bi;
            normVec1+=ai*ai;
            normVec2+=bi*bi;
        }

	    normVec1=Math.max(normVec1, 1e-15f);
	    normVec2=Math.max(normVec2, 1e-15f);
        return (float)(dotProduct/(Math.sqrt(normVec1)*Math.sqrt(normVec2)));
    }

    /**
     * Computes GC-compensated cosine similarity using default GC mapping.
     *
     * @param a First k-mer frequency array
     * @param b Second k-mer frequency array
     * @param k K-mer length for GC content mapping
     * @return GC-compensated cosine similarity
     */
    public static float cosineSimilarityCompensated(int[] a, int[] b, int k) {
    	return cosineSimilarityCompensated(a, b, k, BinObject.gcmapMatrix[k]);
    }
    
    /**
     * Applies GC bias compensation to frequency array.
     * Groups k-mers by GC content and normalizes within groups to reduce bias.
     *
     * @param a Frequency array to compensate
     * @param k K-mer length
     * @param gcmap GC content mapping for each k-mer position
     * @return GC-compensated frequency array
     */
    public static float[] compensate(int[] a, int k, int[] gcmap) {
    	float[] aSum=new float[k+1];
    	
    	for(int i=0; i<a.length; i++) {
    		int gc=gcmap[i];
    		aSum[gc]+=a[i];
    	}
    	
    	for(int i=0; i<aSum.length; i++) {
    		aSum[i]=1f/Math.max(aSum[i], 1);
    	}
    	assert(Tools.sum(aSum)==1);

    	float[] comp=new float[a.length];
    	for(int i=0; i<a.length; i++) {
        	int gc=gcmap[i];
    		comp[i]=a[i]*aSum[gc];
    	}
    	return comp;
    }
    
    /**
     * Applies GC bias compensation to long frequency array using vectorized operations.
     * @param a Frequency array to compensate
     * @param k K-mer length for GC content mapping
     * @return GC-compensated frequency array
     */
    public static float[] compensate(long[] a, int k) {
    	final int[] gcmap=BinObject.gcmapMatrix[k];
    	return Vector.compensate(a, k, gcmap);
    }
    
    /**
     * Computes GC-compensated cosine similarity with explicit GC mapping.
     * Groups k-mers by GC content, normalizes within groups, then calculates similarity.
     * Returns 0 for invalid or negative results.
     *
     * @param a First k-mer frequency array
     * @param b Second k-mer frequency array
     * @param k K-mer length
     * @param gcmap GC content for each k-mer position
     * @return GC-compensated cosine similarity
     */
    public static float cosineSimilarityCompensated(int[] a, int[] b, int k, int[] gcmap) {
    	
    	float[] aSum=new float[k+1];
    	float[] bSum=new float[k+1];
    	
    	for(int i=0; i<a.length; i++) {
    		int gc=gcmap[i];
    		aSum[gc]+=a[i];
    		bSum[gc]+=b[i];
    	}
    	
    	for(int i=0; i<aSum.length; i++) {
    		aSum[i]=1f/Math.max(aSum[i], 1);
    		bSum[i]=1f/Math.max(bSum[i], 1);
    	}
    	
        float dotProduct=0f;
        float normVec1=0f;
        float normVec2=0f;

        for (int i=0; i<a.length; i++) {
        	int gc=gcmap[i];
        	float ai=a[i]*aSum[gc], bi=b[i]*bSum[gc];
            dotProduct+=ai*bi;
            normVec1+=ai*ai;
            normVec2+=bi*bi;
        }

        float ret=(float)(dotProduct/(Math.sqrt(normVec1)*Math.sqrt(normVec2)));
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

    /**
     * Computes Euclidean distance between two float arrays.
     * Calculates square root of sum of squared differences.
     *
     * @param a First vector
     * @param b Second vector
     * @return Euclidean distance
     */
    public static float euclideanDistance(float[] a, float[] b) {
        float sumSquaredDifferences=0f;

        for (int i=0; i<a.length; i++) {
        	float ai=a[i], bi=b[i];
        	float d=ai-bi;
            sumSquaredDifferences+=d*d;
        }

        return (float)Math.sqrt(sumSquaredDifferences);
    }
    

    /**
     * Computes Euclidean distance for normalized integer arrays.
     * @param a First frequency array
     * @param b Second frequency array
     * @return Euclidean distance, or 0 if result is invalid
     */
    public static float euclideanDistance(int[] a, int[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=euclideanDistance(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

    /**
     * Computes Euclidean distance with precomputed normalization factors.
     *
     * @param a First frequency array
     * @param b Second frequency array
     * @param inva Inverse of sum of array a
     * @param invb Inverse of sum of array b
     * @return Euclidean distance
     */
    public static float euclideanDistance(int[] a, int[] b, float inva, float invb) {
        float sumSquaredDifferences=0f;

        for (int i=0; i<a.length; i++) {
        	float ai=a[i]*inva, bi=b[i]*invb;
        	float d=ai-bi;
            sumSquaredDifferences+=d*d;
        }

        return (float)Math.sqrt(sumSquaredDifferences);
    }
    

    /**
     * Computes Euclidean distance for normalized long arrays.
     * @param a First frequency array
     * @param b Second frequency array
     * @return Euclidean distance, or 0 if result is invalid
     */
    public static float euclideanDistance(long[] a, long[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=euclideanDistance(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

    /**
     * Computes Euclidean distance for long arrays with precomputed normalization.
     *
     * @param a First frequency array
     * @param b Second frequency array
     * @param inva Inverse of sum of array a
     * @param invb Inverse of sum of array b
     * @return Euclidean distance
     */
    public static float euclideanDistance(long[] a, long[] b, float inva, float invb) {
        float sumSquaredDifferences=0f;

        for (int i=0; i<a.length; i++) {
        	float ai=a[i]*inva, bi=b[i]*invb;
        	float d=ai-bi;
            sumSquaredDifferences+=d*d;
        }

        return (float)Math.sqrt(sumSquaredDifferences);
    }
	
	/**
	 * @param a Contig kmer frequencies
	 * @param b Cluster kmer frequencies
	 * @return Score
	 */
	static final float absDif(float[] a, float[] b){
		assert(a.length==b.length);
		double sum=0;
		for(int i=0; i<a.length; i++){
			sum+=Math.abs(a[i]-b[i]);
		}

		return (float)sum;
	}
	
	/**
	 * @param a Contig kmer frequencies
	 * @param b Cluster kmer frequencies
	 * @return Score
	 */
	static final float absDifFloat(float[] a, float[] b){
    	if(Shared.SIMD) {return Vector.absDifFloat(a, b);}
		assert(a.length==b.length);
		float sum=0;
		for(int i=0; i<a.length; i++){
			sum+=Math.abs(a[i]-b[i]);
		}
		return (float)sum;
	}
    
    /**
     * Computes absolute difference for normalized integer arrays.
     * @param a First frequency array
     * @param b Second frequency array
     * @return Sum of absolute differences, or 0 if result is invalid
     */
    public static float absDif(int[] a, int[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=absDif(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }
    
	/**
	 * Computes sum of absolute differences with precomputed normalization.
	 *
	 * @param a First frequency array
	 * @param b Second frequency array
	 * @param inva Inverse of sum of array a
	 * @param invb Inverse of sum of array b
	 * @return Sum of absolute differences
	 */
	static final float absDif(int[] a, int[] b, float inva, float invb){
		assert(a.length==b.length);
		float sum=0;
		for(int i=0; i<a.length; i++){
			float ai=a[i]*inva, bi=b[i]*invb;
			sum+=Math.abs(ai-bi);
		}
		return sum;
	}
    
    /**
     * Computes GC-compensated absolute difference using vectorized operations.
     * Applies GC bias correction before calculating absolute differences.
     *
     * @param a First k-mer frequency array
     * @param b Second k-mer frequency array
     * @param k K-mer length for GC compensation
     * @return GC-compensated absolute difference, clamped to [0,1]
     */
    public static float absDifComp(long[] a, long[] b, int k) {
    	float[] af=compensate(a, k);
    	float[] bf=compensate(b, k);
    	float ret=Vector.absDifFloat(af, bf);
    	return Tools.mid(0, 1, (Float.isFinite(ret) && ret>0 ? ret : 0));
    }
    
    /**
     * Computes absolute difference for normalized long arrays.
     * @param a First frequency array
     * @param b Second frequency array
     * @return Sum of absolute differences, or 0 if result is invalid
     */
    public static float absDif(long[] a, long[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=absDif(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }
    
	/**
	 * Computes sum of absolute differences for long arrays with normalization.
	 *
	 * @param a First frequency array
	 * @param b Second frequency array
	 * @param inva Inverse of sum of array a
	 * @param invb Inverse of sum of array b
	 * @return Sum of absolute differences
	 */
	private static final float absDif(long[] a, long[] b, float inva, float invb){
		assert(a.length==b.length);
		float sum=0;
		for(int i=0; i<a.length; i++){
			float ai=a[i]*inva, bi=b[i]*invb;
			sum+=Math.abs(ai-bi);
		}
		return sum;
	}

    /**
     * Computes Jensen-Shannon divergence between two probability distributions.
     * Adds small epsilon (0.0005) to prevent log(0) errors.
     *
     * @param a First probability distribution
     * @param b Second probability distribution
     * @return Jensen-Shannon divergence value
     */
    public static float jensenShannonDivergence(float[] a, float[] b) {
        float kldSumA=0, kldSumB=0;
        for (int i=0; i<a.length; i++) {
        	float ai=a[i]+0.0005f, bi=b[i]+0.0005f;//Prevents zero values
        	float avgi=(ai+bi)*0.5f;
            kldSumA+=ai*Math.log(ai/avgi);
            kldSumA+=bi*Math.log(bi/avgi);
        }
        return (kldSumA+kldSumB)*invLog2*0.5f;
    }
    

    /**
     * Computes Jensen-Shannon divergence for normalized integer arrays.
     * @param a First frequency array
     * @param b Second frequency array
     * @return Jensen-Shannon divergence, or 0 if result is invalid
     */
    public static float jensenShannonDivergence(int[] a, int[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=jensenShannonDivergence(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

    /**
     * Computes Jensen-Shannon divergence with precomputed normalization factors.
     * Adds epsilon to normalized values to prevent log(0) errors.
     *
     * @param a First frequency array
     * @param b Second frequency array
     * @param inva Inverse of sum of array a
     * @param invb Inverse of sum of array b
     * @return Jensen-Shannon divergence value
     */
    public static float jensenShannonDivergence(int[] a, int[] b, float inva, float invb) {
        float kldSumA=0, kldSumB=0;
        for (int i=0; i<a.length; i++) {
        	float ai=a[i]*inva+0.0005f, bi=b[i]*invb+0.0005f;//Prevents zero values
        	float avgi=(ai+bi)*0.5f;
            kldSumA+=ai*Math.log(ai/avgi);
            kldSumA+=bi*Math.log(bi/avgi);
        }
        return (kldSumA+kldSumB)*invLog2*0.5f;
    }

//    public static float jensenShannonDivergenceSlow(float[] a, float[] b) {
//        float[] avg=new float[a.length];
//        for (int i=0; i<a.length; i++) {
//        	float ai=a[i], bi=b[i];
//            avg[i]=(ai+bi)*0.5f;
//        }
//
//        return (kullbackLeiblerDivergence(a, avg)+kullbackLeiblerDivergence(b, avg))*0.5f;
//    }
//
//    public static float kullbackLeiblerDivergence(float[] p, float[] q) {
//        float sum=0f;
//        for (int i=0; i<p.length; i++) {
//        	float pi=p[i], qi=q[i];
//            if (p[i]!=0) {
//                sum+=p[i]*Math.log(pi/qi);
//            }
//        }
//        return sum*invLog2;
//    }

    /**
     * Computes Hellinger distance between two probability distributions.
     * Calculates sqrt(sum((sqrt(ai) - sqrt(bi))^2)) / sqrt(2).
     *
     * @param a First probability distribution
     * @param b Second probability distribution
     * @return Hellinger distance value
     */
    public static float hellingerDistance(float[] a, float[] b) {
        float sum=0f;
        for (int i=0; i<a.length; i++) {
        	float ai=a[i], bi=b[i];
        	float d=(float)(Math.sqrt(ai)-Math.sqrt(bi));
            sum+=d*d;
        }
        return (float)Math.sqrt(sum)*invRoot2;
    }
    
    /**
     * Computes Hellinger distance for normalized integer arrays.
     * @param a First frequency array
     * @param b Second frequency array
     * @return Hellinger distance, or 0 if result is invalid
     */
    public static float hellingerDistance(int[] a, int[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=hellingerDistance(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

    /**
     * Computes Hellinger distance with precomputed normalization factors.
     *
     * @param a First frequency array
     * @param b Second frequency array
     * @param inva Inverse of sum of array a
     * @param invb Inverse of sum of array b
     * @return Hellinger distance value
     */
    public static float hellingerDistance(int[] a, int[] b, float inva, float invb) {
        float sum=0f;
        for (int i=0; i<a.length; i++) {
        	float ai=a[i]*inva, bi=b[i]*invb;
        	float d=(float)(Math.sqrt(ai)-Math.sqrt(bi));
            sum+=d*d;
        }
        return (float)Math.sqrt(sum)*invRoot2;
    }
    
    /**
     * Computes Hellinger distance for normalized long arrays.
     * @param a First frequency array
     * @param b Second frequency array
     * @return Hellinger distance, or 0 if result is invalid
     */
    public static float hellingerDistance(long[] a, long[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=hellingerDistance(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

    /**
     * Computes Hellinger distance for long arrays with precomputed normalization.
     *
     * @param a First frequency array
     * @param b Second frequency array
     * @param inva Inverse of sum of array a
     * @param invb Inverse of sum of array b
     * @return Hellinger distance value
     */
    public static float hellingerDistance(long[] a, long[] b, float inva, float invb) {
        float sum=0f;
        for (int i=0; i<a.length; i++) {
        	float ai=a[i]*inva, bi=b[i]*invb;
        	float d=(float)(Math.sqrt(ai)-Math.sqrt(bi));
            sum+=d*d;
        }
        return (float)Math.sqrt(sum)*invRoot2;
    }
    
    /** This is a KS test for binned histograms, not raw values */
    public static float ksTest(float[] histogram1, float[] histogram2) {
        // Ensure both histograms have the same length
        if (histogram1.length!=histogram2.length) {
            throw new IllegalArgumentException("Histograms must have the same number of bins");
        }

        float cd1=0, cd2=0, dMax=0;

        // Compute the KS statistic (maximum absolute difference between the two CDFs)
        for (int i=0; i<histogram1.length; i++) {
        	cd1+=histogram1[i];
        	cd2+=histogram2[i];
            dMax=(float)Math.max(dMax, Math.abs(cd1-cd2));
        }

        return dMax;
    }
    
    /** This is a KS test for binned histograms, not raw values */
    public static float ksTest(int[] a, int[] b, float inva, float invb) {
        // Ensure both histograms have the same length
        if (a.length!=b.length) {
            throw new IllegalArgumentException("Histograms must have the same number of bins");
        }

        float cda=0, cdb=0, dMax=0;

        // Compute the KS statistic (maximum absolute difference between the two CDFs)
        for (int i=0; i<a.length; i++) {
        	float ai=a[i]*inva, bi=b[i]*invb;
        	cda+=ai;
        	cdb+=bi;
            dMax=(float)Math.max(dMax, Math.abs(cda-cdb));
        }

        return dMax;
    }

    /** Square root of 2 constant for Hellinger distance normalization */
    private static final float root2=(float)Math.sqrt(2);
    /** Natural logarithm of 2 for Jensen-Shannon divergence conversion */
    private static final float log2=(float)Math.log(2);
    /** Inverse of square root of 2 for efficient Hellinger distance calculation */
    private static final float invRoot2=1/root2;
    /** Inverse of natural logarithm of 2 for Jensen-Shannon divergence */
    private static final float invLog2=1/log2;


    /** Enable GC bias compensation in similarity calculations */
    public static boolean GC_COMPENSATED=false;
    
    //2531 kcps (times include contig loading)
    //26 clusters
//    Completeness Score:             60.278
//    Contamination Score:            2.1925
    /** Enable cosine similarity calculation in composite measures */
    public static boolean COSINE=true;
    //2796 kcps
    //21 clusters
    //Completeness Score:             60.909
    //Contamination Score:            2.3108
    /** Enable Euclidean distance calculation in composite measures */
    public static boolean EUCLID=false;//0.008
    //2636 kcps
    //23 clusters at 4x threshold of cosine
//    Completeness Score:             60.947
//    Contamination Score:            1.7679
    /** Enable absolute difference calculation in composite measures */
    public static boolean ABSOLUTE=false; //Best at 0.089
    //183 kcps
    //22 clusters
//  Completeness Score:             59.169
//  Contamination Score:            2.1959
    /** Enable Jensen-Shannon divergence calculation in composite measures */
    public static boolean JSD=false;
    //953 kcps
    //~22 at 2x threshold of cosine
//    Completeness Score:             61.072
//    Contamination Score:            1.9358
    /** Enable Hellinger distance calculation in composite measures */
    public static boolean HELLINGER=false;//0.0425
    //1859 kcps
    //20 clusters
//    Completeness Score:             26.380
//    Contamination Score:            3.1930
    /** Enable Kolmogorov-Smirnov test calculation in composite measures */
    public static boolean KST=false;
    
}
