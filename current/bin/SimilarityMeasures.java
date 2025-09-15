package bin;

import java.util.Arrays;

import shared.Parse;
import shared.Shared;
import shared.Tools;
import shared.Vector;

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

    public static float cosineDifference(float[] a, float[] b) {
    	return 1-cosineSimilarity(a, b);
    }
    
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

    public static float cosineDifference(int[] a, int[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=1-cosineSimilarity(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

    public static float cosineDifference(int[] a, int[] b, float inva, float invb) {
    	return 1-cosineSimilarity(a, b, inva, invb);
    }

    public static float cosineDifferenceCompensated(int[] a, int[] b, int k) {
    	return 1-cosineSimilarityCompensated(a, b, k, BinObject.gcmapMatrix[k]);
    }
    
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

    public static float cosineDifference(long[] a, long[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=1-cosineSimilarity(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

    public static float cosineDifference(long[] a, long[] b, float inva, float invb) {
    	return 1-cosineSimilarity(a, b, inva, invb);
    }
    
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

    public static float cosineSimilarityCompensated(int[] a, int[] b, int k) {
    	return cosineSimilarityCompensated(a, b, k, BinObject.gcmapMatrix[k]);
    }
    
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
    
    public static float[] compensate(long[] a, int k) {
    	final int[] gcmap=BinObject.gcmapMatrix[k];
    	return Vector.compensate(a, k, gcmap);
    }
    
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

    public static float euclideanDistance(float[] a, float[] b) {
        float sumSquaredDifferences=0f;

        for (int i=0; i<a.length; i++) {
        	float ai=a[i], bi=b[i];
        	float d=ai-bi;
            sumSquaredDifferences+=d*d;
        }

        return (float)Math.sqrt(sumSquaredDifferences);
    }
    

    public static float euclideanDistance(int[] a, int[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=euclideanDistance(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

    public static float euclideanDistance(int[] a, int[] b, float inva, float invb) {
        float sumSquaredDifferences=0f;

        for (int i=0; i<a.length; i++) {
        	float ai=a[i]*inva, bi=b[i]*invb;
        	float d=ai-bi;
            sumSquaredDifferences+=d*d;
        }

        return (float)Math.sqrt(sumSquaredDifferences);
    }
    

    public static float euclideanDistance(long[] a, long[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=euclideanDistance(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

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
    
    public static float absDif(int[] a, int[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=absDif(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }
    
	static final float absDif(int[] a, int[] b, float inva, float invb){
		assert(a.length==b.length);
		float sum=0;
		for(int i=0; i<a.length; i++){
			float ai=a[i]*inva, bi=b[i]*invb;
			sum+=Math.abs(ai-bi);
		}
		return sum;
	}
    
    public static float absDifComp(long[] a, long[] b, int k) {
    	float[] af=compensate(a, k);
    	float[] bf=compensate(b, k);
    	float ret=Vector.absDifFloat(af, bf);
    	return Tools.mid(0, 1, (Float.isFinite(ret) && ret>0 ? ret : 0));
    }
    
    public static float absDif(long[] a, long[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=absDif(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }
    
	private static final float absDif(long[] a, long[] b, float inva, float invb){
		assert(a.length==b.length);
		float sum=0;
		for(int i=0; i<a.length; i++){
			float ai=a[i]*inva, bi=b[i]*invb;
			sum+=Math.abs(ai-bi);
		}
		return sum;
	}

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
    

    public static float jensenShannonDivergence(int[] a, int[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=jensenShannonDivergence(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

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

    public static float hellingerDistance(float[] a, float[] b) {
        float sum=0f;
        for (int i=0; i<a.length; i++) {
        	float ai=a[i], bi=b[i];
        	float d=(float)(Math.sqrt(ai)-Math.sqrt(bi));
            sum+=d*d;
        }
        return (float)Math.sqrt(sum)*invRoot2;
    }
    
    public static float hellingerDistance(int[] a, int[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=hellingerDistance(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

    public static float hellingerDistance(int[] a, int[] b, float inva, float invb) {
        float sum=0f;
        for (int i=0; i<a.length; i++) {
        	float ai=a[i]*inva, bi=b[i]*invb;
        	float d=(float)(Math.sqrt(ai)-Math.sqrt(bi));
            sum+=d*d;
        }
        return (float)Math.sqrt(sum)*invRoot2;
    }
    
    public static float hellingerDistance(long[] a, long[] b) {
    	float inva=1f/Math.max(1, Tools.sum(a));
    	float invb=1f/Math.max(1, Tools.sum(b));
    	float ret=hellingerDistance(a, b, inva, invb);
    	return (Float.isFinite(ret) && ret>0 ? ret : 0);
    }

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

    private static final float root2=(float)Math.sqrt(2);
    private static final float log2=(float)Math.log(2);
    private static final float invRoot2=1/root2;
    private static final float invLog2=1/log2;


    public static boolean GC_COMPENSATED=false;
    
    //2531 kcps (times include contig loading)
    //26 clusters
//    Completeness Score:             60.278
//    Contamination Score:            2.1925
    public static boolean COSINE=true;
    //2796 kcps
    //21 clusters
    //Completeness Score:             60.909
    //Contamination Score:            2.3108
    public static boolean EUCLID=false;//0.008
    //2636 kcps
    //23 clusters at 4x threshold of cosine
//    Completeness Score:             60.947
//    Contamination Score:            1.7679
    public static boolean ABSOLUTE=false; //Best at 0.089
    //183 kcps
    //22 clusters
//  Completeness Score:             59.169
//  Contamination Score:            2.1959
    public static boolean JSD=false;
    //953 kcps
    //~22 at 2x threshold of cosine
//    Completeness Score:             61.072
//    Contamination Score:            1.9358
    public static boolean HELLINGER=false;//0.0425
    //1859 kcps
    //20 clusters
//    Completeness Score:             26.380
//    Contamination Score:            3.1930
    public static boolean KST=false;
    
}
