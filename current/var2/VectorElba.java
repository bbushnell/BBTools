package var2;

import shared.Tools;

/**
 * Feature vector generator for Elba's quality score prediction model.
 * Creates normalized feature vectors from variant data for machine learning classification.
 * Transforms raw variant metrics into standardized features optimized for quality prediction.
 * @author Brian Bushnell
 */
public class VectorElba {

	/**
	 * Elba's feature vector for quality score prediction.
	 * Focus on quality metrics and statistical measures.
	 * 
	 * @param v Variant to convert
	 * @param pairingRate Overall proper pairing rate
	 * @param totalQualityAvg Average base quality from dataset
	 * @param totalMapqAvg Average mapping quality from dataset 
	 * @param readLengthAvg Average read length 
	 * @param ploidy Sample ploidy 
	 * @param map Scaffold mapping
	 * @return Feature vector for Elba's quality prediction model
	 */
	public static float[] makeElbaVector(Var v, double pairingRate, double totalQualityAvg,
			double totalMapqAvg, double readLengthAvg, int ploidy, ScafMap map){

		ElbaMinMaxValues minMaxVals = new ElbaMinMaxValues();

		float[] vec = new float[11]; 
		int idx = 0;

		final int alleleCount = v.alleleCount();
		final int safeCount = Tools.max(1, alleleCount);


		// raw values
		float raw_sb = (float)v.strandRatio();
		float raw_avg_read_length = (alleleCount > 0 ? (float)(v.lengthSum / (double)safeCount) : 0f);
		float raw_avg_identity = (alleleCount > 0 ? (float)(v.idSum / (double)safeCount) : 0f);
		float raw_avg_mapq = (alleleCount > 0 ? (float)(v.mapQSum / (double)safeCount) : 0f);
		float raw_avg_baseq = (alleleCount > 0 ? (float)(v.baseQSum / (double)safeCount) : 0f);
		float raw_avg_enddist = (alleleCount > 0 ? (float)(v.endDistSum / (double)safeCount) : 0f);
		float raw_avg_depth_per_allele = (float)(v.coverage() / (double)Tools.max(1, alleleCount));


		float transformed_sb = elbaTransformLeftSkewed(raw_sb);
		vec[idx++] = elbaNormalize(transformed_sb, minMaxVals.min_transformed_sb, minMaxVals.max_transformed_sb);

		float transformed_avg_read_length = elbaTransformLeftSkewed(raw_avg_read_length);
		vec[idx++] = elbaNormalize(transformed_avg_read_length, minMaxVals.min_transformed_avg_read_length, minMaxVals.max_transformed_avg_read_length);

		float transformed_avg_identity = elbaTransformLeftSkewed(raw_avg_identity);
		vec[idx++] = elbaNormalize(transformed_avg_identity, minMaxVals.min_transformed_avg_identity, minMaxVals.max_transformed_avg_identity);

		vec[idx++] = elbaNormalize(raw_avg_mapq, minMaxVals.min_avg_mapq, minMaxVals.max_avg_mapq);

		vec[idx++] = (float)(1.0 - v.alleleFraction());

		vec[idx++] = (float)v.alleleFraction(); //AF

		vec[idx++] = elbaNormalize(raw_avg_baseq, minMaxVals.min_avg_baseq, minMaxVals.max_avg_baseq); //Avg_BaseQ

		vec[idx++] = elbaNormalize(raw_avg_enddist, minMaxVals.min_avg_enddist, minMaxVals.max_avg_enddist); //Avg_EndDist

		float transformed_avg_depth_per_allele = elbaTransformRightSkewed(raw_avg_depth_per_allele); //Avg_Depth_Per_Allele
		vec[idx++] = elbaNormalize(transformed_avg_depth_per_allele, minMaxVals.min_transformed_avg_depth_per_allele, minMaxVals.max_transformed_avg_depth_per_allele);

		vec[idx++] = (float)elbaEncodeTyp(v.type()); //TYP

		int len = elbaCalculateLen(v); //LEN
		vec[idx++] = (float)len;

		return vec;
	}

	/** Minimum and maximum values for feature normalization.
	 * Contains predefined bounds for scaling variant metrics to [0,1] range. */
	private static class ElbaMinMaxValues {
		/** Maximum average mapping quality for normalization */
		/** Minimum average mapping quality for normalization */
		public float min_avg_mapq = 0.0f, max_avg_mapq = 45.0f;    
		/** Maximum average base quality for normalization */
		/** Minimum average base quality for normalization */
		public float min_avg_baseq = 0.0f, max_avg_baseq = 45.0f;       
		/** Maximum average end distance for normalization */
		/** Minimum average end distance for normalization */
		public float min_avg_enddist = 0.0f, max_avg_enddist = 75.0f;  
		/** Maximum transformed average depth per allele for normalization */
		/** Minimum transformed average depth per allele for normalization */
		public float min_transformed_avg_depth_per_allele = 0.0f, max_transformed_avg_depth_per_allele = 3.1f;  
		/** Maximum transformed average read length for normalization */
		/** Minimum transformed average read length for normalization */
		public float min_transformed_avg_read_length = 0.0f, max_transformed_avg_read_length = 2.5e21f;       
		/** Maximum transformed average identity for normalization */
		/** Minimum transformed average identity for normalization */
		public float min_transformed_avg_identity = 0.0f, max_transformed_avg_identity = 9.5e29f;    
		/** Maximum transformed strand bias for normalization */
		/** Minimum transformed strand bias for normalization */
		public float min_transformed_sb = 0.0f, max_transformed_sb = 1.0f;               

	}


	/**
	 * Applies log transformation for right-skewed distributions.
	 * Uses log(value + 1) to handle zero values and compress high outliers.
	 * @param value Raw value to transform
	 * @return Log-transformed value
	 */
	private static float elbaTransformRightSkewed(float value) {
		if (value <= 0) return 0.0f;
		return (float)Math.log(value + 1);
	}


	/**
	 * Applies power transformation for left-skewed distributions.
	 * Uses value^10 to expand small differences and compress large ones.
	 * @param value Raw value to transform
	 * @return Power-transformed value
	 */
	private static float elbaTransformLeftSkewed(float value) {
		return (float)Math.pow(value, 10);
	}


	/**
	 * Normalizes value to [0,1] range using min-max scaling.
	 * Returns 0 if min equals max to avoid division by zero.
	 *
	 * @param value Value to normalize
	 * @param min_val Minimum value for scaling
	 * @param max_val Maximum value for scaling
	 * @return Normalized value in [0,1] range
	 */
	private static float elbaNormalize(float value, float min_val, float max_val) {
		if (max_val == min_val) return 0.0f;
		return (value - min_val) / (max_val - min_val);
	}


	/**
	 * Scales quality scores using piecewise linear transformation.
	 * Maps 0-20 to 0-0.25 and 20-100 to 0.75-1.0 with different slopes.
	 * @param qual Quality score to scale
	 * @return Scaled quality value
	 */
	private static float elbaScaleQualNew(float qual) {
		if (qual >= 20.0f) {
			return 0.75f + (qual - 20.0f) / 80.0f * 0.25f; 
		} else {
			return qual / 20.0f * 0.25f; 
		}
	}


	/**
	 * Encodes variant type into categorical values for feature vector.
	 * Maps substitution=100, insertion=10, deletion=1, other=0.
	 * @param varType Variant type code
	 * @return Encoded categorical value
	 */
	private static int elbaEncodeTyp(int varType) {
		switch(varType) {
		case 0: return 100; // sub
		case 1: return 10; //insert
		case 2: return 1; //del
		default: return 0;
		}
	}


	/**
	 * Calculates variant length as absolute difference between read and reference.
	 * @param v Variant to measure
	 * @return Absolute difference between read length and reference length
	 */
	private static int elbaCalculateLen(Var v) {
		return Math.abs(v.readlen() - v.reflen());
	} 


}
