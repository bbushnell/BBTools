package var2;

import shared.Tools;

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

	private static class ElbaMinMaxValues {
		public float min_avg_mapq = 0.0f, max_avg_mapq = 45.0f;    
		public float min_avg_baseq = 0.0f, max_avg_baseq = 45.0f;       
		public float min_avg_enddist = 0.0f, max_avg_enddist = 75.0f;  
		public float min_transformed_avg_depth_per_allele = 0.0f, max_transformed_avg_depth_per_allele = 3.1f;  
		public float min_transformed_avg_read_length = 0.0f, max_transformed_avg_read_length = 2.5e21f;       
		public float min_transformed_avg_identity = 0.0f, max_transformed_avg_identity = 9.5e29f;    
		public float min_transformed_sb = 0.0f, max_transformed_sb = 1.0f;               

	}


	private static float elbaTransformRightSkewed(float value) {
		if (value <= 0) return 0.0f;
		return (float)Math.log(value + 1);
	}


	private static float elbaTransformLeftSkewed(float value) {
		return (float)Math.pow(value, 10);
	}


	private static float elbaNormalize(float value, float min_val, float max_val) {
		if (max_val == min_val) return 0.0f;
		return (value - min_val) / (max_val - min_val);
	}


	private static float elbaScaleQualNew(float qual) {
		if (qual >= 20.0f) {
			return 0.75f + (qual - 20.0f) / 80.0f * 0.25f; 
		} else {
			return qual / 20.0f * 0.25f; 
		}
	}


	private static int elbaEncodeTyp(int varType) {
		switch(varType) {
		case 0: return 100; // sub
		case 1: return 10; //insert
		case 2: return 1; //del
		default: return 0;
		}
	}


	private static int elbaCalculateLen(Var v) {
		return Math.abs(v.readlen() - v.reflen());
	} 


}
