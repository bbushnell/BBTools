package var2;

import shared.Tools;

/**
 * Legacy 11-feature extractor retained from Elba's model work.
 * Applies power/log transforms and fixed min-max scaling without clamping or
 * finite-value validation. This layout differs from the 32-slot ELBA branch
 * in FeatureVectorMaker. Current Java scoring callers use VectorUMP45 directly.
 * @author Brian Bushnell
 */
//TODO: Possible bug [var2/VectorElba#001] (LOW, dead) - no callers of this class found in current Java sources
//(superseded intern scaffolding; the live feature vector is VectorUMP45). elbaScaleQualNew below is also
//never called even within this class. Cleanup candidate; bundle with the FeatureVectorMaker dead-dispatch.
public class VectorElba{
	/**
	 * Allocates an 11-slot vector and a fresh set of scaling bounds.
	 * Slots contain transformed strand ratio, mean read length and mean identity;
	 * mean map/base qualities and end distance; allele fraction and its complement;
	 * transformed coverage per supporting read; type code; and absolute length change.
	 * Identity is taken in stored per-mille units. The guarded means are zero for
	 * nonpositive allele count. The variant is read without caching a revised AF.
	 *
	 * @param v Non-null variant with accumulated statistics
	 * @param pairingRate Unused compatibility parameter
	 * @param totalQualityAvg Unused compatibility parameter
	 * @param totalMapqAvg Unused compatibility parameter
	 * @param readLengthAvg Unused compatibility parameter
	 * @param ploidy Unused compatibility parameter
	 * @param map Unused compatibility parameter; may be null
	 * @return Newly allocated vector; entries are not restricted to [0,1]
	 */
	public static float[] makeElbaVector(Var v, double pairingRate, double totalQualityAvg,
			double totalMapqAvg, double readLengthAvg, int ploidy, ScafMap map){
		ElbaMinMaxValues minMaxVals=new ElbaMinMaxValues();

		float[] vec=new float[11];
		int idx=0;

		final int alleleCount=v.alleleCount();
		final int safeCount=Tools.max(1, alleleCount);

		// raw values
		float raw_sb=(float)v.strandRatio();
		float raw_avg_read_length=(alleleCount>0 ? (float)(v.lengthSum/(double)safeCount) : 0f);
		float raw_avg_identity=(alleleCount>0 ? (float)(v.idSum/(double)safeCount) : 0f);
		float raw_avg_mapq=(alleleCount>0 ? (float)(v.mapQSum/(double)safeCount) : 0f);
		float raw_avg_baseq=(alleleCount>0 ? (float)(v.baseQSum/(double)safeCount) : 0f);
		float raw_avg_enddist=(alleleCount>0 ? (float)(v.endDistSum/(double)safeCount) : 0f);
		float raw_avg_depth_per_allele=(float)(v.coverage()/(double)Tools.max(1, alleleCount));

		float transformed_sb=elbaTransformLeftSkewed(raw_sb);
		vec[idx++]=elbaNormalize(transformed_sb, minMaxVals.min_transformed_sb, minMaxVals.max_transformed_sb);

		float transformed_avg_read_length=elbaTransformLeftSkewed(raw_avg_read_length);
		vec[idx++]=elbaNormalize(transformed_avg_read_length, minMaxVals.min_transformed_avg_read_length, minMaxVals.max_transformed_avg_read_length);

		float transformed_avg_identity=elbaTransformLeftSkewed(raw_avg_identity);
		vec[idx++]=elbaNormalize(transformed_avg_identity, minMaxVals.min_transformed_avg_identity, minMaxVals.max_transformed_avg_identity);

		vec[idx++]=elbaNormalize(raw_avg_mapq, minMaxVals.min_avg_mapq, minMaxVals.max_avg_mapq);

		vec[idx++]=(float)(1.0-v.alleleFraction());

		vec[idx++]=(float)v.alleleFraction(); //AF

		vec[idx++]=elbaNormalize(raw_avg_baseq, minMaxVals.min_avg_baseq, minMaxVals.max_avg_baseq); //Avg_BaseQ

		vec[idx++]=elbaNormalize(raw_avg_enddist, minMaxVals.min_avg_enddist, minMaxVals.max_avg_enddist); //Avg_EndDist

		float transformed_avg_depth_per_allele=elbaTransformRightSkewed(raw_avg_depth_per_allele); //Avg_Depth_Per_Allele
		vec[idx++]=elbaNormalize(transformed_avg_depth_per_allele, minMaxVals.min_transformed_avg_depth_per_allele, minMaxVals.max_transformed_avg_depth_per_allele);

		vec[idx++]=(float)elbaEncodeTyp(v.type()); //TYP

		int len=elbaCalculateLen(v); //LEN
		vec[idx++]=(float)len;
		return vec;
	}

	/** Fixed scaling endpoints, allocated for each vector.
	 * Values outside the endpoints are extrapolated, not clipped. */
	private static class ElbaMinMaxValues{
		/** Scaling endpoints for average mapping quality */
		public float min_avg_mapq=0.0f, max_avg_mapq=45.0f;
		/** Scaling endpoints for average base quality */
		public float min_avg_baseq=0.0f, max_avg_baseq=45.0f;
		/** Scaling endpoints for average end distance */
		public float min_avg_enddist=0.0f, max_avg_enddist=75.0f;
		/** Scaling endpoints for transformed average depth per allele */
		public float min_transformed_avg_depth_per_allele=0.0f, max_transformed_avg_depth_per_allele=3.1f;
		/** Scaling endpoints for transformed average read length */
		public float min_transformed_avg_read_length=0.0f, max_transformed_avg_read_length=2.5e21f;
		/** Scaling endpoints for transformed average identity */
		public float min_transformed_avg_identity=0.0f, max_transformed_avg_identity=9.5e29f;
		/** Scaling endpoints for transformed strand bias */
		public float min_transformed_sb=0.0f, max_transformed_sb=1.0f;
	}

	/**
	 * Returns zero for nonpositive input, otherwise natural log of value+1.
	 * NaN passes through the logarithm.
	 * @param value Raw feature
	 * @return Log-transformed feature
	 */
	private static float elbaTransformRightSkewed(float value){
		if(value<=0){return 0.0f;}
		return (float)Math.log(value+1);
	}

	/**
	 * Raises the input to the tenth power and narrows the result to float.
	 * No overflow or finite-value check is performed.
	 * @param value Raw feature
	 * @return Power-transformed feature
	 */
	private static float elbaTransformLeftSkewed(float value){
		return (float)Math.pow(value, 10);
	}

	/**
	 * Returns (value-min)/(max-min), or zero when the endpoints are equal.
	 * There is no clipping or finite-value validation.
	 * @param value Feature to scale
	 * @param min_val Lower scaling endpoint
	 * @param max_val Upper scaling endpoint
	 * @return Scaled feature, which may be outside [0,1]
	 */
	private static float elbaNormalize(float value, float min_val, float max_val){
		if(max_val==min_val){return 0.0f;}
		return (value-min_val)/(max_val-min_val);
	}

	/**
	 * Unused piecewise quality transform with a discontinuity at 20.
	 * Values below 20 use qual/20*0.25; values at least 20 use
	 * 0.75+(qual-20)/80*0.25. Neither branch clips its result.
	 * @param qual Quality value
	 * @return Piecewise-scaled value
	 */
	private static float elbaScaleQualNew(float qual){
		if(qual>=20.0f){
			return 0.75f+(qual-20.0f)/80.0f*0.25f;
		}else{
			return qual/20.0f*0.25f;
		}
	}

	/**
	 * Retains the historical literal mapping 0 to 100, 1 to 10 and 2 to 1.
	 * Other inputs map to zero. These codes do not match the original biological
	 * labels under the current Var type constants; see VAR2-014 below.
	 * @param varType Value returned by Var.type()
	 * @return Historical encoded value
	 */
	private static int elbaEncodeTyp(int varType){
		//TODO: Probable bug [VAR2-014] - original labels assumed SUB=0, INS=1, DEL=2;
		//Var currently uses INS=0, NOCALL=1, SUB=2, DEL=3. No current Java caller found;
		//preserve the historical mapping until its model/schema contract is established.
		switch(varType){
		case 0: return 100; //INS under current Var constants
		case 1: return 10; //NOCALL under current Var constants
		case 2: return 1; //SUB under current Var constants
		default: return 0;
		}
	}

	/** Absolute difference between allele length and reference span. */
	private static int elbaCalculateLen(Var v){
		return Math.abs(v.readlen()-v.reflen());
	}
}
