package var2;

import shared.Tools;

/**
 * Builds 33 scaled features for variant scoring. Features are not all bounded to
 * [0,1]; type/platform indicators, logarithms, ratios and affine transforms use
 * the layout shown in makeVector. The composite-score feature is optional.
 * Performs no file I/O itself, but reference-derived features depend on the supplied
 * scaffold map and its bases. May populate the variant's revised-allele-fraction cache.
 * Callers must coordinate shared settings, variants and reusable output arrays.
 * @author UMP45
 * @contributor Sucrose
 */
public class VectorUMP45{

	/*--------------------------------------------------------------*/
	/*----------------        Configuration        ----------------*/
	/*--------------------------------------------------------------*/

	/** Number of feature slots populated by this implementation. */
	public static final int DIMS=33;

	/** Indices for the four platform indicators in slots4–7. */
	public static final int PLATFORM_ILLUMINA=0, PLATFORM_PACBIO=1,
			PLATFORM_NANOPORE=2, PLATFORM_ROCHE=3;

	/** Process-wide platform setting; vector encoding clamps it to0–3. */
	public static int platform=PLATFORM_ILLUMINA;
	/** Enables the non-NN score in slot29; otherwise that slot remains zero. */
	public static boolean includeScore=false;

	/*--------------------------------------------------------------*/
	/*----------------     Vector Construction     ----------------*/
	/*--------------------------------------------------------------*/

	/** Allocates a DIMS-length array and delegates to the reusable-buffer overload.
	 * The returned array is owned by the caller. */
	public static float[] makeVector(Var v, double pairingRate, double totalQualityAvg,
			double totalMapqAvg, double readLengthAvg, int ploidy, ScafMap map){
		return makeVector(v, pairingRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, map, new float[DIMS]);
	}

	/**
	 * Clears the entire supplied array, then fills slots0–32 and returns that same array.
	 * Any extra slots remain zero. Only SUB/INS/DEL set a variant-type indicator;
	 * other types leave all three type slots zero. Platform encoding uses the current
	 * static platform, clamped to0–3. Slot31 is reserved and stays zero.
	 *
	 * Revised AF may be cached into v; an existing cache value is reused. Reference
	 * features use map and its bases. With null map, homopolymer count is zero and
	 * the distance feature uses v.start. A metadata-only map can produce different
	 * features from one with loaded bases. Callers supply their model's reference context.
	 * The final assertion checks finiteness, not feature ranges or model compatibility.
	 * Exceptions can leave a partially filled array and updated variant cache.
	 *
	 * @param v Nonnull variant with the coverage/statistics needed by its accessors
	 * @param pairingRate Dataset pairing rate, used only by the optional score feature
	 * @param totalQualityAvg Dataset base-quality mean, used only by the optional score
	 * @param totalMapqAvg Dataset mapping-quality mean, used only by the optional score
	 * @param readLengthAvg Used by the revised-AF calculation and optional score
	 * @param ploidy Positive sample ploidy, encoded as its inverse and a greater-than-one flag
	 * @param map Reference registry for sequence-derived features, or null for the fallbacks above
	 * @param vec Nonnull caller-owned array with at least DIMS entries; all entries are cleared
	 * @return The supplied vec reference
	 */
	public static float[] makeVector(Var v, double pairingRate, double totalQualityAvg,
			double totalMapqAvg, double readLengthAvg, int ploidy, ScafMap map, float[] vec){

		java.util.Arrays.fill(vec, 0f);
		final int type=v.type();
		final int count=v.alleleCount();
		final double af=v.alleleFraction();

		// 0: Inverse ploidy
		vec[0]=1f/ploidy;

		// 1-3: Variant type 1-hot
		if(type==Var.SUB){vec[1]=1;}
		else if(type==Var.INS){vec[2]=1;}
		else if(type==Var.DEL){vec[3]=1;}

		// 4-7: Platform 1-hot
		vec[4+Tools.mid(0, platform, 3)]=1;

		// 8: Total depth
		vec[8]=log2p1(v.coverage())*INV8;

		// 9: Allelic depth
		vec[9]=log2p1(count)*INV8;

		// 10: Allele fraction
		vec[10]=(float)af;

		// 11: Revised allele fraction
		vec[11]=(float)v.revisedAlleleFraction(af, readLengthAvg);

		// 12-13: Mapping quality
		vec[12]=(count>0) ? (float)(v.mapQAvg()*INV40) : 0f;
		vec[13]=v.mapQMax*INV40;

		// 14-15: Base quality
		vec[14]=(count>0) ? (float)(v.baseQAvg()*INV40) : 0f;
		vec[15]=v.baseQMax*INV40;

		// 16-17: Identity (2x-1 maps 0.5-1.0 to 0-1)
		vec[16]=(count>0) ? (float)(2*v.identityAvg()*0.001-1) : 0f;
		vec[17]=(float)(2*v.idMax*0.001-1);

		// 18-19: End distance
		vec[18]=(count>0) ? log2p1(v.edistAvg())*INV4 : 0f;
		vec[19]=log2p1(v.endDistMax)*INV4;

		// 20: Read length average
		vec[20]=log2p1(count>0 ? v.lengthAvg() : 0)*INV4;

		// 21: Event length
		int eventLen=Tools.max(v.reflen(), v.readlen());
		vec[21]=log2p1(eventLen)*INV8;

		// 22: Strand ratio (raw)
		vec[22]=(float)v.strandRatio();

		// 23: Adjusted strand-balance score
		vec[23]=(float)VarProb.eventProb(v.allelePlusCount(), v.alleleMinusCount());

		// 24: Read bias ratio (raw)
		{
			int r1=v.r1AlleleCount(), r2=v.r2AlleleCount();
			vec[24]=(r1+r2==0) ? 1f : (float)(Tools.min(r1, r2)+1)/(float)Tools.max(r1, r2);
		}

		// 25: Adjusted read-balance score
		vec[25]=(float)VarProb.eventProb(v.r1AlleleCount(), v.r2AlleleCount());

		// 26: Nearby variant count
		vec[26]=1f/(Math.max(v.nearbyVarCount, 0)+1);

		// 27: Proper pair rate
		vec[27]=(count==0) ? 0f : (float)v.properPairRate();

		// 28: Homopolymer count
		int hpc=v.homopolymerCount(map);
		vec[28]=1f/(hpc+1);

		// 29: Non-NN composite score (optional — disabled by default to force NN to learn its own scoring)
		if(includeScore){
			vec[29]=(float)v.score(pairingRate, totalQualityAvg, totalMapqAvg,
					readLengthAvg, 1, ploidy, map, null);
		}

		// 30: Contig end distance
		int ced=(map==null ? v.start : v.contigEndDist(map));
		vec[30]=log2p1(ced)*INV8;

		// 31: Reserved (AF probability — needs copy-count-aware binomial model)
		vec[31]=0;

		// 32: Ploidy>1 flag (0 if haploid, 1 if ploidy>1) — pairs with dim 0 (1/ploidy)
		vec[32]=(ploidy>1)?1f:0f;

		assert(validateVector(vec, v, count, af, pairingRate, totalQualityAvg,
				totalMapqAvg, readLengthAvg, ploidy, map));

		return vec;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Validation          ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks every supplied slot for NaN/infinity, printing context for the first failure.
	 * Invoked by an assertion, so this scan is skipped when assertions are disabled.
	 * Diagnostic accessors can also read the reference or populate the revised-AF cache. */
	private static boolean validateVector(float[] vec, Var v, int count, double af,
			double pairingRate, double totalQualityAvg, double totalMapqAvg,
			double readLengthAvg, int ploidy, ScafMap map){
		for(int i=0; i<vec.length; i++){
			if(Float.isNaN(vec[i]) || Float.isInfinite(vec[i])){
				StringBuilder sb=new StringBuilder();
				sb.append("VectorUMP45: vec[").append(i).append("]=").append(vec[i]);
				sb.append("\n  variant: ").append(v);
				sb.append("\n  type=").append(v.type()).append(" count=").append(count);
				sb.append(" af=").append(af).append(" coverage=").append(v.coverage());
				sb.append("\n  mapQAvg=").append(v.mapQAvg()).append(" baseQAvg=").append(v.baseQAvg());
				sb.append(" idAvg=").append(v.identityAvg()).append(" idMax=").append(v.idMax);
				sb.append("\n  edistAvg=").append(v.edistAvg()).append(" edistMax=").append(v.endDistMax);
				sb.append(" lengthAvg=").append(count>0 ? v.lengthAvg() : 0);
				sb.append("\n  strandRatio=").append(v.strandRatio());
				sb.append(" plusCount=").append(v.allelePlusCount()).append(" minusCount=").append(v.alleleMinusCount());
				sb.append("\n  r1=").append(v.r1AlleleCount()).append(" r2=").append(v.r2AlleleCount());
				sb.append(" properPairRate=").append(count==0 ? "N/A" : ""+v.properPairRate());
				sb.append(" nearbyVarCount=").append(v.nearbyVarCount);
				sb.append("\n  hpc=").append(v.homopolymerCount(map));
				sb.append(" contigEndDist=").append(map==null ? v.start : v.contigEndDist(map));
				sb.append("\n  revisedAF=").append(v.revisedAlleleFraction(af, readLengthAvg));
				sb.append("\n  ploidy=").append(ploidy).append(" pairingRate=").append(pairingRate);
				sb.append(" totalQAvg=").append(totalQualityAvg).append(" totalMQAvg=").append(totalMapqAvg);
				sb.append(" readLenAvg=").append(readLengthAvg);
				System.err.println(sb);
				return false;
			}
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Utility Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns log(x+1) scaled by the stored float approximation to 1/ln(2).
	 * Performs no domain check; invalid inputs propagate non-finite results. */
	private static float log2p1(double x){
		return (float)(Math.log(x+1)*INV_LN2);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Base conversion for logarithmic features. */
	private static final float INV_LN2=(float)(1.0/Math.log(2));
	/** Scales for logarithmic lengths/distances and quality features. */
	private static final float INV4=1f/4;
	private static final float INV8=1f/8;
	private static final float INV40=1f/40;
}
