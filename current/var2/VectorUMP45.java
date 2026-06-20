package var2;

import shared.Tools;

/**
 * Feature vector generator for variant quality prediction.
 * Produces a 33-dimensional normalized vector with deterministic transforms
 * (no external data files required). Includes ploidy, platform, variant type,
 * and the non-NN composite score as inputs.
 * @author UMP45
 */
public class VectorUMP45 {

	public static final int DIMS=33;

	public static final int PLATFORM_ILLUMINA=0, PLATFORM_PACBIO=1,
			PLATFORM_NANOPORE=2, PLATFORM_RESERVED=3;

	public static int platform=PLATFORM_ILLUMINA;
	public static boolean includeScore=false;

	public static float[] makeVector(Var v, double pairingRate, double totalQualityAvg,
			double totalMapqAvg, double readLengthAvg, int ploidy, ScafMap map){
		return makeVector(v, pairingRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, map, new float[DIMS]);
	}

	/**
	 * As makeVector(...), but fills the caller-provided vec (length DIMS) instead of allocating one,
	 * for per-thread reuse in hot loops. The array is zeroed first because several dims (type/platform
	 * 1-hot, reserved) rely on zero-initialization.
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
		vec[12]=(float)(v.mapQAvg()*INV40);
		vec[13]=v.mapQMax*INV40;

		// 14-15: Base quality
		vec[14]=(float)(v.baseQAvg()*INV40);
		vec[15]=v.baseQMax*INV40;

		// 16-17: Identity (2x-1 maps 0.5-1.0 to 0-1)
		vec[16]=(float)(2*v.identityAvg()*0.001-1);
		vec[17]=(float)(2*v.idMax*0.001-1);

		// 18-19: End distance
		vec[18]=log2p1(v.edistAvg())*INV4;
		vec[19]=log2p1(v.endDistMax)*INV4;

		// 20: Read length average
		vec[20]=log2p1(count>0 ? v.lengthAvg() : 0)*INV4;

		// 21: Event length
		int eventLen=Tools.max(v.reflen(), v.readlen());
		vec[21]=log2p1(eventLen)*INV8;

		// 22: Strand ratio (raw)
		vec[22]=(float)v.strandRatio();

		// 23: Strand bias probability
		vec[23]=(float)VarProb.eventProb(v.allelePlusCount(), v.alleleMinusCount());

		// 24: Read bias ratio (raw)
		{
			int r1=v.r1AlleleCount(), r2=v.r2AlleleCount();
			vec[24]=(r1+r2==0) ? 1f : (float)(Tools.min(r1, r2)+1)/(float)Tools.max(r1, r2);
		}

		// 25: Read bias probability
		vec[25]=(float)VarProb.eventProb(v.r1AlleleCount(), v.r2AlleleCount());

		// 26: Nearby variant count
		vec[26]=1f/(v.nearbyVarCount+1);

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

	private static float log2p1(double x){
		return (float)(Math.log(x+1)*INV_LN2);
	}

	private static final float INV_LN2=(float)(1.0/Math.log(2));
	private static final float INV4=1f/4;
	private static final float INV8=1f/8;
	private static final float INV40=1f/40;
}
