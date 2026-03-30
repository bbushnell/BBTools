package cardinality;

import shared.Tools;

/**
 * LDLC calibration v3: outputs DLC, LDLC, FGRA, HC, HLL, and
 * history-corrected Hybrid from UDLL6 registers.
 * Focus: 512 buckets. All estimates computed inline — no pipeline changes needed.
 */
public class LDLCCalibrate3 {

	public static void main(String[] args){
		int buckets=512;
		int numDDLs=4000;
		for(String arg : args){
			String[] ab=arg.split("=");
			if(ab[0].equals("buckets")){buckets=Integer.parseInt(ab[1]);}
			else if(ab[0].equals("ddls")){numDDLs=Integer.parseInt(ab[1]);}
		}
		final int maxCard=buckets*256;
		final int p=Integer.numberOfTrailingZeros(buckets); // bucketBits

		// Build per-tier state multipliers (same as UltraLogLog8.summarize)
		final double[][] tierMult=buildTierMult();

		// Build checkpoint array
		int nc=0;
		for(int c=1; c<=maxCard; c=(int)(c*1.01)+1) nc++;
		final int[] checks=new int[nc];
		int ci=0;
		for(int c=1; c<=maxCard; c=(int)(c*1.01)+1) checks[ci++]=c;

		// Accumulators
		final double[] fgraAbsSum=new double[nc];
		final double[] hllAbsSum=new double[nc];
		final double[] dlcAbsSum=new double[nc];
		final double[] histHybAbsSum=new double[nc];
		// Signed error too
		final double[] fgraErrSum=new double[nc];
		final double[] hllErrSum=new double[nc];
		final double[] dlcErrSum=new double[nc];
		final double[] histHybErrSum=new double[nc];

		for(int inst=0; inst<numDDLs; inst++){
			UltraDynamicLogLog6 udll=new UltraDynamicLogLog6(buckets, 31, -1, 0);
			long seed=inst*999999937L+42;
			int nextCheck=0;
			for(int card=1; card<=maxCard && nextCheck<nc; card++){
				seed=Tools.hash64shift(seed);
				udll.hashAndStore(seed);
				if(card==checks[nextCheck]){
					final int minZeros=udll.getMinZeros();

					// FGRA
					final double fgra=udll.fgraEstimate();

					// Build NLZ counts + uncorrected HLL sum + state-corrected HLL sum
					final int[] nlzCounts=new int[64];
					double hllSum=0;        // uncorrected
					double hllSumCorr=0;    // state-CF corrected
					int count=0;
					for(int i=0; i<buckets; i++){
						int reg=udll.getRegPublic(i);
						if(reg==0) continue;
						int nlzPart=reg>>>2;
						int state=reg&3;
						int absNlz=nlzPart-UltraDynamicLogLog6.HISTORY_MARGIN+minZeros;
						if(absNlz<0) absNlz=0;
						if(absNlz<64) nlzCounts[absNlz]++;
						count++;

						// Uncorrected HLL contribution
						hllSum+=Math.pow(2.0, -absNlz);

						// State-corrected contribution
						int tier=Math.min(absNlz, 63);
						double m=tierMult[tier][state];
						hllSumCorr+=Math.pow(2.0, -absNlz)*m;
					}

					// Uncorrected HLL: alpha_m * count^2 / hllSum
					final double alpha_m=0.7213/(1.0+1.079/buckets);
					final double hll=(count>0 && hllSum>0) ? alpha_m*count*(double)count/hllSum : 0;

					// State-corrected Mean: (count+B)/(2B) * 2*count^2/hllSumCorr
					final double corrMean=(count>0 && hllSumCorr>0) ?
						((double)(count+buckets)/(2.0*buckets)) * 2.0*count*(double)count/hllSumCorr : 0;

					// True DLC from CardinalityStats
					int emptyCount=buckets-count;
					CardinalityStats cs=CardinalityStats.fromNlzCounts(
						nlzCounts, buckets, udll.microIndex, null, 0, null, null);
					final double trueDLC=cs.dlcEst;
					final double lcMin=cs.lcMin;

					// History-corrected Hybrid: blend LCmin → corrMean using DLC as switch
					final double hb0=0.20*buckets, hb1=5.0*buckets;
					double histHyb;
					if(trueDLC<=hb0){
						histHyb=lcMin;
					}else if(trueDLC<hb1){
						double t=Math.log(trueDLC/hb0)/Math.log(hb1/hb0);
						histHyb=(1-t)*lcMin+t*corrMean;
					}else{
						histHyb=corrMean;
					}

					// Accumulate
					double c=(double)card;
					fgraAbsSum[nextCheck]+=Math.abs(fgra-c)/c;
					hllAbsSum[nextCheck]+=Math.abs(hll-c)/c;
					dlcAbsSum[nextCheck]+=Math.abs(trueDLC-c)/c;
					histHybAbsSum[nextCheck]+=Math.abs(histHyb-c)/c;
					fgraErrSum[nextCheck]+=(fgra-c)/c;
					hllErrSum[nextCheck]+=(hll-c)/c;
					dlcErrSum[nextCheck]+=(trueDLC-c)/c;
					histHybErrSum[nextCheck]+=(histHyb-c)/c;

					nextCheck++;
				}
			}
		}

		System.out.println("TrueCard\tFGRA_abs\tHLL_abs\tDLC_abs\tHistHyb_abs\tFGRA_err\tHLL_err\tDLC_err\tHistHyb_err");
		for(int i=0; i<nc; i++){
			System.out.printf("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f%n",
				checks[i],
				fgraAbsSum[i]/numDDLs,
				hllAbsSum[i]/numDDLs,
				dlcAbsSum[i]/numDDLs,
				histHybAbsSum[i]/numDDLs,
				fgraErrSum[i]/numDDLs,
				hllErrSum[i]/numDDLs,
				dlcErrSum[i]/numDDLs,
				histHybErrSum[i]/numDDLs);
		}
	}

	/** Build per-tier, per-state multipliers matching UltraLogLog8. */
	private static double[][] buildTierMult(){
		final double[] baseCF=UltraLogLog8.STATE_CF_ACTIVE;
		final double offset=UltraLogLog8.STATE_CF_OFFSET;
		final double power=UltraLogLog8.STATE_POWER;
		final double[][] tierMult=new double[64][4];
		for(int t=0; t<64; t++){
			double[] cf;
			double off;
			if(t==0){cf=UltraLogLog8.TIER0_CF; off=UltraLogLog8.TIER0_CF_OFFSET;}
			else if(t==1){cf=UltraLogLog8.TIER1_CF; off=UltraLogLog8.TIER1_CF_OFFSET;}
			else if(t==2){cf=UltraLogLog8.TIER2_CF; off=UltraLogLog8.TIER2_CF_OFFSET;}
			else{cf=baseCF; off=offset;}
			for(int s=0; s<4; s++){
				tierMult[t][s]=Math.pow(2.0, -(cf[s]*power+off));
			}
		}
		return tierMult;
	}
}
