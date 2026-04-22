package cardinality;

import shared.Tools;

/**
 * LDLC calibration v2: tries multiple blend strategies.
 * Key insight: blend in RESIDUAL space, not estimate space.
 * LDLC = DLC + alpha*(L2_normalized - DLC)
 */
public class LDLCCalibrate2 {

	/*--------------------------------------------------------------*/
	/*----------------        Static Main          ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		int buckets=512;
		int numDDLs=4000;
		for(String arg : args){
			final String[] ab=arg.split("=");
			if(ab[0].equals("buckets")){buckets=Integer.parseInt(ab[1]);}
			else if(ab[0].equals("ddls")){numDDLs=Integer.parseInt(ab[1]);}
		}
		final int maxCard=buckets*256;

		// Build checkpoint array
		int nc=0;
		for(int c=1; c<=maxCard; c=(int)(c*1.01)+1){nc++;}
		final int[] checks=new int[nc];
		int ci=0;
		for(int c=1; c<=maxCard; c=(int)(c*1.01)+1){checks[ci++]=c;}

		// First pass: compute the average L2/DLC ratio for normalization
		double ratioSum=0; int ratioCount=0;
		{
			for(int inst=0; inst<Math.min(numDDLs, 500); inst++){
				final UltraDynamicLogLog6 udll=new UltraDynamicLogLog6(buckets, 31, -1, 0);
				long seed=inst*999999937L+42;
				int nextCheck=0;
				for(int card=1; card<=maxCard && nextCheck<nc; card++){
					seed=Tools.hash64shift(seed);
					udll.hashAndStore(seed);
					if(card==checks[nextCheck]){
						if(card>8*buckets){
							final int minZeros=udll.getMinZeros();
							final int[] nlzCounts=buildNlzCounts(udll, buckets, minZeros);
							final CardinalityStats cs=CardinalityStats.fromNlzCounts(
								nlzCounts, buckets, 0L, null, 0, null, null);
							final double trueDLC=cs.dlcEst;
							final double[] layers=udll.udlcEstimate();
							final double l2=layers[1];
							if(l2>0 && trueDLC>0){
								ratioSum+=trueDLC/l2;
								ratioCount++;
							}
						}
						nextCheck++;
					}
				}
			}
		}
		final double l2Scale=ratioSum/ratioCount;
		System.err.printf("L2 normalization scale: %.4f (from %d samples)%n", l2Scale, ratioCount);

		// Accumulators for multiple blend weights
		final double[] alphas={0.0, 0.05, 0.10, 0.15, 0.18, 0.20, 0.25, 0.30, 0.40, 0.50};
		final double[][] absSum=new double[alphas.length][nc];
		final double[][] errSum=new double[alphas.length][nc];
		final double[] dlcAbsSum=new double[nc];
		final double[] fgraAbsSum=new double[nc];
		final double[] dlcErrSum=new double[nc];
		final double[] fgraErrSum=new double[nc];

		// Second pass: full run
		for(int inst=0; inst<numDDLs; inst++){
			final UltraDynamicLogLog6 udll=new UltraDynamicLogLog6(buckets, 31, -1, 0);
			long seed=inst*999999937L+42;
			int nextCheck=0;
			for(int card=1; card<=maxCard && nextCheck<nc; card++){
				seed=Tools.hash64shift(seed);
				udll.hashAndStore(seed);
				if(card==checks[nextCheck]){
					final int minZeros=udll.getMinZeros();
					final int[] nlzCounts=buildNlzCounts(udll, buckets, minZeros);
					final CardinalityStats cs=CardinalityStats.fromNlzCounts(
						nlzCounts, buckets, 0L, null, 0, null, null);
					final double trueDLC=cs.dlcEst;
					final double fgra=udll.fgraEstimate();

					final double[] layers=udll.udlcEstimate();
					final double l2=layers[1];
					final int b2=countHistoryBuckets(udll, buckets);
					final double coverage=(double)b2/buckets;

					// DLC and FGRA
					dlcAbsSum[nextCheck]+=Math.abs(trueDLC-card)/(double)card;
					fgraAbsSum[nextCheck]+=Math.abs(fgra-card)/(double)card;
					dlcErrSum[nextCheck]+=(trueDLC-card)/(double)card;
					fgraErrSum[nextCheck]+=(fgra-card)/(double)card;

					// Try each alpha: LDLC = DLC + alpha*(L2_scaled - DLC)
					// With coverage ramp: no blend below 30%, full at 60%
					for(int ai=0; ai<alphas.length; ai++){
						double ldlc;
						if(l2<=0 || coverage<0.30){
							ldlc=trueDLC;
						}else{
							final double l2Norm=l2*l2Scale;
							double effectiveAlpha=alphas[ai];
							if(coverage<0.60){
								effectiveAlpha*=(coverage-0.30)/0.30;
							}
							ldlc=trueDLC+effectiveAlpha*(l2Norm-trueDLC);
						}
						absSum[ai][nextCheck]+=Math.abs(ldlc-card)/(double)card;
						errSum[ai][nextCheck]+=(ldlc-card)/(double)card;
					}
					nextCheck++;
				}
			}
		}

		// Header
		final StringBuilder sb=new StringBuilder("TrueCard\tDLC_abs\tFGRA_abs\tDLC_err\tFGRA_err");
		for(double a : alphas){sb.append(String.format("\ta%.2f_abs\ta%.2f_err", a, a));}
		System.out.println(sb);

		for(int i=0; i<nc; i++){
			final StringBuilder line=new StringBuilder();
			line.append(String.format("%d\t%.6f\t%.6f\t%.6f\t%.6f",
				checks[i],
				dlcAbsSum[i]/numDDLs,
				fgraAbsSum[i]/numDDLs,
				dlcErrSum[i]/numDDLs,
				fgraErrSum[i]/numDDLs));
			for(int ai=0; ai<alphas.length; ai++){
				line.append(String.format("\t%.6f\t%.6f",
					absSum[ai][i]/numDDLs,
					errSum[ai][i]/numDDLs));
			}
			System.out.println(line);
		}

		// Summary: mean abs error for card > 4*B
		System.err.println("\n=== Summary (card > "+4*buckets+") ===");
		int hiStart=0;
		for(int i=0; i<nc; i++){if(checks[i]>4*buckets){hiStart=i; break;}}
		final int hiCount=nc-hiStart;

		double dlcMean=0, fgraMean=0;
		for(int i=hiStart; i<nc; i++){
			dlcMean+=dlcAbsSum[i]/numDDLs;
			fgraMean+=fgraAbsSum[i]/numDDLs;
		}
		dlcMean/=hiCount; fgraMean/=hiCount;
		System.err.printf("DLC:  %.4f%%  FGRA: %.4f%%%n", dlcMean*100, fgraMean*100);

		for(int ai=0; ai<alphas.length; ai++){
			double mean=0;
			for(int i=hiStart; i<nc; i++){
				mean+=absSum[ai][i]/numDDLs;
			}
			mean/=hiCount;
			System.err.printf("alpha=%.2f: %.4f%% (%.1f%% vs DLC, %.1f%% vs FGRA)%n",
				alphas[ai], mean*100,
				(1-mean/dlcMean)*100,
				(1-mean/fgraMean)*100);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods       ----------------*/
	/*--------------------------------------------------------------*/

	/** Build NLZ histogram from UDLL6 registers. */
	private static int[] buildNlzCounts(UltraDynamicLogLog6 udll, int buckets, int minZeros){
		final int[] nlzCounts=new int[64];
		for(int i=0; i<buckets; i++){
			final int reg=udll.getRegPublic(i);
			if(reg==0){continue;}
			final int nlzPart=reg>>>2;
			final int absNlz=nlzPart-UltraDynamicLogLog6.HISTORY_MARGIN+minZeros;
			if(absNlz>=0 && absNlz<64){nlzCounts[absNlz]++;}
		}
		return nlzCounts;
	}

	/** Count buckets with at least one history bit set. */
	private static int countHistoryBuckets(UltraDynamicLogLog6 udll, int buckets){
		int count=0;
		for(int i=0; i<buckets; i++){
			final int reg=udll.getRegPublic(i);
			if(reg>0 && (reg&3)!=0){count++;}
		}
		return count;
	}
}
