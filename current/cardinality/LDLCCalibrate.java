package cardinality;

import shared.Tools;

/**
 * LDLC (Layered Dynamic Linear Counting) calibration driver.
 * Outputs true DLC (from CardinalityStats), LDLC blend, L2, FGRA, and LC
 * at every checkpoint. Focus: 512 buckets (384 bytes = 3 bits/bucket equivalent).
 */
public class LDLCCalibrate {

	public static void main(String[] args){
		int buckets=512;
		int numDDLs=4000;
		float blendL2=0.18f;
		for(String arg : args){
			String[] ab=arg.split("=");
			if(ab[0].equals("buckets")){buckets=Integer.parseInt(ab[1]);}
			else if(ab[0].equals("ddls")){numDDLs=Integer.parseInt(ab[1]);}
			else if(ab[0].equals("blend")){blendL2=Float.parseFloat(ab[1]);}
		}
		final int maxCard=buckets*256;

		// Build checkpoint array (~1% spacing)
		int nc=0;
		for(int c=1; c<=maxCard; c=(int)(c*1.01)+1) nc++;
		final int[] checks=new int[nc];
		int ci=0;
		for(int c=1; c<=maxCard; c=(int)(c*1.01)+1) checks[ci++]=c;

		// Accumulators: per-checkpoint sums across instances
		final double[] dlcSum=new double[nc], l2Sum=new double[nc];
		final double[] fgraSum=new double[nc], lcSum=new double[nc];
		final double[] ldlcSum=new double[nc];
		final double[] dlcAbsSum=new double[nc], fgraAbsSum=new double[nc];
		final double[] lcAbsSum=new double[nc], ldlcAbsSum=new double[nc];
		final double[] l2AbsSum=new double[nc];
		final int[] l2CountSum=new int[nc]; // how many buckets have history

		final float blendL1=1.0f-blendL2;

		for(int inst=0; inst<numDDLs; inst++){
			UltraDynamicLogLog6 udll=new UltraDynamicLogLog6(buckets, 31, -1, 0);
			long seed=inst*999999937L+42;
			int nextCheck=0;
			for(int card=1; card<=maxCard && nextCheck<nc; card++){
				seed=Tools.hash64shift(seed);
				udll.hashAndStore(seed);
				if(card==checks[nextCheck]){
					// Get FGRA
					final double fgra=udll.fgraEstimate();

					// Build NLZ histogram from registers for true DLC
					final int minZeros=udll.getMinZeros();
					final int[] nlzCounts=new int[64];
					int emptyCount=0;
					for(int i=0; i<buckets; i++){
						int reg=udll.getRegPublic(i);
						if(reg==0){emptyCount++; continue;}
						int nlzPart=reg>>>2;
						int absNlz=nlzPart-UltraDynamicLogLog6.HISTORY_MARGIN+minZeros;
						if(absNlz>=0 && absNlz<64){nlzCounts[absNlz]++;}
					}

					// Build CardinalityStats for true DLC
					CardinalityStats cs=CardinalityStats.fromNlzCounts(
						nlzCounts, buckets, 0L, null, 0, null, null);
					final double trueDLC=cs.dlcEst;
					final double lc=(emptyCount>0 && emptyCount<buckets) ?
						(double)buckets*Math.log((double)buckets/emptyCount) : 0;

					// Get L2 from udlcEstimate
					final double[] layers=udll.udlcEstimate();
					final double l2=layers[1];
					// Count how many buckets have history (for coverage)
					int b2=countHistoryBuckets(udll, buckets);

					// Compute L2 scale factor (global constant from prior analysis)
					// At 512 buckets: ~1.293, at 1024: ~1.290
					// We use a fixed 1.29 for now
					final double l2Scaled=l2*1.29;

					// LDLC blend: use true DLC as L1, blend with scaled L2
					// Ramp: only blend L2 when coverage > 30%, full blend at > 60%
					double ldlc;
					double coverage=(double)b2/buckets;
					if(coverage<0.30 || l2<=0){
						ldlc=trueDLC;
					}else if(coverage<0.60){
						double ramp=(coverage-0.30)/0.30;
						ldlc=trueDLC*(1.0-ramp*blendL2) + l2Scaled*(ramp*blendL2);
					}else{
						ldlc=trueDLC*blendL1 + l2Scaled*blendL2;
					}

					dlcSum[nextCheck]+=trueDLC;
					l2Sum[nextCheck]+=l2;
					fgraSum[nextCheck]+=fgra;
					lcSum[nextCheck]+=lc;
					ldlcSum[nextCheck]+=ldlc;
					dlcAbsSum[nextCheck]+=Math.abs(trueDLC-card)/(double)card;
					l2AbsSum[nextCheck]+=Math.abs(l2-card)/(double)card;
					fgraAbsSum[nextCheck]+=Math.abs(fgra-card)/(double)card;
					lcAbsSum[nextCheck]+=Math.abs(lc-card)/(double)card;
					ldlcAbsSum[nextCheck]+=Math.abs(ldlc-card)/(double)card;
					l2CountSum[nextCheck]+=b2;
					nextCheck++;
				}
			}
		}

		System.out.println("TrueCard\tDLC\tLDLC\tL2\tFGRA\tLC\tDLC_abs\tLDLC_abs\tL2_abs\tFGRA_abs\tLC_abs\tL2_coverage");
		for(int i=0; i<nc; i++){
			System.out.printf("%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.4f%n",
				checks[i],
				dlcSum[i]/numDDLs,
				ldlcSum[i]/numDDLs,
				l2Sum[i]/numDDLs,
				fgraSum[i]/numDDLs,
				lcSum[i]/numDDLs,
				dlcAbsSum[i]/numDDLs,
				ldlcAbsSum[i]/numDDLs,
				l2AbsSum[i]/numDDLs,
				fgraAbsSum[i]/numDDLs,
				lcAbsSum[i]/numDDLs,
				(double)l2CountSum[i]/(numDDLs*buckets));
		}
	}

	/** Count buckets with at least one history bit set. */
	private static int countHistoryBuckets(UltraDynamicLogLog6 udll, int buckets){
		int count=0;
		for(int i=0; i<buckets; i++){
			int reg=udll.getRegPublic(i);
			if(reg>0 && (reg&3)!=0) count++;
		}
		return count;
	}
}
