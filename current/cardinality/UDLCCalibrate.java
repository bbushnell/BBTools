package cardinality;

import shared.Tools;

/**
 * UDLC calibration: outputs raw L1, L2, L3 estimates at every checkpoint.
 * Blending and normalization done in Python.
 */
public class UDLCCalibrate {

	public static void main(String[] args){
		int buckets=2048;
		int numDDLs=4000;
		for(String arg : args){
			String[] ab=arg.split("=");
			if(ab[0].equals("buckets")){buckets=Integer.parseInt(ab[1]);}
			else if(ab[0].equals("ddls")){numDDLs=Integer.parseInt(ab[1]);}
		}
		final int maxCard=buckets*256;

		// Build checkpoint array
		int nc=0;
		for(int c=1; c<=maxCard; c=(int)(c*1.01)+1) nc++;
		final int[] checks=new int[nc];
		int ci=0;
		for(int c=1; c<=maxCard; c=(int)(c*1.01)+1) checks[ci++]=c;

		// Accumulators
		final double[] fgraSum=new double[nc];
		final double[] l1Sum=new double[nc], l2Sum=new double[nc], l3Sum=new double[nc];
		final double[] fgraAbsSum=new double[nc];
		final double[] l1AbsSum=new double[nc];

		for(int inst=0; inst<numDDLs; inst++){
			UltraDynamicLogLog6 udll=new UltraDynamicLogLog6(buckets, 31, -1, 0);
			long seed=inst*999999937L+42;
			int nextCheck=0;
			for(int card=1; card<=maxCard && nextCheck<nc; card++){
				seed=Tools.hash64shift(seed);
				udll.hashAndStore(seed);
				if(card==checks[nextCheck]){
					final double fgra=udll.fgraEstimate();
					final double[] layers=udll.udlcEstimate();
					fgraSum[nextCheck]+=fgra;
					fgraAbsSum[nextCheck]+=Math.abs(fgra-card)/(double)card;
					l1Sum[nextCheck]+=layers[0];
					l1AbsSum[nextCheck]+=Math.abs(layers[0]-card)/(double)card;
					l2Sum[nextCheck]+=layers[1];
					l3Sum[nextCheck]+=layers[2];
					nextCheck++;
				}
			}
		}

		System.out.println("TrueCard\tFGRA\tL1\tL2\tL3\tFGRA_abs\tL1_abs");
		for(int i=0; i<nc; i++){
			System.out.printf("%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.6f\t%.6f%n",
				checks[i],
				fgraSum[i]/numDDLs,
				l1Sum[i]/numDDLs,
				l2Sum[i]/numDDLs,
				l3Sum[i]/numDDLs,
				fgraAbsSum[i]/numDDLs,
				l1AbsSum[i]/numDDLs);
		}
	}
}
