package cardinality;

import shared.Tools;

/**
 * Dense UDLC test: 1% checkpoints from 8B to 16B (16384-32768 at 2048 buckets).
 * Prints per-checkpoint signed error for L1, L2, and blends to show sine wave phase.
 */
public class UDLCDense {

	public static void main(String[] args){
		final int buckets=2048;
		final int numInstances=2000;
		final int startCard=buckets*8;   // 16384
		final int endCard=buckets*16;    // 32768

		// Dense checkpoints: ~1% increments
		int numChecks=0;
		for(int card=startCard; card<=endCard; card=(int)(card*1.01)+1) numChecks++;
		final int[] checkpoints=new int[numChecks];
		int ci=0;
		for(int card=startCard; card<=endCard; card=(int)(card*1.01)+1) checkpoints[ci++]=card;

		// Accumulators
		final double[] fgraErr=new double[numChecks];
		final double[] l1Err=new double[numChecks];
		final double[] l2Err=new double[numChecks];
		final double[] l3Err=new double[numChecks];
		final double[] fgraAbs=new double[numChecks];
		final double[] l1Abs=new double[numChecks];
		final double[] l2Abs=new double[numChecks];
		final double[] l3Abs=new double[numChecks];

		for(int inst=0; inst<numInstances; inst++){
			UltraDynamicLogLog6 udll=new UltraDynamicLogLog6(buckets, 31, -1, 0);
			long seed=inst*999999937L+42;
			int nextCheck=0;
			for(int card=1; card<=endCard && nextCheck<numChecks; card++){
				seed=Tools.hash64shift(seed);
				udll.hashAndStore(seed);
				if(card==checkpoints[nextCheck]){
					final double fgra=udll.fgraEstimate();
					final double[] udlc=udll.udlcEstimate();
					final double dlc1=udlc[0], dlc2=udlc[1], dlc3=udlc[2];
					fgraErr[nextCheck]+=(fgra-card)/card;
					l1Err[nextCheck]+=(dlc1-card)/card;
					if(dlc2>0) l2Err[nextCheck]+=(dlc2-card)/card;
					if(dlc3>0) l3Err[nextCheck]+=(dlc3-card)/card;
					fgraAbs[nextCheck]+=Math.abs(fgra-card)/card;
					l1Abs[nextCheck]+=Math.abs(dlc1-card)/card;
					if(dlc2>0) l2Abs[nextCheck]+=Math.abs(dlc2-card)/card;
					if(dlc3>0) l3Abs[nextCheck]+=Math.abs(dlc3-card)/card;
					nextCheck++;
				}
			}
		}

		// Print header and data
		System.out.println("TrueCard\tFGRA_signed\tL1_signed\tL2_signed\tL3_signed\tFGRA_abs\tL1_abs\tL2_abs\tL3_abs");
		for(int i=0; i<numChecks; i++){
			System.out.printf("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f%n",
				checkpoints[i],
				fgraErr[i]/numInstances, l1Err[i]/numInstances,
				l2Err[i]/numInstances, l3Err[i]/numInstances,
				fgraAbs[i]/numInstances, l1Abs[i]/numInstances,
				l2Abs[i]/numInstances, l3Abs[i]/numInstances);
		}
	}
}
