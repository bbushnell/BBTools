package cardinality;

import shared.Tools;

/**
 * Debug: 8-bucket UDLL, print registers and DLC at each layer.
 *
 * @author Brian Bushnell
 */
public class UDLCDebug {

	/*--------------------------------------------------------------*/
	/*----------------        Static Main          ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		final int buckets=8;

		// Print averages at checkpoints
		final int[] checkpoints={4, 8, 16, 32, 64, 128, 256};
		System.out.println("=== Averaged over 10000 instances, 8 buckets ===");
		System.out.println("Card\tFGRA\tDLC_L1\tDLC_L2\tDLC_L3");
		for(int card : checkpoints){
			final int N=10000;
			double fgraSum=0, dlc1Sum=0, dlc2Sum=0, dlc3Sum=0;
			for(int inst=0; inst<N; inst++){
				final UltraDynamicLogLog6 udll=new UltraDynamicLogLog6(buckets, 31, -1, 0);
				long seed=inst*999999937L+card*12345L;
				for(int c=0; c<card; c++){
					seed=Tools.hash64shift(seed);
					udll.hashAndStore(seed);
				}
				fgraSum+=udll.fgraEstimate();
				final double[] udlc=udll.udlcEstimate();
				dlc1Sum+=udlc[0];
				dlc2Sum+=udlc[1];
				dlc3Sum+=udlc[2];
			}
			System.out.printf("%d\t%.1f\t%.1f\t%.1f\t%.1f%n",
				card, fgraSum/N, dlc1Sum/N, dlc2Sum/N, dlc3Sum/N);
		}

		// Print detailed single instances
		for(int card : new int[]{16, 32, 64}){
			System.out.printf("%n=== Single instance, card=%d, 8 buckets ===%n", card);
			final UltraDynamicLogLog6 udll=new UltraDynamicLogLog6(buckets, 31, -1, 0);
			long seed=42;
			for(int c=0; c<card; c++){
				seed=Tools.hash64shift(seed);
				udll.hashAndStore(seed);
			}
			udll.printRegisters();
			final double[] udlc=udll.udlcEstimate();
			System.out.printf("FGRA=%.2f  DLC_L1=%.2f  DLC_L2=%.2f  DLC_L3=%.2f%n",
				udll.fgraEstimate(), udlc[0], udlc[1], udlc[2]);
		}
	}
}
