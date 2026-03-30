package cardinality;

import shared.Tools;

/**
 * Quick test: UDLC layer-peeling estimates vs FGRA at various cardinalities.
 * Prints TrueCard, FGRA, DLC_L1, DLC_L2, DLC_L3, UDLC_avg for each checkpoint.
 */
public class UDLCTest {

	public static void main(String[] args){
		final int buckets=2048;
		final int numInstances=1000;
		final int maxCard=buckets*16; // 32768, well past the blend zone

		System.out.println("TrueCard\tFGRA\tDLC_L1\tDLC_L2\tDLC_L3\tUDLC_avg\tFGRA_err\tL1_err\tUDLC_err");

		// Checkpoints at ~10% intervals
		int nextCheck=1;
		final double[] fgraSum=new double[maxCard+1];
		final double[] l1Sum=new double[maxCard+1];
		final double[] l2Sum=new double[maxCard+1];
		final double[] l3Sum=new double[maxCard+1];
		final double[] avgSum=new double[maxCard+1];
		final int[] counts=new int[maxCard+1];

		for(int inst=0; inst<numInstances; inst++){
			UltraDynamicLogLog6 udll=new UltraDynamicLogLog6(buckets, 31, -1, 0);
			long seed=inst*12345678901L+42;
			for(int card=1; card<=maxCard; card++){
				seed=Tools.hash64shift(seed);
				udll.hashAndStore(seed);
				// Report at checkpoints
				if(card==1 || card==2 || card==5 || card==10 || card==20 || card==50 ||
				   card==100 || card==200 || card==500 || card==1000 || card==2000 ||
				   card==4000 || card==8000 || card==16000 || card==maxCard){
					final double fgra=udll.fgraEstimate();
					final double[] udlc=udll.udlcEstimate();
					fgraSum[card]+=fgra;
					l1Sum[card]+=udlc[0];
					l2Sum[card]+=udlc[1];
					l3Sum[card]+=udlc[2];
					avgSum[card]+=udlc[3];
					counts[card]++;
				}
			}
		}

		// Print averages
		for(int card=1; card<=maxCard; card++){
			if(counts[card]==0){continue;}
			final double n=counts[card];
			final double fgra=fgraSum[card]/n;
			final double l1=l1Sum[card]/n;
			final double l2=l2Sum[card]/n;
			final double l3=l3Sum[card]/n;
			final double avg=avgSum[card]/n;
			final double fgraErr=Math.abs(fgra-card)/card;
			final double l1Err=Math.abs(l1-card)/card;
			final double avgErr=Math.abs(avg-card)/card;
			System.out.printf("%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.5f\t%.5f\t%.5f%n",
				card, fgra, l1, l2, l3, avg, fgraErr, l1Err, avgErr);
		}
	}
}
