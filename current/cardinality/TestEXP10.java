package cardinality;

import shared.Tools;

/**
 * Stress-tests UltraDynamicLogLog6 UDLC estimation at cardinality=8000
 * across 1000 instances, reporting average L1/L2 estimates and shrinkage.
 *
 * @author Brian Bushnell
 */
public class TestEXP10 {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		final int buckets=2048;
		final int numInstances=1000;
		final int testCard=8000;

		double l1Sum=0, l2Sum=0, avgSum=0;
		int count=0;

		for(int inst=0; inst<numInstances; inst++){
			final UltraDynamicLogLog6 udll=new UltraDynamicLogLog6(buckets, 31, -1, 0);
			long seed=inst*12345678901L+42;
			for(int card=1; card<=testCard; card++){
				seed=Tools.hash64shift(seed);
				udll.hashAndStore(seed);
			}
			final double[] udlc=udll.udlcEstimate();
			l1Sum+=udlc[0];
			l2Sum+=udlc[1];
			avgSum+=udlc[3];
			count++;
		}

		System.out.printf("After %d instances at card=%d:%n", count, testCard);
		System.out.printf("Avg L1 = %.2f%n", l1Sum/count);
		System.out.printf("Avg L2 = %.2f%n", l2Sum/count);
		System.out.printf("Avg result = %.2f%n", avgSum/count);
		System.out.printf("Ratio L1/L2 = %.4f%n", (l1Sum/count)/(l2Sum/count));

		// Manual calculation
		final double dlc1=l1Sum/count;
		final double dlc2=l2Sum/count;
		double expectedResult=dlc1;
		if(dlc2>0 && dlc1>0){
			final double ratio=dlc1/dlc2;
			final double shrinkage=Math.min(2.0, Math.max(0.5, 2.0/ratio));
			expectedResult=dlc1*Math.pow(shrinkage, 0.3);
			System.out.printf("Shrinkage factor = %.4f, expected result = %.2f%n", shrinkage, expectedResult);
		}
	}

}
