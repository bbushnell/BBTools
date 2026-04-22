package cardinality;

import rand.FastRandomXoshiro;

/**
 * Tests DynamicLogLog2 DLCBest estimator at cardinality=1
 * across 10 seeds, comparing LC, LCmin, DLCBest, and DLC estimates.
 *
 * @author Brian Bushnell
 */
public class TestDLCBest1 {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		DynamicLogLog2.IGNORE_OVERFLOW=true;
		for(int seed=0; seed<10; seed++){
			final DynamicLogLog2 dll=new DynamicLogLog2(2048, 31, seed, 0);
			dll.add(new FastRandomXoshiro(seed).nextLong());
			final double[] est=dll.rawEstimates();
			// est[13]=DLCBest, est[5]=LC, est[9]=LCmin
			System.out.printf("seed=%d  LC=%.4f  LCmin=%.4f  DLCBest=%.4f  DLC=%.4f%n",
				seed, est[5], est[9], est[13], est[11]);
		}
	}

}
