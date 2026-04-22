package cardinality;

import java.util.Random;

/**
 * Speed and branch-rate comparison across estimator types.
 * Runs many estimator instances to steady state and reports elapsed time,
 * average branch rates, and final cardinality estimate.
 *
 * @author Brian Bushnell
 */
public class SpeedTest {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		final int buckets=512;
		final int k=31;
		final long seed=12345L;
		final int numEstimators=1000;
		final long maxCard=buckets*1024L; // 524288

		System.out.println("Type\tTime_ms\tBranch1Rate\tBranch2Rate\tFinalCard");

		runTest("DLL4", numEstimators, buckets, k, seed, maxCard, "dll4");
		runTest("LL6", numEstimators, buckets, k, seed, maxCard, "ll6");
		runTest("DDL", numEstimators, buckets, k, seed, maxCard, "ddl");
		runTest("DDL8", numEstimators, buckets, k, seed, maxCard, "ddl8");
		runTest("ErtlULLb", numEstimators, buckets, k, seed, maxCard, "ertlb");
		runTest("ULLc", numEstimators, buckets, k, seed, maxCard, "ullc");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Runs a single estimator-type benchmark.
	 * Creates numEst instances, feeds maxCard random longs into each,
	 * and accumulates branch rates for types that expose them.
	 */
	static void runTest(String label, int numEst, int buckets, int k, long seed, long maxCard, String type){
		final long t0=System.currentTimeMillis();
		double b1sum=0, b2sum=0;
		long lastCard=0;
		for(int e=0; e<numEst; e++){
			final CardinalityTracker ct=DDLCalibrationDriver.makeInstance(type, buckets, k, seed+e, 0);
			final Random rng=new Random(seed+e+1000000);
			for(long c=0; c<maxCard; c++){
				ct.add(rng.nextLong());
			}
			ct.lastCardinality=-1;
			lastCard=ct.cardinality();
			if(ct instanceof DynamicLogLog4){
				final DynamicLogLog4 d=(DynamicLogLog4)ct;
				b1sum+=d.branch1Rate(); b2sum+=d.branch2Rate();
			}else if(ct instanceof DynamicDemiLog){
				final DynamicDemiLog d=(DynamicDemiLog)ct;
				b1sum+=d.branch1Rate(); b2sum+=d.branch2Rate();
			}else if(ct instanceof DynamicDemiLog8){
				final DynamicDemiLog8 d=(DynamicDemiLog8)ct;
				b1sum+=d.branch1Rate(); b2sum+=d.branch2Rate();
			}
		}
		final long t1=System.currentTimeMillis();
		final double b1avg=b1sum/numEst;
		final double b2avg=b2sum/numEst;
		System.out.printf("%s\t%d\t%.6f\t%.6f\t%d%n", label, (t1-t0), b1avg, b2avg, lastCard);
	}

}
