package cardinality;

import java.util.Random;

/**
 * Speed and branch-rate comparison across estimator types.
 */
public class SpeedTest {
    public static void main(String[] args){
        final int buckets=512;
        final int k=31;
        final long seed=12345L;
        final int numEstimators=1000;
        final long maxCard=buckets*1024L; // 524288

        System.out.println("Type\tTime_ms\tBranch1Rate\tBranch2Rate\tFinalCard");

        // DLL4
        runTest("DLL4", numEstimators, buckets, k, seed, maxCard, "dll4");
        // LL6
        runTest("LL6", numEstimators, buckets, k, seed, maxCard, "ll6");
        // DDL
        runTest("DDL", numEstimators, buckets, k, seed, maxCard, "ddl");
        // DDL8
        runTest("DDL8", numEstimators, buckets, k, seed, maxCard, "ddl8");
        // ErtlULLb
        runTest("ErtlULLb", numEstimators, buckets, k, seed, maxCard, "ertlb");
        // ULLc
        runTest("ULLc", numEstimators, buckets, k, seed, maxCard, "ullc");
    }

    static void runTest(String label, int numEst, int buckets, int k, long seed, long maxCard, String type){
        long t0=System.currentTimeMillis();
        double b1sum=0, b2sum=0;
        long lastCard=0;
        for(int e=0; e<numEst; e++){
            CardinalityTracker ct=DDLCalibrationDriver.makeInstance(type, buckets, k, seed+e, 0);
            Random rng=new Random(seed+e+1000000);
            for(long c=0; c<maxCard; c++){
                ct.add(rng.nextLong());
            }
            ct.lastCardinality=-1;
            lastCard=ct.cardinality();
            // Extract branch rates if available
            if(ct instanceof DynamicLogLog4){
                DynamicLogLog4 d=(DynamicLogLog4)ct;
                b1sum+=d.branch1Rate(); b2sum+=d.branch2Rate();
            }else if(ct instanceof DynamicDemiLog){
                DynamicDemiLog d=(DynamicDemiLog)ct;
                b1sum+=d.branch1Rate(); b2sum+=d.branch2Rate();
            }else if(ct instanceof DynamicDemiLog8){
                DynamicDemiLog8 d=(DynamicDemiLog8)ct;
                b1sum+=d.branch1Rate(); b2sum+=d.branch2Rate();
            }else if(ct instanceof ULLc){
                ULLc d=(ULLc)ct;
                b1sum+=d.branch1Rate(); b2sum+=d.branch2Rate();
            }
        }
        long t1=System.currentTimeMillis();
        double b1avg=b1sum/numEst;
        double b2avg=b2sum/numEst;
        System.out.printf("%s\t%d\t%.6f\t%.6f\t%d%n", label, (t1-t0), b1avg, b2avg, lastCard);
    }
}
