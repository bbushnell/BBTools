package cardinality;
import rand.FastRandomXoshiro;
public class TestCard1 {
    public static void main(String[] args){
        DynamicLogLog2.IGNORE_OVERFLOW=true;
        for(int seed=0; seed<20; seed++){
            DynamicLogLog2 dll=new DynamicLogLog2(2048, 31, seed, 0);
            FastRandomXoshiro rng=new FastRandomXoshiro(seed);
            long key=rng.nextLong();
            dll.add(key);
            dll.rawEstimates();
            int[] nlz=dll.nlzCounts;
            int V=nlz[0];
            int filled=2048-V;
            int nlzVal=Long.numberOfLeadingZeros(key);
            // Top stored tier
            int topCount=0, topTier=-1;
            for(int k=65; k>=0; k--){
                if(k+1<nlz.length && nlz[k+1]>0){topCount=nlz[k+1]; topTier=k; break;}
            }
            // LC compensation
            double missingFills=topCount*((double)V/2048);
            int lcV=Math.max(0, (int)Math.round(V-missingFills));
            double lcComp=2048.0*Math.log(2048.0/Math.max(lcV, 0.5));
            double lcPlain=2048.0*Math.log(2048.0/Math.max(V, 0.5));
            System.out.printf("seed=%2d NLZ=%2d filled=%d V=%d topTier=%d topCnt=%d mFills=%.2f lcV=%d LC_plain=%.3f LC_comp=%.3f err=%+.3f%n",
                seed, nlzVal, filled, V, topTier, topCount, missingFills, lcV, lcPlain, lcComp, (lcComp-1));
        }
    }
}
