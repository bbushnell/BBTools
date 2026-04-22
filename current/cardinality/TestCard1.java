package cardinality;

import rand.FastRandomXoshiro;

/**
 * Tests DynamicLogLog2 single-key insertion across 20 seeds,
 * printing NLZ distribution and linear-counting compensation values.
 *
 * @author Brian Bushnell
 */
public class TestCard1 {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		DynamicLogLog2.IGNORE_OVERFLOW=true;
		final int buckets=2048;
		for(int seed=0; seed<20; seed++){
			final DynamicLogLog2 dll=new DynamicLogLog2(buckets, 31, seed, 0);
			final FastRandomXoshiro rng=new FastRandomXoshiro(seed);
			final long key=rng.nextLong();
			dll.add(key);
			dll.rawEstimates();
			final int[] nlz=dll.nlzCounts;
			final int V=nlz[0];
			final int filled=buckets-V;
			final int nlzVal=Long.numberOfLeadingZeros(key);
			// Top stored tier
			int topCount=0, topTier=-1;
			for(int k=65; k>=0; k--){
				if(k+1<nlz.length && nlz[k+1]>0){topCount=nlz[k+1]; topTier=k; break;}
			}
			// LC compensation
			final double missingFills=topCount*((double)V/buckets);
			final int lcV=Math.max(0, (int)Math.round(V-missingFills));
			final double lcComp=buckets*Math.log(buckets/(double)Math.max(lcV, 0.5));
			final double lcPlain=buckets*Math.log(buckets/(double)Math.max(V, 0.5));
			System.out.printf("seed=%2d NLZ=%2d filled=%d V=%d topTier=%d topCnt=%d mFills=%.2f lcV=%d LC_plain=%.3f LC_comp=%.3f err=%+.3f%n",
				seed, nlzVal, filled, V, topTier, topCount, missingFills, lcV, lcPlain, lcComp, (lcComp-1));
		}
	}

}
