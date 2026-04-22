package cardinality;

import shared.Tools;

/**
 * Standalone HC diagnostic for CDLL4: compute per-tier HC estimates.
 * Feeds random hashes into a CompressedDynamicLogLog4, then reads bucket
 * nibbles to build per-tier history counts and DLC reference estimates.
 */
public class TestCDLL4HC {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		final int B=512;
		CompressedDynamicLogLog4.DEBUG_HIST=false;
		final CompressedDynamicLogLog4 c=new CompressedDynamicLogLog4(B, 31, 12345L, 0);
		final java.util.Random r=new java.util.Random(42);
		final int trueCard=200000;
		for(int i=0; i<trueCard; i++){c.add(r.nextLong());}

		// Extract bucket data: absTier and hist for each bucket
		// Force summarize to get floor level (via getMinZeros compat accessor)
		c.rawEstimates();
		final int minZeros=c.getMinZeros();

		// Read nibbles directly (we're in the same package)
		final int[] nlzBucketCount=new int[64];
		final int[][] nlzHbitSet=new int[1][64]; // 1-bit history
		final int[] histAtTier=new int[64]; // count of hist=1 per tier
		final int[] totalAtTier=new int[64];

		for(int i=0; i<B; i++){
			final int nib=c.readNibble(i);
			final int tp=nib>>>1;
			final int hist=nib&1;
			int absTier;
			if(tp>0){
				absTier=(tp-1)+minZeros;
			}else if(minZeros>0){
				absTier=minZeros-1; // floor-level bucket
			}else{
				continue; // empty
			}
			if(absTier>=0 && absTier<64){
				nlzBucketCount[absTier]++;
				totalAtTier[absTier]++;
				if(hist==1){nlzHbitSet[0][absTier]++;}
				histAtTier[absTier]+=hist;
			}
		}

		System.err.println("trueCard="+trueCard+" minZeros="+minZeros+" B="+B);
		System.err.println();
		System.err.println("Tier  Buckets  Hist=1  Unseen  HC_est       Ratio");

		// HC: for each tier t, look at tier t+1
		int maxTier=0;
		for(int t=63; t>=0; t--){if(nlzBucketCount[t]>0){maxTier=t; break;}}

		for(int t=0; t<=maxTier+1; t++){
			final int sourceTier=t+1;
			if(sourceTier>=64){continue;}
			final int beff=nlzBucketCount[sourceTier];
			final int unseen=beff-nlzHbitSet[0][sourceTier];
			if(beff<2 || unseen<1 || unseen>=beff){continue;}

			final double est=Math.pow(2.0, (t+1)*TIER_SCALE)*B*Math.log((double)beff/unseen);
			final double ratio=est/trueCard;

			System.err.printf("%-5d %-8d %-7d %-7d %-12.0f %.3f%n",
				t, beff, nlzHbitSet[0][sourceTier], unseen, est, ratio);
		}

		// Also compute DLC estimate for reference
		System.err.println();
		System.err.println("DLC tiers for reference:");
		int vk=0;
		for(int k=0; k<64; k++){
			vk+=nlzBucketCount[k]; // floor-level + tier k buckets
			if(vk>0 && vk<B){
				final double dlcEst=Math.pow(2.0, k*TIER_SCALE)*B*Math.log((double)B/vk);
				System.err.printf("  DLC tier %d: vk=%d est=%.0f ratio=%.3f%n",
					k, vk, dlcEst, dlcEst/trueCard);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	static final int MANTISSA_THRESHOLD=(int)Math.round((2.0-Math.sqrt(2.0))*65536);
	static final double TIER_SCALE=1.5;
}
