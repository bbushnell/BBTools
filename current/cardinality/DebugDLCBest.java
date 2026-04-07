package cardinality;

import rand.FastRandomXoshiro;

/**
 * Debug tool: creates one DLL2 with io=t, adds elements to ~600K,
 * then prints full NLZ state and DLCBest calculation step by step.
 */
public class DebugDLCBest {
	public static void main(String[] args){
		DynamicLogLog2.IGNORE_OVERFLOW=true;
		long targetCard=600000;
		for(String arg : args){
			String[] s=arg.split("=");
			if(s[0].equals("card")){targetCard=Long.parseLong(s[1]);}
		}

		DynamicLogLog2 dll=new DynamicLogLog2(2048, 31, 42, 0);
		FastRandomXoshiro rng=new FastRandomXoshiro(42);
		for(long i=0; i<targetCard; i++){
			dll.add(rng.nextLong());
		}

		// Force summarize to populate nlzCounts
		dll.rawEstimates();
		int[] nlz=dll.nlzCounts;

		System.out.println("=== DLL2 Debug at trueCard="+targetCard+" ===");
		System.out.println("B=2048, IGNORE_OVERFLOW="+DynamicLogLog2.IGNORE_OVERFLOW);
		System.out.println();

		// Print NLZ histogram
		int empties=nlz[0];
		System.out.println("Empties (nlz[0]): "+empties);
		int total=empties;
		int lowestActive=-1;
		int highestActive=-1;
		for(int k=0; k<65; k++){
			if(nlz[k+1]>0){
				System.out.println("  absNlz="+k+": "+nlz[k+1]+" buckets");
				total+=nlz[k+1];
				if(lowestActive<0) lowestActive=k;
				highestActive=k;
			}
		}
		System.out.println("Total buckets: "+total);
		System.out.println("Lowest active tier (=minZeros): "+lowestActive);
		System.out.println("Highest active tier: "+highestActive);
		System.out.println();

		// Now trace DLCBest calculation
		int B=2048;
		int V=empties;
		double target=B*0.25; // =512
		System.out.println("DLCBest target V_k = "+target);
		System.out.println();

		int startTier=lowestActive;
		int vk=V;
		int bestTier=startTier;
		double bestDist=Math.abs(vk-target);

		System.out.println("Tier-by-tier V_k and occupancy:");
		System.out.printf("  tier=%d: V_k=%d, occ=%d, dist=%.1f%n", startTier, vk, B-vk, bestDist);

		int[] vks=new int[66];
		vks[startTier]=vk;
		for(int tier=startTier+1; tier<=highestActive+1 && tier<66; tier++){
			if(tier<nlz.length) vk+=nlz[tier]; // nlz[tier] = buckets at absNlz=tier-1
			vks[tier]=vk;
			double dist=Math.abs(vk-target);
			String marker="";
			if(dist<bestDist && vk>=1 && vk<B){
				bestDist=dist;
				bestTier=tier;
				marker=" <-- new best";
			}
			System.out.printf("  tier=%d: V_k=%d, occ=%d, dist=%.1f%s%n",
				tier, vk, B-vk, dist, marker);
			if(vk>=B) break;
		}

		System.out.println();
		System.out.println("Selected bestTier="+bestTier+", V_k="+vks[bestTier]+", occ="+(B-vks[bestTier]));

		// Now show what dlcFromVk does
		int bestVk=vks[bestTier];
		double rawEst=(1L<<bestTier)*(double)B*Math.log((double)B/Math.max(bestVk, 0.5));
		System.out.println();
		System.out.println("Raw LC estimate: 2^"+bestTier+" * "+B+" * ln("+B+"/"+bestVk+") = "+rawEst);
		System.out.println("Signed error (raw): "+(rawEst-targetCard)/(double)targetCard);

		// Now check correction
		int nativeRelTiers=3;
		int topStoredTier=startTier+nativeRelTiers-1;
		System.out.println();
		System.out.println("nativeRelTiers="+nativeRelTiers);
		System.out.println("topStoredTier = startTier+2 = "+topStoredTier);
		System.out.println("bestTier="+bestTier+" == topStoredTier? "+(bestTier==topStoredTier));

		// The correction check in dlcFromVk:
		// cfTier = topStoredTier (only correct when tier == topStoredTier)
		// if(tier < overflowTier) return raw  --> tier < topStoredTier means NO correction
		// So correction fires when tier >= topStoredTier, i.e., tier == topStoredTier
		boolean wouldCorrect = (bestTier >= topStoredTier);
		System.out.println("Would correction fire? "+wouldCorrect);

		if(wouldCorrect){
			double x=Math.log((double)B/bestVk)/Math.log(2);
			double err=CorrectionFactor.dlcTierSignedErr(x, CorrectionFactor.DLC_TIER_ERR_DLL2);
			double cf=1.0/(1.0+err);
			double corrected=rawEst*cf;
			int occ=B-bestVk;
			System.out.println("x = log2(B/V_k) = "+x);
			System.out.println("Polynomial signed err = "+err);
			System.out.println("CF = 1/(1+err) = "+cf);
			System.out.println("Corrected estimate = "+corrected);
			System.out.println("Corrected signed error = "+(corrected-targetCard)/(double)targetCard);
			System.out.println("Occupancy = "+occ+", min occ threshold = "+(int)(B*0.70));
			System.out.println("Passes occ guard? "+(occ >= B*0.70));
		}

		// Also show what CardStats actually computes
		System.out.println();
		System.out.println("=== Actual CardStats output ===");
		double[] est=dll.rawEstimates();
		// DLCBest is at index 13 in the legacy array
		CardStats s=new CardStats(null, nlz, 0, 0, 0, 0,
			B, 0, targetCard, null, 2048, 0,
			3, CorrectionFactor.DLC_TIER_ERR_DLL2);
		System.out.println("dlcBestF = "+s.dlcBest());
		System.out.println("dlcBest signed err = "+(s.dlcBest()-targetCard)/(double)targetCard);
		System.out.println("lcMinF = "+s.lcMin());
		System.out.println("dlcRawF = "+s.dlcRaw());
	}
}
