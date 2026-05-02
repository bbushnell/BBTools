package cardinality;

import shared.Tools;

public class TTLLDiag {

	public static void main(String[] args){
		CardinalityParser.initializeAll("ttll", 2048, 31, null, null);

		// Test multiple seeds to find one that produces extreme Mean+H at cardinality=1
		System.err.println("--- Scanning seeds for extreme Mean+H at cardinality=1 ---");
		double worstErr=0;
		long worstSeed=0;
		for(long seed=0; seed<10000; seed++){
			TwinTailLogLog t=new TwinTailLogLog(2048, 31, seed, 0);
			t.add(1);
			double[] ldlc=t.ldlcEstimate();
			double meanH=ldlc[6];
			double err=Math.abs(meanH-1.0);
			if(err>worstErr){
				worstErr=err;
				worstSeed=seed;
			}
		}
		System.err.println("Worst seed: "+worstSeed+" (Mean+H error="+worstErr+")");
		System.err.println();

		// Now trace the worst seed in detail
		long[] seeds={12345, worstSeed};
		for(long seed : seeds){
			TwinTailLogLog t=new TwinTailLogLog(2048, 31, seed, 0);
			t.add(1);

			System.err.println("=== Seed "+seed+", cardinality=1 ===");

			// Find the filled bucket
			// We can't access regs directly, but we can re-hash to find what bucket it went to
			long key=Tools.hash64shift(1L^t.hashXor);
			int bucket=(int)(key&(2048-1));
			int hashNLZ=Long.numberOfLeadingZeros(key);
			System.err.println("  hash key="+Long.toHexString(key));
			System.err.println("  bucket="+bucket+" hashNLZ="+hashNLZ);

			double[] ldlc=t.ldlcEstimate();
			System.err.println("  Mean+H = "+ldlc[6]);
			System.err.println("  LDLC   = "+ldlc[0]);
			System.err.println();
		}
	}
}
