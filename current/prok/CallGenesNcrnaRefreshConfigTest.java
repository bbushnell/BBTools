package prok;

import map.LongHashSet;

/** Focused parser-to-scavenger plumbing checks for the production claimed-window refresh flag. */
public class CallGenesNcrnaRefreshConfigTest {

	public static void main(String[] args){
		final boolean original=CallGenes.NCRNA_REFRESH_CLAIMED_WINDOWS;
		try{
			if(original){throw new AssertionError("v2refresh must default to false to preserve master's claimed-window policy");}
			if(CallGenes.parseNcrnaRefreshFlag("notv2refresh", "f")){
				throw new AssertionError("Unknown flag was consumed by v2refresh parser");
			}
			if(CallGenes.NCRNA_REFRESH_CLAIMED_WINDOWS){throw new AssertionError("Unknown flag changed the default");}

			final NcrnaScavenger scavenger=new NcrnaScavenger(new byte[][]{{'A','C','G','T','A','C','G','T'}}, null, null,
				new LongHashSet(), 17, 1, 0);
			if(!CallGenes.parseNcrnaRefreshFlag("v2refresh", "f")){
				throw new AssertionError("v2refresh=f was not parsed");
			}
			GeneCaller.applyNcrnaProductionControls(scavenger);
			if(scavenger.refreshClaimedWindows){throw new AssertionError("v2refresh=f did not reach NcrnaScavenger");}

			if(!CallGenes.parseNcrnaRefreshFlag("V2REFRESH", "t")){
				throw new AssertionError("v2refresh=t was not parsed case-insensitively");
			}
			GeneCaller.applyNcrnaProductionControls(scavenger);
			if(!scavenger.refreshClaimedWindows){throw new AssertionError("v2refresh=t did not reach NcrnaScavenger");}
		}finally{
			CallGenes.NCRNA_REFRESH_CLAIMED_WINDOWS=original;
		}
		System.out.println("PASS CallGenesNcrnaRefreshConfigTest");
	}
}
