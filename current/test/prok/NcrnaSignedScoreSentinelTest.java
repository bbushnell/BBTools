package test.prok;

import prok.CallGenes;

/** Signed score overrides are selected only by explicit flags, so -1 and 0 are data. */
public final class NcrnaSignedScoreSentinelTest {
	public static void main(String[] args){
		final String old=CallGenes.NCRNA_FAMILY_FILTER;
		try{
			CallGenes.NCRNA_FAMILY_FILTER="tmrna";
			check(CallGenes.resolveFamilyScore("tmrna", -1f, true, 7f, true, 5f)==-1f,
				"Explicit family score -1 must outrank generic and shipped scores");
			check(CallGenes.resolveFamilyScore("tmrna", 0f, true, 7f, true, 5f)==0f,
				"Explicit family score 0 must outrank generic and shipped scores");
			check(CallGenes.resolveFamilyScore("tmrna", 0f, false, -1f, true, 5f)==-1f,
				"Explicit generic score -1 must be a real signed score for a targeted family");
			check(CallGenes.resolveFamilyScore("tmrna", 0f, false, 0f, true, 5f)==0f,
				"Explicit generic score 0 must be a real signed score for a targeted family");
			check(CallGenes.resolveFamilyScore("tmrna", 0f, false, 7f, false, 5f)==5f,
				"Unset family and generic scores must inherit the registry default");
			check(CallGenes.resolveFamilyScore("sixs", 0f, false, 7f, true, 5f)==5f,
				"Generic sweep scores must not affect untargeted families");
			boolean failed=false;
			try{CallGenes.resolveFamilyScore("tmrna", Float.NaN, false, 7f, false, 5f);}
			catch(IllegalArgumentException expected){failed=true;}
			check(failed, "Nonfinite signed score carriers must fail even when their flag is unset");
		}finally{
			CallGenes.NCRNA_FAMILY_FILTER=old;
		}
		System.out.println("PASS NcrnaSignedScoreSentinelTest explicit_minus1 explicit_zero precedence");
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
