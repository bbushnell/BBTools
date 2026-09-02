package prok;

import java.util.Arrays;

/** Focused configuration checks for the default-off paired 6S ncRNA bundle. */
public class CallGenesSixsConfigTest {

	public static void main(String[] args){
		if(CallGenes.SIXS_ENABLED){throw new AssertionError("sixS caller is unexpectedly enabled by default");}
		check("sixs_rf01685", CallGenes.parseNcrnaFamily("sixs_rf01685"));
		assertBareSixsAliasRejected("sixs");
		assertBareSixsAliasRejected("ssrs");
		checkFloat(0.70f, CallGenes.defaultNcrnaIdPass("sixs_rf00013"));
		checkFloat(0.60f, CallGenes.defaultNcrnaIdBorderline("sixs_rf00013"));
		checkFloat(0.70f, CallGenes.defaultNcrnaIdBorderline("sixs_rf01685"));
		checkOffsets("sixs_rf00013");
		checkOffsets("sixs_rf01685");
		testSubordinateGate();
		testSpecificOverridesRequireSixsGate();
		testNcrnaGateDoesNotRegisterSixsByDefault();
		// This is intentionally resource-backed: the sealed caller values live at the
		// addNcrnaFamily construction site, not only in the sweep-default helpers.
		testShippedBundle();
		testMissingRf01685KmersFailLoud();
		System.out.println("PASS CallGenesSixsConfigTest");
	}

	private static void testSubordinateGate(){
		final boolean oldNcrna=CallGenes.NCRNA_FAMILIES_ENABLED, oldSixs=CallGenes.SIXS_ENABLED;
		try{
			CallGenes.NCRNA_FAMILIES_ENABLED=false;
			CallGenes.SIXS_ENABLED=true;
			try{CallGenes.validateNcrnaGateCombo(); throw new AssertionError("sixs=t without ncrna=t was accepted");}
			catch(AssertionError expected){/* pass */}
			CallGenes.NCRNA_FAMILIES_ENABLED=true;
			CallGenes.validateNcrnaGateCombo();
		}finally{CallGenes.NCRNA_FAMILIES_ENABLED=oldNcrna; CallGenes.SIXS_ENABLED=oldSixs;}
	}

	private static void testSpecificOverridesRequireSixsGate(){
		final boolean oldNcrna=CallGenes.NCRNA_FAMILIES_ENABLED, oldSixs=CallGenes.SIXS_ENABLED;
		final String oldKmers=CallGenes.SIXS_RF00013_KMERS_OVERRIDE;
		final int oldPad=CallGenes.SIXS_RF01685_PAD_OVERRIDE;
		try{
			CallGenes.NCRNA_FAMILIES_ENABLED=true;
			CallGenes.SIXS_ENABLED=false;
			CallGenes.SIXS_RF00013_KMERS_OVERRIDE="test.fa";
			assertOverrideRejected();
			CallGenes.SIXS_RF00013_KMERS_OVERRIDE=null;
			CallGenes.SIXS_RF01685_PAD_OVERRIDE=1;
			assertOverrideRejected();
			CallGenes.SIXS_ENABLED=true;
			CallGenes.validateNcrnaSweepOverrides();
		}finally{
			CallGenes.NCRNA_FAMILIES_ENABLED=oldNcrna; CallGenes.SIXS_ENABLED=oldSixs;
			CallGenes.SIXS_RF00013_KMERS_OVERRIDE=oldKmers; CallGenes.SIXS_RF01685_PAD_OVERRIDE=oldPad;
		}
	}

	private static void testShippedBundle(){
		final boolean oldNcrna=CallGenes.NCRNA_FAMILIES_ENABLED, oldSixs=CallGenes.SIXS_ENABLED;
		try{
			GeneCaller.ncrnaFamilies.clear();
			CallGenes.NCRNA_FAMILIES_ENABLED=true;
			CallGenes.SIXS_ENABLED=true;
			CallGenes.loadNcrnaResources();
			if(countSixs()!=2){throw new AssertionError("sixs=t did not register exactly its paired two-family bundle");}
			checkBundle(find("sixs_rf00013"), 16, 250, 590, 4000, 0.70f, 0.60f, 0.64f);
			checkBundle(find("sixs_rf01685"), 17, 125, 19, 437, 0.70f, 0.70f, 0.75f);
		}finally{
			GeneCaller.ncrnaFamilies.clear();
			CallGenes.NCRNA_FAMILIES_ENABLED=oldNcrna; CallGenes.SIXS_ENABLED=oldSixs;
		}
	}

	/** The sub-gate is load-bearing: ordinary ncrna=t must not promote unfinished sixS. */
	private static void testNcrnaGateDoesNotRegisterSixsByDefault(){
		final boolean oldNcrna=CallGenes.NCRNA_FAMILIES_ENABLED, oldSixs=CallGenes.SIXS_ENABLED;
		try{
			GeneCaller.ncrnaFamilies.clear();
			CallGenes.NCRNA_FAMILIES_ENABLED=true;
			CallGenes.SIXS_ENABLED=false;
			CallGenes.loadNcrnaResources();
			if(countSixs()!=0){throw new AssertionError("ncrna=t registered sixS despite sixs=f");}
		}finally{
			GeneCaller.ncrnaFamilies.clear();
			CallGenes.NCRNA_FAMILIES_ENABLED=oldNcrna; CallGenes.SIXS_ENABLED=oldSixs;
		}
	}

	/** Injects an absent path through the per-family seam; real resources stay untouched. */
	private static void testMissingRf01685KmersFailLoud(){
		final boolean oldNcrna=CallGenes.NCRNA_FAMILIES_ENABLED, oldSixs=CallGenes.SIXS_ENABLED;
		final String oldKmers=CallGenes.SIXS_RF01685_KMERS_OVERRIDE;
		final String missing=new java.io.File(System.getProperty("java.io.tmpdir"),
				"CallGenesSixsConfigTest_missing_rf01685_"+System.nanoTime()+".fa").getAbsolutePath();
		try{
			GeneCaller.ncrnaFamilies.clear();
			CallGenes.NCRNA_FAMILIES_ENABLED=true;
			CallGenes.SIXS_ENABLED=true;
			CallGenes.SIXS_RF01685_KMERS_OVERRIDE=missing;
			try{
				CallGenes.loadNcrnaResources();
				throw new AssertionError("Missing RF01685 kmer override did not fail loud");
			}catch(IllegalArgumentException expected){
				if(expected.getMessage()==null || !expected.getMessage().contains("sixs_rf01685")
						|| !expected.getMessage().contains(missing)){
					throw new AssertionError("Missing RF01685 kmer failure did not identify the subtype and path", expected);
				}
			}
		}finally{
			GeneCaller.ncrnaFamilies.clear();
			CallGenes.NCRNA_FAMILIES_ENABLED=oldNcrna; CallGenes.SIXS_ENABLED=oldSixs;
			CallGenes.SIXS_RF01685_KMERS_OVERRIDE=oldKmers;
		}
	}

	private static NcrnaFamily find(String name){
		for(NcrnaFamily f : GeneCaller.ncrnaFamilies){if(name.equals(f.name)){return f;}}
		throw new AssertionError("Missing expected sixS family "+name);
	}

	private static int countSixs(){
		int count=0;
		for(NcrnaFamily f : GeneCaller.ncrnaFamilies){if(f.name.startsWith("sixs_")){count++;}}
		return count;
	}

	private static void checkBundle(NcrnaFamily f, int kLong, int pad, int models, int kmers,
			float idPass, float idBorderline, float hbmPass){
		if(f.kLong!=kLong || f.windowPad!=pad || f.library.length!=models || f.models==null
				|| f.models.length!=f.library.length || f.modelNames==null || f.modelNames.length!=f.library.length || f.kmerSet.size()!=kmers
				|| f.minLen!=60 || f.indexK!=7 || f.indexTopN!=60 || f.adaptive || f.fixedMinHits!=1){
			throw new AssertionError("Unexpected sealed sixS bundle shape for "+f.name);
		}
		checkFloat(0f, f.scoreA); checkFloat(20f, f.scoreB); checkFloat(0.85f, f.collapseFrac);
		checkFloat(idPass, f.idPass); checkFloat(idBorderline, f.idBorderline); checkFloat(hbmPass, f.hbmPass);
	}

	private static void assertOverrideRejected(){
		try{CallGenes.validateNcrnaSweepOverrides(); throw new AssertionError("sixS override without sixs=t was accepted");}
		catch(IllegalArgumentException expected){/* pass */}
	}

	private static void checkOffsets(String family){
		final int[] expected={-3,-2,-1,0,1,2};
		if(!Arrays.equals(expected, CallGenes.boundaryStartOffsets(family))
				|| !Arrays.equals(expected, CallGenes.boundaryStopOffsets(family))){
			throw new AssertionError("Unexpected provisional 6S boundary offsets for "+family);
		}
	}

	private static void assertBareSixsAliasRejected(String value){
		try{CallGenes.parseNcrnaFamily(value); throw new AssertionError("Ambiguous ncrnafamily="+value+" was accepted");}
		catch(IllegalArgumentException expected){/* pass */}
	}

	private static void check(String expected, String observed){
		if(!expected.equals(observed)){throw new AssertionError("Expected "+expected+", observed "+observed);}
	}

	private static void checkFloat(float expected, float observed){
		if(Float.floatToIntBits(expected)!=Float.floatToIntBits(observed)){
			throw new AssertionError("Expected "+expected+", observed "+observed);
		}
	}
}
