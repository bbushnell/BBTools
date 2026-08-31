package prok;

import java.util.Arrays;

/** Focused checks for the experimental tmRNA ncRNA-family gate and sweep aliases. */
public class CallGenesTmrnaConfigTest {

	public static void main(String[] args){
		testShippedOff();
		check("tmrna", CallGenes.parseNcrnaFamily("tmrna"));
		check("tmrna", CallGenes.parseNcrnaFamily("tm-rna"));
		check("tmrna", CallGenes.parseNcrnaFamily("ssra"));
		checkFloat(0.62f, CallGenes.defaultNcrnaIdPass("tmrna"));
		checkFloat(0.60f, CallGenes.defaultNcrnaIdBorderline("tmrna"));
		testSubordinateGate();
		testSweepRequiresTmrnaGate();
		testSpecificOverridesRequireTmrnaGate();
		testCompressedKmerResource();
		testDefaultBundle();
		testDefaultBoundaryBundle();
		System.out.println("PASS CallGenesTmrnaConfigTest");
	}

	private static void testShippedOff(){
		if(CallGenes.TMRNA_ENABLED){throw new AssertionError("tmRNA caller is unexpectedly enabled by default");}
	}

	private static void testDefaultBundle(){
		final boolean oldNcrna=CallGenes.NCRNA_FAMILIES_ENABLED;
		final boolean oldTmrna=CallGenes.TMRNA_ENABLED;
		try{
			GeneCaller.ncrnaFamilies.clear();
			CallGenes.NCRNA_FAMILIES_ENABLED=true;
			CallGenes.TMRNA_ENABLED=true;
			CallGenes.loadNcrnaResources();
			NcrnaFamily observed=null;
			for(NcrnaFamily family : GeneCaller.ncrnaFamilies){
				if(family.name.equals("tmrna")){observed=family; break;}
			}
			if(observed==null){throw new AssertionError("Default tmRNA bundle was not loaded");}
			checkFloat(5f, observed.scoreA);
			checkFloat(3f, observed.scoreB);
			checkFloat(0.62f, observed.idPass);
			checkFloat(0.60f, observed.idBorderline);
			checkFloat(0.62f, observed.hbmPass);
			if(observed.windowPad!=370 || observed.library.length!=1128 || observed.kmerSet.size()!=1011){
				throw new AssertionError("Unexpected default bundle: pad="+observed.windowPad
					+" models="+observed.library.length+" kmers="+observed.kmerSet.size());
			}
		}finally{
			GeneCaller.ncrnaFamilies.clear();
			CallGenes.NCRNA_FAMILIES_ENABLED=oldNcrna;
			CallGenes.TMRNA_ENABLED=oldTmrna;
		}
	}

	private static void testDefaultBoundaryBundle(){
		final boolean oldNcrna=CallGenes.NCRNA_FAMILIES_ENABLED;
		final boolean oldTmrna=CallGenes.TMRNA_ENABLED;
		final boolean oldBoundary=CallGenes.NCRNA_BOUNDARY_NN_ENABLED;
		try{
			GeneCaller.ncrnaFamilies.clear();
			CallGenes.NCRNA_FAMILIES_ENABLED=true;
			CallGenes.TMRNA_ENABLED=true;
			CallGenes.NCRNA_BOUNDARY_NN_ENABLED=true;
			CallGenes.loadNcrnaResources();
			NcrnaFamily observed=null;
			for(NcrnaFamily family : GeneCaller.ncrnaFamilies){
				if(family.name.equals("tmrna")){observed=family; break;}
			}
			if(observed==null || observed.boundary5NetTemplate==null || observed.boundary3NetTemplate==null
					|| observed.boundaryStartTable==null || observed.boundaryStopTable==null){
				throw new AssertionError("Default tmRNA boundary bundle was not loaded completely");
			}
			if(observed.boundaryStartInside!=7 || observed.boundaryStartOutside!=0
					|| observed.boundaryStopInside!=11 || observed.boundaryStopOutside!=0
					|| observed.boundaryMeanLen<=0f){
				throw new AssertionError("Unexpected tmRNA boundary metadata");
			}
			if(!Arrays.equals(new int[]{-3,-2,-1,0,1,2}, observed.boundaryStartOffsets)
					|| !Arrays.equals(new int[]{-3,-2,-1,0,1,2}, observed.boundaryStopOffsets)){
				throw new AssertionError("Unexpected tmRNA boundary offsets");
			}
		}finally{
			GeneCaller.ncrnaFamilies.clear();
			CallGenes.NCRNA_FAMILIES_ENABLED=oldNcrna;
			CallGenes.TMRNA_ENABLED=oldTmrna;
			CallGenes.NCRNA_BOUNDARY_NN_ENABLED=oldBoundary;
		}
	}

	private static void testCompressedKmerResource(){
		final map.LongHashSet kmers=ProkObject.loadLongKmersByType(17, "tmrna");
		if(kmers==null || kmers.size()!=1011){
			throw new AssertionError("Expected 1011 unique tmRNA 17-mers from compressed resource, observed "
				+(kmers==null ? "null" : kmers.size()));
		}
	}

	private static void testSpecificOverridesRequireTmrnaGate(){
		final boolean oldNcrna=CallGenes.NCRNA_FAMILIES_ENABLED;
		final boolean oldTmrna=CallGenes.TMRNA_ENABLED;
		final String oldConsensus=CallGenes.TMRNA_CONSENSUS_OVERRIDE;
		final String oldKmers=CallGenes.TMRNA_KMERS_OVERRIDE;
		final int oldPad=CallGenes.TMRNA_PAD_OVERRIDE;
		final float oldScoreA=CallGenes.TMRNA_SCORE_A_OVERRIDE;
		try{
			CallGenes.NCRNA_FAMILIES_ENABLED=true;
			CallGenes.TMRNA_ENABLED=false;
			assertSpecificOverrideRejected("model", new Runnable(){
				@Override public void run(){CallGenes.TMRNA_CONSENSUS_OVERRIDE="test.fa";}
			});
			CallGenes.TMRNA_CONSENSUS_OVERRIDE=null;
			assertSpecificOverrideRejected("kmer", new Runnable(){
				@Override public void run(){CallGenes.TMRNA_KMERS_OVERRIDE="test.fa";}
			});
			CallGenes.TMRNA_KMERS_OVERRIDE=null;
			assertSpecificOverrideRejected("padding", new Runnable(){
				@Override public void run(){CallGenes.TMRNA_PAD_OVERRIDE=1;}
			});
			CallGenes.TMRNA_PAD_OVERRIDE=-1;
			assertSpecificOverrideRejected("score", new Runnable(){
				@Override public void run(){CallGenes.TMRNA_SCORE_A_OVERRIDE=1f;}
			});
			CallGenes.TMRNA_ENABLED=true;
			CallGenes.validateNcrnaSweepOverrides();
		}finally{
			CallGenes.NCRNA_FAMILIES_ENABLED=oldNcrna;
			CallGenes.TMRNA_ENABLED=oldTmrna;
			CallGenes.TMRNA_CONSENSUS_OVERRIDE=oldConsensus;
			CallGenes.TMRNA_KMERS_OVERRIDE=oldKmers;
			CallGenes.TMRNA_PAD_OVERRIDE=oldPad;
			CallGenes.TMRNA_SCORE_A_OVERRIDE=oldScoreA;
		}
	}

	private static void assertSpecificOverrideRejected(String label, Runnable setter){
		setter.run();
		boolean failed=false;
		try{CallGenes.validateNcrnaSweepOverrides();}
		catch(IllegalArgumentException expected){failed=true;}
		if(!failed){throw new AssertionError("tmRNA "+label+" override without tmrna=t was accepted");}
	}

	private static void testSubordinateGate(){
		final boolean oldNcrna=CallGenes.NCRNA_FAMILIES_ENABLED;
		final boolean oldTmrna=CallGenes.TMRNA_ENABLED;
		try{
			CallGenes.NCRNA_FAMILIES_ENABLED=false;
			CallGenes.TMRNA_ENABLED=true;
			boolean failed=false;
			try{CallGenes.validateNcrnaGateCombo();}
			catch(AssertionError expected){failed=true;}
			if(!failed){throw new AssertionError("tmrna=t without ncrna=t was accepted");}
			CallGenes.NCRNA_FAMILIES_ENABLED=true;
			CallGenes.validateNcrnaGateCombo();
		}finally{
			CallGenes.NCRNA_FAMILIES_ENABLED=oldNcrna;
			CallGenes.TMRNA_ENABLED=oldTmrna;
		}
	}

	private static void testSweepRequiresTmrnaGate(){
		final boolean oldNcrna=CallGenes.NCRNA_FAMILIES_ENABLED;
		final boolean oldTmrna=CallGenes.TMRNA_ENABLED;
		final String oldFamily=CallGenes.NCRNA_FAMILY_FILTER;
		try{
			CallGenes.NCRNA_FAMILIES_ENABLED=true;
			CallGenes.TMRNA_ENABLED=false;
			CallGenes.NCRNA_FAMILY_FILTER="tmrna";
			boolean failed=false;
			try{CallGenes.validateNcrnaSweepOverrides();}
			catch(IllegalArgumentException expected){failed=true;}
			if(!failed){throw new AssertionError("tmRNA sweep without tmrna=t was accepted");}
			CallGenes.TMRNA_ENABLED=true;
			CallGenes.validateNcrnaSweepOverrides();
		}finally{
			CallGenes.NCRNA_FAMILIES_ENABLED=oldNcrna;
			CallGenes.TMRNA_ENABLED=oldTmrna;
			CallGenes.NCRNA_FAMILY_FILTER=oldFamily;
		}
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
