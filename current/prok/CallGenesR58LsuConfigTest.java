package prok;

import java.util.ArrayList;

/** Focused checks for the packaged-resource, default-off 5.8S/LSU gate and explicit overrides. */
public class CallGenesR58LsuConfigTest {

	private static final String ROOT="/mnt/c/codex-lbl/Citan/workspace/ribo_new_members_20260901/"
		+"component_split_k2p5_target010_20260902/";
	private static final String R58_K=ROOT+"development_r58_17mers.fa";
	private static final String R58_C=ROOT+"split_baseline_r58_cid076/consensus.fa";
	private static final String R58_M=ROOT+"split_baseline_r58_cid076/models.hbm";
	private static final String LSU_K=ROOT+"development_lsu_17mers.fa";
	private static final String LSU_C=ROOT+"split_baseline_lsu_cid088/consensus.fa";
	private static final String LSU_M=ROOT+"split_baseline_lsu_cid088/models.hbm";

	public static void main(String[] args){
		if(CallGenes.R58LSU_ENABLED){throw new AssertionError("R58/LSU unexpectedly enabled by default");}
		if(!"r58".equals(CallGenes.parseNcrnaFamily("r58")) || !"lsu".equals(CallGenes.parseNcrnaFamily("lsu"))){
			throw new AssertionError("R58/LSU family names are not accepted");
		}
		if(CallGenes.defaultNcrnaIdPass("r58")!=0.60f || CallGenes.defaultNcrnaIdBorderline("lsu")!=0.55f){
			throw new AssertionError("Unexpected R58/LSU default thresholds");
		}
		testGateImplication();
		testBundledPair();
		testExplicitPair();
		testMissingLsuRollsBack();
		System.out.println("PASS CallGenesR58LsuConfigTest");
	}

	private static void testGateImplication(){
		final Snapshot s=new Snapshot();
		try{
			CallGenes.NCRNA_FAMILIES_ENABLED=false;
			CallGenes.R58LSU_ENABLED=false;
			CallGenes.setR58LsuEnabled(true);
			if(!CallGenes.R58LSU_ENABLED || !CallGenes.NCRNA_FAMILIES_ENABLED){
				throw new AssertionError("r58lsu=t did not enable the generic ncRNA machinery");
			}
		}finally{s.restore();}
	}

	private static void testBundledPair(){
		final Snapshot s=new Snapshot();
		try{
			GeneCaller.ncrnaFamilies.clear();
			CallGenes.NCRNA_FAMILIES_ENABLED=true; CallGenes.R58LSU_ENABLED=true;
			setPaths(null, null, null, null, null, null);
			CallGenes.loadNcrnaResources();
			if(count("r58")!=1 || count("lsu")!=1){throw new AssertionError("Bundled R58/LSU pair did not register exactly once");}
			check(find("r58"), 1, 500, 140, 135);
			check(find("lsu"), 14, 500, 60, 3500);
		}finally{s.restore();}
	}

	private static void testExplicitPair(){
		final Snapshot s=new Snapshot();
		try{
			GeneCaller.ncrnaFamilies.clear();
			CallGenes.NCRNA_FAMILIES_ENABLED=true; CallGenes.R58LSU_ENABLED=true;
			setPaths(R58_K, R58_C, R58_M, LSU_K, LSU_C, LSU_M);
			CallGenes.loadNcrnaResources();
			if(count("r58")!=1 || count("lsu")!=1){throw new AssertionError("Explicit R58/LSU pair did not register exactly once");}
			check(find("r58"), 2, 500, 140, 135);
			check(find("lsu"), 11, 1000, 60, 3500);
		}finally{s.restore();}
	}

	private static void testMissingLsuRollsBack(){
		final Snapshot s=new Snapshot();
		try{
			GeneCaller.ncrnaFamilies.clear();
			CallGenes.NCRNA_FAMILIES_ENABLED=true; CallGenes.R58LSU_ENABLED=true;
			setPaths(R58_K, R58_C, R58_M, "/tmp/no_such_lsu_kmers.fa", LSU_C, LSU_M);
			try{CallGenes.loadNcrnaResources(); throw new AssertionError("Missing LSU resource was accepted");}
			catch(IllegalArgumentException expected){/* pass */}
			if(findOrNull("r58")!=null || findOrNull("lsu")!=null){throw new AssertionError("Failed paired load left a partial family registration");}
		}finally{s.restore();}
	}

	private static void setPaths(String rk, String rc, String rm, String lk, String lc, String lm){
		CallGenes.R58_KMERS_OVERRIDE=rk; CallGenes.R58_CONSENSUS_OVERRIDE=rc; CallGenes.R58_MODELS_OVERRIDE=rm;
		CallGenes.LSU_KMERS_OVERRIDE=lk; CallGenes.LSU_CONSENSUS_OVERRIDE=lc; CallGenes.LSU_MODELS_OVERRIDE=lm;
	}

	private static NcrnaFamily find(String name){
		for(NcrnaFamily f : GeneCaller.ncrnaFamilies){if(name.equals(f.name)){return f;}}
		throw new AssertionError("Missing "+name);
	}

	private static NcrnaFamily findOrNull(String name){
		for(NcrnaFamily f : GeneCaller.ncrnaFamilies){if(name.equals(f.name)){return f;}}
		return null;
	}

	private static int count(String name){
		int count=0;
		for(NcrnaFamily f : GeneCaller.ncrnaFamilies){if(name.equals(f.name)){count++;}}
		return count;
	}

	private static void check(NcrnaFamily f, int models, int kmers, int minLen, int pad){
		if(f.library.length!=models || f.models.length!=models || f.modelNames.length!=models || f.kmerSet.size()!=kmers
				|| f.kLong!=17 || f.minLen!=minLen || f.windowPad!=pad || f.indexK!=7 || f.indexTopN!=100
				|| f.adaptive || f.fixedMinHits!=1 || f.idPass!=0.60f || f.idBorderline!=0.55f || f.hbmPass!=0.60f){
			throw new AssertionError("Unexpected "+f.name+" bundle shape");
		}
	}

	private static class Snapshot {
		final boolean ncrna=CallGenes.NCRNA_FAMILIES_ENABLED, r58lsu=CallGenes.R58LSU_ENABLED;
		final String rk=CallGenes.R58_KMERS_OVERRIDE, rc=CallGenes.R58_CONSENSUS_OVERRIDE, rm=CallGenes.R58_MODELS_OVERRIDE;
		final String lk=CallGenes.LSU_KMERS_OVERRIDE, lc=CallGenes.LSU_CONSENSUS_OVERRIDE, lm=CallGenes.LSU_MODELS_OVERRIDE;
		final ArrayList<NcrnaFamily> families=new ArrayList<NcrnaFamily>(GeneCaller.ncrnaFamilies);
		void restore(){
			GeneCaller.ncrnaFamilies.clear(); GeneCaller.ncrnaFamilies.addAll(families);
			CallGenes.NCRNA_FAMILIES_ENABLED=ncrna; CallGenes.R58LSU_ENABLED=r58lsu;
			CallGenes.R58_KMERS_OVERRIDE=rk; CallGenes.R58_CONSENSUS_OVERRIDE=rc; CallGenes.R58_MODELS_OVERRIDE=rm;
			CallGenes.LSU_KMERS_OVERRIDE=lk; CallGenes.LSU_CONSENSUS_OVERRIDE=lc; CallGenes.LSU_MODELS_OVERRIDE=lm;
		}
	}
}
