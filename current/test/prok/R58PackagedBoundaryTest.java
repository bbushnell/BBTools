package test.prok;

import java.util.ArrayList;
import prok.CallGenes;
import prok.FamilyResourceResolver;
import prok.GeneCaller;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.R58BoundaryConfig;

/** Packaged 5.8S inference and explicit opt-out must follow real family loading.
 * @author Raiden */
public final class R58PackagedBoundaryTest {
	public static void main(String[] args){
		check(args.length<=1 && (args.length==0 || args[0].startsWith("barealias=")),"Optional barealias names a producer-created installed test resource");
		final boolean pair=CallGenes.R58LSU_ENABLED,nn=CallGenes.R58_BOUNDARY_NN_ENABLED,explicit=CallGenes.R58_BOUNDARY_NN_EXPLICIT;
		final boolean umbrella=CallGenes.NCRNA_FAMILIES_ENABLED;
		final String oldNet=CallGenes.R58_BOUNDARY_NET_OVERRIDE;
		final ArrayList<NcrnaFamily> old=new ArrayList<NcrnaFamily>(GeneCaller.ncrnaFamilies);
		try{
			CallGenes.R58_BOUNDARY_CONFIG.clear();CallGenes.R58_BOUNDARY_NN_EXPLICIT=false;
			CallGenes.setR58LsuEnabled(false);check(!CallGenes.R58_BOUNDARY_NN_ENABLED,"A disabled eukaryotic family must not request boundary resources");
			CallGenes.setR58LsuEnabled(true);CallGenes.validateNcrnaGateCombo();
			check(CallGenes.R58_BOUNDARY_NN_ENABLED,"The selected 5.8S family enables packaged refinement by default");
			GeneCaller.ncrnaFamilies.clear();CallGenes.loadNcrnaResources();
			final NcrnaFamily on=family("r58"),lsu=family("lsu");
			check(on.boundary5NetTemplate!=null && on.boundary3NetTemplate!=null && on.boundary5NetTemplate.numInputs()==10
				&& on.boundary5NetTemplate.numOutputs()==1,"Manifest must load the trained shared ten-input scorer");
			check(on.boundaryStartInside==7 && on.boundaryStopInside==7 && on.boundaryStartOutside==0 && on.boundaryStopOutside==0
				&& on.boundaryMeanLen==153f,"Inference must retain all-inside7-mer tables and the trained length normalization");
			check(on.boundaryStartOffsets.length==68 && on.boundaryStartOffsets[0]==-11 && on.boundaryStartOffsets[67]==56
				&& on.boundaryStopOffsets.length==84 && on.boundaryStopOffsets[0]==-46 && on.boundaryStopOffsets[83]==37
				&& on.boundaryOnRawEndpoints,"Selected ranges and retained-anchor eligibility must reach the loaded family");
			check(lsu.boundary5NetTemplate==null && lsu.boundary3NetTemplate==null,"5.8S packaging must not enable LSU refinement");
			CallGenes.setR58BoundaryEnabled(false);CallGenes.setR58LsuEnabled(false);CallGenes.setR58LsuEnabled(true);
			CallGenes.validateNcrnaGateCombo();check(!CallGenes.R58_BOUNDARY_NN_ENABLED,"Explicit boundary opt-out wins before family selection");
			GeneCaller.ncrnaFamilies.clear();CallGenes.loadNcrnaResources();
			final NcrnaFamily off=family("r58");check(off.boundary5NetTemplate==null && off.boundaryStartTable==null && !off.boundaryOnRawEndpoints,"Disabled refinement must load no net/table and no raw-anchor path");
			CallGenes.R58_BOUNDARY_NN_EXPLICIT=false;CallGenes.setR58LsuEnabled(true);CallGenes.setR58BoundaryEnabled(false);
			check(!CallGenes.R58_BOUNDARY_NN_ENABLED,"Explicit boundary opt-out wins after family selection");
			final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();options.boundary=true;options.resources.boundary="/nonexistent/r58_fixture_missing.bbnet";
			boolean rejected=false;try{NcrnaFamilyLoader.load("r58",options,new FamilyResourceResolver());}catch(IllegalArgumentException e){rejected=e.getMessage().contains("Missing required r58 resource boundary");}
			check(rejected,"An explicit missing network must fail rather than fall back to the packaged network");
			final String packaged=new FamilyResourceResolver().resolve("r58","boundary",null).path;
			check(CallGenes.explicitNcrnaResource("r58","boundary",null)==null,"Omitted overrides must delegate to the manifest");
			check(CallGenes.explicitNcrnaResource("r58","boundary","?callgenes_5.8S/boundary.bbnet").equals(packaged),"Resource-prefix overrides must resolve to the same trained payload");
			check(CallGenes.explicitNcrnaResource("r58","fixture","model.pgm.gz").equals(dna.Data.findPath("?model.pgm.gz",false)),"A bare installed resource name must retain lookup semantics");
			rejected=false;try{CallGenes.explicitNcrnaResource("r58","boundary","r58_fixture_absent_boundary.bbnet");}catch(IllegalArgumentException e){rejected=true;}
			check(rejected,"A missing bare override must not become null and silently select the packaged network");
			if(args.length==1){
				final String alias=args[0].substring(10);check(!alias.isEmpty() && alias.indexOf('/')<0,"Bare-name test must not collapse to an explicit path");
				CallGenes.setR58BoundaryEnabled(true);CallGenes.R58_BOUNDARY_NET_OVERRIDE=alias;
				GeneCaller.ncrnaFamilies.clear();CallGenes.loadNcrnaResources();
				check(family("r58").boundary5NetTemplate!=null,"The actual family factory must load a bare installed network name");
				CallGenes.R58_BOUNDARY_NET_OVERRIDE="r58_fixture_absent_boundary.bbnet";GeneCaller.ncrnaFamilies.clear();
				rejected=false;try{CallGenes.loadNcrnaResources();}catch(IllegalArgumentException e){rejected=e.getMessage().contains("Missing required r58 resource boundary");}
				check(rejected,"The actual family factory must reject a missing bare network name instead of packaged fallback");
			}
			System.out.println("PASS packaged 5.8S network/tables/ranges, raw-anchor default, argument-order opt-out, missing-resource failure and LSU isolation");
		}finally{
			GeneCaller.ncrnaFamilies.clear();GeneCaller.ncrnaFamilies.addAll(old);
			CallGenes.R58LSU_ENABLED=pair;CallGenes.R58_BOUNDARY_NN_ENABLED=nn;CallGenes.R58_BOUNDARY_NN_EXPLICIT=explicit;
			CallGenes.NCRNA_FAMILIES_ENABLED=umbrella;CallGenes.R58_BOUNDARY_CONFIG.clear();
			CallGenes.R58_BOUNDARY_NET_OVERRIDE=oldNet;
		}
	}
	private static NcrnaFamily family(String name){for(NcrnaFamily f:GeneCaller.ncrnaFamilies){if(f.name.equals(name)){return f;}}throw new AssertionError("Missing family "+name);}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
