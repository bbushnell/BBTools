package test.prok;

import java.util.Arrays;
import prok.CallGenes;
import prok.GeneCaller;
import prok.NcrnaFamily;
import prok.Orf;
import prok.ProkObject;
import prok.R58BoundaryConfig;

/** Tests scoped SAM-I boundary controls against core-only and boundary fixture roots.
 * Run in a fresh JVM; the supplied boundary network is a loading fixture, not a
 * trained SAM-I model or scientific accuracy test.
 * @author Brian Bushnell, Raiden
 */
public final class CallGenesSam1BoundaryTest {

	public static void main(String[] args){
		check(args.length==2, "Supply core-only and boundary-fixture family roots");
		check(!CallGenes.SAM1_ENABLED && CallGenes.SAM1_BOUNDARY_NN_OVERRIDE==null
			&& !CallGenes.NCRNA_BOUNDARY_NN_ENABLED && GeneCaller.ncrnaFamilies.isEmpty(),
			"Configuration checks require fresh startup state");
		testAliases();
		testOutputLabel();
		testDefaults(args[1]);
		testGates();
		testLoading(args[0], args[1]);
		System.out.println("SAM1_BOUNDARY_CONFIG_PASS gates=PASS ranges=PASS loader=PASS off_no_boundary_io=PASS");
	}
	private static void testAliases(){
		check(!CallGenes.SAM1_ENABLED, "Neither activation flag is present at fresh startup");
		for(String key:new String[]{"sam1", "sami", "SaM1", "SaMi"}){
			check(CallGenes.parseSam1Flag(key, "t") && CallGenes.SAM1_ENABLED,
				"Both case-insensitive spellings must activate the same SAM-I field: "+key);
			check(CallGenes.parseSam1Flag(key, "f") && !CallGenes.SAM1_ENABLED,
				"Both case-insensitive spellings must also disable the family: "+key);
		}
	}
	private static void testOutputLabel(){
		final byte[] bases=new byte[128];Arrays.fill(bases, (byte)'A');
		for(int strand=0; strand<2; strand++){
			final Orf orf=new Orf("fixture", 10, 100, strand, 0, bases, false, ProkObject.RNA);
			orf.ncrnaFamily="sam1";
			check(orf.toGff().endsWith(",family:SAM-I") && "sam1".equals(orf.ncrnaFamily),
				"GFF must display SAM-I without changing the internal key used by caller accounting");
			orf.ncrnaFamily="fmn";
			check(orf.toGff().endsWith(",family:FMN"), "Another family retains its own selected scientific display name");
		}
	}
	private static void testDefaults(String boundaryRoot){
		final boolean oldFamilies=CallGenes.NCRNA_FAMILIES_ENABLED;
		CallGenes.parseSam1Flag("SaMi", "t");CallGenes.NCRNA_FAMILIES_ENABLED=true;
		CallGenes.SAM1_RESOURCE_ROOT=boundaryRoot;
		CallGenes.validateSam1BoundaryConfig();CallGenes.loadSam1Resources();
		final NcrnaFamily family=onlyFamily();
		check(family.boundary5NetTemplate!=null && family.boundary3NetTemplate!=null,
			"Enabling SAM-I alone must load the accepted endpoint networks without a global boundary flag");
		check(family.indexTopN==238, "Bare SAM-I must use the selected238-model shortlist cap");
		checkRanges(family, 12);GeneCaller.ncrnaFamilies.clear();
		family.boundaryStartOffsets[0]=0;
		CallGenes.loadSam1Resources();checkRanges(onlyFamily(), 12);GeneCaller.ncrnaFamilies.clear();
		CallGenes.SAM1_ENABLED=false;CallGenes.NCRNA_FAMILIES_ENABLED=oldFamilies;
		CallGenes.SAM1_RESOURCE_ROOT=null;
		CallGenes.validateSam1BoundaryConfig();
	}

	private static void testGates(){
		check(CallGenes.parseSam1Flag("sam1boundarynn", "f"), "SAM-I parser must consume the explicit disable flag");
		CallGenes.parseSam1Flag("sam1root", "/absent/sam1_boundary_fixture");
		CallGenes.loadSam1Resources();
		check(GeneCaller.ncrnaFamilies.isEmpty(), "Disabled SAM-I must not open its missing resource root");
		CallGenes.parseSam1Flag("sam1boundarynn", "t");
		expectFailure(()->CallGenes.validateSam1BoundaryConfig(), "sam1=t and ncrna=t");
		CallGenes.parseSam1Flag("sam1", "t");
		CallGenes.NCRNA_FAMILIES_ENABLED=true;
		CallGenes.validateSam1BoundaryConfig();
		CallGenes.parseSam1Flag("SAM1BOUNDARYSTARTRANGE", "-8,8");
		CallGenes.validateSam1BoundaryConfig();
		for(String invalid:new String[]{"1,8", "-8,-1", "-257,0", "-8,8,9", "bad"}){
			expectFailure(()->CallGenes.parseSam1Flag("sam1boundarystoprange", invalid), "Boundary range");
		}
		CallGenes.parseSam1Flag("sam1boundarystoprange", "-12,12");
		CallGenes.validateSam1BoundaryConfig();
		CallGenes.NCRNA_FAMILIES_ENABLED=false;
		expectFailure(()->CallGenes.validateSam1BoundaryConfig(), "sam1=t and ncrna=t");
		CallGenes.NCRNA_FAMILIES_ENABLED=true;
		CallGenes.parseSam1Flag("sam1boundarynn", "f");
		CallGenes.NCRNA_BOUNDARY_NN_ENABLED=true;
		check(!CallGenes.ncrnaBoundaryEnabled("sam1") && CallGenes.ncrnaBoundaryEnabled("fmn"),
			"SAM-I explicit false must override the global switch without changing another family");
		CallGenes.SAM1_BOUNDARY_NN_OVERRIDE=null;
		check(CallGenes.ncrnaBoundaryEnabled("sam1"), "An absent SAM-I override retains selected default-on refinement");
		CallGenes.NCRNA_BOUNDARY_NN_ENABLED=false;
		check(CallGenes.ncrnaBoundaryEnabled("sam1"), "SAM-I refinement defaults on independently of the global switch");
		CallGenes.parseSam1Flag("sam1boundarynn", "t");
		check(CallGenes.ncrnaBoundaryEnabled("sam1") && !CallGenes.ncrnaBoundaryEnabled("fmn"),
			"SAM-I explicit true must not enable boundary loading in unrelated families");
	}

	private static void testLoading(String coreRoot, String boundaryRoot){
		CallGenes.parseSam1Flag("sam1root", coreRoot);
		CallGenes.parseSam1Flag("sam1boundarynn", "f");
		CallGenes.loadSam1Resources();
		final NcrnaFamily off=onlyFamily();
		check(off.boundary5NetTemplate==null && off.boundary3NetTemplate==null,
			"The core-only root has no boundary resources; explicit off must load successfully");
		checkRanges(off, 12);
		GeneCaller.ncrnaFamilies.clear();
		CallGenes.parseSam1Flag("sam1boundarynn", "t");
		expectFailure(()->CallGenes.loadSam1Resources(), "Missing required family manifest role");
		check(GeneCaller.ncrnaFamilies.isEmpty(), "Missing boundary resources must not publish a partial family");
		CallGenes.parseSam1Flag("sam1root", boundaryRoot);
		CallGenes.loadSam1Resources();
		final NcrnaFamily on=onlyFamily();
		check(on.boundary5NetTemplate!=null && on.boundary3NetTemplate!=null,
			"Enabled refinement must load both boundary network templates");
		checkRanges(on, 12);
		check(!on.boundaryOnRawEndpoints && on.trimAlignmentExtent,
			"New scoped flags must preserve existing trim and raw-endpoint invocation policy");
		GeneCaller.ncrnaFamilies.clear();
		on.boundaryStartOffsets[0]=0;
		CallGenes.parseSam1Flag("sam1boundarystoprange", "-20,20");
		CallGenes.loadSam1Resources();
		checkRanges(onlyFamily(), 20);
		check(Arrays.deepEquals(off.library, on.library) && Arrays.equals(off.modelNames, on.modelNames)
			&& off.kmerSet.size()==on.kmerSet.size(), "Boundary loading must preserve the fixture's core model/seed library");
	}

	private static NcrnaFamily onlyFamily(){
		check(GeneCaller.ncrnaFamilies.size()==1, "Successful SAM-I load must publish exactly one family");
		final NcrnaFamily family=GeneCaller.ncrnaFamilies.get(0);
		check("sam1".equals(family.name), "Boundary controls must attach only to the selected family");
		return family;
	}

	private static void checkRanges(NcrnaFamily family, int stopRadius){
		check(Arrays.equals(family.boundaryStartOffsets, R58BoundaryConfig.range("-8,8"))
			&& Arrays.equals(family.boundaryStopOffsets, R58BoundaryConfig.range(-stopRadius+","+stopRadius)),
			"Explicit inclusive ranges must reach the loaded family without exposing mutable config arrays");
	}

	private static void expectFailure(Runnable action, String message){
		assert(action!=null && message!=null) : "Negative fixtures must exercise the real parser or loader and identify the expected failure";
		try{action.run();}catch(IllegalArgumentException expected){
			check(expected.getMessage()!=null && expected.getMessage().contains(message), "Unexpected failure: "+expected);
			return;
		}
		throw new AssertionError("Expected failure containing: "+message);
	}

	private static void check(boolean condition, String why){if(!condition){throw new AssertionError(why);}}
}
