package test.prok;

import java.util.Arrays;

import prok.CallGenes;
import prok.ConservedRnaSeedIndex;
import prok.FamilyResourceResolver;
import prok.GeneCaller;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaFamilyRegistry;

/** Checks TPP loader/caller agreement against a supplied tiny resource fixture.
 * Each invocation needs a fresh JVM; no scientific accuracy claim is made.
 * @author Brian, Raiden
 */
public final class CallGenesTppTest {

	public static void main(String[] args){
		check(args.length==1, "Supply one TPP loading-fixture directory");
		check(!CallGenes.TPP_ENABLED && GeneCaller.ncrnaFamilies.isEmpty(), "Fresh startup must leave TPP disabled");
		check(!CallGenes.TPP_CDS_STOP, "The CDS stop ceiling must remain opt-in");
		CallGenes.parseTppFlag("tpproot", args[0]+"/absent");
		CallGenes.loadTppResources();
		check(GeneCaller.ncrnaFamilies.isEmpty(), "Disabled TPP must not open even an explicitly missing root");
		check("tpp".equals(CallGenes.parseNcrnaFamily("TPP")), "Family selectors must normalize to the registered key");
		check("TPP".equals(NcrnaFamilyRegistry.displayName("tpp")), "GFF attribution must use the biological TPP label");
		check("tpp".equals(NcrnaFamilyRegistry.modelOverrideMember("tpp").id), "The existing reducer must accept TPP model overrides");

		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.endpoints=false;
		options.resources.root=args[0];
		final NcrnaFamily direct=NcrnaFamilyLoader.load("tpp", options, new FamilyResourceResolver());
		check(direct.library.length==3 && direct.indexTopN==361 && direct.windowPad==232,
			"A replacement library retains the selected production cap and padding");
		check(direct.adaptive && direct.adaptFloor==5f && direct.adaptTopFrac==.1875f && direct.adaptQFrac==0f,
			"The selected adaptive mapping policy must reach direct loading");
		check(direct.kLong==ConservedRnaSeedIndex.K && direct.seedMinHits==2 && direct.seedDistinct,
			"Direct loading must retain the shared seed K and distinct two-key gate");
		CallGenes.NCRNA_FAMILIES_ENABLED=true;
		CallGenes.parseTppFlag("TPP", "t");
		CallGenes.parseTppFlag("tpproot", args[0]);
		CallGenes.loadTppResources();
		final NcrnaFamily called=onlyFamily();
		check(called.cdsStopConstraintFactory==null, "Omitting tppcdsstop must leave the endpoint constraint absent");
		check(called.windowPad==direct.windowPad && called.indexTopN==direct.indexTopN
			&& called.adaptive==direct.adaptive && called.adaptFloor==direct.adaptFloor
			&& called.adaptTopFrac==direct.adaptTopFrac && called.adaptQFrac==direct.adaptQFrac
			&& called.seedMinHits==direct.seedMinHits && called.seedDistinct==direct.seedDistinct
			&& Arrays.deepEquals(called.library, direct.library) && Arrays.equals(called.modelNames, direct.modelNames),
			"Direct loader and CallGenes must bind the same TPP models and controls");
		GeneCaller.ncrnaFamilies.clear();

		options.windowPad=123;
		options.topN=1;
		options.seedMinHits=1;
		final NcrnaFamily overridden=NcrnaFamilyLoader.load("tpp", options, new FamilyResourceResolver());
		check(overridden.windowPad==123 && overridden.indexTopN==1
			&& overridden.seedMinHits==1 && !overridden.seedDistinct,
			"Explicit controls must override derived defaults in direct family loading");
		System.out.println("TPP_LOADER_CALLER_AGREEMENT_PASS");
	}

	private static NcrnaFamily onlyFamily(){
		check(GeneCaller.ncrnaFamilies.size()==1, "A selected TPP factory must append exactly one family");
		final NcrnaFamily family=GeneCaller.ncrnaFamilies.get(0);
		check("tpp".equals(family.name), "TPP must retain a distinct registered caller identity");
		return family;
	}

	private static void check(boolean condition, String why){
		if(!condition){throw new AssertionError(why);}
	}
}
