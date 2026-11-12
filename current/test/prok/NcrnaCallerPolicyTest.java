package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.FamilyResourceResolver;
import prok.GeneCaller;
import prok.NcrnaCallerPolicy;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaFamilyRegistry;
import prok.NcrnaJointPolicy;
import prok.NcrnaScavenger;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Exact family isolation and legacy/patience propagation through native owners.
 * The final loader check uses the installed SAM-I package, not genomic inputs.
 * No calling, grading, alignment kernel or scientific policy selection is run.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCallerPolicyTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaCallerPolicyTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown test option: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && parser.in1==null,
				"Supply out=fresh_fixture_directory with no input or overwrite");
			new NcrnaCallerPolicyTest(new File(parser.out1)).run();
			pp.outstream.println("NCRNA_CALLER_POLICY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaCallerPolicyTest(final File root_){
		root=root_;
		check(root.mkdir(), "Policy fixture requires a fresh directory under an existing parent");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	private void run(){
		final String rows="candidate\tcobalamin_rf00174\tadaptive\t384\t2\t1\t0\t0\t8\n"
			+"candidate\tcobalamin_rf01482\tadaptive\t39\t1\t1\t0.25\t0.00390625\t2147483647\n"
			+"candidate\tcobalamin_rf01689\tadaptive\t8\t1\t1\t0\t0\t2147483647\n"
			+"candidate\tsam1\tadaptive\t2\t2\t1\t0.25\t0.0625\t0\n";
		final NcrnaCallerPolicy policy=load("valid.tsv", HEADER+rows);
		final ArrayList<NcrnaFamily> families=new ArrayList<NcrnaFamily>();
		final String[] names={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689", "sam1"};
		final int[] caps={384, 39, 8, 2}, patience={8, Integer.MAX_VALUE, Integer.MAX_VALUE, 0};
		for(int i=0; i<names.length; i++){
			final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
			policy.apply(names[i], options);
			check(options.topN==caps[i] && options.modelPatience==patience[i] && options.minHits==1,
				"Each canonical member must consume its own literal tuple: "+names[i]);
			check(options.indexK==-1 && options.seedK==-1 && Float.isNaN(options.idPass)
				&& Float.isNaN(options.scoreB), "Policy application must not supply mapping K, seeds, identity or path weights");
			families.add(family(names[i], options));
		}
		policy.validateLoaded(families);
		check(families.get(0).adaptFloor==2 && families.get(1).adaptTopFrac==.25f
			&& families.get(1).adaptQFrac==.00390625f && families.get(2).adaptTopFrac==0,
			"Minority-member fractions must not inherit the main member's controls");
		final NcrnaFamilyLoader.Options unrelated=new NcrnaFamilyLoader.Options();
		unrelated.topN=17;
		policy.apply("tpp", unrelated);
		check(unrelated.topN==17 && unrelated.modelPatience==NcrnaFamilyLoader.Options.UNSET_PATIENCE,
			"An omitted family must retain caller-owned settings and registry patience");
		for(final NcrnaFamily family:families){
			final NcrnaScavenger worker=worker(family);
			check(worker.modelPatience()==family.modelPatience,
				"The actual GeneCaller configuration path must install verifier stopping: "+family.name);
		}
		final NcrnaFamilyLoader.Options legacyOptions=new NcrnaFamilyLoader.Options();
		legacyOptions.topN=3;
		legacyOptions.adaptive=false;
		legacyOptions.floor=0;
		legacyOptions.minHits=1;
		legacyOptions.topFrac=0;
		legacyOptions.qFrac=0;
		legacyOptions.modelPatience=-1;
		check(worker(family("sam1", legacyOptions)).modelPatience()==-1, "An explicit legacy policy must retain native stopping");
		final ArrayList<NcrnaFamily> missing=new ArrayList<NcrnaFamily>(families);
		missing.remove(2);
		rejected(()->policy.validateLoaded(missing), "not enabled/loaded");
		final ArrayList<NcrnaFamily> repeated=new ArrayList<NcrnaFamily>(families);
		repeated.add(families.get(0));
		rejected(()->policy.validateLoaded(repeated), "Repeated loaded");
		families.get(0).modelPatience=7;
		rejected(()->policy.validateLoaded(families), "complete explicit policy");
		families.get(0).modelPatience=8;
		rejected(()->load("two.tsv", HEADER+rows+rows.replace("candidate", "second")), "exactly one policy");
		rejected(()->load("group.tsv", HEADER+rows.replace("cobalamin_rf00174", "cobalamin")), "canonical caller IDs");
		rejected(()->load("alias.tsv", HEADER+"candidate\tSAMI\tfixed\t2\t0\t1\t0\t0\t-1\n"), "No unified ncRNA resource registration");
		rejected(()->load("trna.tsv", HEADER+"candidate\ttRNA\tfixed\t2\t0\t1\t0\t0\t-1\n"), "generic shared-boundary caller");
		rejected(()->load("negative.tsv", HEADER+rows.replace("\t8\n", "\t-2\n")), "patience must be");
		rejected(()->NcrnaCallerPolicy.load(null), "requires a joint-policy");
		loader(policy);
		packagedPolicies();
	}

	/** Literal owner-selected tuples are independent of the integration fixture's package-derived root. */
	private static void packagedPolicies(){
		final NcrnaJointPolicy.Member[] expected={
			new NcrnaJointPolicy.Member("sixs_rf00013", true, 56, 1f, 1, 0f, 0f, 0),
			new NcrnaJointPolicy.Member("sixs_rf01685", true, 5, 16f, 1, 0f, 0f, Integer.MAX_VALUE),
			new NcrnaJointPolicy.Member("cobalamin_rf00174", true, 384, 2f, 1, 0f, 0f, 8),
			new NcrnaJointPolicy.Member("cobalamin_rf01482", false, 100, 0f, 1, 0f, 0f, -1),
			new NcrnaJointPolicy.Member("cobalamin_rf01689", false, 100, 0f, 1, 0f, 0f, -1),
			new NcrnaJointPolicy.Member("sam1", false, 238, 0f, 1, 0f, 0f, 16),
			new NcrnaJointPolicy.Member("tpp", true, 361, 5f, 1, .1875f, 0f, -1)
		};
		final int[] mappingK={7, 7, 9, 9, 9, 7, 7};
		for(int i=0; i<expected.length; i++){
			final NcrnaJointPolicy.Member policy=expected[i];
			final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
			options.endpoints=false;
			options.boundaryOverride=Boolean.FALSE;
			final NcrnaFamily family=NcrnaFamilyLoader.load(policy.name, options, new FamilyResourceResolver());
			final String view=NcrnaFamilyRegistry.get(policy.name).group==null ? "target" : "union";
			policy.validateBaseline(family, view);
			check(family.indexK==mappingK[i] && worker(family).modelPatience()==policy.patience,
				"Literal package controls must reach the actual verifier without an explicit caller policy: "+policy.name);
		}
	}

	/** Exercises the production loader and its sentinel without opening any genomes. */
	private void loader(final NcrnaCallerPolicy policy){
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.endpoints=false;
		options.boundaryOverride=Boolean.FALSE;
		policy.apply("sam1", options);
		final NcrnaFamily loaded=NcrnaFamilyLoader.load("sam1", options, new FamilyResourceResolver());
		check(loaded.modelPatience==0 && loaded.indexTopN==2 && loaded.adaptive
			&& loaded.adaptFloor==2 && loaded.adaptTopFrac==.25f && loaded.adaptQFrac==.0625f,
			"Loaded immutable model data must carry the exact explicit mapping/verifier tuple");
		check(worker(loaded).modelPatience()==0, "The real loaded family must reach the actual worker configuration");
		options.modelPatience=NcrnaFamilyLoader.Options.UNSET_PATIENCE;
		check(NcrnaFamilyLoader.load("sam1", options, new FamilyResourceResolver()).modelPatience==16,
			"Omitted patience must inherit the packaged SAM-I budget despite explicit mapping overrides");
		final NcrnaFamilyLoader.Options packagedOptions=new NcrnaFamilyLoader.Options();
		packagedOptions.endpoints=false;
		packagedOptions.boundaryOverride=Boolean.FALSE;
		final NcrnaFamily packaged=NcrnaFamilyLoader.load("sam1", packagedOptions, new FamilyResourceResolver());
		final NcrnaJointPolicy.Member baseline=new NcrnaJointPolicy.Member("sam1", false, 238, 0f, 1, 0f, 0f, 16);
		baseline.validateBaseline(packaged, "target");
		check(worker(packaged).modelPatience()==16,
			"A packaged nonlegacy baseline must reach the production worker and pass sweep validation");
		final NcrnaJointPolicy.Member stale=new NcrnaJointPolicy.Member("sam1", false, 238, 0f, 1, 0f, 0f, -1);
		rejected(()->stale.validateBaseline(packaged, "target"), "preserve every packaged control");
		packagedOptions.modelPatience=0;
		final NcrnaFamily overridden=NcrnaFamilyLoader.load("sam1", packagedOptions, new FamilyResourceResolver());
		final NcrnaJointPolicy.Member disguised=new NcrnaJointPolicy.Member("sam1", false, 238, 0f, 1, 0f, 0f, 0);
		rejected(()->disguised.validateBaseline(overridden, "target"), "packaged mapping mode");
		options.modelPatience=-2;
		rejected(()->NcrnaFamilyLoader.load("sam1", options, new FamilyResourceResolver()), "patience must be");
		options.modelPatience=0;
		rejected(()->NcrnaFamilyLoader.load("tRNA", options, new FamilyResourceResolver()), "generic shared-boundary caller");
	}

	private static NcrnaFamily family(final String name, final NcrnaFamilyLoader.Options options){
		final NcrnaFamily family=new NcrnaFamily(name,
			new byte[][]{bytes("ACGT"), bytes("AAAA"), bytes("CCCC")}, null, null, null,
			17, 17, 30, 1, options.topN, options.adaptive, options.floor, options.topFrac, options.qFrac, options.minHits);
		family.modelPatience=options.modelPatience;
		return family;
	}

	/** Calls the same public production configuration seam as lazy worker initialization. */
	private static NcrnaScavenger worker(final NcrnaFamily family){
		final NcrnaScavenger worker=new NcrnaScavenger(family.library, family.models, family.modelNames,
			family.kmerSet, family.kLong, family.minLen, family.windowPad,
			family.indexK, family.indexTopN, family.adaptive, family.adaptFloor, family.adaptTopFrac,
			family.adaptQFrac, family.fixedMinHits, family.scoreA, family.scoreB, family.idPass, family.idBorderline,
			family.boundary5NetTemplate, family.boundary3NetTemplate, family.boundaryStartTable, family.boundaryStopTable,
			family.boundaryStartInside, family.boundaryStartOutside, family.boundaryStopInside, family.boundaryStopOutside,
			family.boundaryMeanLen, family.boundaryStartOffsets, family.boundaryStopOffsets,
			family.boundaryMarginStart, family.boundaryMarginStop);
		GeneCaller.applyNcrnaFamilyControls(worker, family);
		return worker;
	}

	private NcrnaCallerPolicy load(final String file, final String text){
		final String path=new File(root, file).toString();
		try(Output out=new Output(path)){out.writer.print(bytes(text));}
		return NcrnaCallerPolicy.load(path);
	}
	private static void rejected(final Runnable action, final String diagnostic){
		try{action.run();}
		catch(IllegalArgumentException | IllegalStateException expected){
			check(expected.getMessage().contains(diagnostic), "Unexpected policy failure: "+expected);
			return;
		}
		throw new IllegalStateException("Expected rejection containing: "+diagnostic);
	}
	private static byte[] bytes(final String text){return text.getBytes(StandardCharsets.UTF_8);}
	private static void check(final boolean valid, final String why){
		if(!valid){throw new IllegalStateException(why);}
	}
	private final File root;
	private static final String HEADER="policy\tmember\tmode\ttopn\tfloor\tminhits\ttopfrac\tqfrac\tpatience\n";
}
