package prok;

import parse.LineParser1;
import prok.NcrnaTextIO.Input;

/** Immutable serving policy and shared resources for the three-member caller.
 * Loading is cold; each GeneCaller creates its own ranking, CM and acceptance state.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminCallerResources{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public CobalaminCallerResources(final Policy policy_, final NcrnaFamily[] families_){
		require(policy_!=null && families_!=null && families_.length==3,
			"The packaged caller binds one declared policy and all three loaded members");
		policy=policy_;
		families=families_.clone();
		require(families[1]!=null && families[1].boundaryRefinerFactory instanceof NcrnaPositionalBoundaryRefiner.Resources,
			"The packaged second member requires its selected positional endpoint resources");
		final NcrnaPositionalBoundaryRefiner.Resources tips=(NcrnaPositionalBoundaryRefiner.Resources)families[1].boundaryRefinerFactory;
		require(tips.fullErrors() && tips.radius()==10 && tips.wordLength()==9 && tips.modelOffsets()!=null && !tips.gated(),
			"The selected second-member route uses centered K9 full-error tips without an experimental movement gate");
		require(families[2]!=null && families[2].boundaryRefinerFactory==null,
			"The selected third member retains its native endpoint route");
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final FamilyResourceResolver.Overrides overrides=new FamilyResourceResolver.Overrides();
		overrides.root=policy.root;
		final String[] paths=new String[3];
		final String[][] names=new String[3][];
		final byte[][][] sequences=new byte[3][][];
		for(int member=0; member<3; member++){
			final NcrnaFamily family=families[member];
			require(family!=null && FAMILIES[member].equals(family.name) && family.idPass==FLOORS[member]
				&& family.idBorderline==FLOORS[member] && family.seedMinHits==2 && family.seedDistinct,
				"The packaged policy requires its ordered member axes, native identity floors and two distinct discovery keys");
			names[member]=family.modelNames;
			sequences[member]=family.library;
			paths[member]=resolver.resolve(family.name, "ddl_profiles", overrides).path;
			caps[member]=family.keyBounds!=null ? family.keyBounds : NcrnaOverlapWindowBuilder.Caps.load(
				resolver.resolve(family.name, "key_bounds", overrides).path, MEMBERS[member], family.kmerSet.size(), policy.percentile);
			requireLegacyBounds(caps[member]);
			caps[member].validateKeys(family.kmerSet.toArray());
			require(policy.pads[member]<=family.windowPad,
				"Member proposal padding cannot exceed the loaded caller's outer discovery padding: "+family.name);
		}
		profiles=NcrnaDdlProfileResource.loadClanWithConsensusReplacement(paths, names, sequences);
		require(profiles.memberOnly.orderOnly() && profiles.memberOnly.buckets()==policy.bins
			&& profiles.memberOnly.seed()==12345 && !profiles.memberOnly.consensus()
			&& !profiles.memberOnly.suppressed() && profiles.memberOnly.orderingBonus()==0
			&& profiles.replacement.consensusReplaced(),
			"Serving uses the member-only index and its nonempty-consensus replacement, without bonus votes or suppression");
		final NcrnaPathScoreConfig scores=new NcrnaPathScoreConfig(resolver.resolve(FAMILIES[0], "path_policy", overrides).path);
		scores.bind(families);
		scoreCell=scores.cell(policy.pathCell);
		scoreCell.validateRuntime(families);
		acceptance=CobalaminCmAcceptance.Resources.serving(families, policy.root,
			resolver.resolve(FAMILIES[0], "acceptance_policy", overrides).path, profiles.memberOnly, policy.depth);
		for(final float upper:acceptance.upper){require(upper==.75f, "The packaged alignment-gated policy must bypass at identity0.75");}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** The merged cobalamin route still uses its reviewed external ceilings and skips unmeasured keys. */
	public static void requireLegacyBounds(final NcrnaOverlapWindowBuilder.Caps caps){
		require(caps!=null && caps.chainCeiling==0 && caps.constructionPopulation==null,
			"Merged cobalamin requires its existing legacy bounds and external policy; new generic ceiling/population metadata is not active on this route");
	}

	/** Called once when the fresh worker's clan is built, before its first query. */
	void configure(final NcrnaClanScavenger clan){
		require(clan!=null && clan.nativeRescue()==null, "Packaged serving cannot inherit an experimental rescue walk");
		int model=0;
		for(int member=0; member<3; member++){
			for(final String name:families[member].modelNames){
				require(model<clan.modelCount() && name.equals(clan.modelName(model)) && clan.member(model)==member,
					"The worker must retain the exact loaded model order and member ownership");
				model++;
			}
		}
		require(model==clan.modelCount(), "The packaged worker cannot add or omit a bound model");
		clan.setModelPatience(policy.patience);
		clan.setWinnerRule(NcrnaMatchedBaseWinner.MODEL_SCORE, .58);
		clan.setScorePatience(policy.shortlist);
		clan.setConsensusPolicy(profiles, false);
		clan.setDdlOnlyAdmission(policy.minimumHits);
		clan.setPathScoreCell(scoreCell);
		clan.setMemberProposalPads(policy.pads);
		clan.setOverlapWindows(true, caps);
		clan.setWindowCeilings(policy.ceilings);
		clan.setCandidateAcceptance(new CobalaminCmAcceptance(acceptance, clan, 0));
	}
	NcrnaDdlProfileResource memberProfiles(){return profiles.memberOnly;}

	/** A caller's explicitly reset/replaced family roster cannot inherit old assets. */
	boolean matches(final java.util.List<NcrnaFamily> roster){
		require(roster!=null, "Serving resources must bind an actual loaded family roster");
		int mask=0;
		for(final NcrnaFamily family:roster){
			for(int member=0; member<families.length; member++){
				if(FAMILIES[member].equals(family.name)){
					final int bit=1<<member;
					if(family!=families[member] || (mask&bit)!=0){return false;}
					mask|=bit;
				}
			}
		}
		return mask==(1<<families.length)-1;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Packaged Policy        ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Policy{
		private Policy(final String path, final String root_){
			root=root_;
			final LineParser1 p=new LineParser1('\t');
			try(Input input=new Input(path)){
				RrnaResourceIO.header(input.reader, "#BBTools cobalamin caller v1", path);
				RrnaResourceIO.header(input.reader, HEADER, path);
				final byte[] row=input.reader.nextLine();
				require(row!=null, "A declared caller policy needs one complete operating point: "+path);
				p.set(row);
				require(p.terms()==16, "The caller policy must declare all ranking, window, path and acceptance fields: "+path);
				bins=p.parseInt(0);
				shortlist=p.parseInt(1);
				patience=p.parseInt(2);
				minimumHits=p.parseInt(3);
				percentile=p.parseInt(4);
				for(int member=0; member<3; member++){
					pads[member]=p.parseInt(5+member);
					ceilings[member]=p.parseInt(8+member);
					require(pads[member]>=0 && ceilings[member]>0, "Window padding and chain ceilings are base counts");
				}
				pathCell=p.parseString(11);
				depth=p.parseInt(13);
				require(bins==1024 && shortlist==100 && patience==4 && minimumHits>0 && minimumHits<=bins
					&& !pathCell.isEmpty() && p.parseFloat(12)==.60f && depth>=1 && depth<=shortlist
					&& p.termEquals("learned", 14) && p.termEquals("true", 15),
					"The packaged merged caller requires supported ranking, aligned acceptance and learned main ends with GA hold");
				NcrnaOverlapWindowBuilder.capRule(percentile);
				require(input.reader.nextLine()==null, "A caller package selects one operating point, never an implicit sweep");
			}
		}

		/** An absent policy retains complete legacy packages used by development tools.
		 * A partial new package must fail rather than mix new and old serving resources. */
		public static Policy find(final FamilyResourceResolver resolver, final FamilyResourceResolver.Overrides overrides){
			require(resolver!=null && overrides!=null, "Caller policy lookup needs its actual group origin");
			final FamilyResourceResolver.Handle handle=resolver.optional(FAMILIES[0], "caller_policy", overrides);
			if(handle!=null){return new Policy(handle.path, resolver.root(FAMILIES[0], overrides));}
			for(final String family:FAMILIES){
				for(final String role:POLICY_ROLES){
					require(resolver.optional(family, role, overrides)==null,
						"New cobalamin resources require a complete caller_policy binding: "+family+'/'+role);
				}
			}
			return null;
		}

		/** Runs before main tip lookup; the covariance choice determines which payloads open. */
		void apply(final String family, final NcrnaFamilyLoader.Options options, final FamilyResourceResolver resolver){
			require(NcrnaFamilyRegistry.isCobalamin(family) && root.equals(resolver.root(family, options.resources)),
				"One serving policy requires all members from its own complete resource group");
			final FamilyResourceResolver.Overrides o=options.resources;
			require(o.consensus==null && o.models==null && o.seeds==null && o.cutoffs==null
				&& o.boundary==null && o.boundary3==null && o.boundaryStart==null && o.boundaryStop==null
				&& o.tables==null && o.endpoints==null && o.voteOffsets==null && o.joinedModelPositions==null,
				"Use a complete cobalamin resource group, not individual payloads mixed with a packaged caller policy");
			// Keep legacy working-directory names from shadowing this selected group.
			o.root=root;
			if(FAMILIES[0].equals(family)){
				resolver.resolve(family, "cm_input_padding", options.resources);
				require(options.cobalaminCmPolicy==CobalaminCmEndpointRefiner.OFF
					|| options.cobalaminCmPolicy==CobalaminCmEndpointRefiner.LEARNED && options.cobalaminCmGa,
					"The packaged main member uses learned covariance ends with GA hold");
				options.cobalaminCmPolicy=CobalaminCmEndpointRefiner.LEARNED;
				options.cobalaminCmGa=true;
			}else{
				require(options.cobalaminCmPolicy==CobalaminCmEndpointRefiner.OFF && !options.cobalaminCmGa,
					"Minor covariance models score acceptance only; their native endpoint routes remain active");
				if(FAMILIES[1].equals(family)){resolver.resolve(family, "positional_config", options.resources);}
			}
		}
		/** The unified loader validates the same percentile as the legacy member consumer. */
		int keyPercentile(){return percentile;}
		private final String root, pathCell;
		private final int bins, shortlist, patience, minimumHits, percentile, depth;
		private final int[] pads=new int[3], ceilings=new int[3];
		private static final String HEADER="bins\tshortlist\tpatience\tddl_min_hits\tkey_percentile\tmain_pad\tsecond_pad\tthird_pad\tmain_ceiling\tsecond_ceiling\tthird_ceiling\tpath_cell\tacceptance_floor\tacceptance_depth\tmain_cm_policy\tmain_cm_ga";
	}

	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private final Policy policy;
	private final NcrnaFamily[] families;
	private final NcrnaDdlProfileResource.ConsensusPair profiles;
	private final NcrnaOverlapWindowBuilder.Caps[] caps=new NcrnaOverlapWindowBuilder.Caps[3];
	private final NcrnaPathScoreConfig.Cell scoreCell;
	private final CobalaminCmAcceptance.Resources acceptance;
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String[] FAMILIES={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
	private static final float[] FLOORS={.63f, .66f, .66f};
	private static final String[] POLICY_ROLES={"ddl_profiles", "key_bounds", "path_policy", "acceptance_policy"};
}
