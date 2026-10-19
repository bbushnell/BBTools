package prok;

/** Stable resource identities and calibrated policies for migrated families.
 * Other families retain their existing factories until individually migrated.
 * @author Brian, Raiden */
public final class NcrnaFamilyRegistry {

	/** Preserve existing family floors; introducing edge rescue does not tune length. */
	public static void enablePhysicalEdges(NcrnaFamily family){
		final Entry entry=get(family.name);
		if("tRNA".equals(entry.id) || family.modelEndClipping){
			throw new IllegalArgumentException("Physical-edge policy requires a generic registered caller without primary model clipping: "+family.name);
		}
		family.physicalEdgeMinLen=family.prok16sPartialMinLen>0 ? family.prok16sPartialMinLen : family.minLen;
		// Physical18S rescue retains its measured cutoff when unrestricted interior
		// rescue is disabled. Primary model thresholds remain unchanged.
		if("euk18S".equals(entry.id) && Float.isNaN(family.modelClipRescueId)){
			family.modelClipRescueId=.66f;
		}
		assert(family.physicalEdgeMinLen>0 && family.physicalEdgeMinLen<=family.minLen)
			: "Existing partial floors or unchanged family minLen bound physical-edge rescue";
		System.err.println("Physical-edge policy: family="+family.name+" minLen="+family.physicalEdgeMinLen+" interiorModelRescue=false");
	}

	public static Entry get(String id){
		if("5S".equalsIgnoreCase(id) || "prok5S".equalsIgnoreCase(id)){return PROK5S;}
		if("23S".equalsIgnoreCase(id)){return PROK23S;}
		if("16S".equalsIgnoreCase(id)){return PROK16S;}
		if("euk5S".equalsIgnoreCase(id)){return EUK5S;}
		if("tRNA".equalsIgnoreCase(id)){return TRNA;}
		if("euk18S".equalsIgnoreCase(id)){return EUK18S;}
		if("rnasep".equalsIgnoreCase(id)){return RNASEP;}
		if("tmrna".equalsIgnoreCase(id)){return TMRNA;}
		if("r58".equalsIgnoreCase(id) || "5.8S".equalsIgnoreCase(id)){return R58;}
		if("lsu".equalsIgnoreCase(id) || "eukLSU".equalsIgnoreCase(id)){return EUKLSU;}
		if("sixs_rf00013".equalsIgnoreCase(id)){return SIXS_RF00013;}
		if("sixs_rf01685".equalsIgnoreCase(id)){return SIXS_RF01685;}
		if("srp_small".equalsIgnoreCase(id)){return SRP_SMALL;}
		if("srp_large".equalsIgnoreCase(id)){return SRP_LARGE;}
		throw new IllegalArgumentException("No unified ncRNA resource registration for family: "+id);
	}

	public static final class Entry {
		private Entry(String id_, String directory_, String gate_, String[] models_,
				int minLen_, int windowPad_, int indexK_, int topN_, int minHits_,
				float weight_, float idPass_, float idBorderline_, float hbmPass_, boolean hasHbm_,
				int outputType_, boolean pass2_, boolean ranked_, boolean strict_,
				int partialMinLen_, int cap_, int voteSlack_, int boundaryRadius_,
				boolean preserveContaining_, String[] rawEndpointModels_,
				float collapse_, boolean modelCutoffs_, boolean dynamicModels_, int minSeedKeys_, String privateLegacyConsensus_, BoundaryPolicy boundary_, Group group_, String member_, String legacyRole_){
			this(id_, directory_, gate_, models_, minLen_, windowPad_, indexK_, topN_, minHits_,
				weight_, idPass_, idBorderline_, hbmPass_, hasHbm_, outputType_, pass2_, ranked_, strict_,
				partialMinLen_, cap_, voteSlack_, boundaryRadius_, preserveContaining_, rawEndpointModels_,
				collapse_, modelCutoffs_, dynamicModels_, minSeedKeys_, privateLegacyConsensus_, boundary_, group_, member_, legacyRole_, null);
		}
		private Entry(String id_, String directory_, String gate_, String[] models_,
				int minLen_, int windowPad_, int indexK_, int topN_, int minHits_,
				float weight_, float idPass_, float idBorderline_, float hbmPass_, boolean hasHbm_,
				int outputType_, boolean pass2_, boolean ranked_, boolean strict_,
				int partialMinLen_, int cap_, int voteSlack_, int boundaryRadius_,
				boolean preserveContaining_, String[] rawEndpointModels_,
				float collapse_, boolean modelCutoffs_, boolean dynamicModels_, int minSeedKeys_, String privateLegacyConsensus_, BoundaryPolicy boundary_, Group group_, String member_, String legacyRole_, ResourcePolicy resources_){
			assert(id_!=null && directory_!=null && models_!=null && (models_.length>0 || boundary_!=null))
				: "Fixed rosters bind rRNA resources; shared-boundary families derive their aligned roster from the library";
			id=id_;directory=directory_;gate=gate_;models=models_;
			assert(member_!=null && (group_==null ? member_.isEmpty() : group_.hasMember(member_) && directory_.equals(group_.directory)))
				: "Grouped entries must name a declared member and the group's one resource root";
			group=group_;member=member_;
			minLen=minLen_;windowPad=windowPad_;indexK=indexK_;topN=topN_;minHits=minHits_;
			weight=weight_;idPass=idPass_;idBorderline=idBorderline_;hbmPass=hbmPass_;hasHbm=hasHbm_;
			outputType=outputType_;pass2=pass2_;rankedFallback=ranked_;strictIndexCutoff=strict_;
			partialMinLen=partialMinLen_;quantumCap=cap_;voteSlack=voteSlack_;boundaryRadius=boundaryRadius_;
			preserveContaining=preserveContaining_;rawEndpointModels=rawEndpointModels_;
			collapse=collapse_;modelCutoffs=modelCutoffs_;dynamicModels=dynamicModels_;
			minSeedKeys=minSeedKeys_;
			privateLegacyConsensus=privateLegacyConsensus_;
			assert(legacyRole_!=null) : "Legacy consumers must explicitly select their original reference role";
			legacyConsensusRole=legacyRole_;
			resources=resources_;packagedEndpoints=resources_==null ? boundary_==null : resources_.packagedEndpoints;
			boundary=boundary_;seedMinHits=boundary==null ? 2 : 1;scoreA=boundary==null ? 0f : boundary.scoreA;
		}
		public final String id, directory, gate;
		final String privateLegacyConsensus;
		final String legacyConsensusRole;
		final ResourcePolicy resources;
		final boolean packagedEndpoints;
		public final String member;
		public final Group group;
		final String[] models;
		final String[] rawEndpointModels;
		final String[] coreRoles={"consensus", "seeds"};
		final String[] endpointRoles={"tables", "endpoints"};
		final int seedK=17, seedMinHits;
		final int minSeedKeys;
		final int minLen, windowPad, indexK, topN, minHits, outputType;
		final int partialMinLen, quantumCap, voteSlack, boundaryRadius;
		final float scoreA;
		final BoundaryPolicy boundary;
		final float collapse;
		final float weight, idPass, idBorderline, hbmPass;
		final boolean hasHbm, pass2, rankedFallback, strictIndexCutoff;
		final boolean preserveContaining;
		/** An explicit custom consensus may change the roster and default shortlist size. */
		final boolean modelCutoffs, dynamicModels;
		int[] boundaryOffsets(boolean start){
			if(id.equals("r58")){return R58BoundaryConfig.defaultOffsets(start);}
			if(boundary!=null && boundary.trna){return start ? TrnaBoundaryVectorGen.START_OFFSETS : TrnaBoundaryVectorGen.STOP_OFFSETS;}
			if(boundary!=null){return start ? NcrnaFamily.LEGACY_START_OFFSETS : boundary.stopOffsets;}
			if(boundaryRadius<0){return start ? NcrnaFamily.LEGACY_START_OFFSETS : NcrnaFamily.LEGACY_STOP_OFFSETS;}
			assert(boundaryRadius>0) : "Centered boundary offsets must bracket the accepted alignment";
			final int[] offsets=new int[2*boundaryRadius+1];
			for(int i=0; i<offsets.length; i++){offsets[i]=i-boundaryRadius;}return offsets;
		}
	}
	private static final Entry PROK5S=new Entry("prok5S", "callgenes_5S", "prok5S seeded selection",
		new String[]{"universal", "bacteria", "archaea", "plastid"},
		90, 150, 7, 3, 2, 3f, .80f, .60f, .60f, true, ProkObject.r5S,
		true, false, false, 0, 1, 60, -1, false, new String[0], .85f, false, false, 1, null, null, null, "", "consensus");
	private static final Entry PROK23S=new Entry("23S", "callgenes_23S", "23S seeded selection",
		new String[]{"universal", "bacteria", "archaea", "plastid"},
		2600, 3500, 9, 3, 128, 1f, .60f, .60f, 1.01f, false, ProkObject.r23S,
		false, true, true, 1800, 2, 150, 16, false, new String[0], .85f, false, false, 1, null, null, null, "", "consensus");
	private static final Entry PROK16S=new Entry("16S", "callgenes_16S", "16S seeded selection",
		new String[]{"universal", "bacteria", "archaea", "plastid", "mitochondria"},
		1400, 1800, 9, 1, 1, 1f, .66f, .66f, 1.01f, false, ProkObject.r16S,
		true, false, false, 1200, 1, 60, -1, true, new String[]{"mitochondria"}, .85f, false, false, 1, null, null, null, "", "consensus");
	private static final Entry EUK5S=new Entry("euk5S", "callgenes_euk5S", "euk5S selected calling",
		new String[]{"euk5S_universal", "euk5S_fungi", "euk5S_plant", "euk5S_animal", "euk5S_protist",
			"euk5S_residual_protozoa", "euk5S_residual_invertebrate"},
		90, 150, 7, 7, 2, 3f, .68f, .68f, 1.01f, false, ProkObject.r5S,
		false, false, false, 0, 1, 60, -1, false, new String[0], .75f, true, true, 2, "euk5S_consensus_sequence.fa", null, null, "", "consensus");
	private static final Entry TRNA=new Entry("tRNA", "callgenes_tRNA", "independent tRNA caller",
		new String[0], 40, 83, 7, 30, 12, 1f, .75f, .65f, .75f, true, ProkObject.tRNA,
		true, false, false, 0, 1, 60, -1, false, new String[0], .9f, false, true, 1, null,
		new BoundaryPolicy(null, null, 0f, true, 8f, .48f, .060f, TrnaBoundaryVectorGen.STOP_OFFSETS, true),
		null, "", "legacy_consensus");
	private static final Entry EUK18S=new Entry("euk18S", "callgenes_18S", "euk18s=t",
		new String[]{"euk18S_universal", "euk18S_fungi", "euk18S_plant", "euk18S_invertebrate",
			"euk18S_vertebrate_other", "euk18S_vertebrate_mammalian", "euk18S_mitochondria",
			"euk18S_protist", "euk18S_protozoa_component_001", "euk18S_protozoa_component_003"},
		916, 3109, 9, 10, 96, 1f, 1f, 1f, 1.01f, false, ProkObject.r18S,
		false, true, true, 0, 1, 60, -1, false, new String[0], .9f, true, true, 2, null, null, null, "", "legacy_consensus",
		new ResourcePolicy(false,
			new String[]{"consensus", "seeds", "cutoffs", "vote_offsets", "joined_model_positions"},
			new String[]{"euk18S_consensus_sequence.fa", "euk18S_17mers.fa.gz", "euk18S_model_cutoffs.tsv",
				"euk18S_seed_offsets.tsv", "euk18S_seed_positions.tsv"}));
	private static final Entry RNASEP=new Entry("rnasep", "callgenes_RNaseP", "ncrna=t",
		new String[0], 200, 320, 7, 60, 12, 1.6f, .65f, .58f, .62f, true, ProkObject.RNA,
		true, false, false, 0, 1, 60, -1, false, new String[0], .85f, false, true, 1, null,
		new BoundaryPolicy("rnasep_consensus.fa", "rnasep_models.hbm", 0f, true, 100f, .74f, .072f, new int[]{-1,0,1,2,3,4}), null, "", "consensus");
	private static final Entry TMRNA=new Entry("tmrna", "callgenes_tmRNA", "ncrna=t tmrna=t",
		new String[0], 180, 370, 7, 60, 12, 3f, .62f, .60f, .62f, true, ProkObject.RNA,
		true, false, false, 0, 1, 60, -1, false, new String[0], .85f, false, true, 1, null,
		new BoundaryPolicy("tmrna_consensus.fa", "tmrna_models.hbm", 5f, true, 11f, .48f, .072f, new int[]{-3,-2,-1,0,1,2}), null, "", "consensus");
	private static final Entry R58=new Entry("r58", "callgenes_5.8S", "ncrna=t r58lsu=t",
		new String[0], 140, 135, 7, 1, 2, 1f, .60f, .60f, .60f, true, ProkObject.RNA,
		true, false, false, 0, 1, 60, -1, false, new String[0], .85f, false, true, 1, null,
		new BoundaryPolicy("r58_consensus.fa", "r58_models.hbm", 0f, false, 0f, 0f, 0f, new int[]{-3,-2,-1,0,1,2}), null, "", "consensus");
	private static final Entry EUKLSU=new Entry("lsu", "callgenes_eukLSU", "ncrna=t r58lsu=t",
		new String[0], 60, 3500, 9, 4, 448, 1f, .60f, .54f, .54f, true, ProkObject.RNA,
		true, false, false, 0, 1, 150, -1, false, new String[0], .85f, false, true, 1, null,
		new BoundaryPolicy("lsu_consensus.fa", "lsu_models.hbm", 0f, false, 0f, 0f, 0f, new int[]{-4,-3,-2,-1,0,1,2,3,4}), null, "", "consensus");
	private static final Group SIXS_GROUP=new Group("6S", "callgenes_6S", new String[]{"RF00013", "RF01685"});
	private static final Entry SIXS_RF00013=new Entry("sixs_rf00013", "callgenes_6S", "ncrna=t sixs=t",
		new String[0], 60, 250, 7, 60, 1, 20f, .70f, .60f, .64f, true, ProkObject.RNA,
		true, false, false, 0, 1, 60, -1, false, new String[0], .85f, false, true, 1, null,
		new BoundaryPolicy("sixs_rf00013_consensus.fa", "sixs_rf00013_models.hbm", 0f, false, 0f, 0f, 0f, new int[]{-3,-2,-1,0,1,2}), SIXS_GROUP, "RF00013", "consensus");
	private static final Entry SIXS_RF01685=new Entry("sixs_rf01685", "callgenes_6S", "ncrna=t sixs=t",
		new String[0], 60, 125, 7, 60, 1, 20f, .70f, .70f, .75f, true, ProkObject.RNA,
		true, false, false, 0, 1, 60, -1, false, new String[0], .85f, false, true, 1, null,
		new BoundaryPolicy("sixs_rf01685_consensus.fa", "sixs_rf01685_models.hbm", 0f, false, 0f, 0f, 0f, new int[]{-3,-2,-1,0,1,2}), SIXS_GROUP, "RF01685", "consensus");
	private static final Group SRP_GROUP=new Group("SRP", "callgenes_SRP", new String[]{"small", "large"},
		new String[]{"consensus", "models"}, new String[]{"seeds"}, new String[]{"boundary", "boundary_start", "boundary_stop"});
	private static final Entry SRP_SMALL=new Entry("srp_small", "callgenes_SRP", "ncrna=t",
		new String[0], 80, 75, 7, 60, 12, 1.8f, .70f, .64f, .70f, true, ProkObject.RNA,
		true, false, false, 0, 1, 60, -1, false, new String[0], .85f, false, true, 1, null,
		new BoundaryPolicy("srp_small_consensus.fa", "srp_small_models.hbm", 0f, true, 23f, .68f, .072f, new int[]{-2,-1,0,1,2,3}), SRP_GROUP, "small", "consensus");
	private static final Entry SRP_LARGE=new Entry("srp_large", "callgenes_SRP", "ncrna=t",
		new String[0], 150, 250, 7, 60, 12, 1f, .70f, .64f, .66f, true, ProkObject.RNA,
		true, false, false, 0, 1, 60, -1, false, new String[0], .85f, false, true, 1, null,
		new BoundaryPolicy("srp_large_consensus.fa", "srp_large_models.hbm", 0f, true, 45f, .75f, .072f, new int[]{-2,-1,0,1,2,3}), SRP_GROUP, "large", "consensus");
	/** Optional lookup names and endpoint availability; legacy role belongs only to Entry. */
	static final class ResourcePolicy {
		ResourcePolicy(boolean endpoints_, String[] roles_, String[] names_){
			assert(roles_!=null && names_!=null && roles_.length==names_.length)
				: "Resolver compatibility lookup requires a name for every declared role";
			packagedEndpoints=endpoints_;roles=roles_.clone();names=names_.clone();
		}
		String cwdName(String role){
			assert(role!=null) : "Compatibility lookup must select one explicit resource role";
			for(int i=0; i<roles.length; i++){if(roles[i].equals(role)){return names[i];}}return null;
		}
		final boolean packagedEndpoints;
		private final String[] roles, names;
	}
	/** One manifest binds every required member and any explicitly shared root payload. */
	static final class Group {
		Group(String id_, String directory_, String[] members_){
			this(id_, directory_, members_, new String[]{"consensus", "models", "seeds"}, new String[0], new String[0]);
		}
		Group(String id_, String directory_, String[] members_, String[] required_, String[] shared_, String[] optional_){
			assert(id_!=null && directory_!=null && members_!=null && members_.length>1)
				: "A group is a declared resource root containing multiple distinct caller members";
			assert(required_!=null && required_.length>0 && shared_!=null && optional_!=null)
				: "Manifest validation needs explicit required, shared and optional role lists";
			id=id_;directory=directory_;members=members_.clone();
			requiredRoles=required_.clone();sharedRoles=shared_.clone();optionalRoles=optional_.clone();
			for(String role:sharedRoles){
				assert(!contains(requiredRoles, role) && !contains(optionalRoles, role))
					: "A shared role has exactly one root binding, never independent member defaults: "+role;
			}
			for(int i=0; i<members.length; i++){
				FamilyResourceManifest.validateRelative(members[i]);
				assert(members[i].indexOf('/')<0) : "Resolver member paths must be one directory component";
				for(int j=0; j<i; j++){assert(!members[j].equals(members[i])) : "Group members must identify distinct subdirectories";}
			}
		}
		boolean hasMember(String member){
			assert(member!=null) : "Manifest member matching needs an explicit label";
			for(String value:members){if(value.equals(member)){return true;}}return false;
		}
		boolean shared(String role){return contains(sharedRoles, role);}
		boolean permits(String member, String role){
			assert(role!=null) : "Grouped roles are constrained by the package's declared contract";
			return member.equals(".") ? shared(role) : hasMember(member) && (contains(requiredRoles, role) || contains(optionalRoles, role));
		}
		private static boolean contains(String[] values, String role){
			assert(values!=null && role!=null) : "Role matching requires a complete registration contract";
			for(String value:values){if(value.equals(role)){return true;}}return false;
		}
		final String id, directory;
		final String[] members;
		final String[] requiredRoles, sharedRoles, optionalRoles;
	}
	/** Shared-end families retain their measured shortlist policy and legacy CWD overrides. */
	static final class BoundaryPolicy {
		BoundaryPolicy(String consensus_, String models_, float scoreA_, boolean adaptive_, float floor_, float topFrac_, float qFrac_, int[] stopOffsets_){
			this(consensus_, models_, scoreA_, adaptive_, floor_, topFrac_, qFrac_, stopOffsets_, false);
		}
		BoundaryPolicy(String consensus_, String models_, float scoreA_, boolean adaptive_, float floor_, float topFrac_, float qFrac_, int[] stopOffsets_, boolean trna_){
			assert((trna_ || consensus_!=null && models_!=null) && stopOffsets_!=null)
				: "Generic factories honor CWD names; tRNA uses resource roots only and retains its own boundary adapter";
			consensus=consensus_;models=models_;scoreA=scoreA_;floor=floor_;topFrac=topFrac_;stopOffsets=stopOffsets_;
			adaptive=adaptive_;qFrac=qFrac_;trna=trna_;
		}
		final String consensus, models;
		final boolean adaptive, trna;
		final float scoreA, floor, topFrac, qFrac;
		final int[] stopOffsets;
	}
	private NcrnaFamilyRegistry(){}
}
