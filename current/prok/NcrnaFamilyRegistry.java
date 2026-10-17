package prok;

/** Stable resource identities and calibrated policies for migrated families.
 * Other families retain their existing factories until individually migrated.
 * @author Brian, Raiden */
public final class NcrnaFamilyRegistry {

	public static Entry get(String id){
		if("5S".equalsIgnoreCase(id) || "prok5S".equalsIgnoreCase(id)){return PROK5S;}
		if("23S".equalsIgnoreCase(id)){return PROK23S;}
		if("16S".equalsIgnoreCase(id)){return PROK16S;}
		if("euk5S".equalsIgnoreCase(id)){return EUK5S;}
		throw new IllegalArgumentException("No unified ncRNA resource registration for family: "+id);
	}

	public static final class Entry {
		private Entry(String id_, String directory_, String gate_, String[] models_,
				int minLen_, int windowPad_, int indexK_, int topN_, int minHits_,
				float weight_, float idPass_, float idBorderline_, float hbmPass_, boolean hasHbm_,
				int outputType_, boolean pass2_, boolean ranked_, boolean strict_,
				int partialMinLen_, int cap_, int voteSlack_, int boundaryRadius_,
				boolean preserveContaining_, String[] rawEndpointModels_,
				float collapse_, boolean modelCutoffs_, boolean dynamicModels_, int minSeedKeys_, String privateLegacyConsensus_){
			assert(id_!=null && directory_!=null && models_.length>0) : "Resource identity and ordered models bind every loaded family";
			id=id_;directory=directory_;gate=gate_;models=models_;
			minLen=minLen_;windowPad=windowPad_;indexK=indexK_;topN=topN_;minHits=minHits_;
			weight=weight_;idPass=idPass_;idBorderline=idBorderline_;hbmPass=hbmPass_;hasHbm=hasHbm_;
			outputType=outputType_;pass2=pass2_;rankedFallback=ranked_;strictIndexCutoff=strict_;
			partialMinLen=partialMinLen_;quantumCap=cap_;voteSlack=voteSlack_;boundaryRadius=boundaryRadius_;
			preserveContaining=preserveContaining_;rawEndpointModels=rawEndpointModels_;
			collapse=collapse_;modelCutoffs=modelCutoffs_;dynamicModels=dynamicModels_;
			minSeedKeys=minSeedKeys_;
			privateLegacyConsensus=privateLegacyConsensus_;
		}
		public final String id, directory, gate;
		final String privateLegacyConsensus;
		public final String member="";
		final String[] models;
		final String[] rawEndpointModels;
		final String[] coreRoles={"consensus", "seeds"};
		final String[] endpointRoles={"tables", "endpoints"};
		final int seedK=17, seedMinHits=2;
		final int minSeedKeys;
		final int minLen, windowPad, indexK, topN, minHits, outputType;
		final int partialMinLen, quantumCap, voteSlack, boundaryRadius;
		final float scoreA=0f;
		final float collapse;
		final float weight, idPass, idBorderline, hbmPass;
		final boolean hasHbm, pass2, rankedFallback, strictIndexCutoff;
		final boolean preserveContaining;
		/** An explicit custom consensus may change the roster and default shortlist size. */
		final boolean modelCutoffs, dynamicModels;
		int[] boundaryOffsets(boolean start){
			if(boundaryRadius<0){return start ? NcrnaFamily.LEGACY_START_OFFSETS : NcrnaFamily.LEGACY_STOP_OFFSETS;}
			assert(boundaryRadius>0) : "Centered boundary offsets must bracket the accepted alignment";
			final int[] offsets=new int[2*boundaryRadius+1];
			for(int i=0; i<offsets.length; i++){offsets[i]=i-boundaryRadius;}return offsets;
		}
	}
	private static final Entry PROK5S=new Entry("prok5S", "callgenes_5S", "prok5S seeded selection",
		new String[]{"universal", "bacteria", "archaea", "plastid"},
		90, 150, 7, 3, 2, 3f, .80f, .60f, .60f, true, ProkObject.r5S,
		true, false, false, 0, 1, 60, -1, false, new String[0], .85f, false, false, 1, null);
	private static final Entry PROK23S=new Entry("23S", "callgenes_23S", "23S seeded selection",
		new String[]{"universal", "bacteria", "archaea", "plastid"},
		2600, 3500, 9, 3, 128, 1f, .60f, .60f, 1.01f, false, ProkObject.r23S,
		false, true, true, 1800, 2, 150, 16, false, new String[0], .85f, false, false, 1, null);
	private static final Entry PROK16S=new Entry("16S", "callgenes_16S", "16S seeded selection",
		new String[]{"universal", "bacteria", "archaea", "plastid", "mitochondria"},
		1400, 1800, 9, 1, 1, 1f, .66f, .66f, 1.01f, false, ProkObject.r16S,
		true, false, false, 1200, 1, 60, -1, true, new String[]{"mitochondria"}, .85f, false, false, 1, null);
	private static final Entry EUK5S=new Entry("euk5S", "callgenes_euk5S", "euk5S selected calling",
		new String[]{"euk5S_universal", "euk5S_fungi", "euk5S_plant", "euk5S_animal", "euk5S_protist",
			"euk5S_residual_protozoa", "euk5S_residual_invertebrate"},
		90, 150, 7, 7, 2, 3f, .68f, .68f, 1.01f, false, ProkObject.r5S,
		false, false, false, 0, 1, 60, -1, false, new String[0], .75f, true, true, 2, "euk5S_consensus_sequence.fa");
	private NcrnaFamilyRegistry(){}
}
