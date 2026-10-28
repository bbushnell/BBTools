package prok;

import java.util.ArrayList;
import consensus.BaseGraph;
import fileIO.ByteFile;
import map.LongHashSet;
import parse.LineParser1;

/** Constructs the existing caller container without changing verification algorithms.
 * CallGenes uses core-only loading then its established endpoint configuration phase.
 * @author Brian, Raiden */
public final class NcrnaFamilyLoader {

	public static NcrnaFamily load(String familyId){return load(familyId, new Options(), new FamilyResourceResolver());}
	public static NcrnaFamily load(String familyId, Options options, FamilyResourceResolver resolver){
		require(options!=null && resolver!=null, "Explicit load options and resolver required");
		if(!options.enabled){return null;}//A disabled selection never resolves or opens a resource.
		final NcrnaFamilyRegistry.Entry p=NcrnaFamilyRegistry.get(familyId);
		require(options.collapsedWindowPad>=-1, "Collapsed-window padding uses -1 only for the registry default");
		final int seedK=options.seedK<0 ? p.seedK : options.seedK;
		require(seedK==p.seedK || p.boundary!=null && p.boundary.trna,
			"Only explicit legacy tRNA calling may select a different seed K: "+p.id);
		final String lib=resolver.resolve(p.id, "consensus", options.resources).path;
		final byte[][] library=TrnaConsensusBuilder.loadLibrary(lib);
		final String[] names=TrnaConsensusBuilder.lastLoadedNames;
		final boolean expandedMito=MitoSsuModelPolicy.extended(p.id, names);
		final boolean ordered=p.models.length>0 && (!p.dynamicModels || options.resources.consensus==null);
		require(library!=null && names!=null && library.length>0 && names.length==library.length
			&& (!ordered || library.length==p.models.length || expandedMito),
			"Incompatible "+p.id+" consensus resource: require its complete ordered model roster: "+lib);
		final map.ObjectSet<String> seen=new map.ObjectSet<String>(String.class);
		for(int i=0; i<names.length; i++){
			final String token=RrnaResourceIO.first(names[i]);
			require(!token.isEmpty() && seen.add(token) && (!ordered || (expandedMito && i>=p.models.length) || p.models[i].equals(token))
				&& library[i]!=null && library[i].length>=(p.dynamicModels ? p.seedK : 1),
				"Incompatible "+p.id+" consensus resource: model order, unique names or sequence length at "+i+": "+lib);
		}
		final BaseGraph[] models;
		require(p.boundary==null || options.hbm, "Selected "+p.id+" requires its index-aligned HBM library");
		final boolean useHbm=p.hasHbm && options.hbm;
		if(useHbm){
			final String hbm=resolver.resolve(p.id, "models", options.resources).path;
			models=TrnaConsensusBuilder.loadModels(hbm);
			CallGenes.requireNcrnaLibraryModelAlignment(p.id, library, models, names, lib, hbm);
		}else{models=null;}
		require(options.seeds || p.boundary!=null && p.boundary.trna,
			"Only the independent tRNA compatibility path may disable seed loading");
		final String seeds=options.seeds ? resolver.seeds(p.id, seedK, options.resources).path : null;
		final LongHashSet keys=options.seeds ? ProkObject.loadLongKmers(seeds, seedK) : null;
		require(!options.seeds || keys!=null && keys.size()>=p.minSeedKeys,
			"Incompatible "+p.id+" seed resource: require at least "+p.minSeedKeys+" distinct forward k-mers: "+seeds);
		final int k=options.indexK<0 ? p.indexK : options.indexK;
		final int top=options.topN<0 ? (p.dynamicModels && p.models.length>0 ? library.length : p.topN) : options.topN;
		final int hits=options.minHits<0 ? p.minHits : options.minHits;
		final boolean adaptive=options.adaptive==null ? p.boundary!=null && p.boundary.adaptive : options.adaptive;
		final float floor=Float.isNaN(options.floor) ? (p.boundary==null ? 0f : p.boundary.floor) : options.floor;
		final float topFrac=Float.isNaN(options.topFrac) ? (p.boundary==null ? 0f : p.boundary.topFrac) : options.topFrac;
		final float qFrac=Float.isNaN(options.qFrac) ? (p.boundary==null ? 0f : p.boundary.qFrac) : options.qFrac;
		TrnaKmerIndex.validateConfiguration(k, floor, topFrac, qFrac, hits);
		TrnaKmerIndex.validateTopN(top);
		final float weight=Float.isNaN(options.weight) ? p.weight : options.weight;
		final boolean memberIdentity=NcrnaFamilyRegistry.isCobalamin(p.id);
		final FamilyResourceResolver.Handle cutoff=memberIdentity && options.cutoffs ? resolver.optional(p.id, "cutoffs", options.resources) : null;
		final MemberCutoffs member=cutoff==null ? null : loadMemberCutoffs(cutoff.path, p.id);
		final float idPass=Float.isNaN(options.idPass) ? (member==null ? p.idPass : member.pass) : options.idPass;
		final float idBorderline=Float.isNaN(options.idBorderline) ? (member==null ? p.idBorderline : member.borderline) : options.idBorderline;
		final int seedMinHits=options.seedMinHits<0 ? p.seedMinHits : options.seedMinHits;
		require(Float.isFinite(idPass) && idPass>=0 && idPass<=1 && Float.isFinite(idBorderline)
			&& idBorderline>=0 && idBorderline<=1 && seedMinHits>0, "Invalid identity or seed gate for "+p.id);
		require(!memberIdentity || idBorderline<=idPass,
			"B12 borderline identity must not exceed pass identity after resource/CLI precedence: "+p.id);
		require(Float.isFinite(weight) && weight>0, "Path weight must be finite and positive for "+p.id);
		int pad=options.windowPad<0 ? p.windowPad : options.windowPad;
		if(options.longestPadFlank>=0){
			int longest=0;for(byte[] ref:library){longest=Math.max(longest, ref.length);}
			require(longest<=Integer.MAX_VALUE-options.longestPadFlank-p.seedK,
				"Longest-model padding must leave room for the flank and seed span: "+p.id);
			pad=longest+options.longestPadFlank;
		}
		final int minLen=options.minLen<0 ? p.minLen : options.minLen;
		final float scoreA=Float.isNaN(options.scoreA) ? p.scoreA : options.scoreA;
		final float scoreB=Float.isNaN(options.scoreB) ? weight : options.scoreB;
		final float collapse=Float.isNaN(options.collapse) ? p.collapse : options.collapse;
		final float hbmPass=useHbm ? (Float.isNaN(options.hbmPass) ? p.hbmPass : options.hbmPass) : 1.01f;
		require(!options.boundary || p.boundary!=null, "Shared boundary loading requires a registered shared-end adapter: "+p.id);
		require(!options.boundary || !p.boundary.trna,
			"Independent tRNA boundary templates must use loadTrnaBoundary, never the generic10-input adapter");
		final NcrnaBoundaryResources boundary=options.boundary ? new NcrnaBoundaryResources(p.id, library, options.resources, resolver) : null;
		final NcrnaFamily family=new NcrnaFamily(p.id, library, models, names, keys, seedK, minLen, pad,
			k, top, adaptive, floor, topFrac, qFrac, hits,
			scoreA, scoreB, idPass, idBorderline, hbmPass, collapse,
			boundary==null ? null : boundary.net, boundary==null ? null : boundary.net,
			boundary==null ? null : boundary.start.table, boundary==null ? null : boundary.stop.table,
			boundary==null ? -1 : boundary.start.insideCount, boundary==null ? -1 : boundary.start.outsideCount,
			boundary==null ? -1 : boundary.stop.insideCount, boundary==null ? -1 : boundary.stop.outsideCount,
			boundary==null ? 0f : boundary.meanLen,
			options.boundaryStartOffsets==null ? p.boundaryOffsets(true) : options.boundaryStartOffsets,
			options.boundaryStopOffsets==null ? p.boundaryOffsets(false) : options.boundaryStopOffsets,
			options.boundaryMarginStart, options.boundaryMarginStop);
		family.outputType=p.outputType;family.seedMinHits=seedMinHits;family.scavengePass2=p.pass2;
		if("fmn".equals(p.id)){family.seedDistinct=seedMinHits>1;}//The direct loader also enforces FMN's distinct-key policy.
		family.rankedModelFallback=p.rankedFallback;family.strictIndexCutoff=p.strictIndexCutoff;family.trimAlignmentExtent=p.boundary!=null;
		family.reuseConsensusAlignment=p.boundary==null;family.voteSlack=p.voteSlack;
		family.prok16sPartialMinLen=p.partialMinLen;family.quantumCapMultiplier=p.quantumCap;
		family.preserveContainingWindows=p.preserveContaining;
		family.collapsedWindowPad=options.collapsedWindowPad<0 ? p.collapsedWindowPad : options.collapsedWindowPad;
		require(family.collapsedWindowPad>=0 && (family.collapsedWindowPad==0
			|| 2L*family.collapsedWindowPad+seedK>=minLen),
			"Collapsed-window padding must be disabled or span the unchanged gene-length floor: "+p.id);
		if(p.modelCutoffs && options.cutoffs){loadCutoffs(family, options.resources, resolver);}
		if(MitoSsuModelPolicy.coverage(p.id,names)){
			NcrnaModelThresholds.loadLengths(resolver.resolve(p.id,"length_bounds",options.resources).path,family);
		}else if(options.resources.consensus==null || options.resources.root!=null || options.resources.lengthBounds!=null){
			final FamilyResourceResolver.Handle lengths=resolver.optional(p.id, "length_bounds", options.resources);
			if(lengths!=null){
				require(!"tRNA".equals(p.id), "The independent tRNA verifier does not consume generic model length bounds");
				NcrnaModelThresholds.loadLengths(lengths.path, family);
			}
		}
		if(options.endpoints && p.packagedEndpoints){family.setRrnaEndpointInference(loadEndpoints(family, options.resources, resolver), null);}
		return family;
	}
	static LongHashSet loadSeedKeys(String family, int k, FamilyResourceResolver resolver){
		require(resolver!=null, "Seed-only compatibility loading still requires the unified path resolver");
		return ProkObject.loadLongKmers(resolver.seeds(family, k, null).path, k);
	}
	/** Core and boundary selections are independent. Null core arrays are valid for
	 * the public boundary-only API; this carrier must never enter the generic caller list. */
	public static NcrnaFamily loadTrnaBoundary(byte[][] library, BaseGraph[] models, String[] names, LongHashSet keys, int seedK,
			FamilyResourceResolver.Overrides resources, FamilyResourceResolver resolver){
		require(resources!=null && resolver!=null, "tRNA boundary selection requires explicit override state and a resolver");
		final NcrnaFamilyRegistry.Entry p=NcrnaFamilyRegistry.get("tRNA");
		final String five=resolver.resolve(p.id, "boundary", resources).path;
		final boolean shared=resources.boundary3==null || resources.boundary3.equals(resources.boundary);
		String three=five;
		if(!shared){
			final FamilyResourceResolver.Overrides separate=new FamilyResourceResolver.Overrides();
			separate.boundary=resources.boundary3;
			three=resolver.resolve(p.id, "boundary", separate).path;
		}
		final TrnaBoundaryResources b=new TrnaBoundaryResources(five, three,
			resolver.resolve(p.id, "boundary_start", resources).path,
			resolver.resolve(p.id, "boundary_stop", resources).path, shared);
		return new NcrnaFamily(p.id, library, models, names, keys, seedK, p.minLen, p.windowPad,
			p.indexK, p.topN, p.boundary.adaptive, p.boundary.floor, p.boundary.topFrac, p.boundary.qFrac, p.minHits,
			p.scoreA, p.weight, p.idPass, p.idBorderline, p.hbmPass, p.collapse,
			b.net5, b.net3, b.start.table, b.stop.table, b.start.insideCount, b.start.outsideCount,
			b.stop.insideCount, b.stop.outsideCount, 0f, p.boundaryOffsets(true), p.boundaryOffsets(false), 0f, 0f);
	}
	static void loadVoteOffsets(NcrnaFamily family, FamilyResourceResolver.Overrides overrides,
			FamilyResourceResolver resolver, String expectedSeedSha80){
		require(family!=null && resolver!=null, "Vote offsets must bind a loaded family's exact seed set");
		family.voteTable=SeedOffsetTable.load(resolver.resolve(family.name, "vote_offsets", overrides).path,
			family.kLong, family.kmerSet, expectedSeedSha80);
	}
	static void loadJoinedPositions(NcrnaFamily family, FamilyResourceResolver.Overrides overrides,
			FamilyResourceResolver resolver){
		require(family!=null && resolver!=null, "Joined positions must bind loaded model names, sequences and seeds");
		family.joinedPositions=SeedModelPositionResource.load(resolver.resolve(family.name, "joined_model_positions", overrides).path,
			family.kmerSet, family.modelNames, family.library);
	}
	static RrnaEndpointCallerFeatures.Resources loadEndpoints(NcrnaFamily family,
			FamilyResourceResolver.Overrides overrides, FamilyResourceResolver resolver){
		require(family!=null, "Endpoint resources require an initialized registered family");
		final NcrnaFamilyRegistry.Entry p=NcrnaFamilyRegistry.get(family.name);
		final FamilyResourceResolver.Handle tables=resolver.tables(family.name, overrides);
		final RrnaEndpointCallerFeatures.Resources loaded=RrnaEndpointResourceLoader.load(tables, family.modelNames, family.library,
			MitoSsuModelPolicy.rawNames(family, p.rawEndpointModels));
		final RrnaEndpointCallerFeatures.Resources inference=RrnaEndpointNetworkLoader.load(resolver.resolve(family.name, "endpoints", overrides), loaded);
		require(inference.acceptedVerifierIdentity==RrnaEndpointCallerFeatures.Resources.supportsAcceptedIdentity(p.id),
			p.id+" nets must preserve their trained identity feature; euk5S uses candidate-span identity");
		return inference;
	}
	static void loadCutoffs(NcrnaFamily family, FamilyResourceResolver.Overrides overrides, FamilyResourceResolver resolver){
		require(family!=null && resolver!=null, "Model cutoffs require a loaded family and its resolver");
		final ArrayList<NcrnaFamily> one=new ArrayList<NcrnaFamily>(1);one.add(family);
		NcrnaModelThresholds.load(resolver.resolve(family.name, "cutoffs", overrides).path, one);
	}
	/** A generic member uses one scalar identity pair throughout its primary and
	 * rescue paths, unlike the existing single-alignment per-model cutoff tables. */
	private static MemberCutoffs loadMemberCutoffs(final String path, final String family){
		final ByteFile in=RrnaResourceIO.open(path);
		final LineParser1 parser=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, "family\tidpass\tidborderline", path);
			final byte[] row=in.nextLine();
			require(row!=null, "A declared member cutoff file requires one identity row: "+path);
			parser.set(row);
			require(parser.terms()==3 && parser.termEquals(family, 0),
				"Member cutoffs must bind exactly the loaded family: "+family+" in "+path);
			final float pass=(float)RrnaPositionalKmerTable.savedDouble(parser, 1);
			final float borderline=(float)RrnaPositionalKmerTable.savedDouble(parser, 2);
			require(Float.isFinite(pass) && Float.isFinite(borderline) && borderline>=0 && borderline<=pass && pass<=1,
				"Member cutoffs require finite 0<=idborderline<=idpass<=1: "+path);
			require(in.nextLine()==null, "One member cutoff file cannot silently override another family or repeat rows: "+path);
			return new MemberCutoffs(pass, borderline);
		}finally{require(!in.close(), "Member cutoff read failed: "+path);}
	}
	private static final class MemberCutoffs {
		MemberCutoffs(float pass_, float borderline_){pass=pass_;borderline=borderline_;}
		final float pass, borderline;
	}
	public static final class Options {
		public boolean enabled=true, hbm=true, seeds=true, endpoints=true, cutoffs=true, boundary=false;
		public Boolean adaptive=null;
		public int indexK=-1, topN=-1, minHits=-1, seedMinHits=-1, windowPad=-1, minLen=-1;
		public int seedK=-1;
		public int longestPadFlank=-1;
		public int collapsedWindowPad=-1;
		public int[] boundaryStartOffsets, boundaryStopOffsets;
		public float boundaryMarginStart=0f, boundaryMarginStop=0f;
		public float floor=Float.NaN, topFrac=Float.NaN, qFrac=Float.NaN, weight=Float.NaN;
		public float idPass=Float.NaN, idBorderline=Float.NaN, hbmPass=Float.NaN;
		public float scoreA=Float.NaN, scoreB=Float.NaN, collapse=Float.NaN;
		public final FamilyResourceResolver.Overrides resources=new FamilyResourceResolver.Overrides();
	}
	private static void require(boolean ok, String why){FamilyResourceManifest.require(ok, why);}
	private NcrnaFamilyLoader(){}
}
