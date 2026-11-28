package prok;

import java.util.ArrayList;
import consensus.BaseGraph;
import fileIO.ByteFile;
import fileIO.FileFormat;
import map.LongHashSet;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;

/** Constructs the existing caller container without changing verification algorithms.
 * CallGenes uses core-only loading then its established endpoint configuration phase.
 * @author Brian, Raiden
 */
public final class NcrnaFamilyLoader {

	public static NcrnaFamily load(String familyId){return load(familyId, new Options(), new FamilyResourceResolver());}
	public static NcrnaFamily load(String familyId, Options options, FamilyResourceResolver resolver){
		require(options!=null && resolver!=null, "Explicit load options and resolver required");
		if(!options.enabled){return null;}//A disabled selection never resolves or opens a resource.
		final NcrnaFamilyRegistry.Entry p=NcrnaFamilyRegistry.get(familyId);
		CobalaminCallerResources.Policy callerPolicy=null;
		if(NcrnaFamilyRegistry.isCobalamin(p.id)){
			callerPolicy=options.cobalaminCallerPolicy!=null ? options.cobalaminCallerPolicy
				: CobalaminCallerResources.Policy.find(resolver, options.resources);
			if(callerPolicy!=null){callerPolicy.apply(p.id, options, resolver);}
		}
		final boolean cm=options.cobalaminCmPolicy!=CobalaminCmEndpointRefiner.OFF;
		require(options.cobalaminCmPolicy>=CobalaminCmEndpointRefiner.OFF
			&& options.cobalaminCmPolicy<=CobalaminCmEndpointRefiner.POOLED && (!cm || NcrnaFamilyRegistry.isCobalamin(p.id)),
			"CM endpoints require an explicitly selected cobalamin member and known operating point");
		require(!options.cobalaminCmGa || cm, "A CM score hold requires an explicit CM endpoint mode");
		if(cm){
			require(!Boolean.TRUE.equals(options.boundaryOverride) && options.resources.boundary==null
				&& options.resources.boundary3==null && options.resources.boundaryStart==null && options.resources.boundaryStop==null
				&& options.resources.voteOffsets==null && options.resources.tables==null && options.resources.endpoints==null
				&& options.boundaryStartOffsets==null && options.boundaryStopOffsets==null
				&& !options.boundaryMarginStartSet && !options.boundaryMarginStopSet,
				"Select covariance endpoints or explicit tip/vote resources, never both");
			options.boundary=false;
			options.endpoints=false;
		}
		final int patience=options.modelPatience==Options.UNSET_PATIENCE ? p.modelPatience : options.modelPatience;
		NcrnaPatience.validate(patience);
		require(patience==-1 || p.boundary!=null && !p.boundary.trna,
			"Explicit ncRNA model patience requires a generic shared-boundary caller: "+p.id);
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
			require(!token.isEmpty() && seen.add(token)
				&& (!ordered || (expandedMito && i>=p.models.length) || p.models[i].equals(token))
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
		final int top=options.topN<0 ?
			(p.topN<0 || p.dynamicModels && p.models.length>0 ? library.length : p.topN) : options.topN;
		final int hits=options.minHits<0 ? p.minHits : options.minHits;
		final boolean adaptive=options.adaptive==null ? p.boundary!=null && p.boundary.adaptive : options.adaptive;
		final float floor=setNonnegative(options.floor, "floor") ? options.floor : (p.boundary==null ? 0f : p.boundary.floor);
		final float topFrac=setNonnegative(options.topFrac, "topFrac") ? options.topFrac : (p.boundary==null ? 0f : p.boundary.topFrac);
		final float qFrac=setNonnegative(options.qFrac, "qFrac") ? options.qFrac : (p.boundary==null ? 0f : p.boundary.qFrac);
		TrnaKmerIndex.validateConfiguration(k, floor, topFrac, qFrac, hits);
		TrnaKmerIndex.validateTopN(top);
		final float weight=setNonnegative(options.weight, "weight") ? options.weight : p.weight;
		final boolean memberIdentity=NcrnaFamilyRegistry.isCobalamin(p.id);
		final FamilyResourceResolver.Handle cutoff=memberIdentity && options.cutoffs ?
			resolver.optional(p.id, "cutoffs", options.resources) : null;
		final MemberCutoffs member=cutoff==null ? null : loadMemberCutoffs(cutoff.path, p.id);
		final float idPass=setNonnegative(options.idPass, "idPass") ? options.idPass : (member==null ? p.idPass : member.pass);
		final float idBorderline=setNonnegative(options.idBorderline, "idBorderline") ?
			options.idBorderline : (member==null ? p.idBorderline : member.borderline);
		final int seedMinHits=options.seedMinHits<0 ? p.seedMinHits : options.seedMinHits;
		require(Float.isFinite(idPass) && idPass>=0 && idPass<=1 && Float.isFinite(idBorderline)
			&& idBorderline>=0 && idBorderline<=1 && seedMinHits>0, "Invalid identity or seed gate for "+p.id);
		require(!memberIdentity || idBorderline<=idPass,
			"B12 borderline identity must not exceed pass identity after resource/CLI precedence: "+p.id);
		require(Float.isFinite(weight) && weight>0, "Path weight must be finite and positive for "+p.id);
		int pad=options.windowPad<0 ? p.windowPad : options.windowPad;
		// A negative registered pad derives the longest model length. Explicit
		// window padding still wins unless longestPadFlank was separately requested.
		final int longestFlank=options.longestPadFlank>=0 ? options.longestPadFlank : pad<0 ? 0 : -1;
		if(longestFlank>=0){
			int longest=0;for(byte[] ref:library){longest=Math.max(longest, ref.length);}
			require(longest<=Integer.MAX_VALUE-longestFlank-p.seedK,
				"Longest-model padding must leave room for the flank and seed span: "+p.id);
			pad=longest+longestFlank;
		}
		final int minLen=options.minLen<0 ? p.minLen : options.minLen;
		final float scoreA=options.scoreASet ? options.scoreA : p.scoreA;
		final float scoreB=options.scoreBSet ? options.scoreB : weight;
		require(Float.isFinite(scoreA) && Float.isFinite(scoreB),
			"Signed path scores must be finite for "+p.id);
		final float collapse=setNonnegative(options.collapse, "collapse") ? options.collapse : p.collapse;
		final float hbmPass=useHbm ? (setNonnegative(options.hbmPass, "hbmPass") ? options.hbmPass : p.hbmPass) : 1.01f;
		require(collapse>=0f && collapse<=1f, "Collapse fraction must be in [0,1] for "+p.id);
		PackagedBoundaryConfig packagedBoundary=null;
		PackagedPositionalBoundaryConfig positionalBoundary=null;
		if(!cm && p.boundary!=null && !p.boundary.trna){
			final FamilyResourceResolver.Handle config=resolver.optional(p.id, "boundary_config", options.resources);
			if(config!=null){
				packagedBoundary=PackagedBoundaryConfig.read(config.path);
				packagedBoundary.applyDefaults(p.id, options, library);
			}
			final FamilyResourceResolver.Handle positional=resolver.optional(p.id, "positional_config", options.resources);
			if(positional!=null){
				require(packagedBoundary==null, "One family cannot declare both sliding and positional boundary configurations: "+p.id);
				positionalBoundary=PackagedPositionalBoundaryConfig.read(positional.path);
				positionalBoundary.applyDefaults(p.id, options, library);
			}
		}
		require(!options.boundary || p.boundary!=null,
			"Shared boundary loading requires a registered shared-end adapter: "+p.id);
		require(!options.boundary || !p.boundary.trna,
			"Independent tRNA boundary templates must use loadTrnaBoundary, never the generic10-input adapter");
		final NcrnaBoundaryResources boundary=options.boundary ?
			new NcrnaBoundaryResources(p.id, library, options.resources, resolver) : null;
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
		family.modelPatience=patience;
		if("fmn".equals(p.id) || "sam1".equals(p.id) || "tpp".equals(p.id) || "glycine".equals(p.id)){
			family.seedDistinct=seedMinHits>1;//Direct loads preserve riboswitch distinct-key support.
		}
		family.rankedModelFallback=p.rankedFallback; family.strictIndexCutoff=p.strictIndexCutoff;
		family.trimAlignmentExtent=p.boundary!=null;
		family.reuseConsensusAlignment=p.boundary==null;family.voteSlack=p.voteSlack;
		if("glycine".equals(p.id)){
			// Glycine verifies every window with Quantum and retains its endpoints;
			// the legacy short-window and second endpoint alignments use Scrabble.
			family.quantumThresh=0;
			family.reuseConsensusAlignment=true;
			family.trimAlignmentExtent=false;
			// Shared endpoint nets must run on that retained Quantum anchor;
			// finishAcceptedBoundary otherwise requires a successful second trim.
			family.boundaryOnRawEndpoints=true;
		}
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
		if(options.endpoints && p.packagedEndpoints){
			family.setRrnaEndpointInference(loadEndpoints(family, options.resources, resolver), null);
		}
		if(packagedBoundary!=null){packagedBoundary.loaded(family);}
		if(positionalBoundary!=null){positionalBoundary.loaded(family, resolver);}
		if(cm){
			family.boundaryRefinerFactory=CobalaminCmEndpointRefiner.load(family, resolver, options.resources, options.cobalaminCmPolicy, options.cobalaminCmGa);
			family.boundaryOnRawEndpoints=true;
		}
		// Library-only consumers intentionally load no seed-dependent serving resources.
		final FamilyResourceResolver.Handle bounds=options.seeds ? resolver.optionalKeyBounds(p.id, options.resources) : null;
		if(bounds!=null){
			require(seedK==17 && keys!=null, "Per-key windows require the actual loaded forward K17 seed axis: "+p.id);
			final String boundsMember=p.group==null ? p.id : p.member;
			final NcrnaOverlapWindowBuilder.Caps caps=memberIdentity ?
				NcrnaOverlapWindowBuilder.Caps.load(bounds.path, boundsMember, keys.size(), callerPolicy==null ? 99 : callerPolicy.keyPercentile())
				: NcrnaOverlapWindowBuilder.Caps.loadGeneric(bounds.path, boundsMember, keys.size());
			caps.validateKeys(keys.toArray());
			if(NcrnaFamilyRegistry.isCobalamin(p.id)){
				CobalaminCallerResources.requireLegacyBounds(caps);
				require(callerPolicy!=null, "Legacy cobalamin bounds require their existing external caller policy: "+p.id);
			}
			family.keyBounds=caps;
		}
		return family;
	}
	/** Bind sibling observations once before workers share these read-only families.
	 * Standalone families and the independent legacy cobalamin route stay unchanged. */
	public static void attributeMemberKeys(final ArrayList<NcrnaFamily> families){
		require(families!=null, "Member attribution must bind the actual loaded family population before workers start");
		for(final NcrnaFamily family:families){
			require(family!=null, "A loaded family list cannot contain an absent member");
			if(family.keyBounds==null){continue;}
			final NcrnaFamilyRegistry.Entry entry=NcrnaFamilyRegistry.get(family.name);
			if(entry.group==null || NcrnaFamilyRegistry.isCobalamin(entry.id)){continue;}
			NcrnaOverlapWindowBuilder.Caps caps=family.keyBounds;
			for(final NcrnaFamily sibling:families){
				require(sibling!=null, "A loaded sibling list cannot contain an absent member");
				if(sibling==family || sibling.keyBounds==null){continue;}
				if(NcrnaFamilyRegistry.get(sibling.name).group==entry.group){
					caps=caps.withoutSiblingOnly(sibling.keyBounds);
				}
			}
			assert(caps.unmeasuredKeys==family.keyBounds.unmeasuredKeys
				&& caps.chainCeiling==family.keyBounds.chainCeiling)
				: "Sibling attribution must preserve the construction population and chain ceiling used by window geometry";
			family.keyBounds=caps;
		}
	}
	static LongHashSet loadSeedKeys(String family, int k, FamilyResourceResolver resolver){
		require(resolver!=null, "Seed-only compatibility loading still requires the unified path resolver");
		return ProkObject.loadLongKmers(resolver.seeds(family, k, null).path, k);
	}
	/** Core and boundary selections are independent. Null core arrays are valid for
	 * the public boundary-only API; this carrier must never enter the generic caller list. */
	public static NcrnaFamily loadTrnaBoundary(byte[][] library, BaseGraph[] models, String[] names,
			LongHashSet keys, int seedK,
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
		family.joinedPositions=SeedModelPositionResource.load(
			resolver.resolve(family.name, "joined_model_positions", overrides).path,
			family.kmerSet, family.modelNames, family.library);
	}
	static RrnaEndpointCallerFeatures.Resources loadEndpoints(NcrnaFamily family,
			FamilyResourceResolver.Overrides overrides, FamilyResourceResolver resolver){
		require(family!=null, "Endpoint resources require an initialized registered family");
		final NcrnaFamilyRegistry.Entry p=NcrnaFamilyRegistry.get(family.name);
		final FamilyResourceResolver.Handle tables=resolver.tables(family.name, overrides);
		final RrnaEndpointCallerFeatures.Resources loaded=RrnaEndpointResourceLoader.load(tables,
			family.modelNames, family.library,
			MitoSsuModelPolicy.rawNames(family, p.rawEndpointModels));
		final RrnaEndpointCallerFeatures.Resources inference=RrnaEndpointNetworkLoader.load(
			resolver.resolve(family.name, "endpoints", overrides), loaded);
		require(inference.acceptedVerifierIdentity==RrnaEndpointCallerFeatures.Resources.supportsAcceptedIdentity(p.id),
			p.id+" nets must preserve their trained identity feature; euk5S uses candidate-span identity");
		return inference;
	}
	static void loadCutoffs(NcrnaFamily family, FamilyResourceResolver.Overrides overrides,
			FamilyResourceResolver resolver){
		require(family!=null && resolver!=null, "Model cutoffs require a loaded family and its resolver");
		final ArrayList<NcrnaFamily> one=new ArrayList<NcrnaFamily>(1);one.add(family);
		NcrnaModelThresholds.load(resolver.resolve(family.name, "cutoffs", overrides).path, one);
	}
	/** A generic member uses one scalar identity pair throughout its primary and
	 * rescue paths, unlike the existing single-alignment per-model cutoff tables. */
	private static MemberCutoffs loadMemberCutoffs(final String path, final String family){
		final LineParser1 parser=new LineParser1('\t');
		try(Input input=new Input(FileFormat.testInput(path, FileFormat.TEXT, null, true, true))){
			final ByteFile in=input.reader;
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
		}
	}
	private static final class MemberCutoffs {
		MemberCutoffs(float pass_, float borderline_){pass=pass_;borderline=borderline_;}
		final float pass, borderline;
	}
	public static final class Options {
		/** Distinguishes an omitted registry choice from an explicit legacy value -1. */
		public static final int UNSET_PATIENCE=Integer.MIN_VALUE;
		public int modelPatience=UNSET_PATIENCE;
		/** Zero leaves the entire ordinary endpoint path and loading unchanged. */
		public int cobalaminCmPolicy=CobalaminCmEndpointRefiner.OFF;
		public boolean cobalaminCmGa=false;
		/** The group loader shares one cold policy parse among its three members. */
		public CobalaminCallerResources.Policy cobalaminCallerPolicy;
		/** Explicit enable selection; null permits a packaged default. Absent packages keep boundary unchanged. */
		public Boolean boundaryOverride;
		public boolean enabled=true, hbm=true, seeds=true, endpoints=true, cutoffs=true, boundary=false;
		public Boolean adaptive=null;
		public int indexK=-1, topN=-1, minHits=-1, seedMinHits=-1, windowPad=-1, minLen=-1;
		public int seedK=-1;
		public int longestPadFlank=-1;
		public int collapsedWindowPad=-1;
		public int[] boundaryStartOffsets, boundaryStopOffsets;
		/** Effective margins; unset ends inherit packaged defaults, or zero without a package. */
		public float boundaryMarginStart=0f, boundaryMarginStop=0f;
		/** Explicit per-end selections, including zero, take precedence over packaged defaults. */
		public boolean boundaryMarginStartSet=false, boundaryMarginStopSet=false;
		public float floor=-1f, topFrac=-1f, qFrac=-1f, weight=-1f;
		public float idPass=-1f, idBorderline=-1f, hbmPass=-1f;
		public float scoreA=0f, scoreB=0f, collapse=-1f;
		public boolean scoreASet=false, scoreBSet=false;
		public final FamilyResourceResolver.Overrides resources=new FamilyResourceResolver.Overrides();
	}
	private static boolean setNonnegative(float value, String label){
		require(Float.isFinite(value) && (value==-1f || value>=0f),
			label+" uses exact -1 for unset and finite nonnegative values for explicit overrides");
		return value>=0f;
	}
	private static void require(boolean ok, String why){FamilyResourceManifest.require(ok, why);}
	private NcrnaFamilyLoader(){}
}
