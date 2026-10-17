package prok;

import java.util.ArrayList;
import consensus.BaseGraph;
import map.LongHashSet;

/** Constructs the existing caller container without changing verification algorithms.
 * CallGenes uses core-only loading then its established endpoint configuration phase.
 * @author Brian, Raiden */
public final class NcrnaFamilyLoader {

	public static NcrnaFamily load(String familyId){return load(familyId, new Options(), new FamilyResourceResolver());}
	public static NcrnaFamily load(String familyId, Options options, FamilyResourceResolver resolver){
		require(options!=null && resolver!=null, "Explicit load options and resolver required");
		if(!options.enabled){return null;}//A disabled selection never resolves or opens a resource.
		final NcrnaFamilyRegistry.Entry p=NcrnaFamilyRegistry.get(familyId);
		final String lib=resolver.resolve(p.id, "consensus", options.resources).path;
		final byte[][] library=TrnaConsensusBuilder.loadLibrary(lib);
		final String[] names=TrnaConsensusBuilder.lastLoadedNames;
		final boolean ordered=!p.dynamicModels || options.resources.consensus==null;
		require(library!=null && names!=null && library.length>0 && names.length==library.length
			&& (!ordered || library.length==p.models.length),
			"Incompatible "+p.id+" consensus resource: require its complete ordered model roster: "+lib);
		final map.ObjectSet<String> seen=new map.ObjectSet<String>(String.class);
		for(int i=0; i<names.length; i++){
			final String token=RrnaResourceIO.first(names[i]);
			require(!token.isEmpty() && seen.add(token) && (!ordered || p.models[i].equals(token))
				&& library[i]!=null && library[i].length>=(p.dynamicModels ? p.seedK : 1),
				"Incompatible "+p.id+" consensus resource: model order, unique names or sequence length at "+i+": "+lib);
		}
		final BaseGraph[] models;
		final boolean useHbm=p.hasHbm && options.hbm;
		if(useHbm){
			final String hbm=resolver.resolve(p.id, "models", options.resources).path;
			models=TrnaConsensusBuilder.loadModels(hbm);
			CallGenes.requireNcrnaLibraryModelAlignment(p.id, library, models, names, lib, hbm);
		}else{models=null;}
		final String seeds=resolver.resolve(p.id, "seeds", options.resources).path;
		final LongHashSet keys=ProkObject.loadLongKmers(seeds, p.seedK);
		require(keys!=null && keys.size()>=p.minSeedKeys,
			"Incompatible "+p.id+" seed resource: require at least "+p.minSeedKeys+" distinct forward k-mers: "+seeds);
		final int k=options.indexK<0 ? p.indexK : options.indexK;
		final int top=options.topN<0 ? (p.dynamicModels ? library.length : p.topN) : options.topN;
		final int hits=options.minHits<0 ? p.minHits : options.minHits;
		TrnaKmerIndex.validateConfiguration(k, options.floor, options.topFrac, options.qFrac, hits);
		TrnaKmerIndex.validateTopN(top);
		final float weight=Float.isNaN(options.weight) ? p.weight : options.weight;
		final float idPass=Float.isNaN(options.idPass) ? p.idPass : options.idPass;
		final float idBorderline=Float.isNaN(options.idBorderline) ? p.idBorderline : options.idBorderline;
		final int seedMinHits=options.seedMinHits<0 ? p.seedMinHits : options.seedMinHits;
		require(Float.isFinite(idPass) && idPass>=0 && idPass<=1 && Float.isFinite(idBorderline)
			&& idBorderline>=0 && idBorderline<=1 && seedMinHits>0, "Invalid identity or seed gate for "+p.id);
		require(Float.isFinite(weight) && weight>0, "Path weight must be finite and positive for "+p.id);
		final int pad=options.windowPad<0 ? p.windowPad : options.windowPad;
		final float scoreA=Float.isNaN(options.scoreA) ? p.scoreA : options.scoreA;
		final float scoreB=Float.isNaN(options.scoreB) ? weight : options.scoreB;
		final float collapse=Float.isNaN(options.collapse) ? p.collapse : options.collapse;
		final NcrnaFamily family=new NcrnaFamily(p.id, library, models, names, keys, p.seedK, p.minLen, pad,
			k, top, options.adaptive, options.floor, options.topFrac, options.qFrac, hits,
			scoreA, scoreB, idPass, idBorderline, useHbm ? p.hbmPass : 1.01f, collapse,
			null, null, null, null, -1, -1, -1, -1, 0f, p.boundaryOffsets(true), p.boundaryOffsets(false));
		family.outputType=p.outputType;family.seedMinHits=seedMinHits;family.scavengePass2=p.pass2;
		family.rankedModelFallback=p.rankedFallback;family.strictIndexCutoff=p.strictIndexCutoff;family.trimAlignmentExtent=false;
		family.reuseConsensusAlignment=true;family.voteSlack=p.voteSlack;
		family.prok16sPartialMinLen=p.partialMinLen;family.quantumCapMultiplier=p.quantumCap;
		family.preserveContainingWindows=p.preserveContaining;
		if(p.modelCutoffs && options.cutoffs){loadCutoffs(family, options.resources, resolver);}
		if(options.endpoints){family.setRrnaEndpointInference(loadEndpoints(family, options.resources, resolver), null);}
		return family;
	}
	static RrnaEndpointCallerFeatures.Resources loadEndpoints(NcrnaFamily family,
			FamilyResourceResolver.Overrides overrides, FamilyResourceResolver resolver){
		require(family!=null, "Endpoint resources require an initialized registered family");
		final NcrnaFamilyRegistry.Entry p=NcrnaFamilyRegistry.get(family.name);
		final FamilyResourceResolver.Handle tables=resolver.tables(family.name, overrides);
		final RrnaEndpointCallerFeatures.Resources loaded=RrnaEndpointResourceLoader.load(tables, family.modelNames, family.library, p.rawEndpointModels);
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
	public static final class Options {
		public boolean enabled=true, hbm=true, endpoints=true, cutoffs=true, adaptive=false;
		public int indexK=-1, topN=-1, minHits=-1, seedMinHits=-1, windowPad=-1;
		public float floor=0f, topFrac=0f, qFrac=0f, weight=Float.NaN;
		public float idPass=Float.NaN, idBorderline=Float.NaN;
		public float scoreA=Float.NaN, scoreB=Float.NaN, collapse=Float.NaN;
		public final FamilyResourceResolver.Overrides resources=new FamilyResourceResolver.Overrides();
	}
	private static void require(boolean ok, String why){FamilyResourceManifest.require(ok, why);}
	private NcrnaFamilyLoader(){}
}
