package prok;

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
		require(library!=null && names!=null && library.length==p.models.length && names.length==library.length,
			"Consensus must retain the complete ordered "+p.id+" model roster: "+lib);
		for(int i=0; i<names.length; i++){
			require(p.models[i].equals(RrnaResourceIO.first(names[i])) && library[i]!=null && library[i].length>0,
				"Consensus model order differs for "+p.id+" at "+i+": "+names[i]);
		}
		final BaseGraph[] models;
		if(options.hbm){
			final String hbm=resolver.resolve(p.id, "models", options.resources).path;
			models=TrnaConsensusBuilder.loadModels(hbm);
			CallGenes.requireNcrnaLibraryModelAlignment(p.id, library, models, names, lib, hbm);
		}else{models=null;}
		final String seeds=resolver.resolve(p.id, "seeds", options.resources).path;
		final LongHashSet keys=ProkObject.loadLongKmers(seeds, p.seedK);
		require(keys!=null && keys.size()>0, "Empty required "+p.id+" seed resource: "+seeds);
		final int k=options.indexK<0 ? p.indexK : options.indexK;
		final int top=options.topN<0 ? p.topN : options.topN, hits=options.minHits<0 ? p.minHits : options.minHits;
		TrnaKmerIndex.validateConfiguration(k, options.floor, options.topFrac, options.qFrac, hits);
		TrnaKmerIndex.validateTopN(top);
		final float weight=Float.isNaN(options.weight) ? p.weight : options.weight;
		require(Float.isFinite(weight) && weight>0, "Path weight must be finite and positive for "+p.id);
		final NcrnaFamily family=new NcrnaFamily(p.id, library, models, names, keys, p.seedK, p.minLen, p.windowPad,
			k, top, options.adaptive, options.floor, options.topFrac, options.qFrac, hits,
			p.scoreA, weight, p.idPass, p.idBorderline, options.hbm ? p.hbmPass : 1.01f, p.collapse,
			null, null, null, null, -1, -1, -1, -1, 0f, NcrnaFamily.LEGACY_START_OFFSETS, NcrnaFamily.LEGACY_STOP_OFFSETS);
		family.outputType=p.outputType;family.seedMinHits=p.seedMinHits;family.scavengePass2=true;
		family.rankedModelFallback=false;family.strictIndexCutoff=false;family.trimAlignmentExtent=false;
		family.reuseConsensusAlignment=true;family.voteSlack=60;
		if(options.endpoints){family.setRrnaEndpointInference(loadEndpoints(family, options.resources, resolver), null);}
		return family;
	}
	static RrnaEndpointCallerFeatures.Resources loadEndpoints(NcrnaFamily family,
			FamilyResourceResolver.Overrides overrides, FamilyResourceResolver resolver){
		require(family!=null && family.name.equals("prok5S"), "Only the registered 5S pilot may use this endpoint adapter");
		final FamilyResourceResolver.Handle tables=resolver.tables(family.name, overrides);
		final RrnaEndpointCallerFeatures.Resources loaded=RrnaEndpointResourceLoader.load(tables, family.modelNames, family.library);
		final RrnaEndpointCallerFeatures.Resources inference=RrnaEndpointNetworkLoader.load(resolver.resolve(family.name, "endpoints", overrides), loaded);
		require(inference.acceptedVerifierIdentity, "5S nets must preserve the trained accepted-verifier-window identity feature");
		return inference;
	}
	public static final class Options {
		public boolean enabled=true, hbm=true, endpoints=true, adaptive=false;
		public int indexK=-1, topN=-1, minHits=-1;
		public float floor=0f, topFrac=0f, qFrac=0f, weight=Float.NaN;
		public final FamilyResourceResolver.Overrides resources=new FamilyResourceResolver.Overrides();
	}
	private static void require(boolean ok, String why){FamilyResourceManifest.require(ok, why);}
	private NcrnaFamilyLoader(){}
}
