package prok;

import java.io.File;
import java.util.Arrays;
import consensus.BaseGraph;
import fileIO.FileFormat;
import map.LongHashSet;
import ml.CellNet;
import shared.Shared;
import stream.Read;
import stream.ReadInputStream;

/** Resource migration parity against an explicit, unchanged pre-migration resources tree.
 * @author Brian, Raiden */
public final class NcrnaFamilyLoaderTest {
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("baseline="), "Supply baseline=/absolute/pre-migration/resources");
		Shared.setThreads(1);Read.TO_UPPER_CASE=true;
		final String old=args[0].substring(9);
		require(new File(old, "5S_consensus_sequence.fa").isFile(), "Baseline must retain original resources");
		final NcrnaFamily family=NcrnaFamilyLoader.load("5S");
		final byte[][] library=TrnaConsensusBuilder.loadLibrary(old+"/5S_consensus_sequence.fa");
		final String[] names=TrnaConsensusBuilder.lastLoadedNames;
		require(Arrays.deepEquals(library, family.library) && Arrays.equals(names, family.modelNames), "Model order and consensus bytes changed");
		final BaseGraph[] models=TrnaConsensusBuilder.loadModels(old+"/5S.g5.models.hbm");
		require(models.length==family.models.length, "HBM model count changed");
		for(int i=0; i<models.length; i++){
			require(models[i].toText().toString().equals(family.models[i].toText().toString()), "HBM counts changed for model "+i);
		}
		final LongHashSet seeds=ProkObject.loadLongKmers(old+"/5S_17mers.fa.gz", 17);
		compareKeys(seeds, family.kmerSet);compareKeys(seeds, ProkObject.loadLongKmersByType(17, "5S"));
		final RrnaEndpointCallerFeatures.Resources oldTables=RrnaEndpointResourceLoader.load("prok5S", old+"/5S_endpoint_tables", names, library, new String[0]);
		final RrnaEndpointCallerFeatures.Resources endpoints=RrnaEndpointNetworkLoader.load(old+"/5S_endpoints.bbnets.gz", oldTables);
		compareEndpoints(endpoints, family.rrnaEndpointFeatures);
		final Read[] original=ReadInputStream.toReadArray(FileFormat.testInput(old+"/5S_consensus_sequence.fa", FileFormat.FA, null, false, false), -1);
		for(boolean mito:new boolean[]{false, true}){for(boolean plastid:new boolean[]{false, true}){
			final Read[] loaded=ProkObject.loadConsensusSequenceType("5S", mito, plastid);int n=0;
			for(Read r:original){
				if(mito && r.id.startsWith("mito") || plastid && r.id.startsWith("plastid")){continue;}
				require(n<loaded.length && r.id.equals(loaded[n].id) && Arrays.equals(r.bases, loaded[n].bases), "Legacy first/all reference contract changed");n++;
			}require(n==loaded.length, "Legacy organelle-filtered count changed");
		}}
		controls(family, old);
		deferredLegacySlots();
		System.out.println("PASS family loader: models="+family.library.length+" seeds="+seeds.size()+" HBM="+models.length+" networks=8 tables=8; ordered references, exact tables and predictions, defaults and failure controls");
	}
	private static void controls(NcrnaFamily family, String old){
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final RrnaEndpointCallerFeatures.Resources genericTables=RrnaEndpointResourceLoader.load(
			resolver.resolve("5S", "tables", null), family.modelNames, family.library);
		compareTable(genericTables.five[0], family.rrnaEndpointFeatures.five[0]);
		final NcrnaFamilyLoader.Options o=new NcrnaFamilyLoader.Options();
		o.resources.root=old+"/nonexistent-family";o.enabled=false;
		require(NcrnaFamilyLoader.load("5S", o, new FamilyResourceResolver())==null, "Disabled family must return before I/O");
		o.enabled=true;fails(()->NcrnaFamilyLoader.load("5S", o, new FamilyResourceResolver()), "manifest");
		o.resources.root=null;o.resources.models=old+"/absent.hbm";
		fails(()->NcrnaFamilyLoader.load("5S", o, new FamilyResourceResolver()), "models");
		o.hbm=false;o.endpoints=false;o.resources.endpoints=old+"/absent.bbnets";o.resources.tables=old+"/absent-tables";
		final NcrnaFamily disabled=NcrnaFamilyLoader.load("5S", o, new FamilyResourceResolver());
		require(disabled.models==null && disabled.hbmPass==1.01f && disabled.rrnaEndpointFeatures==null,
			"Disabled HBM/endpoint modes must not open their explicit missing resources");
		o.resources.models=null;o.hbm=true;o.resources.seeds=old+"/missing-seeds.fa";
		fails(()->NcrnaFamilyLoader.load("5S", o, new FamilyResourceResolver()), "seeds");
		o.resources.seeds=null;o.endpoints=true;
		fails(()->NcrnaFamilyLoader.load("5S", o, new FamilyResourceResolver()), "tables");
		o.resources.tables=old+"/5S_endpoint_tables";
		fails(()->NcrnaFamilyLoader.load("5S", o, new FamilyResourceResolver()), "endpoints");
		o.resources.endpoints=old+"/5S_endpoints.bbnets.gz";
		o.resources.consensus=old+"/5S_consensus_sequence.fa";o.resources.models=old+"/5S.g5.models.hbm";o.resources.seeds=old+"/5S_17mers.fa.gz";
		o.indexK=7;o.topN=3;o.minHits=2;o.weight=3f;
		final NcrnaFamily explicit=NcrnaFamilyLoader.load("prok5S", o, new FamilyResourceResolver());
		require(explicit.indexK==family.indexK && explicit.indexTopN==family.indexTopN && explicit.fixedMinHits==family.fixedMinHits
			&& explicit.scoreB==family.scoreB && explicit.seedMinHits==family.seedMinHits && explicit.voteTable==null
			&& !explicit.voteWindows && !explicit.voteEnds, "Explicit defaults or absent offset policy changed");
		compareEndpoints(family.rrnaEndpointFeatures, explicit.rrnaEndpointFeatures);
		o.endpoints=false;o.indexK=8;o.topN=2;o.minHits=4;o.weight=2f;
		final NcrnaFamily override=NcrnaFamilyLoader.load("5S", o, new FamilyResourceResolver());
		require(override.indexK==8 && override.indexTopN==2 && override.fixedMinHits==4 && override.scoreB==2f, "Explicit policy overrides lost precedence");
		fails(()->NcrnaFamilyLoader.load("unknown"), "registration");
	}
	private static void deferredLegacySlots(){
		ProkObject.load5Skmers=ProkObject.load5SSequence=false;
		ProkObject.loadLongKmers();ProkObject.loadConsensusSequenceFromFile(false, false);
		require(ProkObject.r5SKmers==null && ProkObject.r5SSequence==null, "Disabled 5S slots must remain unopened");
		final LongHashSet ssu=ProkObject.ssuKmers;final Read[] refs=ProkObject.r16SSequence;
		ProkObject.load5Skmers=ProkObject.load5SSequence=true;
		ProkObject.loadLongKmers();ProkObject.loadConsensusSequenceFromFile(false, false);
		require(ProkObject.r5SKmers!=null && ProkObject.r5SSequence!=null && ProkObject.r5SSequence.length==4,
			"An in-process legacy 5S consumer must complete the deferred resource slots");
		require(ProkObject.ssuKmers==ssu && ProkObject.r16SSequence==refs, "Deferred 5S must not reload unrelated family slots");
	}
	private static void compareKeys(LongHashSet a, LongHashSet b){
		final long[] x=a.toArray(), y=b.toArray();Arrays.sort(x);Arrays.sort(y);
		require(Arrays.equals(x, y), "Caller-visible forward17-mer sets changed");
	}
	private static void compareEndpoints(RrnaEndpointCallerFeatures.Resources a, RrnaEndpointCallerFeatures.Resources b){
		require(a!=null && b!=null && a.radius==b.radius && a.inputs==b.inputs && a.sites==b.sites
			&& a.acceptedVerifierIdentity==b.acceptedVerifierIdentity, "Endpoint identity/geometry changed");
		for(int model=0; model<a.names.length; model++){
			compareTable(a.five[model], b.five[model]);compareTable(a.three[model], b.three[model]);
			compareNet(a.fiveNets[model], b.fiveNets[model]);compareNet(a.threeNets[model], b.threeNets[model]);
		}
	}
	private static void compareTable(RrnaPositionalKmerTable.Table a, RrnaPositionalKmerTable.Table b){
		require(a.k==b.k && a.offset==b.offset && a.radius==b.radius && a.alpha==b.alpha && a.records==b.records
			&& Arrays.equals(a.scores, b.scores) && Arrays.equals(a.valid, b.valid) && Arrays.equals(a.clipped, b.clipped)
			&& Arrays.equals(a.ambiguous, b.ambiguous), "Endpoint table counts, scores or geometry changed: "+a.model+"/"+a.end);
	}
	private static void compareNet(CellNet a, CellNet b){
		require(a.numInputs()==b.numInputs() && a.numOutputs()==b.numOutputs(), "Network dimensions changed");
		final CellNet x=a.copy(false), y=b.copy(false);final float[] input=new float[a.numInputs()];
		for(int row=0; row<8; row++){
			for(int i=0; i<input.length; i++){input[i]=(float)Math.sin((row+1)*(i+1));}
			x.applyInput(input);y.applyInput(input);x.feedForward();y.feedForward();
			require(Arrays.equals(x.getOutput(), y.getOutput()), "Network prediction bits changed");
		}
	}
	private static void fails(Runnable action, String text){
		try{action.run();}catch(IllegalArgumentException e){require(e.getMessage().contains(text), "Failure did not explain missing "+text+": "+e);return;}
		throw new AssertionError("Invalid selected-family resources did not fail loudly: "+text);
	}
	private static void require(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
