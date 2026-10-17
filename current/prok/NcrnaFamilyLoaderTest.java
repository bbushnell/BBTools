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
		require((args.length==1 || args.length==2) && args[0].startsWith("baseline="), "Supply baseline=/absolute/pre-migration/resources [family=16S|23S|euk5S]");
		Shared.setThreads(1);Read.TO_UPPER_CASE=true;
		final String old=args[0].substring(9);
		if(args.length==2){
			if(args[1].equals("family=euk5S")){euk5s(old);return;}
			require(args[1].equals("family=23S") || args[1].equals("family=16S"), "Unknown migrated family fixture");
			longRrna(old, args[1].substring(7));return;
		}
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
	private static void euk5s(String old){
		final String consensus=old+"/euk5S_consensus_sequence.fa.gz", seeds=old+"/euk5S_17mers.fa.gz";
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final NcrnaFamily family=NcrnaFamilyLoader.load("euk5S", new NcrnaFamilyLoader.Options(), resolver);
		final byte[][] library=TrnaConsensusBuilder.loadLibrary(consensus);
		final String[] names=TrnaConsensusBuilder.lastLoadedNames;
		require(Arrays.deepEquals(library, family.library) && Arrays.equals(names, family.modelNames), "euk5S model order and consensus bytes changed");
		compareKeys(ProkObject.loadLongKmers(seeds, 17), family.kmerSet);
		compareKeys(family.kmerSet, ProkObject.loadLongKmersByType(17, "euk5S"));
		final String legacy=new File(old, "euk5S_consensus_sequence.fa").isFile() ? old+"/euk5S_consensus_sequence.fa" : consensus;
		final Read[] original=ReadInputStream.toReadArray(FileFormat.testInput(legacy, FileFormat.FA, null, false, false), -1);
		final Read[] loadedLegacy=ProkObject.loadConsensusSequenceType("euk5S", false, false);
		require(original.length==loadedLegacy.length, "Private legacy euk5S reference count changed");
		for(int i=0; i<original.length; i++){
			require(original[i].id.equals(loadedLegacy[i].id) && Arrays.equals(original[i].bases, loadedLegacy[i].bases),
				"Legacy euk5S reference order or bytes changed");
		}
		final RrnaEndpointCallerFeatures.Resources oldTables=RrnaEndpointResourceLoader.load("euk5S", old+"/euk5S_endpoint", names, library);
		compareEndpoints(RrnaEndpointNetworkLoader.load(old+"/euk5S_endpoints.bbnets.gz", oldTables), family.rrnaEndpointFeatures);
		require(!family.rrnaEndpointFeatures.acceptedVerifierIdentity, "euk5S must retain candidate-span identity");
		for(int i=0; i<names.length; i++){
			require(family.rrnaEndpointFeatures.five[i].counts!=null && family.rrnaEndpointFeatures.three[i].counts!=null,
				"euk5S uses full positional counts rather than the prokaryotic scores-only adapter");
		}
		require(family.library.length==7 && family.models==null && family.indexTopN==7 && family.indexK==7 && family.fixedMinHits==2
			&& family.kLong==17 && family.seedMinHits==2 && !family.seedDistinct && family.minLen==90 && family.windowPad==150
			&& family.scoreA==0f && family.scoreB==3f && family.idPass==.68f && family.idBorderline==.68f && family.hbmPass==1.01f
			&& family.collapseFrac==.75f && family.outputType==ProkObject.r5S && !family.scavengePass2 && !family.rankedModelFallback
			&& !family.strictIndexCutoff && !family.trimAlignmentExtent && family.reuseConsensusAlignment && family.quantumCapMultiplier==1,
			"euk5S migration changed the selected calling policy");
		final NcrnaFamilyLoader.Options o=new NcrnaFamilyLoader.Options();
		o.resources.consensus=consensus;o.resources.seeds=seeds;o.resources.cutoffs=old+"/euk5S_model_cutoffs.tsv";
		o.resources.tables=old+"/euk5S_endpoint";o.resources.endpoints=old+"/euk5S_endpoints.bbnets.gz";
		final NcrnaFamily explicit=NcrnaFamilyLoader.load("euk5S", o, resolver);
		compareEndpoints(family.rrnaEndpointFeatures, explicit.rrnaEndpointFeatures);
		for(int i=0; i<names.length; i++){
			require(family.modelThresholds.pass(i)==explicit.modelThresholds.pass(i)
				&& family.modelThresholds.borderline(i)==explicit.modelThresholds.borderline(i), "Model cutoff binding changed");
		}
		o.resources.cutoffs=old+"/absent-cutoffs.tsv";
		fails(()->NcrnaFamilyLoader.load("euk5S", o, resolver), "cutoffs");
		o.cutoffs=false;o.endpoints=false;o.resources.tables=old+"/absent-tables";o.resources.endpoints=old+"/absent-nets";
		o.scoreA=2f;o.scoreB=4f;o.idPass=.73f;o.idBorderline=.72f;o.collapse=.8f;o.windowPad=123;
		final NcrnaFamily scalar=NcrnaFamilyLoader.load("euk5S", o, resolver);
		require(scalar.modelThresholds==null && scalar.rrnaEndpointFeatures==null && scalar.scoreA==2f && scalar.scoreB==4f
			&& scalar.idPass==.73f && scalar.idBorderline==.72f && scalar.collapseFrac==.8f && scalar.windowPad==123,
			"Scalar controls or disabled cutoff/endpoint resource I/O changed");
		final NcrnaFamilyLoader.Options off=new NcrnaFamilyLoader.Options();off.resources.root=old+"/absent-package";off.enabled=false;
		require(NcrnaFamilyLoader.load("euk5S", off, resolver)==null, "Disabled euk5S must perform no I/O");
		off.enabled=true;fails(()->NcrnaFamilyLoader.load("euk5S", off, resolver), "manifest");
		final NcrnaFamilyLoader.Options core=new NcrnaFamilyLoader.Options();core.endpoints=false;core.cutoffs=false;
		NcrnaFamilyLoader.load("16S", core, resolver);
		final NcrnaFamily again=NcrnaFamilyLoader.load("euk5S", core, resolver);
		require(Arrays.deepEquals(family.library, again.library) && Arrays.equals(family.modelNames, again.modelNames),
			"Switching families on the same resolver must not reuse another manifest");
		compareKeys(family.kmerSet, again.kmerSet);
		System.out.println("PASS euk5S family loader: seven models, exact seeds/cutoffs/full tables/network predictions, candidate-span identity, overrides and same-JVM compatibility");
	}
	private static void longRrna(String old, String id){
		final boolean s16=id.equals("16S");
		final String consensus=old+"/"+id+"_consensus_sequence.fa", seedsPath=old+"/"+id+"_17mers.fa.gz";
		final String[] raw=s16 ? new String[]{"mitochondria"} : new String[0];
		require(new File(consensus).isFile(), "Long-rRNA baseline must retain pre-migration resources");
		final NcrnaFamily family=NcrnaFamilyLoader.load(id);
		final byte[][] library=TrnaConsensusBuilder.loadLibrary(consensus);
		final String[] names=TrnaConsensusBuilder.lastLoadedNames;
		require(Arrays.deepEquals(library, family.library) && Arrays.equals(names, family.modelNames), "Long-rRNA ordered consensus bytes changed");
		final LongHashSet keys=ProkObject.loadLongKmers(seedsPath, 17);
		compareKeys(keys, family.kmerSet);compareKeys(keys, ProkObject.loadLongKmersByType(17, id));
		final RrnaEndpointCallerFeatures.Resources tables=RrnaEndpointResourceLoader.load(id, old+"/"+id+"_endpoint_tables", names, library, raw);
		compareEndpoints(RrnaEndpointNetworkLoader.load(old+"/"+id+"_endpoints.bbnets.gz", tables), family.rrnaEndpointFeatures);
		require(family.models==null && family.outputType==(s16 ? ProkObject.r16S : ProkObject.r23S) && family.kLong==17 && family.minLen==(s16 ? 1400 : 2600)
			&& family.windowPad==(s16 ? 1800 : 3500) && family.indexK==9 && family.indexTopN==(s16 ? 1 : 3) && family.fixedMinHits==(s16 ? 1 : 128)
			&& !family.adaptive && family.adaptFloor==0 && family.adaptTopFrac==0 && family.adaptQFrac==0
			&& family.scoreA==0 && family.scoreB==1 && family.idPass==(s16 ? .66f : .60f) && family.idBorderline==family.idPass && family.hbmPass==1.01f
			&& family.collapseFrac==.85f && family.seedMinHits==2 && !family.seedDistinct && family.scavengePass2==s16
			&& family.rankedModelFallback==!s16 && family.strictIndexCutoff==!s16 && !family.trimAlignmentExtent && family.reuseConsensusAlignment
			&& family.prok16sPartialMinLen==(s16 ? 1200 : 1800) && family.quantumCapMultiplier==(s16 ? 1 : 2) && family.voteSlack==(s16 ? 60 : 150)
			&& family.preserveContainingWindows==s16 && !family.voteWindows && !family.voteEnds && family.voteTable==null,
			"Long-rRNA migration changed the calibrated verification/window policy");
		if(s16){
			require(Arrays.equals(family.boundaryStartOffsets, NcrnaFamily.LEGACY_START_OFFSETS)
				&& Arrays.equals(family.boundaryStopOffsets, NcrnaFamily.LEGACY_STOP_OFFSETS), "16S offset geometry changed");
		}else{require(family.boundaryStartOffsets.length==33 && family.boundaryStartOffsets[0]==-16
			&& family.boundaryStartOffsets[32]==16 && Arrays.equals(family.boundaryStartOffsets, family.boundaryStopOffsets),
			"Historical centered offset geometry must survive core loading");}
		final Read[] original=ReadInputStream.toReadArray(FileFormat.testInput(consensus, FileFormat.FA, null, false, false), -1);
		for(boolean mito:new boolean[]{false, true}){for(boolean plastid:new boolean[]{false, true}){
			final Read[] loaded=ProkObject.loadConsensusSequenceType(id, mito, plastid);int n=0;
			for(Read r:original){
				if(mito && r.id.startsWith("mito") || plastid && r.id.startsWith("plastid")){continue;}
				require(n<loaded.length && r.id.equals(loaded[n].id) && Arrays.equals(r.bases, loaded[n].bases), "Long-rRNA legacy reference order/filter changed");n++;
			}require(n==loaded.length, "Long-rRNA legacy filtered count changed");
		}}
		final NcrnaFamilyLoader.Options o=new NcrnaFamilyLoader.Options();
		o.resources.root=old+"/absent-package";o.enabled=false;
		require(NcrnaFamilyLoader.load(id, o, new FamilyResourceResolver())==null, "Disabled family must perform no resource I/O");
		o.enabled=true;fails(()->NcrnaFamilyLoader.load(id, o, new FamilyResourceResolver()), "manifest");
		o.resources.root=null;o.resources.seeds=old+"/absent-seeds.fa";
		fails(()->NcrnaFamilyLoader.load(id, o, new FamilyResourceResolver()), "seeds");
		o.resources.seeds=seedsPath;o.resources.consensus=consensus;o.resources.tables=old+"/"+id+"_endpoint_tables";
		o.resources.endpoints=old+"/"+id+"_endpoints.bbnets.gz";
		compareEndpoints(family.rrnaEndpointFeatures, NcrnaFamilyLoader.load(id, o, new FamilyResourceResolver()).rrnaEndpointFeatures);
		o.resources.tables=old+"/absent-tables";
		fails(()->NcrnaFamilyLoader.load(id, o, new FamilyResourceResolver()), "tables");
		o.endpoints=false;o.idPass=.65f;o.idBorderline=.64f;o.seedMinHits=3;o.indexK=8;o.topN=2;o.minHits=10;
		final NcrnaFamily override=NcrnaFamilyLoader.load(id, o, new FamilyResourceResolver());
		require(override.idPass==.65f && override.idBorderline==.64f && override.seedMinHits==3 && override.indexK==8
			&& override.indexTopN==2 && override.fixedMinHits==10 && override.rrnaEndpointFeatures==null,
			"Explicit family scalar overrides and endpoint-off must retain precedence");
		if(s16){ProkObject.load16SSequence=false;}else{ProkObject.load23SSequence=false;}
		ProkObject.loadConsensusSequenceFromFile(false, false);
		require((s16 ? ProkObject.r16SSequence : ProkObject.r23SSequence)==null, "Disabled family consensus cache must remain unopened");
		final Read[] refs=ProkObject.r5SSequence;
		if(s16){ProkObject.load16SSequence=true;}else{ProkObject.load23SSequence=true;}
		ProkObject.loadConsensusSequenceFromFile(false, false);
		final Read[] loaded=s16 ? ProkObject.r16SSequence : ProkObject.r23SSequence;
		require(loaded!=null && loaded.length==(s16 ? 5 : 4) && ProkObject.r5SSequence==refs,
			"A later in-process consumer must complete only its deferred family slot");
		System.out.println("PASS "+id+" family loader: exact ordered consensus, forward seeds, endpoint tables/predictions, policy, overrides and deferred cache");
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
			&& a.acceptedVerifierIdentity==b.acceptedVerifierIdentity && Arrays.equals(a.names, b.names)
			&& Arrays.equals(a.raw, b.raw), "Endpoint identity/geometry or raw-model policy changed");
		for(int model=0; model<a.names.length; model++){
			if(a.raw[model]){
				require(a.five[model]==null && a.three[model]==null && a.fiveNets[model]==null && a.threeNets[model]==null
					&& b.five[model]==null && b.three[model]==null && b.fiveNets[model]==null && b.threeNets[model]==null,
					"An explicitly raw model must not acquire endpoint tables or networks");continue;
			}
			compareTable(a.five[model], b.five[model]);compareTable(a.three[model], b.three[model]);
			compareNet(a.fiveNets[model], b.fiveNets[model]);compareNet(a.threeNets[model], b.threeNets[model]);
		}
	}
	private static void compareTable(RrnaPositionalKmerTable.Table a, RrnaPositionalKmerTable.Table b){
		require(a.k==b.k && a.offset==b.offset && a.radius==b.radius && a.alpha==b.alpha && a.records==b.records
			&& Arrays.equals(a.scores, b.scores) && Arrays.equals(a.valid, b.valid) && Arrays.equals(a.clipped, b.clipped)
			&& Arrays.equals(a.ambiguous, b.ambiguous) && Arrays.deepEquals(a.counts, b.counts), "Endpoint table counts, scores or geometry changed: "+a.model+"/"+a.end);
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
