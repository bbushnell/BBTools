package benchmark;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.io.OutputStream;
import java.io.BufferedOutputStream;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

import benchmark.CallGenesLengthScorecard.SavedBands;
import benchmark.CallGenesLengthScorecard.ScoreCounts;
import benchmark.CallGenesStopAudit.ReferenceSet;
import dna.GeneticCode;
import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import gff.GffLine;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import prok.CallGenes;
import prok.GeneCaller;
import prok.GeneModel;
import prok.Orf;
import prok.ProkObject;
import prok.PGMTools;
import prok.StatsContainer;
import prok.StopTailCascadeResources;
import prok.StopTailCascadeSupport;
import prot.FamilyBonusFormula;
import prot.FamilyMatchInformation;
import prot.FamilyTailCascade;
import shared.Shared;
import stream.Read;
import stream.ReadInputStream;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Everyday RNA-on caller measurement using previously accepted locus/input bindings.
 * This is not an adoption/integrity validator. One whole-cohort invocation loads
 * resources once, preselects PGMs serially, calls genomes in parallel, then grades
 * after every caller has joined. No per-genome hash walks or GFF parity passes.
 * @author Sayu
 */
public final class CallGenesLeanReplayBenchmark{

	public static void main(String[] args)throws Exception{
		final PreParser pp=new PreParser(args, CallGenesLeanReplayBenchmark.class, false);
		try{new CallGenesLeanReplayBenchmark(pp.args).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/** Require explicit parallelism and inherited identities rather than choosing a new recipe. */
	private CallGenesLeanReplayBenchmark(String[] args){
		invocation=args.clone();
		final HashSet<String> seen=new HashSet<String>();
		final FamilyBonusFormula.Builder builder=new FamilyBonusFormula.Builder();
		for(String arg:args){
			final String[] kv=Parse.splitOnFirst(arg, '=');
			require(kv.length==2 && !kv[1].isEmpty(), "Require nonempty key=value arguments");
			final String key=kv[0].toLowerCase(Locale.ROOT), value=kv[1];
			require(seen.add(key), "Duplicate lean option: "+key);
			if(key.equals("in")){input=path(value);}
			else if(key.equals("loci")){loci=path(value);}
			else if(key.equals("locisha80")){lociPin=pinLabel(value);}
			else if(key.equals("pgm")){pgm=path(value);}
			else if(key.equals("pgmpolicy")){policy=CallGenesReplayRnaPolicy.parseModelPolicy(value);}
			else if(key.equals("selectioncache")){selectionCachePath=path(value);}
			else if(key.equals("selectioncachesha80")){selectionCachePin=measuredPin(value);}
			else if(key.equals("modelinventory")){modelInventoryPath=path(value);}
			else if(key.equals("modelinventorysha80")){modelInventoryPin=measuredPin(value);}
			else if(key.equals("modelrelocationmap")){modelRelocationMap=path(value);}
			else if(key.equals("modelrelocationroot")){modelRelocationRoot=path(value);}
			else if(key.equals("retentionbgzip")){retentionBgzip=path(value);}
			else if(key.equals("retentionbgzipsha80")){/* Legacy metadata only; not verified. */}
			else if(key.equals("retentionpigz")){retentionPigz=path(value);}
			else if(key.equals("retentionpigzsha80")){/* Legacy metadata only; not verified. */}
			else if(key.equals("alignmentwork")){
				require(value.equals("true") || value.equals("false"), "alignmentwork requires true or false");
				recordAlignmentWork=Boolean.parseBoolean(value);
			}
			else if(key.equals("retaingffs")){
				require(value.equals("true") || value.equals("false"), "retaingffs requires true or false");
				retainGffs=Boolean.parseBoolean(value);
			}
			else if(key.equals("retainrna")){
				require(value.equals("true") || value.equals("false"), "retainrna requires true or false");
				retainRna=Boolean.parseBoolean(value);
			}
			else if(key.equals("retentionworkers")){retentionWorkers=Parse.parseIntKMG(value);}
			else if(key.equals("verifyselectioncache")){
				require(value.equals("true") || value.equals("false"), "verifyselectioncache requires true or false");
				verifySelectionCache=Boolean.parseBoolean(value);
			}
			else if(key.equals("tailreference")){resource=path(value);}
			else if(key.equals("tailreferencesha80")){resourcePin=pinLabel(value);}
			else if(key.equals("tailinformation")){information=path(value);}
			else if(key.equals("fmnroot")){fmn=path(value);}
			else if(key.equals("out")){output=path(value);}
			else if(key.equals("genomes")){expected=Parse.parseIntKMG(value);}
			else if(key.equals("workers")){workers=Parse.parseIntKMG(value);}
			else if(key.equals("searchthreads")){searchThreads=Parse.parseIntKMG(value);}
			else if(key.equals("passes")){passes=Integer.parseInt(value);}
			else if(key.equals("evidencecache")){
				require(value.equals("true") || value.equals("false"), "evidencecache requires true or false");
				reusePassEvidence=Boolean.parseBoolean(value);
			}
			else if(key.equals("cpus")){cpus=Parse.parseIntKMG(value);}
			else if(key.equals("top")){
				top=Parse.parseIntKMG(value);
				require(top>0, "Explicit top must be positive");
			}
			else if(key.equals("minhits")){
				minHits=Parse.parseIntKMG(value);
				require(minHits>=0, "minhits must be nonnegative");
			}
			else if(key.equals("identitylowerreduction")){identityReduction=Float.parseFloat(value);}
			else if(key.equals("tailrawscorei0override")){rawScoreI0Override=Parse.parseBoolean(value);}
			else if(builder.parse(key, value)){}
			else{throw new IllegalArgumentException("Unknown lean option: "+key);}
		}
		formula=builder.build();
		require(input!=null && loci!=null && lociPin!=null && pgm!=null && resource!=null && resourcePin!=null
			&& fmn!=null && output!=null && policy!=null && formula!=null,
			"Explicit inputs, locus identity, PGM policy, RNA root and formula required");
		require(expected>0 && workers>0 && workers<=expected && searchThreads>0 && cpus>0
			&& (long)workers+searchThreads<=cpus, "Genome and shared search workers must fit the explicit CPU budget");
		require(passes==1 || passes==2, "Lean passes must be 1 or 2; each policy trains from its own first-pass calls");
		require(passes!=2 || !recordAlignmentWork,
			"Two-pass work requires separate per-pass counters; existing alignmentwork counters are one-pass only");
		require(Float.isFinite(identityReduction) && identityReduction>=0 && identityReduction<=100,
			"Invalid identity reduction");
		require(formula.requiresInformation()==(information!=null),
			"Information input must match the explicit formula requirement");
		require((selectionCachePath==null)==(selectionCachePin==null)
			&& (modelInventoryPath==null)==(modelInventoryPin==null),
			"Selection/model manifest paths require paired sha80 pins");
		require(selectionCachePath==null || modelInventoryPath!=null,
			"Selection cache requires its complete model inventory");
		require(!verifySelectionCache || selectionCachePath!=null, "Fresh verification requires a selection cache");
		require(modelInventoryPath==null || policy==CallGenesReplayRnaPolicy.ModelPolicy.PER_PHYLUM,
			"Selection cache/export requires per-phylum policy");
		require((modelRelocationMap==null)==(modelRelocationRoot==null),
			"Model relocation requires both map and root");
		require(modelRelocationMap==null || selectionCachePath!=null && modelInventoryPath!=null && !verifySelectionCache,
			"Model relocation requires the original inventory/cache and cached taxonomy without fresh verification");
		require(!(retainGffs || recordAlignmentWork) || retentionBgzip!=null && retentionPigz!=null,
				"Retained output requires explicit canonical compression tool paths");
		require(!retainRna || retainGffs, "RNA retention requires retaingffs=true and never enters CDS grading");
		require(retentionWorkers>0 && retentionWorkers<=cpus && (retainGffs || recordAlignmentWork || retentionWorkers==1),
			"Retention workers must fit the post-caller CPU budget and require retained output");
		CallGenes.validateFormulaBonusPolicy(true, top!=0 || minHits>=0 || identityReduction!=0, 0);
	}

	/** One owner scopes all initialization, worker joins and cleanup, including failures. */
	private void run()throws Exception{
		final long start=System.nanoTime(), cpu=CallGenesBenchmarkSupport.processCpuNanos();
		require(!Files.exists(output) && output.getParent()!=null && Files.isDirectory(output.getParent()),
			"Require fresh output under an existing parent");
		final ArrayList<Genome> genomes=readInputs();
		final StopTailCascadeResources resources=StopTailCascadeResources.readStructural(resource.toString(), resourcePin);
		require(resources.minQueryAa==100 && resources.maxQueryAa==100 && resources.firstStart && resources.allFamilySets
			&& resources.hasShortBins(), "Lean comparison currently binds the accepted cap100/SHORT40..99 scanner only");
		checkScannerRecipes(genomes, resources);
		resources.protectOutputs(new String[]{output.resolve("scorecard.tsv").toString()});
		Files.createDirectory(output);
		if(recordAlignmentWork){
			Files.createDirectory(output.resolve("alignment_work"));
			Files.write(output.resolve("alignment_work_notes.txt"), (
				"One primitive row per stop/view with alignment_calls>0 or alignment_accepts>0. Zero-event eligible/shortlist counts remain complete in aggregate table. No reference classification in caller. query_length_aa is the actual per-view peptide length, before alignment.\n"+
				"stop_left_zero_based is the scanner stop codon genomic left edge; strand0 plus,1 minus. TrueStop1based=left+(plus?3:1).\n"+
				"SHORT and LONG are actual family-set views, not peptide length bins. Queries counted once before four replay arms.\n"+
				"shortlist means ranked top cap before scalar gates. alignment_calls counts real aligner invocations, including failed/short-span DP, excluding short-input skips, disabled calibration, cached-provider calls and shortcuts.\n"+
				"alignment_accepts counts passes at the real identity/BLOSUM accept sites. Ordinary search stops at its first accepted family, so at most one per query.\n"+
				"No additional alignment or traceback is performed. Per-view CPU/wall stage timing is unavailable; existing whole-run phases and process CPU retained. Final tally files use BGZF level9 with one compressor thread; raw temporary TSV files removed only after successful compression.\n").getBytes(StandardCharsets.US_ASCII));
		}
		Files.copy(input, output.resolve("accepted_inputs.tsv"));
		Files.write(output.resolve("benchmark.args"),
			(String.join("\n", invocation)+"\n").getBytes(StandardCharsets.UTF_8));
		final String[] args=CallGenesReplayRnaPolicy.controlArgs(genomes.get(0).fasta.toString(),
			output.resolve("unused_control.gff").toString(), pgm.toString(),
			CallGenesReplayRnaPolicy.Mode.NORMAL, fmn, policy);
		Files.write(output.resolve("startup.args"), (String.join("\n", args)+"\n").getBytes(StandardCharsets.UTF_8));
		final long initialization=System.nanoTime();
		final CallGenes owner=new CallGenes(args);
		owner.withReplayResources(()->{
			try{
				phase("caller_initialization", System.nanoTime()-initialization);
				CallGenesReplayRnaPolicy.writeRuntime(output.resolve("rna_runtime.tsv"),
					CallGenesReplayRnaPolicy.Mode.NORMAL);
				final String[] defaults=owner.configuredModelPaths();
				require(defaults.length==1 && Paths.get(defaults[0]).toRealPath().equals(pgm.toRealPath()),
					"PGM must bind the actual selected fallback; per-phylum mode must not override ordinary defaults");
				execute(owner, genomes, resources);
			}catch(Exception problem){throw new IllegalStateException("Lean benchmark failed; no PASS", problem);}
		});
		require(!Thread.currentThread().isInterrupted(), "Interrupted execution cannot publish PASS");
		phase("measured_run", System.nanoTime()-start);
		final String timing="metric\tvalue\nprocess_cpu_seconds\t"+(CallGenesBenchmarkSupport.processCpuNanos()-cpu)*1e-9+
			"\njava_peak_threads\t"+ManagementFactory.getThreadMXBean().getPeakThreadCount()+"\nheap_limit_bytes\t"+Runtime.getRuntime().maxMemory()+"\n";
		Files.write(output.resolve("runtime.tsv"), timing.getBytes(StandardCharsets.UTF_8));
		Files.write(output.resolve("phases.tsv"), phases.toBytes());
		Files.write(output.resolve("recipe.txt"), ("mode=lean_measurement_only\nin="+input+"\nloci="+loci+
			"\nlocisha80_declared="+lociPin+"\ntailreference="+resource+"\ntailreferencesha80_declared="+resourcePin+
			"\npgm="+pgm+"\npgmpolicy="+policy+"\nrnacalling=normal\nfmnroot="+fmn+"\nworkers="+workers+
			"\nmodel_selection_mode="+(selectionCachePath==null ? "fresh" : verifySelectionCache ? "fresh_verify" : "cached")+
			"\nselectioncache="+selectionCachePath+"\nselectioncachesha80="+selectionCachePin+
			"\nmodelinventory="+modelInventoryPath+"\nmodelinventorysha80="+modelInventoryPin+
			(modelRelocationMap==null ? "" : "\nmodelrelocationmap="+modelRelocationMap+"\nmodelrelocationroot="+modelRelocationRoot+
				"\nmodel_relocation_identity=original_manifests_inherited;external_byte_comparison_required;no_new_digests")+
			"\nselection_content_identity=inherited_from_packet_guards;no_per_genome_hash_walk"+
			"\nsearchthreads="+searchThreads+"\ncpus="+cpus+"\ngenomes="+expected+
			"\npasses="+passes+"\nrefinement="+(passes==1 ? "none" : "own_policy_calls_GeneModel_process_PGMTools_mix_0.50_0.10_0.10")+
			"\nfamily_evidence_scope="+(passes==1 ? "single_pass" : "both_passes")+
			"\nevidencecache="+(passes==2 && reusePassEvidence)+"\nmodel_selection_scope=frozen_initial_model_before_refinement"+
			"\naccepted_input_binding=accepted_inputs.tsv;identities_declared_not_recomputed"+
			"\nbenchmark_arguments=benchmark.args\neffective_top="+(top==0 ? "manifest" : Integer.toString(top))+
			"\neffective_minhits="+(minHits<0 ? "manifest" : Integer.toString(minHits))+
			"\nidentity_lower_reduction_pp="+identityReduction+"\ntailinformation="+(information==null ? "none" : information)+
			"\ninput_hashes_recomputed=false\nresource_loader=structural\nscanner_table=11\nquery_cap=100\n"+
			"zero_control=whole_cohort_aggregate_16_counts_all_bands;bonus0;protection_disabled\n"+
			"zero_control_limit=cancelling_per_genome_changes_can_escape;guarded_mode_required_for_adoption\n"+
			"grading=after_all_caller_workers_join;frozen_reference_bands\n"+
			"retaingffs="+retainGffs+"\nretainrna="+retainRna+"\nretentionworkers="+retentionWorkers+"\n"+
			"timing=whole_process_CPU_once;worker_CPU_excludes_shared_search_and_IO_threads\n").getBytes(StandardCharsets.UTF_8));
		Files.write(output.resolve("PASS"),
			("CALLGENES_LEAN_MEASUREMENT_PASS genomes="+expected+" adoption=false\n").getBytes(StandardCharsets.US_ASCII));
	}

	/** Publish every model before workers, then use one shared classifier and saved band map. */
	private void execute(CallGenes owner, ArrayList<Genome> genomes, StopTailCascadeResources resources)throws Exception{
		long start=System.nanoTime();
		final CallGenesReplaySelectionCache.Inventory inventory=modelInventoryPath==null ? null :
			modelRelocationMap==null ? CallGenesReplaySelectionCache.readInventory(modelInventoryPath, modelInventoryPin) :
			CallGenesReplaySelectionCache.readRelocatedInventory(modelInventoryPath, modelInventoryPin,
				modelRelocationMap, modelRelocationRoot);
		if(inventory!=null){inventory.requireFallback(pgm.toString());}
		final CallGenesReplaySelectionCache cache=selectionCachePath==null ? null :
			modelRelocationMap==null ? CallGenesReplaySelectionCache.read(selectionCachePath, selectionCachePin, inventory) :
			CallGenesReplaySelectionCache.readRelocatedCache(selectionCachePath, selectionCachePin, inventory);
		if(cache!=null){
			require(cache.size()==genomes.size(), "Selection cache must cover the complete ordered cohort");
			for(int i=0;i<genomes.size();i++){
				final Genome genome=genomes.get(i);
				cache.requireGenome(i, genome.id, genome.fasta.toString(), genome.fastaPin);
			}
		}
		final ByteBuilder cacheOutput=inventory==null ? null : CallGenesReplaySelectionCache.begin(inventory);
		final ByteBuilder resolution=modelRelocationMap==null ? null : CallGenesReplaySelectionCache.beginResolution();
		final CallGenes.FileModelSelection[] selections=new CallGenes.FileModelSelection[genomes.size()];
		final ByteBuilder selected=new ByteBuilder("genome\tpolicy\tdomain\tphylum\tsource\tfallback\tcache_hit\tmodel_path\n");
		for(int i=0;i<genomes.size();i++){
			final Genome genome=genomes.get(i);
			final CallGenes.FileModelSelection selection=cache==null ? owner.selectModelForReplay(genome.fasta.toString()) :
				cache.select(owner, i, verifySelectionCache);
			selections[i]=selection;
			if(resolution!=null){
				final String[] paths=selection.modelPaths();
				require(paths.length==1, "Relocation evidence requires one actual selected model");
				cache.appendResolved(cacheOutput, resolution, i, paths[0], selection.domain, selection.phylum,
					selection.classificationSource, selection.fallback, selection.cacheHit);
			}else if(cacheOutput!=null){
				CallGenesReplaySelectionCache.append(cacheOutput, inventory, genome.id,
					genome.fasta.toString(), genome.fastaPin, selection);
			}
			for(String name:selection.modelPaths()){
				selected.append(genomes.get(i).id).tab().append(policy.toString()).tab().append(text(selection.domain)).tab()
					.append(text(selection.phylum)).tab().append(text(selection.classificationSource)).tab().append(selection.fallback)
					.tab().append(selection.cacheHit ? 1 : 0).tab().append(name).nl();
			}
		}
		for(CallGenes.FileModelSelection selection:selections){
			CallGenes.validateReplayModelGeometry(selections[0].model, selection.model);
		}
		if(resolution!=null){
			cache.verifyInverse(cacheOutput);
			Files.write(output.resolve("selection_cache.inverse.tsv"), cacheOutput.toBytes());
			Files.write(output.resolve("model_resolution.tsv"), resolution.toBytes());
			Files.copy(selectionCachePath, output.resolve("selection_cache.tsv"));
		}else if(cacheOutput!=null){Files.write(output.resolve("selection_cache.tsv"), cacheOutput.toBytes());}
		Files.write(output.resolve("pgm_selections.tsv"), selected.toBytes());
		phase("model_selection", System.nanoTime()-start);
		start=System.nanoTime();
		final HashSet<String> ids=new HashSet<String>();
		for(Genome genome:genomes){ids.add(genome.id);}
		final SavedBands bands=CallGenesLengthScorecard.readSavedBands(loci, ids);
		phase("saved_bands_load", System.nanoTime()-start);
		final ThreadMXBean clock=ManagementFactory.getThreadMXBean();
		require(clock.isCurrentThreadCpuTimeSupported(), "Worker CPU clock unavailable");
		if(!clock.isThreadCpuTimeEnabled()){clock.setThreadCpuTimeEnabled(true);}
		start=System.nanoTime();
		final FamilyMatchInformation info=information==null ? null : FamilyMatchInformation.read(information.toString());
		try(StopTailCascadeSupport support=resources.loadSupportWithOverrides(searchThreads, top, minHits,
			identityReduction, info)){
			final FamilyBonusFormula boundFormula=support.bindFormulaEvidence(formula, rawScoreI0Override);
			Files.write(output.resolve("formula.tsv"), boundFormula.recipe().toBytes());
			Shared.AMINO_IN=false;
			require(ProkObject.callCDS, "Lean measurement requires CDS calling alongside RNA");
			phase("classifier_load", System.nanoTime()-start);
			start=System.nanoTime();
			final ArrayList<Result> results=CallGenesLeanWorkers.run(genomes.size(), workers,
				ordinal->callGenome(genomes.get(ordinal), selections[ordinal].model, support, boundFormula, clock));
			phase("parallel_calling", System.nanoTime()-start);
			start=System.nanoTime();
			grade(genomes, results, bands);
			phase("grading_and_tables", System.nanoTime()-start);
		}
		if(inventory!=null){inventory.verify();}
	}

	/** Every reader, sequence copy, replay and GeneCaller belongs to this task alone. */
	private Result callGenome(Genome genome, GeneModel selectedModel, StopTailCascadeSupport support,
			FamilyBonusFormula boundFormula, ThreadMXBean clock){
		if(passes==2){return callGenomeTwoPass(genome, selectedModel, support, boundFormula, clock);}
		final long start=System.nanoTime(), cpu=clock.getCurrentThreadCpuTime();
		require(cpu>=0, "Worker CPU clock returned no measurement");
		final Result result=new Result(genome.id);
		final GeneticCode code=GeneticCode.forTable(11);
		final GeneModel model=selectedModel;
		final GeneCaller[] callers=new GeneCaller[ARMS.length];
		for(int a=0;a<ARMS.length;a++){
			callers[a]=CallGenes.makeGeneCaller(model);
			if(a==3){callers[a].setTailCascadeFormulaSupport(support, boundFormula, true, resourcePin, true,
				rawScoreI0Override);}
			else if(a>0){callers[a].setTailCascadeSupportStructural(support, a==1 ? 0 : 400, a>1,
				resourcePin, false);}
		}
		final Streamer reader=StreamerFactory.makeStreamer(FileFormat.testInput(genome.fasta.toString(),
			FileFormat.FASTA, null, false, true), null, true, -1, false, true, 1);
		final HashSet<String> contigs=new HashSet<String>();
		Throwable failure=null;
		final ByteBuilder workBuffer=recordAlignmentWork ? new ByteBuilder(65536) : null;
		try(OutputStream work=recordAlignmentWork ? new BufferedOutputStream(Files.newOutputStream(
				output.resolve("alignment_work").resolve(genome.id+".tsv.tmp")), 65536) : null){
			if(work!=null){work.write(WORK_STOP_HEADER.getBytes(StandardCharsets.US_ASCII));}
			reader.start();
			for(ListNum<Read> batch=reader.nextList();batch!=null;batch=reader.nextList()){
				for(Read read:batch){
						require(read!=null && read.id!=null && read.bases!=null && read.mate==null,
							"Require named unpaired nucleotide contigs");
					if(!read.validated()){read.validate(true);}
						final String name=contigName(read.id);
						require(contigs.add(name), "Duplicate GFF contig key in "+genome.id);
						final StopTailCascadeSupport.Replay replay=support.prepareReplay(read, code);
						result.queries+=replay.queryCount();
						result.accepted+=replay.acceptedCount();
						if(recordAlignmentWork){countSearchWork(result, replay.evidence(), name, work, workBuffer);}
					for(int a=0;a<ARMS.length;a++){
							final ArrayList<Orf> calls=replay.call(callers[a], model);
							captureCalls(result, a, name, calls);
					}
				}
			}
		}catch(IOException problem){
			failure=problem;
			throw new UncheckedIOException("Alignment-work output failed: "+genome.id, problem);
		}
		catch(RuntimeException | Error problem){
			failure=problem;
			throw problem;
		}
		finally{
			if(ReadWrite.closeStream(reader)){
				final IllegalStateException problem=new IllegalStateException("Genome reader failed: "+genome.id);
				if(failure!=null){failure.addSuppressed(problem);}else{throw problem;}
			}
		}
		for(int a=1;a<ARMS.length;a++){
			final GeneCaller.TailCascadeStats stats=callers[a].tailCascadeStats;
			result.searched[a]=stats.searched;
			result.reused[a]=stats.reused;
			result.accepts[a]=stats.accepted;
		}
		result.wall=System.nanoTime()-start;
		result.cpu=clock.getCurrentThreadCpuTime()-cpu;
		return result;
	}

	/** Whole-genome training follows ordinary runOnePass/makeMultipassModel weights,
	 * but each policy supplies its own calls. Sequence-derived evidence may be shared;
	 * start choice, native pruning, formula bonus and DP are always recomputed. */
	private Result callGenomeTwoPass(Genome genome, GeneModel selectedModel, StopTailCascadeSupport support,
			FamilyBonusFormula boundFormula, ThreadMXBean clock){
		final long start=System.nanoTime(), cpu=clock.getCurrentThreadCpuTime();
		require(cpu>=0, "Worker CPU clock returned no measurement");
		final Result result=new Result(genome.id);
		final GeneticCode code=GeneticCode.forTable(11);
		final ArrayList<Read> reads=ReadInputStream.toReads(genome.fasta.toString(), FileFormat.FASTA, -1);
		final GeneCaller[] first=new GeneCaller[ARMS.length];
		final ArrayList<ArrayList<Orf>> trainingCalls=new ArrayList<ArrayList<Orf>>(ARMS.length);
		for(int a=0;a<ARMS.length;a++){
			first[a]=makePolicyCaller(selectedModel, support, boundFormula, a);
			trainingCalls.add(new ArrayList<Orf>());
		}
		final StopTailCascadeSupport.Replay[] cached=new StopTailCascadeSupport.Replay[reads.size()];
		final HashSet<String> names=new HashSet<String>();
		for(int i=0;i<reads.size();i++){
			final Read read=reads.get(i);
			require(read!=null && read.id!=null && read.bases!=null && read.mate==null,
				"Require named unpaired nucleotide contigs");
			if(!read.validated()){read.validate(true);}
			require(names.add(contigName(read.id)), "Duplicate GFF contig key in "+genome.id);
			final StopTailCascadeSupport.Replay replay=support.prepareReplay(read, code);
			result.firstQueries+=replay.queryCount();
			result.firstAccepted+=replay.acceptedCount();
			if(reusePassEvidence){cached[i]=replay;}
			for(int a=0;a<ARMS.length;a++){trainingCalls.get(a).addAll(replay.call(first[a], selectedModel));}
		}
		final GeneModel[] models=new GeneModel[ARMS.length];
		final GeneCaller[] callers=new GeneCaller[ARMS.length];
		for(int a=0;a<ARMS.length;a++){
			final ArrayList<GffLine>[] lines=Orf.toGffLinesByType(trainingCalls.get(a), "CDS,rRNA,tRNA");
			result.firstCds[a]=lines[0].size();
			// makeGeneCaller(model) preserves the benchmark's legacy caller codon
			// policy. Explicit table11 above belongs only to family-query scanning.
			final GeneModel trained=new GeneModel(true);
			require(!trained.process(reads, lines[0], lines[1], lines[2]),
				"First-pass model training failed for "+genome.id+" policy="+ARMS[a]);
			for(StatsContainer container:trained.allContainers){container.calculate();}
			// Same whole-genome blend as CallGenes.makeMultipassModel(...,passes=2).
			models[a]=PGMTools.mix(0.50, 0.10, 0.10, true, selectedModel, trained);
			CallGenes.validateReplayModelGeometry(selectedModel, models[a]);
			trainingCalls.get(a).clear();
			callers[a]=makePolicyCaller(models[a], support, boundFormula, a);
		}
		for(int i=0;i<reads.size();i++){
			final Read read=reads.get(i);
			StopTailCascadeSupport.Replay replay=cached[i];
			if(replay==null){
				// A missing contig entry computes every eligible stop, including any stop
				// newly reached by the caller. The fixed-DNA scanner has no model input.
				replay=support.prepareReplay(read, code);
				result.recomputedQueries+=replay.queryCount();
			}else{result.cachedQueries+=replay.queryCount();}
			result.queries+=replay.queryCount();
			result.accepted+=replay.acceptedCount();
			for(int a=0;a<ARMS.length;a++){
				captureCalls(result, a, contigName(read.id), replay.call(callers[a], models[a]));
			}
			cached[i]=null;
		}
		require(result.firstQueries==result.queries && result.firstAccepted==result.accepted,
			"Identical genome/code/resources must preserve the sequence-derived stop evidence across refinement");
		for(int a=1;a<ARMS.length;a++){
			final GeneCaller.TailCascadeStats stats=callers[a].tailCascadeStats;
			result.searched[a]=stats.searched;
			result.reused[a]=stats.reused;
			result.accepts[a]=stats.accepted;
		}
		result.wall=System.nanoTime()-start;
		result.cpu=clock.getCurrentThreadCpuTime()-cpu;
		return result;
	}

	private GeneCaller makePolicyCaller(GeneModel model, StopTailCascadeSupport support,
		FamilyBonusFormula formula, int arm){
		final GeneCaller caller=CallGenes.makeGeneCaller(model);
		if(arm==3){caller.setTailCascadeFormulaSupport(support, formula, true, resourcePin, true,
			rawScoreI0Override);}
		else if(arm>0){caller.setTailCascadeSupportStructural(support, arm==1 ? 0 : 400, arm>1,
			resourcePin, false);}
		return caller;
	}

	/** Primitive work events only, after search workers join. No reference labels enter
	 * the caller. Only stop/views with a real call or accept need a trace row; totals cover all eligible queries.
	 * Replaying the four caller arms does not write or count another search.
	 */
	private static void countSearchWork(Result result, StopTailCascadeSupport.ContigEvidence evidence, String contig,
			OutputStream out, ByteBuilder buffer)throws IOException{
		require(evidence.searchesAllFamilySets(), "Per-view counters require explicit all-family-set evidence");
		for(int i=0;i<evidence.size();i++){
			final long key=evidence.key(i);
			final int strand=(int)(key&1), left=(int)(key>>>1);
			for(int view=0;view<2;view++){
				final FamilyTailCascade.Assignment a=evidence.familySetAssignment(i, view==0);
				if(a==null){continue;}
				final long[] n=result.work[view];
				n[0]++;
				if(a.shortlistCandidates>0){n[1]++;}
				n[2]+=a.alignmentCalls;
				n[3]+=a.alignmentAccepts;
				n[4]+=a.alignmentInputSkips;
				n[5]+=a.alignmentUnavailable;
				n[6]+=a.alignmentProviderCalls;
				final boolean shortcut=a.assigned() && a.alignmentCalls==0 && a.alignmentProviderCalls==0;
				if(shortcut){n[7]++;}
				if(a.alignmentCalls>0){n[8]++;}
				if(a.alignmentAccepts>0){n[9]++;}
				n[10]+=a.dpAttempts;
				require(a.alignmentAccepts<=a.alignmentCalls && a.alignmentAccepts<=1 &&
					a.dpAttempts==a.alignmentCalls+a.alignmentInputSkips+a.alignmentUnavailable+a.alignmentProviderCalls,
					"Actual alignment calls/skips/providers must partition alignment-stage visits");
				if(a.alignmentCalls==0 && a.alignmentAccepts==0){continue;}
				buffer.append(result.id).tab().append(contig).tab().append(left).tab().append(strand).tab()
					.append(view==0 ? "SHORT" : "LONG").tab().append(a.shortlistCandidates>0 ? 1 : 0).tab()
					.append(a.alignmentCalls).tab().append(a.alignmentAccepts).tab().append(a.alignmentInputSkips).tab()
					.append(a.alignmentUnavailable).tab().append(a.alignmentProviderCalls).tab().append(shortcut ? 1 : 0).tab()
					.append(a.dpAttempts).tab().append(evidence.familySetQueryLength(i, view==0)).nl();
				if(buffer.length()>=65536){
					out.write(buffer.array, 0, buffer.length());
					buffer.clear();
				}
			}
		}
		if(buffer.length()>0){
			out.write(buffer.array, 0, buffer.length());
			buffer.clear();
		}
	}
	private static void appendWork(ByteBuilder out, String genome, long[][] work){
		for(int v=0;v<2;v++){
			out.append(genome).tab().append(v==0 ? "SHORT" : "LONG");
			for(long n:work[v]){out.tab().append(n);}
			out.nl();
		}
	}

	/** Copy native RNA bytes before compact CDS projection discards non-CDS fields. */
	private void captureCalls(Result result, int arm, String contig, ArrayList<Orf> calls){
		require(result!=null && arm>=0 && arm<ARMS.length, "Capture must bind one completed caller arm");
		if(retainRna && (arm==0 || arm==3)){
			if(result.rna==null){
				result.rna=new CallGenesRetainedExport.Rna[]{
					new CallGenesRetainedExport.Rna(), new CallGenesRetainedExport.Rna()};
			}
			result.rna[arm==0 ? 0 : 1].append(contig, calls);
		}
		result.calls[arm].append(contig, calls);
	}

	/** All grading and the single conservation/control check occur after the worker barrier. */
	private void grade(ArrayList<Genome> genomes, ArrayList<Result> results, SavedBands bands){
		require(results.size()==genomes.size(), "Completed cohort differs from ordered input roster");
		final ArrayList<ArrayList<ScoreCounts>> totals=new ArrayList<ArrayList<ScoreCounts>>();
		for(String arm:ARMS){totals.add(new ArrayList<ScoreCounts>());}
		final long[] cds=new long[ARMS.length], searched=new long[ARMS.length],
			reused=new long[ARMS.length], accepted=new long[ARMS.length];
		long queries=0, accepts=0;
		final ByteBuilder counts=new ByteBuilder(CallGenesLengthScorecard.scoreHeader());
		final ByteBuilder timing=new ByteBuilder("genome\tworker_wall_seconds\tworker_cpu_seconds\n");
		final ByteBuilder features=new ByteBuilder("genome\tarm\ttotal_features\tCDS\tnon_CDS\n");
		final ByteBuilder workRows=new ByteBuilder(WORK_HEADER);
		final ByteBuilder refinement=new ByteBuilder("genome\tpasses\tevidencecache\tfirst_pass_queries\tfirst_pass_accepts\tfinal_pass_queries\tcached_queries\trecomputed_queries\tnative_first_CDS\tformula_first_CDS\n");
		final long[][] workTotal=new long[2][WORK_FIELDS];
		final Path retained=retainGffs ? output.resolve("retained_gffs") : null;
		final ByteBuilder inventory=retainGffs ? new ByteBuilder("genome\tnative_gff\tformula_gff\n") : null;
		final ByteBuilder exportCosts=retainGffs ? new ByteBuilder("genome\tpolicy\traw_payload_bytes\tstored_compressed_bytes\texport_wall_seconds\tcompression\tcommand\n") : null;
		final long[] exportTotals=retainGffs ? new long[5] : null;
		if(retainGffs){
			validateRetentionCompression();
			try{Files.createDirectory(retained);}
			catch(IOException failure){
				throw new UncheckedIOException("Retained GFF directory must be fresh: "+retained, failure);
			}
		}
		for(int i=0;i<genomes.size();i++){
			final Genome genome=genomes.get(i);
			final Result result=results.get(i);
			require(result!=null && genome.id.equals(result.id), "Missing, repeated or reordered genome result");
			final ReferenceSet refs=CallGenesStopAudit.loadReference(genome.reference.toString());
			final CallGenesLengthScorecard.ReferenceBands reference=bands.bind(genome.id, refs);
			queries+=result.queries;
			accepts+=result.accepted;
			refinement.append(genome.id).tab().append(passes).tab().append(passes==2 && reusePassEvidence)
				.tab().append(result.firstQueries).tab().append(result.firstAccepted).tab()
				.append(result.queries).tab().append(result.cachedQueries).tab()
				.append(result.recomputedQueries).tab().append(result.firstCds[0]).tab()
				.append(result.firstCds[3]).nl();
			for(int a=0;a<ARMS.length;a++){
				final CallGenesLeanCalls calls=result.calls[a];
				final ArrayList<GffLine> lines=calls.lines();
				final ScoreCounts score=reference.grade(lines);
				//Persist the already graded compact records only after caller workers join.
				if(retainGffs && retentionWorkers==1 && (a==0 || a==3)){
					retainCalls(retained, genome.id, ARMS[a], lines, exportCosts, exportTotals);
				}
				totals.get(a).add(score);
				counts.append(score.rows(ARMS[a], genome.id));
				cds[a]+=calls.size();
				searched[a]+=result.searched[a];
				reused[a]+=result.reused[a];
				accepted[a]+=result.accepts[a];
				features.append(genome.id).tab().append(ARMS[a]).tab().append(calls.features()).tab().append(calls.size())
					.tab().append(calls.features()-calls.size()).nl();
			}
			appendWork(workRows, genome.id, result.work);
			for(int v=0;v<2;v++){
				for(int f=0;f<WORK_FIELDS;f++){workTotal[v][f]+=result.work[v][f];}
			}
			timing.append(genome.id).tab().appendSlow(result.wall*1e-9).tab().appendSlow(result.cpu*1e-9).nl();
			if(retainGffs){
				inventory.append(genome.id).tab().append(retained.resolve(genome.id+".native.gff.gz").toString())
					.tab().append(retained.resolve(genome.id+".formula.gff.gz").toString()).nl();
			}
			if(!retainGffs || retentionWorkers==1 && !retainRna){results.set(i, null);}
		}
		appendWork(workRows, "ALL", workTotal);
		final ScoreCounts[] aggregate=new ScoreCounts[ARMS.length];
		final ByteBuilder summary=new ByteBuilder(CallGenesLengthScorecard.scoreHeader());
		for(int a=0;a<ARMS.length;a++){
			aggregate[a]=CallGenesLengthScorecard.combineScores(totals.get(a));
			final long[] grade=aggregate[a].grade();
			require(grade[5]+grade[7]==cds[a], "Whole-cohort query counts are not conserved");
			if(a>0){
				require(searched[a]==queries && reused[a]==queries && accepted[a]==accepts,
					"Whole-cohort replay evidence is not conserved");
			}
			summary.append(aggregate[a].rows(ARMS[a], "ALL"));
		}
		require(Arrays.equals(aggregate[0].rows("control", "ALL").toBytes(),
			aggregate[1].rows("control", "ALL").toBytes()),
			"Whole-cohort zero-effect control changed aggregate16-count vectors");
		write(output.resolve("by_genome.tsv"), counts);
		write(output.resolve("scorecard.tsv"), summary);
		write(output.resolve("worker_times.tsv"), timing);
		write(output.resolve("features.tsv"), features);
		write(output.resolve("refinement.tsv"), refinement);
		if(recordAlignmentWork){
			write(output.resolve("alignment_work.tsv"), workRows);
			compressWork(genomes);
		}
		if(retainGffs){
			if(retentionWorkers>1){
				retainParallel(genomes, results, retained, exportCosts, exportTotals);
			}
			else if(retainRna){retainRna(genomes, results);}
			write(output.resolve("retained_gff_costs.tsv"), exportCosts);
			write(output.resolve("retained_gff_cost.tsv"), new ByteBuilder("key\tvalue\nretained_files\t")
				.append(2L*genomes.size()).append("\nraw_payload_bytes\t").append(exportTotals[0])
				.append("\nstored_compressed_bytes\t").append(exportTotals[1]).append("\nexport_wall_seconds\t")
				.appendSlow(exportTotals[2]*1e-9).append("\nbgzf_files\t").append(exportTotals[3])
				.append("\npigz_files\t").append(exportTotals[4])
				.append("\ncompression_threshold_raw_bytes\t200000\ncompression_threshold_rule\traw>200000:BGZF9;raw<=200000:pigz-11\ncompression_threads\t1\nraw_files_retained\tfalse\n")
				.append("bgzip\t").append(retentionBgzip.toString())
				.append("\nbgzip_version\t").append(Data.BGZIP_VERSION)
				.append("\npigz\t").append(retentionPigz.toString())
				.append("\npigz_version\t").append(Data.PIGZ_VERSION)
				.append("\nexport_workers\t").append(retentionWorkers)
				.append("\nexport_wall_semantics\tsum_of_file_task_seconds").nl());
			validateRetentionCompression();
			write(output.resolve("retained_gffs.tsv"), inventory);
			for(int i=0; i<results.size(); i++){
				results.set(i, null);
			}
		}
	}

	/** Final tally compression is bounded and occurs after calling, never on search threads. */
	private void compressWork(ArrayList<Genome> genomes){
		new CallGenesRetainedExport.Tools(retentionBgzip, retentionPigz);
		final long start=System.nanoTime();
		final ArrayList<String> receipts=CallGenesLeanWorkers.run(genomes.size(),
			Math.min(retentionWorkers, genomes.size()), index->{
				final String id=genomes.get(index).id;
				final Path dir=output.resolve("alignment_work");
				final Path raw=dir.resolve(id+".tsv.tmp"), file=dir.resolve(id+".tsv.gz"),
					err=dir.resolve(id+".compression.stderr");
				final long began=System.nanoTime();
				try{
					final long bytes=Files.size(raw);
					final String[] command={retentionBgzip.toString(), "-c", "-f", "-l", "9", "-@", "1",
						raw.toString()};
					final Process process=new ProcessBuilder(command).redirectOutput(file.toFile())
						.redirectError(err.toFile()).start();
					final int exit=process.waitFor();
					require(exit==0 && Files.size(file)>0,
						"Alignment-work compression failed: "+id+" exit="+exit);
					//TODO: Verify decompressed bytes before retiring raw one-pass alignment-work output; exit0/nonempty is insufficient. Deferred packaging gap; accepted two-pass baseline disables this mode.
					Files.delete(raw);
					return id+"\t"+file+"\t"+bytes+"\t"+Files.size(file)+"\t"+
						(System.nanoTime()-began)*1e-9+"\t"+exit+"\t"+String.join(" ", command)+"\n";
				}catch(IOException problem){
					throw new UncheckedIOException("Alignment-work compression failed: "+id, problem);
				}
				catch(InterruptedException problem){
					Thread.currentThread().interrupt();
					throw new IllegalStateException("Interrupted alignment-work compression: "+id, problem);
				}
			});
		final ByteBuilder rows=new ByteBuilder("genome\toutput\traw_bytes\tstored_bytes\twall_seconds\texit\tcommand\n");
		for(String receipt:receipts){rows.append(receipt);}
		write(output.resolve("alignment_work_compression.tsv"), rows);
		phase("alignment_work_compression", System.nanoTime()-start);
	}

	/** One bounded pool serves CDS and RNA files without serializing independent kinds. */
	private void retainParallel(ArrayList<Genome> genomes, ArrayList<Result> results, Path retained,
		ByteBuilder costs, long[] totals){
		final Path rnaDirectory=retainRna ? prepareRna(results) : null;
		final int perGenome=retainRna ? 4 : 2, count=perGenome*genomes.size();
		final long start=System.nanoTime();
		final ArrayList<CallGenesRetainedExport.Receipt> receipts=CallGenesRetainedExport.write(count,
			Math.min(retentionWorkers, count), index->{
				final int ordinal=index/perGenome, slot=index%perGenome;
				final String arm=(slot&1)==0 ? "native" : "formula";
				if(slot<2){
					final ArrayList<GffLine> lines=results.get(ordinal).calls[slot==0 ? 0 : 3].lines();
					return new CallGenesRetainedExport.Payload(retained, genomes.get(ordinal).id, arm,
						CallGenesRetainedExport.cdsText(lines), lines.size());
				}
				final CallGenesRetainedExport.Rna rna=results.get(ordinal).rna[slot-2];
				return new CallGenesRetainedExport.Payload(rnaDirectory, genomes.get(ordinal).id, arm,
					rna.text(), rna.rows());
			}, new CallGenesRetainedExport.Tools(retentionBgzip, retentionPigz));
		phase("parallel_retained_export", System.nanoTime()-start);
		final ArrayList<CallGenesRetainedExport.Receipt> rnaReceipts=retainRna ?
			new ArrayList<CallGenesRetainedExport.Receipt>() : null;
		final ByteBuilder processes=new ByteBuilder("genome\tpolicy\tpopulation\toutput\texit\tcommand\n");
		for(CallGenesRetainedExport.Receipt receipt:receipts){
			if(receipt.file.getParent().equals(retained)){
				receipt.appendProcess(processes, "CDS");
				receipt.appendCost(costs);
				totals[0]=Math.addExact(totals[0], receipt.raw);
				totals[1]=Math.addExact(totals[1], receipt.stored);
				totals[2]=Math.addExact(totals[2], receipt.nanos);
				totals[receipt.bgzf ? 3 : 4]++;
			}else{
				require(retainRna && receipt.file.getParent().equals(rnaDirectory),
					"Unexpected retained feature output");
				rnaReceipts.add(receipt);
				receipt.appendProcess(processes, "non_CDS");
			}
		}
		write(output.resolve("retained_export_processes.tsv"), processes);
		if(retainRna){writeRnaTables(genomes, rnaReceipts);}
	}

	/** RNA is emitted from original Orfs and kept out of every CDS reconstruction/grade. */
	private Path prepareRna(ArrayList<Result> results){
		final Path directory=output.resolve("retained_rna_gffs");
		try{Files.createDirectory(directory);}
		catch(IOException failure){throw new UncheckedIOException("RNA retention directory must be fresh: "+directory,
			failure);}
		for(Result result:results){
			require(result!=null && result.rna!=null && result.rna.length==2,
				"Missing pre-compaction RNA capture");
			for(int slot=0;slot<2;slot++){
				final CallGenesLeanCalls calls=result.calls[slot==0 ? 0 : 3];
				require(result.rna[slot].rows()==calls.features()-calls.size(),
					"RNA retention must conserve every emitted non-CDS feature");
			}
		}
		return directory;
	}

	/** Serial retention still uses the isolated per-process backend for its added RNA files. */
	private void retainRna(ArrayList<Genome> genomes, ArrayList<Result> results){
		final Path directory=prepareRna(results);
		final long start=System.nanoTime();
		final ArrayList<CallGenesRetainedExport.Receipt> receipts=CallGenesRetainedExport.write(2*genomes.size(),
			Math.min(retentionWorkers, 2*genomes.size()), index->{
				final int ordinal=index/2, slot=index%2;
				final CallGenesRetainedExport.Rna rna=results.get(ordinal).rna[slot];
				return new CallGenesRetainedExport.Payload(directory, genomes.get(ordinal).id,
					slot==0 ? "native" : "formula", rna.text(), rna.rows());
			}, new CallGenesRetainedExport.Tools(retentionBgzip, retentionPigz));
		phase("retained_RNA_export", System.nanoTime()-start);
		writeRnaTables(genomes, receipts);
	}

	/** Publish separate RNA inventory, native-formatter provenance, costs and process receipts.
	 * Receipts arrive in genome/arm input order only after every active writer joins. */
	private void writeRnaTables(ArrayList<Genome> genomes, ArrayList<CallGenesRetainedExport.Receipt> receipts){
		require(receipts.size()==2*genomes.size(), "RNA inventory must cover both retained arms of every genome");
		final ByteBuilder inventory=new ByteBuilder("genome\tnative_rna_gff\tformula_rna_gff\n");
		final ByteBuilder provenance=new ByteBuilder("genome\tpolicy\tfeature_rows\traw_payload_bytes\tstored_compressed_bytes\tformatter\tpopulation\n");
		final ByteBuilder costs=new ByteBuilder("genome\tpolicy\traw_payload_bytes\tstored_compressed_bytes\texport_wall_seconds\tcompression\tcommand\n");
		final ByteBuilder processes=new ByteBuilder("genome\tpolicy\tpopulation\toutput\texit\tcommand\n");
		for(int i=0;i<genomes.size();i++){
			inventory.append(genomes.get(i).id).tab().append(receipts.get(2*i).file.toString()).tab()
				.append(receipts.get(2*i+1).file.toString()).nl();
		}
		for(CallGenesRetainedExport.Receipt receipt:receipts){
			receipt.appendCost(costs);
			receipt.appendProcess(processes, "non_CDS");
			provenance.append(receipt.genome).tab().append(receipt.arm).tab().append(receipt.rows)
				.tab().append(receipt.raw).tab().append(receipt.stored)
				.append("\tOrf.appendGff_before_LeanCalls\tall_emitted_non_CDS\n");
		}
		write(output.resolve("retained_rna_costs.tsv"), costs);
		write(output.resolve("retained_rna_provenance.tsv"), provenance);
		write(output.resolve("retained_rna_processes.tsv"), processes);
		write(output.resolve("retained_rna_gffs.tsv"), inventory);
	}

	/** Bounded diagnostic output under this run, retaining exactly the graded CDS order. */
	private static void retainCalls(Path directory, String genome, String arm, ArrayList<GffLine> lines,
			ByteBuilder costs, long[] totals){
		require(directory!=null && directory.isAbsolute() && genome.matches("g[0-9]{4}")
			&& (arm.equals("native") || arm.equals("formula")) && costs!=null && totals!=null && totals.length==5,
			"Retained GFF paths and costs require an accepted genome and policy");
		final long start=System.nanoTime();
		final Path file=directory.resolve(genome+"."+arm+".gff.gz");
		final ByteBuilder text=new ByteBuilder("##gff-version 3\n");
		for(GffLine line:lines){
			line.appendTo(text);
			text.nl();
		}
		final boolean bgzf=text.length()>200000;
		final int zip=ReadWrite.ZIPLEVEL, maxThreads=ReadWrite.MAX_ZIP_THREADS(),
			block=ReadWrite.PIGZ_BLOCKSIZE, iterations=ReadWrite.PIGZ_ITERATIONS;
		final boolean change=ReadWrite.ALLOW_ZIPLEVEL_CHANGE, useBgzf=ReadWrite.USE_BGZF,
			useBgzip=ReadWrite.USE_BGZIP, forceBgzip=ReadWrite.FORCE_BGZIP,
			preferBgzip=ReadWrite.PREFER_BGZIP, nativeBgzf=ReadWrite.ALLOW_NATIVE_BGZF,
			usePigz=ReadWrite.USE_PIGZ, forcePigz=ReadWrite.FORCE_PIGZ, raw=ReadWrite.RAWMODE,
			setThreads=ReadWrite.SET_ZIP_THREADS;
		final long compressed;
		try{
			require(!Files.exists(file) && !Files.isSymbolicLink(file), "Retained output must be fresh: "+file);
			//Only this post-worker writer is active; restore every changed compression setting.
			ReadWrite.ZIPLEVEL=bgzf ? 9 : 11;
			ReadWrite.ALLOW_ZIPLEVEL_CHANGE=false;
			ReadWrite.USE_BGZF=ReadWrite.FORCE_BGZIP=ReadWrite.PREFER_BGZIP=bgzf;
			ReadWrite.USE_BGZIP=ReadWrite.USE_PIGZ=true;
			ReadWrite.ALLOW_NATIVE_BGZF=false;
			ReadWrite.FORCE_PIGZ=!bgzf;
			ReadWrite.RAWMODE=false;
			ReadWrite.setZipThreads(1);
			ReadWrite.PIGZ_BLOCKSIZE=128;
			ReadWrite.PIGZ_ITERATIONS=-1;
			final ByteStreamWriter writer=new ByteStreamWriter(file.toString(), false, false, true);
			writer.start();
			try{writer.print(text);}
			finally{require(!writer.poisonAndWait(), "Retained compression writer/process failed: "+file);}
			compressed=Files.size(file);
		}catch(IOException failure){throw new UncheckedIOException("Retained GFF write failed: "+file, failure);}
		finally{
			ReadWrite.ZIPLEVEL=zip;
			ReadWrite.ALLOW_ZIPLEVEL_CHANGE=change;
			ReadWrite.USE_BGZF=useBgzf;
			ReadWrite.USE_BGZIP=useBgzip;
			ReadWrite.FORCE_BGZIP=forceBgzip;
			ReadWrite.PREFER_BGZIP=preferBgzip;
			ReadWrite.ALLOW_NATIVE_BGZF=nativeBgzf;
			ReadWrite.USE_PIGZ=usePigz;
			ReadWrite.FORCE_PIGZ=forcePigz;
			ReadWrite.RAWMODE=raw;
			ReadWrite.setZipThreads(maxThreads);
			ReadWrite.SET_ZIP_THREADS=setThreads;
			ReadWrite.PIGZ_BLOCKSIZE=block;
			ReadWrite.PIGZ_ITERATIONS=iterations;
		}
		final long nanos=System.nanoTime()-start;
		totals[0]=Math.addExact(totals[0], text.length());
		totals[1]=Math.addExact(totals[1], compressed);
		totals[2]=Math.addExact(totals[2], nanos);
		totals[bgzf ? 3 : 4]++;
		costs.append(genome).tab().append(arm).tab().append(text.length()).tab().append(compressed)
				.tab().appendSlow(nanos*1e-9).tab().append(bgzf ? "BGZF_level9" : "pigz_11")
				.tab().append(bgzf ? "bgzip -c -f -l 9 -@ 1" : "pigz -c -p 1 -11").nl();
	}

	/** Match each PATH-selected external compressor to its explicit path, without stamps. */
	private void validateRetentionCompression(){
		try{
			retentionTool(retentionBgzip, "bgzip");
			retentionTool(retentionPigz, "pigz");
			require(Data.BGZIP() && Data.BGZIP_VERSION_levelFlag && Data.BGZIP_VERSION_threadsFlag
				&& Data.PIGZ() && Data.PIGZ_VERSION_23plus,
				"Compressors must support BGZF level9/one thread and pigz -11; fallback is forbidden");
		}catch(IOException failure){throw new UncheckedIOException("Retention compressor preflight failed", failure);}
	}
	private static void retentionTool(Path binary, String name)throws IOException{
		require(binary!=null && binary.equals(binary.toRealPath()) && Files.isRegularFile(binary)
			&& Files.isExecutable(binary), "Explicit canonical retention compressor required: "+name);
		final String search=System.getenv("PATH");
		require(search!=null, "PATH is required for the existing compression launcher");
		for(String directory:search.split(java.io.File.pathSeparator, -1)){
			final Path candidate=Paths.get(directory.isEmpty() ? "." : directory).resolve(name);
			if(Files.isRegularFile(candidate) && Files.isExecutable(candidate)){
				require(candidate.toRealPath().equals(binary),
					"PATH-selected compressor differs from its binding: "+name);
				return;
			}
		}
		throw new IllegalArgumentException("Missing PATH-selected compressor: "+name);
	}

	/** Use the accepted StartAccounting input table directly; do not remap or sample the cohort. */
	private ArrayList<Genome> readInputs(){
		final ArrayList<Genome> rows=new ArrayList<Genome>();
		final HashSet<String> ids=new HashSet<String>();
		final ByteFile reader=ByteFile.makeByteFile1(input.toString(), false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			require(Arrays.equals(reader.nextLine(), CallGenesBenchmarkSupport.INPUT_HEADER.getBytes(StandardCharsets.US_ASCII)),
				"Require accepted13-column StartAccounting inputs");
			for(byte[] line=reader.nextLine();line!=null;line=reader.nextLine()){
				lp.set(line);
				require(lp.terms()==13, "Malformed saved input row");
				final String id=lp.parseString(0);
				require(id.matches("g[0-9]{4}") && ids.add(id), "Unsafe or repeated genome key");
				rows.add(new Genome(id, path(lp.parseString(1)), pinLabel(lp.parseString(2)),
					path(lp.parseString(3)), pinLabel(lp.parseString(4)), path(lp.parseString(9))));
			}
		}finally{require(!reader.close(), "Saved input table failed to read");}
		require(rows.size()==expected, "Input table must cover the exact requested cohort");
		return rows;
	}

	/** The accepted diagnostic layout binds each saved scanner status to its caller recipe. */
	private void checkScannerRecipes(ArrayList<Genome> genomes, StopTailCascadeResources resources)throws Exception{
		final HashSet<Path> seen=new HashSet<Path>();
		for(Genome genome:genomes){
			require(genome.annotations.getFileName().toString().equals("diagnostic_annotations.tsv")
					&& genome.annotations.getParent().getFileName().toString().equals(genome.id),
					"Require retained per-genome diagnostic binding");
			final Path recipe=genome.annotations.getParent().getParent().resolve("recipe.txt");
			if(!seen.add(recipe)){continue;}
			final HashMap<String, String> values=new HashMap<String, String>();
			for(String line:Files.readAllLines(recipe, StandardCharsets.UTF_8)){
				final String[] kv=Parse.splitOnFirst(line, '=');
				require(kv.length==2 && values.put(kv[0], kv[1])==null, "Invalid saved scanner recipe");
			}
			require("11".equals(values.get("scanner_table")) && "100".equals(values.get("query_cap"))
				&& "true".equals(values.get("firststart")) && "true".equals(values.get("searchall"))
				&& "false".equals(values.get("fixed100")) && "true".equals(values.get("diagnostics"))
				&& resources.queryLengths().equals(values.get("query_lengths")),
				"Saved locus scanner contract differs from lean caller");
		}
	}
	private static String contigName(String name){
		int end=0;
		while(end<name.length() && name.charAt(end)!=' ' && name.charAt(end)!='\t'){end++;}
		require(end>0, "Empty emitted contig name");
		return name.substring(0, end);
	}
	private void phase(String name, long nanos){phases.append(name).tab().appendSlow(nanos*1e-9).nl();}
	private static String text(String value){return value==null ? "none" : value;}
	private static String pinLabel(String value){
		require(value.matches("[0-9a-f]{20}|unmeasured"), "Require sha80 or explicit unmeasured identity label");
		return value;
	}
	private static String measuredPin(String value){
		require(value.matches("[0-9a-f]{20}"), "Selection manifests require measured sha80");
		return value;
	}
	private static Path path(String value){
		final Path p=Paths.get(value);
		require(p.isAbsolute() && p.normalize().equals(p), "Require normalized absolute paths");
		return p;
	}
	/** Tables are already buffered; synchronous output has no writer thread to outlive failure. */
	private static void write(Path path, ByteBuilder text){
		try{Files.write(path, text.toBytes());}
		catch(IOException failure){throw new UncheckedIOException("Lean table write failed: "+path, failure);}
	}
	private static void require(boolean value, String message){if(!value){throw new IllegalArgumentException(message);}}
	private static final class Genome{
		Genome(String id_, Path fasta_, String fp_, Path reference_, String rp_, Path annotations_){
			id=id_;
			fasta=fasta_;
			fastaPin=fp_;
			reference=reference_;
			referencePin=rp_;
			annotations=annotations_;
		}
		final String id, fastaPin, referencePin;
		final Path fasta, reference, annotations;
	}
	private static final class Result{
		Result(String id_){
			id=id_;
			for(int a=0;a<calls.length;a++){calls[a]=new CallGenesLeanCalls();}
		}
		final String id;
		final CallGenesLeanCalls[] calls=new CallGenesLeanCalls[ARMS.length];
		CallGenesRetainedExport.Rna[] rna;
		final long[][] work=new long[2][WORK_FIELDS];
		final long[] searched=new long[ARMS.length], reused=new long[ARMS.length], accepts=new long[ARMS.length];
		long queries, accepted, wall, cpu;
		long firstQueries, firstAccepted, cachedQueries, recomputedQueries;
		final long[] firstCds=new long[ARMS.length];
	}
	private static final int WORK_FIELDS=11;
	private static final String WORK_STOP_HEADER="genome\tcontig\tstop_left_zero_based\tstrand_0plus_1minus\tview\thas_ranked_shortlist\talignment_calls\talignment_accepts\tshort_input_skips\tunavailable_calibration\tcached_provider_calls\tshortcut_accepts\talignment_stage_visits\tquery_length_aa\n";
	private static final String WORK_HEADER="genome\tview\teligible_queries\tqueries_with_shortlist\talignment_calls\talignment_accepts\tshort_input_skips\tunavailable_calibration\tcached_provider_calls\tshortcut_accepts\tqueries_with_alignment\tqueries_with_alignment_accept\talignment_stage_visits\n";
	private static final String[] ARMS={"native", "zero_effect", "bonus400", "formula"};
	private final ByteBuilder phases=new ByteBuilder("phase\twall_seconds\n");
	private Path input, loci, pgm, resource, information, fmn, output;
	private Path selectionCachePath, modelInventoryPath;
	private Path modelRelocationMap, modelRelocationRoot;
	private String selectionCachePin, modelInventoryPin;
	private boolean verifySelectionCache;
	private boolean retainGffs, retainRna, recordAlignmentWork;
	private int retentionWorkers=1;
	private Path retentionBgzip, retentionPigz;
	private boolean rawScoreI0Override;
	private String lociPin, resourcePin;
	private CallGenesReplayRnaPolicy.ModelPolicy policy;
	private final FamilyBonusFormula formula;
	private final String[] invocation;
	private int expected, workers, searchThreads, cpus, top, minHits=-1;
	private int passes=1;
	private boolean reusePassEvidence=true;
	private float identityReduction;
}
