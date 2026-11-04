package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.Data;
import map.ObjectIntMap;
import parse.Parse;
import parse.Parser;
import prok.NcrnaBoundaryComparison.Change;
import prok.NcrnaBoundaryComparison.Result;
import prok.NcrnaSupport.Graph;
import prok.NcrnaSupport.Library;
import prok.NcrnaTextIO.Output;
import prok.NucleotideHbmCodec.Model;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/** Serial regeneration with atomic member-support acceptance and caller rollback.
 * The reducer owns CLI parsing and the complete corpus. This stage owns candidate
 * pairs, assignments and output files; all alignment uses NcrnaSupport, all locus
 * grading uses NcrnaBoundaryComparison, and the public consensus builder owns
 * consensus construction. The stage never invokes a development driver.
 * Each instance runs once; pairs, graph and status are read after process returns.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaModelRegeneration{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Expands a scientific consensus recipe; input, grouping and budgets belong to this stage. */
	NcrnaModelRegeneration(final String config, final long heapMiB_, final int maxCalls_,
			final int threads_){
		if(config==null || heapMiB_<1 || maxCalls_<1 || threads_<1){
			throw new IllegalArgumentException("Regeneration needs consensusconfig, positive consensusheap, maxregencalls and threads");
		}
		heapMiB=heapMiB_;
		maxCalls=maxCalls_;
		threads=threads_;
		consensusArgs=Parser.parseConfig(new String[]{"config="+NcrnaCallerProcess.configurationPath(new File(config))});
		for(final String arg:consensusArgs){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			if(split.length!=2 || !(a.equals("passes") || a.equals("minid") || a.equals("minidentity")
					|| a.equals("endtrim") || a.equals("endtrimfrac") || a.equals("stochastic")
					|| a.equals("stochasticseed"))){
				throw new IllegalArgumentException("consensusconfig accepts passes, minid, endtrim, stochastic and stochasticseed only: "+arg);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns verified pairs/graph, retaining the stage-entry baseline on a failed stage.
	 * Original greedy coverage and stage-entry retention are separate populations.
	 * The exact existing caller matching is borrowed, never reconstructed or rerun.
	 */
	void process(final String input, final NcrnaModelPairs baseline, final Graph initial,
			final NcrnaCallerProtection caller, final File directory, final PrintStream diagnostics){
		root=directory;
		if(!root.mkdir()){throw new IllegalArgumentException("Regeneration directory must be fresh: "+root);}
		final Library genes=NcrnaSupport.loadLibrary(input);
		if(!Arrays.equals(initial.modelNames, baseline.names) || initial.queryNames.length!=genes.names.length){
			throw new IllegalArgumentException("Regeneration baseline must use the complete corpus and selected model axis");
		}
		final ObjectIntMap<String> geneIndex=index(genes.names);
		final int[] initialOwners=assignments(initial, geneIndex);
		final int[] owners=initialOwners.clone();
		final NcrnaModelPairs current=baseline.copy();
		serial(current, genes, geneIndex, owners, initial.cutoff, diagnostics);
		final Graph series=scan(genes, null, current.library(), initial.cutoff);
		final int[] seriesOwners=assignments(series, geneIndex);
		for(int gene=0; gene<owners.length; gene++){
			if(owners[gene]!=seriesOwners[gene]){finalMigrations++;}
		}
		writeAssignments("series.assignments.tsv", genes.names, baseline.names, initialOwners, owners, seriesOwners);
		final Result baselineCalls=caller.selectedResult();
		if(baselineCalls==null){throw new IllegalStateException("Regeneration requires the verified caller result from reduction");}
		final boolean callerPass=protect(current, baseline, caller, baselineCalls);
		// Restored pairs can change support: scan once after restoration, otherwise reuse the series pass.
		final Graph candidate=restoredPairs==0 ? series : scan(genes, null, current.library(), initial.cutoff);
		final int[] finalOwners=assignments(candidate, geneIndex);
		writeAssignments("candidate.assignments.tsv", genes.names, baseline.names, initialOwners, seriesOwners, finalOwners);
		final Retention retention=selectResult(baseline, current, initial, candidate, baselineCalls,
			callerPass, initialOwners, finalOwners);
		final int[] exportedOwners=assignments(graph, geneIndex);
		writeAssignments("final.assignments.tsv", genes.names, baseline.names, initialOwners, initialOwners, exportedOwners);
		try(Output output=new Output(path("exported.caller.tsv"))){
			output.writer.print("new_call\tmatched\tunmatched\tmissed\texact_5prime\texact_3prime\tboth_exact\tlost\tgained\n");
			writeCall(output, new ByteBuilder(), pairs==baseline ? 0 : newCalls, baselineCalls, selectedCalls);
		}
		try(Output output=new Output(path("exported.models.tsv"))){
			output.writer.print("model\tfinal_disposition\n");
			final ByteBuilder row=new ByteBuilder();
			for(int model=0; model<pairs.names.length; model++){
				row.clear();
				final String decision=pairs==baseline ? "BASELINE_FALLBACK" : pairs.same(baseline, model) ? "ORIGINAL_PAIR" : "REGENERATED_PAIR";
				output.writer.print(row.append(pairs.names[model]).tab().append(decision).nl());
			}
		}
		try(Output output=new Output(path("summary.tsv"))){
			final ByteBuilder row=new ByteBuilder("status\tmodels\tkept\treverted\tempty\tmember_migrations\tfinal_migrations\tnew_calls\tcaller_restorations\tcorpus_lost\tcorpus_gained\talignments\tinitial_supported\tinitial_unsupported\tcandidate_supported\n");
			row.append(status).tab().append(baseline.names.length).tab().append(kept).tab().append(reverted).tab().append(empty).tab();
			row.append(memberMigrations).tab().append(finalMigrations).tab().append(newCalls).tab().append(restoredPairs).tab();
			row.append(retention.lost).tab().append(retention.gained).tab().append(alignments).tab().append(retention.before).tab();
			output.writer.print(row.append(initialOwners.length-retention.before).tab().append(retention.after).nl());
		}
		diagnostics.println(status+"; new caller passes="+newCalls+", lost support="+retention.lost);
	}

	/** Selects pairs, graph and caller result together after both retention checks.
	 * Owner arrays were joined to the same named gene axis by assignments; counts
	 * alone cannot substitute for retaining each initially supported gene.
	 */
	private Retention selectResult(final NcrnaModelPairs baseline, final NcrnaModelPairs current,
			final Graph initial, final Graph candidate, final Result baselineCalls,
			final boolean callerPass, final int[] initialOwners, final int[] finalOwners){
		final Retention retention=new Retention(initialOwners, finalOwners);
		if(!callerPass || retention.lost>0){
			pairs=baseline;
			graph=initial;
			selectedCalls=baselineCalls;
			if(retention.lost>0){status=(status==null ? "Baseline retained" : status)+"; regeneration lost corpus support";}
		}else{
			pairs=current;
			final boolean[] all=new boolean[baseline.names.length];
			Arrays.fill(all, true);
			graph=candidate.subset(all);
			status="Regeneration passed caller and corpus retention";
		}
		return retention;
	}

	/** Named-gene support conservation after caller-triggered paired restorations. */
	private static final class Retention{
		Retention(final int[] initialOwners, final int[] finalOwners){
			assert(initialOwners.length==finalOwners.length)
				: "Retention compares one unchanged named gene axis: before="+initialOwners.length+", after="+finalOwners.length;
			int lost_=0, gained_=0, before_=0, after_=0;
			for(int gene=0; gene<initialOwners.length; gene++){
				if(initialOwners[gene]>=0){before_++;}
				if(finalOwners[gene]>=0){after_++;}
				if(initialOwners[gene]>=0 && finalOwners[gene]<0){lost_++;}
				if(initialOwners[gene]<0 && finalOwners[gene]>=0){gained_++;}
			}
			lost=lost_;
			gained=gained_;
			before=before_;
			after=after_;
		}
		/** Genes losing/gaining support and total supported genes before/after the candidate. */
		final int lost, gained, before, after;
	}

	/** Visits fixed model indices once; a failed candidate never mutates pairs or assignments. */
	private void serial(final NcrnaModelPairs current, final Library genes,
			final ObjectIntMap<String> geneIndex, final int[] owners, final float cutoff, final PrintStream diagnostics){
		try(Output output=new Output(path("models.tsv"))){
			output.writer.print("model\tdecision\tmembers\tunsupported\tmigrations\n");
			final ByteBuilder row=new ByteBuilder();
			for(int model=0; model<current.names.length; model++){
				final IntList members=new IntList();
				for(int gene=0; gene<owners.length; gene++){if(owners[gene]==model){members.add(gene);}}
				String decision="EMPTY_UNCHANGED";
				int unsupported=0, moved=0;
				if(members.isEmpty()){
					empty++;
				}else{
					final File work=fresh("model_"+model);
					final String memberPath=new File(work, "members.fa").getPath();
					writeMembers(memberPath, genes, members);
					final Model candidate=build(work, current.names[model]);
					if(candidate==null){
						reverted++;
						decision="NO_CANDIDATE_REVERTED";
					}else{
						final NcrnaModelPairs proposal=current.copy();
						proposal.models[model]=candidate;
						proposal.headers[model]=NcrnaModelPairs.header(candidate);
						final Graph score=scan(genes, members, proposal.library(), cutoff);
						if(score.queryNames.length!=members.size){throw new IllegalStateException("Member scoring omitted genes");}
						final int[] ids=new int[score.queryNames.length];
						for(int q=0; q<ids.length; q++){
							ids[q]=geneIndex.get(score.queryNames[q]);
							if(ids[q]<0 || owners[ids[q]]!=model){throw new IllegalStateException("Member score has an unexpected gene");}
							if(score.starts[q]==score.starts[q+1]){unsupported++;}
						}
						if(unsupported>0){
							reverted++;
							decision="UNSUPPORTED_REVERTED";
						}else{
							// Commit only after every member passes. Nothing asynchronous observes this state.
							current.models[model]=candidate;
							current.headers[model]=proposal.headers[model];
							for(int q=0; q<ids.length; q++){
								if(owners[ids[q]]!=score.bestModels[q]){moved++;}
								owners[ids[q]]=score.bestModels[q];
							}
							kept++;
							memberMigrations+=moved;
							decision="KEPT";
						}
					}
				}
				row.clear();
				output.writer.print(row.append(current.names[model]).tab().append(decision).tab().append(members.size).tab()
					.append(unsupported).tab().append(moved).nl());
				diagnostics.println("Regenerated model "+(model+1)+"/"+current.names.length+": "+decision);
			}
		}
		assert(kept+reverted+empty==current.names.length)
			: "Every fixed model index receives exactly one serial decision: models="+current.names.length+
			", kept="+kept+", reverted="+reverted+", empty="+empty;
	}

	/** Runs only new calls, restoring distinct named baseline winners between proposals. */
	private boolean protect(final NcrnaModelPairs current, final NcrnaModelPairs baseline,
			final NcrnaCallerProtection caller, final Result original){
		final ObjectIntMap<String> indices=index(baseline.names);
		try(Output output=new Output(path("caller_results.tsv")); Output restorations=new Output(path("caller_restorations.tsv"))){
			output.writer.print("new_call\tmatched\tunmatched\tmissed\texact_5prime\texact_3prime\tboth_exact\tlost\tgained\n");
			restorations.writer.print("after_call\tmodel\n");
			final ByteBuilder row=new ByteBuilder();
			writeCall(output, row, 0, original, original);
			for(int pass=1; pass<=maxCalls; pass++){
				final File work=fresh("caller_"+pass);
				current.write(work, "models.hbmc");
				final Result result=caller.evaluate(work, "models.hbmc", original);
				newCalls++;
				writeCall(output, row, pass, original, result);
				NcrnaBoundaryComparison.writeTransitions(new File(work, "transitions.tsv").getPath(), original, result);
				final Change change=new Change(original, result);
				if(change.lost==0){
					selectedCalls=result;
					return true;
				}
				if(pass==maxCalls){
					status="Baseline retained: regeneration caller budget exhausted";
					return false;
				}
				int restored=0;
				for(int reference=0; reference<original.matches.length; reference++){
					if(original.matches[reference]<0 || result.matches[reference]>=0){continue;}
					final String name=original.call(reference).model;
					final int model=indices.get(name);
					if(model<0){throw new IllegalStateException("Verified baseline winner is absent from its library: "+name);}
					if(current.restore(baseline, model)){
						restored++;
						row.clear();
						restorations.writer.print(row.append(pass).tab().append(name).nl());
					}
				}
				restoredPairs+=restored;
				if(restored==0){
					status="Baseline retained: lost caller winners were already original; no progress";
					return false;
				}
			}
		}
		throw new AssertionError("Positive finite caller budget must return success or a reported fallback");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Native Operations      ----------------*/
	/*--------------------------------------------------------------*/

	/** Invokes the public builder in one owned child with one model group and fixed input/output routes. */
	private Model build(final File work, final String name){
		final String script=Data.findPath("?trnaconsensus.sh", false);
		if(script==null){throw new IllegalStateException("Cannot locate trnaconsensus.sh");}
		final String fasta=new File(work, "candidate.fa").getPath();
		final String hbm=new File(work, "candidate.hbmc").getPath();
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(script);
		command.add("simd="+(Shared.SIMD ? "t" : "f"));
		for(final String arg:consensusArgs){command.add(arg);}
		command.add("in="+new File(work, "members.fa").getPath());
		command.add("out="+fasta);
		command.add("outmodel="+hbm);
		command.add("cluster=f");
		command.add("mingroup=1");
		command.add("mincluster=1");
		command.add("recruit=f");
		command.add("prune=0");
		command.add("reassignrounds=0");
		command.add("t=1");
		command.add("ow=f");
		NcrnaProcess.addJvmOptions(command, heapMiB);
		NcrnaProcess.run(command, null, "Consensus regeneration for "+name);
		// HBM is the builder's paired consensus plus counts; no need to reread its redundant FASTA output.
		try(NucleotideHbmCodec.Reader reader=new NucleotideHbmCodec.Reader(hbm)){
			final Model model=reader.next();
			if(reader.next()!=null){throw new IllegalStateException("One member group produced more than one consensus model");}
			if(model==null){return null;}
			Tools.toUpperCase(model.bases);
			return NcrnaModelPairs.rename(model, name);
		}
	}

	/** Scores retained genes directly; the written member FASTA belongs only to the child builder. */
	private Graph scan(final Library genes, final IntList members, final Library library, final float cutoff){
		final NcrnaSupport.ScanResult result=NcrnaSupport.scanWithStats(genes, members, library, cutoff, threads);
		alignments+=result.alignments;
		return result.graph;
	}

	/** Retains stable gene IDs for the engine join and rejects unintended anticodon subgroups. */
	private static void writeMembers(final String path, final Library genes, final IntList members){
		final ByteBuilder row=new ByteBuilder();
		try(Output output=new Output(path)){
			for(int i=0; i<members.size; i++){
				final int gene=members.get(i);
				if(TrnaConsensusBuilder.parseAnticodon(genes.names[gene])!=null){
					throw new IllegalArgumentException("Regeneration gene IDs must not create anticodon subgroups: "+genes.names[gene]);
				}
				row.clear();
				output.writer.print(row.append('>').append(genes.names[gene]).nl().append(genes.bases[gene]).nl());
			}
		}
	}

	private static ObjectIntMap<String> index(final String[] names){
		final ObjectIntMap<String> result=new ObjectIntMap<String>(String.class);
		for(int i=0; i<names.length; i++){result.put(names[i], i);}
		return result;
	}

	/** Joins exact gene identities; equal row counts alone cannot establish the population. */
	private static int[] assignments(final Graph graph, final ObjectIntMap<String> indices){
		final int[] owners=new int[graph.queryNames.length];
		Arrays.fill(owners, -2);
		for(int i=0; i<owners.length; i++){
			final int gene=indices.get(graph.queryNames[i]);
			if(gene<0 || gene>=owners.length || owners[gene]!=-2){throw new IllegalArgumentException("Support graph gene identities differ from the corpus");}
			owners[gene]=graph.starts[i]<graph.starts[i+1] ? graph.bestModels[i] : -1;
		}
		for(final int owner:owners){if(owner==-2){throw new IllegalArgumentException("Support graph omitted a corpus gene");}}
		return owners;
	}

	private void writeAssignments(final String filename, final String[] genes, final String[] models,
			final int[] initial, final int[] before, final int[] after){
		try(Output output=new Output(path(filename))){
			output.writer.print("gene\tinitial_model\tprevious_model\tcurrent_model\n");
			final ByteBuilder row=new ByteBuilder();
			for(int i=0; i<genes.length; i++){
				row.clear();
				output.writer.print(row.append(genes[i]).tab().append(modelName(models, initial[i])).tab()
					.append(modelName(models, before[i])).tab().append(modelName(models, after[i])).nl());
			}
		}
	}

	private static void writeCall(final Output output, final ByteBuilder row, final int pass, final Result before, final Result after){
		final Change change=new Change(before, after);
		row.clear();
		output.writer.print(row.append(pass).tab().append(after.matched).tab().append(after.unmatched).tab().append(after.missed).tab()
			.append(after.exactStart).tab().append(after.exactStop).tab().append(after.exactBoth).tab().append(change.lost).tab().append(change.gained).nl());
	}

	private File fresh(final String name){
		final File directory=new File(root, name);
		if(!directory.mkdir()){throw new IllegalArgumentException("Stage output directory must be fresh: "+directory);}
		return directory;
	}

	private String path(final String name){return new File(root, name).getPath();}
	private static String modelName(final String[] names, final int model){return model<0 ? "NONE" : names[model];}
	/** All member and whole-corpus alignment work performed by this stage. */
	long alignments(){return alignments;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Expanded scientific consensus options; input, grouping and budgets are supplied by this stage. */
	final String[] consensusArgs;
	/** Heap budget for the one consensus-builder child, in MiB. */
	final long heapMiB;
	/** New caller-pass budget, excluding the reused protected baseline. */
	final int maxCalls;
	/** Alignment-worker budget; consensus rebuilding itself is single-threaded. */
	final int threads;
	/** Fresh directory owned by this one regeneration invocation. */
	private File root;
	/** Serial-stage decisions, counted before caller rollback can restore accepted pairs. */
	private int kept, reverted, empty;
	/** New whole-genome caller passes and baseline pairs restored after caller losses. */
	private int newCalls, restoredPairs;
	/** Serial assignment-change events, then changes from serial ownership to the first whole-corpus rescan. */
	private long memberMigrations, finalMigrations;
	/** Total gene-to-model alignments across member and whole-corpus scans. */
	private long alignments;
	/** Verified caller matching corresponding to the selected exported pairs. */
	private Result selectedCalls;
	/** Final paired library: accepted candidate or the unchanged protected baseline. */
	NcrnaModelPairs pairs;
	/** Support graph belonging to pairs on the same gene/model axes. */
	Graph graph;
	/** Human-readable final selection or fallback reason. */
	String status;
}
