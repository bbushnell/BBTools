package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import map.ObjectIntMap;
import parse.Parse;
import parse.Parser;
import prok.NcrnaBoundaryComparison.Result;
import prok.NcrnaLocusComparison.Locus;
import prok.NcrnaSupport.Library;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Tests model proposals with ordinary CallGenes and restores original winners.
 * One child JVM runs at a time through the installed public caller script.
 * NcrnaProcess owns the child process group and stops it on interruption.
 * No in-process replay cache is used. Files retain the
 * executed caller configuration, full GFF and reference-level losses.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCallerProtection{

	/** Loads the fixed caller policy, genome paths and reference loci for one reduction.
	 * The caller validates positive thread and call budgets before construction.
	 * @param family_ One registered member with paired-model override support
	 * @param truthFamily Label selecting this family's rows in the supplied truth
	 * @param heap_ Native binary-size JVM heap value, such as 1g
	 * @param maxCalls_ Total initial-protection caller budget, including its baseline
	 */
	NcrnaCallerProtection(final String config, final String genomes, final String truth,
			final String family_, final String truthFamily, final int threads_,
			final String heap_, final int maxCalls_){
		final NcrnaFamilyRegistry.Entry entry=NcrnaFamilyRegistry.modelOverrideMember(family_);
		family=entry.id;
		displayFamily=NcrnaFamilyRegistry.displayName(family);
		baseArgs=Parser.parseConfig(new String[]{"config="+config});
		validatePolicy(baseArgs);
		genomePaths=NcrnaCallerEvaluation.readGenomes(genomes);
		references=NcrnaLocusComparison.readTruth(truth, truthFamily);
		threads=threads_;
		heap=heap_;
		maxCalls=maxCalls_;
		if(Parse.parseKMG(heap)<=0){throw new IllegalArgumentException("callerheap must be a positive JVM heap size");}
	}

	/** Returns the first protected reduced set, or the verified original fallback. */
	boolean[] protect(final Library library, final String ref, final String hbm,
			final boolean[] proposal, final File directory){
		assert(proposal.length==library.names.length) : "Proposal and original winning-model lookup share one library axis";
		final ObjectIntMap<String> modelIndex=new ObjectIntMap<String>(String.class);
		for(int model=0; model<library.names.length; model++){
			modelIndex.put(library.names[model], model);
		}
		final File baselineDir=fresh(directory, "baseline");
		final ArrayList<Locus> baseline=call(ref, hbm, baselineDir);
		final Result baselineResult=NcrnaBoundaryComparison.grade(references, baseline, null);
		selectedResult=baselineResult;
		final int[] originalMatches=baselineResult.matches;
		if(matched(originalMatches)==0){throw new IllegalArgumentException("No baseline references were detected; check family activation and truthfamily");}
		final boolean[] kept=proposal.clone();
		try(Output output=new Output(new File(directory, "caller_results.tsv").getPath())){
			final ByteStreamWriter report=output.writer;
			final ByteBuilder row=new ByteBuilder();
			report.print("iteration\tmodels\tcalls\tmatched\tunmatched_calls\tmissed_references\tprotected_lost\tnewly_matched\tboth_ends_exact\texact_5prime\texact_3prime\n");
			writeResult(report, row, 0, library.names.length, baseline, originalMatches, originalMatches);
			for(int iteration=1; iteration<maxCalls; iteration++){
				final File run=fresh(directory, "candidate_"+iteration);
				ReduceNcrnaModels.writeLibrary(library, hbm, kept, run);
				final ArrayList<Locus> calls=call(new File(run, "consensus.fa").getPath(),
					new File(run, ReduceNcrnaModels.modelFilename(hbm)).getPath(), run);
				final Result graded=NcrnaBoundaryComparison.grade(references, calls, baselineResult);
				final int[] current=graded.matches;
				writeResult(report, row, iteration, count(kept), calls, current, originalMatches);
				int lost=0;
				int restored=0;
				try(Output lossOutput=new Output(new File(run, "lost.tsv").getPath())){
					final ByteStreamWriter losses=lossOutput.writer;
					losses.print("reference_index\tseqid\tstart\tstop\tstrand\toriginal_model\tadded_back\n");
					for(int reference=0; reference<current.length; reference++){
						if(originalMatches[reference]<0 || current[reference]>=0){continue;}
						lost++;
						final String winner=baseline.get(originalMatches[reference]).model;
						final int model=modelIndex.get(winner);
						if(model<0){throw new IllegalStateException("Baseline winning model is absent from original consensus: "+winner);}
						final boolean added=!kept[model];
						if(added){
							kept[model]=true;
							restored++;
						}
						final Locus locus=references.get(reference);
						row.clear();
						row.append(reference).tab().append(locus.seqid).tab().append(locus.start).tab().append(locus.stop).tab();
						row.append(locus.strand==0 ? '+' : '-').tab().append(winner).tab().append(added ? 1 : 0).nl();
						losses.print(row);
					}
				}
				if(lost==0){
					selectedResult=graded;
					resultStatus="Caller protection passed";
					return kept;
				}
				if(restored==0){
					resultStatus="Original library retained: lost loci already had their original winning models";
					return original(library.names.length);
				}
			}
			resultStatus="Original library retained: caller budget exhausted before a reduced set passed";
			return original(library.names.length);
		}
	}

	String status(){return resultStatus;}

	/** The exact verified matching for the returned subset, reused by regeneration. */
	Result selectedResult(){return selectedResult;}

	/** Grades a new paired-library proposal on the same reference objects and caller policy. */
	Result evaluate(final File directory, final String modelFilename, final Result baseline){
		return NcrnaBoundaryComparison.grade(references, call(new File(directory, "consensus.fa").getPath(),
			new File(directory, modelFilename).getPath(), directory), baseline);
	}

	/** Invokes CallGenes in a fresh JVM, preserving ordinary genome-level behavior. */
	private ArrayList<Locus> call(final String ref, final String hbm, final File directory){
		final File config=new File(directory, "caller.config");
		final File gff=new File(directory, "calls.gff");
		NcrnaCallerProcess.writeConfiguration(config, baseArgs, genomePaths, gff, family, threads,
			"ncrnaconsensus="+new File(ref).getAbsolutePath(), "ncrnamodels="+new File(hbm).getAbsolutePath());
		NcrnaCallerProcess.runInstalled(config, Parse.parseKMGBinary(heap)/(1024L*1024));
		return NcrnaLocusComparison.readCalls(gff.getPath(), displayFamily);
	}

	/** Keeps the caller policy scientific; this driver owns inputs, outputs and budgets. */
	public static void validatePolicy(final String[] args){
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			if(split.length!=2 || arg.indexOf('\n')>=0 || arg.indexOf('\r')>=0 || a.startsWith("-") ||
					a.startsWith("xmx") || a.startsWith("xms") || a.startsWith("xmn") || a.startsWith("xss") ||
					a.equals("ea") || a.equals("da") || a.equals("eoom") || a.equals("exitonoutofmemoryerror") ||
					a.equals("in") || a.equals("ref") || a.equals("infna") || a.equals("fnain") || a.equals("fna") ||
					a.equals("in1") || a.equals("in2") || a.equals("input") || a.equals("input1") || a.equals("input2") ||
					a.startsWith("out") || a.endsWith("out") || a.endsWith("log") || a.endsWith("diag") ||
					a.endsWith("hist") || a.equals("amino") || a.equals("json") ||
					a.equals("t") || a.equals("threads") || a.equals("workers") ||
					a.equals("reads") || a.equals("maxreads") || a.equals("samplerate") ||
					a.equals("stats") || a.equals("genehist") || a.equals("ncrnaconsensus") || a.equals("ncrnamodels")){
				throw new IllegalArgumentException("callerconfig must contain scientific options only; driver owns "+a);
			}
		}
	}

	private void writeResult(final ByteStreamWriter output, final ByteBuilder bb, final int iteration,
			final int modelCount, final ArrayList<Locus> calls, final int[] matches, final int[] baseline){
		int tp=0, lost=0, gained=0, exact=0, exactStart=0, exactStop=0;
		for(int i=0; i<matches.length; i++){
			if(matches[i]>=0){
				tp++;
				final Locus call=calls.get(matches[i]), reference=references.get(i);
				if(call.start==reference.start && call.stop==reference.stop){exact++;}
				if(NcrnaBoundaryComparison.offset(call, reference, true)==0){exactStart++;}
				if(NcrnaBoundaryComparison.offset(call, reference, false)==0){exactStop++;}
				if(baseline[i]<0){gained++;}
			}else if(baseline[i]>=0){lost++;}
		}
		bb.clear();
		bb.append(iteration).tab().append(modelCount).tab().append(calls.size()).tab().append(tp).tab();
		bb.append(calls.size()-tp).tab().append(references.size()-tp).tab().append(lost).tab().append(gained).tab().append(exact).tab();
		bb.append(exactStart).tab().append(exactStop).nl();
		output.print(bb);
	}

	private static File fresh(final File parent, final String name){
		final File directory=new File(parent, name);
		if(!directory.mkdir()){throw new IllegalArgumentException("Caller run directory must be new: "+directory);}
		return directory;
	}

	private static int count(final boolean[] kept){
		int count=0;
		for(final boolean value:kept){if(value){count++;}}
		return count;
	}

	private static int matched(final int[] matches){
		int count=0;
		for(final int call:matches){if(call>=0){count++;}}
		return count;
	}

	private static boolean[] original(final int size){
		final boolean[] kept=new boolean[size];
		Arrays.fill(kept, true);
		return kept;
	}

	/** Registry identity used in the ordinary caller's generated configuration. */
	private final String family;
	/** Scientific GFF label used to select this member's emitted calls. */
	private final String displayFamily;
	/** Native JVM heap value, converted to MiB only when launching the child. */
	private final String heap;
	/** Fixed scientific options; input, output and process budgets belong to this driver. */
	private final String[] baseArgs;
	/** Ordered genome files shared by baseline and every proposal. */
	private final ArrayList<String> genomePaths;
	/** Stable reference order is the axis of every saved matching array. */
	private final ArrayList<Locus> references;
	/** Threads passed to each sequential caller subprocess. */
	private final int threads;
	/** Initial-protection budget including baseline; regeneration has its own budget. */
	private final int maxCalls;
	/** Human-readable disposition of the latest protect operation. */
	private String resultStatus="Caller protection not completed";
	/** Verified matching belonging to the subset returned by protect. */
	private Result selectedResult;
}
