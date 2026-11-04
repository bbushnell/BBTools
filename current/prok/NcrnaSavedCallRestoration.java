package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;

import fileIO.FileFormat;
import map.ObjectIntMap;
import prok.NcrnaBoundaryComparison.Change;
import prok.NcrnaBoundaryComparison.Result;
import prok.NcrnaLocusComparison.Locus;
import prok.NcrnaSupport.Graph;
import prok.NcrnaSupport.Library;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** One paired-library restoration from supplied caller evidence, without running a caller.
 * Reuses the development restoration rule and the native baseline-priority matcher.
 * Reports describe the supplied GFFs; the exported proposal is not caller-verified.
 * The driver owns parsing and input selection. This helper owns fresh output files
 * and mutates only the candidate's paired-record selection, on the calling thread.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaSavedCallRestoration{

	/*--------------------------------------------------------------*/
	/*----------------          Restoration         ----------------*/
	/*--------------------------------------------------------------*/

	/** Restores each distinct lost baseline winner once and returns the saved-call loss count.
	 * An already-original lost winner cannot make progress; no proposal is exported in that case.
	 */
	static int restore(final NcrnaModelPairs baseline, final NcrnaModelPairs current,
			final Result before, final Result after, final File directory,
			final String modelFilename, final PrintStream diagnostics){
		assert(Arrays.equals(baseline.names, current.names) && before.truth==after.truth)
			: "Restoration must share one ordered model axis and one reference-locus population";
		if(before.matched==0){throw new IllegalArgumentException("Baseline has no detected references; verify family and truth inputs");}
		final ObjectIntMap<String> indices=new ObjectIntMap<String>(String.class);
		int originalBefore=0;
		for(int model=0; model<baseline.names.length; model++){
			indices.put(baseline.names[model], model);
			if(current.same(baseline, model)){originalBefore++;}
		}
		for(Locus call:before.calls){
			if(!indices.contains(call.model)){
				throw new IllegalArgumentException("Baseline winner is absent from the original library: "+call.model);
			}
		}
		if(!directory.mkdir()){throw new IllegalArgumentException("Restoration output must be a fresh directory: "+directory);}
		final Change change=new Change(before, after);
		int restored=0;
		try(Output output=new Output(new File(directory, "lost.tsv").getPath())){
			output.writer.print("reference_index\tseqid\tstart\tstop\tstrand\tbaseline_model\tnewly_restored\n");
			final ByteBuilder row=new ByteBuilder();
			for(int reference=0; reference<before.matches.length; reference++){
				if(before.matches[reference]<0 || after.matches[reference]>=0){continue;}
				final Locus locus=before.truth.get(reference);
				final String winner=before.call(reference).model;
				final int model=indices.get(winner);
				assert(model>=0) : "Every baseline model was checked before paired restoration: "+winner;
				final boolean added=current.restore(baseline, model);
				if(added){restored++;}
				row.clear();
				row.append(reference).tab().append(locus.seqid).tab().append(locus.start).tab().append(locus.stop).tab();
				output.writer.print(row.append(locus.strand==0 ? '+' : '-').tab().append(winner).tab().append(added ? 1 : 0).nl());
			}
		}
		writeSummary(directory, baseline.names.length, before, after, change, originalBefore, restored);
		if(change.lost>0 && restored==0){
			throw new IllegalStateException("Saved-call restoration made no progress: all lost baseline winners already have their original pairs; inspect lost.tsv");
		}
		current.write(directory, modelFilename);
		diagnostics.println("Saved calls lost "+change.lost+" baseline detections; restored "+restored
			+" original pairs. Proposal is not caller-verified; run CallGenes before use.");
		return change.lost;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Optional Support       ----------------*/
	/*--------------------------------------------------------------*/

	/** Measures named corpus retention after zero-loss saved calls, without changing any pair.
	 * Equal supported totals do not establish retention; losses and gains are joined by gene ID.
	 */
	static void measureSupport(final FileFormat genes, final Graph before, final Library library,
			final File directory, final int threads, final int threadsIn, final PrintStream diagnostics){
		assert(Arrays.equals(before.modelNames, library.names))
			: "Final support measurement must retain the restored library's ordered model axis";
		final NcrnaSupport.ScanResult scan=NcrnaSupport.scanWithStats(genes, library, before.cutoff, threads, threadsIn);
		final Graph after=scan.graph;
		if(after.queryNames.length!=before.queryNames.length){
			throw new IllegalArgumentException("Final support scan omitted or added corpus genes");
		}
		final ObjectIntMap<String> indices=new ObjectIntMap<String>(String.class);
		for(int gene=0; gene<after.queryNames.length; gene++){indices.put(after.queryNames[gene], gene);}
		int originalSupported=0, finalSupported=0, lost=0, gained=0;
		try(Output output=new Output(new File(directory, "support.changes.tsv").getPath())){
			output.writer.print("gene\toriginal_supported\tfinal_supported\n");
			final ByteBuilder row=new ByteBuilder();
			for(int gene=0; gene<before.queryNames.length; gene++){
				final int next=indices.get(before.queryNames[gene]);
				if(next<0){throw new IllegalArgumentException("Final support scan has a different gene identity axis");}
				final boolean a=before.starts[gene]<before.starts[gene+1];
				final boolean b=after.starts[next]<after.starts[next+1];
				if(a){originalSupported++;}
				if(b){finalSupported++;}
				if(a && !b){lost++;}
				if(!a && b){gained++;}
				if(a!=b){
					row.clear();
					output.writer.print(row.append(before.queryNames[gene]).tab().append(a ? 1 : 0).tab().append(b ? 1 : 0).nl());
				}
			}
		}
		assert(finalSupported-originalSupported==gained-lost)
			: "Named support gains and losses must reconcile both complete populations";
		NcrnaSupportIO.write(after, new File(directory, "support.tsv").getPath());
		final ByteBuilder row=new ByteBuilder("genes\toriginal_supported\tfinal_supported\tlost\tgained\tinitial_unsupported\talignments\tretention_pass\n");
		row.append(before.queryNames.length).tab().append(originalSupported).tab().append(finalSupported).tab().append(lost).tab().append(gained).tab();
		row.append(before.queryNames.length-originalSupported).tab().append(scan.alignments).tab().append(lost==0 ? 1 : 0).nl();
		try(Output output=new Output(new File(directory, "support.summary.tsv").getPath())){output.writer.print(row);}
		diagnostics.println("Final corpus support: "+lost+" previously supported genes lost; "+gained
			+" newly supported. This measurement does not change the restoration policy or verify the caller.");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Reporting          ----------------*/
	/*--------------------------------------------------------------*/

	/** Grades supplied calls only; the final column explicitly labels the exported proposal. */
	private static void writeSummary(final File directory, final int models, final Result baseline,
			final Result current, final Change change, final int originalBefore, final int restored){
		final ByteBuilder row=new ByteBuilder("models\tbaseline_tp\tcurrent_tp\tunmatched\tfn\tprecision\trecall\tboth_exact\tlost\tgained\toriginal_pairs_before\tnewly_restored\toriginal_pairs_after\tstatus\n");
		row.append(models).tab().append(baseline.matched).tab().append(current.matched).tab().append(current.unmatched).tab();
		row.append(current.missed).tab();
		if(current.calls.isEmpty()){row.append("NA");}
		else{row.append(current.precision(), 8, true);}
		row.tab().append(current.recall(), 8, true).tab().append(current.exactBoth).tab().append(change.lost).tab().append(change.gained).tab();
		row.append(originalBefore).tab().append(restored).tab().append(originalBefore+restored);
		row.append("\tnot caller-verified; metrics describe supplied calls\n");
		try(Output output=new Output(new File(directory, "summary.tsv").getPath())){output.writer.print(row);}
	}

	private NcrnaSavedCallRestoration(){}
}
