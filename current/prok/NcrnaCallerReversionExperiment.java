package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.FileFormat;
import map.ObjectIntMap;
import parse.Parse;
import parse.Parser;
import prok.NcrnaBoundaryComparison.Change;
import prok.NcrnaBoundaryComparison.Result;
import prok.NcrnaLocusComparison.Locus;
import prok.NcrnaSupport.Library;
import prok.NcrnaSupport.Graph;
import prok.NcrnaTextIO.Output;
import prok.NucleotideHbmCodec.Model;
import stream.Read;
import structures.ByteBuilder;

/** One development-only caller-protection step for a fixed model identity axis.
 * Uses saved ordinary GFFs and the existing baseline-priority matcher. Lost
 * references restore their baseline winning consensus/HBM pair from the original
 * library. Other pairs remain current. The shell recipe runs the next ordinary
 * caller pass; this step never predicts that the proposed changes recover a gene.
 * Preserved standalone experiment: created to restore pairs using saved GFFs
 * without launching another caller. The new reducer covers the normal workflow
 * but does not expose this caller-free mode. Retained pending disposition for
 * future saved-call restoration; no new accuracy or cost claim is implied.
 * Native a/b parsing and processing follow template.A_SampleBasic. The enclosing
 * experiment entry point owns PreParser and its diagnostic stream.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCallerReversionExperiment{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Called only after the shared development entry point expands native config files. */
	static void run(final String[] args, final PrintStream diagnostics){
		new NcrnaCallerReversionExperiment(args).process(diagnostics);
	}

	/** Reads a complete saved-call comparison and two index-aligned paired libraries. */
	private NcrnaCallerReversionExperiment(final String[] args){
		final Parser parser=new Parser();
		for(String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			if(a.equals("mode")){
				if(!"protect".equalsIgnoreCase(b)){throw new IllegalArgumentException("This development step requires mode=protect");}
			}else if(a.equals("installation")){
				// Shared launcher metadata; this step consumes saved calls and launches no child.
			}else if(a.equals("ref")){
				originalFasta=b;
			}else if(a.equals("hbm")){
				currentHbm=b;
			}else if(a.equals("originalhbm")){
				originalHbm=b;
			}else if(a.equals("baseline")){
				baselineGff=b;
			}else if(a.equals("calls")){
				currentGff=b;
			}else if(a.equals("truth")){
				truthPath=b;
			}else if(a.equals("family")){
				family=b;
			}else if(a.equals("truthfamily")){
				truthFamily=b;
			}else if(a.equals("expectedmodels")){
				expectedModels=Parse.parseIntKMG(b);
			}else if(a.equals("genes")){
				genes=b;
			}else if(a.equals("graph")){
				originalGraph=b;
			}else if(a.equals("expectedgenes")){
				expectedGenes=Parse.parseIntKMG(b);
			}else if(a.equals("cutoff")){
				cutoff=Float.parseFloat(b);
			}else if(!parser.parse(arg, a, b)){
				throw new IllegalArgumentException("Unknown reversion parameter: "+arg);
			}
		}
		currentFasta=parser.in1;
		threads=parser.workers()>0 ? parser.workers() : shared.Shared.threads();
		if(genes==null ? originalGraph!=null || expectedGenes!=-1 || !Float.isNaN(cutoff)
				: originalGraph==null || expectedGenes<1 || !(cutoff>=0 && cutoff<=1) || threads<1){
			throw new IllegalArgumentException("Final support measurement requires genes/graph/expectedgenes/cutoff together and positive threads");
		}
		if(currentFasta==null || originalFasta==null || currentHbm==null || originalHbm==null
				|| baselineGff==null || currentGff==null || truthPath==null || family==null
				|| truthFamily==null || parser.out1==null || expectedModels<1){
			throw new IllegalArgumentException("Supply in/ref/hbm/originalhbm/baseline/calls/truth/family/truthfamily/out and positive expectedmodels");
		}
		if(parser.in2!=null || parser.out2!=null || parser.maxReads>=0 || parser.samplerate!=1 || parser.append){
			throw new IllegalArgumentException("Caller protection uses complete saved inputs and a fresh output directory");
		}
		root=new File(parser.out1).getAbsoluteFile();
		if(root.exists() || !root.getParentFile().isDirectory()){
			throw new IllegalArgumentException("out= must be a fresh directory under an existing parent");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Reversion          ----------------*/
	/*--------------------------------------------------------------*/

	/** Restores each distinct lost baseline winner once, retaining all other pairs. */
	private void process(final PrintStream diagnostics){
		final Library original=NcrnaSupport.loadLibrary(originalFasta);
		final Library current=NcrnaSupport.loadLibrary(currentFasta);
		if(original.names.length!=expectedModels || !Arrays.equals(original.names, current.names)){
			throw new IllegalArgumentException("Caller reversion requires the same fixed ordered model identities in both libraries");
		}
		final Model[] originalModels=NcrnaSequentialRegenerationExperiment.loadModels(original, originalHbm);
		final Model[] currentModels=NcrnaSequentialRegenerationExperiment.loadModels(current, currentHbm);
		final ObjectIntMap<String> index=new ObjectIntMap<String>(String.class);
		final boolean[] restore=new boolean[expectedModels];
		int originalBefore=0;
		for(int model=0; model<expectedModels; model++){
			index.put(original.names[model], model);
			restore[model]=samePair(original, originalModels, current, currentModels, model);
			if(restore[model]){originalBefore++;}
		}
		final ArrayList<Locus> truth=NcrnaLocusComparison.readTruth(truthPath, truthFamily);
		final Result baseline=NcrnaBoundaryComparison.grade(truth, NcrnaLocusComparison.readCalls(baselineGff, family), null);
		final Result candidate=NcrnaBoundaryComparison.grade(truth, NcrnaLocusComparison.readCalls(currentGff, family), baseline);
		if(baseline.matched==0){throw new IllegalArgumentException("Baseline has no detected references; verify family and truth inputs");}
		for(Locus call:baseline.calls){
			if(!index.contains(call.model)){throw new IllegalArgumentException("Baseline winner is absent from the original library: "+call.model);}
		}
		if(!root.mkdir()){throw new IllegalStateException("Cannot create fresh reversion directory");}
		final Change change=new Change(baseline, candidate);
		int restored=0;
		try(Output output=new Output(path("lost.tsv"))){
			output.writer.print("reference_index\tseqid\tstart\tstop\tstrand\tbaseline_model\tnewly_restored\n");
			final ByteBuilder row=new ByteBuilder();
			for(int reference=0; reference<truth.size(); reference++){
				if(baseline.matches[reference]<0 || candidate.matches[reference]>=0){continue;}
				final Locus locus=truth.get(reference);
				final String winner=baseline.call(reference).model;
				final int model=index.get(winner);
				assert(model>=0) : "Every baseline winner was checked against the fixed library before restoration";
				final boolean added=!restore[model];
				if(added){
					restore[model]=true;
					restored++;
				}
				row.clear();
				row.append(reference).tab().append(locus.seqid).tab().append(locus.start).tab().append(locus.stop).tab();
				row.append(locus.strand==0 ? '+' : '-').tab().append(winner).tab().append(added ? 1 : 0).nl();
				output.writer.print(row);
			}
		}
		writeSummary(baseline, candidate, change, originalBefore, restored);
		if(change.lost>0 && restored==0){
			throw new IllegalStateException("Caller protection made no progress: all lost baseline winners already have their original pairs; inspect lost.tsv");
		}
		final ArrayList<Read> records=new ArrayList<Read>(expectedModels);
		try(NucleotideHbmCodec.Writer writer=new NucleotideHbmCodec.Writer(path("models.hbmc"), true, false)){
			for(int model=0; model<expectedModels; model++){
				final Library source=restore[model] ? original : current;
				final Model pair=restore[model] ? originalModels[model] : currentModels[model];
				assert(Arrays.equals(source.bases[model], pair.bases)) : "Restoration must keep consensus and HBM in the same coordinate frame";
				records.add(new Read(source.bases[model], null, source.headers[model], model));
				writer.add(pair);
			}
			writer.finish();
		}
		NcrnaSequentialRegenerationExperiment.writeFasta(path("consensus.fa"), records);
		if(change.lost==0 && genes!=null){measureSupport(current, diagnostics);}
		diagnostics.println("Caller reversion step completed: "+change.lost+" baseline detections lost in supplied calls; "+restored+" original pairs restored.");
	}

	/** Measures final corpus retention once after caller convergence, without altering the proposal.
	 * Equal supported totals do not imply retention: losses and gains are joined by gene ID.
	 */
	private void measureSupport(final Library library, final PrintStream diagnostics){
		final Graph before=NcrnaSupportIO.read(originalGraph);
		if(before.queryNames.length!=expectedGenes || before.cutoff!=cutoff || !Arrays.equals(before.modelNames, library.names)){
			throw new IllegalArgumentException("Original support graph must match the declared complete gene/model axes and cutoff");
		}
		final FileFormat fasta=FileFormat.testInput(genes, FileFormat.FASTA, null, true, true);
		final NcrnaSupport.ScanResult scan=NcrnaSupport.scanWithStats(fasta, library, cutoff, threads, 1);
		final Graph after=scan.graph;
		if(after.queryNames.length!=expectedGenes){throw new IllegalArgumentException("Final support scan omitted or added corpus genes");}
		final ObjectIntMap<String> indices=new ObjectIntMap<String>(String.class);
		for(int gene=0; gene<after.queryNames.length; gene++){indices.put(after.queryNames[gene], gene);}
		int originalSupported=0, finalSupported=0, lost=0, gained=0;
		try(Output output=new Output(path("support.changes.tsv"))){
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
		assert(finalSupported-originalSupported==gained-lost) : "Named support gains and losses must reconcile both complete populations";
		NcrnaSupportIO.write(after, path("support.tsv"));
		final ByteBuilder row=new ByteBuilder("genes\toriginal_supported\tfinal_supported\tlost\tgained\tinitial_unsupported\talignments\tretention_pass\n");
		row.append(expectedGenes).tab().append(originalSupported).tab().append(finalSupported).tab().append(lost).tab().append(gained).tab();
		row.append(expectedGenes-originalSupported).tab().append(scan.alignments).tab().append(lost==0 ? 1 : 0).nl();
		try(Output output=new Output(path("support.summary.tsv"))){output.writer.print(row);}
		diagnostics.println("Final corpus support: "+lost+" previously supported genes lost; "+gained+" newly supported. This measurement does not change the reversion policy.");
	}

	/** Whole paired-record equality distinguishes already-original winners from useful restorations. */
	private static boolean samePair(final Library a, final Model[] am, final Library b, final Model[] bm, final int model){
		return a.headers[model].equals(b.headers[model]) && Arrays.equals(am[model].header, bm[model].header)
			&& Arrays.equals(am[model].bases, bm[model].bases) && Arrays.equals(am[model].counts, bm[model].counts);
	}

	/** Counts grade the supplied calls, never the untested proposal emitted by this step. */
	private void writeSummary(final Result baseline, final Result current, final Change change,
			final int originalBefore, final int restored){
		final ByteBuilder row=new ByteBuilder("models\tbaseline_tp\tcurrent_tp\tunmatched\tfn\tprecision\trecall\tboth_exact\tlost\tgained\toriginal_pairs_before\tnewly_restored\toriginal_pairs_after\n");
		row.append(expectedModels).tab().append(baseline.matched).tab().append(current.matched).tab().append(current.unmatched).tab();
		row.append(current.missed).tab();
		if(current.calls.isEmpty()){row.append("NA");}
		else{row.append(current.precision(), 8);}
		row.tab().append(current.recall(), 8).tab().append(current.exactBoth).tab().append(change.lost).tab().append(change.gained).tab();
		row.append(originalBefore).tab().append(restored).tab().append(originalBefore+restored).nl();
		try(Output output=new Output(path("summary.tsv"))){output.writer.print(row);}
	}

	private String path(final String name){return new File(root, name).getPath();}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String currentFasta, originalFasta, currentHbm, originalHbm;
	private String baselineGff, currentGff, truthPath, family, truthFamily;
	private int expectedModels=-1;
	private String genes, originalGraph;
	private int expectedGenes=-1, threads;
	private float cutoff=Float.NaN;
	private final File root;
}
