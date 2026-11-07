package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectIntMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaSupport.Graph;
import prok.NcrnaSupport.Library;
import prok.NcrnaBoundaryComparison.Result;
import prok.NcrnaLocusComparison.Locus;
import prok.NcrnaTextIO.Output;
import shared.KillSwitch;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;

/** Reduces or regenerates an ncRNA library, with optional caller protection.
 * Forked from template.A_SampleStreamerMT. The shared alignment worker stage is
 * NcrnaSupport; this driver owns CLI parsing, proposal selection and output.
 * Original resources are never changed. A proposal without full-genome calling
 * is explicitly labeled as such and is not a validated caller replacement.
 * @author Brian Bushnell, Raiden
 */
public final class ReduceNcrnaModels{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs reduction and optional caller protection, then closes the diagnostic stream. */
	public static void main(final String[] args){
		final Timer t=new Timer();
		final ReduceNcrnaModels tool=new ReduceNcrnaModels(args);
		try{
			tool.process(t);
		}finally{
			Shared.closeStream(tool.outstream);
		}
	}

	/** Parses native key=value options and validates full-corpus unpaired mode and input/output routes.
	 * Failure after preprocessing closes the selected diagnostic log; successful
	 * construction leaves it open for processing and main's final cleanup.
	 * @param args Command-line arguments, including any config= files
	 */
	public ReduceNcrnaModels(String[] args){
		{
			final PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		// Main cannot own cleanup until this constructor returns successfully.
		try{
			ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
			final Parser parser=parse(args);
			in=parser.in1;
			out=parser.out1;
			threads=Tools.max(1, parser.workers()>0 ? parser.workers() : Shared.threads());
			threadsIn=parser.threadsIn;
			if(parser.append || parser.in2!=null || parser.out2!=null || parser.maxReads>=0 || parser.samplerate!=1){
				throw new IllegalArgumentException("Reduction requires the complete unpaired corpus and a fresh output directory");
			}
			ReadWrite.setZipThreads(threads);
			validateParams();
			in=Tools.fixExtension(in);
			ref=Tools.fixExtension(ref);
			hbm=Tools.fixExtension(hbm);
			// Saved comparisons may legitimately compare a file with itself.
			if(!Tools.testInputFiles(restore || compareOnly, true, in, ref, hbm, genomes, truth, callerConfig,
					graphPath, candidate, candidateHbm, baselineGff, callsGff)){
				throw new IllegalArgumentException("Cannot read a required reducer input");
			}
			ffin=in==null ? null : FileFormat.testInput(in, FileFormat.FASTA, null, true, true);
			if(ffin!=null && !ffin.fasta()){throw new IllegalArgumentException("in= must contain sense-oriented gene FASTA");}
		}catch(RuntimeException | Error failure){
			try{Shared.closeStream(outstream);}catch(RuntimeException | Error cleanup){
				if(cleanup!=failure){failure.addSuppressed(cleanup);}
			}
			throw failure;
		}
	}

	/** Parses reducer-specific options while retaining the standard BBTools option handling. */
	private Parser parse(final String[] args){
		final Parser parser=new Parser();
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			if(scanConflict(a)){incompatibleScanOption=a;}
			if(compareConflict(a)){incompatibleCompareOption=a;}
			if(a.equals("mode")){
				if(!"restore".equalsIgnoreCase(b) && !"scan".equalsIgnoreCase(b) && !"compare".equalsIgnoreCase(b)){
					throw new IllegalArgumentException("The explicit mode must be restore, scan or compare");
				}
				restore="restore".equalsIgnoreCase(b);
				scanOnly="scan".equalsIgnoreCase(b);
				compareOnly="compare".equalsIgnoreCase(b);
			}else if(a.equals("rcomp")){
				rcomp=Parse.parseBoolean(b);
				rcompSet=true;
			}else if(a.equals("reportcuts")){
				final String[] values=NcrnaNumbers.value(a, b).split(",", -1);
				reportCuts=new float[values.length];
				for(int i=0; i<values.length; i++){
					reportCuts[i]=NcrnaNumbers.floatValue(a, values[i]);
					if(!(reportCuts[i]>=0 && reportCuts[i]<=1)){
						throw new IllegalArgumentException("reportcuts requires finite identities in [0,1]");
					}
				}
			}else if(a.equals("candidate")){
				candidate=b;
			}else if(a.equals("candidatehbm")){
				candidateHbm=b;
			}else if(a.equals("baseline")){
				baselineGff=b;
			}else if(a.equals("calls")){
				callsGff=b;
			}else if(a.equals("graph")){
				graphPath=b;
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("hbm")){
				hbm=b;
			}else if(a.equals("cutoff")){
				cutoff=NcrnaNumbers.floatValue(a, b);
				requestedCutoff=NcrnaNumbers.doubleValue(a, b);
				cutoffSet=true;
			}else if(a.equals("target")){
				target=NcrnaNumbers.doubleValue(a, b);
			}else if(a.equals("maxmodelloss")){
				maxModelLoss=NcrnaNumbers.intKMG(a, b);
				if(maxModelLoss<0){throw new IllegalArgumentException("maxmodelloss must be nonnegative when supplied");}
			}else if(a.equals("minmodels")){
				minModels=NcrnaNumbers.intKMG(a, b);
				if(minModels<0){throw new IllegalArgumentException("minmodels must be nonnegative when supplied");}
			}else if(a.equals("genomes")){
				genomes=b;
			}else if(a.equals("truth")){
				truth=b;
			}else if(a.equals("family")){
				family=b;
			}else if(a.equals("truthfamily")){
				truthFamily=b;
			}else if(a.equals("callerconfig")){
				callerConfig=b;
			}else if(a.equals("callerheap")){
				callerHeap=b;
			}else if(a.equals("maxcalls")){
				maxCalls=NcrnaNumbers.intKMG(a, b);
			}else if(a.equals("regenerate")){
				regenerate=Parse.parseBoolean(b);
			}else if(a.equals("reduce")){
				reduce=Parse.parseBoolean(b);
			}else if(a.equals("callerverify")){
				callerVerify=Parse.parseBoolean(b);
			}else if(a.equals("consensusconfig")){
				consensusConfig=b;
			}else if(a.equals("consensusheap")){
				consensusHeap=NcrnaNumbers.binaryValue(a, b)/(1024L*1024);
			}else if(a.equals("maxregencalls")){
				maxRegenCalls=NcrnaNumbers.intKMG(a, b);
			}else if(parser.parse(arg, a, b)){
				// Standard BBTools options.
			}else{
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		return parser;
	}

	/** Identifies explicit non-scan requests; ordinary reducer defaults are not conflicts. */
	private static boolean scanConflict(final String a){
		switch(a){
			case "candidate": case "candidatehbm": case "baseline": case "calls":
			case "graph": case "target": case "maxmodelloss": case "minmodels": case "genomes": case "truth":
			case "family": case "truthfamily": case "callerconfig": case "callerheap":
			case "maxcalls": case "regenerate": case "reduce": case "callerverify":
			case "consensusconfig": case "consensusheap": case "maxregencalls":
				return true;
			default: return false;
		}
	}

	/** Comparison consumes only saved calls and truth; reject explicit work from other modes. */
	private static boolean compareConflict(final String a){
		switch(a){
			case "baseline": case "calls": case "truth": case "family": case "truthfamily":
				return false;
			case "in": case "in1": case "ref": case "hbm": case "cutoff": case "rcomp": case "reportcuts":
				return true;
			default: return scanConflict(a);
		}
	}

	/** Requires a fresh destination and the complete option set when caller protection is requested. */
	private void validateParams(){
		if(compareOnly){
			if(incompatibleCompareOption!=null || in!=null){
				throw new IllegalArgumentException("mode=compare accepts saved calls and truth only; incompatible option: "+incompatibleCompareOption);
			}
			if(baselineGff==null || callsGff==null || truth==null || family==null || truthFamily==null || out==null){
				throw new IllegalArgumentException("mode=compare requires baseline/calls/truth/family/truthfamily/out");
			}
			if(new File(out).exists()){throw new IllegalArgumentException("out= must name a fresh directory: "+out);}
			return;
		}
		if(scanOnly){
			if(incompatibleScanOption!=null){
				throw new IllegalArgumentException("mode=scan does not accept "+incompatibleScanOption+"; no reduction, regeneration or caller runs in this mode");
			}
			if(in==null || ref==null || out==null || !(cutoff>=0 && cutoff<=1)){
				throw new IllegalArgumentException("mode=scan requires in=queries.fa ref=consensus.fa out=freshdir cutoff=[0,1]");
			}
			if(reportCuts==null){reportCuts=new float[]{cutoff};}
			if(new File(out).exists()){throw new IllegalArgumentException("out= must name a fresh directory: "+out);}
			return;
		}
		if(rcompSet || reportCuts!=null){throw new IllegalArgumentException("rcomp and reportcuts require mode=scan");}
		if(maxModelLoss>=0 && (restore || !reduce)){
			throw new IllegalArgumentException("maxmodelloss applies only to model reduction, not restoration or fixed-library regeneration");
		}
		if(minModels>=0 && (restore || !reduce)){
			throw new IllegalArgumentException("minmodels applies only to model reduction, not restoration or fixed-library regeneration");
		}
		if(restore){
			if(ref==null || hbm==null || candidate==null || candidateHbm==null || baselineGff==null
					|| callsGff==null || truth==null || family==null || truthFamily==null || out==null){
				throw new IllegalArgumentException("mode=restore requires ref/hbm/candidate/candidatehbm/baseline/calls/truth/family/truthfamily/out");
			}
			if(callerVerify || genomes!=null || callerConfig!=null || regenerate || !reduce || consensusConfig!=null){
				throw new IllegalArgumentException("mode=restore requires explicit callerverify=f and no live-caller or regeneration request; its proposal is not caller-verified");
			}
			if((in==null)!=(graphPath==null) || (in==null && cutoffSet) || (in!=null && !(cutoff>=0 && cutoff<=1))){
				throw new IllegalArgumentException("Optional restoration support checking requires in=genes.fa graph=original.support.tsv cutoff=[0,1] together");
			}
			if(new File(out).exists()){throw new IllegalArgumentException("out= must name a fresh directory: "+out);}
			return;
		}
		if(candidate!=null || candidateHbm!=null || baselineGff!=null || callsGff!=null){
			throw new IllegalArgumentException("Saved calls and candidate library inputs require mode=restore");
		}
		if(in==null || ref==null || out==null || !(cutoff>=0 && cutoff<=1) || !(target>0 && target<=1)){
			throw new IllegalArgumentException("Require in=genes.fa ref=consensus.fa out=freshdir cutoff=[0,1] target=(0,1]");
		}
		if((genomes==null)!=(truth==null)){
			throw new IllegalArgumentException("genomes= and truth= must be supplied together for caller protection");
		}
		if(!reduce && !regenerate){
			throw new IllegalArgumentException("reduce=f requires regenerate=t to rebuild the supplied fixed library");
		}
		if(!callerVerify && (reduce || !regenerate || genomes!=null || truth!=null || family!=null
				|| truthFamily!=null || callerConfig!=null)){
			throw new IllegalArgumentException("callerverify=f requires reduce=f regenerate=t and no caller-protection inputs; its output is not caller-verified");
		}
		if(genomes!=null && (hbm==null || family==null || truthFamily==null || callerConfig==null || maxCalls<(reduce ? 2 : 1))){
			throw new IllegalArgumentException("Caller protection requires hbm, family, truthfamily, callerconfig and maxcalls>=2 for reduction or >=1 for a fixed library");
		}
		if(regenerate && (hbm==null || consensusConfig==null || consensusHeap<1 || maxRegenCalls<1)){
			throw new IllegalArgumentException("regenerate=t requires hbm, consensusconfig, positive consensusheap and maxregencalls");
		}
		if(regenerate && callerVerify && genomes==null){
			throw new IllegalArgumentException("Regeneration requires caller protection by default; only reduce=f callerverify=f explicitly opts out");
		}
		if(new File(out).exists()){throw new IllegalArgumentException("out= must name a fresh directory: "+out);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	/** Scans once, selects all models or a subset, then optionally protects and regenerates them. */
	private void process(final Timer timer){
		if(compareOnly){
			compareSavedCalls();
			timer.stop();
			outstream.println("Time: "+timer);
			return;
		}
		if(scanOnly){
			scanWindows();
			timer.stop();
			outstream.println("Time: "+timer);
			return;
		}
		if(restore){
			restoreSavedCalls();
			timer.stop();
			outstream.println("Time: "+timer);
			return;
		}
		if(!callerVerify){
			outstream.println("Regeneration proposal only: not caller-verified. Verify the resulting library with CallGenes before use.");
		}
		final NcrnaModelRegeneration regeneration=regenerate ? new NcrnaModelRegeneration(consensusConfig,
			consensusHeap, maxRegenCalls, threads) : null;
		final NcrnaCallerProtection protection=genomes==null ? null : new NcrnaCallerProtection(callerConfig,
			genomes, truth, family, truthFamily, threads, callerHeap, maxCalls);
		final Library library=NcrnaSupport.loadLibrary(ref);
		final NcrnaSupport.ScanResult scan=graphPath==null
			? NcrnaSupport.scanWithStats(ffin, library, cutoff, threads, threadsIn) : null;
		final Graph original=scan==null ? readSupportGraph(library) : scan.graph;
		if(scan==null){outstream.println("Initial support loaded from graph; no fresh alignments performed for the initial scan.");}
		final boolean[] proposal;
		if(reduce){
			proposal=NcrnaModelSubset.reduce(library.names.length, original.starts, original.models, target,
				maxModelLoss<0 ? Integer.MAX_VALUE : maxModelLoss, minModels<0 ? 0 : minModels).kept;
		}else{
			proposal=new boolean[library.names.length];
			Arrays.fill(proposal, true);
		}
		final File directory=new File(out).getAbsoluteFile();
		if(!directory.mkdir()){throw new IllegalArgumentException("Could not create fresh output directory: "+directory);}
		NcrnaSupportIO.write(original, new File(directory, "original.support.tsv").getPath());
		boolean[] kept=proposal;
		String status="Proposal only; caller protection not run";
		if(protection!=null){
			if(reduce){kept=protection.protect(library, ref, hbm, kept, directory);}
			else{protection.verifyOriginal(library, ref, hbm, directory);}
			status=protection.status();
		}
		Graph reduced=original.subset(kept);
		if(regeneration==null){
			writeLibrary(library, hbm, kept, directory);
		}else{
			final NcrnaModelPairs baseline=new NcrnaModelPairs(library, hbm, kept);
			regeneration.process(in, baseline, reduced, protection, new File(directory, "regeneration"), outstream);
			regeneration.pairs.write(directory, modelFilename(hbm));
			reduced=regeneration.graph;
			status=status+"; "+regeneration.status;
		}
		NcrnaSupportIO.write(reduced, new File(directory, "members.tsv").getPath());
		writeSummary(directory, original, reduced, status, scan, regeneration==null ? 0 : regeneration.alignments());
		timer.stop();
		outstream.println("Models retained: "+reduced.modelNames.length+" of "+original.modelNames.length);
		outstream.println("Result: "+status);
		outstream.println("Time: "+timer);
	}

	/** Scores every requested orientation, without proposing or exporting a model library. */
	private void scanWindows(){
		assert(scanOnly && reportCuts!=null) : "Scan validation must establish threshold reporting before any alignment";
		final Library library=NcrnaSupport.loadLibrary(ref);
		if(hbm!=null){
			final boolean[] all=new boolean[library.names.length];
			Arrays.fill(all, true);
			new NcrnaModelPairs(library, hbm, all);
		}
		final NcrnaSupport.ScanResult scan=NcrnaSupport.scanWithStats(ffin, library, cutoff, threads, threadsIn, rcomp);
		final File directory=new File(out).getAbsoluteFile();
		if(!directory.mkdir()){throw new IllegalArgumentException("Could not create fresh output directory: "+directory);}
		NcrnaSupportIO.write(scan.graph, new File(directory, "original.support.tsv").getPath());
		writeScanReports(directory, scan);
		outstream.println("Exhaustive Quantum scan: "+scan.alignments+" pair alignments; identity only, no caller or HBM rescue.");
	}

	/** Maps every oriented query to its physical input and counts each window once per threshold. */
	private void writeScanReports(final File directory, final NcrnaSupport.ScanResult scan){
		final Graph graph=scan.graph;
		final int orientations=rcomp ? 2 : 1;
		if(graph.reduced || graph.queryNames.length%orientations!=0){
			throw new IllegalStateException("Scan reports require a complete ORIGINAL orientation population");
		}
		final int windows=graph.queryNames.length/orientations;
		final long[] passed=new long[reportCuts.length];
		final ByteBuilder bb=new ByteBuilder();
		try(Output output=new Output(new File(directory, "window_queries.tsv").getPath())){
			output.writer.print("window_id\tquery_id\tstrand\tbest_model_index\tbest_identity_a48\n");
			for(int first=0; first<graph.queryNames.length; first+=orientations){
				final int nameLength=windowNameLength(graph, first);
				float best=graph.bestIdentities[first];
				for(int strand=0; strand<orientations; strand++){
					final int query=first+strand;
					best=Math.max(best, graph.bestIdentities[query]);
					bb.clear();
					// Copy only the original identifier prefix into the reused byte buffer.
					for(int i=0; i<nameLength; i++){bb.append(graph.queryNames[first].charAt(i));}
					bb.tab().append(graph.queryNames[query]).tab().append(strand==0 ? '+' : '-').tab();
					bb.append(graph.bestModels[query]).tab().appendFloatA48(graph.bestIdentities[query]).nl();
					output.writer.print(bb);
				}
				for(int i=0; i<reportCuts.length; i++){
					if(best>=reportCuts[i]){passed[i]++;}
				}
			}
		}
		assert(windows>0 && scan.alignments==(long)graph.queryNames.length*graph.modelNames.length)
			: "FPR denominator counts physical inputs; exhaustive work counts every orientation-model pair";
		try(Output output=new Output(new File(directory, "scan.summary.tsv").getPath())){
			output.writer.print("cutoff\tcutoff_a48\twindows\toriented_queries\tmodels\tpassing_windows\tpassing_fraction\tpair_alignments\tactive_workers\tinput_batches\trcomp\n");
			for(int i=0; i<reportCuts.length; i++){
				bb.clear();
				bb.append(reportCuts[i], 9, true).tab().appendFloatA48(reportCuts[i]).tab();
				bb.append(windows).tab().append(graph.queryNames.length).tab().append(graph.modelNames.length).tab();
				bb.append(passed[i]).tab().append(passed[i]/(double)windows, 9, true).tab();
				bb.append(scan.alignments).tab().append(scan.activeWorkers).tab().append(scan.batches).tab().append(rcomp).nl();
				output.writer.print(bb);
			}
		}
	}

	/** Pair adjacency follows NcrnaSupport.ProcessThread; suffixes keep source IDs collision-safe. */
	private int windowNameLength(final Graph graph, final int first){
		final String forward=graph.queryNames[first];
		if(!rcomp){return forward.length();}
		final String reverse=graph.queryNames[first+1];
		final int length=forward.length()-2;
		if(length<1 || !forward.endsWith("/+") || !reverse.endsWith("/-")
				|| reverse.length()!=forward.length() || !forward.regionMatches(0, reverse, 0, length)){
			throw new IllegalStateException("Missing or inconsistent strand pair for physical window: "+forward);
		}
		return length;
	}

	/** Compares saved detections on one reference axis without touching models or invoking callers. */
	private void compareSavedCalls(){
		assert(compareOnly) : "Saved-call comparison dispatch must bypass alignment, restoration and regeneration";
		final ArrayList<Locus> references=NcrnaLocusComparison.readTruth(truth, truthFamily);
		final Result before=NcrnaBoundaryComparison.grade(references, NcrnaLocusComparison.readCalls(baselineGff, family), null);
		final Result after=NcrnaBoundaryComparison.grade(references, NcrnaLocusComparison.readCalls(callsGff, family), before);
		final File directory=new File(out).getAbsoluteFile();
		if(!directory.mkdir()){throw new IllegalArgumentException("Could not create fresh output directory: "+directory);}
		NcrnaBoundaryComparison.writeSavedSummary(new File(directory, "half_shorter.summary.tsv").getPath(), before, after);
		NcrnaBoundaryComparison.writeTransitions(new File(directory, "half_shorter.transitions.tsv").getPath(), before, after);
		outstream.println("Saved calls compared with half-shorter overlap; no caller run or library change. Unmatched calls are not adjudicated false positives.");
	}

	/** Restores paired records from saved evidence without constructing or invoking a caller. */
	private void restoreSavedCalls(){
		final Library original=NcrnaSupport.loadLibrary(ref);
		final Library proposed=NcrnaSupport.loadLibrary(candidate);
		if(!Arrays.equals(original.names, proposed.names)){
			throw new IllegalArgumentException("Saved-call restoration requires identical ordered model identities");
		}
		final boolean[] all=new boolean[original.names.length];
		Arrays.fill(all, true);
		final NcrnaModelPairs baseline=new NcrnaModelPairs(original, hbm, all);
		final NcrnaModelPairs current=new NcrnaModelPairs(proposed, candidateHbm, all);
		final ArrayList<Locus> references=NcrnaLocusComparison.readTruth(truth, truthFamily);
		final Result before=NcrnaBoundaryComparison.grade(references, NcrnaLocusComparison.readCalls(baselineGff, family), null);
		final Result after=NcrnaBoundaryComparison.grade(references, NcrnaLocusComparison.readCalls(callsGff, family), before);
		final File directory=new File(out).getAbsoluteFile();
		outstream.println("Saved-call restoration proposal: not caller-verified. Run CallGenes on the proposal before use.");
		final int lost=NcrnaSavedCallRestoration.restore(baseline, current, before, after, directory,
			modelFilename(candidateHbm), outstream);
		if(lost==0 && in!=null){
			NcrnaSavedCallRestoration.measureSupport(ffin, readSupportGraph(original), current.library(),
				directory, threads, threadsIn, outstream);
		}
	}

	/** Checks named axes without recomputing scores or retaining another full sequence corpus. */
	private Graph readSupportGraph(final Library library){
		final Graph graph=NcrnaSupportIO.read(graphPath);
		if(graph.cutoff!=cutoff || !Arrays.equals(graph.modelNames, library.names)){
			throw new IllegalArgumentException("Saved support graph cutoff or ordered model identities differ from ref=");
		}
		assert(ffin!=null) : "Saved graph use retains the complete in= corpus for its named-axis check";
		final ObjectIntMap<String> indices=new ObjectIntMap<String>(String.class);
		for(int i=0; i<graph.queryNames.length; i++){indices.put(graph.queryNames[i], i);}
		final boolean[] seen=new boolean[graph.queryNames.length];
		final Streamer input=StreamerFactory.makeStreamer(ffin, null, true, -1, false, true, threadsIn);
		int found=0;
		String error=null;
		try{
			input.start();
			for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
				for(Read read:batch){
					try{
						if(read==null || read.bases==null || read.bases.length==0 || read.mate!=null){
							throw new IllegalArgumentException("Saved graph input requires nonempty unpaired genes");
						}
						final String name=NcrnaSupport.recordName(read.id);
						final int index=indices.get(name);
						if(index<0 || seen[index]){
							throw new IllegalArgumentException("Saved graph has a different or repeated corpus gene: "+name);
						}
						seen[index]=true;
						found++;
					}catch(IllegalArgumentException invalid){
						// Drain the native producer before throwing a normal input error;
						// close is not cancellation and must not strand unread batches.
						if(error==null){error=invalid.getMessage();}
					}
				}
			}
		}catch(Throwable failure){KillSwitch.exceptionKill(failure);}
		if(ReadWrite.closeStream(input)){throw new IllegalStateException("Saved graph corpus read failed: "+in);}
		if(error!=null){throw new IllegalArgumentException(error);}
		if(found!=seen.length){throw new IllegalArgumentException("Saved graph contains genes absent from the complete input corpus");}
		return graph;
	}

	/** Exports a paired subset with original headers and HBM count data, and upper-case sequences.
	 * Each HBM is checked against the same-position consensus before selection.
	 * This is a library subset, not an installable family resource directory.
	 */
	static void writeLibrary(final Library library, final String hbm,
			final boolean[] kept, final File directory){
		new NcrnaModelPairs(library, hbm, kept).write(directory, hbm==null ? null : modelFilename(hbm));
	}

	/** Preserves the source HBM encoding and gzip suffix in the output filename. */
	static String modelFilename(final String hbm){
		final String suffix=NucleotideHbmCodec.isCompact(hbm) ? ".hbmc" : ".hbm";
		return "models"+suffix+(hbm.endsWith(".gz") ? ".gz" : "");
	}

	/** Reports full-corpus coverage and measured work without implying genomic false-positive validation. */
	private void writeSummary(final File directory, final Graph original, final Graph reduced,
			final String status, final NcrnaSupport.ScanResult scan, final long regenerationAlignments){
		final ByteBuilder bb=new ByteBuilder();
		bb.append("status\t").append(status).nl();
		bb.append("models_original\t").append(original.modelNames.length).nl();
		bb.append("models_retained\t").append(reduced.modelNames.length).nl();
		bb.append("genes\t").append(original.queryNames.length).nl();
		bb.append("coverage_original\t").append(coverage(original), 8, true).nl();
		bb.append("coverage_retained\t").append(coverage(reduced), 8, true).nl();
		bb.append("identity_cutoff\t").append(requestedCutoff, 8, true).nl();
		bb.append("coverage_target\t");
		if(reduce){bb.append(target, 8, true);}
		else{bb.append("NA");}
		bb.nl();
		if(maxModelLoss>=0){
			bb.append("max_model_loss\t").append(maxModelLoss).nl();
			bb.append("exported_support_queries_lost\t").append(writeLostQueries(directory, original, reduced)).nl();
		}
		if(minModels>=0){bb.append("min_models\t").append(minModels).nl();}
		final long initialAlignments=scan==null ? 0 : scan.alignments;
		bb.append("pair_alignments\t").append(initialAlignments+regenerationAlignments).nl();
		bb.append("initial_pair_alignments\t").append(initialAlignments).nl();
		bb.append("regeneration_pair_alignments\t").append(regenerationAlignments).nl();
		bb.append("input_batches\t").append(scan==null ? 0 : scan.batches).nl();
		bb.append("active_workers\t").append(scan==null ? 0 : scan.activeWorkers).nl();
		if(scan==null){bb.append("initial_support\tSaved graph; no fresh alignments performed for the initial scan.\n");}
		bb.append("false_positive_rate\tUNMEASURED\n");
		bb.append(reduce ? "method\tGreedy subset proposal; no global minimum claim.\n"
			: "method\tFixed-library regeneration; no model reduction.\n");
		if(genomes==null){bb.append("limitation\tNo full-genome caller protection was performed.\n");}
		try(Output output=new Output(new File(directory, "summary.txt").getPath())){
			output.writer.print(bb);
		}
	}

	/** Lists input-supported genes with no supporter in the exported graph.
	 * Caller restoration and optional regeneration may alter the original proposal.
	 * These are support losses, not claims about whole-genome caller detections.
	 */
	private static int writeLostQueries(final File directory, final Graph original, final Graph reduced){
		assert(Arrays.equals(original.queryNames, reduced.queryNames))
			: "Support-loss reporting compares the same complete gene axis before and after export";
		int lost=0;
		final ByteBuilder row=new ByteBuilder();
		try(Output output=new Output(new File(directory, "lost_queries.tsv").getPath())){
			output.writer.print("query_id\n");
			for(int gene=0; gene<original.queryNames.length; gene++){
				if(original.starts[gene+1]>original.starts[gene] && reduced.starts[gene+1]==reduced.starts[gene]){
					output.writer.print(row.clear().append(original.queryNames[gene]).nl());
					lost++;
				}
			}
		}
		return lost;
	}

	/** Fraction of all input genes having at least one retained model meeting the alignment cutoff. */
	private static double coverage(final Graph graph){
		assert(graph.queryNames.length>0) : "Coverage denominator includes the complete nonempty query population; genes="+graph.queryNames.length;
		int covered=0;
		for(int gene=0; gene<graph.queryNames.length; gene++){
			if(graph.starts[gene+1]>graph.starts[gene]){covered++;}
		}
		return covered/(double)graph.queryNames.length;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Gene FASTA, consensus FASTA, optional paired HBM, and fresh output directory. */
	private String in, ref, hbm, out;
	/** Optional scored graph replacing only the initial exhaustive alignment scan. */
	private String graphPath;
	/** Explicit saved-call mode and its candidate paired library plus two saved GFFs. */
	private boolean restore=false;
	/** Read-only saved-GFF comparison using the existing half-shorter matcher. */
	private boolean compareOnly=false;
	private String incompatibleCompareOption;
	/** Pure exhaustive window assay; ordinary reducer/regenerator behavior is unchanged. */
	private boolean scanOnly=false;
	private boolean rcomp=true;
	private boolean rcompSet=false;
	private float[] reportCuts;
	/** Explicit irrelevant option, independent of default values in other modes. */
	private String incompatibleScanOption;
	private String candidate, candidateHbm, baselineGff, callsGff;
	/** Optional protection inputs: genome list, truth table, caller/truth labels and scientific config. */
	private String genomes, truth, family, truthFamily, callerConfig;
	/** Native JVM heap value for each sequential caller child, such as 2g. */
	private String callerHeap="2g";
	/** Minimum native float alignment identity for a gene-to-model support edge. */
	private float cutoff=Float.NaN;
	/** Original double value for human reports; alignments retain native float semantics. */
	private double requestedCutoff=Double.NaN;
	/** Distinguishes absent optional restoration support settings from an explicit cutoff. */
	private boolean cutoffSet=false;
	/** Required supported fraction of ALL input genes, including initially unsupported genes. */
	private double target=0.98;
	/** Unset (-1) preserves the legacy total-target-only policy; explicit zero forbids per-drop loss. */
	private int maxModelLoss=-1;
	/** Unset (-1) preserves legacy reduction; an explicit floor may stop zero-loss drops early. */
	private int minModels=-1;
	/** Alignment-worker budget shared by support and caller stages. */
	private int threads=1;
	/** Native input-reader worker setting, separate from alignment workers. */
	private int threadsIn=1;
	/** Initial caller-protection budget, including its baseline call. */
	private int maxCalls=4;
	/** Whether to rebuild retained models after the initial caller protection passes. */
	private boolean regenerate=false;
	/** Whether to select a smaller library before optional protection/regeneration. */
	private boolean reduce=true;
	/** Regeneration verifies caller retention unless fixed-library mode explicitly opts out. */
	private boolean callerVerify=true;
	/** Scientific consensus-builder recipe used only for regeneration. */
	private String consensusConfig;
	/** Consensus-builder child heap budget in MiB. */
	private long consensusHeap=1024;
	/** Additional regeneration caller budget; the protected baseline is reused. */
	private int maxRegenCalls=5;
	/** Validated complete unpaired gene FASTA input descriptor. */
	private FileFormat ffin;
	/** Selected diagnostic stream; constructor failure or main closes it. */
	private PrintStream outstream=System.err;
}
