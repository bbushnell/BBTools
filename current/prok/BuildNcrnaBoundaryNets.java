package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.Locale;
import java.util.concurrent.ArrayBlockingQueue;

import dna.Data;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaBoundaryComparison.Change;
import prok.NcrnaBoundaryComparison.Result;
import prok.NcrnaFamilyExport.Boundary;
import prok.NcrnaNativeTrainer.Method;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Learns, compares and refits ordinary shared-end ncRNA boundary networks.
 * CLI/lifecycle adapted from template.A_SampleByteFileMT; streaming work lives
 * in the Streamer-template learning components. The main thread owns stages,
 * immutable populations and reports. Only independent caller trials overlap.
 * out= is the final loadable family directory. Intermediate learning products,
 * caller outputs and reports remain in the separate out.work directory.
 * Neither directory may exist on entry. Installed resources are never changed.
 * @author Brian Bushnell, Raiden
 */
public final class BuildNcrnaBoundaryNets{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs learning and evaluation, closing the diagnostic stream when processing finishes or fails. */
	public static void main(final String[] args){
		final Timer timer=new Timer();
		final BuildNcrnaBoundaryNets tool=new BuildNcrnaBoundaryNets(args);
		try{
			tool.process(timer);
		}finally{
			Shared.closeStream(tool.outstream);
		}
	}

	/** Parses native key=value options and validates routes and geometry before creating outputs.
	 * Failure after preprocessing closes the selected diagnostic log; successful
	 * construction leaves it open for processing and main's final cleanup.
	 * @param args Command-line arguments, including any config= files
	 */
	public BuildNcrnaBoundaryNets(String[] args){
		{
			final PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		// Main cannot own cleanup until this constructor returns successfully.
		try{
			final Parser parser=parse(args);
			in=parser.in1;
			out=parser.out1;
			threads=Tools.max(1, parser.workers()>0 ? parser.workers() : Shared.threads());
			if(parser.append || parser.in2!=null || parser.out2!=null || parser.maxReads>=0 || parser.samplerate!=1){
				throw new IllegalArgumentException("Boundary learning requires the complete unpaired context corpus and fresh output directories");
			}
			ReadWrite.setZipThreads(threads);
			validateParams();
		}catch(RuntimeException | Error failure){
			try{Shared.closeStream(outstream);}catch(RuntimeException | Error cleanup){
				if(cleanup!=failure){failure.addSuppressed(cleanup);}
			}
			throw failure;
		}
	}

	/** Standard a/b parsing; scientific choices remain explicit rather than hidden in a recipe wrapper. */
	private Parser parse(final String[] args){
		final Parser parser=new Parser();
		for(String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			if(a.equals("roster")){
				roster=b;
			}else if(a.equals("loci")){
				loci=b;
			}else if(a.equals("exclude")){
				exclude=b;
			}else if(a.equals("graph")){
				graph=b;
			}else if(a.equals("original")){
				original=b;
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("hbm")){
				hbm=b;
			}else if(a.equals("seeds")){
				seeds=b;
			}else if(a.equals("family")){
				family=b;
			}else if(a.equals("truthfamily")){
				truthFamily=b;
			}else if(a.equals("genomes")){
				genomes=b;
			}else if(a.equals("truth")){
				truth=b;
			}else if(a.equals("callerconfig")){
				callerConfig=b;
			}else if(a.equals("classifierconfig")){
				classifierConfig=b;
			}else if(a.equals("regressionconfig")){
				regressionConfig=b;
			}else if(a.equals("splitseed")){
				splitSeed=NcrnaNumbers.longValue(a, b);
			}else if(a.equals("shards")){
				shards=NcrnaNumbers.intKMG(a, b);
			}else if(a.equals("features")){
				features=b;
			}else if(a.equals("startradius")){
				startRadius=NcrnaNumbers.intValue(a, b);
			}else if(a.equals("stopradius")){
				stopRadius=NcrnaNumbers.intValue(a, b);
			}else if(a.equals("startinside")){
				startInside=NcrnaNumbers.intValue(a, b);
			}else if(a.equals("startoutside")){
				startOutside=NcrnaNumbers.intValue(a, b);
			}else if(a.equals("stopinside")){
				stopInside=NcrnaNumbers.intValue(a, b);
			}else if(a.equals("stopoutside")){
				stopOutside=NcrnaNumbers.intValue(a, b);
			}else if(a.equals("marginstart")){
				startMargin=NcrnaNumbers.floatValue(a, b);
			}else if(a.equals("marginstop")){
				stopMargin=NcrnaNumbers.floatValue(a, b);
			}else if(a.equals("cutoff")){
				cutoff=NcrnaNumbers.floatValue(a, b);
			}else if(a.equals("maxlost")){
				maxLost=NcrnaNumbers.intKMG(a, b);
			}else if(a.equals("minprecision")){
				minPrecision=NcrnaNumbers.doubleValue(a, b);
			}else if(a.equals("tiepreference")){
				try{
					tiePreference=Method.valueOf(NcrnaNumbers.value(a, b).toUpperCase(Locale.ROOT));
				}catch(IllegalArgumentException failure){
					throw new IllegalArgumentException("Invalid "+a+"="+b+"; use classifier or regression", failure);
				}
			}else if(parser.parse(arg, a, b)){
				// Standard BBTools options.
			}else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		return parser;
	}

	/** Rejects malformed routes and geometry before creating files or generating vectors. */
	private void validateParams(){
		final String[] required={in, roster, loci, graph, original, ref, hbm, seeds, family,
			truthFamily, genomes, truth, callerConfig, classifierConfig, regressionConfig, out};
		for(String value:required){
			if(value==null || value.isEmpty()){throw new IllegalArgumentException("Supply all required learning, library, evaluation, recipe and output arguments; see launcher help");}
		}
		family=CallGenes.parseNcrnaFamily(family);
		NcrnaFamilyRegistry.sharedBoundaryMember(family);
		if(!"v1".equals(features) || splitSeed<0 || shards<1 || startRadius<0 || stopRadius<0
				|| maxLost<0 || !Double.isFinite(minPrecision) || minPrecision<.95 || minPrecision>1){
			throw new IllegalArgumentException("Require features=v1, splitseed>=0, shards>=1, nonnegative radii/maxlost, minprecision in [0.95,1]");
		}
		geometry(startInside, startOutside);
		geometry(stopInside, stopOutside);
		if(startRadius==0 && stopRadius==0){
			throw new IllegalArgumentException("startradius and stopradius cannot both be zero: training requires negative examples");
		}
		boundary(null, null);
		if(!Tools.testInputFiles(true, false, in, roster, loci, graph, original, ref, hbm, seeds,
				exclude, genomes, truth, callerConfig, classifierConfig, regressionConfig)){
			throw new IllegalArgumentException("Cannot read a required boundary-builder input");
		}
		NcrnaCallerProcess.configurationPath(new File(callerConfig));
		NcrnaCallerProcess.configurationPath(new File(classifierConfig));
		NcrnaCallerProcess.configurationPath(new File(regressionConfig));
		destination=new File(out).getAbsoluteFile();
		work=new File(destination.getPath()+".work");
		final String path=destination.getPath();
		if(!path.equals(path.trim()) || path.indexOf(',')>=0 || path.indexOf('=')>=0
				|| path.indexOf('\n')>=0 || path.indexOf('\r')>=0){
			throw new IllegalArgumentException("out= must be representable in native config and shard-list filenames: "+path);
		}
		if(destination.exists() || work.exists() || !destination.getParentFile().isDirectory()){
			throw new IllegalArgumentException("out= and out.work must be fresh directories beneath an existing parent");
		}
	}

	/** Requires nonnegative endpoint-window base counts whose sum is a supported k-mer length. */
	private static void geometry(final int inside, final int outside){
		if(inside<0 || outside<0 || (long)inside+outside<1 || (long)inside+outside>15){
			throw new IllegalArgumentException("Each endpoint needs nonnegative inside/outside counts summing to 1..15");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	/** Keeps the original coverage axis separate from the accepted model/feature axis. */
	private void process(final Timer timer){
		final NcrnaCallerEvaluation evaluation=new NcrnaCallerEvaluation(callerConfig, genomes, truth, family, truthFamily);
		requireRegisteredMapping(Parser.parseConfig(new String[]{"config="+callerConfig}));
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.resources.consensus=ref;
		options.resources.models=hbm;
		options.resources.seeds=seeds;
		options.boundaryOverride=false;
		options.endpoints=false;
		// Ordinary CallGenes explicitly chooses registry topN, even for dynamic libraries.
		options.topN=NcrnaFamilyRegistry.get(family).topN;
		final NcrnaFamily accepted=NcrnaFamilyLoader.load(family, options, new FamilyResourceResolver());
		final int median=PackagedBoundaryConfig.upperMedianLength(accepted.library);
		if(startRadius>median || stopRadius>median){throw new IllegalArgumentException("Radii exceed accepted-library upper-median length: "+median);}
		final NcrnaLearningPopulation population=NcrnaLearningPopulation.load(NcrnaSupportIO.read(graph),
			NcrnaSupport.loadLibrary(original).names, roster, loci, exclude, splitSeed);
		final String launcher=Data.findPath("?train.sh", false);
		if(launcher==null){throw new IllegalStateException("Cannot locate the native train.sh launcher");}
		final NcrnaNativeTrainer trainer=new NcrnaNativeTrainer(new File(launcher).getParentFile(), threads);
		if(!work.mkdir()){throw new IllegalArgumentException("Could not create fresh work directory: "+work);}
		writePopulation(population);
		final File tuning=fresh("tuning");
		final NcrnaLearningVectors vectors=learn(population, accepted, tuning, false);
		final File classifier=new File(tuning, "classifier.bbnet");
		final File regression=new File(tuning, "regression.bbnet");
		trainer.train(Method.CLASSIFIER, classifierConfig, vectors.trainPaths, vectors.validationPath, classifier);
		trainer.train(Method.REGRESSION, regressionConfig, vectors.trainPaths, vectors.validationPath, regression);
		final Boundary off=boundary(null, null);
		final Boundary c=boundary(classifier, tuning);
		final Boundary r=boundary(regression, tuning);
		final File[] roots={new File(work, "off_family"), new File(work, "classifier_family"), new File(work, "regression_family")};
		export(accepted, off, roots[0]);
		export(accepted, c, roots[1]);
		export(accepted, r, roots[2]);
		final File[] calls=callTrials(evaluation, roots);
		final Result baseline=evaluation.grade(calls[0], null);
		final Result classifierResult=evaluation.grade(calls[1], baseline);
		final Result regressionResult=evaluation.grade(calls[2], baseline);
		NcrnaBoundaryComparison.writeSummary(new File(work, "comparison.tsv").getPath(), baseline, classifierResult, regressionResult);
		NcrnaBoundaryComparison.writeTransitions(new File(work, "classifier_transitions.tsv").getPath(), baseline, classifierResult);
		NcrnaBoundaryComparison.writeTransitions(new File(work, "regression_transitions.tsv").getPath(), baseline, regressionResult);
		final Method selected=NcrnaBoundaryComparison.choose(baseline, classifierResult, regressionResult, maxLost, minPrecision, tiePreference);
		Boundary product=off;
		String status="OFF";
		Result finalResult=baseline;
		if(selected!=null){
			final File refit=fresh("refit");
			final NcrnaLearningVectors all=learn(population, accepted, refit, true);
			assert(all.trainPaths.length==1 && all.validationPath==null)
				: "Final refit includes every eligible locus and monitors that same single input, never a holdout";
			final File net=new File(refit, "boundary.bbnet");
			trainer.train(selected, selected==Method.CLASSIFIER ? classifierConfig : regressionConfig,
				all.trainPaths, all.trainPaths[0], net);
			final Boundary refitted=boundary(net, refit);
			final File root=new File(refit, "family");
			export(accepted, refitted, root);
			final File gff=evaluation.run(root, new File(refit, "caller"), threads, childHeap());
			final Result result=evaluation.grade(gff, baseline);
			NcrnaBoundaryComparison.writeTransitions(new File(work, "refit_transitions.tsv").getPath(), baseline, result);
			writeRefit(baseline, result);
			if(new Change(baseline, result).eligible(maxLost, minPrecision)){
				product=refitted;
				finalResult=result;
				status="ENABLED_"+selected.name();
			}else{status="OFF_REFIT_DID_NOT_PASS";}
		}
		export(accepted, product, destination);
		writeDecision(status, selected, baseline, finalResult);
		timer.stop();
		outstream.println("Boundary result: "+status);
		outstream.println("Family directory: "+destination);
		outstream.println("Learning and evaluation records: "+work);
		outstream.println("Time: "+timer);
	}

	/** Rejects selected-family mapping overrides to keep learning and caller shortlists consistent.
	 * These generic flags target the selected family. Historical r58/lsu and bare
	 * tRNA index flags affect other callers and remain ordinary evaluation policy.
	 * Supporting selected overrides later requires a shared public options parser.
	 * @param arguments Configuration-expanded native arguments
	 * @throws IllegalArgumentException If an argument overrides the registered mapping policy
	 */
	public static void requireRegisteredMapping(final String[] arguments){
		for(String arg:arguments){
			final String a=Parse.splitOnFirst(arg, '=')[0].toLowerCase(Locale.ROOT);
			switch(a){
				case "ncrnaindexk": case "ncrnamapk": case "ncrnaindextopn":
				case "ncrnaindexminhits": case "ncrnafixedminhits":
				case "ncrnaadaptiveminhits": case "ncrnaadaptiveshortlist":
				case "ncrnaadaptfloor": case "ncrnaadapttopfrac": case "ncrnaadaptqfrac":
				case "ncrnastrictindexcutoff": case "ncrnarankedfallback":
					throw new IllegalArgumentException("Boundary learning uses registered family mapping settings; remove selected-family override: "+a);
				default: break;
			}
		}
	}

	/** Two complete input passes: count true endpoint words, then emit candidate features. */
	private NcrnaLearningVectors learn(final NcrnaLearningPopulation population, final NcrnaFamily accepted,
			final File directory, final boolean refit){
		final NcrnaLearningTables tables=NcrnaLearningTables.build(population, in,
			startInside, startOutside, stopInside, stopOutside, refit, threads);
		tables.start.write(new File(directory, "start.tsv").getPath());
		tables.stop.write(new File(directory, "stop.tsv").getPath());
		final NcrnaLearningVectors vectors=NcrnaLearningVectors.write(population, in, accepted, tables,
			startRadius, stopRadius, refit ? 1 : shards, threads, new File(directory, "vectors").getPath());
		final ByteBuilder report=new ByteBuilder("records\tfiltered\tunavailable\tno_model\tprocessed\tcandidates\n");
		report.append(vectors.records).tab().append(vectors.filtered).tab().append(vectors.unavailable).tab();
		report.append(vectors.noModel).tab().append(vectors.processed).tab().append(vectors.candidates).nl();
		report.append("shard\trows\tpositive_labels\n");
		for(int i=0; i<vectors.rows.length; i++){
			report.append(i).tab().append(vectors.rows[i]).tab().append(vectors.positives[i]).nl();
			if(vectors.rows[i]==0 || vectors.positives[i]==0 || vectors.positives[i]==vectors.rows[i]){
				throw new IllegalArgumentException("Each native training/monitoring shard needs both target classes: "+i);
			}
		}
		write(new File(directory, "learning.tsv"), report);
		return vectors;
	}

	/** Describes the requested network and endpoint tables; a null network represents OFF. */
	private Boundary boundary(final File net, final File directory){
		return new Boundary(net==null ? null : net.getPath(), directory==null ? null : new File(directory, "start.tsv").getPath(),
			directory==null ? null : new File(directory, "stop.tsv").getPath(), startRadius, stopRadius, startMargin, stopMargin, cutoff);
	}

	/** Copies the accepted library and boundary policy into a fresh loadable family directory. */
	private void export(final NcrnaFamily accepted, final Boundary policy, final File root){
		NcrnaFamilyExport.write(accepted, ref, hbm, seeds, policy, root);
	}

	/** Creates one stage directory beneath this invocation's work directory. */
	private File fresh(final String name){
		final File directory=new File(work, name);
		if(!directory.mkdir()){throw new IllegalStateException("Could not create fresh stage directory: "+directory);}
		return directory;
	}

	/** The aggregate child heap budget equals the invoking -Xmx, in addition to the parent heap. */
	private static long childHeap(){return Math.max(1, Runtime.getRuntime().maxMemory()/(1024L*1024));}

	/*--------------------------------------------------------------*/
	/*----------------         Caller Trials        ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts at most min(t,3) coarse tasks, dividing worker and heap budgets.
	 * ThreadWaiter deliberately retries interrupts; child ownership instead requires
	 * cancelling and joining every started worker before interruption propagates.
	 */
	private File[] callTrials(final NcrnaCallerEvaluation evaluation, final File[] roots){
		final int concurrency=Math.min(threads, roots.length);
		final File[] results=new File[roots.length];
		for(int first=0; first<roots.length; first+=concurrency){
			final int count=Math.min(concurrency, roots.length-first);
			final ArrayBlockingQueue<Trial> completed=new ArrayBlockingQueue<Trial>(count);
			final Trial[] workers=new Trial[count];
			int started=0;
			boolean interrupted=false;
			try{
				if(Thread.currentThread().isInterrupted()){throw new InterruptedException("Caller trials interrupted before launch");}
				for(int i=0; i<count; i++){
					final int id=first+i;
					workers[i]=new Trial(evaluation, roots[id], new File(work, "caller_"+id),
						threads/count+(i<threads%count ? 1 : 0), Math.max(1, childHeap()/count), completed);
					workers[i].start();
					started++;
				}
				for(int i=0; i<count; i++){
					final Trial done=completed.take();
					if(done.failure!=null){throw new IllegalStateException("Boundary caller trial failed: "+done.directory, done.failure);}
				}
				for(int i=0; i<count; i++){results[first+i]=workers[i].result;}
			}catch(InterruptedException failure){
				interrupted=true;
				throw new IllegalStateException("Boundary caller trials interrupted", failure);
			}finally{
				for(int i=0; i<started; i++){
					if(workers[i].isAlive()){workers[i].interrupt();}
				}
				for(int i=0; i<started; i++){
					while(workers[i].isAlive()){
						try{workers[i].join();}catch(InterruptedException again){interrupted=true;}
					}
				}
				if(interrupted){Thread.currentThread().interrupt();}
			}
		}
		return results;
	}

	/** One ordinary caller; queue handoff publishes result/failure to the owner. */
	private static final class Trial extends Thread{
		Trial(final NcrnaCallerEvaluation evaluation_, final File root_, final File directory_,
				final int threads_, final long heap_, final ArrayBlockingQueue<Trial> done_){
			evaluation=evaluation_;
			root=root_;
			directory=directory_;
			threads=threads_;
			heap=heap_;
			done=done_;
		}
		/** Publishes the call result or failure before notifying the owning thread. */
		@Override
		public void run(){
			try{result=evaluation.run(root, directory, threads, heap);}catch(Throwable caught){failure=caught;}finally{done.add(this);}
		}
		final NcrnaCallerEvaluation evaluation;
		final File root, directory;
		final int threads;
		final long heap;
		final ArrayBlockingQueue<Trial> done;
		File result;
		Throwable failure;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Reports           ----------------*/
	/*--------------------------------------------------------------*/

	/** Records original-model assignments, source-group splits and locus eligibility. */
	private void writePopulation(final NcrnaLearningPopulation population){
		try(Output output=new Output(new File(work, "population.tsv").getPath())){
			final ByteStreamWriter writer=output.writer;
			writer.print("locus\tsource\ttaxid\toriginal_model\tvalidation\teligible\n");
			final ByteBuilder row=new ByteBuilder();
			final NcrnaLearningSplit.Result split=population.split();
			for(int gene=0; gene<population.ids.length; gene++){
				row.clear();
				row.append(population.ids[gene]).tab().append(population.sourceName(gene)).tab();
				row.append(split.taxids[split.groupByGene[gene]]).tab().append(population.models[population.originalModel[gene]]).tab();
				row.append(split.validationByGene[gene] ? 1 : 0).tab().append(population.eligible(gene) ? 1 : 0).nl();
				writer.print(row);
			}
		}
	}

	/** Post-refit metrics remain separately visible even when eligibility sends the final product to OFF. */
	private void writeRefit(final Result baseline, final Result refit){
		final Change change=new Change(baseline, refit);
		final ByteBuilder report=new ByteBuilder("tp\tfp\tfn\tprecision\trecall\texact_5prime\texact_3prime\tboth_exact\tlost\tgained\tboth_exact_gain\n");
		report.append(refit.matched).tab().append(refit.unmatched).tab().append(refit.missed).tab();
		if(refit.calls.isEmpty()){report.append("NA");}else{report.append(refit.precision(), 8, true);}
		report.tab().append(refit.recall(), 8, true).tab();
		report.append(refit.exactStart).tab().append(refit.exactStop).tab().append(refit.exactBoth).tab();
		report.append(change.lost).tab().append(change.gained).tab().append(change.exactGain).nl();
		write(new File(work, "refit.tsv"), report);
	}

	/** Records the tuning choice separately from the final post-refit acceptance result. */
	private void writeDecision(final String status, final Method selected, final Result baseline, final Result result){
		final ByteBuilder report=new ByteBuilder("status\t").append(status).nl();
		report.append("tuning_method\t").append(selected==null ? "OFF" : selected.name()).nl();
		report.append("split_seed\t").append(splitSeed).nl();
		report.append("max_lost\t").append(maxLost).nl();
		report.append("min_precision\t").append(minPrecision, 8, true).nl();
		report.append("tie_preference\t").append(tiePreference.name()).nl();
		report.append("baseline_both_exact\t").append(baseline.exactBoth).nl();
		report.append("final_both_exact\t").append(result.exactBoth).nl();
		report.append("refit_evaluation\t").append(selected==null ? "Not performed; neither tuning candidate was eligible."
			: "All-eligible coverage evidence, not independent holdout performance.").nl();
		write(new File(work, "summary.tsv"), report);
	}

	/** Writes a fresh report and propagates writer failures through the output owner. */
	private static void write(final File path, final ByteBuilder text){
		try(Output output=new Output(path.getPath())){
			output.writer.print(text);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Complete source-context FASTA and fresh exported family directory. */
	private String in, out;
	/** Source roster, locus ledger and optional excluded original-locus identifiers. */
	private String roster, loci, exclude;
	/** Original support graph/consensus and the accepted consensus/HBM/seed paths. */
	private String graph, original, ref, hbm, seeds;
	/** Registered caller member and its label in the independent evaluation truth. */
	private String family, truthFamily;
	/** Evaluation genome list, truth table and native caller/trainer recipes. */
	private String genomes, truth, callerConfig, classifierConfig, regressionConfig;
	/** Inference feature contract; only v1 is currently supported. */
	private String features="v1";
	/** Total worker budget and number of whole-locus tuning shards. */
	private int threads=1, shards=1;
	/** Candidate search distances in bases from each annotated endpoint. */
	private int startRadius=-1, stopRadius=-1;
	/** Endpoint word geometry in bases; inside plus outside is the k-mer length. */
	private int startInside=-1, startOutside=-1, stopInside=-1, stopOutside=-1;
	/** Maximum permitted baseline detections lost during method selection/refit. */
	private int maxLost=0;
	/** Explicit seed for assigning entire source TaxID groups to validation. */
	private long splitSeed=-1;
	/** Required movement score gains and optional absolute boundary score floor. */
	private float startMargin=0, stopMargin=0, cutoff=Float.NaN;
	/** Minimum end-to-end detection precision for an eligible candidate. */
	private double minPrecision=.95;
	/** Method selected when eligible classifier and regression gains tie. */
	private Method tiePreference=Method.REGRESSION;
	/** Final export and separate intermediate-work directory, both fresh on entry. */
	private File destination, work;
	/** Diagnostic stream owned by this invocation unless it is standard output/error. */
	private PrintStream outstream=System.err;
}
