package prot;

import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import parse.Parser;
import prok.CallGenes;
import prok.GeneModel;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/**
 * Public assembly-batch engine. Taxonomy is prepared before PGMs and network
 * resources; each gene/inference worker owns its scratch while sharing frozen
 * bindings. A report is published only after every input bin succeeds.
 * @author Yoimiya
 */
public final class MagQCAssemblyBatch {

	/** Runs the native public client; the prok.ProkCC facade preserves its public name. */
	public static void main(String[] args) throws Exception{
		new MagQCAssemblyBatch(args).process();
	}

	/** Parses one optional relocatable release config and ordinary public flags. */
	MagQCAssemblyBatch(String[] args) throws IOException{
		this(args, false);
	}

	/** The package fixture may retain raw rows to compare against the independent harness. */
	MagQCAssemblyBatch(String[] args, boolean retainVectors_) throws IOException{
		retainVectors=retainVectors_;
		options=parseOptions(args);
		sharedComposite=lockedMode(required("compositemode"));
		sharedSubnets=lockedMode(required("subnetmode"));
		parallelLoad=parallelMode(required("loadmode"));
		swapNL=Parse.parseBoolean(required("swapnl"));
		verbose=Parse.parseBoolean(required("verbose"));
		timings=options.containsKey("timings") && Parse.parseBoolean(required("timings"));
		jobs=inputs(required("in"));
		threads=Math.min(jobs.size(), options.containsKey("t") ? Parse.parseIntKMG(required("t")) : Shared.threads());
		if(threads<1){throw new IllegalArgumentException("t must be positive");}
		passes=Integer.parseInt(required("passes"));
		if(passes<1){throw new IllegalArgumentException("passes must be positive");}
		if(!required("policy").equals("BOUNDED_LOOKAHEAD") || !required("lookahead").equals("4")){
			throw new IllegalArgumentException("This release requires policy=BOUNDED_LOOKAHEAD lookahead=4");
		}
		if(!required("pgmmode").equals("taxonomy") && !required("pgmmode").equals("default")){
			throw new IllegalArgumentException("pgmmode must be taxonomy or default");
		}
		taxMode=required("taxmode");
		if(!taxMode.equals("server") && !taxMode.equals("local")){
			throw new IllegalArgumentException("taxmode must be server or local");
		}
		override=MagQCAssemblyInput.overrideTaxonomy(options);
		compMultiplier=multiplier(required("comperrormultiplier"));
		contamMultiplier=multiplier(required("contamerrormultiplier"));
		String output=options.get("out");
		if("null".equalsIgnoreCase(output)){output=null;}
		else if("-".equals(output)){output="stdout";}
		final boolean overwrite=Parse.parseBoolean(required("ow"));
		// This publisher already visits jobs in order; ByteStreamWriter.print uses its ordinary buffer.
		ffout=FileFormat.testOutput(output, FileFormat.TEXT, "txt", true, overwrite, false, false);
		if(ffout!=null){
			if(!ffout.canWrite()){throw new IllegalArgumentException("Cannot write output: "+output+" (overwrite="+overwrite+")");}
			if(!ffout.stdio() && !ffout.devnull()){
				final Path target=Paths.get(ffout.name()).toAbsolutePath().normalize();
				final String[] names=new String[jobs.size()+1];
				for(int i=0; i<jobs.size(); i++){
					names[i]=jobs.get(i).path;
					if(Files.exists(target) && Files.isSameFile(target, Paths.get(names[i]))){
						throw new IllegalArgumentException("Output must not overwrite an input assembly: "+output);
					}
				}
				names[jobs.size()]=target.toString();
				Tools.testForDuplicateFiles(true, names);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Execution          ----------------*/
	/*--------------------------------------------------------------*/

	/** Reports the optional model/sidecar assets before loading or service access. */
	static void requireOptionalAssets(HashMap<String,String> options){
		assert(options!=null) : "Resource preflight requires the parsed release configuration";
		boolean missing=false;
		for(String key:OPTIONAL_ASSETS){
			final String path=options.get(key);
			if(path==null || !Files.isRegularFile(Paths.get(path))){missing=true;}
		}
		if(!missing){return;}
		final StringBuilder message=new StringBuilder("Missing ProkCC model assets. Install the files named by your release config:\n");
		for(String key:OPTIONAL_ASSETS){
			final String path=options.get(key);
			message.append("  ").append(key).append(": ").append(path==null ? "not configured" : path).append('\n');
			if(path!=null){
				final String url=shared.Resources.downloadURL(path);
				if(url!=null){message.append("    Download: ").append(url).append('\n');}
			}
		}
		message.append(shared.Resources.prokccDownloadInstructions())
			.append("Use the matching release config and metadata; custom configs retain their explicit paths.");
		throw new IllegalArgumentException(message.toString());
	}

	/** Configures each global-dependent phase before starting its private workers. */
	void process() throws Exception{
		final long processStart=timestamp();
		prok.ProkObject.verbose=verbose;
		requireOptionalAssets(options);
		if(Parse.parseBoolean(required("deterministic"))){
			Shared.SIMD_FMA=false; Shared.SIMD_FEED_FORWARD=false;
		}
		// Bind operator-supplied resources before any service request or output creation.
		for(String name:PINNED_RESOURCES){MagQCNetworkHarness.pinned(options, name, name+"sha80");}
		if(taxMode.equals("local")){MagQCNetworkHarness.pinned(options, "taxsketch", "taxsketchsha80");}
		long phaseStart=timestamp();
		if(parallelLoad){loadInParallel();}
		else{binding=MagQCAssemblyInput.loadBinding(options);}
		final long bindingNanos=elapsed(phaseStart);
		if(override==null){
			try(MagQCAssemblyInput.SketchSession session=MagQCAssemblyInput.openSketchSession(MagQCAssemblyInput.normalSearch(options))){
				sketches=session;
				runWorkers(true);
			}finally{sketches=null;}
		}else{runWorkers(true);}
		// PGM parsing changes caller geometry: finish all distinct loads before callers start.
		phaseStart=timestamp();
		final ArrayList<GeneModel> models=new ArrayList<GeneModel>();
		for(Job job:jobs){
			job.model=required("pgmmode").equals("default") ? GeneCallAdapter.defaultModel() :
				CallGenes.getPhylumPGM(job.taxonomy.phylum.equals("unknown") ? null : job.taxonomy.phylum);
			if(!models.contains(job.model)){models.add(job.model);}
		}
		long callerSetupNanos,networkNanos;
		try(FastaInCacheRowBuilder.CallerSession session=FastaInCacheRowBuilder.openSession(models.toArray(new GeneModel[models.size()]))){
			callers=session;
			callerSetupNanos=elapsed(phaseStart);
			phaseStart=timestamp();
			if(!parallelLoad){
				master=loadSubnets();
				scorer=sharedComposite ? new MagQCNetworkHarness.Scorer(
					new MagQCCompositeInference.Binding(options, master.preparedInputWidth())) :
					new MagQCNetworkHarness.Scorer(options, master.preparedInputWidth());
				validateLoadedNetworks();
			}
			networkNanos=parallelLoad ? 0 : elapsed(phaseStart);
			runWorkers(false);
		}finally{callers=null;}
		publish();
		for(int i=0; i<jobs.size(); i++){
			if(i>0){System.err.println();}
			System.err.print(jobs.get(i).humanReport);
		}
		if(timings){printTimings(bindingNanos,callerSetupNanos,networkNanos,elapsed(processStart));}
		if(verbose){
			System.err.println("ProkCC completed "+jobs.size()+" bins with "+threads+" workers; error estimates are "+required("errorfitset"));
		}
	}

	/** Three independent resource loads join before taxonomy/calling workers use their state. */
	private void loadInParallel(){
		final ArrayList<ResourceThread> loaders=new ArrayList<ResourceThread>(3);
		for(int part=0; part<3; part++){loaders.add(new ResourceThread(this, part));}
		startAndJoin(loaders);
		for(ResourceThread loader:loaders){
			if(loader.failure!=null){throw new RuntimeException("ProkCC resource load failed: "+loader.getName(), loader.failure);}
			assert(loader.success) : "Joined resource loader must finish or retain its exception";
		}
		validateLoadedNetworks();
	}

	/** The subnet loader owns its formatter construction; no global CellNet mode is changed. */
	private MagQCVectorMaker loadSubnets(){
		return MagQCVectorMaker.initializePrepared(required("bundle"), required("familylist"), required("subnetmanifest"),
			required("subnetmanifestsha80"), required("expectedcopytable"), required("expectedcopytablesha80"),
			required("subnetpopulations"), required("subnetpopulationssha80"), verbose);
	}

	/** Joins establish publication, then independent model and formatter layouts must agree. */
	private void validateLoadedNetworks(){
		if(scorer.input.length!=master.preparedInputWidth()){
			throw new IllegalArgumentException("Composite input width differs from the loaded subnet layout");
		}
		if(!sharedComposite && (scorer.dummy || scorer.worker.getTag("magqc_fixture")!=null)){
			throw new IllegalArgumentException("The public client requires a real six-output model, not fixture weights");
		}
	}

	/** One coarse task per resource kind, with failure retained until every started loader joins. */
	private static final class ResourceThread extends Thread{
		ResourceThread(MagQCAssemblyBatch owner_, int part_){
			super("prokcc-load-"+part_); owner=owner_; part=part_;
			assert(part>=0 && part<3) : "Resource loader index selects HBM, subnets or composite";
		}
		@Override public void run(){
			final long start=owner.timestamp();
			try{
				if(part==0){owner.binding=MagQCAssemblyInput.loadBinding(owner.options);}
				else if(part==1){owner.master=owner.loadSubnets();}
				else if(owner.sharedComposite){
					owner.scorer=new MagQCNetworkHarness.Scorer(prok.ProkCC.preloadCCSynced(owner.options));
				}else{owner.scorer=MagQCNetworkHarness.Scorer.loadConfigured(owner.options, true);}
				success=true;
			}catch(Throwable t){failure=t;}
			finally{owner.resourceNanos[part]=owner.elapsed(start);}
		}
		private final MagQCAssemblyBatch owner;
		private final int part;
		private Throwable failure;
		private boolean success;
	}

	/** Optional monotonic phase clock; disabled runs do not collect timings. */
	private long timestamp(){return timings ? System.nanoTime() : 0;}
	private long elapsed(long start){return timings ? System.nanoTime()-start : 0;}

	/**
	 * Reports wall time at real pipeline boundaries. Worker phases are summed
	 * across bins, so only a single-bin run has an additive wall-time split.
	 * Binding includes HBM/shortlist validation; inference includes feature building
	 * and lazy worker network copies. Unattributed I/O/setup remains explicit.
	 */
	private void printTimings(long bindingNanos,long callerSetupNanos,long networkNanos,long totalNanos){
		assert(timings && totalNanos>=0) : "Phase timings require one enabled monotonic process interval";
		long taxonomy=0,calling=callerSetupNanos,inference=0;
		for(Job job:jobs){taxonomy+=job.taxonomyNanos; calling+=job.callingNanos; inference+=job.inferenceNanos;}
		final ByteBuilder out=new ByteBuilder(1024);
		out.append("#prokcc_timings\tphase\tseconds\n");
		if(parallelLoad){
			appendTiming(out,"parallel_resource_load_wall",bindingNanos);
			appendTiming(out,"parallel_hbm_binding_task",resourceNanos[0]);
			appendTiming(out,"parallel_subnet_task",resourceNanos[1]);
			appendTiming(out,"parallel_composite_task",resourceNanos[2]);
		}else{
			appendTiming(out,"hbm_and_assignment_binding",bindingNanos);
			appendTiming(out,"subnet_and_composite_load",networkNanos);
		}
		appendTiming(out,"quickclade_call_worker_sum",taxonomy);
		appendTiming(out,"gene_call_and_assignment_worker_sum",calling);
		appendTiming(out,"features_and_inference_worker_sum",inference);
		if(jobs.size()==1){appendTiming(out,"other_setup_io_and_report",totalNanos-bindingNanos-networkNanos-taxonomy-calling-inference);}
		appendTiming(out,"java_process",totalNanos);
		System.err.print(out.toString());
	}

	/** One diagnostic row; no timing result is mixed into biological score output. */
	private static void appendTiming(ByteBuilder out,String phase,long nanos){
		assert(nanos>=0) : "Disjoint single-bin timing intervals must not exceed the enclosing process interval: "+phase;
		out.append("prokcc_timing\t").append(phase).tab().appendSlow(nanos/1e9).nl();
	}

	/** Uses coarse workers; all started workers join before any resource session closes. */
	private void runWorkers(boolean prepare){
		final AtomicInteger next=new AtomicInteger();
		final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){workers.add(new ProcessThread(this, next, prepare));}
		startAndJoin(workers);
		for(ProcessThread worker:workers){
			if(worker.failure!=null){throw new RuntimeException("ProkCC bin processing failed", worker.failure);}
			if(!worker.success){throw new IllegalStateException("ProkCC worker ended without success");}
		}
	}

	/** Joins the successfully started prefix even when native thread creation fails. */
	void startAndJoin(ArrayList<? extends Thread> workers){
		int started=0;
		boolean interrupted=false;
		try{
			for(Thread worker:workers){worker.start(); started++;}
		}catch(Throwable t){
			failure=t;
			throw new RuntimeException("Could not start ProkCC workers", t);
		}finally{
			// ThreadWaiter cannot wait on NEW threads: they never reach TERMINATED.
			// The preallocated list retains the started prefix without allocating at OOM.
			for(int i=0; i<started; i++){
				final Thread worker=workers.get(i);
				while(worker.isAlive()){
					try{worker.join();}catch(InterruptedException e){failure=e; interrupted=true;}
				}
			}
			if(interrupted){Thread.currentThread().interrupt();}
		}
		if(failure!=null){throw new RuntimeException("ProkCC workers failed", failure);}
	}

	/** One private formatter, scorer, alignment scratch and one assembly at a time per worker. */
	private static final class ProcessThread extends Thread{
		ProcessThread(MagQCAssemblyBatch owner_, AtomicInteger next_, boolean prepare_){
			owner=owner_; next=next_; prepare=prepare_;
		}
		@Override public void run(){
			try{
				final MagQCVectorMaker formatter=prepare ? null : owner.master.newPreparedWorker(owner.sharedSubnets);
				final MagQCNetworkHarness.Scorer scoring=prepare ? null : owner.scorer.newWorker();
				final ProteinSearcher searcher=prepare ? null : new ProteinSearcher();
				if(searcher!=null){searcher.aligner="d55";}
				final ByteBuilder vector=new ByteBuilder(131072), output=new ByteBuilder(1024);
				for(int index=next.getAndIncrement(); index<owner.jobs.size() && owner.failure==null; index=next.getAndIncrement()){
					final Job job=owner.jobs.get(index);
					final long binStart=System.nanoTime();
					if(prepare){
						job.inputPin=MagQCNetworkHarness.sha80(job.path);
						final ArrayList<Read> contigs=MagQCAssemblyInput.readContigs(job.path);
						final long taxonomyStart=owner.timestamp();
						job.taxonomy=owner.override==null ? owner.sketches.classify(contigs, owner.options) : owner.override;
						job.taxonomyNanos=owner.elapsed(taxonomyStart);
						job.workerNanos=System.nanoTime()-binStart;
					}else{
						final ArrayList<Read> contigs=MagQCAssemblyInput.readContigs(job.path);
						if(!job.inputPin.equals(MagQCNetworkHarness.sha80(job.path))){
							throw new IllegalArgumentException("Assembly changed after taxonomy preparation: "+job.path);
						}
						final MagQCAssemblyInput.Taxonomy tax=job.taxonomy;
						final long callingStart=owner.timestamp();
						final MagQCPreparedBin bin=owner.callers.build(contigs, job.path, tax.status, tax.domain, tax.phylum,
							job.model, owner.passes, owner.binding, searcher, ProteinSearcher.AssignPolicy.BOUNDED_LOOKAHEAD, 4);
						job.callingNanos=owner.elapsed(callingStart);
						final long inferenceStart=owner.timestamp();
						vector.clear(); formatter.formatPreparedBin(vector, bin);
						if(owner.retainVectors){job.vector=vector.toBytes();}
						scoring.evaluate(vector, job.path);
						job.inferenceNanos=owner.elapsed(inferenceStart);
						job.report=new MagQCAssemblyReport(contigs, bin, owner.swapNL);
						job.workerNanos+=System.nanoTime()-binStart;
						output.clear(); owner.appendScore(output, scoring, job, index);
						job.result=output.toBytes();
					}
				}
				success=true;
			}catch(Throwable t){failure=t; owner.failure=t;}
		}
		private final MagQCAssemblyBatch owner;
		private final AtomicInteger next;
		private final boolean prepare;
		private Throwable failure;
		private boolean success;
	}

	/** Keeps raw heads unchanged, appending provenance and the two configured estimates. */
	private void appendScore(ByteBuilder out, MagQCNetworkHarness.Scorer scoring, Job job, int ordinal){
		assert(scoring.names.length==6) : "The validated public model must supply the six frozen output heads";
		out.append(ordinal).tab().append(job.path);
		for(int i=0; i<6; i++){out.tab().appendSlow(scoring.getOutput(i));}
		final double compError=scoring.getOutput(2), contamError=scoring.getOutput(3);
		if(compError<0 || contamError<0){throw new IllegalStateException("Negative predicted error: "+job.path);}
		final double comp=compError*compMultiplier, contam=contamError*contamMultiplier;
		if(!Double.isFinite(comp) || !Double.isFinite(contam)){throw new IllegalArgumentException("Error multiplier overflow: "+job.path);}
		out.tab().appendSlow(comp).tab().appendSlow(contam).tab().append(MagQCAssemblyInput.taxonomySource(options));
		out.tab().append(job.taxonomy.status).tab().append(job.taxonomy.domain).tab().append(job.taxonomy.phylum);
		out.tab().append(job.inputPin);
		job.report.append(out, job.taxonomy, scoring.getOutput(0), scoring.getOutput(1));
		assert(job.workerNanos>=0) : "Bin worker time sums its two completed monotonic processing intervals";
		out.tab().appendSlow(job.workerNanos/1e9);
		out.nl();
		job.humanReport=job.report.human(job.taxonomy, scoring.getOutput(0), scoring.getOutput(1), comp, contam);
	}

	/** Publishes optional TSV data only after all bins succeed; screen reports are separate stderr output. */
	private void publish() throws IOException{
		if(ffout==null){return;}// Omitted out or explicit out=null disables data output.
		if(!ffout.stdio() && !ffout.devnull()){
			Files.createDirectories(Paths.get(ffout.name()).toAbsolutePath().getParent());
		}
		final ByteStreamWriter writer=new ByteStreamWriter(ffout);
		writer.start();
		try{
			final ByteBuilder header=new ByteBuilder(2048);
			header.append("#schema_version\tprokcc_assembly_scores_v4\n");
			MagQCNetworkHarness.appendBindings(header, options);
			header.append("#netsha80\t").append(required("netsha80")).nl();
			appendReportOptions(header);
			header.append("#error_estimates\tpredicted absolute errors in fraction units; not confidence intervals\n");
			header.append("#assembly_statistics\twhole FASTA records; GC=GC/ACGT; coding density=summed CDS bp/assembly bp; quality=RNA-aware BinStats.type\n");
			header.append("#bin_worker_wall_seconds\tsum of per-bin taxonomy and calling/inference worker intervals, including input I/O; excludes shared setup, dispatch queues, phase barriers and report publication\n");
			header.append("#columns\trow_index\tbin_id");
			for(String name:scorer.names){header.tab().append(name);}
			header.append("\tscaled_error_gene_completeness\tscaled_error_gene_contamination\ttaxonomy_source\tqc_status\tqc_domain\tqc_phylum\tinput_sha80");
			header.append(MagQCAssemblyReport.columns(swapNL)).append("\tbin_worker_wall_seconds\n");
			writer.print(header);
			for(Job job:jobs){
				if(job.result==null){throw new IllegalStateException("Missing completed bin: "+job.path);}
				writer.print(job.result);
			}
		}finally{
			if(writer.poisonAndWait()){throw new IOException("ProkCC report writer failed");}
		}
		if(ffout.stdio() && System.out.checkError()){throw new IOException("ProkCC stdout publication failed");}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Input Parsing        ----------------*/
	/*--------------------------------------------------------------*/

	/** Resource paths use the config directory; user paths use the working directory. A bare input file is accepted. */
	static HashMap<String,String> parseOptions(String[] args){
		Path base=Paths.get("").toAbsolutePath();
		boolean configSeen=false;
		for(String arg:args){
			if(arg.toLowerCase(Locale.ROOT).startsWith("config=")){
				final String value=arg.substring(7);
				if(configSeen || value.isEmpty() || value.indexOf(',')>=0){throw new IllegalArgumentException("Use one release config file");}
				base=Paths.get(value).toAbsolutePath().getParent(); configSeen=true;
			}
		}
		final HashMap<String,String> values=new HashMap<String,String>();
		// Expand config with the native reader, but do not run PreParser's output
		// redirection or global setters before this public allowlist is validated.
		for(String arg:Parser.parseConfig(args)){
			// The wrapper prepends config=, so recognize the first bare input even
			// after config expansion. Explicit or repeated inputs still fail below.
			if(arg.indexOf('=')<0 && !values.containsKey("in") && Files.isRegularFile(Paths.get(arg))){arg="in="+arg;}
			final int equals=arg.indexOf('=');
			if(equals<1){throw new IllegalArgumentException("Expected key=value: "+arg);}
			String key=arg.substring(0, equals).toLowerCase(Locale.ROOT);
			if(key.equals("threads")){key="t";}
			if(key.equals("swapln")){key="swapnl";}
			if(key.equals("overwrite")){key="ow";}
			final String value=arg.substring(equals+1);
			if(!ALLOWED.contains(key) || values.put(key, value)!=null || value.isEmpty()){
				throw new IllegalArgumentException("Unknown, duplicate or empty public option: "+key);
			}
			for(int i=0; i<value.length(); i++){
				if(Character.isISOControl(value.charAt(i))){throw new IllegalArgumentException("Option must fit one report field: "+key);}
			}
		}
		for(String key:RESOURCE_PATHS){
			if(values.containsKey(key)){
				final Path path=Paths.get(values.get(key));
				values.put(key, (path.isAbsolute() ? path : base.resolve(path)).normalize().toString());
			}
		}
		defaults(values, "policy", "BOUNDED_LOOKAHEAD", "lookahead", "4", "passes", "1", "pgmmode", "taxonomy",
			"taxaddress", "refseq", "taxmode", "server", "normalsearch", "t", "deterministic", "t", "comperrormultiplier", "1.0", "contamerrormultiplier", "1.0",
			"errorfitset", "UNCALIBRATED", "errorfitdate", "NA", "errorcoverage", "NA",
			"compositemode", "locked", "subnetmode", "locked", "loadmode", "parallel", "swapnl", "f", "verbose", "f", "ow", "t");
		lockedMode(values.get("compositemode")); lockedMode(values.get("subnetmode"));
		parallelMode(values.get("loadmode"));
		MagQCAssemblyInput.normalSearch(values);
		return values;
	}

	/** Selects shared locked inference or explicit private workers for comparison. */
	static boolean lockedMode(String value){
		if("locked".equals(value)){return true;}
		if("worker".equals(value)){return false;}
		throw new IllegalArgumentException("Inference mode must be worker or locked: "+value);
	}

	/** Resource-loading concurrency is independent of the number of bin-processing workers. */
	static boolean parallelMode(String value){
		if("parallel".equals(value)){return true;}
		if("serial".equals(value)){return false;}
		throw new IllegalArgumentException("loadmode must be serial or parallel: "+value);
	}

	/** Enumerates a requested directory once, never recursively; an explicit list retains its order. */
	static ArrayList<Job> inputs(String input) throws IOException{
		final ArrayList<Path> paths=new ArrayList<Path>();
		final Path single=Paths.get(input);
		if(Files.isDirectory(single)){
			try(DirectoryStream<Path> entries=Files.newDirectoryStream(single)){
				for(Path path:entries){if(Files.isRegularFile(path) && fastaName(path)){paths.add(path);}}
			}
			Collections.sort(paths);
		}else{for(String path:input.split(",", -1)){paths.add(Paths.get(path));}}
		if(paths.isEmpty()){throw new IllegalArgumentException("No assembly FASTAs supplied");}
		final HashSet<Path> unique=new HashSet<Path>();
		final ArrayList<Job> jobs=new ArrayList<Job>(paths.size());
		for(Path path:paths){
			if(!Files.isRegularFile(path) || !Files.isReadable(path) || !fastaName(path)){
				throw new IllegalArgumentException("Input must be a readable assembly FASTA file: "+path);
			}
			final Path canonical=path.toRealPath();
			if(!unique.add(canonical)){throw new IllegalArgumentException("Duplicate assembly: "+path);}
			final String name=path.toAbsolutePath().normalize().toString();
			for(int i=0; i<name.length(); i++){
				if(Character.isISOControl(name.charAt(i))){throw new IllegalArgumentException("Assembly path must fit one report field");}
			}
			jobs.add(new Job(name));
		}
		return jobs;
	}

	/** Supported ordinary FASTA names, including native gzip/bzip2/xz compressed inputs. */
	private static boolean fastaName(Path path){
		String name=path.getFileName().toString().toLowerCase(Locale.ROOT);
		for(String suffix:new String[]{".gz", ".bz2", ".xz"}){
			if(name.endsWith(suffix)){name=name.substring(0, name.length()-suffix.length()); break;}
		}
		return name.endsWith(".fa") || name.endsWith(".fna") || name.endsWith(".fasta") || name.endsWith(".fas");
	}

	/** Defaults are explicit reportable release behavior, not hidden data-dependent values. */
	private static void defaults(HashMap<String,String> values, String... pairs){
		assert((pairs.length&1)==0) : "Option defaults are alternating key/value pairs";
		for(int i=0; i<pairs.length; i+=2){if(!values.containsKey(pairs[i])){values.put(pairs[i], pairs[i+1]);}}
	}
	private String required(String key){return MagQCNetworkHarness.required(options, key);}

	/** Records the search policy and, when selected, the pinned local taxonomy input. */
	void appendReportOptions(ByteBuilder header){
		for(String key:REPORT_OPTIONS){header.append('#').append(key).tab().append(required(key)).nl();}
		if(taxMode.equals("local")){
			header.append("#taxmode\tlocal\n#taxsketch\t").append(required("taxsketch"))
				.nl().append("#taxsketchsha80\t").append(required("taxsketchsha80")).nl();
		}
	}

	/** D231 requires positive finite constants and no silent clipping. */
	static double multiplier(String value){
		final double number=Double.parseDouble(value);
		if(!Double.isFinite(number) || number<=0){throw new IllegalArgumentException("Error multipliers must be positive and finite");}
		return number;
	}

	/** Only small durable metadata and final report bytes are retained between phases. */
	static final class Job{
		Job(String path_){path=path_;}
		final String path;
		String inputPin;
		MagQCAssemblyInput.Taxonomy taxonomy;
		GeneModel model;
		byte[] result;
		byte[] vector;
		MagQCAssemblyReport report;
		String humanReport;
		long taxonomyNanos,callingNanos,inferenceNanos;
		long workerNanos;
	}

	private final HashMap<String,String> options;
	final ArrayList<Job> jobs;
	private final boolean retainVectors;
	private final boolean timings, sharedComposite, sharedSubnets, parallelLoad, swapNL, verbose;
	private final long[] resourceNanos=new long[3];
	private final int threads, passes;
	private final String taxMode;
	private final double compMultiplier, contamMultiplier;
	private final FileFormat ffout;
	private final MagQCAssemblyInput.Taxonomy override;
	private volatile Throwable failure;
	private MagQCAssemblyInput.SketchSession sketches;
	private FastaInCacheRowBuilder.CallerSession callers;
	private ProteinSearcher.AssignmentBinding binding;
	private MagQCVectorMaker master;
	private MagQCNetworkHarness.Scorer scorer;
	private static final String[] PINNED_RESOURCES={"bundle", "familylist", "subnetmanifest", "expectedcopytable", "subnetpopulations", "net"};
	private static final String[] OPTIONAL_ASSETS={"net", "bundle", "hbmbundle", "sidecar"};
	private static final String[] RESOURCE_PATHS={"bundle", "familylist", "subnetmanifest", "expectedcopytable", "subnetpopulations", "net",
		"profile", "roster", "ref", "rolemanifest", "core", "coveringsets", "sidecar", "hbmbundle", "hbmprovenance", "taxsketch"};
	private static final String[] REPORT_OPTIONS={"profilesha80", "policy", "lookahead", "passes", "pgmmode", "normalsearch", "deterministic",
		"comperrormultiplier", "contamerrormultiplier", "errorfitset", "errorfitdate", "errorcoverage", "compositemode", "subnetmode", "loadmode", "swapnl"};
	private static final HashSet<String> ALLOWED=new HashSet<String>();
	static{
		ALLOWED.addAll(Arrays.asList(RESOURCE_PATHS)); ALLOWED.addAll(Arrays.asList(REPORT_OPTIONS));
		ALLOWED.addAll(Arrays.asList("in", "out", "t", "taxaddress", "taxmode", "taxsketchsha80", "taxdomain", "taxphylum", "timings", "verbose", "ow"));
		for(String key:PINNED_RESOURCES){ALLOWED.add(key+"sha80");}
	}
}
