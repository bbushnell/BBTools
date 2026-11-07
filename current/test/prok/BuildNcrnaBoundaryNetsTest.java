package test.prok;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Locale;
import java.util.Random;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.ByteFile;
import ml.CellNet;
import ml.Function;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import parse.LineParser1;
import prok.BuildNcrnaBoundaryNets;
import prok.FamilyResourceResolver;
import prok.NcrnaCallerEvaluation;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaLearningPopulation;
import prok.NcrnaSupport;
import prok.NcrnaSupportIO;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Builds twenty exact synthetic teaching loci; checks the complete builder's
 * no-change control independently through its reports and the public loader.
 * Huge movement margins force an OFF result regardless of trainer variability.
 * This is orchestration coverage, not a biological accuracy benchmark.
 * @author Brian Bushnell, Raiden
 */
public final class BuildNcrnaBoundaryNetsTest{

	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, null, false);
		try(AutoCloseable log=()->Shared.closeStream(pp.outstream)){
			run(pp.args);
		}
	}

	/** Main owns its selected log across every requested fixture mode. */
	private static void run(final String[] args) throws Exception{
		final Parser parser=new Parser();
		String mode=null, learned=null, method="regression";
		for(String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("mode")){mode=b;}
			else if(a.equals("learned")){learned=b;}
			else if(a.equals("method")){method=b;}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown fixture argument: "+arg);}
		}
		if(parser.out1==null){throw new IllegalArgumentException("Fixture requires out=");}
		final File root=new File(parser.out1).getAbsoluteFile();
		final File saved=learned==null ? new File(root, "product.work/tuning") : new File(learned).getAbsoluteFile();
		if("prepare".equals(mode)){prepare(root);}
		else if("preparefixed".equals(mode)){prepare(root, true);}
		else if("fixedpolicy".equals(mode)){fixedPolicy(root);}
		else if("checkfixed".equals(mode)){checkFixed(root, method);}
		else if("examples".equals(mode)){examples(root);}
		else if("check".equals(mode)){check(root);}
		else if("mockinputs".equals(mode)){mockInputs(root);}
		else if("checkrefit".equals(mode)){checkRefit(root);}
		else if("lifecycle".equals(mode)){lifecycle(root);}
		else if("precision".equals(mode)){precision(root);}
		else if("checkevaluate".equals(mode)){checkEvaluate(root, saved);}
		else if("reusecontrols".equals(mode)){reuseControls(root, saved);}
		else{throw new IllegalArgumentException("Unknown fixture mode: "+mode);}
	}

	/** Writes small, reproducible teaching genomes and truth, never biological evidence.
	 * Each shipped consensus appears once on the plus strand, with 200 independently
	 * generated bases on each side. Truth follows construction, never caller output.
	 */
	private static void examples(final File root) throws Exception{
		Files.createDirectory(root.toPath());
		final String ref=new FamilyResourceResolver().resolve("sixs_rf01685", "consensus", null).path;
		final NcrnaSupport.Library library=NcrnaSupport.loadLibrary(ref);
		assert(library.names.length>0) : "The installed 6S library supplies the explicit teaching sequences";
		final Random random=new Random(42);
		final int flank=200;
		final ByteBuilder bases=new ByteBuilder(), row=new ByteBuilder();
		try(Output genome=new Output(new File(root, "panel.fa").getPath());
				Output truth=new Output(new File(root, "panel.truth.tsv").getPath())){
			truth.writer.print("seqid\tfamily\tstart\tstop\tstrand\n");
			for(int model=0; model<library.names.length; model++){
				bases.clear().append(">synthetic_").append(model).nl();
				for(int i=0; i<flank; i++){bases.append("ACGT".charAt(random.nextInt(4)));}
				bases.append(library.bases[model]);
				for(int i=0; i<flank; i++){bases.append("ACGT".charAt(random.nextInt(4)));}
				genome.writer.print(bases.nl());
				row.clear().append("synthetic_").append(model).append("\tfixture\t").append(flank+1).tab();
				truth.writer.print(row.append(flank+library.bases[model].length).append("\t+\n"));
			}
		}
		System.out.println("Synthetic teaching panel written: "+library.names.length+" embedded consensuses, 200-base flanks, seed 42");
	}

	private static void prepare(final File root) throws Exception{
		prepare(root, false);
	}

	/** A one-TaxID variant proves fixed refit does not require a feasible holdout. */
	private static void prepare(final File root, final boolean singleTaxid) throws Exception{
		Files.createDirectory(root.toPath());
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final String ref=resolver.resolve("sixs_rf01685", "consensus", null).path;
		final String hbm=resolver.resolve("sixs_rf01685", "models", null).path;
		final String seeds=resolver.resolve("sixs_rf01685", "seeds", null).path;
		final NcrnaSupport.Library library=NcrnaSupport.loadLibrary(ref);
		final byte[] gene=library.bases[0];
		write(root, "original.fa", new ByteBuilder(">original\n").append(gene).nl());
		final ByteBuilder contexts=new ByteBuilder();
		final ByteBuilder roster=new ByteBuilder(NcrnaLearningPopulation.ROSTER_HEADER).nl();
		final ByteBuilder loci=new ByteBuilder(NcrnaLearningPopulation.LOCUS_HEADER).nl();
		final ByteBuilder graph=new ByteBuilder("#ncrna_support_graph\t2\n#scope\tORIGINAL\n#aligner\tquantum\n#cutoff_a48\t");
		graph.appendFloatA48(1).append("\n#models\t1\n#model\t0\toriginal\n#queries\t20\nquery_id\tbest_model_index\tbest_identity_a48\tsupporters\n");
		for(int i=0; i<20; i++){
			final int taxid=singleTaxid ? 1 : i+1;
			contexts.append(">gene").append(i).append(" source=genome").append(i)
				.append(" lflank=20 rflank=20 contig_gc=0.5\nACGTACGTACGTACGTACGT")
				.append(gene).append("ACGTACGTACGTACGTACGT\n");
			roster.append(i).append("\tgenome").append(i).append("\tassembly").append(i).tab().append(taxid)
				.append("\tbacteria\tP\tgff\tfasta\t1\t1000\t1\n");
			loci.append("gene").append(i).append("\tgenome").append(i).append("\tassembly").append(i).tab().append(taxid)
				.append("\tbacteria\tP\tcontig\t+\t101-").append(100+gene.length).tab().append(gene.length)
				.append("\t1\tid\tfalse\tfalse\tANNOTATION_FULL_COPY\t1000\tfalse\t0\t20\t20\t1\n");
			graph.append("gene").append(i).append("\t0\t").appendFloatA48(1).append("\t0:").appendFloatA48(1).nl();
		}
		write(root, "contexts.fa", contexts);
		write(root, "roster.tsv", roster);
		write(root, "loci.tsv", loci);
		write(root, "graph.tsv", graph);
		write(root, "genomes.txt", new ByteBuilder(ref).nl());
		final ByteBuilder truth=new ByteBuilder("seqid\tfamily\tstart\tstop\tstrand\n");
		for(int i=0; i<library.names.length; i++){
			truth.append(library.names[i]).append("\tfixture\t1\t").append(library.bases[i].length).append("\t+\n");
		}
		write(root, "truth.tsv", truth);
		write(root, "caller.config", new ByteBuilder("taxonomy=f\ncallcds=f\n"));
		write(root, "classifier.config", new ByteBuilder("dims=10,4,1\nnetseed=3\nbatches=2\nnets=1\ncycles=1\nscan=0\nvfraction=0\nbalance=0\nfpb=1\ntriage=0\nfinal=rslog\nsimd=f\nt=2\nprintinterval=1\now=f\n"));
		write(root, "regression.config", new ByteBuilder("dims=10,4,1\nseed=3\nepochs=2\nvfraction=0\nfinal=linear\nsimd=f\nt=1\nprintevery=1\now=t\n"));
		final ByteBuilder config=new ByteBuilder();
		for(String[] entry:new String[][]{{"in", "contexts.fa"}, {"roster", "roster.tsv"}, {"loci", "loci.tsv"},
			{"graph", "graph.tsv"}, {"original", "original.fa"}, {"genomes", "genomes.txt"}, {"truth", "truth.tsv"},
			{"callerconfig", "caller.config"}, {"classifierconfig", "classifier.config"}, {"regressionconfig", "regression.config"}}){
			config.append(entry[0]).append('=').append(new File(root, entry[1]).getPath()).nl();
		}
		config.append("ref=").append(ref).nl().append("hbm=").append(hbm).nl().append("seeds=").append(seeds).nl();
		config.append("family=sixs_rf01685\ntruthfamily=fixture\nsplitseed=42\nshards=2\nstartradius=1\nstopradius=1\n");
		config.append("startinside=3\nstartoutside=0\nstopinside=2\nstopoutside=1\nmarginstart=1000000\nmarginstop=1000000\n");
		write(root, "builder.config", config);
		if(singleTaxid){
			for(String method:new String[]{"classifier", "regression"}){
				final ByteBuilder fixed=new ByteBuilder("mode=refit\nmethod=").append(method).nl();
				for(String[] entry:new String[][]{{"in", "contexts.fa"}, {"roster", "roster.tsv"}, {"loci", "loci.tsv"},
					{"graph", "graph.tsv"}, {"original", "original.fa"}, {method+"config", method+".config"}}){
					fixed.append(entry[0]).append('=').append(new File(root, entry[1]).getPath()).nl();
				}
				fixed.append("ref=").append(ref).nl().append("hbm=").append(hbm).nl().append("seeds=").append(seeds).nl();
				fixed.append("family=sixs_rf01685\nstartradius=1\nstopradius=1\nstartinside=3\nstartoutside=0\nstopinside=2\nstopoutside=1\n");
				write(root, "refit_"+method+".config", fixed);
			}
		}
		policy();
	}

	/** Rejects ignored options before any output creation, including explicitly supplied defaults. */
	private static void fixedPolicy(final File root){
		boolean splitRejected=false;
		try{
			NcrnaLearningPopulation.load(NcrnaSupportIO.read(new File(root, "graph.tsv").getPath()),
				NcrnaSupport.loadLibrary(new File(root, "original.fa").getPath()).names,
				new File(root, "roster.tsv").getPath(), new File(root, "loci.tsv").getPath(), null, 42);
		}catch(IllegalArgumentException expected){
			splitRejected=expected.getMessage().contains("at least ten positive TaxID groups");
		}
		assert(splitRejected) : "One-group refit must not relax the original tuning split's minimum-group safeguard";
		final String config="config="+new File(root, "refit_regression.config");
		final String output="out="+new File(root, "fixed_policy");
		new BuildNcrnaBoundaryNets(new String[]{config, output});
		for(String[] entry:new String[][]{{"splitseed", "0"}, {"shards", "1"}, {"maxlost", "0"},
			{"minprecision", "0.95"}, {"tiepreference", "regression"}, {"truthfamily", "fixture"},
			{"genomes", "genomes.txt"}, {"truth", "truth.tsv"}, {"callerconfig", "caller.config"}, {"learned", "missing"}}){
			rejectFixed(new String[]{config, output, entry[0]+"="+entry[1]}, entry[0]);
		}
		rejectFixed(new String[]{config, output, "classifierconfig="+new File(root, "classifier.config")}, "selected method");
		rejectFixed(new String[]{config, output, "method=other"}, "method");
		rejectFixed(new String[]{config, output, "mode=build"}, "only in mode=refit");
		rejectFixed(new String[]{config, output, "mode=evaluate"}, "only in mode=refit");
		final ArrayList<String> missing=new ArrayList<String>();
		for(String arg:Parser.parseConfig(new String[]{config})){
			if(!arg.startsWith("method=")){missing.add(arg);}
		}
		missing.add(output);
		rejectFixed(missing.toArray(new String[missing.size()]), "explicit method");
		assert(!new File(root, "fixed_policy").exists() && !new File(root, "fixed_policy.work").exists())
			: "Invalid fixed-refit options must not create files, train or call genomes";
		System.out.println("Fixed refit policy passed: explicit method, selected config only, no ignored selection controls");
	}

	private static void rejectFixed(final String[] args, final String message){
		boolean rejected=false;
		try{new BuildNcrnaBoundaryNets(args);}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains(message);}
		assert(rejected) : "Fixed refit must reject this invalid contract before execution: "+message;
	}

	/** Checks real native training/export on twenty loci from only one TaxID group. */
	private static void checkFixed(final File root, final String method){
		assert(method.equals("classifier") || method.equals("regression")) : "Fixture checks one explicitly selected native trainer";
		final File output=new File(root, "fixed_"+method), work=new File(output.getPath()+".work");
		assert(line(new File(work, "summary.tsv"), "status\tREFIT_UNEVALUATED")
			&& line(new File(work, "summary.tsv"), "method\t"+method.toUpperCase(Locale.ROOT))
			&& line(new File(work, "summary.tsv"), "trainer_invocations\t1")
			&& line(new File(work, "summary.tsv"), "eligible_loci\t20")
			&& line(new File(work, "summary.tsv"), "withholding\tNONE"))
			: "Fixed refit reports actual all-eligible training without implying caller acceptance";
		assert(line(new File(work, "refit/learning.tsv"), "20\t0\t0\t0\t20\t120")
			&& line(new File(work, "refit/learning.tsv"), "0\t120\t40"))
			: "All twenty loci from a single TaxID teach six candidates/two positives each in one shard";
		for(String absent:new String[]{"tuning", "caller_0", "caller_1", "caller_2", "comparison.tsv", "refit/caller", "refit/vectors.validation.tsv.gz"}){
			assert(!new File(work, absent).exists()) : "Fixed refit must not tune, invoke callers or withhold monitoring data: "+absent;
		}
		int count=0;
		final LineParser1 parser=new LineParser1('\t');
		try(Input input=new Input(new File(work, "population.tsv").getPath())){
			final byte[] header=input.reader.nextLine();
			assert(header!=null && new String(header, StandardCharsets.UTF_8).equals("locus\tsource\ttaxid\toriginal_model\twithheld\teligible"))
				: "Refit population must explicitly describe withholding, not imply an unused validation experiment";
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				parser.set(row);
				assert(parser.terms()==6 && parser.parseInt(2)==1 && parser.parseInt(4)==0 && parser.parseInt(5)==1)
					: "One TaxID's twenty eligible loci must all remain unwithheld";
				count++;
			}
		}
		assert(count==20) : "The complete original locus axis survives refit metadata output";
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.resources.root=output.getPath();
		options.boundary=true;
		final NcrnaFamily loaded=NcrnaFamilyLoader.load("sixs_rf01685", options, new FamilyResourceResolver());
		assert(loaded.library.length==19 && loaded.boundary5NetTemplate!=null && loaded.boundary3NetTemplate!=null
			&& loaded.boundaryStartInside==3 && loaded.boundaryStartOutside==0
			&& loaded.boundaryStopInside==2 && loaded.boundaryStopOutside==1)
			: "Native exported network/table geometry must load for the accepted library without a caller run";
		System.out.println("Fixed refit integration passed: "+method+", twenty loci in one TaxID, one all-data training shard, loadable unevaluated export");
	}

	/** Exercises the public constructor without starting training or writing output. */
	private static void precision(final File root) throws Exception{
		final String config="config="+new File(root, "builder.config");
		final String output="out="+new File(root, "precision_policy");
		final BuildNcrnaBoundaryNets defaults=new BuildNcrnaBoundaryNets(new String[]{config, output});
		final Field field=BuildNcrnaBoundaryNets.class.getDeclaredField("minPrecision");
		field.setAccessible(true);
		assert(field.getDouble(defaults)==.95) : "Widening the accepted explicit range must not change the existing default";
		for(final String value:new String[]{"0", "0.5", "1"}){
			final BuildNcrnaBoundaryNets selected=new BuildNcrnaBoundaryNets(new String[]{config, output, "minprecision="+value});
			assert(field.getDouble(selected)==Double.parseDouble(value))
				: "The explicit precision option must reach the field used by tuning and refit selection: "+value;
		}
		for(final String value:new String[]{"-0.01", "1.01", "NaN", "Infinity"}){
			boolean rejected=false;
			try{new BuildNcrnaBoundaryNets(new String[]{config, output, "minprecision="+value});}
			catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("minprecision");}
			assert(rejected) : "The constructor must reject a precision floor outside the finite [0,1] domain: "+value;
		}
		assert(!new File(root, "precision_policy").exists()) : "Constructor policy checks must not launch learning or create output";
		System.out.println("Boundary precision policy passed: explicit [0,1], unchanged 0.95 default, no training");
	}

	/** Only the selected-family mapping namespace is disallowed. */
	private static void policy(){
		BuildNcrnaBoundaryNets.requireRegisteredMapping(new String[]{"callcds=t", "indexk=7", "r58ncrnaindexk=7", "lsuncrnamapk=9"});
		for(String key:new String[]{"ncrnaindexk", "ncrnamapk", "ncrnaindextopn", "ncrnaindexminhits", "ncrnafixedminhits",
			"ncrnaadaptiveminhits", "ncrnaadaptiveshortlist", "ncrnaadaptfloor", "ncrnaadapttopfrac", "ncrnaadaptqfrac",
			"ncrnastrictindexcutoff", "ncrnarankedfallback"}){
			boolean rejected=false;
			try{BuildNcrnaBoundaryNets.requireRegisteredMapping(new String[]{key.toUpperCase(Locale.ROOT)+"=1"});}
			catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("registered family mapping");}
			assert(rejected) : "Each selected mapping alias must fail explicitly before training: "+key;
		}
	}

	private static void check(final File root){
		final File output=new File(root, "product");
		final File work=new File(root, "product.work");
		assert(line(new File(work, "summary.tsv"), "status\tOFF")) : "Huge margins prohibit corrections; zero gain must produce OFF";
		assert(!new File(work, "refit").exists()) : "An ineligible tuning candidate must not start all-data retraining";
		assert(line(new File(work, "tuning/learning.tsv"), "20\t0\t0\t0\t20\t120"))
			: "Twenty complete contexts contribute six candidates each across the train/validation files";
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.resources.root=output.getPath();
		options.boundary=true;
		final NcrnaFamily loaded=NcrnaFamilyLoader.load("sixs_rf01685", options, new FamilyResourceResolver());
		assert(loaded.library.length==19 && loaded.boundary5NetTemplate==null)
			: "The ordinary loader must accept the complete grouped export and its explicit OFF setting";
		System.out.println("Boundary builder finite integration passed: native training, concurrent ordinary callers, no-gain OFF selection and loadable grouped export");
	}

	/** The existing no-movement fixture remains a no-gain control after reuse. */
	private static void checkEvaluate(final File root, final File learned){
		final File output=new File(root, "evaluated");
		final File work=new File(output.getPath()+".work");
		assert(line(new File(work, "summary.tsv"), "status\tTUNING_ONLY_OFF"))
			: "Evaluation-only mode must label its output as tuning evidence, not a completed all-data refit";
		assert(!new File(work, "tuning").exists() && !new File(work, "refit").exists())
			: "Reusing saved networks must never regenerate vectors or invoke tuning/refit trainers";
		assert(line(new File(work, "summary.tsv"), "learned_directory\t"+learned.getAbsolutePath()))
			: "The report must identify the actual borrowed tuning directory";
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.resources.root=output.getPath();
		options.boundary=true;
		final NcrnaFamily loaded=NcrnaFamilyLoader.load("sixs_rf01685", options, new FamilyResourceResolver());
		assert(loaded.library.length==19 && loaded.boundary5NetTemplate==null)
			: "Evaluation-only OFF still exports a complete loadable group with its nineteen selected models";
		System.out.println("Saved boundary evaluation passed: explicit reuse, no training/refit, labeled tuning-only export");
	}

	/** Rejects wrong saved geometry/dimensions and absent products without starting a caller. */
	private static void reuseControls(final File root, final File learned) throws Exception{
		final Method validate=BuildNcrnaBoundaryNets.class.getDeclaredMethod("validateLearned", File.class);
		validate.setAccessible(true);
		final String config="config="+new File(root, "builder.config");
		final String output="out="+new File(root, "reuse_control");
		final String source="learned="+learned;
		final BuildNcrnaBoundaryNets valid=new BuildNcrnaBoundaryNets(new String[]{config, output, "mode=evaluate", source});
		validate.invoke(valid, learned);
		final BuildNcrnaBoundaryNets wrongGeometry=new BuildNcrnaBoundaryNets(new String[]{config, output,
			"mode=evaluate", source, "startinside=4"});
		rejectLearned(validate, wrongGeometry, learned, "geometry differs");
		final File missing=new File(root, "missing_tuning");
		Files.createDirectory(missing.toPath());
		rejectLearned(validate, valid, missing, "Missing readable saved tuning product");
		// CellNet.randomize requires the native cumulative activation distribution.
		if(Function.TYPE_RATES_CUM==null){Function.normalizeTypeRates();}
		for(int inputs:new int[]{9, 10}){
			final File bad=new File(root, "bad_dimensions_"+inputs);
			Files.createDirectory(bad.toPath());
			for(String name:new String[]{"regression.bbnet", "start.tsv", "stop.tsv"}){
				Files.copy(new File(learned, name).toPath(), new File(bad, name).toPath());
			}
			final CellNet net=new CellNet(new int[]{inputs, 3, inputs==9 ? 1 : 2}, 3, 1, 0, 1, new ArrayList<String>());
			net.randomize();
			write(bad, "classifier.bbnet", net.toBytes());
			rejectLearned(validate, valid, bad, inputs==9 ? "inputs; expected" : "one scalar output");
		}
		assert(!new File(root, "reuse_control").exists() && !new File(root, "reuse_control.work").exists())
			: "Saved-input validation controls must not run callers or create a workflow output";
		System.out.println("Saved tuning controls passed: native reads, geometry, input/output dimensions and missing products");
	}

	private static void rejectLearned(final Method validate, final BuildNcrnaBoundaryNets tool,
			final File learned, final String message) throws Exception{
		boolean rejected=false;
		try{validate.invoke(tool, learned);}
		catch(InvocationTargetException expected){
			rejected=expected.getCause() instanceof IllegalArgumentException && expected.getCause().getMessage().contains(message);
		}
		assert(rejected) : "Saved tuning control must fail for its specific invalid contract: "+message;
	}

	private static boolean line(final File path, final String expected){
		boolean found=false;
		try(Input owner=new Input(path.getPath())){
			final ByteFile file=owner.reader;
			for(byte[] row=file.nextLine(); row!=null; row=file.nextLine()){
				if(expected.equals(new String(row, StandardCharsets.UTF_8))){found=true;}
			}
		}
		return found;
	}

	/** Literal outcomes for orchestration tests; no scientific predictions are inferred from them. */
	private static void mockInputs(final File root){
		final String path=new FamilyResourceResolver().resolve("sixs_rf01685", "consensus", null).path;
		final NcrnaSupport.Library library=NcrnaSupport.loadLibrary(path);
		for(int variant=0; variant<2; variant++){
			final ByteBuilder rows=new ByteBuilder("##gff-version 3\n");
			for(int i=0; i<library.names.length; i++){
				rows.append(library.names[i]).append("\tfixture\tRNA\t").append(variant==0 && i==0 ? 2 : 1)
					.tab().append(library.bases[i].length).append("\t.\t+\t.\tfamily:6S-Flavo;model:")
					.append(library.names[i]).nl();
			}
			write(root, variant==0 ? "mock_off.gff" : "mock_good.gff", rows);
		}
	}

	/** Both branches use actual all-data vectors and trainer products; only the calls are controlled. */
	private static void checkRefit(final File root){
		final boolean rejected="reject".equals(System.getenv("NCRNA_BOUNDARY_TEST_MODE"));
		final File output=new File(root, rejected ? "rejected" : "enabled");
		final File work=new File(output.getPath()+".work");
		assert(line(new File(work, "summary.tsv"), "status\t"+(rejected ? "OFF_REFIT_DID_NOT_PASS" : "ENABLED_CLASSIFIER")))
			: "A one-locus tuning gain selects classifier; identical refit gain enables, zero refit gain restores OFF";
		assert(line(new File(work, "refit/learning.tsv"), "0\t120\t40"))
			: "All twenty eligible loci, including both validation TaxIDs, enter the one refit file with two positives each";
		assert(!new File(work, "refit/vectors.validation.tsv.gz").exists())
			: "All-data monitoring must reuse the training file and must not claim another holdout";
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.resources.root=output.getPath();
		options.boundary=true;
		final NcrnaFamily family=NcrnaFamilyLoader.load("sixs_rf01685", options, new FamilyResourceResolver());
		assert((family.boundary5NetTemplate==null)==rejected)
			: "The native loader must consume the selected refit network or the explicit OFF restoration";
		System.out.println("Controlled caller refit branch passed: "+(rejected ? "rejected" : "enabled"));
	}

	/** Invokes only the new parent's concurrency boundary; its public caller route is unchanged.
	 * Reflection is confined to this test so production gains no injection/configuration API.
	 */
	private static void lifecycle(final File root) throws Exception{
		final boolean interrupt="interrupt".equals(System.getenv("NCRNA_BOUNDARY_TEST_MODE"));
		final File output=new File(root, interrupt ? "interrupted" : "failed");
		final File work=new File(output.getPath()+".work");
		final BuildNcrnaBoundaryNets tool=new BuildNcrnaBoundaryNets(new String[]{
			"config="+new File(root, "builder.config"), "out="+output, "t=3"});
		Files.createDirectory(work.toPath());
		final NcrnaCallerEvaluation evaluation=new NcrnaCallerEvaluation(new File(root, "caller.config").getPath(),
			new File(root, "genomes.txt").getPath(), new File(root, "truth.tsv").getPath(), "sixs_rf01685", "fixture");
		final File family=new File(root, "product");
		final Method method=BuildNcrnaBoundaryNets.class.getDeclaredMethod("callTrials", NcrnaCallerEvaluation.class, File[].class);
		method.setAccessible(true);
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final AtomicBoolean retained=new AtomicBoolean();
		final Thread parent=new Thread(()->{
			try{method.invoke(tool, evaluation, new File[]{family, family, family});}
			catch(InvocationTargetException caught){failure.set(caught.getCause());}
			catch(ReflectiveOperationException caught){failure.set(caught);}
			finally{retained.set(Thread.currentThread().isInterrupted());}
		});
		parent.start();
		try(AutoCloseable cleanup=()->{
			if(parent.isAlive()){
				parent.interrupt();
				parent.join(5000);
			}
		}){
			if(interrupt){
				final long deadline=System.nanoTime()+10_000_000_000L;
				while(parent.isAlive() && !allStarted(work) && System.nanoTime()<deadline){Thread.sleep(10);}
				assert(allStarted(work)) : "All three owned callers must announce startup before the test interrupts their parent";
				parent.interrupt();
			}
			parent.join(15000);
			assert(!parent.isAlive() && failure.get() instanceof IllegalStateException && (!interrupt || retained.get()))
				: "Parent failure/interruption must finish caller cleanup, report failure and preserve interruption";
		}
		for(int i=0; i<3; i++){
			final File directory=new File(work, "caller_"+i);
			final ByteFile file=ByteFile.makeByteFile(new File(directory, "mock.pid").getPath(), false);
			final byte[] row=file.nextLine();
			if(file.close() || row==null){throw new IllegalStateException("Missing mock child identity");}
			final String pid=new String(row, StandardCharsets.US_ASCII);
			assert(pid.matches("[0-9]+")) : "Each mock records only its exact own process identity";
			final Process probe=new ProcessBuilder("kill", "-0", pid).redirectError(new File(directory, "probe.err")).start();
			assert(probe.waitFor()!=0) : "Every owned caller must be gone after parent cleanup: "+pid;
		}
		System.out.println("Concurrent caller lifecycle passed: "+(interrupt ? "parent interruption" : "sibling failure"));
	}

	private static boolean allStarted(final File work){
		for(int i=0; i<3; i++){
			if(new File(work, "caller_"+i+"/mock.pid").length()==0){return false;}
		}
		return true;
	}

	private static void write(final File root, final String name, final ByteBuilder text){
		try(Output owner=new Output(new File(root, name).getPath())){
			owner.writer.print(text);
		}
	}
}
