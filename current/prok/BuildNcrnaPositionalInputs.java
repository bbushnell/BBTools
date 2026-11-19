package prok;

import java.io.File;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Builds all-locus member-pooled positional tables and native training matrices.
 * CLI/lifecycle follows template.A_SampleStreamerMT; native streaming and worker
 * ownership live in NcrnaPositionalTables/Vectors. This tool never trains, calls
 * genomes, changes truth or installs resources. Every output is fresh beneath out.
 * All eligible originals teach the tables; monitoring resamples the same loci.
 * @author Brian Bushnell, Raiden
 */
public final class BuildNcrnaPositionalInputs{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, BuildNcrnaPositionalInputs.class, false);
		try{
			final BuildNcrnaPositionalInputs tool=new BuildNcrnaPositionalInputs(pp.args);
			tool.process();
			timer.stop();
			pp.outstream.println("Positional inputs complete; monitoring uses the same eligible loci.");
			pp.outstream.println("Time: "+timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Parses explicit scientific bindings before creating any output directory. */
	public BuildNcrnaPositionalInputs(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}
			if(a.equals("family")){
				family=b;
			}else if(a.equals("familyroot")){
				familyRoot=b;
			}else if(a.equals("ref")){
				ref=b;
			}else if(a.equals("hbm")){
				hbm=b;
			}else if(a.equals("seeds")){
				seeds=b;
			}else if(a.equals("original")){
				original=b;
			}else if(a.equals("graph")){
				graph=b;
			}else if(a.equals("roster")){
				roster=b;
			}else if(a.equals("loci")){
				loci=b;
			}else if(a.equals("exclude")){
				exclude=b;
			}else if(a.equals("refit")){
				refit=Parse.parseBoolean(b);
			}else if(a.equals("histogram")){
				histogram=b;
			}else if(a.equals("evaluationrows5")){
				evaluationRows5=b;
			}else if(a.equals("evaluationrows3")){
				evaluationRows3=b;
			}else if(a.equals("modeloffsets")){
				modelOffsets=b;
			}else if(a.equals("offsetminimum")){
				offsetMinimum=Parse.parseIntKMG(b);
				if(offsetMinimum<=0){throw new IllegalArgumentException("offsetminimum= must be a positive observed-support floor");}
			}else if(a.equals("evaluationcounts")){
				final String[] terms=NcrnaNumbers.value(a, b).split(",");
				if(terms.length!=2){throw new IllegalArgumentException("evaluationcounts= requires the accepted start,stop trial counts");}
				evaluationCounts=new long[]{Parse.parseKMG(terms[0]), Parse.parseKMG(terms[1])};
			}else if(a.equals("radii")){
				final String[] terms=NcrnaNumbers.value(a, b).split(",");
				radii=new int[terms.length];
				for(int i=0; i<terms.length; i++){radii[i]=Parse.parseIntKMG(terms[i]);}
			}else if(a.equals("samples")){
				samples=Parse.parseIntKMG(b);
			}else if(a.equals("trainseed")){
				trainSeed=Parse.parseKMG(NcrnaNumbers.value(a, b));
			}else if(a.equals("monitorseed")){
				monitorSeed=Parse.parseKMG(NcrnaNumbers.value(a, b));
			}else if(a.equals("evaluationseed")){
				evaluationSeed=Parse.parseKMG(NcrnaNumbers.value(a, b));
			}else if(a.equals("interchange")){
				interchange=Parse.parseBoolean(b);
			}else if(a.equals("fullerrors")){
				fullErrors=Parse.parseBoolean(b);
			}else if(a.equals("prior")){
				prior=Double.parseDouble(b);
			}else if(a.equals("tripletables")){
				tripleTables=Parse.parseBoolean(b);
			}else if(a.equals("startinside")){
				startInside=NcrnaNumbers.intValue(a, b);
			}else if(a.equals("startoutside")){
				startOutside=NcrnaNumbers.intValue(a, b);
			}else if(a.equals("stopinside")){
				stopInside=NcrnaNumbers.intValue(a, b);
			}else if(a.equals("stopoutside")){
				stopOutside=NcrnaNumbers.intValue(a, b);
			}else if(!parser.parse(arg, a, b)){
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		in=parser.in1;
		out=parser.out1;
		workers=Math.max(1, parser.workers()>0 ? parser.workers() : Shared.threads());
		if(in==null || out==null || family==null || familyRoot==null || ref==null || hbm==null || seeds==null
				|| original==null || graph==null || roster==null || loci==null || histogram==null || radii==null
				|| radii.length==0 || samples<1 || trainSeed<0 || monitorSeed<0 || trainSeed==monitorSeed
				|| !(prior>0 && Double.isFinite(prior))){
			throw new IllegalArgumentException("Require in/out/family/familyroot/ref/hbm/seeds/original/graph/roster/loci/histogram/radii/samples/trainseed/monitorseed/prior; see launcher help");
		}
		final int fiveK=NcrnaPositionalTables.wordLength(startInside, startOutside);
		final int threeK=NcrnaPositionalTables.wordLength(stopInside, stopOutside);
		if(tripleTables && (fiveK!=9 || threeK!=9)){
			throw new IllegalArgumentException("tripletables=t requires K9 at both ends; it adds training LOO and HBM coverage features");
		}
		if(interchange && (!tripleTables || radii.length!=1 || radii[0]!=10
				|| !refit && (exclude==null || evaluationSeed<0 || evaluationSeed==trainSeed || evaluationSeed==monitorSeed))){
			throw new IllegalArgumentException("interchange=t requires tripletables=t, radii=10, exclude= and three distinct nonnegative draw seeds");
		}
		if(!interchange && evaluationSeed>=0){throw new IllegalArgumentException("evaluationseed= requires explicit interchange=t");}
		if(fullErrors && !interchange){throw new IllegalArgumentException("fullerrors=t requires explicit interchange=t and its source-excluded 69/22 schema");}
		final boolean fixed=evaluationRows5!=null || evaluationRows3!=null || modelOffsets!=null || evaluationCounts!=null || offsetMinimum>=0;
		if(refit && (!fullErrors || exclude!=null || fixed || evaluationSeed>=0)){
			throw new IllegalArgumentException("All-data refit requires fullerrors=t and no excluded loci, evaluation seed or fixed evaluation trials");
		}
		if(fixed && (!fullErrors || evaluationRows5==null || evaluationRows3==null || modelOffsets==null
				|| evaluationCounts==null || evaluationCounts[0]<=0 || evaluationCounts[1]<=0 || offsetMinimum<=0)){
			throw new IllegalArgumentException("Fixed centered outer trials require fullerrors=t, both evaluationrows, modeloffsets, offsetminimum and exact evaluationcounts");
		}
		if(!(prior/(1<<(2*fiveK))>0 && prior/(1<<(2*threeK))>0)){
			throw new IllegalArgumentException("Total prior must remain positive after per-end division by4^K");
		}
		if(parser.in2!=null || parser.out2!=null || parser.maxReads>=0 || parser.samplerate!=1
				|| parser.overwrite || parser.append){
			throw new IllegalArgumentException("Read every unpaired original context into a fresh output directory; sampling, limits, overwrite and append are unsupported");
		}
		for(int i=0; i<radii.length; i++){
			if(radii[i]<1 || radii[i]>(Integer.MAX_VALUE-4)/2 || (i>0 && radii[i]<=radii[i-1])){
				throw new IllegalArgumentException("Radii must be positive, strictly increasing base counts with representable feature widths");
			}
		}
		entry=NcrnaFamilyRegistry.get(family);
		// Ungrouped registrations use an empty member; their stable ID names the table.
		member=entry.group==null ? entry.id : entry.member;
		root=new File(out).getAbsoluteFile();
		if(root.exists() || root.getParentFile()==null || !root.getParentFile().isDirectory() || !new File(familyRoot).isDirectory()
				|| !Tools.testInputFiles(false, true, in, ref, hbm, seeds, original, graph, roster, loci, histogram, exclude,
					evaluationRows5, evaluationRows3, modelOffsets)){
			throw new IllegalArgumentException("All inputs must be readable; out must be absent beneath an existing parent and familyroot must be a directory");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	/** Resolves the ordinary working library once, then counts and vectorizes all radii. */
	private void process(){
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.resources.consensus=ref;
		options.resources.models=hbm;
		options.resources.seeds=seeds;
		options.resources.root=familyRoot;
		options.boundaryOverride=false;
		options.endpoints=false;
		options.topN=entry.topN;
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		resolver.bindRoot(family, familyRoot);
		final NcrnaFamily accepted=NcrnaFamilyLoader.load(family, options, resolver);
		final NcrnaLearningPopulation population=NcrnaLearningPopulation.loadForRefit(NcrnaSupportIO.read(graph),
			NcrnaSupport.loadLibrary(original).names, roster, loci, exclude);
		final NcrnaFixedEndpointTrials fixed=modelOffsets==null ? null : NcrnaFixedEndpointTrials.load(population,
			accepted, member, evaluationRows5, evaluationRows3, modelOffsets, offsetMinimum, samples, evaluationCounts);
		final NcrnaSignedEndpointErrors errors=fullErrors ? NcrnaSignedEndpointErrors.load(histogram, member) : null;
		final NcrnaBoundaryJitter[] jitter=fullErrors ? null : new NcrnaBoundaryJitter[radii.length];
		if(jitter!=null){
			for(int r=0; r<radii.length; r++){jitter[r]=NcrnaBoundaryJitter.load(histogram, member, radii[r]);}
		}
		if(!root.mkdir()){throw new IllegalArgumentException("Could not create fresh positional output directory: "+root);}
		boolean available=true;
		if(errors!=null){errors.writeBands(path("full_signed_errors.bands.tsv"));}
		else{
			for(int r=0; r<radii.length; r++){
				jitter[r].write(path("r"+radii[r]+".jitter.tsv"));
				available&=jitter[r].available(0) && jitter[r].available(1);
			}
		}
		if(!available){throw new IllegalArgumentException("An observed conditioning window is empty; accounting is saved, but no fallback data will be invented");}
		final NcrnaPositionalTables tables=fullErrors
			? NcrnaPositionalTables.buildFull(population, in, member, radii, startInside, startOutside, stopInside, stopOutside, prior, errors)
			: tripleTables
			? NcrnaPositionalTables.buildTriple(population, in, member, radii, startInside, startOutside, stopInside, stopOutside, prior)
			: NcrnaPositionalTables.build(population, in, member, radii, startInside, startOutside, stopInside, stopOutside, prior);
		for(int r=0; r<radii.length; r++){
			if(fullErrors){
				for(int probe=0; probe<3; probe++){
					tables.fiveWeighted[r][probe].writeScores(path("r"+radii[r]+".5prime.probe"+probe+".weighted.tsv"));
					tables.threeWeighted[r][probe].writeScores(path("r"+radii[r]+".3prime.probe"+probe+".weighted.tsv"));
				}
				continue;
			}
			tables.five[r].write(path("r"+radii[r]+".5prime.table.tsv"));
			tables.three[r].write(path("r"+radii[r]+".3prime.table.tsv"));
			if(tripleTables){
				for(int probe=1; probe<3; probe++){
					final String side=probe==1 ? "left" : "right";
					tables.fiveProbes[r][probe].write(path("r"+radii[r]+".5prime."+side+".table.tsv"));
					tables.threeProbes[r][probe].write(path("r"+radii[r]+".3prime."+side+".table.tsv"));
				}
			}
		}
		final NcrnaPositionalVectors vectors=refit
			? NcrnaPositionalVectors.writeRefit(population, in, accepted, tables, samples, trainSeed, monitorSeed, workers, root)
			: interchange
			? NcrnaPositionalVectors.writeInterchange(population, in, accepted, tables, jitter,
				samples, trainSeed, monitorSeed, evaluationSeed, workers, root, fixed)
			: NcrnaPositionalVectors.write(population, in, accepted, tables, jitter,
				samples, trainSeed, monitorSeed, workers, root);
		assert(vectors.records==tables.records && vectors.filtered==tables.filtered)
			: "Count and vector passes must conserve the same complete original axis and exclusions";
		try(Output output=new Output(path("population.tsv"))){
			output.writer.println("member\toriginals\tfiltered\teligible\ttable_loci\tno_geometry\tusable_context\tmissing_gc\tmissing_flanks\tshort_context\tvector_loci\tno_model");
			output.writer.print(new ByteBuilder().append(member).tab().append(tables.records).tab().append(tables.filtered)
				.tab().append(tables.learning).tab().append(tables.counted).tab().append(tables.noGeometry)
				.tab().append(tables.availability[NcrnaLearningContext.AVAILABLE]).tab().append(tables.availability[NcrnaLearningContext.NO_GC])
				.tab().append(tables.availability[NcrnaLearningContext.NO_FLANK]).tab().append(tables.availability[NcrnaLearningContext.MALFORMED_SPAN])
				.tab().append(vectors.processed).tab().append(vectors.noModel).nl());
		}
	}

	private String path(final String name){return new File(root, name).getPath();}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in, out, family, familyRoot, ref, hbm, seeds, original, graph, roster, loci, histogram, exclude;
	private boolean refit;
	/** Optional frozen outer sidecars and supported offsets; absence preserves the original draw path. */
	private String evaluationRows5, evaluationRows3, modelOffsets;
	private long[] evaluationCounts;
	private int offsetMinimum=-1;
	private String member;
	private File root;
	private NcrnaFamilyRegistry.Entry entry;
	private int[] radii;
	private int samples=-1, workers;
	/** Semantic word placement; defaults preserve the original fully-inside K9 recipe. */
	private int startInside=9, startOutside=0, stopInside=9, stopOutside=0;
	private long trainSeed=-1, monitorSeed=-1, evaluationSeed=-1;
	/** Required positive CLI value; unused zero fails the existing pre-output validation. */
	private double prior;
	/** Optional three-table + leave-one-locus-out + HBM-coverage experiment. */
	private boolean tripleTables;
	/** Plain numeric files plus separate row metadata and held-out evaluation. */
	private boolean interchange;
	/** Opt-in full signed draws, weighted background and an explicit 22nd target. */
	private boolean fullErrors;
}
