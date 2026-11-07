package prok;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.LongHashSet;
import ml.CellNet;
import ml.DataLoader;
import ml.Matrix;
import ml.RegressionTrainer;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.SequenceInput;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Bounded real-contig measurement of boundary-training and inference GC.
 * Uses the Streamer template's native sequence input and ByteFile/LineParser1
 * metadata input. Sampling reads existing contexts, never all source genomes.
 * The original population/split stays intact; exclusions restrict actual vector
 * generation to the selected loci. Diagnostic refit means all selected TRAIN
 * loci, not a release refit. No network is trained and no inference code changes.
 * Recording networks capture pre-normalization inputs, not decision sensitivity.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaGcComparison{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaGcComparison.class, false);
		try{
			final NcrnaGcComparison tool=new NcrnaGcComparison(pp.args);
			tool.process();
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** All scientific routes are explicit; the source-reading sample is bounded by fixed bins. */
	public NcrnaGcComparison(final String[] args){
		final Parser parser=new Parser();
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("graph")){graphPath=b;}
			else if(a.equals("original")){original=b;}
			else if(a.equals("roster")){roster=b;}
			else if(a.equals("loci")){loci=b;}
			else if(a.equals("ref")){ref=b;}
			else if(a.equals("hbm")){hbm=b;}
			else if(a.equals("seeds")){seeds=b;}
			else if(a.equals("family")){familyName=b;}
			else if(a.equals("splitseed")){seed=NcrnaNumbers.longValue(a, b);}
			else if(a.equals("requireexact")){requireExact=Parse.parseBoolean(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		contexts=parser.in1;
		if(contexts==null || parser.out1==null || seed<0 || familyName==null || parser.in2!=null
				|| parser.append || parser.maxReads>=0 || parser.samplerate!=1){
			throw new IllegalArgumentException("Require complete in=contexts, out=fresh_directory, family= and nonnegative splitseed=");
		}
		final String[] paths={contexts, graphPath, original, roster, loci, ref, hbm, seeds};
		for(final String path:paths){if(path==null){throw new IllegalArgumentException("Supply graph/original/roster/loci/ref/hbm/seeds");}}
		if(!Tools.testInputFiles(false, true, paths)){throw new IllegalArgumentException("Unreadable comparison input");}
		root=new File(parser.out1);
		if(root.exists() || !root.getAbsoluteFile().getParentFile().isDirectory()){
			throw new IllegalArgumentException("out= must be fresh beneath an existing parent");
		}
		familyName=CallGenes.parseNcrnaFamily(familyName);
		NcrnaFamilyRegistry.sharedBoundaryMember(familyName);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Selection          ----------------*/
	/*--------------------------------------------------------------*/

	private void process() throws Exception{
		if(Tools.FORCE_JAVA_PARSE_DOUBLE){throw new IllegalArgumentException("Measure the ordinary native decimal reader without a global parser override");}
		Shared.setThreads(1);
		ReadWrite.setZipThreads(1);
		if(!root.mkdir()){throw new IllegalArgumentException("Could not create "+root);}
		final NcrnaSupport.Graph graph=NcrnaSupportIO.read(graphPath);
		final String[] names=NcrnaSupport.loadLibrary(original).names;
		final NcrnaLearningPopulation originalPopulation=NcrnaLearningPopulation.load(graph, names, roster, loci, null, seed);
		loadCoordinates(originalPopulation);
		select(originalPopulation);
		writeSelection(originalPopulation);
		final String exclusions=file("exclusions.txt");
		final boolean[] selected=new boolean[originalPopulation.ids.length];
		for(final Sample sample:samples){selected[sample.gene]=true;}
		try(Output output=new Output(exclusions)){
			for(int gene=0; gene<selected.length; gene++){
				if(!selected[gene]){output.writer.println(originalPopulation.ids[gene]);}
			}
		}
		population=NcrnaLearningPopulation.load(graph, names, roster, loci, exclusions, seed);
		check(Arrays.equals(originalPopulation.split().validationByGene, population.split().validationByGene),
			"Diagnostic exclusions must not change original whole-TaxID split membership");
		loadFamily();
		tables=NcrnaLearningTables.build(population, contexts, 7, 0, 7, 0, true, 1);
		final NcrnaLearningVectors vectors=NcrnaLearningVectors.write(population, contexts, family, tables, 1, 1, 1, 1, file("vectors"));
		check(vectors.validationPath==null && vectors.trainPaths.length==1
			&& vectors.positives[0]==2*vectors.processed && vectors.rows[0]>vectors.positives[0],
			"Native classifier loading requires positive and negative labels; radius one retains both true ends and nearby candidates");
		readMatrices(vectors);
		readSources();
		writeResults(vectors);
		System.out.println("NCRNA_GC_COMPARISON_COMPLETE selected="+samples.length+" matrix_loci="+vectors.processed
			+" source_genomes="+sourceGenomes+" no_training=true no_tolerance_decision=true");
	}

	/** Stores primitive coordinates on the existing original-query axis. */
	private void loadCoordinates(final NcrnaLearningPopulation p){
		lengths=new int[p.ids.length];
		starts=new int[p.ids.length];
		stops=new int[p.ids.length];
		minus=new boolean[p.ids.length];
		final LineParser1 row=new LineParser1('\t');
		try(Input input=new Input(loci)){
			header(input, NcrnaLearningPopulation.LOCUS_HEADER);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				row.set(line);
				check(row.terms()==21, "Expected corrected native 21-column locus metadata");
				final int gene=p.index(row.parseString(0));
				check(gene>=0 && lengths[gene]==0, "Unknown or duplicate original locus");
				lengths[gene]=number(row, 15);
				minus[gene]=row.termEquals("-", 7);
				check(minus[gene] || row.termEquals("+", 7), "Explicit source strand required");
				row.setBounds(8);
				int separator=row.a();
				while(separator<row.b() && line[separator]!='-'){separator++;}
				check(separator<row.b(), "Corrected diagnostic loci require one contiguous interval");
				starts[gene]=NcrnaNumbers.unsignedInt(line, row.a(), separator)-1;
				stops[gene]=NcrnaNumbers.unsignedInt(line, separator+1, row.b())-1;
				check(starts[gene]>=0 && stops[gene]>=starts[gene] && stops[gene]<lengths[gene]
					&& stops[gene]-starts[gene]+1==p.length[gene], "Corrected coordinates must fit the declared contig and core length");
			}
		}
		for(final int length:lengths){check(length>0, "Every original query requires its real source-contig length");}
	}

	/** One lexical source/contig/locus representative per GC-decile/length bin. */
	private void select(final NcrnaLearningPopulation p){
		final Sample[] bins=new Sample[BINS];
		final boolean[] seen=new boolean[p.ids.length];
		final NcrnaLearningContext context=new NcrnaLearningContext();
		final boolean trim=Shared.TRIM_READ_DESCRIPTION, upper=Read.TO_UPPER_CASE;
		Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=true;
		int order=0;
		try(SequenceInput owner=sequenceInput(contexts)){
			final Streamer input=owner.stream;
			try{
				input.start();
				for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
					for(final Read read:batch){
						context.set(read, p);
						final int gene=context.gene;
						check(!seen[gene], "Original context duplicated: "+p.ids[gene]);
						seen[gene]=true;
						final int position=order++;
						if(!p.eligible(gene) || !context.hasFeatures()){
							unavailable++;
							continue;
						}
						final int bin=lengthBin(lengths[gene])*10+Math.min(9, (int)(context.gc*10));
						binCounts[bin]++;
						final String source=p.sourceName(gene);
						final Sample prior=bins[bin];
						if(prior!=null && source.compareTo(prior.source)>0){continue;}
						final String contig=RrnaResourceIO.field(read.id, "contig");
						check(contig!=null && !contig.isEmpty(), "Selected context requires original contig metadata");
						if(prior==null || source.compareTo(prior.source)<0 || contig.compareTo(prior.contig)<0
							|| contig.equals(prior.contig) && p.ids[gene].compareTo(prior.id)<0){
							bins[bin]=new Sample(gene, position, bin, p.ids[gene], source, contig,
								new Read(read.bases.clone(), null, read.id, position, false), context.gc);
						}
					}
				}
			}catch(Throwable failure){KillSwitch.exceptionKill(failure);}
		}finally{
			Shared.TRIM_READ_DESCRIPTION=trim;
			Read.TO_UPPER_CASE=upper;
		}
		check(order==seen.length, "Sampling must account for every original context");
		final ArrayList<Sample> selected=new ArrayList<Sample>();
		for(final Sample sample:bins){if(sample!=null){selected.add(sample);}}
		check(!selected.isEmpty(), "No eligible real source contexts for the diagnostic");
		samples=selected.toArray(new Sample[selected.size()]);
		Arrays.sort(samples, (a, b)->Integer.compare(a.order, b.order));
		for(int i=0; i<samples.length; i++){
			for(int j=0; j<i; j++){
				check(!samples[i].source.equals(samples[j].source) || !samples[i].contig.equals(samples[j].contig),
					"The comparison denominator requires distinct source contigs");
			}
		}
	}

	/** Fixed source-length bins include the integer-to-float exactness boundary. */
	private static int lengthBin(final int length){
		check(length>0, "Contig length must be positive");
		for(int bin=0; bin<LENGTH_LIMITS.length; bin++){if(length<=LENGTH_LIMITS[bin]){return bin;}}
		return LENGTH_LIMITS.length;
	}

	/** Records empty bins and joins selected sources to their original FASTA paths. */
	private void writeSelection(final NcrnaLearningPopulation p){
		try(Output output=new Output(file("bins.tsv"))){
			output.writer.print("length_bin\tgc_decile\teligible_loci\tselected_locus\n");
			final ByteBuilder row=new ByteBuilder();
			for(int bin=0; bin<BINS; bin++){
				String id="NA";
				for(final Sample sample:samples){if(sample.bin==bin){id=sample.id;}}
				output.writer.print(row.clear().append(LENGTH_LABELS[bin/10]).tab().append(bin%10)
					.tab().append(binCounts[bin]).tab().append(id).nl());
			}
		}
		final LineParser1 row=new LineParser1('\t');
		try(Input input=new Input(roster)){
			header(input, NcrnaLearningPopulation.ROSTER_HEADER);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				row.set(line);
				check(row.terms()==11, "Expected native eleven-column source roster");
				final int source=number(row, 0);
				for(final Sample sample:samples){
					if(p.source[sample.gene]==source){sample.fasta=row.parseString(7);}
				}
			}
		}
		for(final Sample sample:samples){check(sample.fasta!=null, "Selected source must join its literal genome path");}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Production Paths      ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads the accepted library without enabling any trained boundary network. */
	private void loadFamily(){
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.resources.consensus=ref;
		options.resources.models=hbm;
		options.resources.seeds=seeds;
		options.boundaryOverride=false;
		options.endpoints=false;
		options.topN=NcrnaFamilyRegistry.get(familyName).topN;
		family=NcrnaFamilyLoader.load(familyName, options, new FamilyResourceResolver());
	}

	/** Loads the unchanged writer's bytes through both real native trainer readers. */
	private void readMatrices(final NcrnaLearningVectors vectors) throws Exception{
		final String path=vectors.trainPaths[0];
		final Matrix classifier=DataLoader.load(path, Integer.MAX_VALUE, false, 0, 0, false, 0)[0].matrix;
		final File unused=new File(root, "NO_TRAINING_MODEL.bbnet");
		final RegressionTrainer regression=new RegressionTrainer(new String[]{"in="+path, "out="+unused,
			"dims=10,1", "epochs=1", "vfraction=0", "threads=1", "simd=f"});
		final Method load=RegressionTrainer.class.getDeclaredMethod("loadData");
		load.setAccessible(true);
		load.invoke(regression);
		final Field dataField=RegressionTrainer.class.getDeclaredField("trainingData");
		dataField.setAccessible(true);
		final Object data=dataField.get(regression);
		final Field inputField=data.getClass().getDeclaredField("inputs");
		inputField.setAccessible(true);
		final float[][] regressionInputs=(float[][])inputField.get(data);
		final Field targetField=data.getClass().getDeclaredField("targets");
		targetField.setAccessible(true);
		final float[][] regressionTargets=(float[][])targetField.get(data);
		check(!unused.exists() && !Tools.FORCE_JAVA_PARSE_DOUBLE, "Only native reading may run; no training or parser-policy switch");
		check(classifier.numInputs==10 && classifier.numOutputs==1 && classifier.inputs.length==vectors.rows[0]
			&& regressionInputs.length==classifier.inputs.length && regressionTargets.length==regressionInputs.length,
			"Both readers must retain every emitted ten-feature row and target");
		final NcrnaBoundaryExamples examples=new NcrnaBoundaryExamples(family, tables.start, tables.stop);
		final NcrnaLearningContext context=new NcrnaLearningContext();
		final NcrnaEndpointTable.LeaveOneOut startCounts=tables.start.leaveOneOut(), stopCounts=tables.stop.leaveOneOut();
		final Capture capture=new Capture();
		final Recorder recorder=new Recorder(capture);
		int row=0;
		for(final Sample sample:samples){
			context.set(sample.context, population);
			if(!examples.prepare(context, true)){continue;}
			sample.model=examples.model();
			startCounts.set(context.read.bases, context.start);
			stopCounts.set(context.read.bases, context.stop);
			final int first=row;
			for(int end=0; end<2; end++){
				for(int offset=-1; offset<=1; offset++){
					if(examples.fill(end==1, offset)!=NcrnaBoundaryExamples.VALID){continue;}
					check(row<classifier.inputs.length && classifier.inputs[row][4]==end && regressionInputs[row][4]==end
						&& classifier.outputs[row][0]==examples.label && regressionTargets[row][0]==examples.label,
						"One worker/shard preserves locus order and each native start/stop candidate and label");
					if(row==first){
						sample.classifierGc=classifier.inputs[row][9];
						sample.regressionGc=regressionInputs[row][9];
					}
					check(bits(sample.classifierGc)==bits(classifier.inputs[row][9])
						&& bits(sample.regressionGc)==bits(regressionInputs[row][9]),
						"Every candidate for both ends of one source locus must carry the same GC");
					capture.reset();
					final int left=context.start+(end==0 ? offset : 0), right=context.stop+(end==1 ? offset : 0);
					final NcrnaEndpointTable table=end==0 ? tables.start : tables.stop;
					NcrnaBoundaryScorer.score(recorder, context.read.bases, left, right, end==1,
						family.library[sample.model], family.models[sample.model], end==0 ? startCounts : stopCounts,
						table.inside, table.outside, context.gc, examples.upperMedianLength);
					check(capture.count==1, "Each explicit scorer candidate must yield one raw network input vector");
					compareFeatures(examples.values, capture.values, classifier.inputs[row], regressionInputs[row], sample.id, row);
					row++;
				}
			}
			check(row-first>=2, "Each prepared locus must retain its two true-boundary rows");
		}
		check(row==classifier.inputs.length, "Missing-model loci and clipped candidates must remain explicit in matrix accounting");
	}

	/** Checks serialization against raw features and inference scoring with the same LOO counts.
	 * Ordinary inference uses full counts; this isolates arithmetic and serialization
	 * without claiming those intentionally different count populations are identical.
	 */
	private void compareFeatures(final float[] raw, final float[] scorer, final float[] classifier,
			final float[] regression, final String locus, final int row){
		check(raw.length==10 && scorer.length==10 && classifier.length==10 && regression.length==10,
			"All parity participants must use the same ten-feature layout");
		for(int feature=0; feature<10; feature++){
			final boolean scorerMatch=bits(raw[feature])==bits(scorer[feature]);
			final boolean classifierMatch=bits(raw[feature])==bits(classifier[feature]);
			final boolean regressionMatch=bits(raw[feature])==bits(regression[feature]);
			if(scorerMatch){scorerExact[feature]++;}
			if(classifierMatch){classifierExact[feature]++;}
			if(regressionMatch){regressionExact[feature]++;}
			if(requireExact){
				check(scorerMatch && classifierMatch && regressionMatch,
					"Boundary feature roundtrip mismatch: locus="+locus+", row="+row+", feature="+feature
					+", raw="+raw[feature]+", scorer="+scorer[feature]+", classifier="+classifier[feature]+", regression="+regression[feature]);
			}
		}
		featureRows++;
	}

	/** Reads each selected original genome once; only chosen contigs enter measurement. */
	private void readSources() throws Exception{
		final Capture capture=new Capture();
		final Recorder recorder=new Recorder(capture);
		final NcrnaScavenger scavenger=new NcrnaScavenger(family.library, family.models, family.modelNames,
			new LongHashSet(), 17, 16, 20, family.indexK, family.indexTopN, family.adaptive,
			family.adaptFloor, family.adaptTopFrac, family.adaptQFrac, family.fixedMinHits,
			0, 1, .70f, .70f, recorder, recorder, tables.start, tables.stop, 7, 0, 7, 0,
			PackagedBoundaryConfig.upperMedianLength(family.library), new int[]{0}, new int[]{0});
		final Method gc=NcrnaScavenger.class.getDeclaredMethod("contigGC", byte[].class);
		final Method refine=NcrnaScavenger.class.getDeclaredMethod("refineBoundaryNN", Orf.class, byte[].class,
			int.class, float.class, boolean.class, int.class);
		gc.setAccessible(true);
		refine.setAccessible(true);
		final map.ObjectSet<String> genomes=new map.ObjectSet<String>(String.class);
		final boolean trim=Shared.TRIM_READ_DESCRIPTION, upper=Read.TO_UPPER_CASE;
		Shared.TRIM_READ_DESCRIPTION=true;
		Read.TO_UPPER_CASE=true;
		try{
			for(final Sample first:samples){
				if(!genomes.add(first.fasta)){continue;}
				sourceGenomes++;
				try(SequenceInput owner=sequenceInput(first.fasta)){
					final Streamer input=owner.stream;
					try{
						input.start();
						for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
							for(final Read read:batch){
								for(final Sample sample:samples){
									if(sample.fasta.equals(first.fasta) && sample.contig.equals(read.id)){
										measure(sample, read.bases, scavenger, gc, refine, capture);
									}
								}
							}
						}
					}catch(Throwable failure){KillSwitch.exceptionKill(failure);}
				}
			}
		}finally{
			Shared.TRIM_READ_DESCRIPTION=trim;
			Read.TO_UPPER_CASE=upper;
		}
		for(final Sample sample:samples){check(sample.measured, "Selected original contig not found: "+sample.source+"/"+sample.contig);}
	}

	/** Checks the source join, then observes the actual private refinement entry point.
	 * Source coordinates are converted to the sense-oriented full contig, as in calling.
	 * Only GC is compared; sampled LOO enrichment is not an all-feature parity claim.
	 */
	private void measure(final Sample sample, final byte[] originalBases, final NcrnaScavenger scavenger,
			final Method gc, final Method refine, final Capture capture) throws Exception{
		check(!sample.measured && originalBases.length==lengths[sample.gene], "Duplicate source record or changed contig length");
		final byte[] bases=minus[sample.gene] ? originalBases.clone() : originalBases;
		if(minus[sample.gene]){dna.AminoAcid.reverseComplementBasesInPlace(bases);}
		final int start=minus[sample.gene] ? bases.length-1-stops[sample.gene] : starts[sample.gene];
		final int stop=minus[sample.gene] ? bases.length-1-starts[sample.gene] : stops[sample.gene];
		final int flank=population.leftFlank[sample.gene];
		for(int i=0; i<=stop-start; i++){
			check(bases[start+i]==sample.context.bases[flank+i], "Original source core changed: "+sample.id);
		}
		sample.geneGc=GeneCaller.calcContigGC(bases);
		check(bits(sample.geneGc)==bits(sample.headerGc), "Existing context GC differs from recomputed GeneCaller source GC: "+sample.id);
		sample.scavengerGc=(Float)gc.invoke(scavenger, (Object)bases);
		check(Float.isFinite(sample.geneGc) && Float.isFinite(sample.scavengerGc), "Selected real source must have defined GC for ULP comparisons");
		if(requireExact){check(bits(sample.geneGc)==bits(sample.scavengerGc), "Sample source GC must match actual Scavenger input for exact whole-path GC parity");}
		if(sample.model>=0){
			capture.reset();
			final Orf orf=new Orf(sample.contig, start, stop, minus[sample.gene] ? Shared.MINUS : Shared.PLUS,
				0, bases, false, ProkObject.RNA);
			refine.invoke(scavenger, orf, bases, sample.model, 1f, true, 0);
			check(capture.count>0 && bits(capture.gc)==bits(sample.scavengerGc), "Actual Scavenger refinement must pass its source GC to the recording net");
			sample.netGc=capture.gc;
			sample.captures=capture.count;
		}
		sample.measured=true;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Reports          ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes complete denominators and differences without selecting a tolerance. */
	private void writeResults(final NcrnaLearningVectors vectors){
		final Stats[] stats={new Stats("helper"), new Stats("classifier_serialization"), new Stats("classifier_combined"),
			new Stats("regression_serialization"), new Stats("regression_combined"), new Stats("reader_agreement")};
		try(Output output=new Output(file("loci.tsv"))){
			output.writer.print("locus\tsource\tcontig\tstrand\tstart1\tstop1\tcontig_length\tsource_fasta\tmodel_status\tnet_captures\tGeneCaller_GC\tScavenger_GC\tnet_GC\tclassifier_GC\tregression_GC\n");
			final ByteBuilder row=new ByteBuilder();
			for(final Sample sample:samples){
				row.clear().append(sample.id).tab().append(sample.source).tab().append(sample.contig).tab()
					.append(minus[sample.gene] ? '-' : '+').tab().append(starts[sample.gene]+1).tab().append(stops[sample.gene]+1)
					.tab().append(lengths[sample.gene]).tab().append(sample.fasta).tab().append(sample.model<0 ? "NO_MODEL" : "MATRIX_AND_NET")
					.tab().append(sample.captures);
				for(final float value:new float[]{sample.geneGc, sample.scavengerGc, sample.netGc, sample.classifierGc, sample.regressionGc}){
					row.tab();
					if(Float.isNaN(value)){row.append("NA");}else{row.append((double)value, 18);}
				}
				output.writer.print(row.nl());
				stats[0].add(sample.geneGc, sample.scavengerGc);
				if(sample.model>=0){
					stats[1].add(sample.classifierGc, sample.geneGc);
					stats[2].add(sample.classifierGc, sample.netGc);
					stats[3].add(sample.regressionGc, sample.geneGc);
					stats[4].add(sample.regressionGc, sample.netGc);
					stats[5].add(sample.classifierGc, sample.regressionGc);
				}
			}
		}
		try(Output output=new Output(file("differences.tsv"))){
			output.writer.print("comparison\tcontigs\texact_bits\tabs_p50\tabs_p95\tabs_max\tulp_p50\tulp_p95\tulp_max\n");
			for(final Stats stat:stats){stat.write(output);}
		}
		try(Output output=new Output(file("accounting.tsv"))){
			output.writer.print("original_contexts\tselected_contigs\tsource_genomes\tselection_unavailable\tfiltered\tvector_unavailable\tno_model\tprocessed\trows\n");
			output.writer.print(new ByteBuilder().append(vectors.records).tab().append(samples.length).tab().append(sourceGenomes)
				.tab().append(unavailable).tab().append(vectors.filtered).tab().append(vectors.unavailable).tab()
				.append(vectors.noModel).tab().append(vectors.processed).tab().append(vectors.rows[0]).nl());
		}
		try(Output output=new Output(file("features.tsv"))){
			output.writer.print("feature\trows\tscorer_matching_LOO_exact\tclassifier_reader_exact\tregression_reader_exact\n");
			final ByteBuilder row=new ByteBuilder();
			for(int feature=0; feature<10; feature++){
				output.writer.print(row.clear().append(FEATURE_NAMES[feature]).tab().append(featureRows).tab()
					.append(scorerExact[feature]).tab().append(classifierExact[feature]).tab().append(regressionExact[feature]).nl());
			}
		}
	}

	/** Float-step distance is defined here only for finite nonnegative GC fractions. */
	private static final class Stats{
		Stats(final String name_){name=name_;}
		void add(final float a, final float b){
			check(a>=0 && a<=1 && b>=0 && b<=1, "GC comparison requires finite fractions in [0,1]");
			check(size<BINS, "At most one selected contig per declared bin");
			errors[size]=Math.abs((double)a-b);
			ulps[size]=Math.abs((long)bits(a)-bits(b));
			if(bits(a)==bits(b)){exact++;}
			size++;
		}
		void write(final Output output){
			Arrays.sort(errors, 0, size);
			Arrays.sort(ulps, 0, size);
			final ByteBuilder row=new ByteBuilder(name).tab().append(size).tab().append(exact);
			if(size==0){row.append("\tNA\tNA\tNA\tNA\tNA\tNA");}
			else{
				final int median=(size*50+99)/100-1, p95=(size*95+99)/100-1;
				row.tab().appendSci(errors[median], 12).tab().appendSci(errors[p95], 12).tab().appendSci(errors[size-1], 12)
					.tab().append(ulps[median]).tab().append(ulps[p95]).tab().append(ulps[size-1]);
			}
			output.writer.print(row.nl());
		}
		final String name;
		final double[] errors=new double[BINS];
		final long[] ulps=new long[BINS];
		int size, exact;
	}

	/** Network instrumentation only: no learned outputs, normalization or parameter changes. */
	private static final class Recorder extends CellNet{
		Recorder(final Capture capture_){capture=capture_;}
		@Override
		public synchronized CellNet copy(final boolean copyDelta){return new Recorder(capture);}
		@Override
		public CellNet applyInput(final float[] values){
			NcrnaGcComparison.check(values.length==10, "Expected the existing ten-feature boundary input");
			capture.add(values);
			return this;
		}
		@Override
		public float feedForward(){return .5f;}
		final Capture capture;
	}

	private static final class Capture{
		void reset(){
			count=0;
			gc=Float.NaN;
		}
		void add(final float[] input){
			if(count>0){check(bits(gc)==bits(input[9]), "All endpoint evaluations for one source contig must receive identical GC");}
			System.arraycopy(input, 0, values, 0, values.length);
			gc=input[9];
			count++;
		}
		final float[] values=new float[10];
		float gc;
		int count;
	}

	/** Only selected records retain their small context bases; source genomes are streamed. */
	private static final class Sample{
		Sample(final int gene_, final int order_, final int bin_, final String id_, final String source_,
				final String contig_, final Read context_, final float gc_){
			gene=gene_;
			order=order_;
			bin=bin_;
			id=id_;
			source=source_;
			contig=contig_;
			context=context_;
			headerGc=gc_;
		}
		final int gene, order, bin;
		final String id, source, contig;
		final Read context;
		final float headerGc;
		String fasta;
		int model=-1, captures;
		float geneGc, scavengerGc, netGc=Float.NaN, classifierGc=Float.NaN, regressionGc=Float.NaN;
		boolean measured;
	}

	private static SequenceInput sequenceInput(final String path){
		return new SequenceInput(StreamerFactory.makeStreamer(FileFormat.testInput(path, FileFormat.FASTA,
			null, true, true), null, false, -1, false, true, 1));
	}
	private static void header(final Input input, final String expected){
		check(Arrays.equals(input.reader.nextLine(), expected.getBytes(StandardCharsets.UTF_8)), "Unexpected native metadata header");
	}
	private static int number(final LineParser1 parser, final int field){
		parser.setBounds(field);
		return NcrnaNumbers.unsignedInt(parser.line(), parser.a(), parser.b());
	}
	private String file(final String name){return new File(root, name).getPath();}
	private static int bits(final float value){return Float.floatToRawIntBits(value);}
	private static void check(final boolean pass, final String reason){
		if(!pass){throw new IllegalArgumentException(reason);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String contexts, graphPath, original, roster, loci, ref, hbm, seeds, familyName;
	private long seed=-1, unavailable;
	private final File root;
	private NcrnaLearningPopulation population;
	private NcrnaFamily family;
	private NcrnaLearningTables tables;
	private int[] lengths, starts, stops;
	private boolean[] minus;
	private Sample[] samples;
	private int sourceGenomes;
	/** Optional strict parity gate; default measurement still reports differences. */
	private boolean requireExact;
	/** One comparison per emitted vector, independently counted for every feature. */
	private long featureRows;
	private final long[] scorerExact=new long[10], classifierExact=new long[10], regressionExact=new long[10];
	private final long[] binCounts=new long[BINS];
	private static final String[] FEATURE_NAMES={"identity", "tip_magnitude", "tip_center", "tip_neighbor",
		"is_stop", "coverage_left", "coverage_center", "coverage_right", "length_ratio", "contig_gc"};
	private static final int[] LENGTH_LIMITS={1000, 10000, 100000, 1000000, 16777216};
	private static final String[] LENGTH_LABELS={"1-1000", "1001-10000", "10001-100000", "100001-1000000", "1000001-16777216", ">16777216"};
	private static final int BINS=60;
}
