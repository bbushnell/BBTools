package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import idaligner.QuantumAligner;
import map.ObjectIntMap;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.Outputs;
import prok.NcrnaTextIO.SequenceInput;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;

/** Streams pooled positional matrices for all radii after one model choice/locus.
 * Adapted from A_SampleStreamerMT and NcrnaLearningVectors. Workers own indexes,
 * random generators and feature scratch. Numbered whole batches cross output
 * queues; empty batches retain their IDs. Sampling keys use original locus index,
 * radius and declared draw seed, so scheduling never changes a locus's trials.
 * Monitoring resamples the SAME eligible loci and is not independent validation.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaPositionalVectors{

	/** Writes four fresh matrices per radius and a complete original-locus ledger.
	 * Tables must have been built from this same original population and FASTA.
	 * All inputs are borrowed read-only; every output remains under directory.
	 */
	public static NcrnaPositionalVectors write(final NcrnaLearningPopulation population, final String fasta,
			final NcrnaFamily family, final NcrnaPositionalTables tables, final NcrnaBoundaryJitter[] jitter,
			final int samples, final long trainSeed, final long monitorSeed, final int workers, final File directory){
		return write(population, fasta, family, tables, jitter, samples, trainSeed, monitorSeed, -1, workers, directory);
	}

	/** Writes plain 69-input/21-target matrices, plus row-aligned identity sidecars.
	 * Explicitly excluded complete loci form the evaluation set and never teach or
	 * subtract from tables. Training and same-locus monitoring retain their LOO view.
	 * Exact assembly overlap fails before any output opens; taxonomy is optional.
	 */
	public static NcrnaPositionalVectors writeInterchange(final NcrnaLearningPopulation population, final String fasta,
			final NcrnaFamily family, final NcrnaPositionalTables tables, final NcrnaBoundaryJitter[] jitter,
			final int samples, final long trainSeed, final long monitorSeed, final long evaluationSeed,
			final int workers, final File directory){
		return writeInterchange(population, fasta, family, tables, jitter, samples, trainSeed, monitorSeed,
			evaluationSeed, workers, directory, null);
	}

	/** Optional frozen outer trials replace only evaluation draws; training and monitor draws are unchanged. */
	public static NcrnaPositionalVectors writeInterchange(final NcrnaLearningPopulation population, final String fasta,
			final NcrnaFamily family, final NcrnaPositionalTables tables, final NcrnaBoundaryJitter[] jitter,
			final int samples, final long trainSeed, final long monitorSeed, final long evaluationSeed,
			final int workers, final File directory, final NcrnaFixedEndpointTrials fixed){
		if(tables==null || !tables.triple() || tables.radii.length!=1 || tables.radii[0]!=10
				|| evaluationSeed<0 || evaluationSeed==trainSeed || evaluationSeed==monitorSeed){
			throw new IllegalArgumentException("Interchange requires three K9 tables, radius10 and three distinct nonnegative draw seeds");
		}
		if(fixed!=null && (!fixed.belongsTo(population, family) || tables.fullErrors==null)){
			throw new IllegalArgumentException("Fixed centered trials require this exact population and full-error weighted tables");
		}
		return write(population, fasta, family, tables, jitter, samples, trainSeed, monitorSeed, evaluationSeed, workers, directory, fixed);
	}

	/** Full-error product refit has only training and same-locus monitoring rows.
	 * It keeps the interchange schema and LOO features without inventing an outer set.
	 */
	public static NcrnaPositionalVectors writeRefit(final NcrnaLearningPopulation population, final String fasta,
			final NcrnaFamily family, final NcrnaPositionalTables tables, final int samples,
			final long trainSeed, final long monitorSeed, final int workers, final File directory){
		if(tables==null || !tables.triple() || tables.fullErrors==null || tables.radii.length!=1 || tables.radii[0]!=10){
			throw new IllegalArgumentException("All-data positional refit requires the full-error K9 three-table radius-ten schema");
		}
		return write(population, fasta, family, tables, null, samples, trainSeed, monitorSeed, -1, workers, directory, null, true);
	}

	private static NcrnaPositionalVectors write(final NcrnaLearningPopulation population, final String fasta,
			final NcrnaFamily family, final NcrnaPositionalTables tables, final NcrnaBoundaryJitter[] jitter,
			final int samples, final long trainSeed, final long monitorSeed, final long evaluationSeed,
			final int workers, final File directory){
		return write(population, fasta, family, tables, jitter, samples, trainSeed, monitorSeed, evaluationSeed, workers, directory, null);
	}

	private static NcrnaPositionalVectors write(final NcrnaLearningPopulation population, final String fasta,
			final NcrnaFamily family, final NcrnaPositionalTables tables, final NcrnaBoundaryJitter[] jitter,
			final int samples, final long trainSeed, final long monitorSeed, final long evaluationSeed,
			final int workers, final File directory, final NcrnaFixedEndpointTrials fixed){
		return write(population, fasta, family, tables, jitter, samples, trainSeed, monitorSeed, evaluationSeed,
			workers, directory, fixed, false);
	}

	private static NcrnaPositionalVectors write(final NcrnaLearningPopulation population, final String fasta,
			final NcrnaFamily family, final NcrnaPositionalTables tables, final NcrnaBoundaryJitter[] jitter,
			final int samples, final long trainSeed, final long monitorSeed, final long evaluationSeed,
			final int workers, final File directory, final NcrnaFixedEndpointTrials fixed, final boolean refit){
		if(population==null || family==null || family.library==null || family.library.length==0
				|| family.modelNames==null || family.modelNames.length!=family.library.length || tables==null
				|| !tables.belongsTo(population) || tables.statusByGene.length!=population.ids.length || tables.records!=population.ids.length
				|| population.split().selected!=0 || samples<1 || workers<1 || trainSeed<0 || monitorSeed<0
				|| trainSeed==monitorSeed || directory==null || !directory.isDirectory()){
			throw new IllegalArgumentException("Positional matrices require the complete original axis, finalized tables, explicit samples/seeds, workers and an existing output directory");
		}
		final NcrnaPositionalVectors result=new NcrnaPositionalVectors(population, family, tables, jitter,
			samples, trainSeed, monitorSeed, evaluationSeed, directory, fixed, refit);
		result.process(fasta, workers);
		result.statistics(new File(directory, "vectors.stats.tsv").getPath());
		return result;
	}

	private NcrnaPositionalVectors(final NcrnaLearningPopulation population_, final NcrnaFamily family_,
			final NcrnaPositionalTables tables_, final NcrnaBoundaryJitter[] jitter_, final int samples_,
			final long trainSeed, final long monitorSeed, final long evaluationSeed, final File directory,
			final NcrnaFixedEndpointTrials fixed_, final boolean refit_){
		population=population_;
		refit=refit_;
		fixed=fixed_;
		family=family_;
		tables=tables_;
		fullErrors=tables.fullErrors;
		plateaus=tables.triple() ? new float[family.library.length] : null;
		if(plateaus!=null){
			if(family.models==null || family.models.length!=family.library.length){
				throw new IllegalArgumentException("Three-probe HBM features require the selected library's paired models");
			}
			for(int model=0; model<plateaus.length; model++){
				if(family.models[model]==null || !Arrays.equals(family.models[model].original, family.library[model])){
					throw new IllegalArgumentException("Three-probe consensus and HBM sequences disagree at model "+model);
				}
				plateaus[model]=TrnaBoundaryFeatures.modelMaxCoverage(family.models[model]);
				if(plateaus[model]<=0){throw new IllegalArgumentException("Three-probe HBM features require positive model coverage");}
			}
		}
		final NcrnaFamilyRegistry.Entry entry=NcrnaFamilyRegistry.get(family.name);
		final String member=entry.group==null ? entry.id : entry.member;
		if(!member.equals(tables.member)){throw new IllegalArgumentException("Positional matrices must use the registered member's own tables and population");}
		samples=samples_;
		interchange=refit || evaluationSeed>=0;
		if(fullErrors!=null && !interchange){throw new IllegalArgumentException("Full-error vectors require the explicit 69-input/22-output interchange mode");}
		seeds=interchange && !refit ? new long[]{trainSeed, monitorSeed, evaluationSeed} : new long[]{trainSeed, monitorSeed};
		if(interchange){validateAssemblies();}
		if(fullErrors==null){
			if(jitter_==null || jitter_.length!=tables.radii.length){throw new IllegalArgumentException("Each positional radius requires its own observed conditioned histograms");}
			jitter=jitter_.clone();
			for(int r=0; r<jitter.length; r++){
				if(jitter[r]==null || jitter[r].radius!=tables.radii[r] || !tables.member.equals(jitter[r].member)
						|| !jitter[r].available(0) || !jitter[r].available(1)){
					throw new IllegalArgumentException("Every member/radius/end needs observed in-window jitter; empty windows cannot produce invented training rows");
				}
			}
		}else{jitter=null;}
		final int count=Math.multiplyExact(tables.radii.length, 2*seeds.length);
		paths=new String[count];
		sidecarPaths=interchange ? new String[count] : null;
		rows=new long[count];
		attempts=new long[count];
		missingWord=new long[count];
		clippedWord=new long[count];
		invalidSpan=new long[count];
		crossedSpan=new long[count];
		outsideContext=new long[count];
		belowMinimumSpan=new long[count];
		zeroAttempts=new long[count];
		zeroWritten=new long[count];
		outOfRangeAttempts=new long[count];
		outOfRangeWritten=new long[count];
		writers=new ByteStreamWriter[(interchange ? 2*count : count)+1+(fixed==null ? 0 : 2)];
		fixedAuditPaths=fixed==null ? null : new String[]{new File(directory, "r10.evaluation.5prime.centered.tsv").getPath(),
			new File(directory, "r10.evaluation.3prime.centered.tsv").getPath()};
		seen=new boolean[population.ids.length];
		modelNames=new String[family.modelNames.length];
		for(int model=0; model<modelNames.length; model++){modelNames[model]=NcrnaSupport.recordName(family.modelNames[model]);}
		for(int r=0; r<tables.radii.length; r++){
			for(int set=0; set<seeds.length; set++){
				for(int end=0; end<2; end++){
					final int file=matrixIndex(r, set, end);
					final String name="r"+tables.radii[r]+'.'+SETS[set]+'.'+ENDS[end];
					paths[file]=new File(directory, name+(interchange ? ".tsv" : ".tsv.gz")).getPath();
					if(interchange){sidecarPaths[file]=new File(directory, name+".rows.tsv").getPath();}
				}
			}
		}
		lociPath=new File(directory, "vectors.loci.tsv").getPath();
	}

	/** Checks identities once on retained metadata, never by a per-vector lookup. */
	private void validateAssemblies(){
		final ObjectIntMap<String> membership=new ObjectIntMap<String>(String.class);
		int evaluation=0;
		for(int gene=0; gene<population.ids.length; gene++){
			if(refit && population.excluded[gene]){
				throw new IllegalArgumentException("All-data refit cannot withhold original locus "+population.ids[gene]);
			}
			if(population.partial[gene]){continue;}
			final String assembly=population.assemblyName(gene);
			if(assembly.equals(".") || assembly.equals("NA")){
				throw new IllegalArgumentException("Interchange split requires an exact assembly identity: "+population.ids[gene]);
			}
			final int mask=population.excluded[gene] ? 2 : 1;
			final int previous=membership.get(assembly);
			if(previous>=0 && previous!=mask){
				throw new IllegalArgumentException("Training and evaluation share assembly "+assembly);
			}
			membership.put(assembly, mask);
			if(mask==2){evaluation++;}
		}
		if(!refit && evaluation==0){throw new IllegalArgumentException("Interchange requires complete explicitly excluded evaluation loci");}
	}

	/** Dense batch IDs bound queued memory and keep row order independent of worker completion. */
	private void process(final String fasta, final int workers){
		final FileFormat format=FileFormat.testInput(fasta, FileFormat.FASTA, null, true, true);
		if(!format.fasta()){throw new IllegalArgumentException("Positional contexts must be FASTA");}
		final boolean trim=Shared.TRIM_READ_DESCRIPTION, upper=Read.TO_UPPER_CASE;
		final boolean bgzf=ReadWrite.USE_BGZF, force=ReadWrite.FORCE_BGZIP, zipExplicit=ReadWrite.SET_ZIP_THREADS;
		final int buffer=Shared.bufferLen(), zipThreads=ReadWrite.MAX_ZIP_THREADS();
		Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=true;
		Shared.setBufferLen(16);
		ReadWrite.USE_BGZF=ReadWrite.FORCE_BGZIP=true;
		// Twelve independent matrices must not each start a full CPU-budget compressor.
		ReadWrite.setZipThreads(1);
		try(NcrnaWorkers cohort=new NcrnaWorkers(); Outputs outputs=new Outputs(writers.length)){
			for(int r=0; r<tables.radii.length; r++){
				for(int set=0; set<seeds.length; set++){
					for(int end=0; end<2; end++){
						final int file=matrixIndex(r, set, end);
						writers[file]=outputs.open(FileFormat.testOutput(paths[file], FileFormat.TEXT, null, true, false, false, true));
						if(interchange){
							writers[file].add(new ByteBuilder().append("#dims\t69\t").append(classes(r)).append("\t0\n"), 0);
							writers[paths.length+file]=outputs.open(FileFormat.testOutput(sidecarPaths[file], FileFormat.TEXT, null, true, false, false, true));
							writers[paths.length+file].add(new ByteBuilder("locus_id\tsource\tassembly\tsplit\tmodel_index\tmodel\tsample\traw_start0\traw_stop0\ttrue_start0\ttrue_stop0\n"), 0);
							continue;
						}
						final ByteBuilder header=new ByteBuilder().append("#dims\t").append(tables.inputs(r)).tab().append(tables.five[r].sites)
							.append("\n#member\t").append(tables.member).append("\tradius\t").append(tables.radii[r]).append("\tend\t").append(ENDS[end])
							.append("\tdraw_seed\t").append(seeds[set]).append("\tsamples\t").append(samples)
							.append("\n#monitoring\tSAME_ELIGIBLE_LOCI_NOT_INDEPENDENT\n");
						if(tables.triple()){
							header.append("#features\t").append(PackagedPositionalBoundaryConfig.TRIPLE_FEATURE_ORDER)
								.append("\n#training_table_view\tLEAVE_ONE_LOCUS_OUT\n");
						}
						writers[file].add(header, 0);
					}
				}
			}
			if(fixed!=null){
				for(int end=0; end<2; end++){
					writers[2*paths.length+end]=outputs.open(FileFormat.testOutput(fixedAuditPaths[end], FileFormat.TEXT, null, true, false, false, true));
					writers[2*paths.length+end].add(new ByteBuilder(NcrnaFixedEndpointTrials.AUDIT_HEADER).nl(), 0);
				}
			}
			writers[writers.length-1]=outputs.open(FileFormat.testOutput(lociPath, FileFormat.TEXT, null, true, false, false, true));
			writers[writers.length-1].add(new ByteBuilder("locus_id\tsource\tmodel_index\tmodel\tcontext_status\tdisposition\n"), 0);
			try(SequenceInput owner=new SequenceInput(StreamerFactory.makeStreamer(format, null, true, -1, false, true, 1))){
				final ArrayList<ProcessThread> threads=new ArrayList<ProcessThread>(workers);
				for(int i=0; i<workers; i++){threads.add(new ProcessThread(owner.stream, this));}
				owner.stream.start();
				cohort.startAndJoin(threads, owner.stream);
				for(final ProcessThread thread:threads){accumulate(thread);}
			}
		}finally{
			Shared.TRIM_READ_DESCRIPTION=trim;
			Read.TO_UPPER_CASE=upper;
			Shared.setBufferLen(buffer);
			ReadWrite.USE_BGZF=bgzf;
			ReadWrite.FORCE_BGZIP=force;
			ReadWrite.setZipThreads(zipThreads);
			ReadWrite.SET_ZIP_THREADS=zipExplicit;
		}
		if(records!=population.ids.length){throw new IllegalArgumentException("Every original locus needs one vector-pass disposition");}
		assert(records==filtered+unavailable+noModel+processed) : "Original contexts must partition without dropping or duplicating loci";
		for(int file=0; file<paths.length; file++){
			final boolean evaluation=interchange && (file/2)%seeds.length==2;
			if(rows[file]==0 && !(evaluation && fixed!=null)){throw new IllegalArgumentException("No usable positional rows in "+paths[file]);}
			final long expected=evaluation && fixed!=null ? fixed.ends[file%2].size() : Math.multiplyExact(evaluation ? evaluationProcessed : processed, samples);
			assert(attempts[file]==expected
					&& attempts[file]==rows[file]+missingWord[file]+invalidSpan[file])
				: "Each requested member/end trial is written or has exactly one explicit limitation: "+paths[file];
		}
		if(fixed!=null){validateFixedOrder();}
	}

	/** Ordered batch emission must preserve the original frozen sidecar sequence, including holds. */
	private void validateFixedOrder(){
		final parse.LineParser1 p=new parse.LineParser1('\t');
		for(int end=0; end<2; end++){
			int count=0;
			try(NcrnaTextIO.Input input=new NcrnaTextIO.Input(fixedAuditPaths[end])){
				RrnaResourceIO.header(input.reader, NcrnaFixedEndpointTrials.AUDIT_HEADER, fixedAuditPaths[end]);
				for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
					p.set(line);
					if(p.terms()!=21 || p.parseInt(11)!=++count){throw new IllegalArgumentException("Fixed outer trial order changed during context rematerialization");}
				}
			}
			if(count!=fixed.ends[end].size()){throw new IllegalArgumentException("Every original trial, including a held span, must remain in the centered ledger");}
		}
	}

	private void accumulate(final ProcessThread worker){
		if(!worker.success){throw new IllegalStateException("Positional vector worker did not finish");}
		for(int i=0; i<worker.genes.size; i++){
			final int gene=worker.genes.get(i);
			if(seen[gene]){throw new IllegalArgumentException("Duplicate original vector context: "+population.ids[gene]);}
			seen[gene]=true;
		}
		records+=worker.genes.size;
		filtered+=worker.filtered;
		unavailable+=worker.unavailable;
		noModel+=worker.noModel;
		processed+=worker.processed;
		evaluationLoci+=worker.evaluationLoci;
		evaluationProcessed+=worker.evaluationProcessed;
		evaluationUnavailable+=worker.evaluationUnavailable;
		evaluationNoModel+=worker.evaluationNoModel;
		evaluationNoTrials+=worker.evaluationNoTrials;
		assert(evaluationLoci==evaluationProcessed+evaluationUnavailable+evaluationNoModel+evaluationNoTrials)
			: "Complete excluded loci must have one explicit evaluation disposition";
		for(int file=0; file<paths.length; file++){
			rows[file]+=worker.rows[file];
			attempts[file]+=worker.attempts[file];
			missingWord[file]+=worker.missing[file];
			clippedWord[file]+=worker.clipped[file];
			invalidSpan[file]+=worker.invalid[file];
			crossedSpan[file]+=worker.crossed[file];
			outsideContext[file]+=worker.outside[file];
			belowMinimumSpan[file]+=worker.shortSpans[file];
			zeroAttempts[file]+=worker.zeroAttempts[file];
			zeroWritten[file]+=worker.zeroWritten[file];
			outOfRangeAttempts[file]+=worker.outOfRangeAttempts[file];
			outOfRangeWritten[file]+=worker.outOfRangeWritten[file];
			assert(fullErrors==null || invalidSpan[file]==crossedSpan[file]+outsideContext[file]+belowMinimumSpan[file])
				: "Full-error invalid trials partition into crossed, out-of-context and optional centered below-minimum spans";
			assert(fullErrors==null || (zeroAttempts[file]>=0 && zeroWritten[file]>=0
				&& outOfRangeAttempts[file]>=0 && outOfRangeWritten[file]>=0
				&& zeroWritten[file]<=zeroAttempts[file]
				&& outOfRangeWritten[file]<=outOfRangeAttempts[file]
				&& zeroAttempts[file]+outOfRangeAttempts[file]<=attempts[file]
				&& zeroWritten[file]+outOfRangeWritten[file]<=rows[file]
				&& rows[file]-zeroWritten[file]-outOfRangeWritten[file]
					<=attempts[file]-zeroAttempts[file]-outOfRangeAttempts[file]))
				: "Correct and out-of-range trial counts are disjoint subsets; unavailable geometry cannot create replacement rows";
		}
	}

	/** One row per real matrix keeps denominators separate from original-locus counts. */
	private void statistics(final String path){
		try(Output output=new Output(path)){
			final ByteBuilder buffer=new ByteBuilder();
			buffer.append("member\tradius_bp\tset\tend\toriginals\tfiltered\tunavailable\tno_model\tprocessed\tsamples\tattempts\twritten\tmissing_word\tinvalid_span\tdraw_seed\tmonitoring");
			if(interchange){buffer.append("\tclipped_word\tambiguous_word\tmatrix\tsidecar\tinputs\toutputs\tweights\tscaling\tfeature_order\tlabel_rule");}
			if(fullErrors!=null){buffer.append("\tcrossed_span\toutside_context\tdraw_distribution\tbackground\tzero_attempts\tzero_written\tout_of_range_attempts\tout_of_range_written");}
			if(fixed!=null){buffer.append("\tfixed_outer_trials\touter_loci_without_frozen_trials\tcentered_below_minimum\tall_trial_ledger");}
			output.writer.print(buffer.nl());
			for(int r=0; r<tables.radii.length; r++){
				for(int set=0; set<seeds.length; set++){
					for(int end=0; end<2; end++){
						final int file=matrixIndex(r, set, end);
						final boolean evaluation=set==2;
						buffer.clear().append(tables.member).tab().append(tables.radii[r]).tab().append(SETS[set]).tab().append(ENDS[end])
							.tab().append(evaluation ? evaluationLoci : records).tab().append(evaluation ? 0 : filtered)
							.tab().append(evaluation ? evaluationUnavailable : unavailable).tab().append(evaluation ? evaluationNoModel : noModel)
							.tab().append(evaluation ? evaluationProcessed : processed)
							.tab().append(samples).tab().append(attempts[file]).tab().append(rows[file]).tab().append(missingWord[file]).tab().append(invalidSpan[file])
							.tab().append(seeds[set]).tab().append(evaluation ? "PANEL_EXCLUDED_NO_SUBTRACTION"
								: refit ? "ALL_DATA_REFIT_SAME_LOCUS_MONITOR_NOT_INDEPENDENT" : "SAME_ELIGIBLE_LOCI_NOT_INDEPENDENT");
						if(interchange){
							buffer.tab().append(clippedWord[file]).tab().append(missingWord[file]-clippedWord[file])
								.tab().append(paths[file]).tab().append(sidecarPaths[file]).append("\t69\t").append(classes(r)).append("\t0\tRAW_NO_TRANSFORM\t")
								.append(PackagedPositionalBoundaryConfig.TRIPLE_FEATURE_ORDER).tab()
								.append(fullErrors==null ? "TRUE_END_MINUS_RAW_END_PLUS_RADIUS" : "TRUE_END_MINUS_RAW_END_PLUS_RADIUS_OR_OUT_OF_RANGE");
						}
						if(fullErrors!=null){
							buffer.tab().append(crossedSpan[file]).tab().append(outsideContext[file])
								.tab().append(evaluation && fixed!=null ? "FIXED_ORIGINAL_OUTER_TRIALS_CENTERED_NO_REDRAW" : "FULL_SIGNED_INCLUDING_ZERO")
								.append("\tNONZERO_SIGNED_WEIGHTED_SEPARATE_DENOMINATORS")
								.tab().append(zeroAttempts[file]).tab().append(zeroWritten[file])
								.tab().append(outOfRangeAttempts[file]).tab().append(outOfRangeWritten[file]);
						}
						if(fixed!=null){buffer.tab().append(evaluation ? fixed.ends[end].size() : 0)
							.tab().append(evaluation ? evaluationNoTrials : 0).tab().append(belowMinimumSpan[file])
							.tab().append(evaluation ? fixedAuditPaths[end] : "NA");}
						output.writer.print(buffer.nl());
					}
				}
			}
		}
	}

	private static final class ProcessThread extends Thread{
		ProcessThread(final Streamer input_, final NcrnaPositionalVectors owner_){
			input=input_;
			owner=owner_;
			final NcrnaFamily family=owner.family;
			index=new TrnaKmerIndex(family.library, family.indexK, family.adaptive,
				family.adaptFloor, family.adaptTopFrac, family.adaptQFrac, family.fixedMinHits);
			features=new NcrnaPositionalBoundaryFeatures[owner.tables.radii.length];
			evaluationFeatures=owner.interchange ? new NcrnaPositionalBoundaryFeatures[features.length] : null;
			for(int r=0; r<features.length; r++){
				features[r]=owner.fullErrors!=null
					? new NcrnaPositionalBoundaryFeatures(owner.tables.fiveProbes[r], owner.tables.threeProbes[r],
						owner.tables.fiveWeighted[r], owner.tables.threeWeighted[r], true)
					: owner.tables.triple()
					? new NcrnaPositionalBoundaryFeatures(owner.tables.fiveProbes[r], owner.tables.threeProbes[r], true)
					: new NcrnaPositionalBoundaryFeatures(owner.tables.five[r], owner.tables.three[r]);
				if(owner.interchange){
					evaluationFeatures[r]=owner.fullErrors==null
						? new NcrnaPositionalBoundaryFeatures(owner.tables.fiveProbes[r], owner.tables.threeProbes[r], false)
						: new NcrnaPositionalBoundaryFeatures(owner.tables.fiveProbes[r], owner.tables.threeProbes[r],
							owner.tables.fiveWeighted[r], owner.tables.threeWeighted[r], false);
				}
			}
			rows=new long[owner.paths.length];
			attempts=new long[rows.length];
			missing=new long[rows.length];
			clipped=new long[rows.length];
			invalid=new long[rows.length];
			crossed=new long[rows.length];
			outside=new long[rows.length];
			shortSpans=new long[rows.length];
			zeroAttempts=new long[rows.length];
			zeroWritten=new long[rows.length];
			outOfRangeAttempts=new long[rows.length];
			outOfRangeWritten=new long[rows.length];
		}
		@Override
		public void run(){
			try{
				for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
					buffers=new ByteBuilder[owner.writers.length];
					for(int file=0; file<buffers.length; file++){buffers[file]=new ByteBuilder(1024);}
					for(final Read read:batch){processRead(read);}
					for(int file=0; file<buffers.length; file++){owner.writers[file].add(buffers[file], batch.id+1);}
					buffers=null;
				}
				success=true;
			}catch(Throwable failure){KillSwitch.exceptionKill(failure);}
		}

		private void processRead(final Read read){
			context.set(read, owner.population);
			final int gene=context.gene;
			genes.add(gene);
			if(context.status!=owner.tables.statusByGene[gene]){throw new IllegalArgumentException("Context availability changed between table/vector passes: "+read.id);}
			final int model;
			final String disposition;
			final boolean evaluation=owner.interchange && owner.population.excluded[gene] && !owner.population.partial[gene];
			if(evaluation){
				// Keep teaching exclusions in the old axis; evaluation counts are separate.
				filtered++;
				evaluationLoci++;
				model=owner.fixed!=null ? owner.fixed.modelByGene[gene] : context.hasFeatures() ? selectModel() : -1;
				if(owner.fixed!=null){
					if(model<0){
						evaluationNoTrials++;
						disposition="EVALUATION_NO_FROZEN_TRIALS";
					}else{
						if(!context.hasFeatures()){throw new IllegalArgumentException("A frozen original evaluation context lost its required features");}
						evaluationProcessed++;
						disposition="EVALUATION_FIXED_TRIALS_REMATERIALIZED";
						emitFixed();
					}
				}else if(!context.hasFeatures()){
					evaluationUnavailable++;
					disposition="EVALUATION_UNAVAILABLE_CONTEXT";
				}else if(model<0){
					evaluationNoModel++;
					disposition="EVALUATION_NO_SHORTLISTED_MODEL";
				}else{
					evaluationProcessed++;
					disposition="EVALUATION_VECTORIZED_ORIGINAL";
					for(int r=0; r<features.length; r++){emit(r, model, true);}
				}
			}else if(!owner.population.eligible(gene)){
				filtered++;
				model=-1;
				disposition=owner.interchange && owner.population.partial[gene] ? "PARTIAL_ORIGINAL" : "FILTERED_ORIGINAL";
			}else if(!context.hasFeatures()){
				unavailable++;
				model=-1;
				disposition="UNAVAILABLE_CONTEXT";
			}else{
				model=selectModel();
				if(model<0){
					noModel++;
					disposition="NO_SHORTLISTED_MODEL";
				}else{
					processed++;
					disposition="VECTORIZED_ORIGINAL";
					for(int r=0; r<features.length; r++){emit(r, model, false);}
				}
			}
			buffers[owner.writers.length-1].append(owner.population.ids[gene]).tab().append(owner.population.sourceName(gene)).tab().append(model)
				.tab().append(model<0 ? "." : owner.modelNames[model]).tab().append(context.status).tab().append(disposition).nl();
		}

		/** Same shortlist/highest-Quantum-identity/first-tie rule as NcrnaBoundaryExamples.prepare. */
		private int selectModel(){
			final byte[] gene=Arrays.copyOfRange(context.read.bases, context.start, context.stop+1);
			final int[] shortlist=index.shortlist(gene, owner.family.indexTopN);
			int selected=-1;
			float best=-1;
			for(final int model:shortlist){
				final float identity=QuantumAligner.alignStatic(gene, owner.family.library[model], null);
				if(!(identity>=0 && identity<=1)){throw new IllegalStateException("Nonfinite model-selection identity cannot define positional examples");}
				if(identity>best){
					best=identity;
					selected=model;
				}
			}
			return selected;
		}

		private void emit(final int r, final int model, final boolean evaluation){
			final int radius=owner.tables.radii[r];
			final NcrnaPositionalBoundaryFeatures selected=evaluation ? evaluationFeatures[r] : features[r];
			if(owner.tables.triple() && !evaluation){selected.setLeaveOneOut(context.read.bases, context.start, context.stop);}
			for(int set=evaluation ? 2 : 0; set<(evaluation ? 3 : 2); set++){
				final long seed=Tools.hash64shift(owner.seeds[set]^((long)context.gene<<32)^radius);
				random[0].setSeed(seed);
				random[1].setSeed(seed^0x5DEECE66DL);
				for(int sample=0; sample<owner.samples; sample++){
					final long start=(long)context.start+draw(r, 0);
					final long stop=(long)context.stop+draw(r, 1);
					final boolean valid=start>=0 && stop>=start && stop<context.read.length();
					final int usable=valid ? (owner.tables.triple()
						? selected.fill(context.read.bases, (int)start, (int)stop, owner.family.library[model],
							owner.family.models[model], owner.plateaus[model], context.gc, 0)
						: selected.fill(context.read.bases, (int)start, (int)stop, owner.family.library[model], context.gc, 0)) : 0;
					for(int end=0; end<2; end++){
						final int file=owner.matrixIndex(r, set, end);
						final long correction=end==0 ? context.start-start : context.stop-stop;
						final boolean outsideCorrection=correction<-radius || correction>radius;
						attempts[file]++;
						if(owner.fullErrors!=null){
							if(correction==0){zeroAttempts[file]++;}
							if(outsideCorrection){outOfRangeAttempts[file]++;}
						}
						if(!valid){
							invalid[file]++;
							if(owner.fullErrors!=null){
								if(stop<start){crossed[file]++;}
								else{outside[file]++;}
							}
						}else if((usable&(1<<end))==0){
							missing[file]++;
							if(owner.interchange && clippedWord(r, end, (int)(end==0 ? start : stop))){clipped[file]++;}
						}else{
							final int truth=end==0 ? context.start : context.stop;
							final int raw=(int)(end==0 ? start : stop);
							final int label=owner.fullErrors!=null && outsideCorrection
								? 2*radius+1 : RrnaEndpointVector.labelIndex(truth, raw, radius);
							final ByteBuilder buffer=buffers[file];
							for(final float value:selected.values[end]){buffer.append(value, 18).tab();}
							for(int site=0; site<owner.classes(r); site++){
								if(site>0){buffer.tab();}
								buffer.append(site==label ? 1 : 0);
							}
							buffer.nl();
							if(owner.interchange){sidecar(file, set, model, sample, (int)start, (int)stop);}
							rows[file]++;
							if(owner.fullErrors!=null){
								if(correction==0){zeroWritten[file]++;}
								if(outsideCorrection){outOfRangeWritten[file]++;}
							}
						}
					}
				}
			}
		}

		private int draw(final int r, final int end){
			return owner.fullErrors==null ? owner.jitter[r].draw(end, random[end])
				: owner.fullErrors.distribution(end).draw(random[end]);
		}

		/** Rebuild features at shifted bounds without redrawing or replacing any outer trial. */
		private void emitFixed(){
			assert(owner.fixed!=null && owner.fullErrors!=null && evaluationFeatures.length==1)
				: "The fixed-trial API is limited to the qualified radius10 full-error feature layout";
			final NcrnaPositionalBoundaryFeatures selected=evaluationFeatures[0];
			for(int end=0; end<2; end++){
				final NcrnaFixedEndpointTrials.End trials=owner.fixed.ends[end];
				final int file=owner.matrixIndex(0, 2, end);
				for(int row=trials.first[context.gene]; row>=0; row=trials.next[row]){
					owner.fixed.center(end, row, context, centered);
					final int model=trials.model[row];
					final long correction=(end==0 ? context.start-centered.start : context.stop-centered.stop);
					attempts[file]++;
					if(correction==0){zeroAttempts[file]++;}
					if(correction<-10 || correction>10){outOfRangeAttempts[file]++;}
					final String disposition;
					if(!centered.valid){
						invalid[file]++;
						if(centered.crossed){crossed[file]++;
						}else if(centered.outside){outside[file]++;
						}else{shortSpans[file]++;}
						disposition=centered.crossed ? "CENTER_CROSSED_HOLD" : centered.outside
							? "CENTER_OUTSIDE_CONTEXT_HOLD" : "CENTER_BELOW_MINIMUM_HOLD";
					}else{
						final int usable=selected.fill(context.read.bases, (int)centered.start, (int)centered.stop,
							owner.family.library[model], owner.family.models[model], owner.plateaus[model], context.gc, 0);
						if((usable&(1<<end))==0){
							missing[file]++;
							if(clippedWord(0, end, (int)(end==0 ? centered.start : centered.stop))){clipped[file]++;}
							disposition="CENTER_MISSING_FEATURE_HOLD";
						}else{
							final int label=correction<-10 || correction>10 ? 21 : (int)correction+10;
							for(final float value:selected.values[end]){buffers[file].append(value, 18).tab();}
							for(int site=0; site<22; site++){
								if(site>0){buffers[file].tab();}
								buffers[file].append(site==label ? 1 : 0);
							}
							buffers[file].nl();
							sidecar(file, 2, model, trials.sample[row], (int)centered.start, (int)centered.stop);
							rows[file]++;
							if(correction==0){zeroWritten[file]++;}
							if(label==21){outOfRangeWritten[file]++;}
							disposition="VECTORIZED";
						}
					}
					buffers[2*owner.paths.length+end].append(trials.original[row]).tab().append(row+1)
						.tab().append(centered.start).tab().append(centered.stop)
						.tab().append(owner.fixed.median5[model]).tab().append(owner.fixed.median3[model]).tab().append(disposition)
						.tab().append(owner.fixed.trainingSupport[model]).tab().append(owner.fixed.supportMinimum)
						.tab().append(owner.fixed.supported[model] ? 1 : 0).tab().append(owner.fixed.minimumLength).nl();
				}
			}
		}

		/** A sidecar row is appended only with its numeric row, in the same batch. */
		private void sidecar(final int file, final int set, final int model, final int sample, final int start, final int stop){
			assert(start>=0 && stop>=start && stop<context.read.length())
				: "Only valid emitted candidate spans receive row metadata";
			final int gene=context.gene;
			buffers[owner.paths.length+file].append(owner.population.ids[gene]).tab()
				.append(owner.population.sourceName(gene)).tab().append(owner.population.assemblyName(gene))
				.tab().append(SETS[set]).tab().append(model).tab().append(owner.modelNames[model]).tab().append(sample)
				.tab().append(start).tab().append(stop).tab().append(context.start).tab().append(context.stop).nl();
		}

		/** Bounds distinguish clipped score windows from ambiguous in-bounds words.
		 * When both occur, the trial is classified as clipped, never counted twice.
		 */
		private boolean clippedWord(final int r, final int end, final int rawEnd){
			final RrnaPositionalKmerTable.Table[] probes=end==0 ? owner.tables.fiveProbes[r] : owner.tables.threeProbes[r];
			for(final RrnaPositionalKmerTable.Table table:probes){
				final long anchor=(long)rawEnd+table.offset;
				if(anchor-table.radius<0 || anchor+table.radius+table.k>context.read.length()){return true;}
			}
			return false;
		}

		private final Streamer input;
		private final NcrnaPositionalVectors owner;
		private final TrnaKmerIndex index;
		private final NcrnaPositionalBoundaryFeatures[] features;
		private final NcrnaPositionalBoundaryFeatures[] evaluationFeatures;
		private final NcrnaLearningContext context=new NcrnaLearningContext();
		private final NcrnaFixedEndpointTrials.Span centered=new NcrnaFixedEndpointTrials.Span();
		private final Random[] random={new Random(0), new Random(0)};
		private final IntList genes=new IntList();
		private final long[] rows, attempts, missing, clipped, invalid;
		private final long[] crossed, outside, shortSpans;
		private final long[] zeroAttempts, zeroWritten, outOfRangeAttempts, outOfRangeWritten;
		private ByteBuilder[] buffers;
		long filtered, unavailable, noModel, processed;
		long evaluationLoci, evaluationUnavailable, evaluationNoModel, evaluationProcessed, evaluationNoTrials;
		boolean success;
	}

	/** Radius, resampling set and endpoint axes map to one native matrix file. */
	public static int index(final int radiusIndex, final int set, final int end){return (radiusIndex<<2)|(set<<1)|end;}
	/** Instance-specific matrix axis; the legacy static index remains two-set only. */
	public int matrixIndex(final int radiusIndex, final int set, final int end){return 2*(radiusIndex*seeds.length+set)+end;}
	private int classes(final int r){return tables.five[r].sites+(fullErrors==null ? 0 : 1);}
	public final String[] paths;
	/** One row-aligned identity file per plain matrix; null in legacy output mode. */
	public final String[] sidecarPaths;
	/** Complete original trial identities, both coordinate frames and explicit holds; optional. */
	public final String[] fixedAuditPaths;
	public final long[] rows, attempts, missingWord, invalidSpan;
	/** Subset of missingWord caused by an out-of-context positional window. */
	public final long[] clippedWord;
	/** Optional full-error invalid-span partition, never silently redrawn. */
	public final long[] crossedSpan, outsideContext;
	/** Opt-in centered trial spans that remain in-context but violate the loaded family's minimum length. */
	public final long[] belowMinimumSpan;
	/** Optional class-mixture accounting before and after explicit geometry losses. */
	public final long[] zeroAttempts, zeroWritten, outOfRangeAttempts, outOfRangeWritten;
	public final String lociPath;
	public long records, filtered, unavailable, noModel, processed;
	/** Complete excluded population and its disjoint evaluation dispositions. */
	public long evaluationLoci, evaluationUnavailable, evaluationNoModel, evaluationProcessed, evaluationNoTrials;
	private final NcrnaFixedEndpointTrials fixed;
	private final NcrnaLearningPopulation population;
	private final NcrnaFamily family;
	private final NcrnaPositionalTables tables;
	private final NcrnaSignedEndpointErrors fullErrors;
	/** Immutable selected-model maximum HBM coverage, optional layout only. */
	private final float[] plateaus;
	private final NcrnaBoundaryJitter[] jitter;
	private final int samples;
	private final long[] seeds;
	private final String[] modelNames;
	private final boolean[] seen;
	private final ByteStreamWriter[] writers;
	private final boolean interchange, refit;
	private static final String[] SETS={"train", "monitor", "evaluation"}, ENDS={"5prime", "3prime"};
}
