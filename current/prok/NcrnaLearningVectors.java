package prok;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import prok.NcrnaTextIO.Outputs;
import prok.NcrnaTextIO.SequenceInput;
import shared.KillSwitch;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;

/** Streams native ten-input/one-label boundary matrices into length-balanced files.
 * Adapted from A_SampleStreamerMT. Workers own feature/index/LOO scratch and
 * transfer whole byte blocks to shared writers, never mutate their print buffers.
 * Whole loci stay in one file. Validation is separate during tuning; refit has
 * only training files and uses leave-one-out for every eligible locus.
 * Header/case/compression settings are process-wide and restored before return;
 * do not run this simultaneously with another pipeline changing these settings.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaLearningVectors{

	/** Writes fresh files at prefix.train.N.tsv.gz and, for tuning, prefix.validation.tsv.gz. */
	public static NcrnaLearningVectors write(final NcrnaLearningPopulation population, final String fasta,
			final NcrnaFamily family, final NcrnaLearningTables tables, final int startRadius,
			final int stopRadius, final int shards, final int workers, final String prefix){
		if(population==null || tables==null || family==null || family.library==null
				|| family.library.length==0 || shards<1 || workers<1 || prefix==null){
			throw new IllegalArgumentException("Vector writing needs population, tables, output prefix and positive shard/worker counts");
		}
		final NcrnaLearningVectors vectors=new NcrnaLearningVectors(population, family, tables, startRadius, stopRadius, shards, prefix);
		vectors.process(fasta, workers);
		return vectors;
	}

	private NcrnaLearningVectors(final NcrnaLearningPopulation population_, final NcrnaFamily family_,
			final NcrnaLearningTables tables_, final int startRadius_, final int stopRadius_, final int shards, final String prefix){
		population=population_;
		family=family_;
		tables=tables_;
		startRadius=startRadius_;
		stopRadius=stopRadius_;
		final int median=PackagedBoundaryConfig.upperMedianLength(family.library);
		if(startRadius<0 || stopRadius<0 || startRadius>median || stopRadius>median){
			throw new IllegalArgumentException("Candidate radii must be nonnegative and no larger than accepted-library upper median");
		}
		if(tables.statusByGene.length!=population.ids.length){throw new IllegalArgumentException("Context status must share the original locus axis");}
		trainPaths=new String[shards];
		for(int shard=0; shard<shards; shard++){trainPaths[shard]=prefix+".train."+shard+".tsv.gz";}
		validationPath=tables.refit ? null : prefix+".validation.tsv.gz";
		writers=new ByteStreamWriter[shards+(tables.refit ? 0 : 1)];
		rows=new long[writers.length];
		positives=new long[writers.length];
		shardBases=new long[shards];
		shardByGene=new int[population.ids.length];
		Arrays.fill(shardByGene, -1);
		seen=new boolean[population.ids.length];
		assignShards();
	}

	/** Greedy longest-first assignment with primitive packed length/index keys. */
	private void assignShards(){
		final long[] order=new long[population.ids.length];
		int count=0, validation=0;
		for(int gene=0; gene<order.length; gene++){
			if(!population.eligible(gene) || tables.statusByGene[gene]!=NcrnaLearningContext.AVAILABLE){continue;}
			if(!tables.refit && population.split().validationByGene[gene]){
				shardByGene[gene]=trainPaths.length;
				validation++;
			}else{order[count++]=((long)population.length[gene]<<32)|(gene&0xffffffffL);}
		}
		if(count<trainPaths.length || (!tables.refit && validation==0)){
			throw new IllegalArgumentException("Not enough usable loci for the requested training shards and validation set");
		}
		Arrays.sort(order, 0, count);
		for(int i=count-1; i>=0; i--){
			final int gene=(int)order[i];
			int shard=0;
			for(int candidate=1; candidate<shardBases.length; candidate++){
				if(shardBases[candidate]<shardBases[shard]){shard=candidate;}
			}
			shardByGene[gene]=shard;
			shardBases[shard]+=population.length[gene];
		}
	}

	/** Starts headers before consumers, then closes readers and writers after all workers finish. */
	private void process(final String fasta, final int workers){
		final FileFormat ff=FileFormat.testInput(fasta, FileFormat.FASTA, null, true, true);
		if(!ff.fasta()){throw new IllegalArgumentException("Boundary contexts must be FASTA");}
		final boolean trim=Shared.TRIM_READ_DESCRIPTION, upper=Read.TO_UPPER_CASE;
		final boolean bgzf=ReadWrite.USE_BGZF, force=ReadWrite.FORCE_BGZIP;
		Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=true;
		ReadWrite.USE_BGZF=true;
		ReadWrite.FORCE_BGZIP=true;
		try(NcrnaWorkers cohort=new NcrnaWorkers(); Outputs outputs=new Outputs(writers.length)){
			for(int shard=0; shard<writers.length; shard++){
				final String path=shard<trainPaths.length ? trainPaths[shard] : validationPath;
				writers[shard]=outputs.open(FileFormat.testOutput(path, FileFormat.TEXT, null, true, false, false, false));
				writers[shard].addJob(new ByteBuilder().append("#dims\t10\t1\n"));
			}
			try(SequenceInput owner=new SequenceInput(StreamerFactory.makeStreamer(ff, null, false, -1, false, true, 1))){
				final Streamer input=owner.stream;
				final ArrayList<ProcessThread> threads=new ArrayList<ProcessThread>(workers);
				for(int i=0; i<workers; i++){threads.add(new ProcessThread(input, this));}
				input.start();
				cohort.startAndJoin(threads, input);
				for(final ProcessThread worker:threads){accumulate(worker);}
			}
		}finally{
			Shared.TRIM_READ_DESCRIPTION=trim;
			Read.TO_UPPER_CASE=upper;
			ReadWrite.USE_BGZF=bgzf;
			ReadWrite.FORCE_BGZIP=force;
		}
		if(errorState){throw new IllegalStateException("Boundary vector stream did not finish cleanly");}
		if(records!=population.ids.length){throw new IllegalArgumentException("Every original locus needs one context in the vector pass");}
		for(int shard=0; shard<rows.length; shard++){
			if(rows[shard]==0){throw new IllegalArgumentException("No vectors in output file "+shard+"; fewer training shards may be needed");}
		}
		assert(records==filtered+unavailable+noModel+processed)
			: "Each context has one disposition: records="+records+", filtered="+filtered+", unavailable="+unavailable+", noModel="+noModel+", processed="+processed;
		assert(candidates==outcomes[0]+outcomes[1]+outcomes[2]+outcomes[3])
			: "Each attempted offset is emitted or rejected once: candidates="+candidates+", outcomes="+Arrays.toString(outcomes);
	}

	/** Merges only completed worker-local counters and detects duplicated source loci. */
	private void accumulate(final ProcessThread worker){
		errorState|=!worker.success;
		for(int i=0; i<worker.genes.size; i++){
			final int gene=worker.genes.get(i);
			if(seen[gene]){throw new IllegalArgumentException("Duplicate vector context: "+population.ids[gene]);}
			seen[gene]=true;
		}
		records+=worker.genes.size;
		filtered+=worker.filtered;
		unavailable+=worker.unavailable;
		noModel+=worker.noModel;
		processed+=worker.processed;
		candidates+=worker.candidates;
		for(int i=0; i<rows.length; i++){
			rows[i]+=worker.rows[i];
			positives[i]+=worker.positives[i];
		}
		for(int i=0; i<outcomes.length; i++){outcomes[i]+=worker.outcomes[i];}
	}

	/** All per-candidate mutation is local; only completed output blocks cross threads. */
	static final class ProcessThread extends Thread{
		ProcessThread(final Streamer input_, final NcrnaLearningVectors owner_){
			input=input_;
			owner=owner_;
			examples=new NcrnaBoundaryExamples(owner.family, owner.tables.start, owner.tables.stop);
			buffers=new ByteBuilder[owner.writers.length];
			rows=new long[buffers.length];
			positives=new long[buffers.length];
		}

		@Override
		public void run(){
			try{
				for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
					for(Read read:batch){processRead(read);}
				}
				for(int shard=0; shard<buffers.length; shard++){flush(shard);}
				success=true;
			}catch(Throwable failure){
				// A failed consumer otherwise strands bounded sequence/output queues.
				KillSwitch.exceptionKill(failure);
			}
		}

		private void processRead(final Read read){
			context.set(read, owner.population);
			genes.add(context.gene);
			if(context.status!=owner.tables.statusByGene[context.gene]){
				throw new IllegalArgumentException("Context availability changed between table and vector passes: "+read.id);
			}
			if(!owner.population.eligible(context.gene)){
				filtered++;
				return;
			}
			if(!context.hasFeatures()){
				unavailable++;
				return;
			}
			final int shard=owner.shardByGene[context.gene];
			assert(shard>=0) : "Every eligible context with features was assigned to one output file: gene="+
				context.gene+", shard="+shard+", outputs="+owner.writers.length;
			if(!examples.prepare(context, shard<owner.trainPaths.length)){
				noModel++;
				return;
			}
			processed++;
			emit(shard, false, owner.startRadius);
			emit(shard, true, owner.stopRadius);
		}

		private void emit(final int shard, final boolean isStop, final int radius){
			for(int offset=-radius; offset<=radius; offset++){
				candidates++;
				final int outcome=examples.fill(isStop, offset);
				outcomes[outcome]++;
				if(outcome!=NcrnaBoundaryExamples.VALID){continue;}
				if(buffers[shard]==null){buffers[shard]=new ByteBuilder(BLOCK_BYTES);}
				final ByteBuilder b=buffers[shard];
				// Preserve the bounded boundary-feature floats through native decimal readers.
				// Six decimals rounded valid GC/enrichment/coverage values before training.
				for(float value:examples.values){b.append(value, 18).tab();}
				b.append(examples.label).nl();
				rows[shard]++;
				positives[shard]+=examples.label;
				if(b.length>=BLOCK_BYTES){flush(shard);}
			}
		}

		private void flush(final int shard){
			final ByteBuilder b=buffers[shard];
			if(b!=null && b.length>0){
				owner.writers[shard].addJob(b);
				buffers[shard]=null;
			}
		}

		final Streamer input;
		final NcrnaLearningVectors owner;
		final NcrnaBoundaryExamples examples;
		final NcrnaLearningContext context=new NcrnaLearningContext();
		final IntList genes=new IntList();
		final ByteBuilder[] buffers;
		final long[] rows, positives;
		final long[] outcomes=new long[4];
		long filtered, unavailable, noModel, processed, candidates;
		boolean success;
	}

	/** One output path per length-balanced training shard. */
	public final String[] trainPaths;
	/** Held-out monitoring output, or null during all-eligible refit. */
	public final String validationPath;
	/** Original locus -> file index; -1 for excluded or unavailable contexts. */
	public final int[] shardByGene;
	/** Eligible gene bases assigned to each training shard before model screening. */
	public final long[] shardBases;
	/** Written vectors by output file; validation, if present, follows the training shards. */
	public final long[] rows;
	/** Positive endpoint labels by the same output-file axis as rows. */
	public final long[] positives;
	/** Candidate counts indexed by NcrnaBoundaryExamples outcome code. */
	public final long[] outcomes=new long[4];
	/** All original contexts read, including excluded or unavailable loci. */
	public long records;
	/** Original loci excluded by the population's eligibility policy. */
	public long filtered;
	/** Eligible loci lacking usable span geometry or whole-contig GC. */
	public long unavailable;
	/** Eligible feature-complete loci with no shortlisted model. */
	public long noModel;
	/** Loci with a selected model, before individual offset clipping. */
	public long processed;
	/** Attempted endpoint offsets, including candidates rejected for span geometry. */
	public long candidates;
	/** Borrowed immutable source identities and partition membership. */
	private final NcrnaLearningPopulation population;
	/** Borrowed accepted consensus/HBM pairs and model-selection policy. */
	private final NcrnaFamily family;
	/** Borrowed learned counts and complete original-context availability. */
	private final NcrnaLearningTables tables;
	/** Inclusive symmetric candidate radii in bases. */
	private final int startRadius, stopRadius;
	/** Original-locus presence, merged only after workers join. */
	private final boolean[] seen;
	/** Writers owned by process scopes and borrowed by workers until they join. */
	private final ByteStreamWriter[] writers;
	/** Set if any joined worker failed to finish. */
	private boolean errorState;
	/** Flush threshold for each worker's owned output block, in bytes. */
	private static final int BLOCK_BYTES=65536;
}
