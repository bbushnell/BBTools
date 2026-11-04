package prok;

import java.util.ArrayList;
import java.util.Arrays;

import fileIO.FileFormat;
import prok.NcrnaTextIO.SequenceInput;
import prok.TrnaBoundaryFeatures.BoundaryType;
import shared.KillSwitch;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;

/** Joins the complete context FASTA and counts both endpoint tables in one pass.
 * Adapted from A_SampleStreamerMT: workers pull native batches, own their maps
 * and parsers, and merge only after all workers join. Context bases are never
 * retained. Excluded/validation records still participate in the complete join.
 * Tuning learns eligible training records; refit learns all eligible records.
 * Changes reader header/case statics during the call; do not invoke concurrently
 * with another pipeline changing those process-wide settings.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaLearningTables{

	/** Reads every original context, failing on unknown, duplicated or missing loci. */
	public static NcrnaLearningTables build(final NcrnaLearningPopulation population, final String fasta,
			final int startInside, final int startOutside, final int stopInside, final int stopOutside,
			final boolean refit, final int workers){
		if(population==null || workers<1){throw new IllegalArgumentException("Context counting needs a population and positive worker count");}
		final NcrnaLearningTables tables=new NcrnaLearningTables(population, startInside, startOutside, stopInside, stopOutside, refit);
		tables.process(fasta, workers);
		return tables;
	}

	private NcrnaLearningTables(final NcrnaLearningPopulation population_, final int si, final int so,
			final int ei, final int eo, final boolean refit_){
		population=population_;
		refit=refit_;
		start=new NcrnaEndpointTable(BoundaryType.START, si, so);
		stop=new NcrnaEndpointTable(BoundaryType.STOP, ei, eo);
		statusByGene=new byte[population.ids.length];
		Arrays.fill(statusByGene, (byte)-1);
	}

	/** Owns reader configuration/lifetime; results become usable only after complete joins. */
	private void process(final String fasta, final int workers){
		final FileFormat ff=FileFormat.testInput(fasta, FileFormat.FASTA, null, true, true);
		if(!ff.fasta()){throw new IllegalArgumentException("Learning contexts must be FASTA: "+fasta);}
		final boolean trim=Shared.TRIM_READ_DESCRIPTION, upper=Read.TO_UPPER_CASE;
		Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=true;
		try(NcrnaWorkers cohort=new NcrnaWorkers();
				SequenceInput owner=new SequenceInput(StreamerFactory.makeStreamer(ff, null, false, -1, false, true, 1))){
			final Streamer input=owner.stream;
			final ArrayList<ProcessThread> threads=new ArrayList<ProcessThread>(workers);
			for(int i=0; i<workers; i++){threads.add(new ProcessThread(input, this));}
			input.start();
			cohort.startAndJoin(threads, input);
			for(final ProcessThread worker:threads){accumulate(worker);}
		}finally{
			Shared.TRIM_READ_DESCRIPTION=trim;
			Read.TO_UPPER_CASE=upper;
		}
		if(errorState){throw new IllegalStateException("Context reading/counting failed: "+fasta);}
		if(records!=statusByGene.length){throw new IllegalArgumentException("Expected "+statusByGene.length+" original contexts, observed "+records);}
		for(int gene=0; gene<statusByGene.length; gene++){
			if(statusByGene[gene]<0){throw new IllegalArgumentException("Missing original context: "+population.ids[gene]);}
		}
		if(start.counted()==0 || stop.counted()==0){throw new IllegalArgumentException("Each learned endpoint table needs at least one valid word");}
		assert(start.counted()+start.missing()+start.ambiguous()==learning
			&& stop.counted()+stop.missing()+stop.ambiguous()==learning)
			: "Each eligible in-scope context supplies one counted word or explicit sequence limitation per end; learning="+learning
			+", start counted/missing/ambiguous="+start.counted()+"/"+start.missing()+"/"+start.ambiguous()
			+", stop counted/missing/ambiguous="+stop.counted()+"/"+stop.missing()+"/"+stop.ambiguous();
	}

	/** Serial accumulation after joins; shared status slots need no per-read synchronization. */
	private void accumulate(final ProcessThread worker){
		errorState|=!worker.success;
		for(int i=0; i<worker.genes.size; i++){
			final int gene=worker.genes.get(i);
			if(statusByGene[gene]>=0){throw new IllegalArgumentException("Duplicate original context: "+population.ids[gene]);}
			statusByGene[gene]=worker.statuses.array[i];
		}
		start.merge(worker.start);
		stop.merge(worker.stop);
		records+=worker.genes.size;
		learning+=worker.learning;
		for(int status=0; status<availability.length; status++){availability[status]+=worker.availability[status];}
	}

	/** Owns native counts, parsed-record scratch and the indices needed for the final join. */
	static final class ProcessThread extends Thread{
		ProcessThread(final Streamer input_, final NcrnaLearningTables owner){
			input=input_;
			population=owner.population;
			refit=owner.refit;
			start=new NcrnaEndpointTable(BoundaryType.START, owner.start.inside, owner.start.outside);
			stop=new NcrnaEndpointTable(BoundaryType.STOP, owner.stop.inside, owner.stop.outside);
		}

		@Override
		public void run(){
			try{
				for(ListNum<Read> batch=input.nextList(); batch!=null; batch=input.nextList()){
					for(Read read:batch){processRead(read);}
				}
				success=true;
			}catch(Throwable failure){
				// A dead consumer could leave the FASTA producer blocked on its bounded queue.
				KillSwitch.exceptionKill(failure);
			}
		}

		private void processRead(final Read read){
			context.set(read, population);
			genes.add(context.gene);
			statuses.append((byte)context.status);
			availability[context.status]++;
			if(!population.eligible(context.gene) || (!refit && population.split().validationByGene[context.gene])){return;}
			learning++;
			if(context.hasGeometry()){
				start.add(read.bases, context.start);
				stop.add(read.bases, context.stop);
			}else{
				start.add(null, 0);
				stop.add(null, 0);
			}
		}

		final Streamer input;
		final NcrnaLearningPopulation population;
		final boolean refit;
		final NcrnaEndpointTable start, stop;
		final NcrnaLearningContext context=new NcrnaLearningContext();
		final IntList genes=new IntList();
		final ByteBuilder statuses=new ByteBuilder();
		final long[] availability=new long[4];
		long learning;
		boolean success;
	}

	/** Completed true-endpoint counts; consumers borrow them without mutation. */
	public final NcrnaEndpointTable start, stop;
	/** True when every eligible locus was counted, including the tuning holdout. */
	public final boolean refit;
	/** Counts over ALL original contexts, indexed by NcrnaLearningContext status. */
	public final long[] availability=new long[4];
	/** Every original context read, including excluded and validation records. */
	public long records;
	/** Eligible records assigned to the training population for this pass. */
	public long learning;
	/** Original-locus availability retained for balancing only vector-eligible loci. */
	public final byte[] statusByGene;
	private final NcrnaLearningPopulation population;
	private boolean errorState;
}
