package prok;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import fileIO.FileFormat;
import map.LongHashSet;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.SequenceInput;
import shared.Shared;
import shared.Timer;
import stream.Read;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Plans whole-genome CM search batches and writes one task's seekable FASTAs.
 * Uses the single-consumer template.A_SampleStreamer reader/writer pattern.
 * Planning reads census metadata only; writing preserves separate contigs and
 * their sequence bytes/orientation, replacing names with input-qualified IDs.
 * Literal target tokens and lengths are retained in a reversible sidecar.
 * No CM search, reference classification or source-file mutation occurs.
 * Output directories must be new. Completion summaries are written last.
 * @author Brian Bushnell, Raiden
 */
public final class PrepareCmSearchBatches{

	/** Runs one metadata-planning or sequence-writing invocation. */
	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PrepareCmSearchBatches tool=new PrepareCmSearchBatches(args);
		try{
			tool.process();
			timer.stop();
			tool.outstream.println("Batches: "+tool.batches.size()+"; time: "+timer);
		}finally{
			Shared.closeStream(tool.outstream);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Accepts standard BBTools options, but never filtering, pairing or overwriting. */
	public PrepareCmSearchBatches(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		args=pp.args;
		outstream=pp.outstream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase();
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("mode")){
				require("plan".equals(b) || "write".equals(b) || "schedule".equals(b), "mode must be plan, write or schedule");
				write=b.equals("write");
				schedule=b.equals("schedule");
			}else if(a.equals("modelcosts")){
				modelCosts=b;
			}else if(a.equals("samplebases")){
				sampleBases=Parse.parseKMG(b);
			}else if(a.equals("workers")){
				workers=Integer.parseInt(b);
			}else if(a.equals("taskperbatch")){
				taskPerBatch=Parse.parseBoolean(b);
			}else if(a.equals("concurrenttasks")){
				concurrentTasks=Integer.parseInt(b);
			}else if(a.equals("extraseconds")){
				extraSeconds=Double.parseDouble(b);
			}else if(a.equals("tasks")){
				tasks=Integer.parseInt(b);
			}else if(a.equals("task")){
				task=Integer.parseInt(b);
			}else if(a.equals("maxgenomes")){
				maxGenomes=Integer.parseInt(b);
			}else if(a.equals("maxbases")){
				maxBases=Parse.parseKMG(b);
			}else if(a.equals("include")){
				include=b;
			}else if(parser.parse(arg, a, b)){
				// Standard native options.
			}else{
				throw new IllegalArgumentException("Unknown parameter: "+arg);
			}
		}
		input=parser.in1;
		output=parser.out1;
		require(input!=null && output!=null, "Required: in=census_or_plan_directory out=new_directory");
		require(!parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1, "Whole genomes require complete unpaired inputs and fresh output");
		require(tasks>0 && maxGenomes>0 && maxBases>0, "Task count and batch limits must be positive");
		require(write ? task>=0 && include==null : task<0, "task= is required only for write; include= is planning-only");
		require(schedule ? modelCosts!=null && (sampleBases==-1 || sampleBases>0) && workers>0 && include==null && concurrentTasks>=0
			&& Double.isFinite(extraSeconds) && extraSeconds>=0
			: modelCosts==null && sampleBases==-1 && workers==8 && extraSeconds==0 && !taskPerBatch && concurrentTasks==0,
			"schedule requires modelcosts= and positive workers=; scheduling options do not apply to other modes");
		require(!new File(output).exists(), "Output directory already exists: "+output);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Planning            ----------------*/
	/*--------------------------------------------------------------*/

	/** Plans before creating output, or validates a complete saved plan before writing. */
	private void process(){
		if(write || schedule){
			loadPlan();
		}else{
			plan();
		}
		require(new File(output).mkdirs(), "Cannot create fresh output directory: "+output);
		if(write){writeTask();}else if(schedule){writeSchedule();}else{writePlan();}
	}

	/** Balances bases across tasks, then packs whole genomes into bounded task-local batches. */
	private void plan(){
		final ObjectSet<String> selected=include==null ? null : readSelection();
		final ArrayList<Genome> genomes=new ArrayList<Genome>();
		final LongHashSet ids=new LongHashSet();
		final ObjectSet<String> paths=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(input)){
			RrnaResourceIO.header(in.reader, NcrnaEvidenceCatalog.CENSUS_HEADER, input);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==5, "Census requires five fields");
				final Genome genome=genome(p, 0);
				require(ids.add(genome.id) && paths.add(genome.path), "Repeated input ID or FASTA path in census");
				if(selected==null || selected.contains(genome.path)){genomes.add(genome);}
			}
		}
		require(!genomes.isEmpty() && (selected==null || selected.size()==genomes.size()), "Selection must match nonempty, unique census paths exactly");
		require(tasks<=genomes.size(), "Do not create empty tasks; reduce tasks= for this population");
		Collections.sort(genomes, new Comparator<Genome>(){
			@Override
			public int compare(final Genome a, final Genome b){
				final int size=Long.compare(b.bases, a.bases);
				return size!=0 ? size : Long.compare(a.id, b.id);
			}
		});
		final long[] load=new long[tasks];
		final Batch[] last=new Batch[tasks];
		for(final Genome genome:genomes){
			require(genome.bases<=maxBases, "A whole genome exceeds maxbases; do not split it silently: "+genome.id);
			int worker=0;
			for(int i=1; i<tasks; i++){if(load[i]<load[worker]){worker=i;}}
			Batch batch=last[worker];
			if(batch==null || batch.genomes.size()>=maxGenomes || genome.bases>maxBases-batch.bases){
				batch=new Batch(batches.size(), worker);
				batches.add(batch);
				last[worker]=batch;
			}
			batch.add(genome);
			load[worker]=Math.addExact(load[worker], genome.bases);
		}
	}

	/** Optional timing subset: first literal TSV field is an absolute FASTA path; no header. */
	private ObjectSet<String> readSelection(){
		final ObjectSet<String> selected=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(include)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				final String path=text(p, 0);
				require(new File(path).isAbsolute() && selected.add(path), "Selection needs distinct absolute FASTA paths");
			}
		}
		return selected;
	}

	/** Writes a single metadata plan and explicit task totals; no sequence file is opened. */
	private void writePlan(){
		final ByteBuilder b=new ByteBuilder();
		final long[] taskBases=new long[tasks], taskGenomes=new long[tasks], taskBatches=new long[tasks];
		long genomes=0, bases=0, records=0;
		try(Output members=new Output(out("members.tsv")); Output sizes=new Output(out("batches.tsv"))){
			members.writer.println(MEMBER_HEADER);
			sizes.writer.println("batch_id\ttask_id\tgenomes\tnucleotide_bases\tfasta_records");
			for(final Batch batch:batches){
				for(final Genome genome:batch.genomes){
					b.clear().append(batch.id).tab().append(batch.task).tab();
					members.writer.print(genome.append(b).nl());
				}
				sizes.writer.print(b.clear().append(batch.id).tab().append(batch.task).tab().append(batch.genomes.size()).tab()
					.append(batch.bases).tab().append(batch.records).nl());
				taskBatches[batch.task]++;
				taskGenomes[batch.task]+=batch.genomes.size();
				taskBases[batch.task]=Math.addExact(taskBases[batch.task], batch.bases);
				genomes+=batch.genomes.size();
				bases=Math.addExact(bases, batch.bases);
				records=Math.addExact(records, batch.records);
			}
		}
		try(Output out=new Output(out("tasks.tsv"))){
			out.writer.println("task_id\tbatches\tgenomes\tnucleotide_bases");
			for(int i=0; i<tasks; i++){
				out.writer.print(b.clear().append(i).tab().append(taskBatches[i]).tab().append(taskGenomes[i]).tab().append(taskBases[i]).nl());
			}
		}
		try(Output out=new Output(out("summary.tsv"))){
			out.writer.println(PLAN_HEADER);
			out.writer.print(b.clear().append(tasks).tab().append(maxGenomes).tab().append(maxBases).tab().append(batches.size()).tab()
				.append(genomes).tab().append(bases).tab().append(records).append("\tMETADATA_ONLY_NO_FASTA_WRITTEN\n"));
		}
	}

	/** Requires the plan's entire membership population, not only the selected task. */
	private void loadPlan(){
		final LineParser1 p=new LineParser1('\t');
		long expectedGenomes, expectedBases, expectedRecords;
		int expectedBatches;
		try(Input in=new Input(new File(input, "summary.tsv").getPath())){
			RrnaResourceIO.header(in.reader, PLAN_HEADER, input);
			final byte[] row=in.reader.nextLine();
			require(row!=null, "Missing completed plan summary");
			p.set(row);
			require(p.terms()==8 && p.termEquals("METADATA_ONLY_NO_FASTA_WRITTEN", 7) && in.reader.nextLine()==null, "Invalid plan summary");
			tasks=p.parseInt(0);
			maxGenomes=p.parseInt(1);
			maxBases=p.parseLong(2);
			expectedBatches=p.parseInt(3);
			expectedGenomes=p.parseLong(4);
			expectedBases=p.parseLong(5);
			expectedRecords=p.parseLong(6);
		}
		require(tasks>0 && task<tasks && maxGenomes>0 && maxBases>0 && expectedBatches>0, "Invalid saved plan dimensions or selected task");
		final LongHashSet ids=new LongHashSet();
		final ObjectSet<String> paths=new ObjectSet<String>(String.class);
		long genomes=0, bases=0, records=0;
		Batch batch=null;
		try(Input in=new Input(new File(input, "members.tsv").getPath())){
			RrnaResourceIO.header(in.reader, MEMBER_HEADER, input);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==7, "Batch membership requires seven fields");
				final int id=p.parseInt(0), worker=p.parseInt(1);
				if(batch==null || batch.id!=id){
					require(id==batches.size() && worker>=0 && worker<tasks, "Plan batches must be dense, grouped and assigned to valid tasks");
					batch=new Batch(id, worker);
					batches.add(batch);
				}
				require(batch.task==worker, "A batch cannot belong to multiple tasks");
				final Genome genome=genome(p, 2);
				require(ids.add(genome.id) && paths.add(genome.path), "Repeated planned input ID or path");
				batch.add(genome);
				require(batch.genomes.size()<=maxGenomes && batch.bases<=maxBases, "Saved batch exceeds its declared whole-genome limits");
				genomes++;
				bases=Math.addExact(bases, genome.bases);
				records=Math.addExact(records, genome.records);
			}
		}
		require(batches.size()==expectedBatches && genomes==expectedGenomes && bases==expectedBases && records==expectedRecords,
			"Membership must conserve the complete saved plan population");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Modeled Scheduling     ----------------*/
	/*--------------------------------------------------------------*/

	/** Predicts longest-processing-time-first queues from measured one-worker CPU.
	 * The model assumes linear base scaling and no scheduler, I/O or contention
	 * delay. Tasks start in ID order when a modeled task slot becomes free.
	 * Extra seconds are an explicit unmeasured sensitivity
	 * allowance beyond the scaled pilot; startup already present in the pilot is
	 * not separately recovered or added again. Jobs are retained output records.
	 */
	private void writeSchedule(){
		final ArrayList<ModelCost> models=readCosts();
		final int taskCount=taskPerBatch ? batches.size() : tasks;
		final int parallel=concurrentTasks==0 ? taskCount : Math.min(concurrentTasks, taskCount);
		final double[] taskSlots=new double[parallel];
		long commonBases=models.get(0).bases;
		for(final ModelCost model:models){if(model.bases!=commonBases){commonBases=0;}}
		final ByteBuilder b=new ByteBuilder();
		double totalCpu=0, fleetWall=0;
		long searches=0;
		try(Output jobsOut=new Output(out("jobs.tsv")); Output tasksOut=new Output(out("tasks.tsv"));
			Output workersOut=new Output(out("workers.tsv")); Output startsOut=new Output(out("task_starts.tsv"))){
			startsOut.writer.println("task_id\tconcurrent_task_slots\tmodeled_slot\tmodeled_start_seconds\tmodeled_end_seconds");
			jobsOut.writer.println("task_id\trank\tbatch_id\tmodel\tinput_bases\tlinear_CPU_seconds\textra_CPU_seconds\tmodeled_CPU_seconds\tmodeled_worker\tmodeled_start_seconds\tmodeled_end_seconds");
			tasksOut.writer.println("task_id\tjobs\tworker_slots\tlinear_CPU_seconds\textra_CPU_seconds\tmodeled_CPU_seconds\tideal_wall_seconds\tmodeled_makespan_seconds\tidle_capacity_seconds\ttail_over_ideal_seconds\tlongest_point_seconds\tstatus");
			workersOut.writer.println("task_id\tworker\tmodeled_busy_seconds\tidle_until_task_end_seconds");
			for(int workerTask=0; workerTask<taskCount; workerTask++){
				final ArrayList<SearchJob> jobs=new ArrayList<SearchJob>();
				for(final Batch batch:batches){
					if((taskPerBatch ? batch.id : batch.task)!=workerTask){continue;}
					for(final ModelCost model:models){
						final double linear=model.cpu*batch.bases/model.bases;
						require(Double.isFinite(linear) && linear>0 && Double.isFinite(linear+extraSeconds), "Modeled search cost must remain finite and positive");
						jobs.add(new SearchJob(batch, model, linear, linear+extraSeconds));
					}
				}
				require(!jobs.isEmpty(), "Every declared task must own modeled searches");
				Collections.sort(jobs, JOB_ORDER);
				final double[] load=new double[workers];
				double linearCpu=0, modeledCpu=0;
				int rank=0;
				for(final SearchJob job:jobs){
					int pick=0;
					for(int i=1; i<workers; i++){if(load[i]<load[pick]){pick=i;}}
					final double start=load[pick], end=start+job.cpu;
					require(Double.isFinite(end), "Modeled worker time overflowed");
					load[pick]=end;
					linearCpu+=job.linear;
					modeledCpu+=job.cpu;
					b.clear().append(workerTask).tab().append(rank++).tab().append(job.batch.id).tab().append(job.model.id).tab().append(job.batch.bases).tab();
					jobsOut.writer.print(b.append(job.linear, 6).tab().append(extraSeconds, 6).tab().append(job.cpu, 6).tab().append(pick).tab()
						.append(start, 6).tab().append(end, 6).nl());
				}
				double wall=0, conserved=0;
				for(final double value:load){
					wall=Math.max(wall, value);
					conserved+=value;
				}
				assert(Math.abs(conserved-modeledCpu)<=1e-8*Math.max(1, modeledCpu)) : "Every modeled search must contribute once to exactly one worker load";
				final double ideal=modeledCpu/workers, idle=Math.max(0, wall*workers-modeledCpu);
				b.clear().append(workerTask).tab().append(jobs.size()).tab().append(workers).tab().append(linearCpu, 6).tab();
				tasksOut.writer.print(b.append(extraSeconds*jobs.size(), 6).tab().append(modeledCpu, 6).tab().append(ideal, 6).tab()
					.append(wall, 6).tab().append(idle, 6).tab().append(Math.max(0, wall-ideal), 6).tab().append(jobs.get(0).cpu, 6)
					.append("\tMODELED_NOT_MEASURED\n"));
				for(int i=0; i<workers; i++){
					workersOut.writer.print(b.clear().append(workerTask).tab().append(i).tab().append(load[i], 6).tab().append(wall-load[i], 6).nl());
				}
				totalCpu+=modeledCpu;
				int slot=0;
				for(int i=1; i<parallel; i++){if(taskSlots[i]<taskSlots[slot]){slot=i;}}
				final double start=taskSlots[slot], end=start+wall;
				require(Double.isFinite(end), "Modeled task-slot time overflowed");
				taskSlots[slot]=end;
				startsOut.writer.print(b.clear().append(workerTask).tab().append(parallel).tab().append(slot).tab()
					.append(start, 6).tab().append(end, 6).nl());
				fleetWall=Math.max(fleetWall, end);
				searches=Math.addExact(searches, jobs.size());
			}
		}
		require(searches==(long)batches.size()*models.size() && Double.isFinite(totalCpu), "Schedule must conserve every batch/model pair without overflow");
		final double slots=(double)parallel*workers;
		try(Output out=new Output(out("summary.tsv"))){
			out.writer.println("tasks\tworkers_per_task\tbatches\tmodels\tsearches\tsample_bases\textra_seconds_per_search\tmodeled_CPU_seconds\tideal_fleet_wall_seconds\tmodeled_fleet_wall_seconds\ttotal_idle_capacity_seconds\tstatus");
			out.writer.print(b.clear().append(taskCount).tab().append(workers).tab().append(batches.size()).tab().append(models.size()).tab()
				.append(searches).tab().append(commonBases).tab().append(extraSeconds, 6).tab().append(totalCpu, 6).tab()
				.append(totalCpu/slots, 6).tab().append(fleetWall, 6).tab().append(Math.max(0, fleetWall*slots-totalCpu), 6)
				.append("\tMODELED_NOT_MEASURED\n"));
		}
	}

	/** Reads legacy paired costs with one supplied sample size, or qualified costs
	 * with a measured base count per model. Ambiguous mixed schemas are rejected.
	 */
	private ArrayList<ModelCost> readCosts(){
		final ArrayList<ModelCost> result=new ArrayList<ModelCost>();
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(modelCosts)){
			final byte[] header=in.reader.nextLine();
			require(header!=null, "Model cost table is empty");
			p.set(header);
			final int fields=p.terms();
			int model=-1, cpu=-1, bases=-1, status=-1;
			boolean legacy=false;
			for(int i=0; i<fields; i++){
				if(p.termEquals("accession", i)){
					require(model<0, "Repeated accession column");
					model=i;
				}else if(p.termEquals("candidate_CPU_seconds", i) || p.termEquals("CPU_seconds", i)){
					require(cpu<0, "Repeated or ambiguous CPU column");
					cpu=i;
					legacy=p.termEquals("candidate_CPU_seconds", i);
				}else if(p.termEquals("sample_bases", i)){
					require(bases<0, "Repeated sample_bases column");
					bases=i;
				}else if(p.termEquals("status", i)){
					require(status<0, "Repeated status column");
					status=i;
				}
			}
			require(model>=0 && cpu>=0, "Costs need accession and a measured CPU column");
			require(legacy ? sampleBases>0 && bases<0 : sampleBases==-1 && bases>=0 && status>=0,
				"Legacy candidate_CPU_seconds needs samplebases=; CPU_seconds needs per-row sample_bases and status, without samplebases=");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==fields, "Model cost row differs from its header");
				final String id=p.parseString(model);
				final double seconds=p.parseDouble(cpu);
				final long measuredBases=legacy ? sampleBases : p.parseLong(bases);
				require(measuredBases>0, "Measured sample bases must be positive for linear scaling");
				require(legacy || p.termEquals("MEASURED", status) || p.termEquals("MEASURED_REUSED", status),
					"Scheduling requires a qualified measured cost for every model");
				require(id.matches("[A-Za-z0-9_.-]+") && ids.add(id), "Model identifiers must be distinct simple tokens");
				require(Double.isFinite(seconds) && seconds>=10, "Scheduling requires measured one-worker baselines of at least ten CPU seconds");
				result.add(new ModelCost(id, seconds, measuredBases));
			}
		}
		require(!result.isEmpty(), "No measured model costs supplied");
		return result;
	}

	private static final class ModelCost{
		ModelCost(final String id_, final double cpu_, final long bases_){
			id=id_;
			cpu=cpu_;
			bases=bases_;
		}
		final String id;
		final double cpu;
		final long bases;
	}
	private static final class SearchJob{
		SearchJob(final Batch batch_, final ModelCost model_, final double linear_, final double cpu_){
			batch=batch_;
			model=model_;
			linear=linear_;
			cpu=cpu_;
		}
		final Batch batch;
		final ModelCost model;
		final double linear, cpu;
	}
	private static final Comparator<SearchJob> JOB_ORDER=(a, b)->{
		int order=Double.compare(b.cpu, a.cpu);
		if(order==0){order=Integer.compare(a.batch.id, b.batch.id);}
		return order==0 ? a.model.id.compareTo(b.model.id) : order;
	};

	/*--------------------------------------------------------------*/
	/*----------------       Sequence Writing       ----------------*/
	/*--------------------------------------------------------------*/

	/** Preserves sequence bytes by disabling Read normalization during this invocation. */
	private void writeTask(){
		final boolean validate=Read.VALIDATE_IN_CONSTRUCTOR;
		final boolean fixHeader=Read.FIX_HEADER;
		Read.VALIDATE_IN_CONSTRUCTOR=false;
		Read.FIX_HEADER=false;
		long bases=0, records=0, genomes=0;
		int written=0;
		try{
			for(final Batch batch:batches){
				if(batch.task!=task){continue;}
				writeBatch(batch);
				written++;
				genomes+=batch.genomes.size();
				bases=Math.addExact(bases, batch.bases);
				records=Math.addExact(records, batch.records);
			}
		}finally{
			Read.VALIDATE_IN_CONSTRUCTOR=validate;
			Read.FIX_HEADER=fixHeader;
		}
		require(written>0, "Selected task has no batch; no preparation completion can be claimed");
		try(Output out=new Output(out("summary.tsv"))){
			out.writer.println("task_id\tbatches\tgenomes\tnucleotide_bases\tfasta_records\tstatus");
			out.writer.print(new ByteBuilder().append(task).tab().append(written).tab().append(genomes).tab()
				.append(bases).tab().append(records).append("\tPREPARED_NOT_SEARCHED\n"));
		}
	}

	/** One FASTA is shared by every model; it is not concatenated across contig boundaries. */
	private void writeBatch(final Batch batch){
		final File dir=new File(output, "batch_"+batch.id);
		require(dir.mkdir(), "Batch directory must be new: "+dir);
		final FileFormat format=FileFormat.testOutput(new File(dir, "batch.fa").getPath(), FileFormat.FASTA, null, true, false, false, true);
		final Writer writer=WriterFactory.makeWriter(format, null, 0, null, false);
		try(SequenceOutput owned=new SequenceOutput(writer);
			Output mapping=new Output(new File(dir, "targets.tsv").getPath());
			Output members=new Output(new File(dir, "members.tsv").getPath())){
			mapping.writer.println(TARGET_HEADER);
			members.writer.println(NcrnaEvidenceCatalog.CENSUS_HEADER);
			final ByteBuilder b=new ByteBuilder();
			long listId=0;
			for(final Genome genome:batch.genomes){
				listId=writeGenome(batch.id, genome, listId, writer, mapping);
				members.writer.print(genome.append(b.clear()).nl());
			}
		}
		require(writer.readsWritten()==batch.records && writer.basesWritten()==batch.bases, "FASTA writer must preserve all planned contigs and bases");
		try(Output out=new Output(new File(dir, "summary.tsv").getPath())){
			out.writer.println("batch_id\ttask_id\tgenomes\tnucleotide_bases\tfasta_records\tstatus");
			out.writer.print(new ByteBuilder().append(batch.id).tab().append(batch.task).tab().append(batch.genomes.size()).tab()
				.append(batch.bases).tab().append(batch.records).append("\tPREPARED_NOT_SEARCHED\n"));
		}
	}

	/** Streams a whole original input once. Any validation failure drains the reader
	 * before propagation; already-submitted batches remain writer-owned and immutable.
	 */
	private long writeGenome(final int batch, final Genome genome, long listId,
		final Writer writer, final Output mapping){
		require(new File(genome.path).length()==genome.bytes, "Source size changed since census: "+genome.path);
		final FileFormat format=FileFormat.testInput(genome.path, FileFormat.FASTA, null, true, true);
		require(format.fasta() && !format.interleaved(), "Batch inputs must be unpaired whole-genome FASTA: "+genome.path);
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		final ByteBuilder b=new ByteBuilder();
		long bases=0, records=0;
		RuntimeException failure=null;
		try(SequenceInput in=new SequenceInput(StreamerFactory.makeStreamer(format, null, true, -1, false, true, 0))){
			in.stream.start();
			for(ListNum<Read> list=in.stream.nextList(); list!=null; list=in.stream.nextList()){
				for(final Read read:list.list){
					if(failure!=null){continue;}
					try{
						if(read==null || read.bases==null || read.length()<1 || read.mate!=null){
							throw new IllegalArgumentException("Expected nonempty unpaired contig in input "+genome.id);
						}
						final String original=targetToken(read.id);
						if(!names.add(original)){throw new IllegalArgumentException("Repeated target token within input "+genome.id+": "+original);}
						records++;
						bases=Math.addExact(bases, read.length());
						read.id=b.clear().append("i").append(genome.id).append("c").append(records).toString();
						mapping.writer.print(b.clear().append(batch).tab().append(read.id).tab().append(genome.id).tab()
							.append(records).tab().append(original).tab().append(read.length()).nl());
					}catch(RuntimeException invalid){failure=invalid;}
				}
				if(failure==null){writer.add(list.list, listId++);}else{in.stream.returnList(list);}
			}
		}
		if(failure!=null){throw failure;}
		require(bases==genome.bases && records==genome.records, "Actual contigs/bases differ from planned input "+genome.id);
		return listId;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Shared Records        ----------------*/
	/*--------------------------------------------------------------*/

	/** Original identifiers use the native ASCII token rule; encoding loss is fatal. */
	static String targetToken(final String header){
		require(header!=null && !header.isEmpty(), "FASTA target identifier is empty");
		int end=0;
		while(end<header.length()){
			final char c=header.charAt(end);
			if(c==' ' || (c>='\t' && c<='\r')){break;}
			if(c<33 || c>126){throw new IllegalArgumentException("FASTA target token is not preserved printable ASCII");}
			end++;
		}
		require(end>0, "FASTA target begins with whitespace");
		return end==header.length() ? header : header.substring(0, end);
	}
	private static Genome genome(final LineParser1 p, final int offset){
		final long id=p.parseLong(offset), bytes=p.parseLong(offset+2), bases=p.parseLong(offset+3), records=p.parseLong(offset+4);
		final String path=text(p, offset+1);
		require(id>0 && p.termEquals(Long.toString(id), offset) && new File(path).isAbsolute()
			&& bytes>0 && bases>0 && records>0, "Census requires canonical positive IDs/dimensions and absolute paths");
		return new Genome(id, path, bytes, bases, records);
	}
	private static String text(final LineParser1 p, final int field){
		p.setBounds(field);
		return new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
	}
	private String out(final String file){return new File(output, file).getPath();}
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}
	private static final class Genome{
		Genome(final long id_, final String path_, final long bytes_, final long bases_, final long records_){
			id=id_;
			path=path_;
			bytes=bytes_;
			bases=bases_;
			records=records_;
		}
		ByteBuilder append(final ByteBuilder b){
			return b.append(id).tab().append(path.getBytes(StandardCharsets.UTF_8)).tab().append(bytes).tab().append(bases).tab().append(records);
		}
		final long id, bytes, bases, records;
		final String path;
	}
	private static final class Batch{
		Batch(final int id_, final int task_){
			id=id_;
			task=task_;
		}
		void add(final Genome genome){
			assert(genome!=null) : "Batch membership must retain each concrete census entry";
			genomes.add(genome);
			bases=Math.addExact(bases, genome.bases);
			records=Math.addExact(records, genome.records);
		}
		final int id, task;
		final ArrayList<Genome> genomes=new ArrayList<Genome>();
		long bases, records;
	}
	/** Owns the factory's caller-driven FASTA writer from startup through closure.
	 * FastqWriterST2.finishError only marks abandonment; it does not close its stream.
	 * Normal finalization is therefore required on validation failures too. The
	 * enclosing resource scope preserves the original failure and suppresses cleanup
	 * failures. A written prefix has no completion summary after a body failure.
	 */
	private static final class SequenceOutput implements AutoCloseable{
		SequenceOutput(final Writer writer_){
			writer=writer_;
			try{
				writer.start();
			}catch(RuntimeException | Error failure){
				try{close();}catch(RuntimeException | Error cleanup){
					if(cleanup!=failure){failure.addSuppressed(cleanup);}
				}
				throw failure;
			}
		}
		@Override
		public void close(){
			if(writer.poisonAndWait()){throw new IllegalStateException("Failed to finalize batch FASTA output");}
		}
		private final Writer writer;
	}
	private String input, output, include, modelCosts;
	private boolean write, schedule, taskPerBatch;
	private int workers=8;
	/** Zero means all tasks may start together; positive values cap modeled slots. */
	private int concurrentTasks=0;
	private long sampleBases=-1;
	private double extraSeconds=0;
	private int tasks=16, task=-1, maxGenomes=64;
	private long maxBases=400000000;
	private final ArrayList<Batch> batches=new ArrayList<Batch>();
	private PrintStream outstream=System.err;
	public static final String MEMBER_HEADER="batch_id\ttask_id\tinput_id\tfasta\tfile_bytes\tnucleotide_bases\tfasta_records";
	public static final String TARGET_HEADER="batch_id\ttarget_id\tinput_id\tcontig_ordinal\toriginal_target\tlength";
	public static final String PLAN_HEADER="tasks\tmax_genomes\tmax_bases\tbatches\tgenomes\tnucleotide_bases\tfasta_records\tstatus";
}
