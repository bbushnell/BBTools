package test.prok;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import parse.LineParser1;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.PrepareCmSearchBatches;
import structures.ByteBuilder;

/** Finite whole-genome batching cases with literal FASTA and map oracles.
 * No CM searches or production genome reads occur.
 * @author Brian Bushnell, Raiden
 */
public final class PrepareCmSearchBatchesTest{

	public static void main(final String[] args){
		String out=null;
		for(final String arg:args){if(arg.startsWith("out=")){out=arg.substring(4);}}
		check(out!=null && new File(out).mkdirs(), "Use out=fresh_fixture_directory");
		run(new File(out).getAbsolutePath());
		System.err.println("PREPARE_CM_SEARCH_BATCHES_TEST_PASS");
	}

	/** Same target in two inputs remains distinct; lowercase/IUPAC and contig boundaries survive. */
	private static void run(final String root){
		final String a=root+"/a.fa", b=root+"/b.fa", c=root+"/c.fa";
		write(a, ">shared description\nacgtNR\n>odd%25\nGGcc\n");
		write(b, ">shared another description\nTTaa\n");
		write(c, ">third\nNY\n");
		final String census=root+"/census.tsv";
		write(census, NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'+row(1, a, 10, 2)+row(2, b, 4, 1)+row(3, c, 2, 1));
		plan(census, root+"/plan", "tasks=2", "maxgenomes=2", "maxbases=12");
		invokeWrite(root+"/plan", root+"/task0", 0);
		invokeWrite(root+"/plan", root+"/task1", 1);
		equal(root+"/task0/batch_0/batch.fa", ">i1c1\nacgtNR\n>i1c2\nGGcc\n");
		equal(root+"/task1/batch_1/batch.fa", ">i2c1\nTTaa\n>i3c1\nNY\n");
		equal(root+"/task0/batch_0/targets.tsv", PrepareCmSearchBatches.TARGET_HEADER+'\n'
			+"0\ti1c1\t1\t1\tshared\t6\n0\ti1c2\t1\t2\todd%25\t4\n");
		equal(root+"/task1/batch_1/targets.tsv", PrepareCmSearchBatches.TARGET_HEADER+'\n'
			+"1\ti2c1\t2\t1\tshared\t4\n1\ti3c1\t3\t1\tthird\t2\n");
		equal(root+"/plan/summary.tsv", PrepareCmSearchBatches.PLAN_HEADER+"\n2\t2\t12\t2\t3\t16\t4\tMETADATA_ONLY_NO_FASTA_WRITTEN\n");
		equal(root+"/task1/summary.tsv", "task_id\tbatches\tgenomes\tnucleotide_bases\tfasta_records\tstatus\n1\t1\t2\t6\t2\tPREPARED_NOT_SEARCHED\n");
		// Repeated original tokens must also coexist inside the very same native batch.
		plan(census, root+"/same_batch", "tasks=1", "maxgenomes=3", "maxbases=16");
		invokeWrite(root+"/same_batch", root+"/same_output", 0);
		equal(root+"/same_output/batch_0/batch.fa", ">i1c1\nacgtNR\n>i1c2\nGGcc\n>i2c1\nTTaa\n>i3c1\nNY\n");
		equal(root+"/same_output/batch_0/targets.tsv", PrepareCmSearchBatches.TARGET_HEADER+'\n'
			+"0\ti1c1\t1\t1\tshared\t6\n0\ti1c2\t1\t2\todd%25\t4\n"
			+"0\ti2c1\t2\t1\tshared\t4\n0\ti3c1\t3\t1\tthird\t2\n");
		// Exercise each packing overflow independently, including writing several batches in one task.
		plan(census, root+"/base_limit", "tasks=1", "maxgenomes=3", "maxbases=12");
		equal(root+"/base_limit/summary.tsv", PrepareCmSearchBatches.PLAN_HEADER+"\n1\t3\t12\t2\t3\t16\t4\tMETADATA_ONLY_NO_FASTA_WRITTEN\n");
		plan(census, root+"/genome_limit", "tasks=1", "maxgenomes=1", "maxbases=100");
		invokeWrite(root+"/genome_limit", root+"/rotated_output", 0);
		equal(root+"/genome_limit/summary.tsv", PrepareCmSearchBatches.PLAN_HEADER+"\n1\t1\t100\t3\t3\t16\t4\tMETADATA_ONLY_NO_FASTA_WRITTEN\n");
		equal(root+"/rotated_output/batch_0/batch.fa", ">i1c1\nacgtNR\n>i1c2\nGGcc\n");
		equal(root+"/rotated_output/batch_1/batch.fa", ">i2c1\nTTaa\n");
		equal(root+"/rotated_output/batch_2/batch.fa", ">i3c1\nNY\n");
		// Selection remains explicit metadata, not an unlabelled shortened full population.
		write(root+"/include.tsv", b+"\tunused_metadata\n");
		plan(census, root+"/selected", "tasks=1", "include="+root+"/include.tsv");
		equal(root+"/selected/summary.tsv", PrepareCmSearchBatches.PLAN_HEADER+"\n1\t64\t400000000\t1\t1\t4\t1\tMETADATA_ONLY_NO_FASTA_WRITTEN\n");
		expectFailure(new Runnable(){public void run(){
			plan(census, root+"/oversize", "tasks=1", "maxbases=9");
		}}, "exceeds maxbases");
		check(!new File(root+"/oversize").exists(), "Oversize planning must fail before creating a success-shaped output");
		write(root+"/bad_ids.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'+row(1, a, 10, 2)+row(1, b, 4, 1));
		expectFailure(new Runnable(){public void run(){
			plan(root+"/bad_ids.tsv", root+"/bad_ids", "tasks=1");
		}}, "Repeated input");
		// Duplicate tokens within one genome remain unsupported, unlike repetition across genomes.
		final String duplicate=root+"/duplicate.fa";
		write(duplicate, ">same\nAC\n>same other description\nGT\n");
		write(root+"/duplicate.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'+row(4, duplicate, 4, 2));
		plan(root+"/duplicate.tsv", root+"/duplicate_plan", "tasks=1");
		expectFailure(new Runnable(){public void run(){
			invokeWrite(root+"/duplicate_plan", root+"/duplicate_output", 0);
		}}, "Repeated target token");
		check(!new File(root+"/duplicate_output/summary.tsv").exists()
			&& !new File(root+"/duplicate_output/batch_0/summary.tsv").exists(), "Failed preparation cannot publish completion");
		closed(root+"/duplicate_output/batch_0/batch.fa");
		// The source census may be internally consistent yet wrong about actual bases.
		write(root+"/wrong_bases.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'+row(2, b, 5, 1));
		plan(root+"/wrong_bases.tsv", root+"/wrong_plan", "tasks=1");
		expectFailure(new Runnable(){public void run(){
			invokeWrite(root+"/wrong_plan", root+"/wrong_output", 0);
		}}, "Actual contigs/bases differ");
		check(!new File(root+"/wrong_output/summary.tsv").exists(), "Wrong source dimensions cannot publish completion");
		closed(root+"/wrong_output/batch_0/batch.fa");
		schedules(root);
	}

	/** Literal synthetic timings expose chunk tails; nonexistent FASTAs must not be opened. */
	private static void schedules(final String root){
		final String census=root+"/schedule_census.tsv", costs=root+"/costs.tsv";
		write(census, NcrnaEvidenceCatalog.CENSUS_HEADER+"\n1\t"+root+"/absent_a.fa\t1\t10\t1\n2\t"
			+root+"/absent_b.fa\t1\t4\t1\n3\t"+root+"/absent_c.fa\t1\t2\t1\n");
		write(costs, "accession\tcandidate_CPU_seconds\nA\t40\nB\t10\n");
		plan(census, root+"/schedule_split", "tasks=1", "maxgenomes=1");
		plan(census, root+"/schedule_coarse", "tasks=1", "maxgenomes=3");
		plan(census, root+"/schedule_two_tasks", "tasks=2", "maxgenomes=3");
		schedule(root+"/schedule_split", root+"/queue_split", costs, 2, 0);
		schedule(root+"/schedule_coarse", root+"/queue_coarse", costs, 2, 0);
		schedule(root+"/schedule_split", root+"/queue_serial", costs, 1, 0);
		schedule(root+"/schedule_split", root+"/queue_extra", costs, 2, 1);
		schedule(root+"/schedule_two_tasks", root+"/queue_two_tasks", costs, 2, 0);
		// Split costs40,16,10,8,4,2 fill two workers40/40; coarse64/16 leaves48 idle.
		metrics(root+"/queue_split/summary.tsv", 6, 80, 40, 0);
		metrics(root+"/queue_coarse/summary.tsv", 2, 80, 64, 48);
		metrics(root+"/queue_serial/summary.tsv", 6, 80, 80, 0);
		// Added one second per invocation gives41,17,11,9,5,3 and loads44/42.
		metrics(root+"/queue_extra/summary.tsv", 6, 86, 44, 2);
		// Independent tasks have loads40/10 and24/6: fleet wall40, four slots,80 idle.
		metrics(root+"/queue_two_tasks/summary.tsv", 4, 80, 40, 80);
		final LineParser1 p=new LineParser1('\t');
		final int[] batchIds={0, 1, 0, 2, 1, 2};
		final String[] modelIds={"A", "A", "B", "A", "B", "B"};
		try(Input in=new Input(root+"/queue_split/jobs.tsv")){
			in.reader.nextLine();
			int rank=0;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				check(rank<6 && p.parseInt(1)==rank && p.parseInt(2)==batchIds[rank] && p.termEquals(modelIds[rank], 3),
					"Queues must retain every batch/model identity in expensive-first order");
				rank++;
			}
			check(rank==6, "Six distinct batch/model searches");
		}
		write(root+"/short_cost.tsv", "accession\tcandidate_CPU_seconds\nA\t9\n");
		expectFailure(new Runnable(){public void run(){
			schedule(root+"/schedule_split", root+"/queue_short", root+"/short_cost.tsv", 2, 0);
		}}, "at least ten CPU seconds");
		check(!new File(root+"/queue_short/summary.tsv").exists(), "Invalid costs cannot publish a schedule summary");
		write(root+"/duplicate_cost.tsv", "accession\tcandidate_CPU_seconds\nA\t40\nA\t10\n");
		expectFailure(new Runnable(){public void run(){
			schedule(root+"/schedule_split", root+"/queue_duplicate", root+"/duplicate_cost.tsv", 2, 0);
		}}, "distinct simple tokens");
		mixedSchedules(root);
	}
	/** Mixed sample sizes must scale independently; batches retain their IDs when
	 * assigned to separate tasks. Literal starts expose a capped-array tail.
	 */
	private static void mixedSchedules(final String root){
		final String costs=root+"/mixed_costs.tsv";
		write(costs, "accession\tCPU_seconds\tsample_bases\tstatus\nA\t40\t10\tMEASURED_REUSED\nB\t20\t20\tMEASURED\n");
		mixedSchedule(root, root+"/queue_mixed", costs, "concurrenttasks=2");
		metrics(root+"/queue_mixed/summary.tsv", 6, 80, 40, 80);
		equal(root+"/queue_mixed/task_starts.tsv", "task_id\tconcurrent_task_slots\tmodeled_slot\tmodeled_start_seconds\tmodeled_end_seconds\n"
			+"0\t2\t0\t0\t40\n1\t2\t1\t0\t16\n2\t2\t1\t16\t24\n");
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(root+"/queue_mixed/jobs.tsv")){
			in.reader.nextLine();
			int rows=0;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				check(p.parseInt(0)==rows/2 && p.parseInt(2)==rows/2 && p.parseInt(1)==rows%2,
					"Each original batch must occupy its own task without renumbering or dropping searches");
				rows++;
			}
			check(rows==6, "Every original batch/model pair must remain present");
		}
		try(Input in=new Input(root+"/queue_mixed/summary.tsv")){
			in.reader.nextLine();
			p.set(in.reader.nextLine());
			check(p.parseInt(0)==3 && p.parseLong(5)==0, "Summary must declare three tasks and mixed sample bases");
		}
		mixedSchedule(root, root+"/queue_capped_serial", costs, "concurrenttasks=1");
		metrics(root+"/queue_capped_serial/summary.tsv", 6, 80, 64, 48);
		write(root+"/unqualified.tsv", "accession\tCPU_seconds\tsample_bases\tstatus\nA\t40\t10\tUNPRICED\n");
		expectFailure(new Runnable(){public void run(){
			mixedSchedule(root, root+"/queue_unqualified", root+"/unqualified.tsv", "concurrenttasks=1");
		}}, "qualified measured cost");
		write(root+"/zero_sample.tsv", "accession\tCPU_seconds\tsample_bases\tstatus\nA\t40\t0\tMEASURED\n");
		expectFailure(new Runnable(){public void run(){
			mixedSchedule(root, root+"/queue_zero_sample", root+"/zero_sample.tsv", "concurrenttasks=1");
		}}, "sample bases must be positive");
		expectFailure(new Runnable(){public void run(){
			mixedSchedule(root, root+"/queue_conflicting_sample", costs, "samplebases=10");
		}}, "without samplebases=");
	}
	private static void mixedSchedule(final String root, final String out, final String costs, final String extra){
		PrepareCmSearchBatches.main(new String[]{"mode=schedule", "in="+root+"/schedule_split", "out="+out,
			"modelcosts="+costs, "workers=2", "taskperbatch=t", extra, "t=1"});
	}
	private static void schedule(final String plan, final String out, final String costs, final int workers, final int extra){
		PrepareCmSearchBatches.main(new String[]{"mode=schedule", "in="+plan, "out="+out, "modelcosts="+costs,
			"samplebases=10", "workers="+workers, "extraseconds="+extra, "t=1"});
	}
	private static void metrics(final String path, final int searches, final double cpu, final double wall, final double idle){
		try(Input in=new Input(path)){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			p.set(in.reader.nextLine());
			check(p.terms()==12 && p.parseLong(4)==searches && p.parseDouble(7)==cpu && p.parseDouble(9)==wall
				&& p.parseDouble(10)==idle && p.termEquals("MODELED_NOT_MEASURED", 11), "Independent literal schedule oracle: "+path);
			check(in.reader.nextLine()==null, "One complete fleet summary");
		}
	}

	private static void plan(final String input, final String output, final String... extra){
		final String[] args=new String[extra.length+4];
		args[0]="mode=plan";
		args[1]="in="+input;
		args[2]="out="+output;
		args[3]="t=1";
		System.arraycopy(extra, 0, args, 4, extra.length);
		PrepareCmSearchBatches.main(args);
	}
	private static void invokeWrite(final String plan, final String output, final int task){
		PrepareCmSearchBatches.main(new String[]{"mode=write", "in="+plan, "out="+output, "task="+task, "t=1"});
	}
	private static String row(final long id, final String fasta, final long bases, final long records){
		return id+"\t"+fasta+"\t"+new File(fasta).length()+"\t"+bases+"\t"+records+'\n';
	}
	private static void write(final String path, final String text){
		try(Output out=new Output(path)){out.writer.print(text);}
	}
	private static void equal(final String path, final String expected){
		final ByteBuilder b=new ByteBuilder();
		try(Input in=new Input(path)){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){b.append(line).nl();}
		}
		check(expected.equals(new String(b.array, 0, b.length(), StandardCharsets.UTF_8)), "Literal fixture differs: "+path);
	}
	private static void expectFailure(final Runnable action, final String fragment){
		try{action.run();}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(fragment), "Wrong rejection: "+failure);
			return;
		}
		throw new AssertionError("Expected rejection: "+fragment);
	}
	/** Linux-only finite check of this JVM's exact output descriptor, not a heap/RSS test. */
	private static void closed(final String output){
		final File directory=new File("/proc/self/fd");
		check(directory.isDirectory(), "Lifecycle fixture requires Linux /proc/self/fd");
		final File[] descriptors=directory.listFiles();
		check(descriptors!=null && descriptors.length<4096, "Unexpected fixture descriptor inventory");
		final String target;
		try{target=new File(output).getCanonicalPath();}
		catch(IOException failure){throw new AssertionError("Cannot resolve fixture output", failure);}
		for(final File descriptor:descriptors){
			try{
				check(!Files.readSymbolicLink(descriptor.toPath()).toString().equals(target), "Failed preparation left its FASTA stream open: "+output);
			}catch(IOException vanished){
				// Another completed native I/O owner may close a descriptor after listing.
				check(!descriptor.exists(), "Could not inspect a still-existing JVM descriptor: "+descriptor);
			}
		}
	}
	private static void check(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
}
