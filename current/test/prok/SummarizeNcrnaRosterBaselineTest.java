package test.prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CaptureNcrnaWindows;
import prok.PrepareNcrnaRosterBaseline;
import prok.PrepareNcrnaBaselineShards;
import prok.SummarizeNcrnaRosterBaseline;
import prok.SummarizeNcrnaBaselineMisses;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Literal saved-event fixture: source-frame identity, neighbor isolation,
 * highest accepted identity, ambiguous proposals, final one-to-one detection,
 * missing/zero-window sources and GC population moments. No DNA or caller.
 * Adapted from template.A_SampleBasic2's output-only parser lifecycle.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaRosterBaselineTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SummarizeNcrnaRosterBaselineTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown fixture argument: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append,
				"Use out=fresh_fixture_directory without input or overwrite");
			new SummarizeNcrnaRosterBaselineTest(new File(parser.out1).getAbsoluteFile()).run();
			pp.outstream.println("NCRNA_ROSTER_BASELINE_SUMMARY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private SummarizeNcrnaRosterBaselineTest(final File root_){
		root=root_;
		check(root.mkdir(), "The saved-baseline fixture requires a fresh directory");
	}

	private void run(){
		final ByteBuilder frames=new ByteBuilder("input_id\tassembly_group\tcoordinate_frame\tsupplied_fasta\tfile_bytes\tnucleotide_bases\tfasta_records\toriginal_genes\tinput_state\n");
		frame(frames, 0, 7, "older", 2, 4, true);
		frame(frames, 1, 7, "live", 2, 3, true);
		frame(frames, 2, 8, "zero", 1, 0, true);
		frame(frames, 3, 9, "missing", 1, 1, false);
		write("population.frames.tsv", frames);
		final ByteBuilder originals=new ByteBuilder(PrepareNcrnaRosterBaseline.ORIGINAL_HEADER)
			.append("\tbaseline_input_id\tbaseline_input_state\n");
		gene(originals, 10, 0, 7, "older", "same%name", 101, 200, 10, false);
		gene(originals, 11, 0, 7, "older", "same%name", 301, 400, 11, false);
		gene(originals, 12, 0, 7, "older", "same%name", 501, 550, 12, false);
		gene(originals, 13, 0, 7, "older", "same%name", 551, 600, 13, false);
		gene(originals, 30, 1, 7, "live", "same%name", 101, 200, 30, false);
		gene(originals, 31, 1, 7, "live", "no_seed", 51, 100, 31, false);
		gene(originals, 32, 1, 7, "live", "no_seed", 201, 250, -1, true);
		gene(originals, 20, 3, 9, "missing", "same%name", 101, 200, -1, true);
		write("population.original_calls.tsv.gz", originals);
		final ByteBuilder assignments=new ByteBuilder("#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus\n");
		assignment(assignments, 10, B);
		assignment(assignments, 11, A);
		assignment(assignments, 12, B);
		assignment(assignments, 13, B);
		assignment(assignments, 30, B);
		assignment(assignments, 31, C);
		write("assignments.tsv", assignments);
		capture("native", 0);
		summarize("native", "summary");
		checkMetric("original_genes", 8);
		checkMetric("matchable_original_genes", 7);
		checkMetric("unresolved_original_genes", 1);
		checkMetric("completed_inputs", 3);
		checkMetric("unresolved_input_genes", 1);
		checkMetric("available_input_genes_without_scheduled_overlap", 2);
		checkMetric("scheduled_window_events", 5);
		checkMetric("actual_alignment_events", 20);
		checkMetric("gene_window_associations_nonadditive_cost", 7);
		checkMetric("captured_windows_containing_original_gene", 5);
		checkMetric("unique_best_members", 3);
		checkMetric("construction_comparable_genes", 3);
		checkMetric("construction_model_changes", 2);
		checkMetric("genes_with_ambiguous_best_acceptance", 2);
		checkMetric("final_one_to_one_detected_genes", 3);
		checkMetric("final_missed_genes", 5);
		checkMetric("exact_duplicate_final_calls", 1);
		checkGenes();
		checkAssociations();
		checkGc();
		checkShards();
		capture("fractional_gc", 0, true);
		summarize("fractional_gc", "fractional_summary");
		checkFractionalGc();
		capture("bad", 1);
		boolean rejected=false;
		try{summarize("bad", "rejected");}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("input identity");}
		check(rejected, "The same contig name in another source must never authorize crossing the exact input path");
		capture("duplicate_name", 2);
		rejected=false;
		try{summarize("duplicate_name", "rejected_duplicate");}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("contig token");}
		check(rejected, "A second ordinal with the same first-token name in one source must not inherit the first contig's gene credit");
		write("unresolved_population.frames.tsv", frames);
		writeUnresolvedOriginals();
		summarize("native", "unresolved_summary", "unresolved_population");
		checkMetric("unresolved_summary", "original_genes", 8);
		checkMetric("unresolved_summary", "unresolved_original_genes", 5);
		checkMetric("unresolved_summary", "final_one_to_one_detected_genes", 1);
		checkMetric("unresolved_summary", "final_missed_genes", 7);
		checkMetric("unresolved_summary", "unique_best_members", 1);
		checkMissClasses();
	}

	/** Add one scheduled-but-uncaptured rejection and one independently marked
	 * incomplete gene. This distinguishes actual capture from scheduling, preserves
	 * a final match without unique membership, and exercises all four classes. */
	private void checkMissClasses(){
		final ByteBuilder bytes=new ByteBuilder();
		try(Input in=new Input(path("population.original_calls.tsv.gz")); Output out=new Output(path("miss_population.original_calls.tsv.gz"))){
			out.writer.print(bytes.clear().append(in.reader.nextLine()).nl());
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(p.parseInt(0)==10){
					p.setBounds(21);
					bytes.clear().append(line, 0, p.a()).append("true").append(line, p.b(), line.length-p.b());
				}else{bytes.clear().append(line);}
				out.writer.print(bytes.nl());
			}
		}
		try(Input in=new Input(path("summary.genes.tsv.gz")); Output out=new Output(path("miss_grade.genes.tsv.gz"))){
			out.writer.print(bytes.clear().append(in.reader.nextLine()).nl());
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				bytes.clear();
				final boolean scheduledOnly=p.parseInt(0)==31;
				for(int column=0; column<26; column++){
					if(column>0){bytes.tab();}
					if(scheduledOnly && column==9){bytes.append(1);}
					else if(scheduledOnly && column==11){bytes.append("PARTIAL_WINDOW_NO_ACCEPTANCE");}
					else{
						p.setBounds(column);
						bytes.append(line, p.a(), p.b()-p.a());
					}
				}
				out.writer.print(bytes.nl());
			}
		}
		try(Input in=new Input(path("summary.gene_windows.tsv.gz")); Output out=new Output(path("miss_grade.gene_windows.tsv.gz"))){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){out.writer.print(bytes.clear().append(line).nl());}
			out.writer.println("31\tframe_1\ti1s1\tNA\ttrue\tfalse\tfalse\t0\ttrue");
		}
		try(Input in=new Input(path("summary.summary.tsv")); Output out=new Output(path("miss_grade.summary.tsv"))){
			out.writer.print(bytes.clear().append(in.reader.nextLine()).nl());
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(p.termEquals("gene_window_associations_nonadditive_cost", 0)){bytes.clear().append("gene_window_associations_nonadditive_cost\t8");}
				else if(p.termEquals("scheduled_window_events", 0)){bytes.clear().append("scheduled_window_events\t6");}
				else{bytes.clear().append(line);}
				out.writer.print(bytes.nl());
			}
		}
		try(Input in=new Input(path("summary.calls.tsv.gz")); Output out=new Output(path("miss_grade.calls.tsv.gz"))){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){out.writer.print(bytes.clear().append(line).nl());}
		}
		write("miss_captures.tsv", new ByteBuilder("capture_id\tcapture_prefix\treference_manifest\n\t")
			.append(path("native")).tab().append(path("unused.reference.tsv")).nl());
		write("miss_manifest.tsv", new ByteBuilder("library\tgrade_prefix\tcapture_manifest\nfixture\t")
			.append(path("miss_grade")).tab().append(path("miss_captures.tsv")).nl());
		SummarizeNcrnaBaselineMisses.main(new String[]{"in="+path("miss_manifest.tsv"), "population="+path("miss_population"),
			"out="+path("miss_report"), "expectedgenes=8", "expectedtruncated=2", "expectedincomplete=1"});
		int rows=0;
		try(Input in=new Input(path("miss_report.classes.tsv"))){
			in.reader.nextLine();
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(p.termEquals("ALL_ORIGINALS", 1)){
					check(p.parseInt(2)==8 && p.parseInt(3)==3 && p.parseInt(4)==1 && p.parseInt(5)==1
						&& p.parseInt(6)==3 && p.parseInt(7)==1 && p.parseInt(9)==1,
						"Every original needs exactly one class; scheduled-only is not captured, and final matching takes precedence over ambiguous membership");
				}else if(p.termEquals("NATIVE_TRUNCATED", 1)){
					check(p.parseInt(2)==2 && p.parseInt(3)==2, "Native truncations keep their own complete denominator");
				}else{
					check(p.termEquals("INCOMPLETE_COVERAGE", 1) && p.parseInt(2)==1 && p.parseInt(6)==1,
						"Incomplete coverage is an independent subset, not inferred from truncation or construction");
				}
				rows++;
			}
		}
		check(rows==3, "The miss report must publish the complete and both requested subset views");
		checkGeometryMiss();
	}

	/** A high-identity accepted/final span over twice the reference length is a
	 * geometric miss, not evidence that every model failed its identity cutoff. */
	private void checkGeometryMiss(){
		final ByteBuilder original=new ByteBuilder(PrepareNcrnaRosterBaseline.ORIGINAL_HEADER)
			.append("\tbaseline_input_id\tbaseline_input_state\n");
		gene(original, 10, 0, 7, "older", "same%name", 101, 200, 10, false);
		write("geometry_population.original_calls.tsv.gz", original);
		write("geometry_grade.genes.tsv.gz", new ByteBuilder(SummarizeNcrnaRosterBaseline.GENES_HEADER).nl()
			.append("10\tframe_0\tsame%25name\t101\t200\t+\t0\t10\t").append(A)
			.append("\t1\t1\tWINDOW_WITHOUT_ACCEPTED_MODEL\t-1\t-1\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tfalse\tNA\t2\tNA\tCONSTRUCTION_REPRESENTATIVE\n"));
		write("geometry_grade.gene_windows.tsv.gz", new ByteBuilder(SummarizeNcrnaRosterBaseline.ASSOCIATION_HEADER).nl()
			.append("10\tframe_0\ti0s0\ti0e0\ttrue\ttrue\tfalse\t2\ttrue\n"));
		write("geometry_grade.calls.tsv.gz", new ByteBuilder(SummarizeNcrnaRosterBaseline.CALLS_HEADER).nl()
			.append("0\tframe_0\tsame%25name\t1\t1000\t+\t").append(A).append("\t0.9\tNA\tfalse\tUNMATCHED_NOT_ADJUDICATED_FP\n"));
		write("geometry_grade.summary.tsv", new ByteBuilder("metric\tcount\noriginal_genes\t1\ngene_window_associations_nonadditive_cost\t1\nfinal_one_to_one_detected_genes\t0\nunique_best_members\t0\ncaptured_window_events\t1\nfinal_call_rows\t1\n"));
		final ByteBuilder events=new ByteBuilder(CaptureNcrnaWindows.RESULT_HEADER).nl();
		event(events, 0, "older", 0, 0, 999, 0, 999, A, .9f, 2, 50);
		write("geometry_capture.results.tsv", events);
		write("geometry_captures.tsv", new ByteBuilder("capture_id\tcapture_prefix\treference_manifest\n\t")
			.append(path("geometry_capture")).tab().append(path("unused.reference.tsv")).nl());
		write("geometry_manifest.tsv", new ByteBuilder("library\tgrade_prefix\tcapture_manifest\nfixture\t")
			.append(path("geometry_grade")).tab().append(path("geometry_captures.tsv")).nl());
		SummarizeNcrnaBaselineMisses.main(new String[]{"in="+path("geometry_manifest.tsv"), "population="+path("geometry_population"),
			"out="+path("geometry_report"), "expectedgenes=1", "expectedtruncated=0", "expectedincomplete=0"});
		try(Input in=new Input(path("geometry_report.classes.tsv"))){
			in.reader.nextLine();
			p.set(in.reader.nextLine());
			check(p.termEquals("ALL_ORIGINALS", 1) && p.parseInt(4)==1 && p.parseInt(11)==1 && p.parseInt(14)==1,
				"Accepted and final overlapping calls failing reciprocal geometry must remain distinct from absent acceptance");
		}
	}

	/** Two captures deliberately reuse local input0/event0 for different global
	 * frames. A third capture has only a zero-gene input; a fourth is unavailable.
	 * Model B's GC values lie on opposite shards, so a partial-SD average fails. */
	private void checkShards(){
		PrepareNcrnaBaselineShards.main(new String[]{"in="+path("population.frames.tsv"),
			"out="+path("planned"), "shards=2", "expectedframes=4", "expectedgenes=8"});
		checkPlan();
		splitCapture(".inputs.tsv", false);
		splitCapture(".results.tsv", true);
		splitCapture(".final_calls.tsv", false);
		write("first.reference.tsv", new ByteBuilder("input_id\tfasta\tgff\nframe_0\t")
			.append(path("older.fa")).tab().append(path("unused0.gff")).nl());
		write("second.reference.tsv", new ByteBuilder("input_id\tfasta\tgff\nframe_1\t")
			.append(path("live.fa")).tab().append(path("unused1.gff")).nl());
		write("zero.reference.tsv", new ByteBuilder("input_id\tfasta\tgff\nframe_2\t")
			.append(path("zero.fa")).tab().append(path("unused2.gff")).nl());
		final ByteBuilder manifest=new ByteBuilder("capture_id\tcapture_prefix\treference_manifest\n");
		manifest.append("first\t").append(path("first")).tab().append(path("first.reference.tsv")).nl();
		write("missing.captures.tsv", new ByteBuilder(manifest));
		manifest.append("second\t").append(path("second")).tab().append(path("second.reference.tsv")).nl();
		write("zero_missing.captures.tsv", new ByteBuilder(manifest));
		manifest.append("zero\t").append(path("zero")).tab().append(path("zero.reference.tsv")).nl();
		write("good.captures.tsv", manifest);
		checkSample();
		summarizeShards("good", "sharded");
		checkMetric("sharded", "original_genes", 8);
		checkMetric("sharded", "completed_inputs", 3);
		checkMetric("sharded", "actual_alignment_events", 20);
		checkMetric("sharded", "final_one_to_one_detected_genes", 3);
		checkMetric("sharded", "final_missed_genes", 5);
		equalOutput("summary.summary.tsv", "sharded.summary.tsv");
		equalOutput("summary.model_gc.tsv", "sharded.model_gc.tsv");
		equalOutput("summary.calls.tsv.gz", "sharded.calls.tsv.gz");
		int rows=0;
		boolean first=false, second=false;
		try(Input in=new Input(path("sharded.genes.tsv.gz"))){
			in.reader.nextLine();
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(p.parseInt(0)==10){first=p.termEquals("first:i0e2", 15);}
				if(p.parseInt(0)==30){second=p.termEquals("second:i0e0", 15);}
				rows++;
			}
		}
		check(rows==8 && first && second, "Different source frames must retain distinct capture namespaces despite repeated local slots/events");
		rejectShards("missing", "Missing inputs, including the zero-call frame, must fail before reporting recall");
		manifest.append("duplicate\t").append(path("first")).tab().append(path("first.reference.tsv")).nl();
		write("duplicate.captures.tsv", manifest);
		rejectShards("duplicate", "The same source frame cannot complete twice under different capture names");
		write("wrong.captures.tsv", new ByteBuilder("capture_id\tcapture_prefix\treference_manifest\nfirst\t")
			.append(path("first")).tab().append(path("second.reference.tsv")).nl());
		rejectShards("wrong", "A valid global input cannot be substituted into another declared shard");
		rejectShards("zero_missing", "Zero-call input bindings must be checked even though they contribute no event rows");
	}

	private void checkPlan(){
		int rows=0, inputs=0, genes=0;
		long bases=0;
		try(Input in=new Input(path("planned/shards.tsv"))){
			in.reader.nextLine();
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				check(p.terms()==7, "Each shard needs counts and exact reference/config bindings");
				rows++;
				inputs+=p.parseInt(2);
				genes+=p.parseInt(3);
				bases+=p.parseLong(4);
			}
		}
		check(rows==2 && inputs==3 && genes==8 && bases==5000,
			"The planner must conserve available inputs, saved bases and all genes including unavailable sources");
		try(Input in=new Input(path("planned/ownership.tsv"))){
			in.reader.nextLine();
			int count=0;
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(p.termEquals("frame_0", 1) || p.termEquals("frame_1", 1)){
					check(p.termEquals("shard_0", 0), "Both recorded frames of accepted assembly group7 must stay together");
				}
				count++;
			}
			check(count==4, "The zero-gene and unavailable frames must both retain ownership rows");
		}
	}

	/** Converts the literal single-capture fixture into two independent native
	 * namespaces plus one zero-gene capture; global input IDs and paths stay intact. */
	private void splitCapture(final String suffix, final boolean events){
		try(Input in=new Input(path("native"+suffix)); Output first=new Output(path("first"+suffix));
			Output second=new Output(path("second"+suffix)); Output zero=new Output(path("zero"+suffix))){
			final ByteBuilder bytes=new ByteBuilder();
			final byte[] header=in.reader.nextLine();
			first.writer.print(bytes.append(header).nl());
			second.writer.print(bytes.clear().append(header).nl());
			zero.writer.print(bytes.clear().append(header).nl());
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				final int old=p.parseInt(0), local=0;
				p.setBounds(0);
				final int end=p.b();
				bytes.clear().append(local);
				if(events){
					p.setBounds(12);
					final int startEvent=p.a(), endEvent=p.b();
					bytes.append(line, end, startEvent-end).append('i').append(local).append('e')
						.append(Parse.parseLong(line, startEvent+3, endEvent));
					bytes.append(line, endEvent, line.length-endEvent);
				}else{bytes.append(line, end, line.length-end);}
				(old==0 ? first : old==1 ? second : zero).writer.print(bytes.nl());
			}
		}
	}

	/** An explicit frame sample retains every selected original, including misses. */
	private void checkSample(){
		SummarizeNcrnaRosterBaseline.main(new String[]{"captures="+path("missing.captures.tsv"),
			"samplemanifest="+path("first.reference.tsv"), "population="+path("population"),
			"assignments="+path("assignments.tsv"), "out="+path("sample"),
			"expectedgenes=4", "expectedframes=1", "expectedmembers=6"});
		checkMetric("sample", "original_genes", 4);
		checkMetric("sample", "completed_inputs", 1);
		checkMetric("sample", "actual_alignment_events", 14);
		checkMetric("sample", "final_one_to_one_detected_genes", 2);
		checkMetric("sample", "final_missed_genes", 2);
		boolean labelled=false;
		try(Input in=new Input(path("sample.summary.tsv"))){
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				if(p.termEquals("scope", 0)){labelled=p.termEquals("EXPLICIT_INPUT_SAMPLE_NOT_FULL_POPULATION", 1);}
			}
		}
		check(labelled, "A one-frame result must never carry an implicit full-population claim");
		write("unknown_sample.tsv", new ByteBuilder("input_id\tfasta\tgff\nunknown\t")
			.append(path("older.fa")).tab().append(path("unused.gff")).nl());
		boolean rejected=false;
		try{
			SummarizeNcrnaRosterBaseline.main(new String[]{"captures="+path("missing.captures.tsv"),
				"samplemanifest="+path("unknown_sample.tsv"), "population="+path("population"),
				"assignments="+path("assignments.tsv"), "out="+path("unknown_sample"),
				"expectedgenes=4", "expectedframes=1", "expectedmembers=6"});
		}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "An unknown sample frame must fail rather than yield a reduced denominator");
	}

	private void summarizeShards(final String manifest, final String output){
		SummarizeNcrnaRosterBaseline.main(new String[]{"captures="+path(manifest+".captures.tsv"),
			"population="+path("population"), "assignments="+path("assignments.tsv"),
			"out="+path(output), "expectedgenes=8", "expectedframes=4", "expectedmembers=6"});
	}
	private void rejectShards(final String manifest, final String why){
		boolean rejected=false;
		try{summarizeShards(manifest, "rejected_"+manifest);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(path("rejected_"+manifest+".summary.tsv")).exists(), why);
	}
	private void equalOutput(final String expected, final String actual){
		try(Input a=new Input(path(expected)); Input b=new Input(path(actual))){
			for(byte[] line=a.reader.nextLine(); line!=null; line=a.reader.nextLine()){
				check(Arrays.equals(line, b.reader.nextLine()), "Global saved grading must preserve the literal fixture result: "+actual);
			}
			check(b.reader.nextLine()==null, "Sharding cannot add result rows: "+actual);
		}
	}

	private void capture(final String prefix, final int failure){
		capture(prefix, failure, false);
	}

	private void capture(final String prefix, final int failure, final boolean fractionalGc){
		final ByteBuilder inputs=new ByteBuilder("input_slot\tinput_id\tfasta\tcontigs\tbases\twindows\tscheduled_windows\talignments\taccepted_proposals\tfinal_target_call_rows\n");
		inputs.append("0\tframe_0\t").append(path("older.fa")).append("\t2\t2000\t4\t4\t14\t4\t3\n")
			.append("1\tframe_1\t").append(path("live.fa")).append("\t2\t2000\t1\t1\t6\t1\t1\n")
			.append("2\tframe_2\t").append(path("zero.fa")).append("\t1\t1000\t0\t0\t0\t0\t0\n");
		write(prefix+".inputs.tsv", inputs);
		final ByteBuilder events=new ByteBuilder(CaptureNcrnaWindows.RESULT_HEADER).nl();
		event(events, 0, "older", 0, 80, 220, 100, 199, A, .8f, 2, 20);
		// This window contains genes10and11; its accepted own span qualifies ONLY11.
		event(events, 0, "older", 1, 80, 420, 300, 399, B, .99f, 3, 40);
		event(events, 0, "older", 2, 80, 220, 100, 199, A, .9f, 4, fractionalGc ? 1 : 60, fractionalGc ? 3 : 100);
		// One100base proposal meets reciprocal1/2 against two adjacent50base genes.
		event(events, 0, "older", 3, 480, 620, 500, 599, B, 1f, 5, 50);
		event(events, 1, failure==1 ? "older" : "live", 0, 80, 220, 100, 199, B, .85f, 6, 30);
		write(prefix+".results.tsv", events);
		final ByteBuilder calls=new ByteBuilder(CaptureNcrnaWindows.FINAL_HEADER).nl();
		call(calls, 0, "older", 100, 199, A, .9f, failure==2 ? 1 : 0);
		call(calls, 0, "older", 100, 199, A, .9f, 0);
		call(calls, 0, "older", 500, 599, B, 1f, 0);
		call(calls, 1, "live", 100, 199, B, .85f, 0);
		write(prefix+".final_calls.tsv", calls);
	}

	private void checkGenes(){
		int rows=0;
		boolean best=false, neighbor=false, missing=false, noSeed=false;
		try(Input in=new Input(path("summary.genes.tsv.gz"))){
			in.reader.nextLine();
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				check(p.terms()==26, "Every gene needs separate proposal/final fields, event work and its original extraction state");
				final int id=p.parseInt(0);
				if(id==10){best=p.termEquals(A, 14) && p.termEquals("i0e2", 15) && p.parseInt(18)==60
					&& p.termEquals("true", 21) && p.parseLong(23)==9 && p.parseLong(24)==4;}
				if(id==11){neighbor=p.termEquals(B, 14) && p.termEquals("false", 21);}
				if(id==20){missing=p.termEquals("UNRESOLVED_INPUT", 11) && p.parseInt(6)==1 && p.termEquals("false", 21);}
				if(id==31){noSeed=p.termEquals("NO_SCHEDULED_WINDOW", 11) && p.termEquals("false", 21);}
				if(id==12 || id==13){check(p.termEquals("AMBIGUOUS_BEST_ACCEPTANCE", 11) && p.termEquals("NA", 14), "One ambiguous acceptance cannot assign two GC members");}
				rows++;
			}
		}
		check(rows==8 && best && neighbor && missing && noSeed, "Actual best-window GC, neighbor isolation, path loss and all-original misses must survive the join");
	}

	private void checkAssociations(){
		int rows=0;
		long associatedWork=0;
		boolean neighbor=false;
		try(Input in=new Input(path("summary.gene_windows.tsv.gz"))){
			in.reader.nextLine();
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				check(p.terms()==9, "Each gene/window association needs containment, own-span qualification and reference-resolution state");
				associatedWork+=p.parseLong(7);
				if(p.parseInt(0)==10 && p.termEquals("i0e1", 3)){
					neighbor=p.termEquals("true", 5) && p.termEquals("false", 6);
				}
				rows++;
			}
		}
		check(rows==7 && associatedWork==28 && neighbor,
			"Gene associations can share20actual alignments: their28associated alignments must never become global cost");
	}

	private void checkGc(){
		boolean a=false, b=false, c=false;
		try(Input in=new Input(path("summary.model_gc.tsv"))){
			in.reader.nextLine();
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(p.termEquals(A, 0)){a=p.parseLong(1)==1 && Math.abs(p.parseDouble(3)-.6)<1e-9 && p.termEquals("ZERO_SD", 7);}
				if(p.termEquals(B, 0)){b=p.parseLong(1)==2 && Math.abs(p.parseDouble(3)-.35)<1e-9 && Math.abs(p.parseDouble(4)-.05)<1e-9;}
				if(p.termEquals(C, 0)){c=p.parseLong(1)==0 && p.termEquals("NO_KNOWN_GC_MEMBERS", 7);}
			}
		}
		check(a && b && c, "All-data profile moments use the selected actual windows, population SD and explicit missing/zero-SD states");
	}

	/** The original one-third window must remain inside its saved zero-width range. */
	private void checkFractionalGc(){
		boolean found=false;
		try(Input in=new Input(path("fractional_summary.model_gc.tsv"))){
			in.reader.nextLine();
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(!p.termEquals(A, 0)){continue;}
				final double gc=1.0/3.0;
				final double mean=Double.parseDouble(p.parseString(3)), sd=Double.parseDouble(p.parseString(4));
				final double minimum=Double.parseDouble(p.parseString(5)), maximum=Double.parseDouble(p.parseString(6));
				check(p.parseLong(1)==1 && p.termEquals("ZERO_SD", 7)
					&& sd==0 && gc==mean && minimum==maximum && gc>=minimum && gc<=maximum,
					"GC profile text must round-trip its original mean/min/max; decimal rounding would reject this original window under an exact observed-range gate");
				found=true;
			}
		}
		check(found, "The non-terminating GC fraction must retain its actual best-model profile");
	}

	private void checkMetric(final String metric, final long expected){
		checkMetric("summary", metric, expected);
	}
	private void checkMetric(final String output, final String metric, final long expected){
		boolean found=false;
		try(Input in=new Input(path(output+".summary.tsv"))){
			in.reader.nextLine();
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				if(p.termEquals(metric, 0)){check(p.parseLong(1)==expected, "Literal baseline metric mismatch: "+metric);found=true;}
			}
		}
		check(found, "A baseline summary cannot omit an expected metric: "+metric);
	}

	private void summarize(final String capture, final String output){
		summarize(capture, output, "population");
	}
	private void summarize(final String capture, final String output, final String population){
		SummarizeNcrnaRosterBaseline.main(new String[]{"in="+path(capture), "population="+path(population),
			"assignments="+path("assignments.tsv"), "out="+path(output), "expectedgenes=8", "expectedframes=4", "expectedmembers=6"});
	}
	private void frame(final ByteBuilder out, final int slot, final int group, final String name, final int contigs, final int genes, final boolean available){
		out.append("frame_").append(slot).tab().append(group).tab().append(path(name+".fa")).tab().append(path(name+".fa"))
			.append("\t0\t").append(1000*contigs).tab().append(contigs).tab().append(genes).tab()
			.append(available ? "AVAILABLE_RECORDED_FRAME" : "UNRESOLVED_SOURCE_UNAVAILABLE").nl();
	}
	private void gene(final ByteBuilder out, final int id, final int slot, final int group, final String file,
			final String contig, final int start, final int stop, final int rep, final boolean truncated){
		out.append(id).tab().append(group).tab().append(path(file+".fa")).tab().append(contig).tab().append(start).tab().append(stop)
			.append("\t+\ttrue\tfalse\t1\t0\t0\tQUALIFIED\t").append(id).tab().append(rep)
			.tab().append(file.equals("missing") ? "SOURCE_UNAVAILABLE" : truncated ? "EXCLUDED_NATIVE_PARTIAL_OR_TRUNCATED" : "CONSTRUCTION_REPRESENTATIVE")
			.append("\t1\t").append(truncated ? 1 : 0).append("\t1\t100\t100\tfalse\tFIXTURE\t0\t0\t0\t0\tFIXTURE\t50.0\tframe_")
			.append(slot).tab().append(file.equals("missing") ? "UNRESOLVED_SOURCE_UNAVAILABLE" : "AVAILABLE_RECORDED_FRAME").nl();
	}
	private void assignment(final ByteBuilder out, final int id, final String model){
		out.append("RF00174_X_").append(id).append(" assembly_group=7\t-1\tunknown\tRF00174_0\t").append(model).append("\tkept\n");
	}
	private void event(final ByteBuilder out, final int slot, final String file, final int event, final int start, final int stop,
			final int acceptedStart, final int acceptedStop, final String model, final float identity, final int work, final int gc){
		event(out, slot, file, event, start, stop, acceptedStart, acceptedStop, model, identity, work, gc, 100);
	}

	private void event(final ByteBuilder out, final int slot, final String file, final int event, final int start, final int stop,
			final int acceptedStart, final int acceptedStop, final String model, final float identity, final int work, final int gc, final int acgt){
		out.append(slot).append("\tframe_").append(slot).tab().append(path(file+".fa"))
			.append("\t0\tsame%25name description\t0\t1000\t0\t").append(start).tab().append(stop).tab().append(start).tab().append(stop)
			.append("\ti").append(slot).append('e').append(event).append("\tACCEPTED_PROPOSAL\t").append(work).append("\t0\t0\t")
			.append(work).tab().append(model).tab().append(identity, 6).append("\t1\t1\t").append(acceptedStart).tab().append(acceptedStop)
			.tab().append(gc).tab().append(acgt).nl();
	}
	private void call(final ByteBuilder out, final int slot, final String file, final int start, final int stop, final String model, final float identity, final int ordinal){
		out.append(slot).append("\tframe_").append(slot).tab().append(path(file+".fa"))
			.tab().append(ordinal).append("\tsame%25name description\t").append(start).tab().append(stop).append("\t0\t").append(model).tab().append(identity, 6).append("\t1\t1\n");
	}
	/** Prior extractor ambiguity must survive even if only one duplicate contig
	 * produced saved windows; missing observations are not proof of uniqueness. */
	private void writeUnresolvedOriginals(){
		try(Input in=new Input(path("population.original_calls.tsv.gz"));
			Output out=new Output(path("unresolved_population.original_calls.tsv.gz"))){
			final ByteBuilder bytes=new ByteBuilder();
			out.writer.print(bytes.append(in.reader.nextLine()).nl());
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				bytes.clear();
				if(p.termEquals("frame_0", 29)){
					p.setBounds(14);
					final int start=p.a();
					p.setBounds(15);
					final int stop=p.b();
					bytes.append(line, 0, start).append("-1\tAMBIGUOUS_CONTIG_ID").append(line, stop, line.length-stop);
				}else{bytes.append(line);}
				out.writer.print(bytes.nl());
			}
		}
	}
	private String path(final String name){return new File(root, name).getPath();}
	private void write(final String name, final ByteBuilder value){try(Output out=new Output(path(name))){out.writer.print(value);}}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final File root;
	private final LineParser1 p=new LineParser1('\t');
	private static final String A="cobalamin_RF00174_A", B="cobalamin_RF00174_B", C="cobalamin_RF00174_C";
}
