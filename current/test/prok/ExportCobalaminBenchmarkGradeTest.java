package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import parse.LineParser1;
import prok.ExportCobalaminBenchmarkGrade;
import prok.CompareCobalaminBenchmarkGrades;
import prok.GradeNcrnaReference;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Finite post-match oracle with repeated contig/reference names across inputs,
 * plus/minus errors, masks, duplicates, zero-call genomes and legacy full-core
 * coordinates different from the source interval. Model-win fixtures use tiny
 * literal consensus FASTAs; no genome DNA or GFF is read.
 * @author Raiden
 */
public final class ExportCobalaminBenchmarkGradeTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, ExportCobalaminBenchmarkGradeTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown benchmark-fixture argument: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && new File(parser.out1).mkdir(), "Supply out=fresh_fixture_directory");
			new ExportCobalaminBenchmarkGradeTest(parser.out1).run();
			pp.outstream.println("COBALAMIN_BENCHMARK_EXPORT_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private ExportCobalaminBenchmarkGradeTest(final String directory){root=new File(directory).getAbsolutePath();}

	private void run() throws Exception{
		write("native.time.tsv", "12.5\t0.5\t2.0\n");
		write("caller.stderr", "Unrelated startup text\n"
			+"ncRNA Alignments (cobalamin_rf00174):\t        11\n"
			+"ncRNA Alignments (cobalamin_rf01482):\t        22\n"
			+"ncRNA Alignments (cobalamin_rf01689):\t        33\n");
		final String inputs=GradeNcrnaReference.INPUT_HEADER+"\n"
			+"a\t/missing/a.fa\t/q/a\t/r/a\t3\t0\t0\t0\t3\t0\n"
			+"b\t/missing/b.fa\t/q/b\t/r/b\t1\t0\t0\t0\t1\t0\n"
			+"empty\t/missing/e.fa\t/q/e\t/r/e\t0\t0\t0\t0\t0\t0\n";
		final String raw=GradeNcrnaReference.RAW_HEADER+"\n"
			+"a\tQUERY\tcobalamin\tshared\t102\t199\t+\tTP\t.\t1\tr\t.\t.\tfalse\ttrue\tcobalamin_RF00174_one\t2\t0\n"
			+"b\tQUERY\tcobalamin\tshared\t301\t395\t-\tTP\t.\t1\tr\t.\t.\tfalse\ttrue\tcobalamin_RF01482_two\t4\t-1\n"
			+"a\tQUERY\tcobalamin\tshared\t500\t550\t+\tFP\t.\t2\t.\t.\t.\tfalse\tfalse\tcobalamin_RF01689_three\t.\t.\n"
			+"a\tQUERY\tcobalamin\tshared\t700\t750\t+\tMASKED\t.\t3\t.\t.\tm\ttrue\tfalse\tcobalamin_RF00174_one\t.\t.\n"
			+"a\tREFERENCE\tcobalamin\tshared\t100\t199\t+\tTP\tr\t.\t.\t1\t.\t.\t.\tRF00174\t.\t.\n"
			+"b\tREFERENCE\tcobalamin\tshared\t300\t399\t-\tTP\tr\t.\t.\t1\t.\t.\t.\tRF00174\t.\t.\n"
			+"a\tREFERENCE\tcobalamin\tshared\t700\t750\t+\tDONT_CARE\tm\t.\t.\t.\t.\t.\t.\tRF00174\t.\t.\n";
		final String summary=GradeNcrnaReference.SUMMARY_HEADER+"\n"
			+"cobalamin\t3\t2\t1\t0\t2\t1\t1\t1\t0\t0\t1\t0\t.666667\t1\t.5\tfixture\n";
		write("canonical.inputs.tsv", inputs);
		write("canonical.raw.tsv", raw);
		write("canonical.summary.tsv", summary);
		invoke("canonical", "canonical.out", false);
		final String matches="#family\tseqid\tstrand\ttruth_start\ttruth_stop\tcall_start\tcall_stop\td5\td3\n"
			+"cobalamin\t1:a/shared\t+\t100\t199\t102\t199\t2\t0\n"
			+"RF00174\t1:a/shared\t+\t100\t199\t102\t199\t2\t0\n"
			+"cobalamin\t1:b/shared\t-\t300\t399\t301\t395\t4\t-1\n"
			+"RF01482\t1:b/shared\t-\t300\t399\t301\t395\t4\t-1\n";
		check(read("canonical.out.matches.tsv").equals(matches), "Exact independent coordinates, gene-sense signs and source-qualified repeated contigs must survive");
		check(read("canonical.out.genomes.tsv").equals("input_id\tcalls\tmatched\tmasked\tunmatched\na\t3\t1\t1\t1\nb\t1\t1\t0\t0\nempty\t0\t0\t0\t0\n"), "Negative inputs remain in the genome denominator");
		check(read("canonical.out.members.tsv").contains("RF00174\t2\t1\t1\t0\t0\t"), "Masked duplicate stays masked in the disjoint member partition");
		check(read("canonical.out.metrics.tsv").contains("\t12.5\t0.5\t2.0\t11\t22\t33\n"), "Copy actual user/system/wall seconds and all three native alignment counters without estimating costs");
		check(read("canonical.out.queries.tsv").contains("a\tshared\t500\t550\t+\tRF01689\tcobalamin_RF01689_three\tUNMATCHED_OTHER\tfalse\t.\n"), "Unmatched short-member calls and original eligibility remain available for length diagnostics");
		check(read("canonical.out.witnesses.tsv").contains("b\tr\tRF01482\tcobalamin_RF01482_two\tshared\t301\t395\t300\t399\t-\t4\t-1\n"), "Winning model and both coordinate intervals must survive for loss attribution");
		legacy(matches);
		inclusive(inputs, raw);
		winnerCounts();
		memberOutcomes();
		genericFamily(inputs, raw, summary);
		write("bad_offset.inputs.tsv", inputs);
		write("bad_offset.summary.tsv", summary);
		write("bad_offset.raw.tsv", raw.replace("two\t4\t-1", "two\t-4\t-1"));
		reject("bad_offset", "offsets");
		write("bad_model.inputs.tsv", inputs);
		write("bad_model.summary.tsv", summary);
		write("bad_model.raw.tsv", raw.replace("cobalamin_RF01689_three", "cobalamin_RF016890_three"));
		reject("bad_model", "Unknown emitted");
		write("bad_reference.inputs.tsv", inputs);
		write("bad_reference.summary.tsv", summary);
		write("bad_reference.raw.tsv", raw.replace("\t1\tr\t.\t.", "\t1\tmissing\t.\t."));
		reject("bad_reference", "declared scored reference");
		write("bad_summary.inputs.tsv", inputs);
		write("bad_summary.raw.tsv", raw);
		write("bad_summary.summary.tsv", summary.replace("cobalamin\t3\t2\t1", "cobalamin\t3\t3\t1"));
		reject("bad_summary", "disagree");
		write("missing_counter.stderr", "ncRNA Alignments (cobalamin_rf00174):\t11\n");
		boolean missingCounter=false;
		try{
			ExportCobalaminBenchmarkGrade.main(new String[]{"in="+path("canonical"), "out="+path("missing_counter"),
				"timing="+path("native.time.tsv"), "callerlog="+path("missing_counter.stderr")});
		}catch(IllegalArgumentException expected){missingCounter=expected.getMessage().contains("lacks");}
		check(missingCounter && !new File(path("missing_counter.metrics.tsv")).exists(), "A partial caller log must not become a completed cost report");
	}

	/** A second family keeps its biological union, literal models and zero-call
	 * genomes, while unrelated raw families cannot enter its normalized report. */
	private void genericFamily(final String inputs, final String raw, final String summary) throws Exception{
		final String selected=raw.replace("\tcobalamin\t", "\tRNaseP\t")
			.replace("cobalamin_RF00174_one", "rnasep_one").replace("cobalamin_RF01482_two", "rnasep_two")
			.replace("cobalamin_RF01689_three", "rnasep_three");
		write("generic.inputs.tsv", inputs);
		write("generic.raw.tsv", selected+raw.substring(raw.indexOf('\n')+1));
		write("generic.summary.tsv", summary.replace("cobalamin\t", "RNaseP\t")+summary.substring(summary.indexOf('\n')+1));
		ExportCobalaminBenchmarkGrade.main(new String[]{"in="+path("generic"), "out="+path("generic.out"), "family=RNaseP"});
		check(read("generic.out.family.tsv").equals("#NCRNA_NORMALIZED_FAMILY_V1\nfamily\tRNaseP\n"),
			"Generic exports must bind the selected family even when no queries exist");
		check(read("generic.out.genomes.tsv").equals(read("canonical.out.genomes.tsv")),
			"The selected family keeps its exact calls and all three input genomes");
		check(read("generic.out.matches.tsv").equals("#family\tseqid\tstrand\ttruth_start\ttruth_stop\tcall_start\tcall_stop\td5\td3\n"
			+"RNaseP\t1:a/shared\t+\t100\t199\t102\t199\t2\t0\n"
			+"RNaseP\t1:b/shared\t-\t300\t399\t301\t395\t4\t-1\n"),
			"Generic matches contain one union view, with actual signed offsets and no invented subtype rows");
		check(read("generic.out.members.tsv").split("\n").length==2
			&& read("generic.out.members.tsv").contains("\nRNaseP\t4\t2\t1\t0\t1\t")
			&& !read("generic.out.metrics.tsv").contains("align_rf")
			&& read("generic.out.queries.tsv").contains("\tfamily_union\t")
			&& !read("generic.out.queries.tsv").contains("emitted_member"),
			"A generic union has one cohort and no unrelated cobalamin telemetry columns");
		for(final String suffix:new String[]{"inputs", "overlap.raw", "overlap.summary"}){
			write("generic_inexact."+suffix+".tsv", read("inclusive."+suffix+".tsv").replace("cobalamin", "RNaseP"));
		}
		ExportCobalaminBenchmarkGrade.main(new String[]{"in="+path("generic_inexact"), "out="+path("generic_inexact.out"), "family=RNaseP", "inexact=t"});
		check(read("generic_inexact.out.matches.tsv").split("\n").length==3
			&& read("generic_inexact.out.matches.tsv").contains("RNaseP\t1:a/shared\t+\t100\t199\t55\t105\t-45\t-94\n"),
			"A second family's recorded inexact witness must survive without rematching or duplicate subtype views");
		CompareCobalaminBenchmarkGrades.main(new String[]{"in1="+path("generic.out"), "in2="+path("generic.out"),
			"out="+path("generic.compare"), "family=RNaseP", "radius=3"});
		check(read("generic.compare.paired_population.tsv").contains("RNaseP\t5prime\t3\t2\t1\n"),
			"The real normalized comparison must consume generic exports and preserve both matched references");
		for(final String suffix:new String[]{"inputs", "overlap.raw", "overlap.summary"}){
			write("generic_empty."+suffix+".tsv", read("axis3."+suffix+".tsv").replace("cobalamin", "RNaseP"));
		}
		ExportCobalaminBenchmarkGrade.main(new String[]{"in="+path("generic_empty"), "out="+path("generic_empty.out"), "family=RNaseP", "inexact=t"});
		check(read("generic_empty.out.family.tsv").equals(read("generic.out.family.tsv"))
			&& read("generic_empty.out.genomes.tsv").contains("empty\t0\t0\t0\t0\n"),
			"An empty family still declares its identity and zero-call input population");
		check(!new File(path("canonical.out.family.tsv")).exists(), "Default cobalamin output file set must remain unchanged");
	}

	/** The main call overlaps six bases of a100-base reference. Its recorded
	 * inclusive match must survive without inventing a reciprocal-half match. */
	private void inclusive(final String inputs, final String raw) throws Exception{
		final String overlapRaw=raw.replace("shared\t102\t199\t+\tTP", "shared\t55\t105\t+\tTP_INEXACT")
			.replace("one\t2\t0\n", "one\t-45\t-94\n")
			.replace("shared\t100\t199\t+\tTP\tr", "shared\t100\t199\t+\tTP_INEXACT\tr");
		final String summary=GradeNcrnaReference.INEXACT_HEADER+'\n'
			+"ALL\tcobalamin\t3\t1\t1\t2\t1\t2\t1\t0\t1\t1\t2\t1\t1\t0\t0\t0\tfixture\n";
		write("inclusive.inputs.tsv", inputs);
		write("inclusive.overlap.raw.tsv", overlapRaw);
		write("inclusive.overlap.summary.tsv", summary);
		invokeInexact("inclusive", "inclusive.out");
		check(read("inclusive.out.matches.tsv").equals("#family\tseqid\tstrand\ttruth_start\ttruth_stop\tcall_start\tcall_stop\td5\td3\n"
			+"cobalamin\t1:a/shared\t+\t100\t199\t55\t105\t-45\t-94\n"
			+"RF00174\t1:a/shared\t+\t100\t199\t55\t105\t-45\t-94\n"
			+"cobalamin\t1:b/shared\t-\t300\t399\t301\t395\t4\t-1\n"
			+"RF01482\t1:b/shared\t-\t300\t399\t301\t395\t4\t-1\n"),
			"Inclusive export retains literal assigned endpoints, input identity and gene-sense errors without rematching");
		check(read("inclusive.out.genomes.tsv").equals("input_id\tcalls\tmatched\tmasked\tunmatched\na\t3\t1\t1\t1\nb\t1\t1\t0\t0\nempty\t0\t0\t0\t0\n"),
			"Inexact matches remain matched while masks, unmatched calls and zero-call inputs keep their denominators");
		check(read("inclusive.out.queries.tsv").contains("a\tshared\t55\t105\t+\tRF00174\tcobalamin_RF00174_one\tTP\ttrue\tr\n"),
			"Existing normalized comparators receive the actual inclusive assignment");
		check(read("inclusive.out.metrics.tsv").contains("3\t2\t1\t4\t2\t1\t0\t1\t1\t0\t0\t0\t")
			&& read("inclusive.out.metrics.tsv").contains("\tCANONICAL_NEW_CM_UNION_INEXACT\t"),
			"Inclusive metrics reconcile assigned matches and explicitly name their matching view");

		write("inexact_without_flag.inputs.tsv", inputs);
		write("inexact_without_flag.raw.tsv", overlapRaw);
		write("inexact_without_flag.summary.tsv", read("canonical.summary.tsv"));
		reject("inexact_without_flag", "unknown recorded outcome");
		write("bad_inexact.inputs.tsv", inputs);
		write("bad_inexact.overlap.raw.tsv", overlapRaw);
		write("bad_inexact.overlap.summary.tsv", GradeNcrnaReference.INEXACT_HEADER+'\n'
			+"ALL\tcobalamin\t3\t1\t1\t2\t1\t2\t0\t0\t1\t1\t2\t1\t1\t0\t0\t0\tfixture\n");
		boolean badCount=false;
		try{invokeInexact("bad_inexact", "bad_inexact.out");}
		catch(IllegalArgumentException expected){badCount=expected.getMessage().contains("inexact match counts");}
		check(badCount && !new File(path("bad_inexact.out.metrics.tsv")).exists(), "An inconsistent recorded inexact count cannot become a normalized result");
		boolean badLegacy=false;
		try{
			ExportCobalaminBenchmarkGrade.main(new String[]{"in="+path("legacy"), "out="+path("legacy_inexact"),
				"legacy=t", "inexact=t", "panel="+path("panel.tsv")});
		}catch(IllegalArgumentException expected){badLegacy=expected.getMessage().contains("legacy grades");}
		check(badLegacy, "Explicit canonical overlap mode must not reinterpret historical legacy matching");
	}
	/** Independent six-count oracle straddles both requested strict thresholds.
	 * An inexact winner counts once; masked and unmatched calls never win. */
	private void winnerCounts() throws Exception{
		final String[] ids={"cobalamin_RF00174_zero", "cobalamin_RF00174_one", "cobalamin_RF01482_fourteen",
			"cobalamin_RF01482_fifteen", "cobalamin_RF01689_fortyfive", "cobalamin_RF01689_fortysix"};
		final int[] counts={0, 1, 14, 15, 45, 46};
		for(int member=0; member<3; member++){
			write("wins.models"+member+".fa", ">"+ids[2*member]+" n=900\nACGT\n>"+ids[2*member+1]+" n=901\nTGCA\n");
		}
		final StringBuilder raw=new StringBuilder(GradeNcrnaReference.RAW_HEADER+'\n');
		int line=0;
		for(int model=0; model<ids.length; model++){
			for(int i=0; i<counts[model]; i++){
				line++;
				final String outcome=line==1 ? "TP_INEXACT" : "TP";
				final int start=line==1 ? 101 : 100;
				raw.append("a\tQUERY\tcobalamin\tc").append(line).append('\t').append(start).append("\t199\t+\t")
					.append(outcome).append("\t.\t").append(line).append("\tr").append(line).append("\t.\t.\tfalse\ttrue\t")
					.append(ids[model]).append(" n=999\t").append(start-100).append("\t0\n");
				raw.append("a\tREFERENCE\tcobalamin\tc").append(line).append("\t100\t199\t+\t").append(outcome)
					.append("\tr").append(line).append("\t.\t.\t").append(line).append("\t.\t.\t.\tRF00174\t.\t.\n");
			}
		}
		check(line==121, "Literal six-model count population has121 matched queries");
		raw.append("a\tQUERY\tcobalamin\textra\t1\t4\t+\tFP\t.\t122\t.\t.\t.\tfalse\tfalse\t").append(ids[0]).append("\t.\t.\n")
			.append("a\tQUERY\tcobalamin\tmask\t1\t4\t+\tMASKED\t.\t123\t.\t.\tm\tfalse\tfalse\t").append(ids[0]).append("\t.\t.\n")
			.append("a\tREFERENCE\tcobalamin\tmask\t1\t4\t+\tDONT_CARE\tm\t.\t.\t.\t.\t.\t.\tRF00174\t.\t.\n");
		write("wins.inputs.tsv", GradeNcrnaReference.INPUT_HEADER+"\na\t/missing/a.fa\t/q/a\t/r/a\t123\t0\t0\t0\t123\t0\nempty\t/missing/e.fa\t/q/e\t/r/e\t0\t0\t0\t0\t0\t0\n");
		write("wins.overlap.raw.tsv", raw.toString());
		write("wins.overlap.summary.tsv", GradeNcrnaReference.INEXACT_HEADER+"\nALL\tcobalamin\t2\t120\t1\t2\t1\t121\t1\t0\t1\t1\t121\t1\t0\t120\t121\t120\tfixture\n");
		ExportCobalaminBenchmarkGrade.main(winnerArguments("wins.out"));
		final LineParser1 p=new LineParser1('\t');
		final String[] rows=read("wins.out.model_wins.tsv").split("\n");
		check(rows.length==7, "Full model roster includes the zero-win model");
		for(int i=0; i<counts.length; i++){
			p.set(rows[i+1].getBytes(StandardCharsets.UTF_8));
			check(p.terms()==5 && p.parseInt(1)==i%2 && p.termEquals(ids[i], 2) && p.parseLong(3)==counts[i],
				"Every original model ordinal and literal assigned-win count must survive");
		}
		final double[][] expected={{6,121,0,4.25,14.5,37.5,46,1,3,5}, {2,1,0,.25,.5,.75,1,1,2,2},
			{2,29,14,14.25,14.5,14.75,15,0,1,2}, {2,91,45,45.25,45.5,45.75,46,0,0,1}};
		final String[] labels={"ALL", "RF00174", "RF01482", "RF01689"};
		final String[] summaries=read("wins.out.model_win_summary.tsv").split("\n");
		check(summaries.length==5, "One ALL plus exactly three emitted-member distributions");
		for(int i=0; i<4; i++){
			p.set(summaries[i+1].getBytes(StandardCharsets.UTF_8));
			check(p.terms()==14 && p.termEquals(labels[i], 0) && p.parseInt(1)==2, "Zero-call genome stays in the declared panel");
			for(int column=2; column<12; column++){
				check(p.parseDouble(column)==expected[i][column-2], "Literal zero/threshold/quartile oracle differs at row="+i+" column="+column);
			}
			check(p.termEquals("LINEAR_INTERPOLATION_AT_(N-1)*P", 12) && p.parseString(13).startsWith("INCLUSIVE_ASSIGNED_QUERY_MATCHES"), "Quantile and matching definitions must remain explicit");
		}
		write("wins.bad.models.fa", ">cobalamin_RF00174_absent\nACGT\n>"+ids[1]+"\nTGCA\n");
		String[] args=winnerArguments("wins.bad_unknown");
		args[3]="ref="+path("wins.bad.models.fa")+','+path("wins.models1.fa")+','+path("wins.models2.fa");
		rejectWinner(args, "Every emitted model");
		args=winnerArguments("wins.bad_count");
		args[4]="expectedmodels=2,2,3";
		rejectWinner(args, "roster count");
		args=winnerArguments("wins.bad_panel");
		args[5]="expectedinputs=3";
		rejectWinner(args, "input population");
		write("wins.foreign.models.fa", ">"+ids[2]+"\nACGT\n>"+ids[3]+"\nTGCA\n");
		args=winnerArguments("wins.bad_member");
		args[3]="ref="+path("wins.foreign.models.fa")+','+path("wins.models1.fa")+','+path("wins.models2.fa");
		rejectWinner(args, "declared member");
		// Existing strict fixture has one winner in each of the first two members,
		// plus a fully zero-win third member; the ordinary export files stay identical.
		write("wins.strict0.fa", ">cobalamin_RF00174_one\nACGT\n>cobalamin_RF00174_unused\nTGCA\n");
		write("wins.strict1.fa", ">cobalamin_RF01482_two\nACGT\n");
		write("wins.strict2.fa", ">cobalamin_RF01689_three\nACGT\n");
		ExportCobalaminBenchmarkGrade.main(new String[]{"in="+path("canonical"), "out="+path("wins.strict"),
			"ref="+path("wins.strict0.fa")+','+path("wins.strict1.fa")+','+path("wins.strict2.fa"), "expectedmodels=2,1,1", "expectedinputs=3"});
		check(read("wins.strict.model_win_summary.tsv").contains("RF01689\t3\t1\t0\t0\t"), "A member with no matched wins must retain its complete model population");
		for(final String suffix:new String[]{".members.tsv", ".genomes.tsv", ".matches.tsv", ".witnesses.tsv", ".queries.tsv", ".references.tsv"}){
			check(read("wins.strict"+suffix).equals(read("canonical.out"+suffix)), "Optional winner reporting must leave existing export bytes unchanged");
		}
	}
	private String[] winnerArguments(final String output){
		return new String[]{"in="+path("wins"), "out="+path(output), "inexact=t", "ref="+path("wins.models0.fa")+','+path("wins.models1.fa")+','+path("wins.models2.fa"), "expectedmodels=2,2,2", "expectedinputs=2"};
	}
	private void rejectWinner(final String[] args, final String reason){
		boolean failed=false;
		try{ExportCobalaminBenchmarkGrade.main(args);}
		catch(IllegalArgumentException error){check(error.getMessage().contains(reason), "Wrong model-win failure: "+error);failed=true;}
		check(failed, "Invalid model-win binding must fail");
		check(!new File(args[1].substring(4)+".model_win_summary.tsv").exists(), "Invalid model roster cannot publish its distribution");
	}
	private void invokeInexact(final String prefix, final String output){
		ExportCobalaminBenchmarkGrade.main(new String[]{"in="+path(prefix), "out="+path(output), "inexact=t"});
	}

	/** Four independent partitions reset query lines and reference names, while
	 * their original input IDs stay disjoint. Expected counts are literal, and the
	 * emitted/representative axes deliberately disagree on every scored match. */
	private void memberOutcomes() throws Exception{
		final String[] ids={"a", "b", "c", "empty"};
		final int[] calls={4, 1, 1, 0};
		final String[] raw={
			"a\tQUERY\tcobalamin\tshared\t102\t199\t+\tTP\t.\t1\tr\t.\t.\tfalse\ttrue\tcobalamin_RF00174_one\t2\t0\n"
			+"a\tQUERY\tcobalamin\tshared\t700\t750\t+\tMASKED\t.\t2\t.\t.\tm\tfalse\tfalse\tcobalamin_RF00174_one\t.\t.\n"
			+"a\tQUERY\tcobalamin\tshared\t500\t550\t+\tFP\t.\t3\t.\t.\t.\tfalse\tfalse\tcobalamin_RF01689_three\t.\t.\n"
			+"a\tQUERY\tcobalamin\tshared\t102\t199\t+\tFP\t.\t4\t.\t.\t.\ttrue\ttrue\tcobalamin_RF01689_three\t.\t.\n"
			+"a\tREFERENCE\tcobalamin\tshared\t100\t199\t+\tTP\tr\t.\t.\t1\t.\t.\t.\tRF01482\t.\t.\n"
			+"a\tREFERENCE\tcobalamin\tshared\t700\t750\t+\tDONT_CARE\tm\t.\t.\t.\t.\t.\t.\tRF00174\t.\t.\n"
			+"a\tREFERENCE\tcobalamin\tshared\t900\t999\t+\tFN\tmiss\t.\t.\t.\t.\t.\t.\tRF01689\t.\t.\n",
			"b\tQUERY\tcobalamin\tshared\t55\t105\t+\tTP_INEXACT\t.\t1\tr\t.\t.\tfalse\ttrue\tcobalamin_RF01482_two\t-45\t-94\n"
			+"b\tREFERENCE\tcobalamin\tshared\t100\t199\t+\tTP_INEXACT\tr\t.\t.\t1\t.\t.\t.\tRF01689\t.\t.\n",
			"c\tQUERY\tcobalamin\tshared\t100\t199\t+\tTP\t.\t1\tr\t.\t.\tfalse\ttrue\tcobalamin_RF01689_three\t0\t0\n"
			+"c\tREFERENCE\tcobalamin\tshared\t100\t199\t+\tTP\tr\t.\t.\t1\t.\t.\t.\tRF00174\t.\t.\n", ""};
		final String[] summaries={
			"ALL\tcobalamin\t1\t1\t1\t2\t1\t1\t0\t1\t2\t1\t2\t1\t1\t0\t1\t0\tfixture\n",
			"ALL\tcobalamin\t1\t0\t1\t1\t0\t1\t1\t0\t0\t0\t1\t0\t0\t0\t0\t0\tfixture\n",
			"ALL\tcobalamin\t1\t1\t0\t0\t0\t1\t0\t0\t0\t0\t1\t0\t0\t1\t1\t1\tfixture\n",
			"ALL\tcobalamin\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tfixture\n"};
		for(int i=0; i<4; i++){
			write("axis"+i+".inputs.tsv", GradeNcrnaReference.INPUT_HEADER+'\n'+ids[i]
				+"\t/missing/"+ids[i]+".fa\t/q\t/r\t"+calls[i]+"\t0\t0\t0\t"+calls[i]+"\t0\n");
			write("axis"+i+".overlap.raw.tsv", GradeNcrnaReference.RAW_HEADER+'\n'+raw[i]);
			write("axis"+i+".overlap.summary.tsv", GradeNcrnaReference.INEXACT_HEADER+'\n'+summaries[i]);
		}
		final String partitions=path("axis0")+','+path("axis1")+','+path("axis2")+','+path("axis3");
		invokeMembers(partitions, "axis.all");
		final String suffix="\tRECORDED_UNION_ASSIGNMENTS_NOT_SUBTYPE_ACCURACY\n";
		final String expected="axis\tmember\tinputs\ttotal\ttp\ttp_inexact\tmatched\tunmatched\tmasked\tmissed\tdont_care\tqualification\n"
			+"QUERY_EMITTED_MEMBER\tRF00174\t4\t2\t1\t0\t1\t0\t1\tNA\tNA"+suffix
			+"QUERY_EMITTED_MEMBER\tRF01482\t4\t1\t0\t1\t1\t0\t0\tNA\tNA"+suffix
			+"QUERY_EMITTED_MEMBER\tRF01689\t4\t3\t1\t0\t1\t2\t0\tNA\tNA"+suffix
			+"REFERENCE_REPRESENTATIVE_MEMBER\tRF00174\t4\t2\t1\t0\t1\tNA\tNA\t0\t1"+suffix
			+"REFERENCE_REPRESENTATIVE_MEMBER\tRF01482\t4\t1\t1\t0\t1\tNA\tNA\t0\t0"+suffix
			+"REFERENCE_REPRESENTATIVE_MEMBER\tRF01689\t4\t2\t0\t1\t1\tNA\tNA\t1\t0"+suffix;
		check(read("axis.all.member_outcomes.tsv").equals(expected),
			"Literal axes must preserve cross-member matches, TP with nonzero offsets, TP_INEXACT, both masks, misses, duplicates and the empty input");
		check(!new File(path("axis.all.metrics.tsv")).exists(), "Member-summary mode must not replay ordinary detailed exports");
		invokeMembers(path("axis3")+','+path("axis2")+','+path("axis1")+','+path("axis0"), "axis.reverse");
		check(read("axis.reverse.member_outcomes.tsv").equals(expected), "Disjoint partition order must not alter global counts");
		rejectMembers(partitions+','+path("axis0"), "axis.repeated", "Repeated input ID");
		for(final String name:new String[]{"class", "line", "member"}){
			write("axis.bad_"+name+".inputs.tsv", read("axis1.inputs.tsv"));
			write("axis.bad_"+name+".overlap.summary.tsv", read("axis1.overlap.summary.tsv"));
		}
		write("axis.bad_class.overlap.raw.tsv", (GradeNcrnaReference.RAW_HEADER+'\n'+raw[1])
			.replace("TP_INEXACT\tr\t.\t.\t1", "TP\tr\t.\t.\t1"));
		rejectMembers(path("axis.bad_class"), "axis.bad_class.out", "match class and reciprocal");
		write("axis.bad_line.overlap.raw.tsv", (GradeNcrnaReference.RAW_HEADER+'\n'+raw[1])
			.replace("TP_INEXACT\tr\t.\t.\t1", "TP_INEXACT\tr\t.\t.\t2"));
		rejectMembers(path("axis.bad_line"), "axis.bad_line.out", "match class and reciprocal");
		write("axis.bad_member.overlap.raw.tsv", (GradeNcrnaReference.RAW_HEADER+'\n'+raw[1]).replace("\tRF01689\t", "\tRF016890\t"));
		rejectMembers(path("axis.bad_member"), "axis.bad_member.out", "representative");
		// Summary-only mode can coexist with existing ordinary outputs. All seven
		// default files, including metrics, must remain byte-identical.
		final String[] ordinary={".members.tsv", ".genomes.tsv", ".matches.tsv", ".witnesses.tsv", ".queries.tsv", ".references.tsv", ".metrics.tsv"};
		final String[] before=new String[ordinary.length];
		for(int i=0; i<ordinary.length; i++){before[i]=read("inclusive.out"+ordinary[i]);}
		check(!new File(path("inclusive.out.member_outcomes.tsv")).exists(), "Default export must not create an optional member summary");
		invokeMembers(path("inclusive"), "inclusive.out");
		check(read("inclusive.out.member_outcomes.tsv").contains("REFERENCE_REPRESENTATIVE_MEMBER\tRF00174\t3\t3\t1\t1\t2\tNA\tNA\t0\t1"),
			"Single-prefix summary preserves three inputs and assigns both matches to the reference member");
		for(int i=0; i<ordinary.length; i++){check(read("inclusive.out"+ordinary[i]).equals(before[i]), "Optional summary must leave existing default output bytes unchanged");}
	}
	private void invokeMembers(final String partitions, final String output){
		ExportCobalaminBenchmarkGrade.main(new String[]{"in="+partitions, "out="+path(output), "mode=members", "inexact=t"});
	}
	private void rejectMembers(final String partitions, final String output, final String reason){
		boolean failed=false;
		try{invokeMembers(partitions, output);}
		catch(IllegalArgumentException error){check(error.getMessage().contains(reason), "Wrong member-summary failure: "+error);failed=true;}
		check(failed && !new File(path(output+".member_outcomes.tsv")).exists(), "Invalid partitions or assignments cannot publish a member summary");
	}

	/** Historical unscored/duplicate behavior is retained, not normalized to new matching. */
	private void legacy(final String matches) throws Exception{
		write("panel.tsv", "locus_id\tgenome\tseqid\tsource_family\tmember_family\tsource_start1\tsource_stop1\tcore_start1\tcore_stop1\tstrand\n"
			+"r\ta\tshared\tRF00174\tRF00174\t90\t210\t100\t199\t+\n"
			+"r\tb\tshared\tRF00174\tRF00174\t290\t410\t300\t399\t-\n"
			+"m\ta\tshared\tRF00174\tUNRESOLVED\t700\t750\t0\t0\t+\n");
		write("legacy.genomes.tsv", "view\tscope\tgenome\ttruth\tcalls\tTP\tFP\tFN\tunscored\texact_duplicates\texact5\texact3\tboth_exact\n"
			+"FULL_CM\tunion\ta\t1\t4\t1\t2\t0\t1\t1\t0\t1\t0\n"
			+"FULL_CM\tunion\tb\t1\t1\t1\t0\t0\t0\t0\t0\t0\t0\n"
			+"FULL_CM\tunion\tempty\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n");
		write("legacy.calls.tsv", "view\tscope\tcall_index\tgenome\tseqid\temitted_family\tstart1\tstop1\tstrand\tmodel\toutcome\ttruth_id\td5\td3\n"
			+"FULL_CM\tunion\t0\ta\tshared\tRF00174\t102\t199\t+\tcobalamin_RF00174_one\tTP\tr\t2\t0\n"
			+"FULL_CM\tunion\t1\tb\tshared\tRF01482\t301\t395\t-\tcobalamin_RF01482_two\tTP\tr\t4\t-1\n"
			+"FULL_CM\tunion\t2\ta\tshared\tRF01689\t500\t550\t+\tcobalamin_RF01689_three\tAPPARENT_FP\tNA\tNA\tNA\n"
			+"FULL_CM\tunion\t3\ta\tshared\tRF00174\t700\t750\t+\tcobalamin_RF00174_one\tUNSCORED\tNA\tNA\tNA\n"
			+"FULL_CM\tunion\t4\ta\tshared\tRF01689\t102\t199\t+\tcobalamin_RF01689_three\tAPPARENT_FP\tNA\tNA\tNA\n");
		write("legacy.summary.tsv", "view\tscope\tgenomes\ttruth\tcalls\tTP\tFP\tFN\tunscored\texact_duplicates\tprecision\trecall\tapparent_FP_per_genome\texact5\texact3\tboth_exact\n"
			+"FULL_CM\tunion\t3\t2\t5\t2\t2\t0\t1\t1\t.5\t1\t.666667\t0\t1\t0\n");
		invoke("legacy", "legacy.out", true);
		check(read("legacy.out.matches.tsv").equals(matches), "Legacy FULL_CM coordinates, not larger source intervals, define matched ends");
		check(read("legacy.out.members.tsv").contains("RF01689\t2\t0\t0\t1\t1\t"), "Cross-member duplicate remains an unmatched duplicate in the legacy union");
		check(read("legacy.out.references.tsv").contains("a\tm\tshared\t700\t750\t+\tfalse\n"), "Historical unresolved source-interval masks must survive for overlap and length diagnostics");
		CompareCobalaminBenchmarkGrades.main(new String[]{"in1="+path("legacy.out"), "in2="+path("legacy.out"), "out="+path("legacy.mask_lengths")});
		final LineParser1 fields=new LineParser1('\t');
		int masked=0;
		for(final String line:read("legacy.mask_lengths.lengths.tsv").split("\n")){
			fields.set(line.getBytes(StandardCharsets.UTF_8));
			if(!fields.termEquals("MASKED", 4)){continue;}
			check(fields.termEquals("m", 5) && fields.termEquals("MASK", 6) && fields.parseInt(7)==51
				&& fields.parseInt(8)==51 && fields.parseInt(9)==51 && fields.parseDouble(10)==1,
				"Legacy mask geometry/length ratio must survive independently of numeric text formatting");
			masked++;
		}
		check(masked==2, "Both comparison sides must retain the same legacy mask witness");
		check(read("legacy.out.metrics.tsv").contains("3\t2\t1\t5\t2\t2\t0\t1\t1\t0\t1\t0\t"), "Java metrics retain mask-reference count separately from masked calls and conserve all five queries");
	}
	private void invoke(final String prefix, final String output, final boolean legacy){
		ExportCobalaminBenchmarkGrade.main(new String[]{"in="+path(prefix), "out="+path(output), "legacy="+legacy,
			"panel="+path("panel.tsv"), "timing="+path("native.time.tsv"), "callerlog="+path("caller.stderr")});
	}
	private void reject(final String prefix, final String reason){
		boolean failed=false;
		try{invoke(prefix, prefix+".rejected", false);}
		catch(IllegalArgumentException expected){failed=expected.getMessage().contains(reason);}
		check(failed, "Corrupt recorded grades must fail for the stated reason: "+prefix);
	}
	private void write(final String name, final String value){try(Output out=new Output(path(name))){out.writer.print(value);}}
	private String read(final String name) throws Exception{return new String(Files.readAllBytes(new File(path(name)).toPath()), StandardCharsets.UTF_8);}
	private String path(final String name){return new File(root, name).getPath();}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final String root;
}
