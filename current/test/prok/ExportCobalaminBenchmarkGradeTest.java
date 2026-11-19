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
 * coordinates different from the source interval. No DNA or GFF exists here.
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
