package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.nio.charset.StandardCharsets;

import prok.GradeNcrnaReference;
import prok.NcrnaDdlCallerPanel;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.SummarizeNcrnaBenchmark;

/** Literal native-grade fixture: exact membership is not source-locus presence.
 * Repeated parent IDs across inputs, an unrankable UNKNOWN parent, a masked
 * multi-supported parent, duplicate/false-positive queries and an empty input
 * retain their original evidence. No matching, sequences or alignment occurs.
 * @author Ganyu
 */
public final class NcrnaExactMembershipTest{

	static void run(final File directory){
		check(directory.mkdir(), "Exact-membership fixture directory must be fresh");
		final Fixture f=new Fixture(directory);
		f.write();
		SummarizeNcrnaBenchmark.main(f.args("plain", null));
		SummarizeNcrnaBenchmark.main(f.args("exact", f.path("strata.tsv")));
		for(final String suffix:new String[]{"cells", "timing", "grades", "comparisons", "changes", "subtype_comparisons", "subtype_changes"}){
			check(rows(f.path("plain."+suffix+".tsv")).equals(rows(f.path("exact."+suffix+".tsv"))),
				"Optional exact partition must not change existing native summary/comparison output: "+suffix);
		}
		final ArrayList<String> summary=rows(f.path("exact.exact_summary.tsv"));
		check(summary.size()==13, "One declared cell retains three states plus total in each of three views");
		expect(summary, "public\t1\t1\tMODEL_SUPPORT_RECOVERY\tSEEN\t1\t0\t0\t1\tNA\t1\t"+QUAL);
		expect(summary, "public\t1\t1\tMODEL_SUPPORT_RECOVERY\tNOT_SEEN\t0\t1\t0\t1\tNA\t0\t"+QUAL);
		expect(summary, "public\t1\t1\tMODEL_SUPPORT_RECOVERY\tUNKNOWN\t1\t0\t0\t1\tNA\t1\t"+QUAL);
		expect(summary, "public\t1\t1\tEXCLUSIVE_SUBTYPE_ACCURACY\tSEEN\t0\t0\t1\t0\tNA\tNA\t"+QUAL);
		expect(summary, "public\t1\t1\tEXCLUSIVE_SUBTYPE_ACCURACY\tALL_LABELLED\t1\t1\t1\t2\tNA\t0.500000\t"+QUAL);
		expect(summary, "public\t1\t1\tCOBALAMIN_UNION_ON_RF00174_PARENTS\tALL_LABELLED\t2\t1\t0\t3\tNA\t0.666667\t"+QUAL);
		final ArrayList<String> references=rows(f.path("exact.exact_references.tsv"));
		check(references.size()==10, "All three labelled parents survive in each view, including exclusive masks");
		expect(references, "public\t1\t1\tMODEL_SUPPORT_RECOVERY\tA\tp%25one\tSEEN\tNOT_SEEN\ttrue\ttrue\tc\t100\t119\t+\tTP\t2");
		expect(references, "public\t1\t1\tMODEL_SUPPORT_RECOVERY\tB\tp%25one\tUNKNOWN\tSEEN\tfalse\tfalse\tc\t100\t119\t+\tTP\t2");
		final ArrayList<String> union=rows(f.path("exact.exact_union_queries.tsv"));
		check(union.size()==6, "All union queries survive, including a different-model TP, duplicate and unmatched FP");
		expect(union, "UNASSIGNED\t"+f.unionQuery("A", 3, 100, "FP", ".", true, MODEL0));
		expect(union, "UNASSIGNED\t"+f.unionQuery("A", 4, 600, "FP", ".", false, MODEL0));
		expect(union, "UNASSIGNED\t"+f.unionQuery("A", 5, 400, "TP", "other", false, MODEL1));
		expect(union, "UNKNOWN\t"+f.unionQuery("B", 2, 100, "TP", "p%25one", false, MODEL0));
		check(rows(f.path("exact.exact_subtype_queries.tsv")).size()==9, "Both RF00174 views retain all four queries, not only TPs");
		final String qqual="QUERY_ROWS_NOT_UNIQUE_LOCI;UNION_INCLUDES_ALL_COBALAMIN_MODELS";
		final ArrayList<String> queries=rows(f.path("exact.exact_query_counts.tsv"));
		expect(queries, "public\t1\t1\tCOBALAMIN_UNION_ON_RF00174_PARENTS\t5\t3\t2\t0\t0\t1\t3\t"+qqual);
		expect(queries, "public\t1\t1\tMODEL_SUPPORT_RECOVERY\t4\t2\t2\t0\t0\t1\t2\t"+qqual);
		expect(queries, "public\t1\t1\tEXCLUSIVE_SUBTYPE_ACCURACY\t4\t1\t1\t2\t0\t1\t3\t"+qqual);

		final ArrayList<String> labels=rows(f.path("strata.tsv"));
		final ArrayList<String> duplicate=new ArrayList<String>(labels);
		duplicate.add(labels.get(1));
		f.writeRows("duplicate.tsv", duplicate);
		reject(f.args("bad_duplicate", f.path("duplicate.tsv")), "Duplicate exact label");
		final ArrayList<String> missing=new ArrayList<String>(labels);
		missing.remove(3);
		f.writeRows("missing.tsv", missing);
		reject(f.args("bad_missing", f.path("missing.tsv")), "requires an exact-sequence label");
		final ArrayList<String> wrong=new ArrayList<String>(labels);
		wrong.set(1, wrong.get(1).replace("A\t", "wrong_input\t"));
		f.writeRows("wrong.tsv", wrong);
		reject(f.args("bad_wrong", f.path("wrong.tsv")), "no original input/parent");
		final ArrayList<String> badUnknown=new ArrayList<String>(labels);
		badUnknown.set(3, badUnknown.get(3).replace("UNKNOWN", "NOT_SEEN"));
		f.writeRows("unknown.tsv", badUnknown);
		reject(f.args("bad_unknown", f.path("unknown.tsv")), "must be UNKNOWN");
		final ArrayList<String> extra=new ArrayList<String>(labels);
		extra.add("A\tother\tNOT_SEEN\tNOT_SEEN\ttrue\tfalse");
		f.writeRows("extra.tsv", extra);
		reject(f.args("bad_extra", f.path("extra.tsv")), "complete scored RF00174 support population");
		reject(f.args("exact", f.path("strata.tsv")), "Output exists");
		System.err.println("NCRNA_EXACT_MEMBERSHIP_TEST_PASS");
	}

	private static final class Fixture{
		Fixture(final File root_){root=root_;}
		String path(final String name){return new File(root, name).toString();}
		String[] args(final String name, final String strata){
			final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+path("times.tsv"),
				"gradeinputs="+path("inputs.tsv"), "grades="+path("grades.tsv"), "raw="+path("raw.tsv"),
				"subtyperaw="+path("subtype_raw.tsv"), "subtypegrades="+path("subtype_summary.tsv"),
				"out="+path(name), "policies=public", "baseline=public", "workers=1", "repeats=1", "t=1"));
			if(strata!=null){args.add("exactstrata="+strata);}
			return args.toArray(new String[args.size()]);
		}
		void write(){
			try(Output out=new Output(path("times.tsv"))){out.writer.println(CONTEXT+"user_seconds\tsystem_seconds\twall_seconds\npublic\t1\t1\t12\t1\t14");}
			try(Output out=new Output(path("strata.tsv"))){
				out.writer.println(NcrnaDdlCallerPanel.STRATA_HEADER);
				out.writer.println("A\tp%25one\tSEEN\tNOT_SEEN\ttrue\ttrue");
				out.writer.println("A\tp2\tNOT_SEEN\tSEEN\ttrue\ttrue");
				out.writer.println("B\tp%25one\tUNKNOWN\tSEEN\tfalse\tfalse");
			}
			try(Output out=new Output(path("inputs.tsv"))){
				out.writer.println(CONTEXT+GradeNcrnaReference.INPUT_HEADER);
				out.writer.println(C+"A\t/A.fa\t/qA\t/rA\t4\t0\t0\t0\t4\t0");
				out.writer.println(C+"B\t/B.fa\t/qB\t/rB\t1\t0\t0\t0\t1\t0");
				out.writer.println(C+"empty\t/empty.fa\t/qE\t/rE\t0\t0\t0\t0\t0\t0");
			}
			try(Output out=new Output(path("grades.tsv"))){
				out.writer.println(CONTEXT+GradeNcrnaReference.SUMMARY_HEADER);
				out.writer.println(C+"cobalamin\t3\t3\t2\t1\t4\t0\t0\t1\t0\t3\t3\t3\t0.6\t0.75\t0.5\tFIXTURE");
			}
			try(Output out=new Output(path("raw.tsv"))){
				out.writer.println(CONTEXT+GradeNcrnaReference.RAW_HEADER);
				unionReference(out, "A", "p%25one", 100, "TP", "2", "RF00174");
				unionReference(out, "A", "p2", 200, "FN", ".", "RF00174");
				unionReference(out, "B", "p%25one", 100, "TP", "2", "RF00174");
				unionReference(out, "A", "other", 400, "TP", "5", "RF01482");
				out.writer.println(unionQuery("A", 2, 100, "TP", "p%25one", false, MODEL0));
				out.writer.println(unionQuery("A", 3, 100, "FP", ".", true, MODEL0));
				out.writer.println(unionQuery("A", 4, 600, "FP", ".", false, MODEL0));
				out.writer.println(unionQuery("A", 5, 400, "TP", "other", false, MODEL1));
				out.writer.println(unionQuery("B", 2, 100, "TP", "p%25one", false, MODEL0));
			}
			try(Output out=new Output(path("subtype_raw.tsv")); Output sums=new Output(path("subtype_summary.tsv"))){
				out.writer.println(CONTEXT+SUBTYPE_RAW);
				sums.writer.println(CONTEXT+SUBTYPE_SUMMARY);
				for(int v=0; v<2; v++){
					final String view=VIEWS[v];
					sums.writer.println(C+view+"\tRF00174\t3\t"+(v==0 ? "2\t2\t1\t3\t0\t0\t1\t2\t2\t2\t0.5\t0.666667" : "1\t1\t1\t2\t1\t2\t1\t1\t1\t1\t0.5\t0.5")+"\t0.5\tFIXTURE");
					sums.writer.println(C+view+"\tRF01482\t3\t"+(v==0 ? "1\t0\t1\t2\t0\t0\t0\t1\t1\t1\t1\t0.5" : "1\t0\t0\t1\t1\t0\t0\t1\t1\t1\t1\t1")+"\t0.5\tFIXTURE");
					sums.writer.println(C+view+"\tRF01689\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tNA\tNA\t0.5\tFIXTURE");
					subReference(out, view, "RF00174", "A", "p%25one", 100, v==0 ? "TP" : "DONT_CARE", v==0 ? "2" : ".", true);
					subReference(out, view, "RF00174", "A", "p2", 200, "FN", ".", false);
					subReference(out, view, "RF00174", "B", "p%25one", 100, "TP", "2", false);
					subReference(out, view, "RF01482", "A", "p%25one", 100, v==0 ? "FN" : "DONT_CARE", ".", true);
					subReference(out, view, "RF01482", "A", "other", 400, "TP", "5", false);
					subQuery(out, view, "RF00174", "A", 2, 100, v==0 ? "TP" : "MASKED", v==0 ? "p%25one" : ".", false, true);
					subQuery(out, view, "RF00174", "A", 3, 100, v==0 ? "FP" : "MASKED", ".", true, true);
					subQuery(out, view, "RF00174", "A", 4, 600, "FP", ".", false, false);
					subQuery(out, view, "RF00174", "B", 2, 100, "TP", "p%25one", false, false);
					subQuery(out, view, "RF01482", "A", 5, 400, "TP", "other", false, false);
				}
			}
		}
		String unionQuery(final String input, final int line, final int start, final String outcome,
				final String parent, final boolean duplicate, final String model){
			return C+input+"\tQUERY\tcobalamin\tc\t"+start+'\t'+(start+19)+"\t+\t"+outcome+"\t.\t"+line+'\t'+parent+
				"\t.\t.\t"+duplicate+"\ttrue\t"+model+'\t'+(outcome.equals("TP") ? "0\t0" : ".\t.");
		}
		void unionReference(final Output out, final String input, final String parent, final int start,
				final String outcome, final String query, final String model){
			out.writer.println(C+input+"\tREFERENCE\tcobalamin\tc\t"+start+'\t'+(start+19)+"\t+\t"+outcome+'\t'+parent+
				"\t.\t.\t"+query+"\t.\t.\t.\t"+model+"\t.\t.");
		}
		void subReference(final Output out, final String view, final String type, final String input,
				final String parent, final int start, final String outcome, final String query, final boolean multi){
			out.writer.println(C+view+'\t'+type+'\t'+input+"\tREFERENCE\t"+parent+"\t.\t"+type+"\tc\t"+start+'\t'+(start+19)+
				"\t+\t"+outcome+"\t.\t"+query+"\t.\t.\t.\t"+(multi ? "RF00174,RF01482" : type)+"\t.\t"+multi+"\t.\t.");
		}
		void subQuery(final Output out, final String view, final String type, final String input, final int line,
				final int start, final String outcome, final String parent, final boolean duplicate, final boolean multi){
			out.writer.println(C+view+'\t'+type+'\t'+input+"\tQUERY\t.\t"+line+"\tcobalamin_"+type+"_m\tc\t"+start+'\t'+(start+19)+
				"\t+\t"+outcome+'\t'+parent+"\t.\t"+(outcome.equals("MASKED") ? "p%25one" : ".")+'\t'+duplicate+
				'\t'+(!outcome.equals("MASKED"))+"\t"+(multi ? "RF00174,RF01482" : type)+"\t.\t"+multi+'\t'+(outcome.equals("TP") ? "0\t0" : ".\t."));
		}
		void writeRows(final String file, final ArrayList<String> rows){try(Output out=new Output(path(file))){for(final String row:rows){out.writer.println(row);}}}
		final File root;
	}

	private static ArrayList<String> rows(final String path){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){result.add(new String(row, StandardCharsets.UTF_8));}
		}
		return result;
	}
	private static void expect(final ArrayList<String> rows, final String expected){check(rows.contains(expected), "Missing literal saved-evidence row: "+expected);}
	private static void reject(final String[] args, final String expected){
		try{SummarizeNcrnaBenchmark.main(args);}
		catch(IllegalArgumentException failure){check(failure.getMessage().contains(expected), "Wrong rejection: "+failure.getMessage()); return;}
		throw new AssertionError("Invalid exact evidence was accepted: "+expected);
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private static final String C="public\t1\t1\t";
	private static final String CONTEXT="policy\tthreads\trepeat\t";
	private static final String MODEL0="cobalamin_RF00174_m";
	private static final String MODEL1="cobalamin_RF01482_m";
	private static final String QUAL="RF00174_LABELLED_PARENTS;PRECISION_UNAVAILABLE_UNMATCHED_QUERIES_UNPARTITIONED";
	private static final String[] VIEWS={"MODEL_SUPPORT_RECOVERY", "EXCLUSIVE_SUBTYPE_ACCURACY"};
	private static final String SUBTYPE_RAW="view\tsubtype\tinput_id\trole\tparent_id\tquery_line\tmodel\tseqid\tstart\tstop\tstrand\toutcome\tmatched_parent_id\tmatched_query_line\tmask_parent_id\tduplicate\tscored_eligible\tscored_support\tmarginal_support\tmulti_supported\td5\td3";
	private static final String SUBTYPE_SUMMARY="view\tsubtype\tinputs\tTP\tFP\tFN\tscored_references\tmasked_references\tmasked_calls\tduplicate_queries\texact_5prime\texact_3prime\texact_both\tprecision\trecall\toverlap\tqualification";
}
