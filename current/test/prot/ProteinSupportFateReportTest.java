package test.prot;

import prot.ProteinSupportFateReport;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

/** Literal fixtures for supported-stop fate joins. @author Yelan, Sayu */
public final class ProteinSupportFateReportTest{
	public static void main(String[] args) throws Exception{
		final Path root=Files.createTempDirectory("protein-support-fate-report-");
		final Path fate=root.resolve("fate.tsv"), ref=root.resolve("reference.tsv"), cand=root.resolve("candidate.tsv");
		final Path manifest=root.resolve("manifest.tsv"), out=root.resolve("report");
		write(fate, FATE_HEADER
			+fate("ctg%2C1", 0, 150, 10, 2, 1, "selected")
			+fate("ctg%2C1", 0, 250, 0, 0, 0, "no_starts")
			+fate("ctg%2C1", 1, 300, 4, 1, 0, "dp_not_selected")
			+fate("ctg%2C1", 1, 400, 5, 1, 1, "selected")
			+fate("ctg%2C1", 1, 500, 2, 0, 0, "soft_pruned"));
		write(ref, ROW_HEADER
			+row("11", "CDS", "ctg,1", 0, 100, 150, "primary", "ID=ref1")
			+row("12", "CDS", "ctg,1", 0, 200, 250, "primary", "ID=ref2")
			+row("13", "CDS", "ctg,1", 1, 300, 350, "pseudo", "ID=ref3")
			+row("14", "CDS", "ctg,1", 2, 700, 750, "unknown_strand", "ID=ref4"));
		write(cand, ROW_HEADER
			+row("21", "CDS", "ctg,1", 0, 100, 150, "primary", "bbgap_status=stop_supported;ID=call1")
			+row("22", "CDS", "ctg,1", 1, 600, 650, "primary", "bbgap_status=no_hit;ID=call2")
			+row("23", "CDS", "ctg,1", 1, 400, 450, "reference_exception_mask", "bbgap_status=stop_supported;ID=call3"));
		write(manifest, "genome\tfate\treference\tcandidate\nfixture\t"+fate+"\t"+ref+"\t"+cand+"\n");
		ProteinSupportFateReport.run(new String[]{"manifest="+manifest, "out="+out});
		final List<String> stops=Files.readAllLines(root.resolve("report.stops.tsv"), StandardCharsets.UTF_8);
		final List<String> losses=Files.readAllLines(root.resolve("report.losses.tsv"), StandardCharsets.UTF_8);
		final List<String> summary=Files.readAllLines(root.resolve("report.summary.tsv"), StandardCharsets.UTF_8);
		check(stops.size()==6 && stops.get(1).contains("ctg%2C1\t0\t150\t1\t1\tprimary"), "Escaped contigs must join to primary reference rows");
		check(stops.get(4).contains("ctg%2C1\t1\t400\t0\t0\tabsent"), "Non-primary selected rows must be conserved");
		check(losses.size()==4 && losses.get(1).contains("\tno_starts\t") && losses.get(2).contains("\tdp_not_selected\t")
			&& losses.get(3).contains("\tsoft_pruned\t"),
			"All unselected supported stops must be retained as losses");
		check(summary.size()==2 && summary.get(1).equals("fixture\t5\t2\t3\t2\t1\t1\t1\t1\t1"), "Summary counts mismatch: "+summary);
		write(fate, FATE_HEADER+fate("ctg%252C1", 0, 150, 1, 1, 1, "selected"));
		reject(manifest, root.resolve("bad-escape"), "Fate/GFF selected-key mismatch");
		write(cand, ROW_HEADER);
		write(fate, FATE_HEADER+fate("ctg%2C1", 0, 150, 1, 1, 0, "selected"));
		reject(manifest, root.resolve("bad-selected"), "Fate class disagrees with counts");
		write(fate, FATE_HEADER+badAccounting("ctg%2C1"));
		reject(manifest, root.resolve("bad-accounting"), "Fate rejected+retained must equal admitted");
		write(fate, FATE_HEADER+badPredicates("ctg%2C1"));
		reject(manifest, root.resolve("bad-predicates"), "Pruning predicates must sum");
		testChanges(root);
		System.err.println("PROTEIN_SUPPORT_FATE_REPORT_TEST_PASS escapes=true conservation=true");
	}
	private static void testChanges(Path root) throws Exception{
		final Path fate=root.resolve("change-fate.tsv"), comparison=root.resolve("comparison.tsv");
		final Path manifest=root.resolve("change-manifest.tsv");
		final String fates=FATE_HEADER
			+fate("ctg%3A%252C%2C1", 0, 100, 2, 1, 1, "selected")
			+fate("ctg%3A%252C%2C1", 1, 200, 0, 0, 0, "no_starts")
			+fate("ctg%3A%252C%2C1", 1, 300, 2, 0, 0, "soft_pruned")
			+fate("ctg%3A%252C%2C1", 0, 400, 2, 1, 0, "dp_not_selected")
			+fate("ctg%3A%252C%2C1", 0, 600, 1, 1, 1, "selected");
		// Raw selected call at100 survives but becomes non-primary: a loss is not always absence.
		final String rows=comparison(0, 100, 2, 2, 1, 1, 1, 0, "lost")
			+comparison(1, 200, 1, 1, 1, 1, 0, 0, "lost")
			+comparison(1, 300, 1, 1, 1, 1, 0, 0, "lost")
			+comparison(0, 400, 1, 1, 1, 1, 0, 0, "lost")
			+comparison(0, 500, 1, 1, 1, 1, 0, 0, "lost")
			+comparison(0, 600, 1, 1, 0, 0, 1, 1, "gained")
			+comparison(1, 700, 1, 1, 0, 0, 1, 1, "gained")
			+comparison(0, 800, 1, 0, 1, 1, 0, 0, "lost")
			+comparison(0, 900, 1, 1, 1, 1, 1, 1, "same_stop_changed_starts")
			// Literal producer shapes: unknown-strand CDS rows survive in all, never primary.
			+"stops-v1\tunknown_dot\t3\t950\t1\t0\t0\t0\t0\t0\texcluded_or_uncalled\t41\t.\t.\n"
			+"stops-v1\tunknown_question\t2\t960\t1\t0\t1\t0\t1\t0\texcluded_or_uncalled\t42\t43\t44\n";
		write(fate, fates); write(comparison, COMPARISON_HEADER+rows);
		write(manifest, "genome\tfate\tcomparison\nfixture\t"+fate+"\t"+comparison+"\n");
		ProteinSupportFateReport.run(new String[]{"mode=changes", "manifest="+manifest, "out="+root.resolve("changes")});
		final List<String> summary=Files.readAllLines(root.resolve("changes.summary.tsv"), StandardCharsets.UTF_8);
		final List<String> changes=Files.readAllLines(root.resolve("changes.changes.tsv"), StandardCharsets.UTF_8);
		check(summary.size()==3 && summary.get(1).equals("fixture\tlost\t5\t4\t1\t1\t1\t1\t1")
			&& summary.get(2).equals("fixture\tgained\t2\t1\t1\t0\t0\t0\t1"), "Unique-key fate partition mismatch: "+summary);
		check(changes.size()==8 && changes.get(1).contains("ctg%3A%252C%2C1:0:100\tctg%3A%252C%2C1\t0\t100\tlost\t2\t2\t1\t1\t1\t0\t1\t1\tselected"),
			"Escaped stop key and raw-presence loss must survive the join");
		check(changes.get(3).endsWith("\t0\t.\tuntraced\t.\t.\t.\t."), "Untraced is unknown, not zero support or a pruning cause");
		// Reordered comparison rows must produce identical sorted output.
		final String[] lines=rows.split("\n");
		final StringBuilder reverse=new StringBuilder(COMPARISON_HEADER);
		for(int i=lines.length-1; i>=0; i--){reverse.append(lines[i]).append('\n');}
		write(comparison, reverse.toString());
		ProteinSupportFateReport.run(new String[]{"mode=changes", "manifest="+manifest, "out="+root.resolve("reverse")});
		check(changes.equals(Files.readAllLines(root.resolve("reverse.changes.tsv"), StandardCharsets.UTF_8))
			&& summary.equals(Files.readAllLines(root.resolve("reverse.summary.tsv"), StandardCharsets.UTF_8)), "Input order must not affect joins");
		write(comparison, COMPARISON_HEADER+rows+lines[0]+"\n");
		rejectChanges(manifest, root.resolve("duplicate-stop"), "Duplicate comparison stop");
		write(comparison, COMPARISON_HEADER+comparison(0, 100, 1, 1, 1, 1, 0, 0, "lost"));
		rejectChanges(manifest, root.resolve("wrong-arm"), "candidate-presence mismatch");
		write(comparison, COMPARISON_HEADER+comparison(0, 100, 1, 1, 0, 0, 1, 1, "lost"));
		rejectChanges(manifest, root.resolve("wrong-change"), "Change class disagrees");
		write(comparison, COMPARISON_HEADER+comparison(0, 100, 1, 2, 0, 0, 1, 1, "gained"));
		rejectChanges(manifest, root.resolve("wrong-count"), "Invalid raw/primary");
		write(comparison, COMPARISON_HEADER+"stops-v1\tunknown_dot\t3\t950\t1\t1\t0\t0\t0\t0\texcluded_or_uncalled\t41\t.\t.\n");
		rejectChanges(manifest, root.resolve("unknown-reference-primary"), "Primary comparison rows require a known strand");
		write(comparison, COMPARISON_HEADER+"stops-v1\tunknown_dot\t3\t950\t1\t0\t1\t1\t0\t0\tlost\t41\t43\t.\n");
		rejectChanges(manifest, root.resolve("unknown-baseline-primary"), "Primary comparison rows require a known strand");
		write(comparison, COMPARISON_HEADER+"stops-v1\tunknown_dot\t3\t950\t1\t0\t0\t0\t1\t1\tgained\t41\t.\t44\n");
		rejectChanges(manifest, root.resolve("unknown-candidate-primary"), "Primary comparison rows require a known strand");
		write(comparison, COMPARISON_HEADER);
		rejectChanges(manifest, root.resolve("missing-selected"), "Selected fate stop missing");
		write(fate, FATE_HEADER);
		ProteinSupportFateReport.run(new String[]{"mode=changes", "manifest="+manifest, "out="+root.resolve("empty")});
		check(Files.readAllLines(root.resolve("empty.summary.tsv"), StandardCharsets.UTF_8).get(1)
			.equals("fixture\tlost\t0\t0\t0\t0\t0\t0\t0"), "Empty valid inputs must emit explicit zero partitions");
		System.err.println("PROTEIN_SUPPORT_FATE_CHANGES_TEST_PASS partitions=true raw_presence=true escaped_keys=true rejection=true order=true unknown_strand_excluded=true");
	}
	private static String comparison(int strand, int terminal, int refRaw, int refPrimary,
		int oldRaw, int oldPrimary, int newRaw, int newPrimary, String change){
		return "stops-v1\tctg:%2C,1\t"+strand+"\t"+terminal+"\t"+refRaw+"\t"+refPrimary+"\t"+oldRaw+"\t"+oldPrimary
			+"\t"+newRaw+"\t"+newPrimary+"\t"+change+"\t11\t21\t31\n";
	}
	private static void rejectChanges(Path manifest, Path out, String diagnostic){
		boolean failed=false;
		try{ProteinSupportFateReport.run(new String[]{"mode=changes", "manifest="+manifest, "out="+out});}
		catch(IllegalArgumentException expected){check(expected.getMessage().contains(diagnostic), "Wrong diagnostic: "+expected); failed=true;}
		check(failed, "Expected changes rejection for "+diagnostic);
		check(!Files.exists(out.resolveSibling(out.getFileName()+".changes.tsv")), "Invalid input must fail before output creation");
	}
	private static String fate(String contig, int strand, int terminal, int admitted, int retained, int selected, String fate){
		final int pruned=admitted-retained;
		return "support-fate-v1\t"+contig+":"+strand+":"+terminal+"\t"+contig+"\t"+strand+"\t"+terminal
			+"\t100\t"+terminal+"\t0\t0\t1\trep\t50.0\t100\t50\t80\t-0.1\t1.25\t"+admitted
			+"\t0\t0\t0\t"+pruned+"\t"+retained+"\t"+selected+"\t"+fate
			+"\t"+(pruned==0 ? "." : "first:"+pruned)+"\t100\t"+terminal+"\t1\t2\t3\t4\t50\t54"
			+"\tbest_retained\t100\t"+terminal+"\t60\t30\t.\n";
	}
	private static String badAccounting(String contig){
		return "support-fate-v1\t"+contig+":0:150\t"+contig+"\t0\t150\t100\t150\t0\t0\t1\trep\t50.0\t100\t50\t80\t-0.1\t1.25"
			+"\t2\t0\t0\t0\t1\t0\t0\tsoft_pruned\tfirst:1\t100\t150\t1\t2\t3\t4\t50\t54"
			+"\tbest_admitted\t100\t150\t60\t30\t.\n";
	}
	private static String badPredicates(String contig){
		return "support-fate-v1\t"+contig+":0:150\t"+contig+"\t0\t150\t100\t150\t0\t0\t1\trep\t50.0\t100\t50\t80\t-0.1\t1.25"
			+"\t2\t0\t0\t0\t2\t0\t0\tsoft_pruned\tfirst:1\t100\t150\t1\t2\t3\t4\t50\t54"
			+"\tbest_admitted\t100\t150\t60\t30\t.\n";
	}
	private static String row(String id, String type, String contig, int strand, int start, int end, String stratum, String attributes){
		return "stops-v1\t"+id+"\t"+type+"\t"+contig+"\t"+strand+"\t"+start+"\t"+end+"\t0\t11\t0\t0\t"+stratum+"\t"+attributes+"\n";
	}
	private static void reject(Path manifest, Path out, String diagnostic){
		boolean failed=false;
		try{ProteinSupportFateReport.run(new String[]{"manifest="+manifest, "out="+out});}
		catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(diagnostic), "Wrong diagnostic: "+expected); failed=true;
		}
		check(failed, "Expected report rejection for "+diagnostic);
	}
	private static void write(Path path, String text) throws Exception{Files.write(path, text.getBytes(StandardCharsets.UTF_8));}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private static final String FATE_HEADER="schema\tstop_id\tcontig\tstrand\tterminal\tlongest_start\tlongest_end\tpack_index\tmodel_rank\tfamily_id\trepresentative\tidentity\tscore64\ttail_length\tmin_length\tmin_start_score\traw_stop_score\tadmitted\trejected_start_score\trejected_length\trejected_limit\tsoft_rejected\tretained\tselected\tfate\tfirst_pruning_predicates\tbest_start\tbest_end\tbest_raw_start\tbest_raw_stop\tbest_inner\tbest_native\tbest_bonus\tbest_total\toverlap_basis\toverlap_start\toverlap_end\tmax_overlap_same\tmax_overlap_opposite\tobserved_selected_overlaps\n";
	private static final String ROW_HEADER="schema\tinput_row\ttype\tcontig\tstrand\tstart\tend\tphase\ttable\tpartial5\tpartial3\tstratum\tattributes\n";
	private static final String COMPARISON_HEADER="schema\tcontig\tstrand\tterminal\treference_raw_rows\treference_primary_rows\tbaseline_raw_rows\tbaseline_primary_rows\tcandidate_raw_rows\tcandidate_primary_rows\tchange\treference_rows\tbaseline_rows\tcandidate_rows\n";
	private ProteinSupportFateReportTest(){}
}
