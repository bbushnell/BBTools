package test.prot;

import prot.ReferenceConcordanceReport;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.TreeMap;

/** Literal round-trip and rejection fixtures for the generic report. @author Sayu */
public final class ReferenceConcordanceReportTest{
	public static void main(String[] args) throws Exception{
		final Path root=Files.createTempDirectory("reference-concordance-test-");
		final Path summary=root.resolve("summary.tsv"), stops=root.resolve("stops.tsv"), ta=root.resolve("a.time"), tb=root.resolve("b.time");
		final Path manifest=root.resolve("manifest.tsv"), output=root.resolve("report");
		final String summaryText="schema\tarm\tmetric\tvalue\n"+metrics("reference", true)+metrics("baseline", false)+metrics("candidate", false)
			+"stops-v1\tdelta\tshared_primary_stops\t2\nstops-v1\tdelta\tgained_primary_stops\t1\nstops-v1\tdelta\tlost_primary_stops\t1\n"
			+"stops-v1\tdelta\tgained_reference_stops\t1\nstops-v1\tdelta\tlost_reference_stops\t1\nstops-v1\tdelta\tsame_stop_changed_starts\t1\n";
		final String stopsText="schema\tcontig\tstrand\tterminal\treference_raw_rows\treference_primary_rows\tbaseline_raw_rows\tbaseline_primary_rows\tcandidate_raw_rows\tcandidate_primary_rows\tchange\treference_rows\tbaseline_rows\tcandidate_rows\n"
			+"stops-v1\tc\t0\t10\t1\t1\t1\t1\t1\t1\tunchanged\t1\t1\t1\n"
			+"stops-v1\tc\t0\t20\t1\t1\t1\t1\t0\t0\tlost\t2\t2\t.\n"
			+"stops-v1\tc\t1\t30\t1\t1\t0\t0\t1\t1\tgained\t3\t.\t2\n"
			+"stops-v1\tc\t1\t40\t1\t0\t1\t0\t1\t0\texcluded_or_uncalled\t4\t3\t3\n"
			+"stops-v1\tc\t1\t50\t0\t0\t1\t1\t1\t1\tsame_stop_changed_starts\t.\t4\t4\n";
		final String manifestHeader="genome\tbaseline\tcandidate\tsummary\tstops\tbaseline_timing\tcandidate_timing\n";
		final String manifestRow="toy\toff\tnew\t"+summary+"\t"+stops+"\t"+ta+"\t"+tb+"\n";
		write(summary, summaryText); write(stops, stopsText); write(ta, "elapsed_seconds\t1.25\n"); write(tb, "elapsed_seconds\t2.5\n");
		write(manifest, manifestHeader+manifestRow);
		ReferenceConcordanceReport.run(new String[]{"in="+manifest, "out="+output});
		final List<String> arms=Files.readAllLines(root.resolve("report.arms.tsv"), StandardCharsets.UTF_8);
		check(arms.size()==3 && arms.get(1).startsWith("toy\toff\t3\t3\t2\t1\t0.6666666666666666\t0.6666666666666666\t1.25\t4\t1\t"), "Literal arm counts/rates/caller time mismatch");
		final List<String> pairs=Files.readAllLines(root.resolve("report.comparisons.tsv"), StandardCharsets.UTF_8);
		check(pairs.size()==2 && pairs.get(1).startsWith("toy\toff\tnew\t2\t1\t1\t1\t1\t1\t"), "Literal paired deltas mismatch");
		final List<String> changes=Files.readAllLines(root.resolve("report.changes.tsv"), StandardCharsets.UTF_8);
		check(changes.size()==4 && changes.get(1).contains("\tlost\t") && changes.get(2).contains("\tgained\t")
			&& changes.get(3).contains("\tsame_stop_changed_starts\t"), "Every changed stop must remain linked to original rows");
		check(Files.readAllLines(root.resolve("report.strata.tsv"), StandardCharsets.UTF_8).size()==7, "Reference and both arm strata must survive");
		rejectExisting(manifest, output); // Existing outputs must never be silently replaced.
		write(summary, summaryText+"stops-v1\tbaseline\tprimary_stops\t3\n"); reject(manifest, root.resolve("duplicate"), "Duplicate/empty summary metric");
		write(summary, summaryText.replace("stops-v1\tbaseline\tmatched_stops\t2\n", "")); reject(manifest, root.resolve("missing"), "Missing metric: matched_stops");
		write(summary, summaryText.replace("0.6666666666666666", "0.8")); reject(manifest, root.resolve("bad-rate"), "Rate disagrees with explicit numerator/denominator");
		write(summary, summaryText);
		write(stops, stopsText.replace("\tlost\t", "\tgained\t")); reject(manifest, root.resolve("bad-change"), "Change label contradicts primary presence");
		write(stops, stopsText+stopsText.substring(stopsText.lastIndexOf("stops-v1"))); reject(manifest, root.resolve("duplicate-stop"), "Duplicate paired-stop key");
		write(stops, stopsText);
		write(manifest, manifestHeader+manifestRow+manifestRow); reject(manifest, root.resolve("duplicate-pair"), "Duplicate comparison");
		write(manifest, manifestHeader+manifestRow); write(ta, "elapsed_seconds\tNaN\n"); reject(manifest, root.resolve("bad-time"), "Invalid caller time");
		write(ta, "elapsed_seconds\t1.25\n");
		// A second comparison may share an arm only when its complete metrics agree.
		final Path other=root.resolve("other.tsv");
		write(other, summaryText.replace("baseline\trows.reference_exception_mask", "baseline\trows.synthetic_stop"));
		write(manifest, manifestHeader+manifestRow+"toy\toff\tthird\t"+other+"\t"+stops+"\t"+ta+"\t"+tb+"\n");
		reject(manifest, root.resolve("arm-conflict"), "Repeated genome/arm has conflicting evidence");
		write(other, summaryText.replace("reference\trows.pseudo", "reference\trows.multipart"));
		reject(manifest, root.resolve("reference-conflict"), "Reference inventory differs between comparisons");
		final TreeMap<String,String> empty=new TreeMap<String,String>();
		for(String key : new String[]{"raw_cds_rows", "primary_rows", "excluded_rows", "primary_stops", "reference_primary_stops", "matched_stops", "unmatched_primary_stops", "missed_reference_stops"}){empty.put(key, "0");}
		empty.put("stop_reference_precision", "NA"); empty.put("stop_reference_recall", "NA");
		ReferenceConcordanceReport.validateMetrics(empty);
		empty.put("stop_reference_precision", "0");
		boolean failed=false; try{ReferenceConcordanceReport.validateMetrics(empty);}catch(IllegalArgumentException expected){
			check(expected.getMessage().contains("Zero denominator requires NA"), "Wrong zero-denominator diagnostic: "+expected); failed=true;
		}
		check(failed, "Zero denominator must not become numeric zero precision");
		System.err.println("REFERENCE_CONCORDANCE_REPORT_TEST_PASS literal=true rejections=11 zero_denominators=true fixtures="+root);
	}
	private static String metrics(String arm, boolean ref){
		final StringBuilder out=new StringBuilder();
		final String[][] pairs={{"raw_cds_rows", "4"}, {"primary_rows", "3"}, {"excluded_rows", "1"}, {"primary_stops", "3"},
			{"reference_primary_stops", "3"}, {"matched_stops", ref ? "3" : "2"}, {"unmatched_primary_stops", ref ? "0" : "1"},
			{"missed_reference_stops", ref ? "0" : "1"}, {"stop_reference_precision", ref ? "1.0" : "0.6666666666666666"},
			{"stop_reference_recall", ref ? "1.0" : "0.6666666666666666"}, {"rows.primary", "3"}, {ref ? "rows.pseudo" : "rows.reference_exception_mask", "1"}};
		for(String[] pair : pairs){out.append("stops-v1\t").append(arm).append('\t').append(pair[0]).append('\t').append(pair[1]).append('\n');}
		return out.toString();
	}
	private static void reject(Path manifest, Path prefix, String diagnostic){
		boolean failed=false;
		try{ReferenceConcordanceReport.run(new String[]{"in="+manifest, "out="+prefix});}catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(diagnostic), "Wrong rejection diagnostic for "+prefix+": "+expected); failed=true;
		}
		check(failed, "Expected malformed/conflicting input rejection: "+prefix);
		for(String suffix : new String[]{".arms.tsv", ".comparisons.tsv", ".strata.tsv", ".changes.tsv"}){
			check(!Files.exists(prefix.resolveSibling(prefix.getFileName()+suffix)), "Invalid input must fail before output: "+prefix);
		}
	}
	private static void rejectExisting(Path manifest, Path prefix) throws Exception{
		final Path arms=prefix.resolveSibling(prefix.getFileName()+".arms.tsv"); final byte[] before=Files.readAllBytes(arms);
		boolean failed=false;
		try{ReferenceConcordanceReport.run(new String[]{"in="+manifest, "out="+prefix});}catch(AssertionError expected){
			check(expected.getMessage().contains("exists and overwrite=false"), "Wrong existing-output assertion: "+expected); failed=true;
		}
		check(failed && java.util.Arrays.equals(before, Files.readAllBytes(arms)), "Existing-output rejection must preserve bytes");
	}
	private static void write(Path path, String text) throws Exception{Files.write(path, text.getBytes(StandardCharsets.UTF_8));}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private ReferenceConcordanceReportTest(){}
}
