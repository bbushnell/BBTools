package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import parse.LineParser1;
import prok.CompareCobalaminBenchmarkGrades;
import prok.ReportNcrnaBoundaryComparison;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Literal five-genome oracle: a retained label switch, short overlapping call,
 * absent call, an eligible call assigned to another reference, and a gained
 * first call. Native matches are declared; the comparison must never rematch.
 * @author Raiden
 */
public final class CompareCobalaminBenchmarkGradesTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, CompareCobalaminBenchmarkGradesTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown comparison-fixture argument: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && new File(parser.out1).mkdir(), "Supply out=fresh_fixture_directory");
			new CompareCobalaminBenchmarkGradesTest(parser.out1).run();
			pp.outstream.println("COBALAMIN_BENCHMARK_COMPARE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	private CompareCobalaminBenchmarkGradesTest(final String directory){root=new File(directory).getAbsolutePath();}

	private void run() throws Exception{
		final String genomes="input_id\tcalls\tmatched\tmasked\tunmatched\n";
		write("before.genomes.tsv", genomes+"a\t1\t1\t0\t0\nb\t1\t1\t0\t0\nc\t1\t1\t0\t0\nd\t1\t1\t0\t0\ne\t0\t0\t0\t0\n");
		write("after.genomes.tsv", genomes+"a\t1\t1\t0\t0\nb\t1\t0\t0\t1\nc\t0\t0\t0\t0\nd\t1\t1\t0\t0\ne\t1\t1\t0\t0\n");
		final String refs="input_id\treference_id\tseqid\tstart\tstop\tstrand\tscored\n"
			+"a\tr\tshared\t100\t199\t+\ttrue\n"
			+"b\tr\tshared\t380107\t380293\t-\ttrue\n"
			+"c\tr\tshared\t100\t199\t+\ttrue\n"
			+"d\tr\tshared\t100\t199\t+\ttrue\n"
			+"d\tr2\tshared\t110\t209\t+\ttrue\n"
			+"e\tr\tshared\t100\t199\t+\ttrue\n";
		write("before.references.tsv", refs);
		write("after.references.tsv", refs);
		final String queries="input_id\tseqid\tstart\tstop\tstrand\temitted_member\tmodel\toutcome\tscored_eligible\tmatched_reference_id\n";
		write("before.queries.tsv", queries
			+"a\tshared\t100\t199\t+\tRF00174\tmain_a\tTP\ttrue\tr\n"
			+"b\tshared\t380162\t380263\t-\tRF01689\tremoved_b n=386\tTP\ttrue\tr\n"
			+"c\tshared\t100\t199\t+\tRF01482\tremoved_c\tTP\ttrue\tr\n"
			+"d\tshared\t100\t199\t+\tRF00174\tmain_d n=42\tTP\ttrue\tr\n");
		final String after=queries
			+"a\tshared\t100\t199\t+\tRF01482\tminor_a\tTP\ttrue\tr\n"
			+"b\tshared\t380175\t380267\t-\tRF01689\tshort_b\tUNMATCHED_OTHER\tfalse\t.\n"
			+"d\tshared\t110\t209\t+\tRF00174\tmain_d\tTP\ttrue\tr2\n"
			+"e\tshared\t100\t199\t+\tRF00174\tmain_e\tTP\ttrue\tr\n";
		write("after.queries.tsv", after);
		final String witness="input_id\treference_id\temitted_member\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3\n";
		write("before.witnesses.tsv", witness
			+"a\tr\tRF00174\tmain_a\tshared\t100\t199\t100\t199\t+\t0\t0\n"
			+"b\tr\tRF01689\tremoved_b n=386\tshared\t380162\t380263\t380107\t380293\t-\t30\t-55\n"
			+"c\tr\tRF01482\tremoved_c\tshared\t100\t199\t100\t199\t+\t0\t0\n"
			+"d\tr\tRF00174\tmain_d n=42\tshared\t100\t199\t100\t199\t+\t0\t0\n");
		write("after.witnesses.tsv", witness
			+"a\tr\tRF01482\tminor_a\tshared\t100\t199\t100\t199\t+\t0\t0\n"
			+"d\tr2\tRF00174\tmain_d\tshared\t110\t209\t110\t209\t+\t0\t0\n"
			+"e\tr\tRF00174\tmain_e\tshared\t100\t199\t100\t199\t+\t0\t0\n");
		write("main.fa", ">main_a\nACGT\n>main_d\nACGT\n>main_e\nACGT\n");
		write("minor.fa", ">minor_a\nACGT\n");
		write("third.fa", ">short_b\nACGT\n");
		invoke("after", "valid");
		final String summary=read("valid.summary.tsv");
		check(summary.startsWith("#COBALAMIN_COMPARISON_COUNTS_V2\nmeasure\tstate\tcount\tqualification\n"), "Comparison counts require an explicit version and long-form schema");
		check(summary.contains("graded_scored_references\tRETAINED\t1\t") && summary.contains("graded_scored_references\tLOST\t3\t")
			&& summary.contains("graded_scored_references\tGAINED\t2\t"), "One retained, three lost and two gained reference IDs must not collapse to a net change");
		check(summary.contains("lost_references_without_any_candidate_overlap\tLOST\t1\t"), "Only the genuinely absent call belongs in the no-overlap add-back class");
		check(summary.contains("genomes_with_any_clan_call\tRETAINED\t3\t") && summary.contains("genomes_with_any_clan_call\tLOST\t1\t")
			&& summary.contains("genomes_with_scored_reference_detection\tRETAINED\t2\t") && summary.contains("genomes_with_scored_reference_detection\tLOST\t2\t"),
			"Raw first/last call and scored detection presence are different quantities");
		check(summary.contains("RF00174_to_RF01482\tRETAINED\t1\t") && summary.contains("genomes_with_any_clan_call\tABSENT_BOTH\t0\t"),
			"Retained locus label switches and true zero populations must survive attribution");
		check(!summary.contains("\tNA") && !summary.contains("\tnull"), "Inapplicable summary states are absent, never numeric placeholders");
		final String loci=read("valid.loci.tsv");
		check(loci.contains("b\tr\tshared\t380107\t380293\t-\tLOST_GRADED_REFERENCE\tOVERLAPPING_CLAN_CALL_PRESENT\tRF01689\tremoved_b n=386\tfalse\t"), "A removed baseline winner with a description and shorter overlapping candidate is not an absent-call loss");
		check(loci.contains("d\tr\tshared\t100\t199\t+\tLOST_GRADED_REFERENCE\tOVERLAPPING_CLAN_CALL_PRESENT\tRF00174\tmain_d n=42\ttrue\t"), "A retained model with description must use its native first-token ID for subset membership");
		check(loci.contains("c\tr\tshared\t100\t199\t+\tLOST_GRADED_REFERENCE\tNO_OVERLAPPING_CLAN_CALL\tRF01482\tremoved_c\tfalse\t"), "Missing-call loss retains baseline member/model and actual subset membership");
		final String overlaps=read("valid.overlaps.tsv");
		boolean shortOverlap=false;
		final LineParser1 fields=new LineParser1('\t');
		for(final String line:overlaps.split("\n")){
			fields.set(line.getBytes(StandardCharsets.UTF_8));
			if(!fields.termEquals("b", 0) || !fields.termEquals("1", 3)){continue;}
			check(fields.parseInt(10)==93 && fields.parseInt(11)==187 && fields.parseInt(12)==93
				&& Math.abs(fields.parseDouble(13)-93.0/187)<5e-9 && fields.parseDouble(14)==1
				&& fields.termEquals("BELOW_RECIPROCAL_HALF", 15), "93/187 must stay below one half regardless of numeric text formatting");
			shortOverlap=true;
		}
		check(shortOverlap, "Expected the short overlapping candidate witness");
		check(overlaps.contains("d\tr\tLOST_GRADED_REFERENCE\t1\tRF00174\tmain_d\t110\t209\tTP\tr2\t90\t100\t100\t"), "Eligible geometry assigned elsewhere must be visible without a second matching pass");
		check(read("valid.lengths.tsv").contains("1\tb\tRF01689\tshort_b\tUNMATCHED_OTHER\tr\tSCORED\t93\t187\t93\t"), "Unmatched third-member lengths must not vanish behind TP-only endpoint statistics");
		changesAudit();
		pairedCases();
		genericFamily();
		unmatchedLengths();
		reportFormats();
		countReportFormats();
		write("bad.queries.tsv", after.replace("UNMATCHED_OTHER", "UNRECOGNIZED"));
		write("bad.genomes.tsv", read("after.genomes.tsv"));
		write("bad.references.tsv", refs);
		write("bad.witnesses.tsv", read("after.witnesses.tsv"));
		boolean rejected=false;
		try{invoke("bad", "bad.out");}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("Unknown normalized");}
		check(rejected && !new File(path("bad.out.summary.tsv")).exists(), "Unknown outcomes must fail before a completed comparison exists");
	}
	private void invoke(final String candidate, final String output){
		CompareCobalaminBenchmarkGrades.main(new String[]{"in1="+path("before"), "in2="+path(candidate), "out="+path(output),
			"subset="+path("main.fa")+','+path("minor.fa")+','+path("third.fa")});
	}

	/** Six candidate unmatched calls exercise all old-outcome dispositions. The
	 * retained TP has an extra nearby call, so overlap must not invent ancestry.
	 * Repeated contig names across inputs and a minus-strand shrink stay separate.
	 */
	private void changesAudit() throws Exception{
		write("audit_before.genomes.tsv", read("before.genomes.tsv")+"f\t1\t0\t1\t0\ng\t1\t0\t0\t1\nh\t0\t0\t0\t0\ni\t2\t0\t1\t1\n");
		write("audit_after.genomes.tsv", read("after.genomes.tsv").replace("a\t1\t1\t0\t0", "a\t2\t1\t0\t1")
			+"f\t1\t0\t0\t1\ng\t1\t0\t0\t1\nh\t1\t0\t0\t1\ni\t1\t0\t0\t1\n");
		for(final String side:new String[]{"before", "after"}){
			write("audit_"+side+".references.tsv", read(side+".references.tsv"));
			write("audit_"+side+".witnesses.tsv", read(side+".witnesses.tsv"));
		}
		write("audit_before.queries.tsv", read("before.queries.tsv")
			+"f\tshared\t100\t199\t+\tRF00174\tf\tMASKED\tfalse\t.\n"
			+"g\tshared\t100\t199\t+\tRF00174\tg\tUNMATCHED_OTHER\tfalse\t.\n"
			+"i\tshared\t100\t149\t+\tRF00174\ti1\tMASKED\tfalse\t.\n"
			+"i\tshared\t200\t249\t+\tRF00174\ti2\tUNMATCHED_OTHER\tfalse\t.\n");
		write("audit_after.queries.tsv", read("after.queries.tsv")
			+"a\tshared\t150\t180\t+\tRF01689\textra\tUNMATCHED_OTHER\tfalse\t.\n"
			+"f\tshared\t90\t209\t+\tRF00174\tf\tUNMATCHED_OTHER\tfalse\t.\n"
			+"g\tshared\t110\t209\t+\tRF00174\tg\tUNMATCHED_OTHER\tfalse\t.\n"
			+"h\tshared\t100\t199\t+\tRF00174\th\tUNMATCHED_OTHER\tfalse\t.\n"
			+"i\tshared\t120\t230\t+\tRF00174\ti\tUNMATCHED_OTHER\tfalse\t.\n");
		CompareCobalaminBenchmarkGrades.main(new String[]{"in1="+path("audit_before"), "in2="+path("audit_after"), "out="+path("audit"), "changesonly=t"});
		final String counts=read("audit.change_summary.tsv");
		for(final String expected:new String[]{
			"candidate_unmatched\tONLY_OLD_LOST_TP\tSHRUNK\t1\t",
			"candidate_unmatched\tONLY_OLD_RETAINED_TP\tSHRUNK\t1\t",
			"candidate_unmatched\tONLY_OLD_MASKED\tGROWN\t1\t",
			"candidate_unmatched\tONLY_OLD_UNMATCHED\tSHIFTED_EQUAL_LENGTH\t1\t",
			"candidate_unmatched\tMIXED_OLD_OUTCOMES\tMULTIPLE_LINKS\t1\t",
			"candidate_unmatched\tNO_OLD_CALL_OVERLAP\tNO_LINK\t1\t",
			"lost_reference\tLOST_SCORED_MATCH\tSHRUNK\t1\t",
			"lost_reference\tLOST_SCORED_MATCH\tSHIFTED_EQUAL_LENGTH\t1\t",
			"lost_reference\tLOST_SCORED_MATCH\tNO_LINK\t1\t",
			"after\tUNMATCHED\tALL\t6\t"}){
			check(counts.contains(expected), "Literal saved-outcome audit disagrees: "+expected);
		}
		check(read("audit.unmatched.tsv").split("\n").length==7 && read("audit.unmatched_links.tsv").split("\n").length==7,
			"Six unmatched calls conserve one row each; the mixed call keeps both links while the new region has none");
		check(read("audit.lost.tsv").split("\n").length==4, "Three lost references remain separate from the six unmatched calls");
		check(read("audit.unmatched_links.tsv").contains("\t31\t-69\t50\t-19\tSHRUNK\t31\t2\t"),
			"An extra call near a retained TP reports its geometry, original reference overlap and reverse degree two");
		check(!new File(path("audit.summary.tsv")).exists(), "changesonly must not repeat the previously completed comparison reports");
	}

	/** Literal coordinates give independent signed errors, populations and sums
	 * of squares. Lost references cannot enter a paired endpoint denominator.
	 */
	private void pairedCases() throws Exception{
		CompareCobalaminBenchmarkGrades.main(new String[]{"in1="+path("before"), "in2="+path("after"), "out="+path("paired"),
			"radius=6", "subset="+path("main.fa")+','+path("minor.fa")+','+path("third.fa")});
		for(final String suffix:new String[]{"summary", "genomes", "loci", "overlaps", "lengths", "call_lengths"}){
			check(read("valid."+suffix+".tsv").equals(read("paired."+suffix+".tsv")), "Optional paired reporting must preserve every historical report byte: "+suffix);
		}
		final LineParser1 fields=pairedRow("paired", "cobalamin", "5prime");
		check(fields.parseLong(3)==1 && fields.parseDouble(4)==0 && fields.parseDouble(5)==0
			&& fields.parseLong(8)==1 && fields.parseLong(6)==0 && fields.parseLong(7)==0, "Unchanged shared references retain zero RMS and explicit better/worse/tied counts");
		final LineParser1 absent=tableRow("paired", "paired_population", "RF01482", "5prime", 5);
		check(absent.parseLong(3)==0 && absent.parseLong(4)==0 && !read("paired.paired_endpoints.tsv").contains("\nRF01482\t")
			&& !read("paired.paired_near_endpoints.tsv").contains("\nRF01482\t") && !read("paired.paired_details.tsv").contains("\nRF01482\t"),
			"ON's switched member cannot create an OFF-cohort observation: population zero remains but nonexistent statistics are omitted");
		check(!read("paired.paired_endpoints.tsv").contains("better_fraction"), "The derived ratio is not part of the version-2 numeric schema");
		final String moreGenomes="f\t1\t1\t0\t0\ng\t1\t1\t0\t0\nh\t1\t1\t0\t0\ni\t1\t1\t0\t0\n";
		final String moreReferences="f\tr\tshared\t100\t199\t+\ttrue\ng\tr\tshared\t100\t199\t-\ttrue\n"
			+"h\tr\tshared\t100\t199\t+\ttrue\ni\tr\tshared\t100\t199\t+\ttrue\n";
		for(final String side:new String[]{"before", "after"}){
			write(side+"_more.genomes.tsv", read(side+".genomes.tsv")+moreGenomes);
			write(side+"_more.references.tsv", read(side+".references.tsv")+moreReferences);
		}
		write("before_more.queries.tsv", read("before.queries.tsv")
			+"f\tshared\t96\t207\t+\tRF00174\tf\tTP\ttrue\tr\n"
			+"g\tshared\t97\t205\t-\tRF01482\tg\tTP\ttrue\tr\n"
			+"h\tshared\t80\t229\t+\tRF01689\th\tTP\ttrue\tr\n"
			+"i\tshared\t90\t211\t+\tRF01689\ti\tTP\ttrue\tr\n");
		write("after_more.queries.tsv", read("after.queries.tsv")
			+"f\tshared\t98\t201\t+\tRF00174\tf\tTP\ttrue\tr\n"
			+"g\tshared\t90\t196\t-\tRF01689\tg\tTP\ttrue\tr\n"
			+"h\tshared\t99\t209\t+\tRF01689\th\tTP\ttrue\tr\n"
			+"i\tshared\t89\t212\t+\tRF01689\ti\tTP\ttrue\tr\n");
		write("before_more.witnesses.tsv", read("before.witnesses.tsv")
			+"f\tr\tRF00174\tf\tshared\t96\t207\t100\t199\t+\t-4\t8\n"
			+"g\tr\tRF01482\tg\tshared\t97\t205\t100\t199\t-\t-6\t3\n"
			+"h\tr\tRF01689\th\tshared\t80\t229\t100\t199\t+\t-20\t30\n"
			+"i\tr\tRF01689\ti\tshared\t90\t211\t100\t199\t+\t-10\t12\n");
		write("after_more.witnesses.tsv", read("after.witnesses.tsv")
			+"f\tr\tRF00174\tf\tshared\t98\t201\t100\t199\t+\t-2\t2\n"
			+"g\tr\tRF01689\tg\tshared\t90\t196\t100\t199\t-\t3\t10\n"
			+"h\tr\tRF01689\th\tshared\t99\t209\t100\t199\t+\t-1\t10\n"
			+"i\tr\tRF01689\ti\tshared\t89\t212\t100\t199\t+\t-11\t13\n");
		CompareCobalaminBenchmarkGrades.main(new String[]{"in1="+path("before_more"), "in2="+path("after_more"), "out="+path("more"), "radius=6"});
		final LineParser1 five=pairedRow("more", "cobalamin", "5prime"), three=pairedRow("more", "cobalamin", "3prime");
		check(five.parseLong(3)==5 && three.parseLong(3)==5 && read("more.summary.tsv").contains("graded_scored_references\tRETAINED\t5\t"),
			"Only five shared reference IDs enter RMS; the three losses and two gains remain separate");
		close(five.parseDouble(4), Math.sqrt(552.0/5));
		close(five.parseDouble(5), Math.sqrt(135.0/5));
		close(three.parseDouble(4), Math.sqrt(1117.0/5));
		close(three.parseDouble(5), Math.sqrt(373.0/5));
		check(five.parseLong(6)==3 && five.parseLong(7)==1 && five.parseLong(8)==1
			&& three.parseLong(6)==2 && three.parseLong(7)==2 && three.parseLong(8)==1, "Absolute-error changes must retain better, worse and unchanged ends separately");
		check(five.parseLong(9)==0 && five.parseLong(10)==1 && three.parseLong(9)==1 && three.parseLong(10)==1,
			"Worsening must split by OFF error inside versus beyond radius, without hiding the tail");
		final LineParser1 near5=tableRow("more", "paired_near_endpoints", "cobalamin", "5prime", 9);
		final LineParser1 near3=tableRow("more", "paired_near_endpoints", "cobalamin", "3prime", 9);
		check(near5.parseLong(4)==3 && near3.parseLong(4)==2, "The separate secondary cohort is selected only from OFF errors");
		close(near5.parseDouble(5), Math.sqrt(52.0/3));
		close(near5.parseDouble(6), Math.sqrt(13.0/3));
		close(near3.parseDouble(5), Math.sqrt(4.5));
		close(near3.parseDouble(6), Math.sqrt(50));
		check(five.parseLong(11)==2 && pairedRow("more", "RF01482", "5prime").parseLong(3)==1,
			"Label switches remain explicit while member cohorts use baseline labels");
		check(read("more.paired_loci.tsv").contains("g\tr\tshared\t-\tRF01482\tRF01689\t-6\t3\t3\t10\ttrue\ttrue\n"),
			"Reverse-strand signed errors must derive from gene-oriented endpoints, not ascending coordinate order");
		final LineParser1 detail5=tableRow("more", "paired_details", "cobalamin", "5prime", 18);
		final LineParser1 detail3=tableRow("more", "paired_details", "cobalamin", "3prime", 18);
		check(detail5.parseLong(3)==5 && detail5.parseLong(4)==1 && detail5.parseLong(5)==1
			&& detail5.parseLong(6)==1 && detail5.parseLong(7)==4
			&& detail3.parseLong(6)==2 && detail3.parseLong(7)==2,
			"Paired exact/within-three counts use the same five references, including error3 and excluding error4");
		close(detail5.parseDouble(8), Math.sqrt(552.0/5));
		close(detail5.parseDouble(9), Math.sqrt(135.0/5));
		close(detail5.parseDouble(10), 8);
		close(detail5.parseDouble(11), 3.4);
		close(detail3.parseDouble(10), 10.6);
		close(detail3.parseDouble(11), 7);
		check(detail5.parseLong(12)==4 && detail5.parseLong(13)==3 && detail5.parseLong(14)==1 && detail5.parseLong(15)==1
			&& detail5.parseLong(16)==2 && tableRow("more", "paired_details", "RF01482", "5prime", 18).parseLong(3)==1,
			"Movement and absolute-error changes conserve their paired N; label switches remain in the baseline member cohort");
		signFlipCase();
		largeErrorCase();
	}

	/** Moving across truth by the same distance changes signed error, not accuracy. */
	private void signFlipCase() throws Exception{
		for(final String side:new String[]{"flip_before", "flip_after"}){
			final boolean before=side.equals("flip_before");
			final int start=before ? 103 : 97, stop=before ? 196 : 202;
			final String errors=before ? "3\t-3" : "-3\t3";
			write(side+".genomes.tsv", "input_id\tcalls\tmatched\tmasked\tunmatched\nz\t1\t1\t0\t0\n");
			write(side+".references.tsv", "input_id\treference_id\tseqid\tstart\tstop\tstrand\tscored\nz\tr\ts\t100\t199\t+\ttrue\n");
			write(side+".queries.tsv", "input_id\tseqid\tstart\tstop\tstrand\temitted_member\tmodel\toutcome\tscored_eligible\tmatched_reference_id\n"
				+"z\ts\t"+start+'\t'+stop+"\t+\tRF00174\tm\tTP\ttrue\tr\n");
			write(side+".witnesses.tsv", "input_id\treference_id\temitted_member\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3\n"
				+"z\tr\tRF00174\tm\ts\t"+start+'\t'+stop+"\t100\t199\t+\t"+errors+'\n');
		}
		CompareCobalaminBenchmarkGrades.main(new String[]{"in1="+path("flip_before"), "in2="+path("flip_after"), "out="+path("flip"), "radius=6"});
		for(final String end:new String[]{"5prime", "3prime"}){
			final LineParser1 row=tableRow("flip", "paired_details", "cobalamin", end, 18);
			check(row.parseLong(3)==1 && row.parseLong(4)==0 && row.parseLong(5)==0 && row.parseLong(6)==1 && row.parseLong(7)==1
				&& row.parseDouble(10)==3 && row.parseDouble(11)==3 && row.parseLong(12)==1
				&& row.parseLong(13)==0 && row.parseLong(14)==0 && row.parseLong(15)==1,
				"A three-base sign flip moves the end while preserving absolute mean, within-three membership and unchanged accuracy");
		}
	}

	/** Reciprocal-half-compatible intervals whose squared offsets exceed int. */
	private void largeErrorCase() throws Exception{
		for(final String side:new String[]{"wide_before", "wide_after"}){
			final boolean before=side.equals("wide_before");
			final String start=before ? "100100" : "100", stop=before ? "300099" : "200099", error=before ? "100000" : "0";
			write(side+".genomes.tsv", "input_id\tcalls\tmatched\tmasked\tunmatched\nz\t1\t1\t0\t0\n");
			write(side+".references.tsv", "input_id\treference_id\tseqid\tstart\tstop\tstrand\tscored\nz\tr\ts\t100\t200099\t+\ttrue\n");
			write(side+".queries.tsv", "input_id\tseqid\tstart\tstop\tstrand\temitted_member\tmodel\toutcome\tscored_eligible\tmatched_reference_id\nz\ts\t"
				+start+'\t'+stop+"\t+\tRF00174\tm\tTP\ttrue\tr\n");
			write(side+".witnesses.tsv", "input_id\treference_id\temitted_member\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3\nz\tr\tRF00174\tm\ts\t"
				+start+'\t'+stop+"\t100\t200099\t+\t"+error+'\t'+error+'\n');
		}
		CompareCobalaminBenchmarkGrades.main(new String[]{"in1="+path("wide_before"), "in2="+path("wide_after"), "out="+path("wide"), "radius=6"});
		final LineParser1 row=pairedRow("wide", "cobalamin", "5prime");
		close(row.parseDouble(4), 100000);
		check(row.parseDouble(5)==0 && tableRow("wide", "paired_population", "cobalamin", "5prime", 5).parseLong(4)==0
			&& !read("wide.paired_near_endpoints.tsv").contains("\ncobalamin\t"),
			"Large errors must not overflow int squares or enter the small-error subset after ON corrects them");
	}
	private LineParser1 pairedRow(final String prefix, final String cohort, final String end) throws Exception{
		return tableRow(prefix, "paired_endpoints", cohort, end, 13);
	}
	private LineParser1 tableRow(final String prefix, final String suffix, final String cohort, final String end, final int width) throws Exception{
		return tableRow(prefix, suffix, cohort, end, width, "COBALAMIN");
	}
	private LineParser1 tableRow(final String prefix, final String suffix, final String cohort, final String end,
			final int width, final String kind) throws Exception{
		final String text=read(prefix+'.'+suffix+".tsv");
		final String version=suffix.equals("paired_details") ? "V1" : "V2";
		check(text.startsWith('#'+kind+'_'+suffix.toUpperCase(Locale.ROOT)+'_'+version+'\n'), "Paired tables require their exact declared version markers");
		final LineParser1 fields=new LineParser1('\t');
		for(final String row:text.split("\n")){
			if(row.startsWith("#")){continue;}
			fields.set(row.getBytes(StandardCharsets.UTF_8));
			if(fields.termEquals(cohort, 0) && fields.termEquals(end, 1)){
				check(fields.terms()==width, "Paired report width must match its declared version");
				for(int i=2; i<width-(suffix.equals("paired_population") ? 0 : 1); i++){
					if(kind.equals("NCRNA") && (suffix.equals("paired_endpoints") && i==11 || suffix.equals("paired_details") && i==16)){
						check(fields.termEquals("NA", i), "A family union cannot claim unmeasured member-label switches");
						continue;
					}
					check(Double.isFinite(fields.parseDouble(i)), "Every emitted paired statistic must be defined and finite");
				}
				return fields;
			}
		}
		throw new AssertionError("Missing paired cohort/end: "+cohort+' '+end);
	}

	/** The same literal before/after loci in another family retain gains, losses
	 * and paired-end harm, with one union cohort rather than cobalamin subtypes. */
	private void genericFamily() throws Exception{
		for(final String side:new String[]{"before", "after"}){
			for(final String suffix:new String[]{"genomes", "queries", "references", "witnesses"}){
				write("generic_"+side+'.'+suffix+".tsv", read(side+"_more."+suffix+".tsv")
					.replace("emitted_member", "family_union").replace("RF00174", "SAM_I").replace("RF01482", "SAM_I").replace("RF01689", "SAM_I"));
			}
			write("generic_"+side+".family.tsv", "#NCRNA_NORMALIZED_FAMILY_V1\nfamily\tSAM_I\n");
		}
		CompareCobalaminBenchmarkGrades.main(new String[]{"in1="+path("generic_before"), "in2="+path("generic_after"),
			"out="+path("generic"), "family=SAM_I", "radius=6"});
		final String summary=read("generic.summary.tsv");
		check(summary.startsWith("#NCRNA_COMPARISON_COUNTS_V2\n")
			&& summary.contains("genomes_with_any_family_call\tLOST\t1\t")
			&& summary.contains("genomes_with_scored_reference_detection\tLOST\t2\t")
			&& summary.contains("genomes_with_scored_reference_detection\tGAINED\t1\t"),
			"Another family must preserve distinct raw-call and scored first/last genome populations");
		final LineParser1 five=tableRow("generic", "paired_endpoints", "SAM_I", "5prime", 13, "NCRNA");
		final LineParser1 three=tableRow("generic", "paired_endpoints", "SAM_I", "3prime", 13, "NCRNA");
		check(five.parseLong(3)==5 && three.parseLong(3)==5 && three.parseLong(7)==2 && three.parseLong(9)==1,
			"The fixed paired cohort must retain both worsening ends, including the previously-near end");
		close(five.parseDouble(4), Math.sqrt(552.0/5));
		close(three.parseDouble(5), Math.sqrt(373.0/5));
		check(read("generic.paired_population.tsv").split("\n").length==4
			&& !summary.contains("RF00174") && !read("generic.paired_endpoints.tsv").contains("cobalamin")
			&& !summary.contains("changed_emitted_member") && !summary.contains("SAM_I_to_SAM_I"),
			"Generic output has exactly one family cohort and no borrowed subtype labels");
		tableRow("generic", "paired_details", "SAM_I", "5prime", 18, "NCRNA");
		for(final String side:new String[]{"empty_before", "empty_after"}){
			write(side+".genomes.tsv", "input_id\tcalls\tmatched\tmasked\tunmatched\nz\t0\t0\t0\t0\n");
			write(side+".references.tsv", "input_id\treference_id\tseqid\tstart\tstop\tstrand\tscored\n");
			write(side+".queries.tsv", "input_id\tseqid\tstart\tstop\tstrand\tfamily_union\tmodel\toutcome\tscored_eligible\tmatched_reference_id\n");
			write(side+".witnesses.tsv", "input_id\treference_id\tfamily_union\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3\n");
			write(side+".family.tsv", "#NCRNA_NORMALIZED_FAMILY_V1\nfamily\t"+(side.equals("empty_before") ? "SAM_I" : "RNaseP")+'\n');
		}
		boolean rejected=false;
		try{
			CompareCobalaminBenchmarkGrades.main(new String[]{"in1="+path("empty_before"), "in2="+path("empty_after"),
				"out="+path("wrong_family"), "family=SAM_I"});
		}catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("requested family");}
		check(rejected && !new File(path("wrong_family.summary.tsv")).exists(),
			"Empty query populations cannot hide a cross-family comparison");
		for(final boolean explicit:new boolean[]{false, true}){
			final String output=path(explicit ? "empty_as_explicit_cobalamin" : "empty_as_default_cobalamin");
			final String[] args=explicit ? new String[]{"in1="+path("empty_before"), "in2="+path("empty_before"), "out="+output, "family=cobalamin"}
				: new String[]{"in1="+path("empty_before"), "in2="+path("empty_before"), "out="+output};
			boolean wrongFamily=false;
			try{CompareCobalaminBenchmarkGrades.main(args);}
			catch(IllegalArgumentException expected){wrongFamily=expected.getMessage().contains("requested family");}
			check(wrongFamily && !new File(output+".summary.tsv").exists(),
				"Default and explicit cobalamin must honor another family's declaration even with no queries");
		}
	}

	/** No reference exists: the physical call remains, while no reference-length
	 * measurement or endpoint distribution may be invented for it.
	 */
	private void unmatchedLengths() throws Exception{
		for(final String side:new String[]{"unmatched_before", "unmatched_after"}){
			write(side+".genomes.tsv", "input_id\tcalls\tmatched\tmasked\tunmatched\nu\t1\t0\t0\t1\n");
			write(side+".references.tsv", "input_id\treference_id\tseqid\tstart\tstop\tstrand\tscored\n");
			write(side+".queries.tsv", "input_id\tseqid\tstart\tstop\tstrand\temitted_member\tmodel\toutcome\tscored_eligible\tmatched_reference_id\n"
				+"u\ts\t10\t29\t+\tRF00174\tm\tUNMATCHED_OTHER\tfalse\t.\n");
			write(side+".witnesses.tsv", "input_id\treference_id\temitted_member\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3\n");
		}
		CompareCobalaminBenchmarkGrades.main(new String[]{"in1="+path("unmatched_before"), "in2="+path("unmatched_after"),
			"out="+path("unmatched"), "radius=6"});
		final String calls=read("unmatched.call_lengths.tsv");
		check(calls.startsWith("#COBALAMIN_CALL_LENGTHS_V2\n") && calls.split("\n").length==4
			&& calls.contains("0\tu\ts\t10\t29\t+\tRF00174\tm\tUNMATCHED_OTHER\t20\t0\t")
			&& calls.contains("1\tu\ts\t10\t29\t+\tRF00174\tm\tUNMATCHED_OTHER\t20\t0\t"),
			"Every unmatched call must retain its length and explicit zero reference-pair count");
		check(read("unmatched.lengths.tsv").startsWith("#COBALAMIN_REFERENCE_LENGTHS_V2\n") && read("unmatched.lengths.tsv").split("\n").length==2,
			"No reference means no reference-length rows, not zero lengths or undefined ratios");
		check(read("unmatched.paired_endpoints.tsv").split("\n").length==2
			&& read("unmatched.paired_near_endpoints.tsv").split("\n").length==2
			&& read("unmatched.paired_details.tsv").split("\n").length==2
			&& read("unmatched.paired_population.tsv").split("\n").length==10,
			"An entirely empty paired population retains eight zero count rows and no statistics");
	}

	/** Real report-reader acceptance/rejection, using literal format fixtures;
	 * no saved biological measurements are rewritten or regenerated.
	 */
	private void reportFormats() throws Exception{
		final String tables=path("format_tables");
		final String replay=path("format_replays");
		check(new File(tables).mkdir() && new File(replay).mkdir(), "Report format fixtures need fresh roots");
		final String metrics="inputs\tscored_references\tmasked_references\tcalls\tmatched\tunmatched\tmissed\tmasked_calls\tduplicate_queries\texact5\texact3\tboth_exact\tprecision\trecall\tqualification\tcpu_user_seconds\tcpu_system_seconds\twall_seconds\talign_rf00174\talign_rf01482\talign_rf01689\n";
		final String memberHeader="member\toff_calls\ton_calls\toff_matched\ton_matched\toff_masked\ton_masked\toff_unmatched\ton_unmatched\tretained_unmatched_off\tretained_unmatched_on\ton_only_unmatched\toff_only_unmatched\tlabel_switch_in\tlabel_switch_out\tunmatched_net_change\n";
		final String spatialHeader="member\ton_unmatched\tretained_exact_unmatched\ton_only_unmatched\ton_only_no_off_overlap\ton_only_off_overlap\ton_only_shares_off_with_another_on\n";
		final String replayHeader="kind\tradius_bp\tmargin\toff_member_cohort\tend\tpaired_total\treplayable_n\tunknown_n\toff_rms_nt\treplay_rms_nt\tbetter\tworse\tunchanged\tbetter_fraction\texact\twithin3\twithin10\tmoved_from_anchor\tunmoved_from_anchor\tworse_off_within_radius\tworse_off_beyond_radius\toff_within_radius_n\toff_within_radius_rms_nt\treplay_same_subset_rms_nt\tqualification\n";
		final String legacyHeader="cohort\tend\tradius_bp\tpaired_n\toff_rms_nt\ton_rms_nt\tbetter\tworse\tunchanged\tbetter_fraction\tworse_off_within_radius\tworse_off_beyond_radius\toff_within_radius_n\toff_within_radius_rms_nt\ton_same_subset_rms_nt\tlabel_switches\tqualification\n";
		String members=memberHeader;
		String spatial=spatialHeader;
		for(final String member:new String[]{"RF00174", "RF01482", "RF01689", "ALL"}){
			members+=member+"\t1\t1\t1\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
			spatial+=member+"\t0\t0\t0\t0\t0\t0\n";
		}
		for(final String kind:new String[]{"sliding", "positional"}){
			for(final int radius:new int[]{6, 10, 16}){
				final String replayCell="format_replays/"+kind+"_r"+radius;
				check(new File(path(replayCell)).mkdir(), "Fresh replay fixture directory required");
				String replayRows=replayHeader;
				for(final String margin:new String[]{"0", "0.005", "0.02", "0.1"}){
					for(final String cohort:new String[]{"cobalamin", "RF00174", "RF01482", "RF01689"}){
						for(final String end:new String[]{"5prime", "3prime"}){
							replayRows+=kind+'\t'+radius+'\t'+margin+'\t'+cohort+'\t'+end+"\t1\t1\t0\t0\t0\t0\t0\t1\tNA\t1\t1\t1\t0\t1\t0\t0\t1\t0\t0\tHISTORICAL_DIAGNOSTIC\n";
						}
					}
				}
				write(replayCell+"/comparison.replay_endpoints.tsv", replayRows);
				for(final int panel:new int[]{1000, 300}){
					final String cell="format_tables/"+kind+"_r"+radius+"_n"+panel;
					check(new File(path(cell)).mkdir(), "Fresh caller format fixture directory required");
					for(final String side:new String[]{"off", "on"}){
						write(cell+'/'+side+".metrics.tsv", metrics+panel+"\t1\t0\t1\t1\t0\t0\t0\t0\t1\t1\t1\t1\t1\tFIXTURE\t1\t0\t1\t1\t0\t0\n");
						write(cell+'/'+side+".cobalamin.radii.tsv", "family\toverlap\tradius_bp\tmatched\treachable5\treachable3\treachable_both\n"
							+"cobalamin\treciprocal_half\t0\t1\t1\t1\t1\ncobalamin\treciprocal_half\t3\t1\t1\t1\t1\n"
							+"cobalamin\treciprocal_half\t10\t1\t1\t1\t1\ncobalamin\treciprocal_half\t100\t1\t1\t1\t1\n");
					}
					write(cell+"/changes.tsv", read("paired.summary.tsv"));
					write(cell+"/members.tsv", members);
					write(cell+"/unmatched_spatial.tsv", spatial);
					for(final String suffix:new String[]{"paired_endpoints", "paired_near_endpoints", "paired_population"}){
						write(cell+'/'+suffix+".tsv", read("paired."+suffix+".tsv").replace("\t6\t", "\t"+radius+"\t"));
					}
				}
			}
		}
		CompareCobalaminBenchmarkGrades.main(new String[]{"mode=report", "in="+tables, "out="+path("report_v2")});
		check(read("report_v2/paired_population.tsv").contains("\tRF01482\t5prime\t6\t0\t0\n"), "Actual report reader must preserve empty cohort populations");
		final String currentHtml=read("report_v2/report.html");
		check(!currentHtml.contains("Historical saved-score diagnostics") && !currentHtml.contains("HISTORICAL_DIAGNOSTIC")
			&& !currentHtml.contains("<td>NA</td>") && !currentHtml.contains("better/(better+worse)"),
			"Default reports must neither read nor display historical replay rows or undefined derived fields");
		boolean replayRejected=false;
		try{CompareCobalaminBenchmarkGrades.main(new String[]{"mode=report", "in="+tables, "in2="+replay, "out="+path("report_replay_rejected")});}
		catch(IllegalArgumentException expected){replayRejected=expected.getMessage().contains("Only explicit historical=t");}
		check(replayRejected && !new File(path("report_replay_rejected")).exists(), "Current mode must reject historical replay input before creating output");
		final String first="format_tables/sliding_r6_n1000/paired_endpoints.tsv";
		replaceFixture(first, legacyHeader, "before_old_header");
		boolean rejected=false;
		try{CompareCobalaminBenchmarkGrades.main(new String[]{"mode=report", "in="+tables, "out="+path("report_old_rejected")});}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(path("report_old_rejected")).exists(), "Default report reader must reject a version-1 paired header before output");
		for(final String kind:new String[]{"sliding", "positional"}){
			for(final int radius:new int[]{6, 10, 16}){
				String old=legacyHeader;
				for(final String cohort:new String[]{"cobalamin", "RF00174", "RF01482", "RF01689"}){
					for(final String end:new String[]{"5prime", "3prime"}){old+=cohort+'\t'+end+'\t'+radius+"\t1\t0\t0\t0\t0\t1\tNA\t0\t0\t1\t0\t0\t0\tHISTORICAL\n";}
				}
				for(final int panel:new int[]{1000, 300}){
					final String cell="format_tables/"+kind+"_r"+radius+"_n"+panel;
					replaceFixture(cell+"/paired_endpoints.tsv", old, "before_historical_mode");
					replaceFixture(cell+"/changes.tsv", "measure\tretained\tlost\tgained\tqualification\ngraded_scored_references\t1\t0\t0\tHISTORICAL\n"
						+"lost_references_without_any_candidate_overlap\tNA\t0\tNA\tHISTORICAL\n", "before_historical_mode");
				}
			}
		}
		CompareCobalaminBenchmarkGrades.main(new String[]{"mode=report", "historical=t", "in="+tables, "in2="+replay, "out="+path("report_historical")});
		final String historicalHtml=read("report_historical/report.html");
		check(historicalHtml.contains("Historical format selected explicitly") && historicalHtml.contains("Historical saved-score diagnostics")
			&& historicalHtml.contains("HISTORICAL_DIAGNOSTIC") && historicalHtml.contains("<td>NA</td>"),
			"Only explicitly selected historical mode may display the unchanged historical replay values");
	}
	private void countReportFormats() throws Exception{
		write("after.metrics.tsv", "inputs\tscored_references\tmasked_references\tcalls\tmatched\tunmatched\tmissed\tmasked_calls\tduplicate_queries\texact5\texact3\tboth_exact\tprecision\trecall\tqualification\tcpu_user_seconds\tcpu_system_seconds\twall_seconds\talign_rf00174\talign_rf01482\talign_rf01689\n"
			+"5\t6\t0\t4\t3\t1\t3\t0\t0\t3\t3\t3\t0.75\t0.5\tFIXTURE\t1\t0\t1\t1\t1\t1\n");
		write("after.members.tsv", "emitted_member\tcalls\tmatched\tmasked\tunmatched_duplicate\tunmatched_other\tqualification\n"
			+"RF00174\t2\t2\t0\t0\t0\tFIXTURE\nRF01482\t1\t1\t0\t0\t0\tFIXTURE\nRF01689\t1\t0\t0\t0\t1\tFIXTURE\n");
		for(final String cohort:new String[]{"cobalamin", "RF00174"}){
			final int n=cohort.equals("cobalamin") ? 3 : 2;
			write("after_"+cohort+".summary.tsv", "family\toverlap\tmatched\texact5\texact3\tboth_exact\n"+cohort+"\treciprocal_half\t"+n+'\t'+n+'\t'+n+'\t'+n+'\n');
			write("after_"+cohort+".radii.tsv", "family\toverlap\tradius_bp\tmatched\treachable5\treachable3\treachable_both\n"+cohort+"\treciprocal_half\t3\t"+n+'\t'+n+'\t'+n+'\t'+n+'\n');
		}
		write("count_manifest.tsv", "configuration\tpanel\tstate\tprefix\tcomparison\nmeasured\t5\tMEASURED\t"+path("after")+'\t'+path("valid")
			+"\nnot_run\t5\tNOT_RUN\t.\t.\n");
		ReportNcrnaBoundaryComparison.main(new String[]{"in="+path("count_manifest.tsv"), "out="+path("count_report"), "family=cobalamin", "mainmember=RF00174"});
		check(!read("count_report/comparison.tsv").contains("not_run") && !read("count_report/comparison.tsv").contains("\tNA")
			&& read("count_report/population.tsv").contains("not_run\t5\tcobalamin\tNOT_RUN\n"), "Unmeasured report configurations belong in state accounting, not numeric measurement rows");
		final String saved=read("valid.summary.tsv");
		replaceFixture("valid.summary.tsv", "measure\tretained\tlost\tgained\tqualification\ngraded_scored_references\t1\t3\t2\tHISTORICAL\n"
			+"lost_references_without_any_candidate_overlap\tNA\t1\tNA\tHISTORICAL\n", "before_historical_mode");
		boolean rejected=false;
		try{ReportNcrnaBoundaryComparison.main(new String[]{"in="+path("count_manifest.tsv"), "out="+path("count_old_rejected"), "family=cobalamin", "mainmember=RF00174"});}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(path("count_old_rejected")).exists(), "The second active comparison reader must reject unselected historical input");
		ReportNcrnaBoundaryComparison.main(new String[]{"historical=t", "in="+path("count_manifest.tsv"), "out="+path("count_historical"), "family=cobalamin", "mainmember=RF00174"});
		check(read("count_historical/report.html").contains("Historical input format selected explicitly"), "Explicit historical count reports must be labeled");
		replaceFixture("valid.summary.tsv", saved, "historical_mode_input");
	}
	/** Keep each deliberately different schema fixture; production writers stay fresh-only. */
	private void replaceFixture(final String name, final String value, final String previousVersion) throws Exception{
		final File previous=new File(path(name+'.'+previousVersion));
		check(!previous.exists(), "A fixture preimage must have a unique path: "+previous);
		Files.move(new File(path(name)).toPath(), previous.toPath());
		write(name, value);
	}
	private static void close(final double actual, final double expected){
		check(Math.abs(actual-expected)<1e-6, "Independent endpoint RMS/fraction oracle differs: actual="+actual+", expected="+expected);
	}
	private void write(final String name, final String value){try(Output out=new Output(path(name))){out.writer.print(value);}}
	private String read(final String name) throws Exception{return new String(Files.readAllBytes(new File(path(name)).toPath()), StandardCharsets.UTF_8);}
	private String path(final String name){return new File(root, name).getPath();}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final String root;
}
