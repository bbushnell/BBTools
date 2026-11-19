package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CobalaminSavedCallAudit;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Literal final-call correspondence oracle, not an implementation of the join.
 * Covers every unmatched category, both directions of label/grade changes,
 * repeated contig names across sources, opposite strands, one-base intersection,
 * duplicate rejection and empty-call genomes. No biological sequence is read.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminSavedCallAuditTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, CobalaminSavedCallAuditTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown fixture argument: "+arg);
			}
			check(parser.out1!=null && new File(parser.out1).mkdir(), "Use one fresh out=fixture directory");
			new CobalaminSavedCallAuditTest(new File(parser.out1).getAbsolutePath()).run();
			pp.outstream.println("COBALAMIN_SAVED_CALL_AUDIT_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private CobalaminSavedCallAuditTest(final String root_){root=root_;}

	private void run() throws Exception{
		final ArrayList<String> before=new ArrayList<String>(), after=new ArrayList<String>();
		before.add(call("g0", 100, 199, '+', "RF00174", "TP", "old1"));
		before.add(call("g0", 300, 399, '+', "RF01482", "MASKED", "old2"));
		before.add(call("g0", 500, 599, '+', "RF01689", "UNMATCHED_OTHER", "old3"));
		before.add(call("g0", 700, 799, '+', "RF00174", "TP", "old4"));
		before.add(call("g0", 750, 849, '+', "RF01482", "MASKED", "old5"));
		before.add(call("g0", 900, 999, '+', "RF00174", "UNMATCHED_OTHER", "old6"));
		before.add(call("g0", 1100, 1199, '+', "RF01482", "UNMATCHED_OTHER", "old7"));
		before.add(call("g0", 1300, 1399, '+', "RF01689", "UNMATCHED_OTHER", "old8"));
		before.add(call("g0", 1800, 1899, '+', "RF01689", "MASKED", "old9"));
		before.add(call("g0", 2100, 2199, '+', "RF00174", "TP", "old10"));
		before.add(call("g1", 1600, 1699, '+', "RF00174", "TP", "other_source"));
		after.add(call("g0", 105, 204, '+', "RF00174", "UNMATCHED_OTHER", "one_matched"));
		after.add(call("g0", 305, 404, '+', "RF01482", "UNMATCHED_OTHER", "one_masked"));
		after.add(call("g0", 505, 604, '+', "RF01689", "UNMATCHED_OTHER", "shift_unmatched"));
		after.add(call("g0", 760, 790, '+', "RF00174", "UNMATCHED_OTHER", "mixed"));
		after.add(call("g0", 900, 999, '+', "RF01689", "UNMATCHED_OTHER", "new_label_model"));
		after.add(call("g0", 1100, 1199, '+', "RF01482", "TP", "now_matched"));
		after.add(call("g0", 1600, 1699, '+', "RF00174", "UNMATCHED_OTHER", "source_isolation"));
		after.add(call("g0", 100, 199, '-', "RF00174", "UNMATCHED_OTHER", "strand_isolation"));
		after.add(call("g0", 1800, 1899, '+', "RF01689", "UNMATCHED_OTHER", "exact_grade_change"));
		after.add(call("g0", 2100, 2140, '+', "RF00174", "UNMATCHED_OTHER", "split1"));
		after.add(call("g0", 2160, 2199, '+', "RF00174", "UNMATCHED_OTHER", "split2"));
		after.add(call("g1", 2400, 2499, '+', "RF00174", "MASKED", "mask_after"));
		final String beforeSources=GENOMES+"g0\t10\t3\t3\t4\ng1\t1\t1\t0\t0\n";
		final String afterSources=GENOMES+"g0\t11\t1\t0\t10\ng1\t1\t0\t1\t0\n";
		input("before", beforeSources, before, false);
		input("after", afterSources, after, false);
		invoke("before", "after", "valid", 2);
		final String members=read("valid/members.tsv"), categories=read("valid/categories.tsv");
		check(read("valid/unmatched_spatial.tsv").contains("ALL\t10\t1\t9\t2\t7\t2\n")
			&& read("valid/unmatched_spatial.tsv").contains("RF00174\t6\t0\t6\t2\t4\t2\n"),
			"Two literal split calls share one OFF span; a single ON overlapping two OFF spans is a different geometry");
		check(members.contains("ALL\t11\t12\t4\t1\t3\t1\t4\t10\t1\t1\t9\t3\t1\t1\t6\n"),
			"Four OFF unmatched become ten ON: one exact retained, nine additions, three departures, net six");
		check(members.contains("RF00174\t5\t7\t4\t0\t0\t1\t1\t6\t1\t0\t6\t0\t0\t1\t5\n")
			&& members.contains("RF01482\t3\t2\t0\t1\t2\t0\t1\t1\t0\t0\t1\t1\t0\t0\t0\n")
			&& members.contains("RF01689\t3\t3\t0\t0\t1\t0\t2\t3\t0\t1\t2\t2\t1\t0\t1\n"),
			"The exact retained call changes label, so member net counts require one outgoing/incoming switch");
		final String[] names={"NOT_UNMATCHED", "RETAINED_EXACT_UNMATCHED", "A_NO_OFF_CALL_AT_LOCUS",
			"B_UNIQUE_OFF_MASKED_OVERLAP", "C_UNIQUE_OFF_MATCHED_OVERLAP", "D_PRIOR_UNMATCHED_SHIFTED",
			"D_MULTIPLE_OR_MIXED_COUNTERPARTS", "D_EXACT_COORDINATES_CHANGED_GRADE"};
		final int[] counts={2, 1, 2, 1, 1, 1, 3, 1};
		for(int i=0; i<counts.length; i++){
			check(categories.contains("ALL\t"+names[i]+'\t'+counts[i]+'\n'), "Independent category count differs: "+names[i]);
		}
		input("reversed", afterSources, after, true);
		invoke("before", "reversed", "reordered", 2);
		check(members.equals(read("reordered/members.tsv")) && categories.equals(read("reordered/categories.tsv")),
			"Final-call order and line IDs cannot change coordinate identity or category totals");
		for(final String prefix:new String[]{"before", "after"}){
			write(prefix+"_unknown.genomes.tsv", read(prefix+".genomes.tsv"));
			write(prefix+"_unknown.queries.tsv", read(prefix+".queries.tsv").replace("\ttrue\t", "\tNA\t").replace("\tfalse\t", "\tNA\t"));
		}
		invoke("before_unknown", "after_unknown", "unknown", 2);
		for(final String table:new String[]{"members.tsv", "categories.tsv", "overlaps.tsv", "complete.tsv", "unmatched_spatial.tsv"}){
			check(read("valid/"+table).equals(read("unknown/"+table)), "Unknown eligibility cannot change saved-grade or spatial accounting: "+table);
		}
		final String unknownCalls=read("unknown/on_calls.tsv");
		check(unknownCalls.contains("\tNA\treference_") && !unknownCalls.contains("\ttrue\t") && !unknownCalls.contains("\tfalse\t"),
			"Legacy unknown eligibility must stay explicitly NA, including already graded matches");
		write("malformed.genomes.tsv", afterSources);
		write("malformed.queries.tsv", read("after.queries.tsv").replace("\ttrue\t", "\tmaybe\t"));
		reject("before", "malformed", "malformed_result", "Invalid scored eligibility");
		final ArrayList<String> duplicate=new ArrayList<String>(after);
		duplicate.add(after.get(0));
		input("duplicate", GENOMES+"g0\t12\t1\t0\t11\ng1\t1\t0\t1\t0\n", duplicate, false);
		reject("before", "duplicate", "bad_duplicate", "Duplicate coordinates");
		input("bad_counts", GENOMES+"g0\t10\t1\t0\t9\ng1\t1\t0\t1\t0\n", after, false);
		reject("before", "bad_counts", "bad_count_result", "do not conserve");
		input("bad_sources", afterSources.replace("g1\t", "g2\t"), after, false);
		reject("before", "bad_sources", "bad_source_result", "ON input");
		input("empty0", GENOMES+"g0\t0\t0\t0\t0\ng1\t0\t0\t0\t0\n", new ArrayList<String>(), false);
		input("empty1", GENOMES+"g0\t0\t0\t0\t0\ng1\t0\t0\t0\t0\n", new ArrayList<String>(), false);
		invoke("empty0", "empty1", "empty_result", 2);
		check(read("empty_result/members.tsv").contains("ALL\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"),
			"Zero-call genomes remain in the complete roster with conserved zero totals");
		final ArrayList<String> edgeBefore=new ArrayList<String>(), edgeAfter=new ArrayList<String>();
		edgeBefore.add(call("edge", Integer.MAX_VALUE-9, Integer.MAX_VALUE, '+', "RF00174", "TP", "edge_before"));
		edgeAfter.add(call("edge", Integer.MAX_VALUE-19, Integer.MAX_VALUE-9, '+', "RF00174", "UNMATCHED_OTHER", "edge_after"));
		input("edge0", GENOMES+"edge\t1\t1\t0\t0\n", edgeBefore, false);
		input("edge1", GENOMES+"edge\t1\t0\t0\t1\n", edgeAfter, false);
		invoke("edge0", "edge1", "edge_result", 1);
		check(read("edge_result/categories.tsv").contains("ALL\tC_UNIQUE_OFF_MATCHED_OVERLAP\t1\n"),
			"One inclusive shared base near the maximum coordinate is a same-locus overlap, without integer overflow");
	}

	private void input(final String prefix, final String sources, final ArrayList<String> calls, final boolean reverse){
		write(prefix+".genomes.tsv", sources);
		final ByteBuilder b=new ByteBuilder(QUERIES);
		for(int i=0; i<calls.size(); i++){b.append(calls.get(reverse ? calls.size()-1-i : i));}
		write(prefix+".queries.tsv", b.toString());
	}
	private static String call(final String source, final int start, final int stop, final char strand,
			final String member, final String outcome, final String model){
		return new ByteBuilder().append(source).append("\tshared\t").append(start).tab().append(stop).tab().append(strand)
			.tab().append(member).tab().append(model).tab().append(outcome).tab().append(outcome.equals("MASKED") ? "false" : "true")
			.tab().append(outcome.equals("TP") ? "reference_"+start : "NA").nl().toString();
	}
	private void invoke(final String before, final String after, final String output, final int inputs){
		CobalaminSavedCallAudit.main(new String[]{"in1="+path(before), "in2="+path(after),
			"out="+path(output), "expectedinputs="+inputs, "t=1"});
	}
	private void reject(final String before, final String after, final String output, final String message){
		boolean rejected=false;
		try{invoke(before, after, output, 2);}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains(message);}
		check(rejected && !new File(path(output+"/complete.tsv")).exists(), "Invalid saved data must fail before an acceptance marker: "+message);
	}
	private void write(final String file, final String text){try(Output output=new Output(path(file))){output.writer.print(text);}}
	private String read(final String file) throws Exception{return new String(Files.readAllBytes(new File(path(file)).toPath()), StandardCharsets.UTF_8);}
	private String path(final String file){return new File(root, file).getPath();}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final String root;
	private static final String GENOMES="input_id\tcalls\tmatched\tmasked\tunmatched\n";
	private static final String QUERIES="input_id\tseqid\tstart\tstop\tstrand\temitted_member\tmodel\toutcome\tscored_eligible\tmatched_reference_id\n";
}
