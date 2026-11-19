package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.GradeNcrnaReference;
import prok.NcrnaTextIO.Output;
import prok.SummarizeNcrnaQueryOutcomes;
import shared.Shared;

/** Literal post-match counting oracle. Duplicate flags deliberately occur in
 * MASKED and FP rows; an unmatched scored-eligible query remains unmatched.
 * A zero-query prefix, excluded family and escaped percent survive unchanged.
 * No sequence, GFF or reference payload exists in this fixture.
 * @author Raiden
 */
public final class SummarizeNcrnaQueryOutcomesTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, SummarizeNcrnaQueryOutcomesTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown outcome-fixture option: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && parser.in1==null, "Supply out=fresh_fixture_directory");
			check(new File(parser.out1).mkdir(), "Fixture output must be fresh");
			new SummarizeNcrnaQueryOutcomesTest(parser.out1).run();
			pp.outstream.println("NCRNA_QUERY_OUTCOMES_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private SummarizeNcrnaQueryOutcomesTest(final String directory){root=new File(directory).getAbsolutePath();}

	private void run() throws Exception{
		final String map=GradeNcrnaReference.RULE_HEADER+"\ncobalamin\talpha\tA\ttrue\n"
			+"cobalamin\tbeta%\tB\tt\ncobalamin\tunused\tC\ttrue\nFMN\tFMN\tFMN\tf\n";
		write("families.tsv", map);
		final String raw=GradeNcrnaReference.RAW_HEADER+'\n'
			+query(2, "cobalamin", "TP", "alpha_0", "r1", ".", "false", "true")
			+query(3, "cobalamin", "MASKED", "alpha_0", ".", "mask1", "true", "false")
			+query(4, "cobalamin", "FP", "alpha_1", ".", ".", "true", "true")
			+query(5, "cobalamin", "FP", "alpha_1", ".", ".", "false", "false")
			+query(6, "cobalamin", "TP", "beta%25_0", "r2", ".", "false", "true")
			+query(7, "cobalamin", "FP", "beta%25_0", ".", ".", "false", "true")
			+query(8, "FMN", "EXCLUDED", "FMN_0", ".", ".", ".", ".")
			+"global\tREFERENCE\tcobalamin\tsame\t2\t9\t+\tTP\tr1\t.\t.\t2\t.\t.\t.\t.\t.\t.\n";
		final String summary=GradeNcrnaReference.SUMMARY_HEADER+"\ncobalamin\t1\t2\t3\t0\t2\t1\t1\t2\t0\t0\t0\t0\t.4\t1\t.5\tfixture\n"
			+"FMN\t1\t0\t0\t0\t0\t0\t0\t0\t1\t0\t0\t0\tNA\tNA\t.5\tfixture\n";
		write("valid.raw.tsv", raw);
		write("valid.summary.tsv", summary);
		invoke("valid", "families.tsv", "out.tsv");
		final String qualification="\tEMITTED_MODEL_ATTRIBUTION_AFTER_JOINT_MATCH_NOT_SUBTYPE_ACCURACY\n";
		check(read("out.tsv").equals(SummarizeNcrnaQueryOutcomes.HEADER+'\n'
			+"cobalamin\talpha\tA\ttrue\t4\t1\t1\t1\t1\t0\t2\t2"+qualification
			+"cobalamin\tbeta%25\tB\ttrue\t2\t1\t0\t0\t1\t0\t0\t2"+qualification
			+"cobalamin\tunused\tC\ttrue\t0\t0\t0\t0\t0\t0\t0\t0"+qualification),
			"Exact prefix buckets, independent flags, excluded/empty rules and escaped text must survive");
		write("bad_summary.raw.tsv", raw);
		write("bad_summary.summary.tsv", summary.replace("cobalamin\t1\t2\t3", "cobalamin\t1\t2\t4"));
		reject("bad_summary", "families.tsv", "disagree");
		write("bad_mask.raw.tsv", raw.replace("mask1\ttrue\tfalse", "mask1\ttrue\ttrue"));
		write("bad_mask.summary.tsv", summary);
		reject("bad_mask", "families.tsv", "Scored-first");
		write("unknown.raw.tsv", raw.replace("alpha_0", "unknown_0"));
		write("unknown.summary.tsv", summary);
		reject("unknown", "families.tsv", "absent");
		write("unknown_outcome.raw.tsv", raw.replace("\tTP\t", "\tUNKNOWN\t"));
		write("unknown_outcome.summary.tsv", summary);
		reject("unknown_outcome", "families.tsv", "requires FP");
		write("overlap.tsv", map+"cobalamin\talpha_\tambiguous\ttrue\n");
		write("ambiguous.raw.tsv", raw.replace("alpha_0", "alpha__0"));
		write("ambiguous.summary.tsv", summary);
		reject("ambiguous", "overlap.tsv", "multiple");
		final String validTp=query(6, "cobalamin", "TP", "beta%25_0", "r2", ".", "false", "true");
		write("duplicate_tp.raw.tsv", raw.replace(validTp, query(6, "cobalamin", "TP", "beta%25_0", "r2", ".", "true", "true")));
		write("duplicate_tp.summary.tsv", summary);
		reject("duplicate_tp", "families.tsv", "Recorded TP");
		write("masked_tp.raw.tsv", raw.replace(validTp, query(6, "cobalamin", "TP", "beta%25_0", "r2", "mask1", "false", "true")));
		write("masked_tp.summary.tsv", summary);
		reject("masked_tp", "families.tsv", "Recorded TP");
		write("boundaries.tsv", GradeNcrnaReference.RULE_HEADER+"\ncobalamin\tA\tA\ttrue\ncobalamin\tAB\tAB\ttrue\n");
		write("boundaries.raw.tsv", GradeNcrnaReference.RAW_HEADER+'\n'
			+query(2, "cobalamin", "FP", "A_0", ".", ".", "false", "false")
			+query(3, "cobalamin", "FP", "AB_0", ".", ".", "false", "false")
			+query(4, "cobalamin", "FP", "A%09tail", ".", ".", "false", "false")
			+query(5, "cobalamin", "FP", "A tail", ".", ".", "false", "false")
			+query(6, "cobalamin", "FP", "A", ".", ".", "false", "false"));
		write("boundaries.summary.tsv", GradeNcrnaReference.SUMMARY_HEADER+"\ncobalamin\t1\t0\t5\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tNA\t.5\tfixture\n");
		invoke("boundaries", "boundaries.tsv", "boundaries.out.tsv", "4,1");
		check(read("boundaries.out.tsv").equals(SummarizeNcrnaQueryOutcomes.HEADER+'\n'
			+"cobalamin\tA\tA\ttrue\t4\t0\t0\t0\t4\t0\t0\t0"+qualification
			+"cobalamin\tAB\tAB\ttrue\t1\t0\t0\t0\t1\t0\t0\t0"+qualification),
			"Whole-prefix boundaries must distinguish A from AB without false ambiguity");
		boolean countFailed=false;
		try{
			SummarizeNcrnaQueryOutcomes.main(new String[]{"in="+path("valid"), "familymap="+path("families.tsv"),
				"family=cobalamin", "expected=4,3,0", "out="+path("bad_count.tsv")});
		}catch(IllegalArgumentException expected){countFailed=expected.getMessage().contains("Included-query count differs");}
		check(countFailed && !new File(path("bad_count.tsv")).exists(), "Per-prefix retained raw counts need an independent oracle");
	}

	/** Eighteen literal grader columns; model strings are already report-escaped. */
	private static String query(final int line, final String family, final String outcome, final String model,
		final String matched, final String mask, final String duplicate, final String eligible){
		return "global\tQUERY\t"+family+"\tsame\t2\t9\t+\t"+outcome+"\t.\t"+line+'\t'+matched
			+"\t.\t"+mask+'\t'+duplicate+'\t'+eligible+'\t'+model+"\t.\t.\n";
	}
	private void invoke(final String prefix, final String map, final String output){
		invoke(prefix, map, output, "4,2,0");
	}
	private void invoke(final String prefix, final String map, final String output, final String expected){
		SummarizeNcrnaQueryOutcomes.main(new String[]{"in="+path(prefix), "familymap="+path(map), "family=cobalamin", "expected="+expected, "out="+path(output)});
	}
	private void reject(final String prefix, final String map, final String reason){
		final String output=prefix+'_'+map+".rejected.tsv";
		boolean failed=false;
		try{invoke(prefix, map, output);}
		catch(IllegalArgumentException expected){failed=expected.getMessage().contains(reason);}
		check(failed && !new File(path(output)).exists(), "Invalid reports must fail before opening output: "+prefix);
	}
	private void write(final String name, final String text){try(Output out=new Output(path(name))){out.writer.print(text);}}
	private String read(final String name) throws Exception{return new String(Files.readAllBytes(new File(path(name)).toPath()), StandardCharsets.UTF_8);}
	private String path(final String name){return new File(root, name).getPath();}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final String root;
}
