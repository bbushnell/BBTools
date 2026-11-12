package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.GradeNcrnaReference;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.SummarizeNcrnaCallOverlap;
import shared.Shared;
import structures.ByteBuilder;

/** Literal saved-grade geometry fixtures; no caller or production-policy changes.
 * Exact half boundaries, duplicates, third candidates, masks, strand separation,
 * repeated contigs across inputs, UTF-8 and shuffled row order have independent
 * expected outputs. Test files belong to the requested fresh directory.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaCallOverlapTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SummarizeNcrnaCallOverlapTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(java.util.Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				require(parser.parse(arg, a, b), "Unknown fixture parameter: "+arg);
			}
			require(parser.out1!=null, "Supply out=/absolute/fresh_directory");
			final File root=new File(parser.out1);
			require(root.isAbsolute() && root.mkdir(), "Fixture output directory must be fresh and absolute");
			run(root);
			pp.outstream.println("NCRNA_CALL_OVERLAP_LITERAL_FIXTURES_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final File root){
		final String[] rows={
			query("1", 100, 199, "+", "TP", "1", false),
			query("1", 200, 279, "+", "FP", "2", false),
			query("1", 100, 199, "+", "FP", "3", true),
			query("1", 300, 379, "+", "FP", "4", false),
			query("1", 350, 449, "+", "FP", "5", false),
			query("1", 650, 749, "+", "FP", "6", false),
			query("1", 50, 149, "+", "FP", "7", false),
			query("1", 40, 139, "+", "FP", "8", false),
			query("1", 200, 279, "-", "FP", "9", false),
			query("1", 110, 209, "+", "FP", "10", false),
			query("1", 610, 619, "+", "MASKED", "11", false),
			reference(100, 299, "TP", "ref%25é"),
			reference(900, 999, "FN", "missed"),
			reference(600, 699, "DONT_CARE", "mask")
		};
		for(int pass=0; pass<2; pass++){
			final String raw=new File(root, "input"+pass+".tsv").getPath();
			try(Output out=new Output(raw)){
				out.writer.println(GradeNcrnaReference.RAW_HEADER);
				for(int i=0; i<rows.length; i++){write(out, rows[pass==0 ? i : rows.length-1-i]);}
				write(out, query("2", 200, 279, "+", "FP", "1", false));
			}
			SummarizeNcrnaCallOverlap.main(new String[]{"in="+raw, "family=glycine", "out="+new File(root, "out"+pass)});
		}
		final String summary=read(new File(root, "out0.summary.tsv"));
		require(summary.endsWith("glycine\t2\t1\t10\t1\t1\t1\t1\t4\t1\t5\t1\t0\tORIGINAL_GRADES_UNCHANGED_IDENTITY_UNAVAILABLE\n"),
			"Literal original counts and exhaustive FP geometry partition disagree: "+summary);
		final String calls=read(new File(root, "out0.calls.tsv"));
		require(calls.contains("inputé\tglycine\tcontigé%25\t200\t279\t+\t2\tglycine_modelé\tfalse\tAT_LEAST_HALF_CALL_OVERLAP\tref%25é\t80\t1\tTP\t5\t.\t0\t1\t0\t100\t199\tglycine_modelé\tNA\tGEOMETRY_ONLY_NOT_RECLASSIFIED\n"),
			"Fully covered second aptamer must retain UTF-8, zero gap, original model and unmodified outcome");
		require(calls.contains("\t7\tglycine_modelé\tfalse\tAT_LEAST_HALF_CALL_OVERLAP\tref%25é\t50\t0.50000000\t"),
			"Exactly half call overlap belongs to the descriptive >= half bin");
		require(calls.contains("\t8\tglycine_modelé\tfalse\tPOSITIVE_OVERLAP_BELOW_HALF\tref%25é\t40\t0.40000000\t"),
			"Partial overlap below half must not become zero or half-covered");
		require(calls.contains("\t6\tglycine_modelé\tfalse\tNO_SCORED_REFERENCE_OVERLAP\t.\t0\t0\t.\t0\tmask\t50\t1\t450\t"),
			"Mask witness and 450 intervening bases must remain separate from scored overlap");
		final String distances=read(new File(root, "out0.distances.tsv"));
		require(distances.contains("glycine\t100\t1\t") && distances.contains("glycine\t200\t2\t")
			&& distances.contains("glycine\t500\t3\t"), "Distance counts must be cumulative, zero-reference-overlap and same-strand only");
		for(final String suffix:new String[]{"calls", "summary", "distances"}){
			require(read(new File(root, "out0."+suffix+".tsv")).equals(read(new File(root, "out1."+suffix+".tsv"))),
				"Shuffling rows within an input changed the descriptive "+suffix+" report");
		}
		final File unicode=new File(root, "unicode-family.tsv");
		try(Output out=new Output(unicode.getPath())){
			write(out, read(new File(root, "input0.tsv")).replace("\tglycine\t", "\tglyciné\t"));
		}
		SummarizeNcrnaCallOverlap.main(new String[]{"in="+unicode, "family=glyciné", "out="+new File(root, "unicode")});
		require(read(new File(root, "unicode.summary.tsv")).equals(summary.replace("glycine\t", "glyciné\t")),
			"UTF-8 family names must select the same literal rows without changing counts");
		final File malformed=new File(root, "noncontiguous.tsv");
		try(Output out=new Output(malformed.getPath())){
			out.writer.println(GradeNcrnaReference.RAW_HEADER);
			write(out, query("1", 10, 20, "+", "FP", "1", false));
			write(out, query("2", 10, 20, "+", "FP", "1", false));
			write(out, query("1", 30, 40, "+", "FP", "2", false));
		}
		boolean failed=false;
		try{
			SummarizeNcrnaCallOverlap.main(new String[]{"in="+malformed, "family=glycine", "out="+new File(root, "bad")});
		}catch(IllegalArgumentException expected){failed=expected.getMessage().contains("contiguous");}
		require(failed && !new File(root, "bad.summary.tsv").exists(), "Noncontiguous input groups must fail before a success summary");
	}

	/** Literal input rows deliberately avoid calling the production geometry code. */
	private static String query(final String input, final int start, final int stop, final String strand,
		final String outcome, final String id, final boolean duplicate){
		return (input.equals("1") ? "inputé" : input)+"\tQUERY\tglycine\tcontigé%25\t"+start+'\t'+stop+'\t'+strand+'\t'+outcome
			+"\t.\t"+id+"\t.\t.\t.\t"+duplicate+"\tfalse\tglycine_modelé\t.\t.\n";
	}
	private static String reference(final int start, final int stop, final String outcome, final String id){
		return "inputé\tREFERENCE\tglycine\tcontigé%25\t"+start+'\t'+stop+"\t+\t"+outcome+'\t'+id
			+"\t.\t.\t.\t.\t.\t.\tRF00504\t.\t.\n";
	}
	private static void write(final Output out, final String row){out.writer.print(row.getBytes(StandardCharsets.UTF_8));}
	private static String read(final File file){
		final ByteBuilder buffer=new ByteBuilder();
		try(Input in=new Input(file.getPath())){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){buffer.append(line).nl();}
		}
		return new String(buffer.toBytes(), StandardCharsets.UTF_8);
	}
	private static void require(final boolean valid, final String why){
		if(!valid){throw new AssertionError(why);}
	}
}
