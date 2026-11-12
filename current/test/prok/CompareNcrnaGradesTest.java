package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import parse.Parser;
import parse.PreParser;
import parse.Parse;
import prok.CompareNcrnaGrades;
import prok.GradeNcrnaReference;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Hand-counted fixed-reference comparisons, including colliding locus IDs in
 * different inputs/families and exchanged detections with unchanged total TP.
 * Test-only literal text is deliberately independent of the comparison algorithm.
 * @author Brian Bushnell, Raiden
 */
public final class CompareNcrnaGradesTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, CompareNcrnaGradesTest.class, false);
		final Parser parser=new Parser();
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(java.util.Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			check(parser.parse(arg, a, b), "Unknown test argument: "+arg);
		}
		check(parser.out1!=null, "Supply out=fresh_fixture_directory");
		try{new CompareNcrnaGradesTest(new File(parser.out1)).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	private CompareNcrnaGradesTest(final File root_){root=root_;}

	private void run(){
		check(root.mkdir(), "Fixture root must be fresh");
		write("families.tsv", GradeNcrnaReference.RULE_HEADER+"\nF\tF\tF\ttrue\nG\tG\tG\ttrue\nX\tX\tX\tfalse\n");
		baselineRaw=GradeNcrnaReference.RAW_HEADER+'\n'
			+query("A", "F", "TP", 1, "R", ".", false, "0", "0")
			+query("A", "F", "TP", 2, "L", ".", false, "0", "0")
			+query("A", "F", "MASKED", 3, ".", "M", false, ".", ".")
			+query("A", "G", "TP", 4, "R", ".", false, "2", "-1").replace("seq\t10\t29", "seq\t11\t27")
			+query("A", "F", "FP", 5, ".", ".", true, ".", ".")
			+query("A", "X", "EXCLUDED", 6, ".", ".", false, ".", ".")
			+reference("A", "F", "R", "TP", "1")
			+reference("A", "F", "L", "TP", "2")
			+reference("A", "F", "M", "DONT_CARE", ".")
			+reference("A", "G", "R", "TP", "4")
			+reference("B", "F", "R", "FN", ".")
			+reference("Z", "F", "Z", "FN", ".");
		candidateRaw=GradeNcrnaReference.RAW_HEADER+'\n'
			+reference("Z", "F", "Z", "FN", ".")
			+reference("B", "F", "R", "TP", "1")
			+query("B", "F", "TP", 1, "R", ".", false, "0", "0")
			+query("A", "F", "TP", 1, "R", ".", false, "1", "0").replace("seq\t10\t29", "seq\t11\t29")
			+query("A", "F", "MASKED", 2, ".", "M", false, ".", ".")
			+query("A", "G", "TP", 3, "R", ".", false, "0", "0")
			+reference("A", "G", "R", "TP", "3")
			+reference("A", "F", "M", "DONT_CARE", ".")
			+reference("A", "F", "L", "FN", ".")
			+reference("A", "F", "R", "TP", "1");
		baselineInputs=GradeNcrnaReference.INPUT_HEADER+'\n'+input("A", 5, 1)+input("B", 0, 0)+input("Z", 0, 0);
		candidateInputs=GradeNcrnaReference.INPUT_HEADER+'\n'+input("Z", 0, 0)+input("B", 1, 0)+input("A", 3, 0);
		baselineSummary=GradeNcrnaReference.SUMMARY_HEADER+'\n'
			+summary("F", "2\t1\t2\t4\t1\t1\t1\t0\t2\t2\t2", "0.666667\t0.500000")
			+summary("G", "1\t0\t0\t1\t0\t0\t0\t0\t0\t0\t0", "1\t1")
			+summary("X", "0\t0\t0\t0\t0\t0\t0\t1\t0\t0\t0", "NA\tNA");
		candidateSummary=GradeNcrnaReference.SUMMARY_HEADER+'\n'
			+summary("F", "2\t0\t2\t4\t1\t1\t0\t0\t1\t2\t1", "1\t0.500000")
			+summary("G", "1\t0\t0\t1\t0\t0\t0\t0\t1\t1\t1", "1\t1")
			+summary("X", "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0", "NA\tNA");
		grade("baseline", baselineRaw, baselineInputs, baselineSummary);
		grade("candidate", candidateRaw, candidateInputs, candidateSummary);
		invoke("candidate", "good", "families.tsv");
		final String[] result=lines("good.summary.tsv");
		check(result.length==4, "Every configured family remains in the comparison");
		check(result[1].equals("F\t1\t1\t1\t1\t1\t1\t1\t0\t1"),
			"One retained, lost and gained TP; one persistent zero-input FN; mask and endpoint loss remain distinct");
		check(result[2].equals("G\t1\t0\t0\t0\t0\t1\t0\t0\t0"), "A same-ID different-family locus retains its independent endpoint improvement");
		check(result[3].equals("X\t0\t0\t0\t0\t0\t0\t0\t0\t0"), "Excluded family stays separate from negatives");
		final String[] counts=lines("good.counts.tsv");
		check(Arrays.asList(counts).contains("baseline\t"+baselineSummary.split("\n")[1]), "Original duplicate/FP counts must survive");
		check(Arrays.asList(counts).contains("candidate\t"+candidateSummary.split("\n")[1]), "Original candidate mask/precision/recall counts must survive");
		final String[] loci=lines("good.loci.tsv");
		check(loci.length==7, "All six reference loci, including masks, must remain");
		check(Arrays.asList(loci).contains("A\tG\tR\tseq\t10\t29\t-\tTP\tTP\t2\t-1\t0\t0\ttrue\tfalse\tfalse"),
			"Minus-strand offsets follow transcript direction and can improve independently of the colliding plus-strand family");
		check(Arrays.asList(loci).contains("B\tF\tR\tseq\t10\t29\t+\tFN\tTP\t.\t.\t0\t0\tfalse\tfalse\tfalse"),
			"Equal contig/reference names in different genomes cannot match each other");
		negative("geometry", candidateRaw.replace(reference("A", "F", "R", "TP", "1"),
			reference("A", "F", "R", "TP", "1").replace("\t10\t29\t", "\t11\t29\t")), candidateInputs, candidateSummary, "geometry");
		negative("disposition", candidateRaw.replace(reference("Z", "F", "Z", "FN", "."),
			reference("Z", "F", "Z", "DONT_CARE", ".")), candidateInputs, candidateSummary, "disposition");
		negative("missing_reference", candidateRaw.replace(reference("Z", "F", "Z", "FN", "."), ""), candidateInputs, candidateSummary, "Missing reference");
		negative("missing_query", candidateRaw.replace(query("B", "F", "TP", 1, "R", ".", false, "0", "0"), ""), candidateInputs, candidateSummary, "query counts");
		negative("different_fasta", candidateRaw, candidateInputs.replace("/A.fa", "/other.fa"), candidateSummary, "FASTA");
		negative("missing_zero_input", candidateRaw, candidateInputs.replace(input("Z", 0, 0), ""), candidateSummary, "populations");
		negative("duplicate_input", candidateRaw, candidateInputs+input("A", 3, 0), candidateSummary, "duplicated");
		negative("changed_overlap", candidateRaw, candidateInputs, candidateSummary.replace("0.500000\tDECLARED", "0.750000\tDECLARED"), "Overlap settings");
		negative("truncated_summary", candidateRaw, candidateInputs, GradeNcrnaReference.SUMMARY_HEADER+'\n', "populations");
		negative("duplicate_tp", candidateRaw+query("A", "F", "TP", 99, "R", ".", false, "0", "0"), candidateInputs, candidateSummary, "Two query TPs");
		final String bQuery=query("B", "F", "TP", 1, "R", ".", false, "0", "0");
		negative("query_contig", candidateRaw.replace(bQuery, bQuery.replace("\tseq\t", "\tother\t")), candidateInputs, candidateSummary, "contig or strand");
		final String gQuery=query("A", "G", "TP", 3, "R", ".", false, "0", "0");
		negative("query_strand", candidateRaw.replace(gQuery, gQuery.replace("\t-\t", "\t+\t")), candidateInputs, candidateSummary, "contig or strand");
		final String shifted=query("A", "F", "TP", 1, "R", ".", false, "1", "0").replace("seq\t10\t29", "seq\t11\t29");
		negative("query_offset", candidateRaw.replace(shifted, shifted.replace("\t1\t0\n", "\t2\t0\n")), candidateInputs, candidateSummary, "signed offsets");
		negative("query_position", candidateRaw.replace(shifted, shifted.replace("seq\t11\t29", "seq\t12\t29")), candidateInputs, candidateSummary, "signed offsets");
		write("other_families.tsv", GradeNcrnaReference.RULE_HEADER+"\nF\tF\tF\tfalse\n");
		expect(()->invoke("candidate", "different_rules", "other_families.tsv"), "family maps differ");
		check(!new File(path("different_rules.summary.tsv")).exists(), "Changed family rules cannot publish a comparison");
		System.err.println("COMPARE_NCRNA_GRADES_TEST_PASS");
	}

	private void negative(final String name, final String raw, final String inputs, final String summary, final String reason){
		grade(name, raw, inputs, summary);
		expect(()->invoke(name, name+"_comparison", "families.tsv"), reason);
		for(final String suffix:new String[]{".summary.tsv", ".counts.tsv", ".loci.tsv"}){
			check(!new File(path(name+"_comparison"+suffix)).exists(), "Rejected grade pairs must fail before opening output");
		}
	}
	private void invoke(final String candidate, final String output, final String rules){
		CompareNcrnaGrades.main(new String[]{"in1="+path("baseline"), "in2="+path(candidate), "out="+path(output),
			"familymap1="+path("families.tsv"), "familymap2="+path(rules), "t=1", "pigz=f", "unpigz=f"});
	}
	private void grade(final String name, final String raw, final String inputs, final String summary){
		write(name+".raw.tsv", raw);
		write(name+".inputs.tsv", inputs);
		write(name+".summary.tsv", summary);
	}
	private String input(final String id, final int included, final int excluded){
		return id+'\t'+path(id+".fa")+'\t'+path(id+".calls.gff")+'\t'+path(id+".ref.gff")+'\t'
			+(included+excluded)+"\t0\t0\t0\t"+included+'\t'+excluded+'\n';
	}
	private static String summary(final String family, final String counts, final String ratios){
		return family+"\t3\t"+counts+'\t'+ratios+"\t0.500000\tDECLARED_INPUT_BINDINGS_EXPLICIT_REFERENCE_STATUS\n";
	}
	private static String reference(final String input, final String family, final String id, final String outcome, final String matched){
		final String span=id.equals("L") ? "50\t69" : id.equals("M") ? "90\t109" : "10\t29";
		final char strand=family.equals("G") ? '-' : '+';
		return input+"\tREFERENCE\t"+family+"\tseq\t"+span+'\t'+strand+'\t'+outcome+'\t'+id+"\t.\t.\t"+matched+"\t.\t.\t.\tmodel\t.\t.\n";
	}
	private static String query(final String input, final String family, final String outcome, final int line,
		final String matched, final String mask, final boolean duplicate, final String d5, final String d3){
		final String span=matched.equals("L") ? "50\t69" : mask.equals("M") ? "90\t109" : "10\t29";
		final boolean eligible=outcome.equals("TP") || duplicate;
		final char strand=family.equals("G") ? '-' : '+';
		return input+"\tQUERY\t"+family+"\tseq\t"+span+'\t'+strand+'\t'+outcome+"\t.\t"+line+'\t'+matched+"\t.\t"+mask+'\t'
			+(outcome.equals("EXCLUDED") ? ".\t." : duplicate+"\t"+eligible)+"\t"+family+"_model\t"+d5+'\t'+d3+'\n';
	}
	private void write(final String name, final String text){
		try(Output out=new Output(path(name))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
	}
	private String[] lines(final String name){
		final ArrayList<String> rows=new ArrayList<String>();
		try(Input in=new Input(path(name))){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				rows.add(new String(line, StandardCharsets.UTF_8));
			}
		}
		return rows.toArray(new String[rows.size()]);
	}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private static void expect(final Runnable action, final String message){
		try{action.run();}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(message), "Unexpected failure: "+failure);
			return;
		}
		throw new AssertionError("Expected failure containing "+message);
	}
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private final File root;
	private String baselineRaw, candidateRaw, baselineInputs, candidateInputs, baselineSummary, candidateSummary;
}
