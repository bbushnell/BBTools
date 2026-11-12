package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import prok.BuildNcrnaReference;
import prok.CompareNcrnaScenarios;
import prok.GradeNcrnaReference;
import prok.WriteNcrnaReference;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Literal synthetic changes to assumed references; no biological policy selected.
 * @author Brian Bushnell, Raiden
 */
final class CompareNcrnaScenariosTest{
	CompareNcrnaScenariosTest(final File root_){root=root_;}

	void run(){
		check(root.mkdir(), "Scenario fixture root must be new");
		final String a=path("A%raw.fa"), b=path("B.fa");
		final String q=write("a.gff", call(10, 29, "F_1", "F")+call(10, 29, "F_duplicate", "F")
			+call(50, 69, "F_new", "F")+call(90, 109, "F_family", "F")+call(130, 149, "F_move", "F")
			+call(170, 189, "F_unknown", "F")+call(210, 229, "X_excluded", "X")
			+"chr\tBBTools\tCDS\t1\t6\t.\t+\t0\t.\nchr\tBBTools\ttRNA\t1\t6\t.\t+\t.\t.\nchr\tBBTools\trRNA\t1\t6\t.\t+\t.\t.\n");
		calls=write("calls.tsv", GradeNcrnaReference.MANIFEST_HEADER+"\nA\t"+a+'\t'+q+"\nB\t"+b+'\t'+write("b.gff", "")+'\n');
		contigs=write("contigs.tsv", WriteNcrnaReference.CONTIG_HEADER+'\n'+fields("A", a, "chr", "300")+fields("B", b, "chr", "300"));
		families=write("families.tsv", GradeNcrnaReference.RULE_HEADER+"\nF\tF\tF\ttrue\nG\tG\tG\ttrue\nX\tX\tX\tfalse\n");
		final String mask=row("mask", "A", 10, 29, "F", "scored", true);
		final String otherMask=row("other_mask", "A", 8, 31, "F", "dont_care", true);
		final String family=row("family", "A", 90, 109, "F", "scored", true);
		final String move=row("move", "A", 130, 149, "F", "scored", true);
		final String missing=row("missing", "A", 250, 269, "F", "scored", true);
		final String unknown=row("unknown", "A", 170, 189, "F", "unresolved", false);
		final String context=row("context", "A", 210, 229, "F", "context", true);
		final String empty=row("empty", "B", 250, 269, "F", "scored", true);
		final String beforeText=CompareNcrnaScenarios.HEADER+'\n'+mask+otherMask+family+move+missing+unknown+context+empty;
		final String afterText=CompareNcrnaScenarios.HEADER+'\n'+row("mask", "A", 10, 29, "F", "dont_care", true)+otherMask
			+row("family", "A", 90, 109, "G", "scored", true)+row("move", "A", 131, 149, "F", "scored", true)
			+row("new", "A", 50, 69, "F", "scored", true)+row("unknown", "A", 170, 189, "F", "removed", true)+context+empty;
		before=write("before.tsv", beforeText);
		after=write("after.tsv", afterText);
		invoke("result", before, after);
		final String qualifier='\t'+CompareNcrnaScenarios.QUALIFICATION;
		final String summaryQualifier="\t0.500000"+qualifier;
		final String[] summary=lines("result/summary.tsv");
		check(summary.length==4, "All three explicitly projected families appear, including excluded calls");
		check(find(summary, "F\t").equals("F\t2\t3\t3\t2\t5\t1\t0\t1\t0\t2\t2\t1\t3\t2\t2\t1\t0\t-1\t-1\t-1\t2"+summaryQualifier),
			"Hypothetical F before=3/3/2 and after=2/2/1; two masked calls explain the denominator change");
		check(find(summary, "G\t").equals("G\t2\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1\t1\t0\t0\t0\t0\t0\t0\t1\t0"+summaryQualifier),
			"Relabeling the reference produces a G miss and does not excuse the F call");
		check(find(summary, "X\t").equals("X\t2\t0\t0\t0\t0\t0\t0\t0\t1\t0\t0\t0\t0\t0\t0\t0\t1\t0\t0\t0\t0"+summaryQualifier), "Excluded queries are retained in both scenarios");
		final String[] counts=lines("result/counts.tsv");
		check(find(counts, "FAMILY\tF\t.\tbefore\t").equals("FAMILY\tF\t.\tbefore\t5\t1\t1\t0\t1"+qualifier), "All before dispositions counted");
		check(find(counts, "FAMILY\tF\t.\tafter\t").equals("FAMILY\tF\t.\tafter\t3\t2\t0\t1\t1"+qualifier), "Unresolved to removed stays unscored");
		check(find(counts, "MODEL\tG\tM\tafter\t").equals("MODEL\tG\tM\tafter\t1\t0\t0\t0\t0"+qualifier), "Model counts retain changed family assignment");
		final String[] changes=lines("result/changes.tsv");
		check(find(changes, "A\tfamily\tPAIRED\t").endsWith("\tfalse\ttrue\tfalse\tfalse"+qualifier), "A family change is paired by declared decision provenance");
		check(find(changes, "A\tmove\tPAIRED\t").endsWith("\ttrue\tfalse\tfalse\tfalse"+qualifier), "Moved coordinates are not inferred additions/removals");
		check(find(changes, "A\tmissing\tREMOVED\t").contains("\tABSENT\t"), "Absent decision differs from disposition=removed");
		check(find(changes, "A\tnew\tADDED\t").contains("\tABSENT\t"), "New declared decision remains an addition");
		final String[] outcomes=lines("result/calls.tsv.gz");
		check(outcomes.length==8, "Every generic call survives, including a duplicate and an exclusion");
		check(find(outcomes, "A\t2\tF\t").contains("\tFP\t.\t.\t.\ttrue\ttrue\tMASKED\t"), "Scored eligibility protects duplicate FP before masking");
		check(find(outcomes, "A\t5\tF\t").endsWith("\ttrue"+qualifier), "Changed physical witness is reported despite unchanged TP outcome");
		check(find(outcomes, "A\t6\tF\t").endsWith("\tfalse"+qualifier), "Unresolved/removed rows provide no scoring authority");
		check(find(outcomes, "A\t7\tX\t").endsWith("\tfalse"+qualifier), "Context does not alter excluded-call outcome");
		check(lines("result/inputs.tsv")[2].endsWith("\t0\t0\t0\t0\t0\t0\t1\t1"+qualifier), "Empty call input remains explicit with one decision on each side");
		final String[] retained=lines("result/dispositions.tsv.gz");
		check(retained.length==17 && retained[1].equals("before\t2\t"+mask.trim()), "Every original row and source line is retained");
		check(!new File(path("result/reference.tsv")).exists() && !new File(path("result/input_1.gff.gz")).exists(), "A diagnostic report never publishes a reference manifest or GFF");
		final String reordered=write("reordered.tsv", CompareNcrnaScenarios.HEADER+'\n'+empty+context+unknown+missing+move+family+otherMask+mask);
		invoke("same", before, reordered);
		for(final String line:lines("same/calls.tsv.gz")){
			if(!line.startsWith("input_id\t")){check(line.endsWith("\tfalse"+qualifier), "Row order alone cannot change chosen witnesses/outcomes");}
		}
		bad("duplicate_key", beforeText+mask, "Duplicate decision_id");
		bad("duplicate_physical", beforeText+mask.replace("mask\tA", "other_id\tA"), "Duplicate/contradictory physical");
		bad("unresolved_mask", beforeText.replace(otherMask, otherMask.replace("RESOLVED", "UNRESOLVED")), "resolved included families");
		bad("outside", beforeText.replace(mask, mask.replace("\t10\t29\t", "\t290\t310\t")), "valid observed bounds");
		bad("missing_provenance", beforeText.replace("fixture:decision", "."), "decision and witness pointers");
		bad("unknown_input", beforeText.replace("\tA\t", "\tC\t"), "absent from the common");
		invoke("half", before, after, "overlap=0.5");
		check(Arrays.equals(summary, lines("half/summary.tsv")) && Arrays.equals(outcomes, lines("half/calls.tsv.gz")),
			"Explicit half overlap must match the omitted default");
		invoke("strict", before, after, "overlap=1");
		check(find(lines("strict/summary.tsv"), "F\t").equals("F\t2\t3\t3\t2\t5\t1\t0\t1\t0\t1\t3\t2\t3\t2\t2\t1\t0\t-2\t0\t0\t2\t1"+qualifier),
			"Strict overlap loses only the moved-after TP to FP+FN; masks and before counts stay unchanged");
		for(final String value:new String[]{"0", "1.01", "NaN"}){
			final String directory="invalid_overlap_"+value;
			try{
				invoke(directory, before, after, "overlap="+value);
				throw new AssertionError("Invalid overlap must fail: "+value);
			}catch(IllegalArgumentException failure){
				check(failure.getMessage().contains("Reciprocal scored overlap"), "Invalid overlap has a specific diagnostic");
				check(!new File(path(directory)).exists(), "Invalid overlap fails before directory creation");
			}
		}
		System.err.println("COMPARE_NCRNA_SCENARIOS_TEST_PASS");
	}

	private void bad(final String name, final String text, final String expected){
		final String file=write(name+".tsv", text);
		try{invoke(name, file, after);}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(expected), "Unexpected scenario failure: "+failure);
			check(!new File(path(name+"/summary.tsv")).exists(), "Failure must not publish the diagnostic completion summary");
			return;
		}
		throw new AssertionError("Expected scenario failure: "+expected);
	}
	private void invoke(final String name, final String first, final String second, final String... extra){
		final String[] base={"in="+first, "in2="+second, "calls="+calls, "contigs="+contigs,
			"familymap="+families, "out="+path(name), "t=1", "bgzip=f", "zipthreads=1"};
		final String[] args=Arrays.copyOf(base, base.length+extra.length);
		System.arraycopy(extra, 0, args, base.length, extra.length);
		CompareNcrnaScenarios.main(args);
	}
	private static String call(final int start, final int stop, final String model, final String family){
		return "chr\tBBTools\tRNA\t"+start+'\t'+stop+"\t.\t+\t.\tmodel:"+model+",family:"+family+'\n';
	}
	private static String row(final String id, final String input, final int start, final int stop, final String family, final String state, final boolean resolved){
		return fields(id, input, "chr", Integer.toString(start), Integer.toString(stop), "+", family,
			resolved ? "RESOLVED" : "UNRESOLVED", state, "proposed_rule", "M", state.equals("scored") ? "." : "uncertain",
			"fixture:decision", "fixture:witness");
	}
	private static String fields(final String... fields){
		final ByteBuilder b=new ByteBuilder();
		for(int i=0; i<fields.length; i++){
			if(i>0){b.tab();}
			BuildNcrnaReference.field(b, fields[i]);
		}
		b.nl();
		return new String(b.array, 0, b.length(), StandardCharsets.UTF_8);
	}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private String write(final String name, final String text){
		final String path=path(name);
		try(Output out=new Output(path)){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
		return path;
	}
	private String[] lines(final String name){
		final ArrayList<String> lines=new ArrayList<String>();
		try(Input in=new Input(path(name))){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				lines.add(new String(line, StandardCharsets.UTF_8));
			}
		}
		final int width=columns(lines.get(0));
		for(final String line:lines){check(columns(line)==width, "Diagnostic table header and rows must have equal widths: "+name);}
		return lines.toArray(new String[lines.size()]);
	}
	private static int columns(final String line){
		int count=1;
		for(int i=0; i<line.length(); i++){if(line.charAt(i)=='\t'){count++;}}
		return count;
	}
	private static String find(final String[] lines, final String prefix){
		for(final String line:lines){if(line.startsWith(prefix)){return line;}}
		throw new AssertionError("Missing row: "+prefix+" in "+Arrays.toString(lines));
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final File root;
	private String before, after, calls, contigs, families;
}
