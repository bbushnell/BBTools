package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import parse.LineParser1;
import prok.GradeNcrnaReference;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceGff;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;

/** Three-input, hand-calculated full grading with raw-row reconstruction.
 * Equal contig names/coordinates occur in different inputs deliberately.
 * @author Brian Bushnell, Raiden
 */
final class GradeNcrnaReferenceTest{

	GradeNcrnaReferenceTest(final File root_){root=root_;}

	void run(){
		check(root.mkdir(), "Grader fixture directory must be new");
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final Row a=codec.create("A", new Locus("seq", "A", 10, 29, '+'), false, null, null, "prior:a");
		final Row am=codec.create("A", new Locus("seq", "A", 90, 109, '+'), true, "partial", null, "prior:am");
		final Row b=codec.create("B", new Locus("seq", "A", 10, 29, '+'), false, null, null, "prior:b");
		final Row bm=codec.create("B", new Locus("seq", "A", 200, 219, '+'), true, "edge", null, "prior:bm");
		final Row c=codec.create("C", new Locus("seq", "A", 10, 29, '-'), false, null, null, "prior:c");
		codec.write(path("a.ref.gff"), "A", Arrays.asList(a, am));
		codec.write(path("b.ref.gff.gz"), "B", Arrays.asList(b, bm));
		codec.write(path("c.ref.gff"), "C", Collections.singletonList(c));
		write("a.calls.gff", "# calls\n"+call("A", 10, 29, '+')+call("A", 10, 29, '+')
			+call("A", 90, 109, '+')+call("B", 90, 109, '+')+call("TPP", 30, 40, '+'));
		write("b.calls.gff", "");
		write("c.calls.gff", call("A", 11, 29, '-'));
		rules=write("families.tsv", GradeNcrnaReference.RULE_HEADER+"\nA\tA\tA\ttrue\nB\tB\tB\ttrue\nTPP\tTPP\tTPP\tfalse\n");
		final String referenceRows=binding("A", "a.ref.gff")+binding("B", "b.ref.gff.gz")+binding("C", "c.ref.gff");
		final String queryRows=binding("C", "c.calls.gff")+binding("B", "b.calls.gff")+binding("A", "a.calls.gff");
		refs=write("references.tsv", GradeNcrnaReference.MANIFEST_HEADER+'\n'+referenceRows);
		queries=write("queries.tsv", GradeNcrnaReference.MANIFEST_HEADER+'\n'+queryRows);
		invoke("good", refs, queries);
		final String[] summary=lines("good.summary.tsv");
		check(summary.length==4 && summary[0].equals(GradeNcrnaReference.SUMMARY_HEADER), "All three explicitly configured families remain in summary");
		check(summary[1].equals("A\t3\t2\t1\t1\t3\t2\t1\t1\t0\t2\t1\t1\t0.666667\t0.666667\t0.500000\tDECLARED_INPUT_BINDINGS_EXPLICIT_REFERENCE_STATUS"),
			"A has two TPs, one duplicate FP, one empty-input FN, two masks and one masked call; reverse endpoint is off by one");
		// ByteBuilder.append(double, decimals) writes exact integers through append(long).
		check(summary[2].equals("B\t3\t0\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tNA\t0.500000\tDECLARED_INPUT_BINDINGS_EXPLICIT_REFERENCE_STATUS"),
			"Wrong-family query is an FP despite overlapping the other family's mask");
		check(summary[3].equals("TPP\t3\t0\t0\t0\t0\t0\t0\t0\t1\t0\t0\t0\tNA\tNA\t0.500000\tDECLARED_INPUT_BINDINGS_EXPLICIT_REFERENCE_STATUS"),
			"Explicit exclusion is visible without pretending to be a scored negative");
		final String[] raw=lines("good.raw.tsv");
		check(raw.length==12 && raw[0].equals(GradeNcrnaReference.RAW_HEADER), "Five reference and six generic query rows must survive");
		final LineParser1 p=new LineParser1('\t');
		int tp=0, fp=0, fn=0, masked=0, excluded=0, maskRefs=0, duplicates=0;
		boolean reverse=false, missingInput=false;
		for(int i=1; i<raw.length; i++){
			p.set(raw[i].getBytes(StandardCharsets.UTF_8));
			check(p.terms()==18, "Every raw row obeys its eighteen-field schema");
			if(p.termEquals("QUERY", 1)){
				if(p.termEquals("TP", 7)){
					tp++;
					check(!p.termEquals(".", 10), "Every TP query identifies its paired reference");
				}
				if(p.termEquals("FP", 7)){fp++;}
				if(p.termEquals("MASKED", 7)){
					masked++;
					check(p.termEquals(am.id, 12), "Masked query points to its same-input/family mask");
				}
				if(p.termEquals("EXCLUDED", 7)){excluded++;}
				if(p.termEquals("true", 13)){duplicates++;}
				if(p.termEquals("C", 0)){
					reverse=p.parseInt(16)==0 && p.parseInt(17)==-1;
				}
			}else{
				if(p.termEquals("FN", 7)){
					fn++;
					missingInput=p.termEquals("B", 0) && p.termEquals(b.id, 8);
				}
				if(p.termEquals("DONT_CARE", 7)){maskRefs++;}
			}
		}
		check(tp==2 && fp==2 && fn==1 && masked==1 && excluded==1 && maskRefs==2 && duplicates==1,
			"Raw outcomes independently reconstruct the complete scored/masked/excluded populations");
		check(reverse && missingInput, "Signed reverse offsets and zero-call input FNs survive exact-input partitioning");
		final String[] inputs=lines("good.inputs.tsv");
		check(inputs.length==4 && inputs[1].endsWith("\t5\t0\t0\t0\t4\t1")
			&& inputs[2].endsWith("\t0\t0\t0\t0\t0\t0") && inputs[3].endsWith("\t1\t0\t0\t0\t1\t0"),
			"Input report preserves explicit empty query files and exclusion counts");
		invoke("explicit_half", refs, queries, "0.5");
		for(final String suffix:new String[]{".raw.tsv", ".inputs.tsv", ".summary.tsv"}){
			check(Arrays.equals(lines("good"+suffix), lines("explicit_half"+suffix)),
				"Explicit overlap=0.5 must reproduce the documented default: "+suffix);
		}
		invoke("strict", refs, queries, "1");
		check(lines("strict.summary.tsv")[1].equals("A\t3\t1\t2\t2\t3\t2\t1\t1\t0\t1\t1\t1\t0.333333\t0.333333\t1\tDECLARED_INPUT_BINDINGS_EXPLICIT_REFERENCE_STATUS"),
			"overlap=1 removes only the reverse off-by-one TP, producing an extra FP+FN with unchanged masks");
		for(final String invalid:new String[]{"0", "1.01", "NaN"}){
			final String name="bad_overlap_"+invalid;
			expect(()->invoke(name, refs, queries, invalid), "reciprocal overlap");
			for(final String suffix:new String[]{".raw.tsv", ".inputs.tsv", ".summary.tsv"}){
				check(!new File(path(name+suffix)).exists(), "Invalid overlap must fail before any output: "+invalid);
			}
		}

		final String duplicate=write("duplicate.tsv", GradeNcrnaReference.MANIFEST_HEADER+'\n'+queryRows+binding("A", "a.calls.gff"));
		expect(()->invoke("duplicate", refs, duplicate), "Duplicate input_id");
		final String missing=write("missing.tsv", GradeNcrnaReference.MANIFEST_HEADER+'\n'+binding("A", "a.calls.gff"));
		expect(()->invoke("missing", refs, missing), "equal input populations");
		final String wrong=write("wrong.tsv", (GradeNcrnaReference.MANIFEST_HEADER+'\n'+queryRows).replace(path("A.fa"), path("wrong.fa")));
		expect(()->invoke("wrong", refs, wrong), "FASTA paths disagree");
		final String repeatedFile=write("repeated_file.tsv", (GradeNcrnaReference.MANIFEST_HEADER+'\n'+queryRows).replace(path("b.calls.gff"), path("a.calls.gff")));
		expect(()->invoke("repeated_file", refs, repeatedFile), "One GFF");
		final String badRef=write("bad.ref.gff", join(lines("a.ref.gff")).replace("reference_status=scored", "reference_status=unknown"));
		final String badManifest=write("bad_refs.tsv", (GradeNcrnaReference.MANIFEST_HEADER+'\n'+referenceRows).replace(path("a.ref.gff"), badRef));
		expect(()->invoke("bad_status", badManifest, queries), "reference_status");
		for(final String name:new String[]{"duplicate", "missing", "wrong", "repeated_file", "bad_status"}){
			check(!new File(path(name+".summary.tsv")).exists(), "A failed grade must not publish a completion summary: "+name);
		}
		checkRestriction(codec, a, am, referenceRows, queryRows);
		System.err.println("GRADE_NCRNA_REFERENCE_TEST_PASS");
	}
	/** Global reference restriction must reproduce a manually selected reference
	 * while keeping empty-query FNs and making each kind of exclusion observable. */
	private void checkRestriction(final NcrnaReferenceGff codec, final Row a, final Row am,
		final String referenceRows, final String queryRows){
		final Row excluded=codec.create("A", new Locus("seq", "TPP", 30, 40, '+'), false, null, null, "prior:excluded");
		final Row unknown=codec.create("A", new Locus("seq", "unsupported", 120, 139, '+'), false, null, null, "prior:unsupported");
		final Row unknownMask=codec.create("A", new Locus("seq", "unsupported", 140, 159, '+'), true, "partial", null, "prior:unsupported_mask");
		codec.write(path("a.global.gff"), "A", Arrays.asList(excluded, a, unknown, am, unknownMask));
		final String expandedRows=referenceRows.replace(path("a.ref.gff"), path("a.global.gff"));
		final String expanded=write("expanded.tsv", GradeNcrnaReference.MANIFEST_HEADER+'\n'+expandedRows);
		final String global=write("global.tsv", GradeNcrnaReference.MANIFEST_HEADER+'\n'+expandedRows+binding("D", "absent_off_panel.gff"));
		check(!new File(path("absent_off_panel.gff")).exists(), "Off-panel payload must be absent so a hidden read cannot pass");
		expect(()->invoke("strict_population", global, queries), "equal input populations");
		expect(()->invoke("strict_family", expanded, queries), "lacks an included query projection");
		invoke("restricted", global, queries, null, true);
		for(final String suffix:new String[]{".raw.tsv", ".summary.tsv"}){
			check(Arrays.equals(lines("good"+suffix), lines("restricted"+suffix)),
				"Explicit restriction must equal the independently authored manual reference, including empty-query FNs: "+suffix);
		}
		check(!new File(path("good.selection.tsv")).exists(), "Default strict mode must not acquire an extra output");
		final String[] selection=lines("restricted.selection.tsv");
		check(selection[0].equals(GradeNcrnaReference.SELECTION_HEADER), "Selection units must have an explicit schema");
		check(selection[1].equals("INPUT_POPULATION\t.\t.\t3\t1\t.\t.\t.\t.\tEXACT_QUERY_INPUT_SET_OFF_PANEL_GFFS_NOT_READ"),
			"One omitted INPUT is not an invented off-panel gene-row count");
		check(Arrays.asList(selection).contains("REFERENCE_FAMILY\tA\tTPP\t.\t.\t0\t0\t1\t0\tEXCLUDED_QUERY_FAMILY"),
			"Explicitly excluded family retains its omitted scored-row count");
		check(Arrays.asList(selection).contains("REFERENCE_FAMILY\tA\tunsupported\t.\t.\t0\t0\t1\t1\tUNMAPPED_REFERENCE_FAMILY"),
			"Unmapped reference families retain separate scored and masked omission counts");
		check(Arrays.asList(selection).contains("REFERENCE_FAMILY\tB\tA\t.\t.\t1\t1\t0\t0\tINCLUDED_PROJECTION"),
			"An empty query input retains both eligible scored and masked references");
		final String noB=write("missing_reference.tsv", GradeNcrnaReference.MANIFEST_HEADER+'\n'+binding("A", "a.global.gff")+binding("C", "c.ref.gff"));
		expect(()->invoke("restricted_missing", noB, queries, null, true), "Missing reference input: B");
		final String wrong=write("restricted_wrong.tsv", (GradeNcrnaReference.MANIFEST_HEADER+'\n'+queryRows).replace(path("B.fa"), path("wrong_B.fa")));
		expect(()->invoke("restricted_wrong", global, wrong, null, true), "FASTA paths disagree");
		final String duplicate=write("restricted_duplicate.tsv", GradeNcrnaReference.MANIFEST_HEADER+'\n'+queryRows+binding("A", "a.calls.gff"));
		expect(()->invoke("restricted_duplicate", global, duplicate, null, true), "Duplicate input_id");
		for(final String name:new String[]{"restricted_missing", "restricted_wrong", "restricted_duplicate"}){
			check(!new File(path(name+".summary.tsv")).exists(), "Invalid restriction cannot publish a completed grade");
		}
		System.err.println("GRADE_NCRNA_REFERENCE_RESTRICTION_TEST_PASS");
	}

	private void invoke(final String name, final String references, final String calls){
		invoke(name, references, calls, null);
	}
	private void invoke(final String name, final String references, final String calls, final String overlap){
		invoke(name, references, calls, overlap, false);
	}
	private void invoke(final String name, final String references, final String calls, final String overlap, final boolean restrict){
		final ArrayList<String> args=new ArrayList<String>(Arrays.asList("in="+calls, "reference="+references,
			"familymap="+rules, "out="+path(name), "t=1", "bgzip=f", "zipthreads=1"));
		if(overlap!=null){args.add("overlap="+overlap);}
		if(restrict){args.add("restrictreference=t");}
		GradeNcrnaReference.main(args.toArray(new String[args.size()]));
	}
	private String binding(final String id, final String gff){return id+'\t'+path(id+".fa")+'\t'+path(gff)+'\n';}
	private static String call(final String family, final int start, final int stop, final char strand){
		return "seq\tBBTools\tRNA\t"+start+'\t'+stop+"\t.\t"+strand+"\t0\tRNA,model:"+family+"_c1,family:"+family+'\n';
	}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private String write(final String name, final String text){
		final String path=path(name);
		try(Output out=new Output(path)){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
		return path;
	}
	private String[] lines(final String name){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path(name))){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				result.add(new String(line, StandardCharsets.UTF_8));
			}
		}
		return result.toArray(new String[result.size()]);
	}
	private static String join(final String[] lines){
		final StringBuilder out=new StringBuilder();
		for(final String line:lines){out.append(line).append('\n');}
		return out.toString();
	}
	private static void expect(final Runnable action, final String fragment){
		try{action.run();}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(fragment), "Unexpected grading failure: "+failure);
			return;
		}
		throw new AssertionError("Expected failure containing "+fragment);
	}
	private static void check(final boolean valid, final String reason){
		if(!valid){throw new AssertionError(reason);}
	}
	private final File root;
	private String rules, refs, queries;
}
