package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.InvalidPathException;
import java.util.ArrayList;
import java.util.Arrays;

import prok.BuildNcrnaReference;
import prok.GradeNcrnaReference;
import prok.NcrnaReferenceGff;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.WriteNcrnaReference;
import structures.ByteBuilder;

/** Synthetic declared-authority serialization, not biological approval or truth.
 * @author Brian Bushnell, Raiden
 */
final class WriteNcrnaReferenceTest{

	WriteNcrnaReferenceTest(final File root_){root=root_;}

	void run(){
		check(root.mkdir(), "Writer fixture directory must be new");
		final String a=path("A%raw.fa"), b=path("B.fa");
		final String target="seq%Ω";
		inputs=write("inputs.tsv", WriteNcrnaReference.INPUT_HEADER+"\nA\t"+a+"\nB\t"+b+'\n');
		final String contigText=WriteNcrnaReference.CONTIG_HEADER+'\n'+fields("A", a, target, "100")+fields("B", b, target, "100");
		contigs=write("contigs.tsv", contigText);
		final String ruleText=WriteNcrnaReference.RULE_HEADER+"\n"
			+fields("keep", "Fam", "scored", "APPROVED", "synthetic_fixture:review")
			+fields("mask", "Fam", "dont_care", "APPROVED", "synthetic_fixture:review")
			+fields("unknown", ".", "unresolved", "APPROVED", "synthetic_fixture:preserve")
			+fields("remove", "Fam", "removed", "APPROVED", "synthetic_fixture:remove")
			+fields("context", "Ctx", "context", "APPROVED", "synthetic_fixture:context");
		rules=write("rules.tsv", ruleText);
		final String keep=decision(target, "10", "29", "+", "Fam", "RESOLVED", "scored", "keep", "RF00050", ".");
		final String mask=decision(target, "40", "49", "-", "Fam", "RESOLVED", "dont_care", "mask", ".", "partial,edge");
		final String unknown=decision(".", ".", ".", ".", ".", "UNRESOLVED", "unresolved", "unknown", ".", "unassessed");
		final String removed=decision(target, "70", "79", "+", "Fam", "RESOLVED", "removed", "remove", ".", "reviewed_bad");
		final String context=decision(target, "60", "69", "+", "Ctx", "RESOLVED", "context", "context", ".", "context_only");
		final String decisionText=WriteNcrnaReference.DECISION_HEADER+'\n'+keep+mask+unknown+removed+context;
		decisions=write("decisions.tsv", decisionText);
		invoke("result", inputs, contigs, rules, decisions);
		check(Arrays.equals(lines("result/dispositions.tsv"), lines("decisions.tsv")), "All five original disposition rows survive unchanged");
		check(Arrays.equals(lines("result/rules.tsv"), lines("rules.tsv")), "Recorded approval/authority claims remain visible, not independently certified");
		final String[] summary=lines("result/summary.tsv");
		check(summary.length==3 && summary[1].equals("A\t5\t1\t1\t1\t1\t1\tDECLARED_APPROVED_RULES_NOT_INDEPENDENTLY_VERIFIED")
			&& summary[2].equals("B\t0\t0\t0\t0\t0\t0\tDECLARED_APPROVED_RULES_NOT_INDEPENDENTLY_VERIFIED"),
			"Every decision has exactly one disposition count; empty input remains explicit");
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final ArrayList<Row> rows=codec.read(path("result/input_1.gff.gz"), "A");
		check(rows.size()==2 && !rows.get(0).mask && rows.get(1).mask && rows.get(0).locus.seqid.equals(target),
			"Only resolved scored/masked rows enter the GFF, with literal Unicode/percent target identity");
		check(rows.get(1).reason.equals("partial,edge") && rows.get(0).evidence.equals(path("result/dispositions.tsv")+":2"),
			"Mask reason and exact companion row pointer survive serialization");
		check(codec.read(path("result/input_2.gff.gz"), "B").isEmpty(), "Zero-disposition input gets an explicitly bound empty GFF");
		final String[] manifest=lines("result/reference.tsv");
		check(manifest.length==3 && manifest[1].equals("A\t"+a+'\t'+path("result/input_1.gff.gz")),
			"Published manifest uses literal paths, not percent-escaped FASTA aliases");
		check(!new File(path("result/reference.partial.tsv")).exists(), "Successful publication leaves only the completed manifest");

		write("a.calls.gff", "");
		write("b.calls.gff", "");
		final String queries=write("queries.tsv", GradeNcrnaReference.MANIFEST_HEADER+"\nB\t"+b+'\t'+path("b.calls.gff")
			+"\nA\t"+a+'\t'+path("a.calls.gff")+'\n');
		final String families=write("families.tsv", GradeNcrnaReference.RULE_HEADER+"\nFam\tF\tFam\ttrue\n");
		GradeNcrnaReference.main(new String[]{"in="+queries, "reference="+path("result/reference.tsv"),
			"familymap="+families, "out="+path("graded"), "t=1"});
		check(lines("graded.summary.tsv")[1].equals("Fam\t2\t0\t0\t1\t1\t1\t0\t0\t0\t0\t0\t0\tNA\t0\t0.500000\tDECLARED_INPUT_BINDINGS_EXPLICIT_REFERENCE_STATUS"),
			"Materialized scored truth yields one FN for empty calls; the mask yields no FN");
		final String reorderedInputs=write("reordered_inputs.tsv", WriteNcrnaReference.INPUT_HEADER+"\nB\t"+b+"\nA\t"+a+'\n');
		final String reordered=write("reordered.tsv", WriteNcrnaReference.DECISION_HEADER+'\n'+context+removed+unknown+mask+keep);
		invoke("reordered", reorderedInputs, contigs, rules, reordered);
		final ArrayList<Row> again=codec.read(path("reordered/input_2.gff.gz"), "A");
		check(again.size()==2 && again.get(0).id.equals(rows.get(0).id) && again.get(1).id.equals(rows.get(1).id),
			"Locus IDs survive input/disposition reordering and different output filenames");

		final String draft=write("draft_rules.tsv", ruleText.replace("APPROVED", "DRAFT_NOT_APPLIED"));
		expect(()->invoke("draft", inputs, contigs, draft, decisions), "APPROVED");
		check(!new File(path("draft")).exists(), "Unapproved rules fail before creating a reference directory");
		bad("unknown_rule", decisionText.replace("\tkeep\t", "\tmissing\t"), "Missing explicit approved rule");
		bad("mask_unresolved", decisionText.replace(mask, mask.replace("RESOLVED", "UNRESOLVED")), "resolved families");
		bad("mask_reason", decisionText.replace(mask, mask.replace("partial,edge", ".")), "explicit reason");
		bad("outside", decisionText.replace(keep, keep.replace("\t10\t29\t", "\t90\t110\t")), "valid observed contig");
		bad("duplicate", decisionText+keep, "Duplicate or contradictory");
		bad("contradictory", decisionText+removed.replace("\t70\t79\t", "\t10\t29\t"), "Duplicate or contradictory");
		bad("mismatch", decisionText.replace("\tscored\tkeep\t", "\tdont_care\tkeep\t"), "conflicts with its approved rule");
		final String wrongContigs=write("wrong_contigs.tsv", contigText.replace(fields("A", a, target, "100"), fields("A", path("wrong.fa"), target, "100")));
		expect(()->invoke("wrong_contigs", inputs, wrongContigs, rules, decisions), "exact roster input/FASTA");
		check(!new File(path("wrong_contigs")).exists(), "Invalid observation binding fails before creating output");
		for(final String badPath:new String[]{"tab\toutput", "cr\routput", "lf\noutput"}){
			expect(()->invoke(badPath, inputs, contigs, rules, decisions), "no tab/CR/LF");
			check(!new File(path(badPath)).exists(), "Invalid manifest paths fail before directory creation");
		}
		boolean unicodePath=true;
		try{new File(path("percent%unicodeΩ")).toPath();}
		catch(InvalidPathException unsupported){unicodePath=false;}
		if(unicodePath){
			invoke("percent%unicodeΩ", inputs, contigs, rules, decisions);
			check(lines("percent%unicodeΩ/reference.tsv")[1].equals("A\t"+a+'\t'+path("percent%unicodeΩ/input_1.gff.gz")),
				"Percent and Unicode output paths remain literal on a supporting filesystem");
			System.err.println("WRITE_NCRNA_REFERENCE_UNICODE_PATH_PASS");
		}else{
			expect(()->invoke("percent%unicodeΩ", inputs, contigs, rules, decisions), "unmappable characters");
			check(!new File(path("percent%unicodeΩ")).exists(), "An unsupported name must not create a replacement-character directory");
			System.err.println("WRITE_NCRNA_REFERENCE_UNICODE_PATH_REJECTED");
		}
		rebaseCases(a, b);
		System.err.println("WRITE_NCRNA_REFERENCE_TEST_PASS");
	}

	/** Literal old-pointer encoding and unknown attributes expose accidental global
	 * replacement or parse/serialize loss. Empty views and masks must survive.
	 */
	private void rebaseCases(final String a, final String b){
		final String old="/old%tree/witness.tsv", oldEncoded="%2Fold%25tree%2Fwitness.tsv";
		final String witness=write("new%witness.tsv", "synthetic witness\n");
		final String header="##gff-version 3\n##reference-input-id A\n";
		final String comment="# The old path remains historical here: "+oldEncoded+"\n\n";
		final String sourceText=header+comment+feature(old, false, oldEncoded)+feature(old, true, oldEncoded);
		final String expectedText=header+comment+feature(witness, false, oldEncoded)+feature(witness, true, oldEncoded);
		final String source=write("rebase_source.gff.gz", sourceText);
		final String empty=write("rebase_empty.gff.gz", "##gff-version 3\n##reference-input-id B\n");
		final String manifest=write("rebase_inputs.tsv", GradeNcrnaReference.MANIFEST_HEADER+"\nA\t"+a+'\t'+source+"\nB\t"+b+'\t'+empty+'\n');
		invokeRebase("rebased", manifest, old, witness);
		check(joined("rebased/input_1.gff.gz").equals(expectedText), "Rebase may change only the evidence value, retaining comments, unknown fields and numeric spelling");
		check(joined("rebased/input_2.gff.gz").equals(joined("rebase_empty.gff.gz")), "Empty reference remains explicitly bound and unchanged");
		final ArrayList<Row> rows=new NcrnaReferenceGff().read(path("rebased/input_1.gff.gz"), "A");
		check(rows.size()==2 && !rows.get(0).mask && rows.get(1).mask && rows.get(1).reason.equals("partial")
			&& rows.get(0).evidence.equals(witness) && rows.get(1).evidence.equals(witness), "Existing reference reader must consume rebased scored and masked rows");
		check(rows.get(0).scores.bitsTenths(0)==100 && rows.get(0).scores.gaBits(0).equals("020.0"),
			"Stored CM measurements retain their meaning and supplied GA text");
		check(lines("rebased/reference.tsv")[1].equals("A\t"+a+'\t'+path("rebased/input_1.gff.gz")), "Only the GFF path changes in each literal manifest row");
		check(joined("rebase_source.gff.gz").equals(sourceText), "Rebasing must never modify its source GFF");
		check(!new File(path("rebased/reference.partial.tsv")).exists(), "Completed rebase publishes the manifest last");
		expect(()->invokeRebase("rebased", manifest, old, witness), "must not exist");
		expect(()->invokeRebase("missing_witness", manifest, old, path("absent.tsv")), "destination witness must already exist");
		check(!new File(path("missing_witness")).exists(), "Missing witness fails before creating output");
		badRebase("wrong_evidence", sourceText.replace("reference_evidence="+oldEncoded, "reference_evidence=elsewhere"), a, old, witness, "old evidence pointer");
		badRebase("missing_evidence", sourceText.replace(";reference_evidence="+oldEncoded, ""), a, old, witness, "Missing reference");
		badRebase("duplicate_evidence", sourceText.replace(";reference_evidence="+oldEncoded,
			";reference_evidence="+oldEncoded+";reference_evidence="+oldEncoded), a, old, witness, "old evidence pointer");
		badRebase("raw_status", sourceText.replace("reference_status=scored", "reference_status=source"), a, old, witness, "scored/dont_care");
		badRebase("wrong_row_input", sourceText.replace(";input_id=A", ";input_id=B"), a, old, witness, "Feature input_id");
		badRebase("wrong_file_input", sourceText.replace("##reference-input-id A", "##reference-input-id B"), a, old, witness, "input directive");
		badRebase("repeated_header", sourceText+header, a, old, witness, "concatenate reference headers");
		System.err.println("WRITE_NCRNA_REFERENCE_REBASE_TEST_PASS");
	}
	private static String feature(final String witness, final boolean mask, final String historicalPointer){
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final Locus locus=new Locus("ctg%A", "Fam", mask ? 40 : 10, mask ? 50 : 20, mask ? '-' : '+');
		final Row row=codec.create("A", locus, mask, mask ? "partial" : null, "RF00050", witness);
		final ByteBuilder encoded=new ByteBuilder();
		codec.append(row, encoded);
		return new String(encoded.array, 0, encoded.length()-1, StandardCharsets.UTF_8)
			+";cm_models=RF00050;cm_bits=010.0;cm_ga=020.0;cm_percent_ga=50.000000"
			+";cm_score_basis=printed_0.1_bit;cm_percent_precision=6_decimal_places"
			+";future_score=01.2300;historical_pointer="+historicalPointer+"\n";
	}
	private void badRebase(final String name, final String gff, final String fasta,
		final String old, final String witness, final String fragment){
		final String source=write(name+".gff", gff);
		final String manifest=write(name+".manifest.tsv", GradeNcrnaReference.MANIFEST_HEADER+"\nA\t"+fasta+'\t'+source+'\n');
		expect(()->invokeRebase(name, manifest, old, witness), fragment);
		check(!new File(path(name+"/reference.tsv")).exists(), "A failed rebase cannot publish a completed manifest");
	}
	private void invokeRebase(final String name, final String manifest, final String old, final String witness){
		WriteNcrnaReference.main(new String[]{"mode=rebase", "in="+manifest, "out="+path(name),
			"evidencefrom="+old, "evidenceto="+witness, "t=1", "bgzip=f", "zipthreads=1"});
	}
	private String joined(final String file){
		final ByteBuilder text=new ByteBuilder();
		try(Input in=new Input(path(file))){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){text.append(line).nl();}
		}
		return new String(text.array, 0, text.length(), StandardCharsets.UTF_8);
	}

	private void bad(final String name, final String text, final String fragment){
		final String file=write(name+".tsv", text);
		expect(()->invoke(name, inputs, contigs, rules, file), fragment);
		check(!new File(path(name+"/reference.tsv")).exists(), "Failed materialization must never publish a completed manifest: "+name);
	}
	private void invoke(final String name, final String roster, final String observations, final String ruleFile, final String data){
		WriteNcrnaReference.main(new String[]{"in="+data, "inputs="+roster, "contigs="+observations,
			"rules="+ruleFile, "out="+path(name), "t=1", "bgzip=f", "zipthreads=1"});
	}
	private static String decision(final String target, final String start, final String stop, final String strand,
		final String family, final String state, final String disposition, final String rule, final String model, final String reason){
		return fields("A", target, start, stop, strand, family, state, disposition, rule, model, reason,
			"synthetic_fixture:decision", "synthetic_source%:witness");
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
		final String file=path(name);
		try(Output out=new Output(file)){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
		return file;
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
	private static void expect(final Runnable action, final String fragment){
		try{action.run();}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(fragment), "Unexpected materializer failure: "+failure);
			return;
		}
		throw new AssertionError("Expected failure containing "+fragment);
	}
	private static void check(final boolean valid, final String reason){
		if(!valid){throw new AssertionError(reason);}
	}
	private final File root;
	private String inputs, contigs, rules, decisions;
}
