package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaReferenceGff;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaReferenceMatcher;
import prok.NcrnaReferenceMatcher.Result;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Literal GFF bytes and matcher expectations, independent of admission policy.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceGffTest{

	static void run(final File root){
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		final Locus core=new Locus("seq", "A", 10, 29, '+');
		final Row scored=codec.create("1", core, false, null, null, "prior:7");
		check(scored.id.equals("ncrna|1:1|3:seq|1:A|10|29|+"), "ID has literal length-prefixed identity, coordinates and strand");
		final String expected="seq\tBBTools\tncRNA\t10\t29\t.\t+\t.\tID=ncrna%7C1%3A1%7C3%3Aseq%7C1%3AA%7C10%7C29%7C%2B"
			+";input_id=1;reference_family=A;reference_status=scored;reference_evidence=prior%3A7\n";
		check(text(codec, scored).equals(expected), "Writer must match independently spelled GFF3 bytes");
		final Row restored=parse(codec, expected, "1");
		check(!restored.mask && restored.id.equals(scored.id) && restored.evidence.equals("prior:7"),
			"Reader retains explicit status, identity and decoded decision pointer");
		final Row mask=codec.create("1", new Locus("seq", "A", 1, 100, '+'), true, "partial,edge", "RF00050", "review:8");
		final Row restoredMask=parse(codec, text(codec, mask), "1");
		check(restoredMask.mask && restoredMask.reason.equals("partial,edge") && restoredMask.model.equals("RF00050"),
			"Escaped mask reason and resolved model survive the GFF row");
		final Locus maskedCall=new Locus("seq", "A", 90, 109, '+');
		final Locus wrongFamily=new Locus("seq", "B", 90, 109, '+');
		final Result result=NcrnaReferenceMatcher.match(Arrays.asList(core, core, maskedCall, wrongFamily),
			new int[]{1, 1, 1, 1}, Arrays.asList(restoredMask.locus, restored.locus), new int[]{1, 1},
			new boolean[]{restoredMask.mask, restored.mask}, .5f);
		check(result.matched==1 && result.falsePositives==2 && result.falseNegatives==0 && result.maskedCalls==1,
			"Parsed GFF status yields one TP, duplicate and wrong-family FPs, and one same-family masked call");
		check(result.callMaskReference[1]<0 && result.callMaskReference[2]==0,
			"Scored eligibility defeats masking even after the scored reference is consumed");
		final String maskId=mask.id;
		codec.create("elsewhere", core, false, null, null, "unrelated");
		check(codec.create("1", mask.locus, true, "partial,edge", "RF00050", "different_evidence").id.equals(maskId),
			"IDs do not depend on creation order, evidence pointer or transient buffers");
		check(!codec.create("2", core, false, null, null, "prior:7").id.equals(scored.id),
			"Identical contig coordinates in different exact inputs have different IDs");
		final Row unicode=codec.create("Ω", new Locus("a%;=,\tΩ", "family|Ω", 1, 1, '-'),
			true, "edge", null, "decision;Ω");
		check(unicode.id.startsWith("ncrna|2:Ω|8:a%;=,\tΩ|9:family|Ω|"), "Identity lengths count UTF-8 bytes");
		final String unicodeText=text(codec, unicode);
		check(unicodeText.startsWith("a%25%3B%3D%2C%09%CE%A9\t"), "Seqid framing preserves literal delimiters and UTF-8");
		final Row unicodeBack=parse(codec, unicodeText, "Ω");
		check(unicodeBack.id.equals(unicode.id) && unicodeBack.locus.seqid.equals(unicode.locus.seqid),
			"Decode once: literal percent signs never become another escape pass");

		expect(()->parse(codec, expected, "2"), "manifest");
		expect(()->parse(codec, expected.replace(";reference_status=scored", ""), "1"), "reference_status");
		expect(()->parse(codec, expected.replace("reference_status=scored", "reference_status=unknown"), "1"), "reference_status");
		expect(()->parse(codec, expected.replace("reference_status=scored", "reference_status=scored;reference_status=scored"), "1"), "Duplicate");
		expect(()->parse(codec, expected.replace("reference_status=scored", "reference_status=dont_care"), "1"), "reason");
		expect(()->parse(codec, expected.replace("reference_status=scored", "reference_status=scored;dont_care_reason=edge"), "1"), "reason");
		expect(()->parse(codec, expected.replace("prior%3A7", "prior%QZ7"), "1"), "percent escape");
		expect(()->parse(codec, expected.replace("ID=ncrna", "ID=other"), "1"), "Reference ID");
		expect(()->parse(codec, expected.replace("\t10\t29\t", "\t10\t2147483648\t"), "1"), "coordinates");
		expect(()->parse(codec, expected.replace("\t+\t", "\t.\t"), "1"), "strand");
		expect(()->codec.create("\uD800", core, false, null, null, "prior"), "Unicode");
		expect(()->codec.create("\uFFFD", core, false, null, null, "prior"), "Unicode");
		files(root, codec, scored, mask, unicode);
		scores(root, codec);
		System.err.println("NCRNA_REFERENCE_GFF_TEST_PASS");
	}

	/** Model-specific decimal evidence survives serialization without inventing precision. */
	private static void scores(final File root, final NcrnaReferenceGff codec){
		final Row band=codec.create("1", new Locus("seq", "A", 10, 29, '+'), true,
			"marginal", "RF00050", "native:3", codec.scores("RF00050", new long[]{180}, new String[]{"20.00"}));
		final String bandText=text(codec, band);
		check(bandText.contains(";cm_bits=18.0;cm_ga=20.00;cm_percent_ga=90.000000;cm_score_basis=printed_0.1_bit;cm_percent_precision=6_decimal_places"),
			"Reported bits and original GA stay exact; display percentage declares its rounding");
		final Row bandBack=parse(codec, bandText, "1");
		check(bandBack.mask && bandBack.scores.bitsTenths(0)==180 && bandBack.scores.gaBits(0).equals("20.00"),
			"Scores are retained measurements, not inferred dispositions");
		final Row multiple=codec.create("1", new Locus("seq", "A", 30, 49, '+'), true,
			"competing_models", "RF00050,RF00162", "native:4,5",
			codec.scores("RF00050,RF00162", new long[]{10, 300}, new String[]{"3.00", "40.0"}));
		final String multiText=text(codec, multiple);
		check(multiText.contains(";cm_bits=1.0,30.0;cm_ga=3.00,40.0;cm_percent_ga=33.333333,75.000000"),
			"Repeating ratios use an explicitly rounded display while both exact source values survive");
		final Row multiBack=parse(codec, multiText, "1");
		check(multiBack.scores.size()==2 && multiBack.scores.bitsTenths(1)==300
			&& multiBack.scores.gaBits(0).equals("3.00") && text(codec, multiBack).equals(multiText),
			"Separate model witnesses cannot flatten to one score; full GFF bytes round-trip");
		final Row prior=codec.create("1", new Locus("seq", "A", 50, 69, '+'), false,
			null, "RF00050", "protected_prior:7",
			codec.scores("RF00050", new long[]{Long.MIN_VALUE}, new String[]{"20.00"}));
		final String priorText=text(codec, prior);
		check(priorText.contains(";cm_bits=.;cm_ga=20.00;cm_percent_ga=.;cm_score_basis=no_current_native_hit"),
			"Hitless protected priors carry present missing fields, never synthetic zero scores");
		final Row priorBack=parse(codec, priorText, "1");
		check(!priorBack.mask && priorBack.scores.missing && priorBack.scores.bitsTenths(0)==Long.MIN_VALUE,
			"Independent prior acceptance survives absent current measurements");
		final Row unknownGa=codec.create("1", new Locus("seq", "A", 70, 89, '+'), false,
			null, null, "protected_prior:8", codec.scores(null, new long[]{Long.MIN_VALUE}, new String[]{"."}));
		check(parse(codec, text(codec, unknownGa), "1").scores.gaBits(0).equals("."),
			"Unavailable prior GA remains unavailable");
		final Row halfway=codec.create("1", new Locus("seq", "A", 90, 109, '+'), true,
			"format_fixture", "RF00050", "native:9", codec.scores("RF00050", new long[]{1}, new String[]{"256.0"}));
		check(text(codec, halfway).contains(";cm_bits=0.1;cm_ga=256.0;cm_percent_ga=0.039063"),
			"Exact rational percentage0.0390625 rounds at the declared sixth place");
		final Row negative=codec.create("1", halfway.locus, true, "format_fixture", "RF00050", "native:10",
			codec.scores("RF00050", new long[]{-1}, new String[]{"256.0"}));
		check(text(codec, parse(codec, text(codec, negative), "1")).contains(";cm_bits=-0.1;cm_ga=256.0;cm_percent_ga=-0.039063"),
			"Negative reported scores retain their sign and round halves away from zero");
		expect(()->parse(codec, bandText.replace(";cm_ga=20.00", ""), "1"), "all six fields");
		expect(()->parse(codec, bandText.replace(";cm_models=RF00050", ""), "1"), "all six fields");
		expect(()->parse(codec, bandText.replace("90.000000", "89.000000"), "1"), "percentage");
		expect(()->parse(codec, bandText.replace("cm_bits=18.0", "cm_bits=18.00"), "1"), "one decimal");
		expect(()->parse(codec, bandText.replace("cm_bits=18.0", "cm_bits=18.0,19.0"), "1"), "equal lengths");
		expect(()->parse(codec, priorText.replace("cm_bits=.", "cm_bits=0.0"), "1"), "explicitly missing");
		expect(()->codec.scores("RF00050", new long[]{100}, new String[]{"0.00"}), "positive");
		expect(()->codec.scores("RF00050,RF00162", new long[]{180}, new String[]{"20.00"}), "corresponding");
		expect(()->codec.scores("RF00050,", new long[]{180, 190}, new String[]{"20.00", "20.00"}), "empty model-list");
		expect(()->codec.scores(",RF00162", new long[]{180, 190}, new String[]{"20.00", "20.00"}), "empty model-list");
		expect(()->codec.scores("RF00050,.", new long[]{180, 190}, new String[]{"20.00", "20.00"}), "missing model marker");
		final Row resolvedPrior=codec.create("1", prior.locus, false, null, "RF00050", "protected_prior:7", multiple.scores);
		final Row resolvedBack=parse(codec, text(codec, resolvedPrior), "1");
		check(resolvedBack.model.equals("RF00050") && resolvedBack.scores.models.equals("RF00050,RF00162"),
			"Current CM witnesses must not overwrite a protected prior's resolved model assignment");
		final String path=new File(root, "scores.gff.gz").getPath();
		codec.write(path, "1", Arrays.asList(band, multiple, prior, unknownGa));
		final ArrayList<Row> reread=codec.read(path, "1");
		check(reread.size()==4 && reread.get(1).scores.size()==2 && reread.get(2).scores.missing,
			"Compressed native I/O retains all score lists and missing-value declarations");
		System.err.println("NCRNA_REFERENCE_GFF_SCORES_TEST_PASS");
	}

	/** Real native I/O checks compressed output, empty-input binding and buffer reuse. */
	private static void files(final File root, final NcrnaReferenceGff codec, final Row scored,
		final Row mask, final Row unicode){
		check(root.mkdir(), "Reference GFF fixture directory must be new");
		final String header="##gff-version 3\n##reference-input-id 1\n";
		final String ordinary=new File(root, "ordinary.gff").getPath();
		codec.write(ordinary, "1", Arrays.asList(scored, mask));
		check(readText(ordinary).equals(header+text(codec, scored)+text(codec, mask)),
			"Native file writer emits the literal version/input directives followed by row bytes");
		final ArrayList<Row> restored=codec.read(ordinary, "1");
		check(restored.size()==2 && !restored.get(0).mask && restored.get(1).mask, "File reader retains scored and masked rows");
		final String empty=new File(root, "empty.gff.gz").getPath();
		codec.write(empty, "1", Collections.<Row>emptyList());
		check(readText(empty).equals(header) && codec.read(empty, "1").isEmpty(), "Empty references retain exact input identity");
		expect(()->codec.read(empty, "2"), "manifest");
		final String utf8=new File(root, "unicode.gff.gz").getPath();
		codec.write(utf8, "Ω", Collections.singletonList(unicode));
		check(codec.read(utf8, "Ω").get(0).id.equals(unicode.id), "Compressed file preserves UTF-8 input and row identity");
		final String invalid=new File(root, "invalid.gff").getPath();
		expect(()->codec.write(invalid, "1", Arrays.asList(scored, scored)), "Duplicate");
		check(!new File(invalid).exists(), "Duplicate rows must fail before opening an output");
		expect(()->codec.write(invalid, "2", Collections.singletonList(scored)), "declared file input");
		check(!new File(invalid).exists(), "Mixed input rows must fail before opening an output");
		final String duplicate=raw(root, "duplicate.gff", header+text(codec, scored)+text(codec, scored));
		expect(()->codec.read(duplicate, "1"), "Duplicate");
		final String noBinding=raw(root, "no_binding.gff", "##gff-version 3\n");
		expect(()->codec.read(noBinding, "1"), "input directive");
		final String concatenated=raw(root, "concatenated.gff", header+header);
		expect(()->codec.read(concatenated, "1"), "Repeated reference header");
		final String fasta=raw(root, "embedded_fasta.gff", header+"##FASTA\n>seq\nA\n");
		expect(()->codec.read(fasta, "1"), "embedded FASTA");
		final ArrayList<Row> many=new ArrayList<Row>();
		for(int i=0; i<200; i++){
			many.add(codec.create("1", new Locus("seq", "A", 2*i+1, 2*i+2, '+'), false, null, null, "fixture"));
		}
		final String buffered=new File(root, "buffered.gff.gz").getPath();
		codec.write(buffered, "1", many);
		final ArrayList<Row> reread=codec.read(buffered, "1");
		check(reread.size()==200, "Buffer flushing must preserve all 200 distinct rows");
		for(int i=0; i<reread.size(); i++){
			check(reread.get(i).locus.start==2*i+1 && reread.get(i).locus.stop==2*i+2,
				"Native writer copies the reusable buffer before it is cleared; wrong row="+i);
		}
		System.err.println("NCRNA_REFERENCE_GFF_FILES_TEST_PASS");
	}

	private static String raw(final File root, final String name, final String content){
		final String path=new File(root, name).getPath();
		try(Output out=new Output(path)){out.writer.print(content.getBytes(StandardCharsets.UTF_8));}
		return path;
	}
	private static String readText(final String path){
		final ByteBuilder bytes=new ByteBuilder();
		try(Input in=new Input(path)){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){bytes.append(line).append('\n');}
		}
		return new String(bytes.array, 0, bytes.length(), StandardCharsets.UTF_8);
	}

	private static String text(final NcrnaReferenceGff codec, final Row row){
		final ByteBuilder bytes=new ByteBuilder();
		codec.append(row, bytes);
		return new String(bytes.array, 0, bytes.length(), StandardCharsets.UTF_8);
	}
	private static Row parse(final NcrnaReferenceGff codec, final String line, final String input){
		return codec.parse(line.substring(0, line.length()-1).getBytes(StandardCharsets.UTF_8), input);
	}
	private static void expect(final Runnable action, final String fragment){
		try{action.run();}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(fragment), "Unexpected failure: "+failure);
			return;
		}
		throw new AssertionError("Expected failure containing "+fragment);
	}
	private static void check(final boolean valid, final String reason){
		if(!valid){throw new AssertionError(reason);}
	}
}
