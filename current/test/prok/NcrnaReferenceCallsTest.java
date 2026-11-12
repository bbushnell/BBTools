package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import prok.NcrnaReferenceCalls;
import prok.NcrnaReferenceCalls.Result;
import prok.NcrnaReferenceCalls.Rule;
import prok.NcrnaTextIO.Output;

/** Literal CallGenes-format queries and independently specified population counts.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceCallsTest{

	static void run(final File root){
		check(root.mkdir(), "Query reader fixture directory must be new");
		final Rule[] rules={new Rule("FMN", "FMN", "FMN", true),
			new Rule("SRP", "SRP_RF00169", "SRP RNA (small)", true),
			new Rule("SRP", "SRP_RF01854", "SRP RNA (large)", true),
			new Rule("TPP", "TPP", "TPP", false)};
		final NcrnaReferenceCalls reader=new NcrnaReferenceCalls(rules);
		final String a="RNA,startScr:0,model:FMN_c1 n=7,family:FMN;partial=10";
		final String legacy="RNA,startScr:0,model:FMN_c1 n=7;partial=10";
		final String small="RNA,model:SRP_RF00169_c2,family:SRP RNA (small)";
		final String large="RNA,model:SRP_RF01854_c3,family:SRP RNA (large)";
		final String tpp="RNA,model:TPP_c1,family:TPP";
		final String rows="# fixture\n"+row("CDS", "CDS,start:ATG")+row("tRNA", "tRNA,model:tRNA_x")
			+row("rRNA", "16S,model:16S_x")+row("RNA", a)+row("RNA", legacy)
			+row("RNA", small)+row("RNA", large)+row("RNA", tpp);
		final String path=write(root, "calls.gff.gz", rows);
		final Result result=reader.read(path, "assembly-A");
		check(result.features==8 && result.cds==1 && result.trna==1 && result.rrna==1
			&& result.included==4 && result.excluded==1 && result.calls.size()==5,
			"Eight features partition into three non-targets, four included calls and one explicitly excluded call");
		check(Arrays.equals(result.ruleCounts, new long[]{2, 1, 1, 1}), "Every generic RNA row is counted under its exact rule");
		check(result.input.equals("assembly-A") && result.calls.get(0).line==5 && result.calls.get(4).line==9,
			"Exact declared input and original physical lines survive comments and ignored feature types");
		check(result.calls.get(0).locus.seqid.equals("seq%20Ω") && result.calls.get(0).model.equals("FMN_c1 n=7"),
			"CallGenes identifiers are literal UTF-8; percent signs and model metadata are not GFF3-decoded");
		check(result.calls.get(0).locus.start==10 && result.calls.get(0).locus.stop==29
			&& result.calls.get(0).locus.strand=='-' && result.calls.get(1).locus.start==10,
			"Coordinates, strand and duplicate calls survive; partial metadata does not trim or mask queries");
		check(result.calls.get(2).locus.family.equals("SRP") && result.calls.get(3).locus.family.equals("SRP")
			&& !result.calls.get(4).included, "Two explicit subtype prefixes project to one union; excluded calls remain visible");
		final Result second=reader.read(path, "assembly-B");
		check(second.input.equals("assembly-B") && result.input.equals("assembly-A") && second.features==8,
			"Repeated reader use neither accumulates counts nor mutates an earlier input's results");
		final String empty=write(root, "empty.gff", "# no calls\n");
		check(reader.read(empty, "assembly-A").features==0, "An empty query file is a counted zero-call input");

		bad(root, reader, "unknown", a.replace("FMN_c1", "FMNX_c1"), "Unknown generic");
		bad(root, reader, "conflict", a.replace("family:FMN", "family:TPP"), "conflicts");
		bad(root, reader, "duplicate_model", a+",model:FMN_c2", "Duplicate model");
		bad(root, reader, "duplicate_family", a+",family:FMN", "Duplicate family");
		bad(root, reader, "missing_model", "RNA,notmodel:FMN_c1,family:FMN", "Missing model");
		bad(root, reader, "empty_model", "RNA,model:,family:FMN", "Missing model");
		bad(root, reader, "excluded_conflict", tpp.replace("family:TPP", "family:FMN"), "conflicts");
		final String unknownType=write(root, "unknown_type.gff", row("ncRNA", a));
		expect(()->reader.read(unknownType, "A"), "Unknown CallGenes feature");
		final String badStrand=write(root, "bad_strand.gff", row("RNA", a).replace("\t-\t", "\t.\t"));
		expect(()->reader.read(badStrand, "A"), "strand");
		final String badRange=write(root, "bad_range.gff", row("RNA", a).replace("\t10\t29\t", "\t10\t2147483648\t"));
		expect(()->reader.read(badRange, "A"), "coordinates");
		final NcrnaReferenceCalls ambiguous=new NcrnaReferenceCalls(new Rule[]{rules[0], new Rule("OTHER", "FMN_c1", "FMN", true)});
		expect(()->ambiguous.read(path, "A"), "Ambiguous");
		expect(()->new NcrnaReferenceCalls(new Rule[]{rules[0], rules[0]}), "Duplicate model prefix");
		System.err.println("NCRNA_REFERENCE_CALLS_TEST_PASS");
	}

	private static void bad(final File root, final NcrnaReferenceCalls reader, final String name,
		final String attributes, final String fragment){
		final String path=write(root, name+".gff", row("RNA", attributes));
		expect(()->reader.read(path, "A"), fragment);
	}
	private static String row(final String type, final String attributes){
		return "seq%20Ω\tBBTools\t"+type+"\t10\t29\t.\t-\t0\t"+attributes+'\n';
	}
	private static String write(final File root, final String name, final String text){
		final String path=new File(root, name).getPath();
		try(Output out=new Output(path)){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
		return path;
	}
	private static void expect(final Runnable action, final String fragment){
		try{action.run();}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(fragment), "Unexpected query failure: "+failure);
			return;
		}
		throw new AssertionError("Expected failure containing "+fragment);
	}
	private static void check(final boolean valid, final String reason){
		if(!valid){throw new AssertionError(reason);}
	}
}
