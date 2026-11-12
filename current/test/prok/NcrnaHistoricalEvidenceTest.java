package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaEvidenceCatalog.Observation;
import prok.NcrnaHistoricalEvidence;
import prok.NcrnaHistoricalEvidence.Feature;
import prok.NcrnaHistoricalEvidence.Source;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Finite source-graph expectations, independent of biological truth admission.
 * Uses the template PreParser/Parser entry point and native byte output.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaHistoricalEvidenceTest{

	/*--------------------------------------------------------------*/
	/*----------------         Initialization       ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaHistoricalEvidenceTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
			}
			check(parser.out1!=null, "Supply out=new_fixture_directory");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Fixture directory must be new: "+root);
			Shared.setThreads(1);
			final Observation observed=observation(root);
			branches(root, observed);
			multipleParents(root, observed);
			multipart(root, observed);
			cycles(root, observed);
			attributes(root, observed);
			biotypes(root, observed);
			coordinates(root, observed);
			sourceStates(root, observed);
			malformed(root, observed);
			pp.outstream.println("NCRNA_HISTORICAL_EVIDENCE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Twenty observed bases; the literal target deliberately has a taxonomy-like suffix. */
	private static Observation observation(final File root){
		final String fasta=write(root, "one.fa", ">seq_tid_7 description\nAAAAAAAAAAAAAAAAAAAA\n");
		final String census=write(root, "census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'
			+"1\t"+fasta+'\t'+new File(fasta).length()+"\t20\t1\n");
		final String members=write(root, "members.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'
			+"fixture\t1\t1\t"+fasta+"\tNA\tNA\tNA\tprokaryotic\tNA\tNA\n");
		final NcrnaEvidenceCatalog catalog=NcrnaEvidenceCatalog.load(census, members);
		final Observation result=catalog.observe(catalog.get("1"));
		check(result.size()==1 && result.length("seq_tid_7")==20, "Fixture contigs must match the literal FASTA token");
		return result;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Structural cases      ----------------*/
	/*--------------------------------------------------------------*/

	/** A selected RNA brings its own descendants and ancestors, never their siblings. */
	private static void branches(final File root, final Observation observed){
		final String rows=row("region", "ID=chromosome")
			+row("gene", "ID=mixed;Parent=chromosome")
			+row("ncRNA", "ID=rna;Parent=mixed")
			+row("exon", "ID=rna_exon;Parent=rna")
			+row("mRNA", "ID=mrna;Parent=mixed")
			+row("CDS", "ID=protein;Parent=mrna")
			+row("gene", "ID=other_gene;Parent=chromosome")
			+row("CDS", "ID=other_protein;Parent=other_gene")
			+row("region", "ID=unused_region");
		final Source source=read(root, "branches", rows, observed);
		check(source.size()==4 && source.featureRows()==9 && source.excludedRows()==5,
			"Only chromosome, mixed gene, RNA and RNA exon survive the nine-row source");
		check(source.excluded("UNSELECTED_CODING")==3 && source.excluded("UNSELECTED_GENERIC_STRUCTURE")==1
			&& source.excluded("UNSELECTED_REGION")==1, "Excluded sibling rows retain their explicit reasons");
		check(source.get(0).id.equals("chromosome") && source.get(1).id.equals("mixed")
			&& source.get(2).id.equals("rna") && source.get(3).id.equals("rna_exon"), "Retained rows stay in source order");
		check(source.get(0).selection()==NcrnaHistoricalEvidence.ANCESTOR
			&& source.get(3).selection()==NcrnaHistoricalEvidence.DESCENDANT,
			"Context ancestors must not seed another downward walk");
	}

	/** A shared child retains every Parent, including outside context, but not that context's other children. */
	private static void multipleParents(final File root, final Observation observed){
		final Source source=read(root, "parents", row("ncRNA", "ID=rna")
			+row("exon", "ID=shared;Parent=rna,outside")
			+row("mRNA", "ID=outside")
			+row("CDS", "ID=sibling;Parent=outside"), observed);
		check(source.size()==3 && source.excludedRows()==1 && source.missingParents()==0,
			"Outside parent is context, its sibling CDS is not");
		check(source.get(1).parents.length==2 && source.get(1).parents[0].equals("rna")
			&& source.get(1).parents[1].equals("outside"), "Every original Parent relationship survives");
	}

	/** Same-ID pieces keep their original ranges; an incompatible reuse is flagged, not silently joined. */
	private static void multipart(final File root, final Observation observed){
		final String first=feature("seq_tid_7", "exon", "1", "3", "+", "ID=part;Parent=rna");
		final String second=feature("seq_tid_7", "exon", "8", "10", "+", "ID=part;Parent=rna");
		final Source source=read(root, "multipart", first+second+row("ncRNA", "ID=rna"), observed);
		check(source.size()==3 && source.missingParents()==0 && source.get(0).stop==3 && source.get(1).start==8,
			"Child-before-parent and multipart rows must not turn into a synthesized span");
		check(!source.get(0).incompatibleId() && !source.get(1).incompatibleId(), "Compatible pieces are not ID conflicts");
		check(new String(source.get(0).raw, StandardCharsets.UTF_8).equals(first.trim()), "Source bytes survive unchanged");
		final Source conflict=read(root, "conflict", row("ncRNA", "ID=reused")
			+feature("other", "CDS", "2", "6", "-", "ID=reused"), observed);
		check(conflict.size()==2 && conflict.get(0).incompatibleId() && conflict.get(1).incompatibleId(),
			"Incompatible target/type/strand reuse remains explicit on every affected row");
		final Source separate=read(root, "separate", row("gene", "ID=rna"), observed);
		check(separate.size()==0 && separate.excludedRows()==1, "ID namespaces do not leak between GFF sources");
	}

	/** Finite walks preserve dangling Parent IDs and actual back-edge witnesses, including self-loops. */
	private static void cycles(final File root, final Observation observed){
		final Source source=read(root, "cycles", row("ncRNA", "ID=a;Parent=b,missing")
			+row("gene", "ID=b;Parent=a")+row("ncRNA", "ID=self;Parent=self"), observed);
		check(source.size()==3 && source.missingParents()==1 && source.missingParent(0).equals("missing"),
			"Dangling selected relationships stay explicit and do not suppress evidence");
		check(source.cycles()==2 && source.cycle(0).child.equals("b") && source.cycle(0).parent.equals("a")
			&& source.cycle(0).line==2 && source.cycle(1).child.equals("self"), "Witnesses identify real source Parent edges");
	}

	/*--------------------------------------------------------------*/
	/*----------------      Parsing and states      ----------------*/
	/*--------------------------------------------------------------*/

	private static void attributes(final File root, final Observation observed){
		final Source source=read(root, "attributes", row("gene", "ID=g%2Cx%3Ay+z")
			+row("ncRNA", "ID=n;Parent=g%2Cx%3Ay+z")
			+row("regulatory", "ID=aptamer;Dbxref=GeneID:7,RFAM%3ARF00504,RFAM:RF99999")
			+row("mystery", "ID=unknown")
			+row("gene", "ID=hint;ncRNA_class=other")
			+row("gene", "ID=biotype;gene_biotype=ncRNA")
			+row("gene", "ID=coding;gene_biotype=protein_coding")
			+row("mRNA", "ID=coding2;gbkey=mRNA"), observed);
		check(source.size()==6 && source.excludedRows()==2 && source.missingParents()==0,
			"Unknown, aptamer and explicit RNA hints remain reviewable without selecting coding-only hints");
		check(source.get(0).id.equals("g,x:y+z") && source.get(1).parents.length==1,
			"Percent-encoded commas are identifier bytes, while plus stays literal");
		final Feature aptamer=source.get(2);
		check(aptamer.rfam.length==2 && aptamer.rfam[0].equals("RF00504") && aptamer.rfam[1].equals("RF99999"),
			"RFAM evidence is retained even outside today's searched model set");
		check(source.get(3).selection()==NcrnaHistoricalEvidence.DIRECT_UNKNOWN,
			"An unknown source type is not reclassified as proven non-coding");
	}

	private static void coordinates(final File root, final Observation observed){
		final Source source=read(root, "coordinates",
			feature("seq_tid_7", "ncRNA", "1", "20", "+", "ID=valid")
			+feature("seq_tid_7", "ncRNA", "1", "21", "-", "ID=outside")
			+feature("seq", "ncRNA", "1", "20", "+", "ID=unresolved")
			+feature("seq_tid_7", "ncRNA", "bad", "20", ".", "ID=bad_start")
			+feature("seq_tid_7", "ncRNA", "8", "3", "?", "ID=reversed")
			+feature("seq_tid_7", "ncRNA", "1", "99999999999999999999", "+", "ID=overflow"), observed);
		check(source.size()==6 && source.get(0).range.equals("VALID_INTERVAL") && source.get(1).range.equals("OUT_OF_RANGE")
			&& source.get(2).range.equals("UNRESOLVED_LITERAL_TARGET"), "Literal resolution and interval validity remain separate");
		check(source.get(3).range.equals("INVALID_COORDINATES") && source.get(4).range.equals("INVALID_COORDINATES")
			&& source.get(5).range.equals("INVALID_COORDINATES"), "Bad, reversed and overflowing coordinates remain explicit evidence");
		check(source.get(1).stop==21 && source.get(3).strand=='.' && source.get(4).strand=='?', "No range clipping or strand invention");
	}

	/** Encoded values are interpreted only after decoding; unknown claims never disappear. */
	private static void biotypes(final File root, final Observation observed){
		final Source source=read(root, "biotypes", row("gene", "ID=rna;gene_biotype=nc%52NA")
			+row("gene", "ID=mrna;gene_biotype=%6DRNA")
			+row("gene", "ID=noncoding;gene_biotype=non_coding")
			+row("exon", "ID=unknown;transcript_biotype=unrecognized")
			+row("gene", "ID=keyrna;gbkey=nc%52NA")
			+row("gene", "ID=keymrna;gbkey=%6DRNA")
			+row("gene", "ID=ordinary;gbkey=Gene;gene_biotype=protein_coding")
			+row("gene", "ID=keygene;gbkey=Gene")
			+row("ncRNA", "."), observed);
		check(source.size()==5 && source.excludedRows()==4, "Decoded mRNA and protein claims stay excluded; unknown biotypes remain reviewable");
		check(source.get(0).id.equals("rna") && source.get(0).selection()==NcrnaHistoricalEvidence.DIRECT_RNA,
			"Encoded ncRNA biotype selects RNA evidence");
		check(source.get(1).id.equals("noncoding") && source.get(1).selection()==NcrnaHistoricalEvidence.DIRECT_UNKNOWN
			&& source.get(2).id.equals("unknown"), "Non-coding and unfamiliar biotypes retain their original claims without family invention");
		check(source.get(3).id.equals("keyrna") && source.get(4).id==null, "Decoded RNA gbkey and missing-attribute sentinel are supported");
		check(new String(source.get(0).raw, StandardCharsets.UTF_8).contains("gene_biotype=nc%52NA"), "Decoding never rewrites original attributes");
	}

	private static void sourceStates(final File root, final Observation observed){
		check(NcrnaHistoricalEvidence.read("NA", observed).state().equals("MISSING_GFF"), "Missing claim is not an empty annotation");
		check(NcrnaHistoricalEvidence.read("relative.gff", observed).state().equals("UNRESOLVED_GFF_PATH"), "Relative path is unresolved");
		check(NcrnaHistoricalEvidence.read(new File(root, "absent.gff").getPath(), observed).state().equals("MISSING_GFF"), "Absent file is missing");
		check(NcrnaHistoricalEvidence.read(root.getPath(), observed).state().equals("NOT_REGULAR_GFF"), "Directory cannot masquerade as a GFF");
		check(read(root, "empty", "# comment\n\n", observed).state().equals("EMPTY_GFF"), "Present empty source is distinct from missing source");
		final Source inline=read(root, "inline", "# comment\n\n"+row("ncRNA", "ID=rna")+"##FASTA\n>seq\nACTG\n", observed);
		check(inline.size()==1 && inline.get(0).line==3 && inline.state().equals("PARSED_GFF"), "Original line numbers survive comments and embedded FASTA");
	}

	private static void malformed(final File root, final Observation observed){
		reject(root, "fields", "seq\tRNA\n", observed, "nine GFF fields");
		reject(root, "missing_equal", row("ncRNA", "ID=rna;Parent"), observed, "missing_equal.gff:1: Malformed GFF attribute");
		reject(root, "duplicate_id", row("ncRNA", "ID=a;ID=b"), observed, "Duplicate GFF ID");
		reject(root, "empty_parent", row("ncRNA", "Parent=a,,b"), observed, "Empty GFF Parent");
		reject(root, "escape", row("ncRNA", "ID=a%xx"), observed, "Invalid GFF percent escape");
		reject(root, "utf8", row("ncRNA", "ID=%FF"), observed, "invalid UTF-8");
		reject(root, "rfam", row("ncRNA", "Dbxref=RFAM:"), observed, "Empty RFAM");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/

	private static Source read(final File root, final String name, final String text, final Observation observed){
		return NcrnaHistoricalEvidence.read(write(root, name+".gff", text), observed);
	}
	private static String row(final String type, final String attributes){return feature("seq_tid_7", type, "1", "20", "+", attributes);}
	private static String feature(final String target, final String type, final String start, final String stop, final String strand, final String attributes){
		return target+"\tfixture\t"+type+'\t'+start+'\t'+stop+"\t.\t"+strand+"\t.\t"+attributes+'\n';
	}
	private static String write(final File root, final String name, final String text){
		final String path=new File(root, name).getPath();
		try(Output output=new Output(path)){output.writer.print(text);}
		return path;
	}
	private static void reject(final File root, final String name, final String text, final Observation observed, final String fragment){
		try{read(root, name, text, observed);}
		catch(IllegalArgumentException expected){
			check(expected.getMessage()!=null && expected.getMessage().contains(fragment), "Unexpected failure: "+expected);
			return;
		}
		throw new AssertionError("Expected failure containing: "+fragment);
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
}
