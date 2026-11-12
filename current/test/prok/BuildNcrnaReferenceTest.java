package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.BuildNcrnaReference;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaReferenceEvidence;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;

/** Exact end-to-end evidence expectations; no biology or performance is inferred.
 * @author Brian Bushnell, Raiden
 */
public final class BuildNcrnaReferenceTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, BuildNcrnaReferenceTest.class, false);
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
			check(root.mkdir(), "Fixture root must be new: "+root);
			new BuildNcrnaReferenceTest(root).run();
			pp.outstream.println("BUILD_NCRNA_REFERENCE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private BuildNcrnaReferenceTest(final File root_){root=root_;}

	private void run(){
		prepare();
		invoke("one", 1, Integer.MAX_VALUE, collector, memberships);
		invoke("two", 2, Integer.MAX_VALUE, collector, memberships);
		for(final String suffix:SUFFIXES){
			final String[] one=lines(path("one")+suffix), two=lines(path("two")+suffix);
			Arrays.sort(one);
			Arrays.sort(two);
			check(Arrays.equals(one, two), "Worker count changed keyed evidence in "+suffix);
		}
		final String[] summary=lines(path("two")+".summary.tsv");
		check(summary.length==2 && summary[1].equals("2\t2\t40\t2\t3\t2\t1\t10\t8\t2\t9\t3\t3\t6\t4\tEVIDENCE_COMPLETE_NOT_ADJUDICATED"),
			"Two genomes, three raw events and eight historical rows are distinct conserved units: "+Arrays.toString(summary));
		final String[] historical=lines(path("two")+".historical.tsv.gz");
		check(historical.length==9, "Retained context and component rows must remain present");
		check(find(historical, "\tregion\t").contains("\tSTRUCTURAL_CONTEXT\t"), "Ancestor region cannot seed gene-overlap support");
		check(find(historical, "\ta\t").endsWith("\t2\t1\t0\t3"), "Every alternative survives; opposite strand never wins support");
		check(find(historical, "\texon\t").endsWith("\t1\t1\t0\t3"), "An exon is a component row, not another inferred gene");
		check(find(historical, "\tunknown_strand\t").contains("\tUNASSESSABLE_STRAND\t"), "Unknown strand is not zero support");
		check(find(historical, "\toutside\t").contains("\tOUT_OF_RANGE\t"), "Bad source coordinates remain explicit");
		check(find(historical, "\tforeign\t").contains("\tUNRESOLVED_LITERAL_TARGET\t"), "No taxonomy-looking alias guessing");
		final String[] sources=lines(path("two")+".sources.tsv.gz");
		check(sources.length==3 && find(sources, "\tPARSED_GFF\t").contains("\t2\t10\t8\t2\t"),
			"Two memberships of one GFF must not duplicate its source rows");
		check(find(sources, "\tMISSING_GFF\t").startsWith("2\tNA\t"), "Missing second-input GFF is not an empty negative source");
		final String[] coverage=lines(path("two")+".coverage.tsv");
		check(coverage.length==4 && contains(coverage, "RF00050\tfixture\t10\t10\t0\t2\t2\t2")
			&& contains(coverage, "RF00162\tfixture\t10\t7\t1\t1\t1\t2")
			&& contains(coverage, "RF00001\tfixture\t10\tNA\tNA\t0\t0\t2"),
			"Model-specific distributions retain different truncations and completed zero-hit models");
		final String[] links=lines(path("two")+".overlaps.tsv.gz");
		check(links.length==13, "Nine historical and three model-model overlaps must survive without pair duplication");
		for(final String link:links){
			if(link.contains("OPPOSITE_STRAND")){check(link.endsWith("\tfalse"), "Opposite-strand geometry cannot win support");}
		}
		final String[] relations=lines(path("two")+".relations.tsv.gz");
		check(relations.length==9, "Five Parent edges, two RFAM claims and one directed cycle witness must survive");
		check(contains(relations, "1\t"+gff+"\t2\tPARENT\ta\troot\tPRESENT_ID"), "A selected RNA retains its actual Parent ID");
		check(contains(relations, "1\t"+gff+"\t2\tRFAM\ta\tRF00050\tSOURCE_CLAIM"), "RFAM is retained as a source claim");
		check(contains(relations, "1\t"+gff+"\t9\tRFAM\tc1\tRF99999\tSOURCE_CLAIM"), "Unrecognized RFAM claims are not discarded");
		check(contains(relations, "1\t"+gff+"\t9\tPARENT\tc1\tmissing\tDANGLING"), "Dangling relationships retain their referring source line");
		check(contains(relations, "1\t"+gff+"\t10\tCYCLE\tc2\tc1\tBACK_EDGE_WITNESS"), "Cycle witness preserves its actual directed edge");
		invoke("sample", 1, 1, collector, memberships);
		check(lines(path("sample")+".summary.tsv")[1].endsWith("\tSAMPLE_EVIDENCE_NOT_ADJUDICATED"), "A bounded prefix must remain labeled a sample");
		final String badCollector=write("bad.collector.tsv", collectorText.replace("2\t"+fasta2, "2\t"+fasta1));
		expect(()->invoke("bad_join", 1, 1, badCollector, memberships), "exact census FASTA");
		check(!new File(path("bad_join")+".summary.tsv").exists(), "Invalid source binding cannot publish completion");
		check(new File(root, "stale").mkdir(), "Stale-directory fixture must be new");
		copy("result1/hits.tbl", "stale/hits.tbl");
		copy("result1/completion.tsv", "stale/completion.tsv");
		write("stale/input.tsv", "input_id\tfasta\tnucleotide_bases\tFASTA_records\n1\t"+fasta2+"\t20\t1\n");
		final String stale=write("stale.collector.tsv", collectorText.replace(path("result1"), path("stale")));
		expect(()->invoke("stale_input", 1, Integer.MAX_VALUE, stale, memberships), "Native search input binding disagrees");
		check(!new File(path("stale_input")+".summary.tsv").exists(), "Matching IDs, dimensions and hits do not excuse a different original FASTA");
		final String malformed=write("malformed.gff", "seq\tfixture\tncRNA\t1\t10\t.\t+\t.\tID=rna;Parent\n");
		final String badMembers=write("bad.members.tsv", membershipText.replace(gff, malformed));
		expect(()->invoke("worker_failure", 2, Integer.MAX_VALUE, collector, badMembers), "Malformed GFF attribute");
		check(!new File(path("worker_failure")+".summary.tsv").exists(), "Worker parse failure must not yield a success summary");
		indexAndEncoding();
		new NcrnaReferenceProfilesTest(new File(root, "profiles")).run();
		new NcrnaReferenceGeometryTest(new File(root, "geometry")).run();
		new NcrnaReferenceAttributesTest(new File(root, "attributes")).run();
		new NcrnaReferenceStrataTest(new File(root, "strata")).run();
		new NcrnaReferencePriorTest(new File(root, "prior")).run();
		NcrnaReferenceMatcherTest.run();
		NcrnaReferenceGffTest.run(new File(root, "reference_gff"));
		NcrnaReferenceCallsTest.run(new File(root, "reference_calls"));
		new GradeNcrnaReferenceTest(new File(root, "reference_grade")).run();
		new NcrnaSubtypeGradeTest(new File(root, "subtype_grade")).run();
		new WriteNcrnaReferenceTest(new File(root, "reference_write")).run();
		new CompareNcrnaScenariosTest(new File(root, "reference_scenarios")).run();
		new BindNcrnaCallsTest(new File(root, "call_binding")).run();
		new DraftNcrnaScenariosTest(new File(root, "draft_scenarios")).run();
		new SelectNcrnaScenariosTest(new File(root, "select_scenarios")).run();
		new NcrnaReferenceCostTest(new File(root, "cost")).run();
		new NcrnaReferencePanelsTest(new File(root, "panels")).run();
		new NcrnaEvaluationPanelTest(new File(root, "evaluation_panel")).run();
		new NcrnaReferenceAnchorsTest(new File(root, "anchors")).run();
	}

	private void indexAndEncoding(){
		Shared.setThreads(1);
		final NcrnaEvidenceCatalog catalog=NcrnaEvidenceCatalog.load(census, memberships);
		final NcrnaReferenceEvidence.Roster roster=NcrnaReferenceEvidence.loadModels(scope, assets, queries);
		final NcrnaReferenceEvidence.Evidence evidence=NcrnaReferenceEvidence.read(path("result1/hits.tbl"), path("result1/completion.tsv"), catalog.observe(catalog.get("1")), roster);
		final BuildNcrnaReference.HitIndex index=new BuildNcrnaReference.HitIndex(evidence);
		final IntList found=new IntList(1);
		index.overlaps("seq", 18, 19, found);
		check(found.size==1 && index.hits[found.get(0)].line==3, "A long earlier interval cannot be hidden by shorter later endpoints");
		index.overlaps("different", 1, 20, found);
		check(found.size==0, "Target namespaces remain separate");
		final ByteBuilder encoded=BuildNcrnaReference.field(new ByteBuilder(), "a%\t\r\nΩ");
		check(Arrays.equals(encoded.toBytes(), "a%25%09%0D%0AΩ".getBytes(StandardCharsets.UTF_8)), "TSV framing must preserve UTF-8 and escaped controls");
	}

	private void prepare(){
		check(new File(root, "result1").mkdir() && new File(root, "result2").mkdir(), "Fixture result directories must be new");
		fasta1=write("one.fa", ">seq\nAAAAAAAAAAAAAAAAAAAA\n");
		fasta2=write("two.fa", ">seq\nCCCCCCCCCCCCCCCCCCCC\n");
		gff=write("source.gff", gff("region", 1, 20, "+", "ID=root")
			+gff("ncRNA", 1, 10, "+", "ID=a;Parent=root;Dbxref=RFAM:RF00050")
			+gff("exon", 2, 3, "+", "ID=part;Parent=a")
			+gff("ncRNA", 11, 15, ".", "ID=unknown_strand")
			+gff("ncRNA", 18, 30, "+", "ID=outside")
			+"seq_tid_7\tfixture\tncRNA\t1\t5\t.\t+\t.\tID=foreign\n"
			+gff("gene", 1, 20, "+", "ID=coding;Parent=root")
			+gff("CDS", 1, 20, "+", "ID=protein;Parent=coding")
			+gff("ncRNA", 1, 2, "+", "ID=c1;Parent=c2,missing;Dbxref=RFAM:RF99999")
			+gff("gene", 1, 2, "+", "ID=c2;Parent=c1"));
		census=write("census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'
			+"1\t"+fasta1+'\t'+new File(fasta1).length()+"\t20\t1\n"
			+"2\t"+fasta2+'\t'+new File(fasta2).length()+"\t20\t1\n");
		membershipText=NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'
			+"first\t1\t1\t"+fasta1+"\tNA\tNA\tNA\tprokaryotic\t"+gff+"\tNA\n"
			+"first\t2\t2\t"+fasta2+"\tNA\tNA\tNA\tprokaryotic\tNA\tNA\n"
			+"second\t1\t1\t"+fasta1+"\tNA\tNA\tNA\tUNKNOWN\t"+gff+"\tNA\n";
		memberships=write("members.tsv", membershipText);
		scope=write("scope.tsv", NcrnaReferenceEvidence.SCOPE_HEADER+'\n'
			+"RF00050\tFMN\tprimary\taptamer_only\tFMN\tTARGET_CLASS\n"
			+"RF00162\tSAM-I\tprimary\taptamer_only\tSAM_I\tTARGET_CLASS\n"
			+"RF00001\t5S\tprimary\tcompartment_unresolved\t5S\tTARGET_CLASS\n");
		assets=write("assets.tsv", NcrnaReferenceEvidence.ASSET_HEADER+'\n'
			+"RF00050\tfixture\t10\t20\t00000000000000000001\t/fixture/fmn.cm\n"
			+"RF00162\tfixture\t10\t30\t00000000000000000002\t/fixture/sam1.cm\n"
			+"RF00001\tfixture\t10\t20\t00000000000000000003\t/fixture/5s.cm\n");
		queries=write("queries.tsv", NcrnaReferenceEvidence.QUERY_HEADER+"\nRF00050\tfmn\t10\nRF00162\tsam1\t10\nRF00001\t5s\t10\n");
		final String header="# Program: cmsearch\n# Version: 1.1.5 (fixture)\n";
		write("result1/hits.tbl", header
			+"seq - fmn RF00050 cm 1 10 1 20 + no 1 0.5 0 50 1e-8 ! fixture\n"
			+"seq - sam1 RF00162 cm 2 8 8 14 + 5' 1 0.5 0 60 1e-9 ! fixture\n"
			+"seq - fmn RF00050 cm 1 10 10 1 - no 1 0.5 0 70 1e-10 ! fixture\n# [ok]\n");
		write("result2/hits.tbl", header+"# [ok]\n");
		write("result1/input.tsv", "input_id\tfasta\tnucleotide_bases\tFASTA_records\n1\t"+fasta1+"\t20\t1\n");
		write("result2/input.tsv", "input_id\tfasta\tnucleotide_bases\tFASTA_records\n2\t"+fasta2+"\t20\t1\n");
		write("result1/completion.tsv", NcrnaReferenceEvidence.COMPLETION_HEADER+"\n"
			+"1\tRF00050\t1\t20\t40\t2\tCOMPLETE\n1\tRF00162\t1\t20\t40\t1\tCOMPLETE\n1\tRF00001\t1\t20\t40\t0\tCOMPLETE\n");
		write("result2/completion.tsv", NcrnaReferenceEvidence.COMPLETION_HEADER+"\n"
			+"2\tRF00050\t1\t20\t40\t0\tCOMPLETE\n2\tRF00162\t1\t20\t40\t0\tCOMPLETE\n2\tRF00001\t1\t20\t40\t0\tCOMPLETE\n");
		collectorText=BuildNcrnaReference.COLLECTOR_HEADER+'\n'
			+"1\t"+fasta1+'\t'+path("result1")+"\t3\t3\t0\t0\t0\n"
			+"2\t"+fasta2+'\t'+path("result2")+"\t3\t0\t0\t0\t0\n";
		collector=write("collector.tsv", collectorText);
	}

	private void invoke(final String name, final int threads, final int maximum, final String collection, final String members){
		BuildNcrnaReference.main(new String[]{"in="+collection, "census="+census, "memberships="+members,
			"scope="+scope, "assets="+assets, "queries="+queries, "out="+path(name), "t="+threads, "maxinputs="+maximum});
	}
	private static String gff(final String type, final int start, final int stop, final String strand, final String attributes){
		return "seq\tfixture\t"+type+'\t'+start+'\t'+stop+"\t.\t"+strand+"\t.\t"+attributes+'\n';
	}
	private String path(final String name){return new File(root, name).getPath();}
	/** Copies only tiny synthetic fixture data, never a program or real genome. */
	private void copy(final String source, final String target){
		try(Input in=new Input(path(source)); Output out=new Output(path(target))){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){out.writer.println(row);}
		}
	}
	private String write(final String name, final String text){
		final String file=path(name);
		try(Output out=new Output(file)){out.writer.print(text);}
		return file;
	}
	private static String[] lines(final String file){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(file)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				result.add(new String(row, StandardCharsets.UTF_8));
			}
		}
		return result.toArray(new String[result.size()]);
	}
	private static String find(final String[] rows, final String fragment){
		String found=null;
		for(final String row:rows){
			if(row.contains(fragment)){
				check(found==null, "Fixture fragment must identify one row: "+fragment);
				found=row;
			}
		}
		check(found!=null, "Missing fixture row: "+fragment);
		return found;
	}
	private static boolean contains(final String[] rows, final String expected){
		for(final String row:rows){if(row.equals(expected)){return true;}}
		return false;
	}
	private static void expect(final Runnable action, final String fragment){
		try{action.run();}
		catch(RuntimeException failure){
			for(Throwable cause=failure; cause!=null; cause=cause.getCause()){
				if(cause.getMessage()!=null && cause.getMessage().contains(fragment)){return;}
			}
			throw new AssertionError("Unexpected failure", failure);
		}
		throw new AssertionError("Expected failure containing: "+fragment);
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final File root;
	private String fasta1, fasta2, gff, census, memberships, scope, assets, queries, collector, collectorText, membershipText;
	private static final String[] SUFFIXES={".contigs.tsv.gz", ".hits.tsv.gz", ".sources.tsv.gz", ".historical.tsv.gz", ".overlaps.tsv.gz", ".relations.tsv.gz", ".coverage.tsv", ".summary.tsv"};
}
