package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import prok.BuildNcrnaReference;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaReferenceEvidence;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;

/** Current-source capture, exact historical-key binding and honest change semantics.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceAttributesTest{

	NcrnaReferenceAttributesTest(final File root_){root=root_;}

	void run(){
		check(root.mkdir(), "Attribute fixture root must be new");
		prepare();
		evidence("full", false, history, relations);
		invoke("full", "one", 1);
		invoke("full", "two", 2);
		final String[] one=lines("one.attributes.tsv.gz"), two=lines("two.attributes.tsv.gz");
		Arrays.sort(one);
		Arrays.sort(two);
		check(Arrays.equals(one, two) && one.length==8, "Both worker counts must retain all seven input/source-row keys");
		check(lines("two.attribute_summary.tsv")[1].equals("3\t3\t7\t7\t2\t4\t7\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tCURRENT_SOURCE_ATTRIBUTES_NOT_ADJUDICATED"),
			"One physical path serves two inputs: seven keyed rows from four selected lines in two paths");
		check(contains(one, "1\t"+encodedA+"\t4\t"+FRAMED_CHILD)
			&& contains(one, "2\t"+encodedA+"\t4\t"+FRAMED_CHILD), "Raw escapes, duplicate parents/RFAMs and both source memberships survive");
		check(contains(one, "1\t"+encodedA+"\t5\tID=child%253B1;Parent=root"), "Repeated/multipart IDs remain separate source lines");
		replaceSource(BASE.replace("pseudo=true", "pseudo=false"), "before_attribute_change");
		invoke("full", "current", 2);
		check(contains(lines("current.attributes.tsv.gz"), "1\t"+encodedA+"\t4\t"+FRAMED_CHILD.replace("pseudo=true", "pseudo=false")),
			"Previously uncaptured attributes may differ; output must report the current bytes honestly");
		check(lines("current.attribute_summary.tsv")[1].endsWith("\tCURRENT_SOURCE_ATTRIBUTES_NOT_ADJUDICATED"), "Identity matches must not certify historical attribute bytes");
		failure("inserted_line", "# inserted\n"+BASE, "Saved line no longer identifies a GFF feature");
		failure("missing_line", BASE.replace(EXON, ""), "ended before all saved retained lines");
		failure("changed_parent", BASE.replace("Parent=root,ghost,root", "Parent=root,other,root"), "Current Parent identity/order differs");
		failure("changed_rfam", BASE.replace("RFAM:RF00050,RFAM:RF00050", "RFAM:RF00050,RFAM:RF00051"), "Current RFAM identity/order differs");
		failure("changed_coordinate", BASE.replace("\tncRNA\t2\t19\t", "\tncRNA\t3\t19\t"), "Current target/type/coordinates/strand/ID differ");
		replaceSource(BASE, "before_valid_sample");
		final String firstRfam="1\t"+encodedA+"\t4\tRFAM\tchild;1\tRF00050\tSOURCE_CLAIM\n";
		evidence("missing_relation", false, history, relations.replace(firstRfam, ""));
		expect(()->invoke("missing_relation", "missing_relation_out", 1), "Missing saved Parent/RFAM claims");
		check(!new File(path("missing_relation_out.attributes.tsv.gz")).exists(), "Incomplete saved relationships fail before writing attributes");
		evidence("duplicate_key", false, history.replace("1\t"+encodedA+"\t4\t", "1\t"+encodedA+"\t3\t"), relations);
		expect(()->invoke("duplicate_key", "duplicate_key_out", 1), "Duplicate input/GFF/line");
		evidence("sample", true, history.substring(0, history.lastIndexOf("3\t"+encodedB)), relations);
		invoke("sample", "sample_out", 2);
		check(lines("sample_out.attribute_summary.tsv")[1].equals("2\t3\t6\t6\t1\t3\t6\tSAMPLE_EVIDENCE_NOT_ADJUDICATED\tCURRENT_SOURCE_ATTRIBUTES_NOT_ADJUDICATED"),
			"A sampled source join retains both observed and catalog denominators");
		check(new File(a).renameTo(new File(root, "missing_source.preimage.gff")), "Preserve the source before testing an unavailable path");
		// ReadWrite.exceptionKill terminates the JVM for absent files. The saved
		// verification producer invokes this full population in a separate process.
		check(new File(b).renameTo(new File(root, "zero_source.preimage.gff")), "The zero-row case must have no readable source payloads");
		zeroRows();
	}

	private void zeroRows(){
		write("zero.summary.tsv", BuildNcrnaReference.SUMMARY_HEADER+"\n3\t3\t300\t3\t0\t3\t3\t0\t0\t0\t0\t0\t1\t3\t3\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\n");
		write("zero.coverage.tsv", BuildNcrnaReference.COVERAGE_HEADER+"\nRF00050\tfixture\t10\tNA\tNA\t0\t0\t3\n");
		final StringBuilder sources=new StringBuilder(BuildNcrnaReference.SOURCE_HEADER).append('\n');
		for(int i=1; i<=3; i++){
			sources.append(i).append('\t').append(i==3 ? encodedB : encodedA).append("\tMISSING_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n");
		}
		write("zero.sources.tsv.gz", sources.toString());
		write("zero.historical.tsv.gz", BuildNcrnaReference.HISTORICAL_HEADER+'\n');
		write("zero.relations.tsv.gz", BuildNcrnaReference.RELATION_HEADER+'\n');
		invoke("zero", "zero_out", 2);
		check(lines("zero_out.attribute_summary.tsv")[1].equals("3\t3\t0\t0\t0\t0\t0\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tCURRENT_SOURCE_ATTRIBUTES_NOT_ADJUDICATED"),
			"An empty retained-row population needs no original payload reads and does not invent annotation absence");
	}

	private void prepare(){
		a=path("percent%source.gff");
		b=path("second.gff");
		encodedA=a.replace("%", "%25");
		encodedB=b.replace("%", "%25");
		writePath(a, BASE);
		writePath(b, "other\tfixture\ttRNA\t1\t10\t.\t-\t.\tID=second;partial=true;exception=RNA%20editing\n");
		final StringBuilder census=new StringBuilder(NcrnaEvidenceCatalog.CENSUS_HEADER).append('\n');
		final StringBuilder members=new StringBuilder(NcrnaEvidenceCatalog.MEMBERSHIP_HEADER).append('\n');
		for(int i=1; i<=3; i++){
			final String fasta=path("absent"+i+".fa");
			check(!new File(fasta).exists(), "Attribute joins must not require FASTA payloads");
			census.append(i).append('\t').append(fasta).append("\t110\t100\t1\n");
			members.append("fixture\t").append(i).append('\t').append(i).append('\t').append(fasta)
				.append("\tNA\tNA\tNA\tUNKNOWN\t").append(i==3 ? b : a).append("\tNA\n");
		}
		write("census.tsv", census.toString());
		write("members.tsv", members.toString());
		write("scope.tsv", NcrnaReferenceEvidence.SCOPE_HEADER+"\nRF00050\tFMN\tprimary\tfixture\tFMN\tTARGET_CLASS\n");
		write("assets.tsv", NcrnaReferenceEvidence.ASSET_HEADER+"\nRF00050\tfixture\t10\t20\t00000000000000000001\t/missing/fmn.cm\n");
		write("queries.tsv", NcrnaReferenceEvidence.QUERY_HEADER+"\nRF00050\tfmn\t10\n");
		final StringBuilder rows=new StringBuilder(), links=new StringBuilder();
		for(int i=1; i<=2; i++){
			rows.append(i).append('\t').append(encodedA).append("\t3\tseq;Ω\tgene\t1\t20\t+\troot\t17\tVALID_INTERVAL\tfalse\tASSESSABLE\t0\t0\t0\t0\t0\t-1\n");
			rows.append(i).append('\t').append(encodedA).append("\t4\tseq;Ω\tncRNA\t2\t19\t+\tchild;1\t3\tVALID_INTERVAL\ttrue\tASSESSABLE\t3\t2\t0\t0\t0\t-1\n");
			rows.append(i).append('\t').append(encodedA).append("\t5\tseq;Ω\texon\t3\t4\t+\tchild;1\t8\tVALID_INTERVAL\ttrue\tASSESSABLE\t1\t0\t0\t0\t0\t-1\n");
			links.append(i).append('\t').append(encodedA).append("\t4\tPARENT\tchild;1\troot\tPRESENT_ID\n");
			links.append(i).append('\t').append(encodedA).append("\t4\tPARENT\tchild;1\tghost\tDANGLING\n");
			links.append(i).append('\t').append(encodedA).append("\t4\tPARENT\tchild;1\troot\tPRESENT_ID\n");
			links.append(i).append('\t').append(encodedA).append("\t4\tRFAM\tchild;1\tRF00050\tSOURCE_CLAIM\n");
			links.append(i).append('\t').append(encodedA).append("\t4\tRFAM\tchild;1\tRF00050\tSOURCE_CLAIM\n");
			links.append(i).append('\t').append(encodedA).append("\t5\tPARENT\tchild;1\troot\tPRESENT_ID\n");
		}
		rows.append("3\t").append(encodedB).append("\t1\tother\ttRNA\t1\t10\t-\tsecond\t1\tVALID_INTERVAL\tfalse\tASSESSABLE\t0\t0\t0\t0\t0\t-1\n");
		history=rows.toString();
		relations=links.toString();
	}

	private void evidence(final String name, final boolean sample, final String rows, final String links){
		final int n=sample ? 2 : 3, retained=sample ? 6 : 7;
		write(name+".summary.tsv", BuildNcrnaReference.SUMMARY_HEADER+'\n'+n+"\t3\t"+(n*100)+'\t'+n+"\t0\t"+n+"\t0\t"+retained+'\t'+retained
			+"\t0\t0\t0\t1\t"+n+'\t'+n+'\t'+(sample ? "SAMPLE_EVIDENCE_NOT_ADJUDICATED" : "EVIDENCE_COMPLETE_NOT_ADJUDICATED")+'\n');
		write(name+".coverage.tsv", BuildNcrnaReference.COVERAGE_HEADER+"\nRF00050\tfixture\t10\tNA\tNA\t0\t0\t"+n+'\n');
		final StringBuilder sources=new StringBuilder(BuildNcrnaReference.SOURCE_HEADER).append('\n');
		for(int i=1; i<=n; i++){
			final int count=i==3 ? 1 : 3;
			sources.append(i).append('\t').append(i==3 ? encodedB : encodedA).append("\tPARSED_GFF\t1\t").append(count).append('\t').append(count)
				.append("\t0\t0\t0\t0\t").append(i==3 ? 0 : 1).append("\t0\n");
		}
		write(name+".sources.tsv.gz", sources.toString());
		write(name+".historical.tsv.gz", BuildNcrnaReference.HISTORICAL_HEADER+'\n'+rows);
		write(name+".relations.tsv.gz", BuildNcrnaReference.RELATION_HEADER+'\n'+links);
	}
	private void failure(final String name, final String current, final String reason){
		replaceSource(current, "before_"+name);
		expect(()->invoke("full", name, 2), reason);
		check(!new File(path(name+".attribute_summary.tsv")).exists(), "A partial/failed source join must not publish completion");
	}
	/** Every overwritten synthetic source keeps its previous bytes under a unique name. */
	private void replaceSource(final String current, final String saved){
		final File previous=new File(path(saved+".gff"));
		check(!previous.exists() && new File(a).renameTo(previous), "Preserve each synthetic source preimage before replacement");
		writePath(a, current);
	}
	private void invoke(final String input, final String output, final int threads){
		BuildNcrnaReference.main(new String[]{"mode=attributes", "in="+path(input), "out="+path(output), "census="+path("census.tsv"),
			"memberships="+path("members.tsv"), "scope="+path("scope.tsv"), "assets="+path("assets.tsv"), "queries="+path("queries.tsv"),
			"t="+threads, "bgzip=f", "zipthreads=1"});
	}
	private void write(final String name, final String text){writePath(path(name), text);}
	private static void writePath(final String file, final String text){try(Output out=new Output(file)){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}}
	private String[] lines(final String name){
		final ArrayList<String> rows=new ArrayList<String>();
		try(Input in=new Input(path(name))){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){rows.add(new String(row, StandardCharsets.UTF_8));}
		}
		return rows.toArray(new String[rows.size()]);
	}
	private static boolean contains(final String[] rows, final String wanted){
		for(final String row:rows){if(row.equals(wanted)){return true;}}
		return false;
	}
	private static void expect(final Runnable action, final String reason){
		try{action.run();}
		catch(RuntimeException expected){
			for(Throwable cause=expected; cause!=null; cause=cause.getCause()){
				if(cause.getMessage()!=null && cause.getMessage().contains(reason)){return;}
			}
			throw new AssertionError("Unexpected fixture failure", expected);
		}
		throw new AssertionError("Expected failure containing: "+reason);
	}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final File root;
	private String a, b, encodedA, encodedB, history, relations;
	private static final String CHILD="ID=child%3B1;Parent=root,ghost,root;Dbxref=RFAM:RF00050,RFAM:RF00050;gene_biotype=pseudogene;pseudo=true;Note=a%3Bb";
	private static final String FRAMED_CHILD="ID=child%253B1;Parent=root,ghost,root;Dbxref=RFAM:RF00050,RFAM:RF00050;gene_biotype=pseudogene;pseudo=true;Note=a%253Bb";
	private static final String EXON="seq%3BΩ\tfixture\texon\t3\t4\t.\t+\t.\tID=child%3B1;Parent=root\n";
	private static final String BASE="# fixture\n\nseq%3BΩ\tfixture\tgene\t1\t20\t.\t+\t.\tID=root;gene_biotype=ncRNA\n"
		+"seq%3BΩ\tfixture\tncRNA\t2\t19\t.\t+\t.\t"+CHILD+'\n'+EXON+"##FASTA\n>seq\nAAAA\n";
}
