package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import prok.BuildNcrnaReference;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaReferenceEvidence;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;

/** Literal source-origin and nonexclusive-count oracles without genomic payloads.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceStrataTest{

	NcrnaReferenceStrataTest(final File root_){root=root_;}

	void run(){
		check(root.mkdir(), "Strata fixture directory must be new");
		prepare();
		invoke("evidence", "capture", "result");
		check(lines("result.strata_summary.tsv")[1].equals("2\t2\t4\t3\t6\t12\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tSTRATA_COMPLETE_NOT_ADJUDICATED"),
			"Four events, three groups and six source rows are separate conserved populations");
		final String[] events=lines("result.event_strata.tsv.gz");
		check(events.length==5, "All four original events must survive");
		final String a=find(events, "1\t/missing/hits1.tbl\t1\t");
		for(final String flag:new String[]{"SOURCE_SAME_ACCESSION", "SOURCE_SAME_FAMILY_OTHER_ACCESSION", "SOURCE_OTHER_KNOWN_FAMILY", "SOURCE_UNKNOWN_RFAM",
			"SOURCE_PSEUDO_FALSE", "SOURCE_PSEUDOGENE_BIOTYPE", "OVERLAPPING_OTHER_FAMILY_MODEL"}){
			check(a.contains(flag), "Retain independent source/model evidence axis: "+flag);
		}
		check(!a.contains("SOURCE_PSEUDO_TRUE") && !a.contains("SOURCE_EXCEPTION_ATTRIBUTE") && !a.contains("SOURCE_PARTIAL_TRUE"),
			"Coding/context, component and unassessable source claims must not transfer to an RNA event");
		check(a.contains("\t4\t0\t2\t1\t1\t0\t"), "Four same-strand rows retain RNA/component/context origin counts");
		final String opposite=find(events, "1\t/missing/hits1.tbl\t3\t");
		check(opposite.contains("OPPOSITE_STRAND_SOURCE_OVERLAP") && opposite.contains("NO_ASSESSABLE_SAME_STRAND_SOURCE_OVERLAP")
			&& !opposite.contains("SOURCE_PSEUDO_FALSE"), "Opposite-strand overlaps are not same-strand RNA claims");
		final String c=find(events, "1\t/missing/hits1.tbl\t2\t");
		check(c.contains("FULL_MODEL_SPAN") && c.contains("NATIVE_TRUNCATED"), "Full model span and native truncation must be independent axes");
		final String other=find(events, "2\t/missing/hits2.tbl\t1\t");
		check(other.contains("NO_ASSESSABLE_SAME_STRAND_SOURCE_OVERLAP") && !other.contains("SOURCE_SAME_ACCESSION"),
			"Equal native-line/contig/coordinates in another genome cannot inherit source claims");
		final String[] counts=lines("result.strata_counts.tsv");
		check(find(counts, "ALL_QUERIED\tALL\t").startsWith("ALL_QUERIED\tALL\t4\t3\t"), "Union groups counted once");
		check(find(counts, "RF00050\tALL\t").startsWith("RF00050\tALL\t3\t3\t"), "One accession can occur in three groups");
		check(find(counts, "RF00001\tALL\t").startsWith("RF00001\tALL\t1\t1\t"), "A shared group appears under both accessions");
		check(find(counts, "RF00051\tALL\t").startsWith("RF00051\tALL\t0\t0\t"), "Queried zero-hit models remain explicit");
		check(find(counts, "RF00050\tSOURCE_PSEUDO_FALSE\t").startsWith("RF00050\tSOURCE_PSEUDO_FALSE\t1\t1\t"),
			"Two false-claim source rows must not count the same event twice");
		final String[] joint=lines("result.joint_strata.tsv");
		check(find(joint, "RF00050\t100_percent\t0\t20_to_below50\tNO_AVAILABLE_GFF\t")
			.endsWith("\t1\t1\tDESCRIPTIVE_BINS_NONADDITIVE_GROUPS"), "Joint bins retain the unannotated full-model event");
		check(find(joint, "RF00050\t50_to_below80_percent\t0\t20_to_below50\tAVAILABLE_GFF_NO_ASSESSABLE_SAME_STRAND_OVERLAP\t")
			.contains("\tfalse\tfalse\tfalse\ttrue\tfalse\t1\t1\t"), "Opposite models and uncorroborated partial span remain separate axes");
		final String[] sources=lines("result.source_strata.tsv.gz");
		check(sources.length==7 && find(sources, "1\t/missing/source.gff\t4\t").contains("UNASSESSABLE_RANGE_OR_SOURCE"),
			"Unknown strand remains a retained unassessable source row");
		check(find(sources, "1\t/missing/source.gff\t5\t").contains("STRUCTURAL_SOURCE_CONTEXT"), "Ancestor-only context never becomes a gene");
		final String[] links=lines("result.source_links.tsv.gz");
		check(links.length==13 && find(links, "1\t/missing/source.gff\t2\t/missing/hits1.tbl\t1\t").contains("SOURCE_PSEUDO_TRUE"),
			"The context pseudo claim remains visible with its exact source/event association");
		copyEvidence("regulatory", overlaps, historyText.replace("\tgene\t", "\tregulatory_region\t"));
		invoke("regulatory", "capture", "regulatory_result");
		final String regulatory=find(lines("regulatory_result.source_links.tsv.gz"), "1\t/missing/source.gff\t2\t/missing/hits1.tbl\t1\t");
		check(regulatory.contains("\tregulatory_region\tOTHER_ROW\t") && regulatory.contains("SOURCE_SAME_FAMILY_OTHER_ACCESSION"),
			"Regulatory RFAM evidence remains associated with the event, including its literal feature type");
		copyEvidence("invalid_hit", overlaps, historyText, hits.replace(hit(2, 1, "RF00050", 10, 19, '+', 1, 10, 0),
			hit(2, 1, "RF00050", 90, 110, '+', 1, 10, 0).replace("VALID_INTERVAL", "OUT_OF_RANGE")));
		invoke("invalid_hit", "capture", "invalid_hit_result");
		final String invalid=find(lines("invalid_hit_result.event_strata.tsv.gz"), "2\t/missing/hits2.tbl\t1\t");
		check(invalid.contains("UNASSESSABLE_RANGE_OR_SOURCE") && !invalid.contains("NO_ASSESSABLE_SAME_STRAND_SOURCE_OVERLAP"),
			"An invalid interval is unassessable, not evidence of absent source support");
		final String longCoordinates=historyText.replace("\t4\tΩ%25seq\tncRNA\t10\t19\t?\tid4\t1\tVALID_INTERVAL\tfalse\tUNASSESSABLE_STRAND",
			"\t4\tΩ%25seq\tncRNA\t-9223372036854775808\t5000000000\t?\tid4\t1\tOUT_OF_RANGE\tfalse\tOUT_OF_RANGE");
		check(!longCoordinates.equals(historyText), "The coordinate fixture must replace the intended literal source row");
		copyEvidence("long_source", overlaps, longCoordinates);
		invoke("long_source", "capture", "long_source_result");
		check(find(lines("long_source_result.source_strata.tsv.gz"), "1\t/missing/source.gff\t4\t")
			.contains("\tncRNA\t-9223372036854775808\t5000000000\t?\t"),
			"Malformed-coordinate sentinel and oversized source coordinates must survive without 32-bit wrapping");
		capture("contradict", attributes.replace("gene_biotype=pseudogene;partial=false", "gene_biotype=pseudogene;partial=false;pseudo=true"));
		invoke("evidence", "contradict", "contradict_result");
		final String conflicting=find(lines("contradict_result.event_strata.tsv.gz"), "1\t/missing/hits1.tbl\t1\t");
		check(conflicting.contains("SOURCE_PSEUDO_TRUE") && conflicting.contains("SOURCE_PSEUDO_FALSE"), "Contradictory RNA-source claims remain visible together");
		capture("duplicate", attributes+attribute(1, "pseudo=true"));
		expect(()->invoke("evidence", "duplicate", "bad_duplicate"), "Repeated captured");
		capture("missing", attributes.substring(0, attributes.lastIndexOf("1\t/missing/source.gff\t6\t")));
		expect(()->invoke("evidence", "missing", "bad_missing"), "cover every historical key");
		capture("malformed", attributes.replace("pseudo=false", "pseudo=%25QZ"));
		expect(()->invoke("evidence", "malformed", "bad_escape"), "Invalid percent escape");
		copyEvidence("duplicate_link", overlaps.replace(link(1, 1, false), link(2, 1, false)));
		expect(()->invoke("duplicate_link", "capture", "bad_link"), "Repeated historical/native overlap");
		for(final String name:new String[]{"bad_duplicate", "bad_missing", "bad_escape", "bad_link"}){
			check(!new File(path(name+".strata_summary.tsv")).exists(), "Failed joins cannot publish completion: "+name);
		}
		sampleAndZero();
		jointCounts();
		new NcrnaReferenceSourceAnchorsTest(root).run();
	}

	private void jointCounts(){
		final String duplicateModelHits=hits.replace(hit(1, 2, "RF00001", 15, 24, '+', 1, 10, 1), hit(1, 2, "RF00050", 15, 24, '+', 1, 10, 0));
		final String duplicateCoverage="RF00050\tfixture\t10\t10\t0\t3\t4\t2\nRF00050\tfixture\t10\t7\t0\t1\t4\t2\n"
			+"RF00051\tfixture\t10\tNA\tNA\t0\t0\t2\nRF00001\tfixture\t10\tNA\tNA\t0\t0\t2\n";
		copyEvidence("shared_bin", overlaps, historyText, duplicateModelHits, duplicateCoverage);
		invoke("shared_bin", "capture", "shared_bin_result");
		final String[] shared=lines("shared_bin_result.joint_strata.tsv");
		check(find(shared, "RF00050\t100_percent\t0\t20_to_below50\tRNA_ROW_OVERLAP\t")
			.contains("\ttrue\tfalse\tfalse\ttrue\tfalse\t2\t1\t"),
			"Two events in one joint bin count as two events and one group, including same-model overlap");
		check(find(shared, "ALL_QUERIED\t100_percent\t0\t20_to_below50\tRNA_ROW_OVERLAP\t")
			.endsWith("\t2\t1\tDESCRIPTIVE_BINS_NONADDITIVE_GROUPS"), "Union bin deduplication is separate from accession bin deduplication");
		final int[] scores={19, 20, 25, 30, 40, 70, 120};
		final String[] bins={"below0", "0_to_below5", "5_to_below10", "10_to_below20", "20_to_below50", "50_to_below100", "100_or_more"};
		for(int i=0; i<scores.length; i++){
			final String name="margin"+scores[i];
			final String original=hit(2, 1, "RF00050", 10, 19, '+', 1, 10, 0);
			final String changed=original.replace("\t50\t1e-8\ttrue", "\t"+scores[i]+"\t1e-8\t"+(scores[i]>=20));
			copyEvidence(name, overlaps, historyText, hits.replace(original, changed));
			invoke(name, "capture", name+"_result");
			check(find(lines(name+"_result.joint_strata.tsv"), "RF00050\t100_percent\t0\t"+bins[i]+"\tNO_AVAILABLE_GFF\t")
				.endsWith("\t1\t1\tDESCRIPTIVE_BINS_NONADDITIVE_GROUPS"), "Exact GA-margin interval boundary: "+scores[i]);
			check(find(lines(name+"_result.event_strata.tsv.gz"), "2\t/missing/hits2.tbl\t1\t")
				.contains("\t"+scores[i]+"\t20\t"+(scores[i]-20)+"\t"), "Exact score, GA and margin remain in event output");
		}
	}

	private void sampleAndZero(){
		write("sample.summary.tsv", BuildNcrnaReference.SUMMARY_HEADER+"\n1\t2\t100\t1\t3\t1\t0\t6\t6\t0\t12\t3\t3\t3\t1\tSAMPLE_EVIDENCE_NOT_ADJUDICATED\n");
		write("sample.coverage.tsv", BuildNcrnaReference.COVERAGE_HEADER+"\nRF00050\tfixture\t10\t10\t0\t1\t2\t1\nRF00050\tfixture\t10\t7\t0\t1\t2\t1\n"
			+"RF00051\tfixture\t10\tNA\tNA\t0\t0\t1\nRF00001\tfixture\t10\t10\t1\t1\t1\t1\n");
		write("sample.sources.tsv.gz", BuildNcrnaReference.SOURCE_HEADER+"\n1\t/missing/source.gff\tPARSED_GFF\t1\t6\t6\t0\t0\t0\t0\t0\t0\n");
		write("sample.hits.tsv.gz", BuildNcrnaReference.HIT_HEADER+'\n'+hits.replace(hit(2, 1, "RF00050", 10, 19, '+', 1, 10, 0), ""));
		write("sample.historical.tsv.gz", historyText);
		write("sample.relations.tsv.gz", relationsText);
		write("sample.overlaps.tsv.gz", BuildNcrnaReference.OVERLAP_HEADER+'\n'+overlaps);
		write("sample_capture.attribute_summary.tsv", ATTRIBUTE_HEADER+"\n1\t2\t6\t6\t1\t6\t6\tSAMPLE_EVIDENCE_NOT_ADJUDICATED\tCURRENT_SOURCE_ATTRIBUTES_NOT_ADJUDICATED\n");
		write("sample_capture.attributes.tsv.gz", "input_id\tgff\tline\traw_attributes\n"+attributes);
		invoke("sample", "sample_capture", "sample_result");
		check(lines("sample_result.strata_summary.tsv")[1].equals("1\t2\t3\t2\t6\t12\tSAMPLE_EVIDENCE_NOT_ADJUDICATED\tSTRATA_COMPLETE_NOT_ADJUDICATED"),
			"Sample strata retain observed and original catalog denominators");
		write("zero.summary.tsv", BuildNcrnaReference.SUMMARY_HEADER+"\n2\t2\t200\t2\t0\t2\t2\t0\t0\t0\t0\t0\t3\t6\t6\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\n");
		write("zero.coverage.tsv", BuildNcrnaReference.COVERAGE_HEADER+"\nRF00050\tfixture\t10\tNA\tNA\t0\t0\t2\n"
			+"RF00051\tfixture\t10\tNA\tNA\t0\t0\t2\nRF00001\tfixture\t10\tNA\tNA\t0\t0\t2\n");
		write("zero.sources.tsv.gz", BuildNcrnaReference.SOURCE_HEADER+"\n1\t/missing/source.gff\tMISSING_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n"
			+"2\tNA\tMISSING_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n");
		write("zero.hits.tsv.gz", BuildNcrnaReference.HIT_HEADER+'\n');
		write("zero.historical.tsv.gz", BuildNcrnaReference.HISTORICAL_HEADER+'\n');
		write("zero.relations.tsv.gz", BuildNcrnaReference.RELATION_HEADER+'\n');
		write("zero.overlaps.tsv.gz", BuildNcrnaReference.OVERLAP_HEADER+'\n');
		write("zero_capture.attribute_summary.tsv", ATTRIBUTE_HEADER+"\n2\t2\t0\t0\t0\t0\t0\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tCURRENT_SOURCE_ATTRIBUTES_NOT_ADJUDICATED\n");
		write("zero_capture.attributes.tsv.gz", "input_id\tgff\tline\traw_attributes\n");
		invoke("zero", "zero_capture", "zero_result");
		check(lines("zero_result.strata_summary.tsv")[1].equals("2\t2\t0\t0\t0\t0\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tSTRATA_COMPLETE_NOT_ADJUDICATED"),
			"Zero-event/source populations remain explicit without any original payload files");
		check(lines("zero_result.joint_strata.tsv").length==1, "Zero events must not invent occupied joint bins");
	}

	private void prepare(){
		write("census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+"\n1\t/missing/one.fa\t110\t100\t1\n2\t/missing/two.fa\t110\t100\t1\n");
		write("members.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+"\nfixture\t1\t1\t/missing/one.fa\tNA\tNA\tNA\tUNKNOWN\t/missing/source.gff\tNA\n"
			+"fixture\t2\t2\t/missing/two.fa\tNA\tNA\tNA\tUNKNOWN\tNA\tNA\n");
		write("scope.tsv", NcrnaReferenceEvidence.SCOPE_HEADER+"\nRF00050\tFMN\tprimary\tfixture\tFMN\tTARGET_CLASS\n"
			+"RF00051\tother_FMN\tprimary\tfixture\tFMN\tTARGET_CLASS\nRF00001\t5S\tprimary\tfixture\t5S\tTARGET_CLASS\n"
			+"RF09999\tunqueried_FMN\trelated\tfixture\tFMN\tRELATED_CLASS\n");
		final StringBuilder assets=new StringBuilder(NcrnaReferenceEvidence.ASSET_HEADER).append('\n');
		for(final String accession:new String[]{"RF00050", "RF00051", "RF00001", "RF09999"}){
			assets.append(accession).append("\tfixture\t10\t20\t00000000000000000001\t/missing/").append(accession).append(".cm\n");
		}
		write("assets.tsv", assets.toString());
		write("queries.tsv", NcrnaReferenceEvidence.QUERY_HEADER+"\nRF00050\ta\t10\nRF00051\tb\t10\nRF00001\tc\t10\n");
		hits=hit(1, 1, "RF00050", 10, 19, '+', 1, 10, 0)
			+hit(1, 2, "RF00001", 15, 24, '+', 1, 10, 1)
			+hit(1, 3, "RF00050", 10, 19, '-', 2, 8, 0)
			+hit(2, 1, "RF00050", 10, 19, '+', 1, 10, 0);
		final StringBuilder history=new StringBuilder(BuildNcrnaReference.HISTORICAL_HEADER).append('\n');
		final StringBuilder relations=new StringBuilder(BuildNcrnaReference.RELATION_HEADER).append('\n');
		final StringBuilder pairs=new StringBuilder();
		for(int line=1; line<=6; line++){
			final String type=line==2 ? "gene" : line==3 ? "exon" : line==5 ? "region" : "ncRNA";
			final boolean linked=line!=4 && line!=5;
			history.append("1\t/missing/source.gff\t").append(line).append("\tΩ%25seq\t").append(type).append("\t10\t19\t")
				.append(line==4 ? '?' : '+').append("\tid").append(line).append(line==5 ? "\t16" : "\t1").append("\tVALID_INTERVAL\tfalse\t")
				.append(line==4 ? "UNASSESSABLE_STRAND" : line==5 ? "STRUCTURAL_CONTEXT" : "ASSESSABLE")
				.append(line==1 ? "\t1" : "\t0").append(linked ? "\t1\t2\t1\t0\t1\n" : "\t0\t0\t0\t0\t-1\n");
			if(linked){
				final String accession=line==1 ? "RF00050" : line==2 ? "RF09999" : line==3 ? "RF99998" : "RF00001";
				relations.append("1\t/missing/source.gff\t").append(line).append("\tRFAM\tid").append(line).append('\t').append(accession).append("\tSOURCE_CLAIM\n");
				pairs.append(link(line, 1, false)).append(link(line, 2, false)).append(link(line, 3, true));
			}
		}
		relations.append("1\t/missing/source.gff\t1\tPARENT\tid1\tid2\tPRESENT_ID\n");
		historyText=history.toString();
		relationsText=relations.toString();
		pairs.append(modelPair(1, 2, 5, false)).append(modelPair(1, 3, 10, true)).append(modelPair(2, 3, 5, true));
		overlaps=pairs.toString();
		copyEvidence("evidence", overlaps);
		attributes=attribute(1, "pseudo=false;Note=do%253Bpseudo%253Dtrue")
			+attribute(2, "pseudo=true;exception=RNA%2520editing")+attribute(3, "partial=true")
			+attribute(4, "pseudo=true")+attribute(5, "partial=true")+attribute(6, "gene_biotype=pseudogene;partial=false;pseudo=false");
		capture("capture", attributes);
	}
	private void copyEvidence(final String name, final String pairs){
		copyEvidence(name, pairs, historyText);
	}
	private void copyEvidence(final String name, final String pairs, final String history){
		copyEvidence(name, pairs, history, hits);
	}
	private void copyEvidence(final String name, final String pairs, final String history, final String nativeHits){
		copyEvidence(name, pairs, history, nativeHits, COVERAGE);
	}
	private void copyEvidence(final String name, final String pairs, final String history, final String nativeHits, final String coverage){
		write(name+".summary.tsv", BuildNcrnaReference.SUMMARY_HEADER+"\n2\t2\t200\t2\t4\t2\t1\t6\t6\t0\t12\t3\t3\t6\t3\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\n");
		write(name+".coverage.tsv", BuildNcrnaReference.COVERAGE_HEADER+'\n'+coverage);
		write(name+".sources.tsv.gz", BuildNcrnaReference.SOURCE_HEADER+"\n1\t/missing/source.gff\tPARSED_GFF\t1\t6\t6\t0\t0\t0\t0\t0\t0\n"
			+"2\tNA\tMISSING_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n");
		write(name+".hits.tsv.gz", BuildNcrnaReference.HIT_HEADER+'\n'+nativeHits);
		write(name+".historical.tsv.gz", history);
		write(name+".relations.tsv.gz", relationsText);
		write(name+".overlaps.tsv.gz", BuildNcrnaReference.OVERLAP_HEADER+'\n'+pairs);
	}
	private void capture(final String name, final String rows){
		write(name+".attribute_summary.tsv", ATTRIBUTE_HEADER+'\n'
			+"2\t2\t6\t6\t1\t6\t6\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tCURRENT_SOURCE_ATTRIBUTES_NOT_ADJUDICATED\n");
		write(name+".attributes.tsv.gz", "input_id\tgff\tline\traw_attributes\n"+rows);
	}
	private static String hit(final int input, final int line, final String accession, final int start, final int stop,
		final char strand, final int from, final int to, final int truncation){
		return input+"\t/missing/hits"+input+".tbl\t"+line+"\tΩ%25seq\t"+start+'\t'+stop+'\t'+strand+'\t'+accession
			+"\tfixture\t"+(accession.equals("RF00001") ? "5S" : "FMN")+"\tTARGET_CLASS\t"+from+'\t'+to+"\t10\t"+truncation
			+"\t50\t1e-8\ttrue\tVALID_INTERVAL\tfalse\tfalse\n";
	}
	private static String link(final int source, final int hit, final boolean opposite){
		return "1\tHISTORICAL_MODEL\t/missing/source.gff\t"+source+"\t/missing/hits1.tbl\t"+hit+'\t'+(hit==2 ? 5 : 10)+'\t'
			+(opposite ? "OPPOSITE_STRAND_GEOMETRY" : "SAME_STRAND_CANDIDATE")+"\tfalse\n";
	}
	private static String modelPair(final int a, final int b, final int overlap, final boolean opposite){
		return "1\tMODEL_MODEL\t/missing/hits1.tbl\t"+a+"\t/missing/hits1.tbl\t"+b+'\t'+overlap+'\t'
			+(opposite ? "OPPOSITE_STRAND_GEOMETRY" : "SAME_STRAND_GEOMETRY")+"\tfalse\n";
	}
	private static String attribute(final int line, final String text){return "1\t/missing/source.gff\t"+line+'\t'+text+'\n';}
	private void invoke(final String evidence, final String capture, final String out){
		BuildNcrnaReference.main(new String[]{"mode=strata", "in="+path(evidence), "attributes="+path(capture), "out="+path(out),
			"census="+path("census.tsv"), "memberships="+path("members.tsv"), "scope="+path("scope.tsv"),
			"assets="+path("assets.tsv"), "queries="+path("queries.tsv"), "t=1", "bgzip=f", "zipthreads=1"});
	}
	private void write(final String name, final String text){try(Output out=new Output(path(name))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}}
	private String[] lines(final String name){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path(name))){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){result.add(new String(row, StandardCharsets.UTF_8));}
		}
		return result.toArray(new String[result.size()]);
	}
	private static String find(final String[] rows, final String start){
		String found=null;
		for(final String row:rows){
			if(row.startsWith(start)){
				check(found==null, "Fixture prefix must identify one row: "+start);
				found=row;
			}
		}
		check(found!=null, "Missing fixture prefix: "+start);
		return found;
	}
	private static void expect(final Runnable action, final String reason){
		try{action.run();}
		catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(reason), "Unexpected fixture failure: "+expected);
			return;
		}
		throw new AssertionError("Expected failure: "+reason);
	}
	private static void check(final boolean condition, final String why){if(!condition){throw new AssertionError(why);}}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private final File root;
	private String hits, historyText, relationsText, overlaps, attributes;
	private static final String ATTRIBUTE_HEADER="inputs_observed\tcatalog_inputs\tretained_source_rows\tmatched_attribute_rows\tsource_paths_read\tunique_selected_lines\tsource_lines_read\tevidence_status\tstatus";
	private static final String COVERAGE="RF00050\tfixture\t10\t10\t0\t2\t3\t2\nRF00050\tfixture\t10\t7\t0\t1\t3\t2\n"
		+"RF00051\tfixture\t10\tNA\tNA\t0\t0\t2\nRF00001\tfixture\t10\t10\t1\t1\t1\t2\n";
}
