package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import prok.BuildNcrnaReference;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;

/** Hand-calculated source concordance over the strata test's literal evidence.
 * No original genome, GFF or model payload exists in this fixture.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceSourceAnchorsTest{
	NcrnaReferenceSourceAnchorsTest(final File root_){root=root_;}
	void run(){
		invoke("evidence", "result", "concordance");
		check(lines("concordance.source_anchor_summary.tsv")[1].equals(
			"2\t2\t4\t4\t12\t3\t2\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tANNOTATION_CONCORDANCE_COMPLETE_NOT_ADMISSION"),
			"Four hits, twelve links, three matching links and two physical anchors are different denominators");
		final String[] loci=lines("concordance.source_anchor_loci.tsv.gz");
		check(loci.length==3, "Opposite-strand and other-input hits must not inherit support");
		final String fmn=find(loci, "1\tΩ%25seq\t+\t10\t19\tRF00050\t");
		check(fmn.contains("\tEXACT_ACCESSION_FLAGGED\t1\t2\t") && fmn.contains("\t2\t5\tNA\tRNA_ROW;GENE_OR_TRANSCRIPT_CONTEXT\t"),
			"FMN has one RNA exact link and one context family-only link, retaining both physical-line witnesses");
		check(fmn.contains("SOURCE_PSEUDO_TRUE") && fmn.contains("SOURCE_RELATED_CLASS")
			&& !fmn.contains("SOURCE_PSEUDOGENE_BIOTYPE") && !fmn.contains("SOURCE_PARTIAL_TRUE"),
			"Warnings must come from a matching source row, not the unrelated overlapping RNA/exon");
		final String[] models=lines("concordance.source_anchor_models.tsv");
		check(find(models, "RF00050\tfixture\t00000000000000000001\t20\t10\tEXACT_ACCESSION_FLAGGED\t")
			.contains("\t1\t1\t1\t2\t0\t0\t1\t1\t0\t1\t0\t1\t0\t30\t30\t30\t30\t30\t"),
			"Multiple source links contribute one full-model score30 aboveGA");
		check(find(models, "RF00001\tfixture\t00000000000000000001\t20\t10\tEXACT_ACCESSION_FLAGGED\t")
			.contains("\t1\t1\t1\t1\t0\t0\t1\t0\t0\t1\t0\t0\t0\t30\t"),
			"A full-span but truncated hit is an observation, not a full-untruncated one");
		check(find(models, "RF00051\tfixture\t00000000000000000001\t20\t10\tEXACT_ACCESSION_CLEAR\t")
			.contains("\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tNA\t"), "Queried models without support remain explicit");
		invoke("sample", "sample_result", "concordance_sample");
		check(lines("concordance_sample.source_anchor_summary.tsv")[1].startsWith("1\t2\t3\t3\t12\t3\t2\tSAMPLE_EVIDENCE_NOT_ADJUDICATED\t"),
			"A saved sample retains its original catalog denominator");
		invoke("zero", "zero_result", "concordance_zero");
		check(lines("concordance_zero.source_anchor_loci.tsv.gz").length==1
			&& lines("concordance_zero.source_anchor_models.tsv").length==19, "Zero hits yield explicit model rows, no invented anchors");
		duplicates(false);
		duplicates(true);
		clearSupport(false);
		clearSupport(true);
		final String original=read("result.source_links.tsv.gz");
		write("bad_native.strata_summary.tsv", read("result.strata_summary.tsv"));
		write("bad_native.source_links.tsv.gz", original.replace("/missing/hits1.tbl", "/missing/wrong.tbl"));
		reject("bad_native", "Source link names no saved native event");
		write("bad_source.strata_summary.tsv", read("result.strata_summary.tsv"));
		write("bad_source.source_links.tsv.gz", original.replace("/missing/source.gff", "/missing/wrong.gff"));
		reject("bad_source", "Source GFF must be declared");
		write("bad_scope.strata_summary.tsv", read("result.strata_summary.tsv").replace("\n2\t2\t", "\n1\t2\t"));
		write("bad_scope.source_links.tsv.gz", original);
		reject("bad_scope", "population must match");
		typedSupport();
		System.err.println("NCRNA_REFERENCE_SOURCE_ANCHORS_TEST_PASS");
	}

	/** Type-only tRNA is family evidence; generic ncRNA and antisense rows are not. */
	private void typedSupport(){
		final NcrnaReferenceSourceAnchorsTest t=new NcrnaReferenceSourceAnchorsTest(new File(root, "typed"));
		check(t.root.mkdir(), "Type-only fixture root must be fresh");
		t.write("census.tsv", prok.NcrnaEvidenceCatalog.CENSUS_HEADER+"\n1\t/missing/typed.fa\t110\t100\t1\n");
		t.write("members.tsv", prok.NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+"\nfixture\t1\t1\t/missing/typed.fa\tNA\tNA\tNA\tUNKNOWN\t/missing/typed.gff\tNA\n");
		t.write("scope.tsv", prok.NcrnaReferenceEvidence.SCOPE_HEADER+"\nRF00005\ttRNA\tprimary\tfixture\ttRNA\tTARGET_CLASS\n");
		t.write("assets.tsv", prok.NcrnaReferenceEvidence.ASSET_HEADER+"\nRF00005\tfixture\t71\t29\t00000000000000000001\t/missing/RF00005.cm\n");
		t.write("queries.tsv", prok.NcrnaReferenceEvidence.QUERY_HEADER+"\nRF00005\ttrna\t71\n");
		t.write("evidence.summary.tsv", BuildNcrnaReference.SUMMARY_HEADER+"\n1\t1\t100\t1\t1\t1\t0\t3\t3\t0\t3\t0\t1\t1\t0\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\n");
		t.write("evidence.sources.tsv.gz", BuildNcrnaReference.SOURCE_HEADER+"\n1\t/missing/typed.gff\tPARSED_GFF\t1\t3\t3\t0\t0\t0\t0\t0\t0\n");
		t.write("evidence.coverage.tsv", BuildNcrnaReference.COVERAGE_HEADER+"\nRF00005\tfixture\t71\t71\t0\t1\t1\t1\n");
		t.write("evidence.hits.tsv.gz", BuildNcrnaReference.HIT_HEADER+"\n1\t/missing/typed.tbl\t1\tchr\t20\t90\t+\tRF00005\tfixture\ttRNA\tTARGET_CLASS\t1\t71\t71\t0\t60\t1e-8\ttrue\tVALID_INTERVAL\tfalse\tfalse\n");
		t.write("evidence.overlaps.tsv.gz", BuildNcrnaReference.OVERLAP_HEADER+"\n"
			+"1\tHISTORICAL_MODEL\t/missing/typed.gff\t1\t/missing/typed.tbl\t1\t71\tSAME_STRAND_CANDIDATE\tfalse\n"
			+"1\tHISTORICAL_MODEL\t/missing/typed.gff\t2\t/missing/typed.tbl\t1\t71\tSAME_STRAND_CANDIDATE\tfalse\n"
			+"1\tHISTORICAL_MODEL\t/missing/typed.gff\t3\t/missing/typed.tbl\t1\t71\tOPPOSITE_STRAND_GEOMETRY\tfalse\n");
		final String summary="inputs_observed\tcatalog_inputs\traw_hit_events\treview_groups\tretained_source_rows\thistorical_model_pairs\tevidence_status\tstatus\n"
			+"1\t1\t1\t1\t3\t3\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tSTRATA_COMPLETE_NOT_ADJUDICATED\n";
		final String links="input_id\tgff\tline\tnative_table\tnative_line\tsource_type\tsource_kind\tassociation\trfam_comparison\tsource_flags\n"
			+"1\t/missing/typed.gff\t1\t/missing/typed.tbl\t1\ttRNA\tRNA_ROW\tSAME_STRAND_OVERLAP\tNO_OVERLAPPING_SOURCE_RFAM\tNONE\n"
			+"1\t/missing/typed.gff\t2\t/missing/typed.tbl\t1\tncRNA\tRNA_ROW\tSAME_STRAND_OVERLAP\tNO_OVERLAPPING_SOURCE_RFAM\tNONE\n"
			+"1\t/missing/typed.gff\t3\t/missing/typed.tbl\t1\ttRNA\tRNA_ROW\tOPPOSITE_STRAND_OVERLAP\tNO_OVERLAPPING_SOURCE_RFAM\tNONE\n";
		t.write("evidence.historical.tsv.gz", BuildNcrnaReference.HISTORICAL_HEADER+"\n"
			+"1\t/missing/typed.gff\t1\tchr\ttRNA\t20\t90\t+\ta\t1\tVALID_INTERVAL\tfalse\tASSESSABLE\t0\t0\t1\t0\t0\t1\n"
			+"1\t/missing/typed.gff\t2\tchr\tncRNA\t20\t90\t+\tb\t1\tVALID_INTERVAL\tfalse\tASSESSABLE\t0\t0\t1\t0\t0\t1\n"
			+"1\t/missing/typed.gff\t3\tchr\ttRNA\t20\t90\t-\tc\t1\tVALID_INTERVAL\tfalse\tASSESSABLE\t0\t0\t0\t1\t0\t-1\n");
		t.write("evidence.relations.tsv.gz", BuildNcrnaReference.RELATION_HEADER+'\n');
		t.write("capture.attribute_summary.tsv", "inputs_observed\tcatalog_inputs\tretained_source_rows\tmatched_attribute_rows\tsource_paths_read\tunique_selected_lines\tsource_lines_read\tevidence_status\tstatus\n"
			+"1\t1\t3\t3\t1\t3\t3\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tCURRENT_SOURCE_ATTRIBUTES_NOT_ADJUDICATED\n");
		t.write("capture.attributes.tsv.gz", "input_id\tgff\tline\traw_attributes\n1\t/missing/typed.gff\t1\tID=a\n1\t/missing/typed.gff\t2\tID=b\n1\t/missing/typed.gff\t3\tID=c\n");
		BuildNcrnaReference.main(new String[]{"mode=strata", "in="+t.path("evidence"), "attributes="+t.path("capture"), "out="+t.path("support"),
			"census="+t.path("census.tsv"), "memberships="+t.path("members.tsv"), "scope="+t.path("scope.tsv"),
			"assets="+t.path("assets.tsv"), "queries="+t.path("queries.tsv"), "t=1", "bgzip=f", "zipthreads=1"});
		check(t.read("support.source_links.tsv.gz").equals(links), "The actual strata producer must emit the independently specified type-only support schema");
		check(t.read("support.strata_summary.tsv").equals(summary), "The produced typed fixture must conserve its declared source population");
		t.invoke("evidence", "support", "result");
		check(t.lines("result.source_anchor_summary.tsv")[1].startsWith("1\t1\t1\t1\t3\t1\t1\t"), "Only the explicit same-strand tRNA row supplies type-only support");
		final String row=t.lines("result.source_anchor_loci.tsv.gz")[1];
		check(row.contains("\tTYPE_FAMILY_ONLY_CLEAR\t1\t1\t")
			&& row.contains("\tNA\tNA\t2\tRNA_ROW\tNONE\t1\t60\t31\t1\t0\tfalse\tfalse\t"),
			"The type witness is explicit, distinct from RFAM and warning bits, with margin31");
		t.write("unknown.strata_summary.tsv", summary);
		t.write("unknown.source_links.tsv.gz", links.replace("NO_OVERLAPPING_SOURCE_RFAM", "SOURCE_UNKNOWN_RFAM"));
		t.invoke("evidence", "unknown", "unknown_result");
		check(t.lines("unknown_result.source_anchor_loci.tsv.gz").length==1, "Unknown RFAM claims cannot be bypassed through the feature type");
		t.write("mixed.strata_summary.tsv", summary);
		t.write("mixed.source_links.tsv.gz", links.replace(
			"\t2\t/missing/typed.tbl\t1\tncRNA\tRNA_ROW\tSAME_STRAND_OVERLAP\tNO_OVERLAPPING_SOURCE_RFAM\t",
			"\t2\t/missing/typed.tbl\t1\tncRNA\tRNA_ROW\tSAME_STRAND_OVERLAP\tSOURCE_SAME_ACCESSION\t"));
		t.invoke("evidence", "mixed", "mixed_result");
		check(t.lines("mixed_result.source_anchor_summary.tsv")[1].startsWith("1\t1\t1\t1\t3\t2\t1\t"),
			"Typed and exact RFAM source links support one physical anchor, not two");
		final String mixed=t.lines("mixed_result.source_anchor_loci.tsv.gz")[1];
		check(mixed.contains("\tEXACT_ACCESSION_CLEAR\t1\t2\t")
			&& mixed.contains("\t3\tNA\t2\tRNA_ROW\tNONE\t1\t60\t31\t"),
			"RFAM basis wins while both source witness lines and the single score remain explicit");
		// Reuse this synthetic71-base/GA29 model shape to test family attribution,
		// not biological tmRNA scores. Its metadata release is explicitly fixture.
		final NcrnaReferenceSourceAnchorsTest m=new NcrnaReferenceSourceAnchorsTest(new File(root, "typed_tmrna"));
		check(m.root.mkdir(), "tmRNA attribution fixture root must be fresh");
		for(final String file:new String[]{"census.tsv", "members.tsv", "scope.tsv", "assets.tsv", "queries.tsv",
			"evidence.summary.tsv", "evidence.sources.tsv.gz", "evidence.coverage.tsv", "evidence.hits.tsv.gz",
			"evidence.overlaps.tsv.gz", "support.strata_summary.tsv", "support.source_links.tsv.gz"}){
			m.write(file, t.read(file).replace("RF00005", "RF00023").replace("\ttRNA\t", "\ttmRNA\t"));
		}
		m.invoke("evidence", "support", "result");
		check(m.lines("result.source_anchor_loci.tsv.gz")[1].contains("\tRF00023\tTYPE_FAMILY_ONLY_CLEAR\t1\t1\t"),
			"Literal tmRNA type provides family-only support to the declared tmRNA model family");
		m.write("mismatch.strata_summary.tsv", summary);
		m.write("mismatch.source_links.tsv.gz", m.read("support.source_links.tsv.gz")
			.replace("\ttmRNA\tRNA_ROW\tSAME_STRAND_OVERLAP", "\ttRNA\tRNA_ROW\tSAME_STRAND_OVERLAP"));
		m.invoke("evidence", "mismatch", "mismatch_result");
		check(m.lines("mismatch_result.source_anchor_loci.tsv.gz").length==1,
			"A tRNA feature cannot supply type-only support to a tmRNA model");
	}

	/** One direct RNA link remains after the flagged context ceases to match. */
	private void clearSupport(final boolean familyOnly){
		final String name=familyOnly ? "family_only" : "clear_exact";
		final StringBuilder selected=new StringBuilder();
		for(String row:lines("result.source_links.tsv.gz")){
			if(row.startsWith("1\t/missing/source.gff\t2\t/missing/hits1.tbl\t1\t")){
				row=row.replace("SOURCE_SAME_FAMILY_OTHER_ACCESSION;SOURCE_RELATED_CLASS", "SOURCE_OTHER_KNOWN_FAMILY");
			}
			if(familyOnly && row.startsWith("1\t/missing/source.gff\t1\t/missing/hits1.tbl\t1\t")){
				row=row.replace("SOURCE_SAME_ACCESSION", "SOURCE_SAME_FAMILY_OTHER_ACCESSION");
			}
			selected.append(row).append('\n');
		}
		write(name+".strata_summary.tsv", read("result.strata_summary.tsv"));
		write(name+".source_links.tsv.gz", selected.toString());
		invoke("evidence", name, name+"_result");
		final String row=find(lines(name+"_result.source_anchor_loci.tsv.gz"), "1\tΩ%25seq\t+\t10\t19\tRF00050\t");
		check(row.contains(familyOnly ? "\tFAMILY_ONLY_CLEAR\t1\t1\t" : "\tEXACT_ACCESSION_CLEAR\t1\t1\t"),
			"A matching RNA annotation keeps exact-accession versus family-only attribution distinct");
		check(row.contains("\tRNA_ROW\tNONE\t1\t50\t30\t1\t0\tfalse\tfalse\t"),
			"Unrelated source warnings cannot contaminate this directly supported score observation");
	}

	/** Same physical model interval appears twice; differing scores may not pick a winner. */
	private void duplicates(final boolean conflicting){
		final String name=conflicting ? "conflicting_native" : "repeated_native";
		write(name+".summary.tsv", read("evidence.summary.tsv").replace("\t4\t2\t1\t6\t6\t0\t12\t3\t", "\t5\t2\t1\t6\t6\t0\t16\t6\t"));
		write(name+".sources.tsv.gz", read("evidence.sources.tsv.gz"));
		write(name+".coverage.tsv", read("evidence.coverage.tsv")
			.replace("RF00050\tfixture\t10\t10\t0\t2\t3\t2", "RF00050\tfixture\t10\t10\t0\t3\t4\t2")
			.replace("RF00050\tfixture\t10\t7\t0\t1\t3\t2", "RF00050\tfixture\t10\t7\t0\t1\t4\t2"));
		final String duplicate="1\t/missing/hits1.tbl\t4\tΩ%25seq\t10\t19\t+\tRF00050\tfixture\tFMN\tTARGET_CLASS\t1\t10\t10\t0\t"
			+(conflicting ? 51 : 50)+"\t1e-8\ttrue\tVALID_INTERVAL\tfalse\tfalse\n";
		write(name+".hits.tsv.gz", read("evidence.hits.tsv.gz")+duplicate);
		final StringBuilder pairs=new StringBuilder(read("evidence.overlaps.tsv.gz"));
		final StringBuilder sourceLinks=new StringBuilder(read("result.source_links.tsv.gz"));
		for(final String row:lines("result.source_links.tsv.gz")){
			if(!row.contains("\t/missing/hits1.tbl\t1\t")){continue;}
			sourceLinks.append(row.replace("\t/missing/hits1.tbl\t1\t", "\t/missing/hits1.tbl\t4\t")).append('\n');
		}
		for(final int source:new int[]{1, 2, 3, 6}){
			pairs.append("1\tHISTORICAL_MODEL\t/missing/source.gff\t").append(source)
				.append("\t/missing/hits1.tbl\t4\t10\tSAME_STRAND_CANDIDATE\tfalse\n");
		}
		pairs.append("1\tMODEL_MODEL\t/missing/hits1.tbl\t1\t/missing/hits1.tbl\t4\t10\tSAME_STRAND_GEOMETRY\tfalse\n")
			.append("1\tMODEL_MODEL\t/missing/hits1.tbl\t2\t/missing/hits1.tbl\t4\t5\tSAME_STRAND_GEOMETRY\tfalse\n")
			.append("1\tMODEL_MODEL\t/missing/hits1.tbl\t3\t/missing/hits1.tbl\t4\t10\tOPPOSITE_STRAND_GEOMETRY\tfalse\n");
		write(name+".overlaps.tsv.gz", pairs.toString());
		write(name+".strata_summary.tsv", read("result.strata_summary.tsv").replace("\n2\t2\t4\t3\t6\t12\t", "\n2\t2\t5\t3\t6\t16\t"));
		write(name+".source_links.tsv.gz", sourceLinks.toString());
		invoke(name, name, name+"_result");
		check(lines(name+"_result.source_anchor_summary.tsv")[1].startsWith("2\t2\t5\t4\t16\t5\t2\t"),
			"Five native events collapse to four physical intervals, only two with annotation support");
		final String fmn=find(lines(name+"_result.source_anchor_models.tsv"), "RF00050\tfixture\t00000000000000000001\t20\t10\tEXACT_ACCESSION_FLAGGED\t");
		check(fmn.contains(conflicting ? "\t1\t1\t2\t4\t1\t1\t0\t0\t0\t" : "\t1\t1\t2\t4\t1\t0\t1\t1\t0\t"),
			"An identical physical interval counts once; disagreeing native observations retain no selected score");
	}
	private void invoke(final String evidence, final String strata, final String out){
		BuildNcrnaReference.main(new String[]{"mode=sourceanchors", "in="+path(evidence), "sources="+path(strata), "out="+path(out),
			"census="+path("census.tsv"), "memberships="+path("members.tsv"), "scope="+path("scope.tsv"),
			"assets="+path("assets.tsv"), "queries="+path("queries.tsv"), "t=1", "bgzip=f", "zipthreads=1"});
	}
	private void reject(final String source, final String reason){
		try{invoke("evidence", source, source+"_result");}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(reason), "Unexpected source-concordance failure: "+failure);
			check(!new File(path(source+"_result.source_anchor_summary.tsv")).exists(), "Rejected metadata cannot publish completion");
			return;
		}
		throw new AssertionError("Expected source-concordance rejection: "+reason);
	}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private void write(final String name, final String value){try(Output out=new Output(path(name))){out.writer.print(value.getBytes(StandardCharsets.UTF_8));}}
	private String[] lines(final String name){
		final ArrayList<String> rows=new ArrayList<String>();
		try(Input in=new Input(path(name))){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){rows.add(new String(row, StandardCharsets.UTF_8));}
		}
		return rows.toArray(new String[rows.size()]);
	}
	private String read(final String name){return String.join("\n", lines(name))+"\n";}
	private static String find(final String[] rows, final String prefix){
		String result=null;
		for(final String row:rows){
			if(!row.startsWith(prefix)){continue;}
			check(result==null, "Expected one row for "+prefix);
			result=row;
		}
		check(result!=null, "Missing row for "+prefix);
		return result;
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final File root;
}
