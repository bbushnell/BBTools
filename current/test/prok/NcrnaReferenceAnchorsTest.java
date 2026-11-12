package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import prok.BuildNcrnaReference;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaReferenceAnchors;
import prok.NcrnaReferenceEvidence;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;

/** Hand-calculated current-hit joins, with prior and construction evidence separate.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceAnchorsTest{
	NcrnaReferenceAnchorsTest(final File root_){root=root_;}
	void run(){
		check(root.mkdir(), "Anchor fixture root must be new");
		prepare(false);
		invoke("manifest.tsv", "good");
		final String[] summary=lines("good.anchor_summary.tsv");
		check(summary.length==3 && summary[1].equals("bound\t8\t1\t0\t1\t6\t2\t2\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tPRIOR_EVIDENCE_JOIN_NOT_ADMISSION")
			&& summary[2].equals("tpp\t3\t1\t1\t0\t1\t2\t2\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tPRIOR_EVIDENCE_JOIN_NOT_ADMISSION"),
			"Every original row has exactly one accounting outcome and retains the complete input-population denominator");
		final String[] rows=lines("good.anchor_rows.tsv.gz");
		check(rows.length==12, "Eight bound decisions and three TPP rows must survive");
		check(find(rows, "bound\t"+path("bound.tsv")+"\t2\t").endsWith("\t1\t1\tSINGLE_CURRENT_OWN_MODEL_HIT"),
			"Own and other models at the same locus remain separate observations");
		check(find(rows, "bound\t"+path("bound.tsv")+"\t3\t").endsWith("\tbound\t2\tNA\tNA\tREPEATED_EXACT_ANCHOR"),
			"A repeated physical anchor is retained with a pointer, not another weighted observation");
		check(find(rows, "bound\t"+path("bound.tsv")+"\t4\t").endsWith("\t0\t0\tNO_CURRENT_OWN_MODEL_OVERLAP"),
			"A same-name contig in another input and an opposite-strand hit cannot supply support");
		check(find(rows, "bound\t"+path("bound.tsv")+"\t5\t").endsWith("\t2\t0\tMULTIPLE_CURRENT_OWN_MODEL_HITS"), "Do not silently choose one hit from a multi-hit locus");
		check(find(rows, "bound\t"+path("bound.tsv")+"\t6\t").endsWith("\t1\t0\tSINGLE_OWN_HIT_WEAK_OVERLAP"), "A one-base touch cannot calibrate a scored match");
		final String fmn=find(rows, "bound\t"+path("bound.tsv")+"\t7\t");
		check(fmn.contains("\tPRIOR_ACCEPTED_WHOLE_INTERVAL\t60\t69\t") && fmn.endsWith("NO_CURRENT_OWN_MODEL_OVERLAP"),
			"Prior accepted FMN survives without a current hit or newly invented aptamer bounds");
		check(find(rows, "bound\t"+path("bound.tsv")+"\t8\t").endsWith("\tNA\tNA\tCURRENT_MODEL_NOT_SEARCHED"), "Unsearched differs from searched-without-support");
		check(find(rows, "tpp\t"+path("tpp.tsv")+"\t2\t").contains("\tΩ%25seq\tRF00059\tCM_SUPPORTED_COMPLETE_CANDIDATE\tREVIEWED_CONSTRUCTION_CORE\t"),
			"Plain TPP target framing must join the escaped UTF-8 current target without changing evidence basis");
		check(find(rows, "tpp\t"+path("tpp.tsv")+"\t3\t").endsWith("\tNA\tNA\tOUTSIDE_DECLARED_INPUT_PATHS"), "No path guessing for prior genomes outside this catalog");
		final String[] models=lines("good.anchor_models.tsv");
		final String cob=find(models, "RF00174\tfixture\t00000000000000000001\t20\t10\ttrue\tPRIOR_ACCEPTED_CORE\t");
		check(cob.contains("\t5\t0\t1\t4\t2\t1\t1\t0\t1\t0\t1\t40\t40\t40\t40\t40\t1\t1\t1\t1\t1\t1\t"),
			"Only the unique qualifying own-model hit contributes margin40/coverage1; weak, duplicate and multi-hit cases stay outside that distribution");
		check(cob.endsWith("\t40\t40\t40\t40\t40\tOBSERVED_HITS_ONLY_NOT_ADMISSION_THRESHOLDS"),
			"The full untruncated subset retains its own five margin quantiles");
		final String tpp=find(models, "RF00059\tfixture\t00000000000000000001\t30\t10\ttrue\tREVIEWED_CONSTRUCTION_CORE\t");
		check(tpp.contains("\t2\t1\t0\t1\t1\t0\t0\t0\t0\t1\t0\t15\t15\t15\t15\t15\t0.700000\t"),
			"Inside TPP has margin15 and coverage7/10; its outside-path row is counted but contributes no score");
		check(tpp.endsWith("\tNA\tNA\tNA\tNA\tNA\tOBSERVED_HITS_ONLY_NOT_ADMISSION_THRESHOLDS"),
			"A partial current hit must not enter the full-model margin distribution");
		check(lines("good.anchor_links.tsv.gz").length==7, "Six actual same-context hit links survive, including weak overlap and both competing models");
		badBound("wrong_length", bound.replace("\t100\tRF00174", "\t101\tRF00174"), "contig length");
		badBound("missing_target", bound.replace("Ω%25seq", "missing"), "literal target");
		badBound("wrong_path", bound.replace("\t/missing/one.fa\tΩ", "\t/missing/two.fa\tΩ"), "exact selected census");
		badBound("unknown_assignment", bound.replace("\tUNREVIEWED\t", "\tINVENTED\t"), "Unknown prior assignment");
		write("bad_tpp.tsv", tppRows.replace("\tg1\ta\t", "\tg1\tb\t"));
		write("bad_tpp.manifest.tsv", NcrnaReferenceAnchors.MANIFEST_HEADER+'\n'+manifest("tpp", "tpp_construction_v1", "bad_tpp.tsv", path("roster.tsv"), 3));
		reject("bad_tpp.manifest.tsv", "bad_tpp", "assembly differs");
		write("wrong_count.manifest.tsv", NcrnaReferenceAnchors.MANIFEST_HEADER+'\n'+manifest("bound", "bound_prior_v1", "bound.tsv", "NA", 9));
		reject("wrong_count.manifest.tsv", "wrong_count", "accounted for exactly once");
		sampledQuantiles();
		System.err.println("NCRNA_REFERENCE_ANCHORS_TEST_PASS");
	}

	/** Three observations in unsorted score order have nearest ranks1/1/2/3/3.
	 * A third unobserved catalog input must leave the output explicitly sampled.
	 */
	private void sampledQuantiles(){
		final NcrnaReferenceAnchorsTest sample=new NcrnaReferenceAnchorsTest(new File(root, "sample"));
		check(sample.root.mkdir(), "Sampled fixture root must be new");
		sample.prepare(true);
		sample.invoke("manifest.tsv", "result");
		check(sample.lines("result.anchor_summary.tsv")[1].equals("bound\t3\t0\t0\t0\t3\t2\t3\tSAMPLE_EVIDENCE_NOT_ADJUDICATED\tPRIOR_EVIDENCE_JOIN_NOT_ADMISSION"),
			"Two observed inputs out of three must remain explicitly sampled, even when every declared anchor has support");
		final String row=find(sample.lines("result.anchor_models.tsv"), "RF00174\tfixture\t00000000000000000001\t20\t10\ttrue\tPRIOR_ACCEPTED_CORE\t");
		check(row.endsWith("\t40\t40\t50\t60\t60\t1\t1\t1\t1\t1\t0\t40\t40\t50\t60\t60\tOBSERVED_HITS_ONLY_NOT_ADMISSION_THRESHOLDS"),
			"Margins60/40/50 sort to40/50/60, with nearest-rank min/p05/median/p95/max40/40/50/60/60 in both full and all-hit distributions");
	}

	private void prepare(final boolean sample){
		write("census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+"\n1\t/missing/one.fa\t101\t100\t1\n2\t/missing/two.fa\t101\t100\t1\n"
			+(sample ? "3\t/missing/three.fa\t101\t100\t1\n" : ""));
		write("members.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+"\nfixture\t1\t1\t/old/one.fa\tNA\tNA\tNA\tprokaryotic\tNA\tNA\n"
			+"fixture\t2\t2\t/missing/two.fa\tNA\tNA\tNA\tprokaryotic\tNA\tNA\n"
			+(sample ? "fixture\t3\t3\t/missing/three.fa\tNA\tNA\tNA\tprokaryotic\tNA\tNA\n" : ""));
		String scope=NcrnaReferenceEvidence.SCOPE_HEADER+'\n', assets=NcrnaReferenceEvidence.ASSET_HEADER+'\n';
		for(final String accession:new String[]{"RF00050", "RF00174", "RF00059", "RF00011"}){
			scope+=accession+'\t'+accession+"\tprimary\tcore\t"+accession+"\tTARGET_CLASS\n";
			final String ga=accession.equals("RF00050") ? "66.2" : accession.equals("RF00059") ? "30" : "20";
			assets+=accession+"\tfixture\t10\t"+ga+"\t00000000000000000001\t/missing/"+accession+".cm\n";
		}
		write("scope.tsv", scope);
		write("assets.tsv", assets);
		write("queries.tsv", NcrnaReferenceEvidence.QUERY_HEADER+"\nRF00050\tfmn\t10\nRF00174\tcob\t10\nRF00059\ttpp\t10\n");
		write("evidence.summary.tsv", BuildNcrnaReference.SUMMARY_HEADER+"\n2\t"+(sample ? 3 : 2)
			+"\t200\t2\t7\t2\t2\t0\t0\t0\t0\t1\t3\t6\t2\t"
			+(sample ? "SAMPLE_EVIDENCE_NOT_ADJUDICATED\n" : "EVIDENCE_COMPLETE_NOT_ADJUDICATED\n"));
		write("evidence.sources.tsv.gz", BuildNcrnaReference.SOURCE_HEADER+"\n1\tNA\tMISSING_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n2\tNA\tMISSING_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n");
		write("evidence.coverage.tsv", BuildNcrnaReference.COVERAGE_HEADER+"\nRF00050\tfixture\t10\t10\t0\t1\t1\t2\nRF00174\tfixture\t10\t10\t0\t5\t5\t2\nRF00059\tfixture\t10\t7\t2\t1\t1\t2\n");
		write("evidence.contigs.tsv.gz", "input_id\tfasta\tliteral_target\tobserved_length\n1\t/missing/one.fa\tΩ%25seq\t100\n2\t/missing/two.fa\tΩ%25seq\t100\n");
		write("evidence.hits.tsv.gz", BuildNcrnaReference.HIT_HEADER+'\n'
			+hit(1, 1, "RF00174", 10, 19, '+', 1, 10, 0, 60)+hit(2, 1, "RF00174", 10, 19, '-', 1, 10, 0, 60)
			+hit(1, 2, "RF00174", 40, 49, '+', 1, 10, 0, 70)+hit(1, 3, "RF00174", 50, 59, '+', 1, 10, 0, 80)
			+hit(1, 4, "RF00059", 70, 79, '+', 2, 8, 2, 45)+hit(1, 5, "RF00050", 10, 19, '+', 1, 10, 0, 90)
			+hit(1, 6, "RF00174", 90, 99, '+', 1, 10, 0, 50));
		write("evidence.overlaps.tsv.gz", BuildNcrnaReference.OVERLAP_HEADER+"\n1\tMODEL_MODEL\t/native1.tbl\t1\t/native1.tbl\t5\t10\tSAME_STRAND_GEOMETRY\tfalse\n");
		bound=NcrnaReferenceAnchors.BOUND_HEADER+'\n'+prior(2, 1, "RF00174", "RF00174", 10, 19)
			+prior(3, 1, "RF00174", "RF00174", 10, 19)+prior(4, 2, "RF00174", "RF00174", 10, 19)
			+prior(5, 1, "RF00174", "RF00174", 40, 59)+prior(6, 1, "RF00174", "RF00174", 80, 90)
			+prior(7, 2, "RF00050", "ACCEPTED_NOVEL_FMN", 60, 69)+prior(8, 1, "RF00011", "RF00011", 20, 29)
			+prior(9, 1, "RF00174", "UNREVIEWED", 30, 39);
		if(sample){
			bound=NcrnaReferenceAnchors.BOUND_HEADER+'\n'+prior(2, 1, "RF00174", "RF00174", 50, 59)
				+prior(3, 1, "RF00174", "RF00174", 10, 19)+prior(4, 1, "RF00174", "RF00174", 40, 49);
		}
		write("bound.tsv", bound);
		tppRows=NcrnaReferenceAnchors.TPP_HEADER+'\n'+tpp("g1", true)+tpp("outside", true)+tpp("g1", false);
		write("tpp.tsv", tppRows);
		write("roster.tsv", "index\tgenome\tassembly\ttaxid\tdomain\tphylum\tgff\tfasta\tloci\tgenome_bases\tcontigs\n"
			+"0\tg1\ta\tNA\tNA\tNA\tNA\t/old/one.fa\t1\t100\t1\n1\toutside\ta\tNA\tNA\tNA\tNA\t/outside.fa\t1\t100\t1\n");
		write("manifest.tsv", NcrnaReferenceAnchors.MANIFEST_HEADER+'\n'
			+manifest("bound", "bound_prior_v1", "bound.tsv", "NA", sample ? 3 : 8)+manifest("tpp", "tpp_construction_v1", "tpp.tsv", path("roster.tsv"), 3));
	}
	private static String hit(final int input, final int line, final String model, final int start, final int stop,
		final char strand, final int from, final int to, final int truncation, final int bits){
		return input+"\t/native"+input+".tbl\t"+line+"\tΩ%25seq\t"+start+'\t'+stop+'\t'+strand+'\t'+model+"\tfixture\t"+model
			+"\tTARGET_CLASS\t"+from+'\t'+to+"\t10\t"+truncation+'\t'+bits+"\t1e-8\ttrue\tVALID_INTERVAL\tfalse\tfalse\n";
	}
	private static String prior(final int line, final int input, final String family, final String assignment, final int start, final int stop){
		final boolean core=family.equals(assignment);
		return "old\t/original.tsv\t"+line+"\tid"+line+"\tg\t/missing/"+(input==1 ? "one" : "two")+".fa\t"+input+"\t/missing/"
			+(input==1 ? "one" : "two")+".fa\tΩ%25seq\t100\t"+family+'\t'+assignment+'\t'+start+'\t'+stop+'\t'
			+(core ? start : 0)+'\t'+(core ? stop : 0)+"\t+\t"+(assignment.equals("ACCEPTED_NOVEL_FMN") ? "PRIOR_REVIEWED_WHOLE_INTERVAL_NOT_NEW_APTAMER_BOUNDS" : "SOURCE_AND_REVIEWED_CORE")
			+"\tfixture_authority\tNA\told_score_42.5\n";
	}
	private static String tpp(final String genome, final boolean eligible){
		final String[] fields=new String[46];
		Arrays.fill(fields, "NA");
		fields[0]=genome+eligible;
		fields[1]=genome;
		fields[2]="a";
		fields[6]="Ω%seq";
		fields[7]="+";
		fields[15]="100";
		fields[39]=eligible ? "70" : "NA";
		fields[40]=eligible ? "79" : "NA";
		fields[41]=Boolean.toString(eligible);
		fields[42]=eligible ? "CM_SUPPORTED_COMPLETE_CANDIDATE" : "UNRESOLVED";
		return String.join("\t", fields)+'\n';
	}
	private String manifest(final String name, final String format, final String file, final String bindings, final int count){
		return name+'\t'+format+'\t'+path(file)+'\t'+bindings+"\tfixture_review\t"+count+'\n';
	}
	private void invoke(final String manifest, final String output){
		BuildNcrnaReference.main(new String[]{"mode=anchors", "in="+path("evidence"), "anchors="+path(manifest),
			"census="+path("census.tsv"), "memberships="+path("members.tsv"), "scope="+path("scope.tsv"),
			"assets="+path("assets.tsv"), "queries="+path("queries.tsv"), "out="+path(output), "t=1", "bgzip=f", "zipthreads=1"});
	}
	private void badBound(final String name, final String content, final String reason){
		write(name+".tsv", content);
		write(name+".manifest.tsv", NcrnaReferenceAnchors.MANIFEST_HEADER+'\n'+manifest("bound", "bound_prior_v1", name+".tsv", "NA", 8));
		reject(name+".manifest.tsv", name, reason);
	}
	private void reject(final String manifest, final String output, final String reason){
		try{invoke(manifest, output);}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(reason), "Unexpected anchor failure: "+failure);
			check(!new File(path(output+".anchor_summary.tsv")).exists(), "Rejected anchor input cannot publish completion");
			return;
		}
		throw new AssertionError("Expected failure containing "+reason);
	}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private void write(final String name, final String text){
		try(Output out=new Output(path(name))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
	}
	private String[] lines(final String name){
		final ArrayList<String> rows=new ArrayList<String>();
		try(Input in=new Input(path(name))){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){rows.add(new String(row, StandardCharsets.UTF_8));}
		}
		return rows.toArray(new String[rows.size()]);
	}
	private static String find(final String[] rows, final String prefix){
		String found=null;
		for(final String row:rows){
			if(!row.startsWith(prefix)){continue;}
			check(found==null, "Expected unique row: "+prefix);
			found=row;
		}
		check(found!=null, "Missing row: "+prefix);
		return found;
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final File root;
	private String bound, tppRows;
}
