package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import prok.BuildNcrnaReference;
import prok.CobalaminPanelPrep;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaReferenceEvidence;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;

/** Exact prior-decision conservation with no sequence files or current hit table.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferencePriorTest{

	NcrnaReferencePriorTest(final File root_){root=root_;}

	void run(){
		check(root.mkdir(), "Prior fixture root must be new");
		prepare();
		invoke("manifest.tsv", "members.tsv", "result");
		check(lines("result.prior_summary.tsv")[1].equals("2\t7\t2\t2\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tPRIOR_DECISIONS_BOUND_NOT_NEWLY_ADJUDICATED"),
			"All seven original rows must survive even when the completed scan has zero hits");
		final String[] rows=lines("result.prior_decisions.tsv.gz");
		check(rows.length==8, "Five panel and two FMN decisions must survive exactly once");
		check(find(rows, "panel\t"+path("panel.tsv")+"\t2\t").contains("\tsame\tg1\t/recorded/one.fa\t1\t/missing/one.fa\tΩ%25seq\t100\tRF00174\tRF00174\t1\t20\t2\t18\t+\t"),
			"The original recorded FASTA and resolved input remain separate, with unchanged source/core bounds");
		check(find(rows, "panel\t"+path("panel.tsv")+"\t3\t").contains("\tsame\tg2\t/missing/two.fa\t2\t/missing/two.fa\tΩ%25seq\t"),
			"Repeated original locus and literal contig names in different inputs must not collide");
		final String accepted=find(rows, "fmn\t"+path("fmn.tsv")+"\t2\t");
		check(accepted.contains("\tRF00050\tACCEPTED_NOVEL_FMN\t50\t80\t0\t0\t-\tPRIOR_REVIEWED_WHOLE_INTERVAL_NOT_NEW_APTAMER_BOUNDS\t"),
			"Accepted FMN survives without a new hit or invented aptamer bounds");
		check(accepted.endsWith("Ω%25seq%0950%0980%09-%09ACCEPTED_NOVEL_FMN%09model%250%090.97%0931%097%0942.5%0966.2%09BELOW_GA%09owner%09reviewer%0900000000000000000001%0900000000000000000002"),
			"The complete original FMN row, including its authority and evidence fields, must remain lossless");
		final String[] counts=lines("result.prior_counts.tsv");
		for(final String expected:new String[]{"panel\tRF00174\t2", "panel\tB12clade\t1", "panel\tUNRESOLVED\t1", "panel\tUNREVIEWED\t1",
			"fmn\tACCEPTED_NOVEL_FMN\t1", "fmn\tUNRESOLVED\t1"}){
			check(contains(counts, expected), "Unchanged assignment count: "+expected);
		}
		badBinding("missing", bindingText.replace("/missing/two.fa", "/not-recorded/two.fa"), "Original FASTA is absent");
		badBinding("changed_coordinate", bindingText.replace("\t50\t80\t", "\t51\t80\t"), "binding fields differ");
		badBinding("changed_decision", bindingText.replace("ACCEPTED_NOVEL_FMN", "UNRESOLVED"), "binding fields differ");
		badBinding("duplicate_binding", bindingText+bindingText.substring(bindingText.indexOf('\n')+1, bindingText.indexOf("\n3\t")+1), "exactly one original FASTA");
		badBinding("absent_binding", FMN_BINDING_HEADER+'\n', "exact source schema and FASTA binding");
		write("ambiguous.members.tsv", membersText.replace("fixture\t2\t2\t/missing/two.fa", "fixture\t2\t2\t/recorded/one.fa"));
		expect(()->invoke("manifest.tsv", "ambiguous.members.tsv", "bad_ambiguous"), "maps to multiple catalog inputs");
		check(!new File(path("bad_ambiguous.prior_summary.tsv")).exists(), "An ambiguous path cannot publish completion");
		badPanel("missing_target", panelText.replace("Ω%seq", "other"), 5, "no observed literal contig");
		badPanel("outside", panelText.replace("\t1\t20\t2\t18\t", "\t1\t101\t2\t18\t"), 5, "interval leaves");
		badPanel("duplicate_key", panelText+"same\tg1\tΩ%seq\tRF00174\tRF00174\t1\t20\t2\t18\t+\n", 6, "Repeated original family/locus key");
		badPanel("wrong_count", panelText, 6, "Prior row count changed");
		badPanel("invented_core", panelText.replace("\tUNRESOLVED\t40\t50\t0\t0\t", "\tUNRESOLVED\t40\t50\t41\t49\t"), 5, "core-presence fields disagree");
		write("empty.panel.tsv", CobalaminPanelPrep.HEADER+'\n');
		write("empty.manifest.tsv", MANIFEST_HEADER+'\n'+manifestRow("panel", "cobalamin_panel_v1", "empty.panel.tsv", "roster.tsv", 0));
		invoke("empty.manifest.tsv", "members.tsv", "empty");
		check(lines("empty.prior_decisions.tsv.gz").length==1 && lines("empty.prior_counts.tsv").length==1
			&& lines("empty.prior_summary.tsv")[1].startsWith("1\t0\t2\t2\t"), "An explicitly empty prior source must not invent decisions");
		System.err.println("NCRNA_REFERENCE_PRIOR_TEST_PASS");
	}

	private void prepare(){
		write("census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+"\n1\t/missing/one.fa\t110\t100\t1\n2\t/missing/two.fa\t110\t100\t1\n");
		membersText=NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+"\nfixture\t1\t1\t/recorded/one.fa\tg1\tNA\tNA\tUNKNOWN\tNA\tNA\n"
			+"fixture\t2\t2\t/missing/two.fa\tg2\tNA\tNA\tUNKNOWN\tNA\tNA\n";
		write("members.tsv", membersText);
		write("scope.tsv", NcrnaReferenceEvidence.SCOPE_HEADER+"\nRF00050\tFMN\tprimary\taptamer_only\tFMN\tTARGET_CLASS\n");
		write("assets.tsv", NcrnaReferenceEvidence.ASSET_HEADER+"\nRF00050\tfixture\t10\t66.2\t00000000000000000001\t/missing/fmn.cm\n");
		write("queries.tsv", NcrnaReferenceEvidence.QUERY_HEADER+"\nRF00050\tfmn\t10\n");
		write("evidence.summary.tsv", BuildNcrnaReference.SUMMARY_HEADER+"\n2\t2\t200\t2\t0\t2\t2\t0\t0\t0\t0\t0\t1\t2\t2\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\n");
		write("evidence.coverage.tsv", BuildNcrnaReference.COVERAGE_HEADER+"\nRF00050\tfixture\t10\tNA\tNA\t0\t0\t2\n");
		write("evidence.sources.tsv.gz", BuildNcrnaReference.SOURCE_HEADER+"\n1\tNA\tMISSING_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n"
			+"2\tNA\tMISSING_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n");
		write("evidence.contigs.tsv.gz", "input_id\tfasta\tliteral_target\tobserved_length\n1\t/missing/one.fa\tΩ%25seq\t100\n2\t/missing/two.fa\tΩ%25seq\t100\n");
		write("roster.tsv", "index\tgenome\tassembly\ttaxid\tdomain\tphylum\tgff\tfasta\tloci\tgenome_bases\tcontigs\n"
			+"0\tg1\ta\tNA\tNA\tNA\tNA\t/recorded/one.fa\t4\t100\t1\n"
			+"1\tg2\tb\tNA\tNA\tNA\tNA\t/missing/two.fa\t1\t100\t1\n");
		panelText=CobalaminPanelPrep.HEADER+"\nsame\tg1\tΩ%seq\tRF00174\tRF00174\t1\t20\t2\t18\t+\n"
			+"same\tg2\tΩ%seq\tRF00174\tRF00174\t1\t20\t2\t18\t+\n"
			+"clan\tg1\tΩ%seq\tRF01482\tB12clade\t22\t35\t23\t34\t-\n"
			+"unresolved\tg1\tΩ%seq\tRF01689\tUNRESOLVED\t40\t50\t0\t0\t+\n"
			+"unreviewed\tg1\tΩ%seq\tRF00174\tUNREVIEWED\t60\t70\t0\t0\t-\n";
		write("panel.tsv", panelText);
		write("fmn.tsv", FMN_HEADER+"\nΩ%seq\t50\t80\t-\tACCEPTED_NOVEL_FMN\tmodel%0\t0.97\t31\t7\t42.5\t66.2\tBELOW_GA\towner\treviewer\t00000000000000000001\t00000000000000000002\n"
			+"Ω%seq\t81\t90\t+\tUNRESOLVED\tmodel\t0.78\tNA\tNA\t42.2\t66.2\tBELOW_GA\tretain_unresolved\treviewer\t00000000000000000001\t00000000000000000002\n");
		bindingText=FMN_BINDING_HEADER+"\n2\t/missing/two.fa\tΩ%seq\t50\t80\t-\tACCEPTED_NOVEL_FMN\n"
			+"3\t/missing/one.fa\tΩ%seq\t81\t90\t+\tUNRESOLVED\n";
		write("bindings.tsv", bindingText);
		write("manifest.tsv", MANIFEST_HEADER+'\n'+manifestRow("panel", "cobalamin_panel_v1", "panel.tsv", "roster.tsv", 5)
			+manifestRow("fmn", "fmn_adjudication_v1", "fmn.tsv", "bindings.tsv", 2));
	}

	private void badBinding(final String name, final String bindings, final String error){
		write(name+".tsv", bindings);
		write(name+".manifest.tsv", MANIFEST_HEADER+'\n'+manifestRow("fmn", "fmn_adjudication_v1", "fmn.tsv", name+".tsv", 2));
		expect(()->invoke(name+".manifest.tsv", "members.tsv", "bad_"+name), error);
		check(!new File(path("bad_"+name+".prior_summary.tsv")).exists(), "Failed binding must not publish completion");
	}
	private void badPanel(final String name, final String rows, final int expected, final String error){
		write(name+".tsv", rows);
		write(name+".manifest.tsv", MANIFEST_HEADER+'\n'+manifestRow("panel", "cobalamin_panel_v1", name+".tsv", "roster.tsv", expected));
		expect(()->invoke(name+".manifest.tsv", "members.tsv", "bad_"+name), error);
		check(!new File(path("bad_"+name+".prior_summary.tsv")).exists(), "Failed import must not publish completion");
	}
	private String manifestRow(final String id, final String format, final String decisions, final String bindings, final int expected){
		return id+'\t'+format+'\t'+path(decisions)+'\t'+path(bindings)+"\tfixture_authority\t"+expected+'\n';
	}
	private void invoke(final String manifest, final String members, final String output){
		BuildNcrnaReference.main(new String[]{"mode=prior", "in="+path("evidence"), "priors="+path(manifest),
			"census="+path("census.tsv"), "memberships="+path(members), "scope="+path("scope.tsv"),
			"assets="+path("assets.tsv"), "queries="+path("queries.tsv"), "out="+path(output), "t=1", "bgzip=f", "zipthreads=1"});
	}
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private void write(final String name, final String text){
		try(Output out=new Output(path(name))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
	}
	private String[] lines(final String name){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path(name))){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){result.add(new String(row, StandardCharsets.UTF_8));}
		}
		return result.toArray(new String[result.size()]);
	}
	private static String find(final String[] rows, final String prefix){
		String result=null;
		for(final String row:rows){
			if(!row.startsWith(prefix)){continue;}
			check(result==null, "Expected exactly one source row: "+prefix);
			result=row;
		}
		check(result!=null, "Missing source row: "+prefix);
		return result;
	}
	private static boolean contains(final String[] rows, final String expected){
		for(final String row:rows){if(row.equals(expected)){return true;}}
		return false;
	}
	private static void expect(final Runnable action, final String fragment){
		try{action.run();}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(fragment), "Unexpected failure: "+failure);
			return;
		}
		throw new AssertionError("Expected failure containing: "+fragment);
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final File root;
	private String membersText, panelText, bindingText;
	private static final String MANIFEST_HEADER="source_id\tformat\tdecisions\tbindings\tauthority\texpected_rows";
	private static final String FMN_BINDING_HEADER="decision_line\tfasta\tseqid\tstart1\tstop1\tstrand\tclassification";
	private static final String FMN_HEADER="seqid\tstart1\tstop1\tstrand\tclassification\twinning_model\tverifier_identity\tmodel_length\toriginal_annotated_members\tbest_CM_bits\tCM_GA_bits\tCM_status\tdecision_source\tindependent_reviewer\ttrace_packet_sha80\tCM_packet_sha80";
}
