package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import prok.BuildNcrnaReference;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaReferenceEvidence;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;

/** Saved-table fixtures with literal arithmetic oracles and nonexistent payloads.
 * Invoked by BuildNcrnaReferenceTest through its ordinary native launcher.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceProfilesTest{

	NcrnaReferenceProfilesTest(final File root_){root=root_;}

	void run(){
		check(root.mkdir(), "Profile fixture directory must be new");
		prepare();
		writeEvidence("good", SUMMARY, COVERAGE, SOURCES);
		invoke("good", "result");
		final String[] models=lines("result.model_profile.tsv");
		check(models.length==4, "All three models, including the zero-hit model, need a profile");
		check(models[1].equals("RF00050\tFMN\tFMN\tTARGET_CLASS\tfixture\t10\t5\t5\t4\t1\t0\t0\t3\t0\t5\t10\t10\t10\taptamer_only"),
			"Weighted spans [5,8,10,10,10] must give p05=5, median=10, p95=10, mode=10");
		check(models[2].equals("RF00162\tSAM-I\tSAM_I\tTARGET_CLASS\tfixture\t6\t5\t3\t2\t0\t0\t1\t2\t0\t3\t6\t6\t6\taptamer_only"),
			"Native both-end truncation stays separate from span and weighted quantiles");
		check(models[3].equals("RF00001\t5S\t5S\tTARGET_CLASS\tfixture\t12\t5\t0\t0\t0\t0\t0\t0\t0\tNA\tNA\tNA\tNA\tcompartment_unresolved"),
			"Zero-hit models have explicit zero counts and no invented quantiles");
		final String[] categories=lines("result.availability.tsv");
		check(categories.length==4 && categories[1].equals("NO_AVAILABLE_GFF\t2\t2\t0\t2\t0\t0")
			&& categories[2].equals("AVAILABLE_ONLY\t2\t2\t2\t0\t1\t1")
			&& categories[3].equals("MIXED_AVAILABLE_UNAVAILABLE\t1\t2\t1\t1\t0\t1"),
			"Three unavailable claims belong to two wholly unavailable inputs plus one mixed input");
		final String[] inputs=lines("result.input_availability.tsv");
		check(inputs.length==6 && inputs[3].equals("3\t1\t1\t1\t0\t8\t2\t6\tMIXED_AVAILABLE_UNAVAILABLE")
			&& inputs[4].equals("4\t1\t0\t0\t1\t0\t0\t0\tAVAILABLE_ONLY"),
			"Duplicate memberships do not duplicate claims; an empty parsed file remains available");
		check(lines("result.profile_summary.tsv")[1].equals("5\t5\t3\t8\t6\t3\t18\t5\t13\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\tPROFILE_COMPLETE_NOT_ADJUDICATED"),
			"Profile accounting must conserve the saved evidence without admitting genes");
		writeEvidence("sample", "3\t5\t300\t3\t8\t4\t2\t18\t5\t13\t0\t0\t3\t9\t6\tSAMPLE_EVIDENCE_NOT_ADJUDICATED\n",
			COVERAGE.replace("\t5\n", "\t3\n"), SOURCES.substring(0, SOURCES.indexOf("4\t/source/empty")));
		invoke("sample", "sample_result");
		check(lines("sample_result.profile_summary.tsv")[1].equals("3\t5\t3\t8\t4\t2\t18\t5\t13\tSAMPLE_EVIDENCE_NOT_ADJUDICATED\tSAMPLE_PROFILE_NOT_ADJUDICATED"),
			"A three-input sample of five catalog inputs must retain both denominator and explicit sample status");
		failure("duplicate_bin", SUMMARY, COVERAGE+FIRST_BIN, SOURCES, "Repeated accession/span/truncation");
		failure("wrong_total", SUMMARY, COVERAGE.replace("\t1\t5\t5\n", "\t1\t6\t5\n"), SOURCES, "Inconsistent accession hit denominator");
		failure("missing_zero", SUMMARY, COVERAGE.replace(ZERO_BIN, ""), SOURCES, "Missing queried-model profile");
		failure("repeat_zero", SUMMARY, COVERAGE+ZERO_BIN, SOURCES, "Only one explicit zero-hit row");
		failure("unknown_claim", SUMMARY, COVERAGE, SOURCES.replace("/source/a", "/source/wrong"), "Source claim is unknown");
		failure("wrong_memberships", SUMMARY, COVERAGE, SOURCES.replace("PARSED_GFF\t2\t8", "PARSED_GFF\t1\t8"), "different membership count");
		failure("missing_claim", SUMMARY, COVERAGE, SOURCES.replace(MIXED_MISSING, ""), "retain all original GFF claims");
		failure("unknown_input", SUMMARY, COVERAGE, SOURCES.replace("5\t/source/bad", "6\t/source/bad"), "unknown catalog input");
		failure("wrong_population", SUMMARY.replace("5\t5\t500", "4\t5\t500"), COVERAGE, SOURCES, "all queried-model completions");
		failure("wrong_bases", SUMMARY.replace("\t500\t", "\t501\t"), COVERAGE, SOURCES, "input dimensions disagree");
		failure("false_empty", SUMMARY, COVERAGE, SOURCES.replace("PARSED_GFF\t1\t10", "EMPTY_GFF\t1\t10"), "empty GFF cannot");
		failure("unknown_state", SUMMARY, COVERAGE, SOURCES.replace("UNREADABLE_GFF", "UNRECOGNIZED_GFF"), "Unavailable GFF state");
		// Equal span counts must select the smaller mode, regardless of row order.
		final String tied=COVERAGE.replace("10\t0\t3\t5\t5", "10\t0\t1\t3\t5")
			.replace("\t1\t5\t5\n", "\t1\t3\t5\n");
		writeEvidence("tied", SUMMARY.replace("500\t5\t8", "500\t5\t6"), tied, SOURCES);
		invoke("tied", "tied_result");
		check(lines("tied_result.model_profile.tsv")[1].endsWith("\t5\t8\t10\t5\taptamer_only"),
			"For spans [5,8,10], nearest-rank median=8 and smallest tied mode=5");
	}

	private void prepare(){
		final StringBuilder census=new StringBuilder(NcrnaEvidenceCatalog.CENSUS_HEADER).append('\n');
		final StringBuilder members=new StringBuilder(NcrnaEvidenceCatalog.MEMBERSHIP_HEADER).append('\n');
		final String[] gffs={"/source/a", "NA", "/source/Ω%file", "/source/empty", "/source/bad"};
		for(int i=1; i<=5; i++){
			final String fasta=path("absent"+i+".fa");
			check(!new File(fasta).exists(), "Profile tests must not have readable sequence payloads");
			census.append(i).append('\t').append(fasta).append("\t110\t100\t1\n");
			members.append("main\t").append(i).append('\t').append(i).append('\t').append(fasta)
				.append("\tNA\tNA\tNA\tUNKNOWN\t").append(gffs[i-1]).append("\tNA\n");
		}
		members.append("extra\t1\t3\t").append(path("absent3.fa")).append("\tNA\tNA\tNA\tUNKNOWN\tNA\tNA\n");
		members.append("extra\t2\t3\t").append(path("absent3.fa")).append("\tNA\tNA\tNA\tUNKNOWN\t/source/Ω%file\tNA\n");
		write("census.tsv", census.toString());
		write("members.tsv", members.toString());
		write("scope.tsv", NcrnaReferenceEvidence.SCOPE_HEADER+'\n'
			+"RF00050\tFMN\tprimary\taptamer_only\tFMN\tTARGET_CLASS\n"
			+"RF00162\tSAM-I\tprimary\taptamer_only\tSAM_I\tTARGET_CLASS\n"
			+"RF00001\t5S\tprimary\tcompartment_unresolved\t5S\tTARGET_CLASS\n");
		write("assets.tsv", NcrnaReferenceEvidence.ASSET_HEADER+'\n'
			+"RF00050\tfixture\t10\t20\t00000000000000000001\t/absent/fmn.cm\n"
			+"RF00162\tfixture\t6\t30\t00000000000000000002\t/absent/sam1.cm\n"
			+"RF00001\tfixture\t12\t20\t00000000000000000003\t/absent/5s.cm\n");
		write("queries.tsv", NcrnaReferenceEvidence.QUERY_HEADER+"\nRF00050\tfmn\t10\nRF00162\tsam1\t6\nRF00001\t5s\t12\n");
	}

	/** Each invalid fixture is independent and must fail before any report is opened. */
	private void failure(final String name, final String summary, final String coverage, final String sources, final String reason){
		writeEvidence(name, summary, coverage, sources);
		boolean failed=false;
		try{invoke(name, name+"_result");}
		catch(IllegalArgumentException expected){
			check(expected.getMessage().contains(reason), "Unexpected failure for "+name+": "+expected);
			failed=true;
		}
		check(failed, "Expected invalid fixture to fail: "+name);
		check(!new File(path(name+"_result.model_profile.tsv")).exists(), "Invalid saved data must not publish profiles");
	}
	private void writeEvidence(final String name, final String summary, final String coverage, final String sources){
		write(name+".summary.tsv", BuildNcrnaReference.SUMMARY_HEADER+'\n'+summary);
		write(name+".coverage.tsv", BuildNcrnaReference.COVERAGE_HEADER+'\n'+coverage);
		write(name+".sources.tsv.gz", BuildNcrnaReference.SOURCE_HEADER+'\n'+sources);
	}
	private void invoke(final String input, final String output){
		BuildNcrnaReference.main(new String[]{"mode=profile", "in="+path(input), "out="+path(output),
			"census="+path("census.tsv"), "memberships="+path("members.tsv"), "scope="+path("scope.tsv"),
			"assets="+path("assets.tsv"), "queries="+path("queries.tsv"), "t=1", "bgzip=f", "zipthreads=1"});
	}
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
	private String path(final String name){return new File(root, name).getAbsolutePath();}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private final File root;
	private static final String SUMMARY="5\t5\t500\t5\t8\t6\t3\t18\t5\t13\t0\t0\t3\t15\t12\tEVIDENCE_COMPLETE_NOT_ADJUDICATED\n";
	private static final String FIRST_BIN="RF00050\tfixture\t10\t5\t0\t1\t5\t5\n";
	private static final String ZERO_BIN="RF00001\tfixture\t12\tNA\tNA\t0\t0\t5\n";
	private static final String COVERAGE=FIRST_BIN+"RF00050\tfixture\t10\t8\t1\t1\t5\t5\n"
		+"RF00050\tfixture\t10\t10\t0\t3\t5\t5\n"
		+"RF00162\tfixture\t6\t3\t3\t1\t3\t5\nRF00162\tfixture\t6\t6\t0\t2\t3\t5\n"+ZERO_BIN;
	private static final String MIXED_MISSING="3\tNA\tMISSING_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n";
	private static final String SOURCES="1\t/source/a\tPARSED_GFF\t1\t10\t3\t7\t7\t0\t0\t0\t0\n"
		+"2\tNA\tMISSING_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n"
		+"3\t/source/Ω%25file\tPARSED_GFF\t2\t8\t2\t6\t6\t0\t0\t0\t0\n"+MIXED_MISSING
		+"4\t/source/empty\tEMPTY_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n"
		+"5\t/source/bad\tUNREADABLE_GFF\t1\t0\t0\t0\t0\t0\t0\t0\t0\n";
}
