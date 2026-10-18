package test.prot;

import prot.HbmTailChangeAudit;

import java.io.File;
import java.nio.file.Files;

/** Literal fixtures for stop-event to normalized-row joins. @author Yelan */
public final class HbmTailChangeAuditTest{
	public static void main(String[] args) throws Exception{
		final File dir=Files.createTempDirectory("hbm-tail-change-audit-test").toFile();
		final File ref=new File(dir, "reference.tsv"), base=new File(dir, "baseline.tsv");
		final File cand=new File(dir, "candidate.tsv"), stops=new File(dir, "stops.tsv");
		Files.write(ref.toPath(), (ROW_HEADER
			+"stops-v1\t10\tCDS\tctg\t0\t100\t150\t0\t11\t0\t0\tprimary\tID=refA\n"
			+"stops-v1\t11\tCDS\tctg\t1\t110\t140\t0\t11\t0\t0\tprimary\tID=lost\n"
			+"stops-v1\t12\tpseudogene\tctg\t0\t50\t150\t-1\t-1\t0\t0\tpseudo\tID=excludedGain\n"
			+"stops-v1\t0\tpseudogene\tctg\t0\t500\t550\t-1\t-1\t0\t0\tpseudo_mask\tID=maskA\n"
			+"stops-v1\t0\tpseudogene\tctg\t0\t600\t650\t-1\t-1\t0\t0\tpseudo_mask\tID=maskB\n").getBytes());
		Files.write(base.toPath(), (ROW_HEADER
			+"stops-v1\t20\tCDS\tctg\t1\t110\t140\t0\t11\t0\t0\tprimary\t"+attrs("20", 100, 0, 40, 10, ".")+"\n"
			+"stops-v1\t21\tCDS\tctg\t0\t300\t350\t0\t11\t0\t0\tprimary\t"+attrs("21", 200, 50, 41, 20, "77")+"\n"
			+"stops-v1\t22\tpseudogene\tctg\t1\t110\t160\t-1\t-1\t0\t0\tpseudo\t"+attrs("22", 100, 0, 40, 10, ".")+"\n").getBytes());
		Files.write(cand.toPath(), (ROW_HEADER
			+"stops-v1\t30\tCDS\tctg\t0\t50\t150\t0\t11\t0\t0\tprimary\t"+attrs("30", 700, 101, 42, 30, "88")+"\n"
			+"stops-v1\t31\tCDS\tctg\t0\t320\t350\t0\t11\t0\t0\tprimary\t"+attrs("31", 300, 31, 43, 31, "99")+"\n").getBytes());
		Files.write(stops.toPath(), (STOP_HEADER
			+"stops-v1\tctg\t0\t150\t2\t1\t0\t0\t1\t1\tgained\t10,12\t.\t30\n"
			+"stops-v1\tctg\t1\t110\t1\t1\t2\t1\t0\t0\tlost\t11\t20,22\t.\n"
			+"stops-v1\tctg\t0\t350\t0\t0\t1\t1\t1\t1\tsame_stop_changed_starts\t.\t21\t31\n").getBytes());
		final File out=new File(dir, "out.tsv");
		HbmTailChangeAudit.run(new String[]{"stops="+stops, "reference="+ref, "baseline="+base, "candidate="+cand, "out="+out});
		final String text=new String(Files.readAllBytes(out.toPath()));
		check(text.contains("\tgained\tctg\t0\t150\t") && text.contains("\tcandidate\t30\t1\t1\t1\t20\t50\t150\t101\t"),
			"Opposite-strand overextension and stop-only gain must be explicit");
		check(text.contains("\tsame_stop_changed_starts\tctg\t0\t350\t") && text.contains("\tcandidate\t31\t0\t0\t0\t.\t320\t350\t31\t"),
			"Start-only row links must retain baseline and candidate rows");
		final File bad=new File(dir, "bad.tsv");
		Files.write(bad.toPath(), (STOP_HEADER+"stops-v1\tctg\t0\t150\t1\t1\t0\t0\t1\t1\tgained\t10\t.\t999\n").getBytes());
		boolean rejected=false;
		try{HbmTailChangeAudit.run(new String[]{"stops="+bad, "reference="+ref, "baseline="+base, "candidate="+cand, "out="+new File(dir, "bad.out")});}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Missing candidate row IDs must fail");
		Files.write(bad.toPath(), ("schema\tinput_row\ttype\tcontig\tstrand\tstart\tend\tphase\ttable\tpartial5\tpartial3\tattributes\n").getBytes());
		rejected=false;
		try{HbmTailChangeAudit.run(new String[]{"stops="+stops, "reference="+bad, "baseline="+base, "candidate="+cand, "out="+new File(dir, "bad2.out")});}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Normalized row schema drift must fail");
		Files.write(bad.toPath(), (STOP_HEADER+"stops-v1\tctg\t0\t150\t1\t1\t0\t0\t1\t1\tgained\t10\t.\t31\n").getBytes());
		rejected=false;
		try{HbmTailChangeAudit.run(new String[]{"stops="+bad, "reference="+ref, "baseline="+base, "candidate="+cand, "out="+new File(dir, "bad3.out")});}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Existing IDs with the wrong event key must fail");
		Files.write(bad.toPath(), (STOP_HEADER+"stops-v1\tctg\t0\t150\t2\t2\t0\t0\t1\t1\tgained\t10,10\t.\t30\n").getBytes());
		rejected=false;
		try{HbmTailChangeAudit.run(new String[]{"stops="+bad, "reference="+ref, "baseline="+base, "candidate="+cand, "out="+new File(dir, "bad4.out")});}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Duplicate row links must fail");
		final File empty=new File(dir, "empty.out");
		Files.write(bad.toPath(), (STOP_HEADER+"stops-v1\tctg\t0\t150\t1\t1\t1\t1\t1\t1\tunchanged\t10\t20\t30\n").getBytes());
		HbmTailChangeAudit.run(new String[]{"stops="+bad, "reference="+ref, "baseline="+base, "candidate="+cand, "out="+empty});
		check(new String(Files.readAllBytes(empty.toPath())).trim().split("\n").length==1, "Empty changed inputs must emit a header");
		System.err.println("HBM_TAIL_CHANGE_AUDIT_TEST_PASS");
	}
	private static String attrs(String fam, int score, int paired, int term, int best, String terminalScore){
		return "bbgap_status=stop_supported;bbgap_native_score=1.5;bbgap_call_bonus=50"
			+";bbgap_rule=hbm_stop_tail_profile_joint_v1;bbgap_family_id="+fam
			+";bbgap_score64="+score+";bbgap_tail_paired="+paired+";bbgap_identity=40.5"
			+";bbgap_tail_terminal_score64="+terminalScore+";bbgap_best_other_score64="+best
			+";bbgap_runner_supported_score64=-1;bbgap_tail_terminal_paired=10"
			+";bbgap_tail_query_suffix=0;bbgap_tail_model_suffix=0;bbgap_tail_query_start=0"
			+";bbgap_tail_model_start=0;bbgap_label_class=fixture;bbgap_label_support=fixture"
			+";bbgap_product=fixture";
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private static final String ROW_HEADER="schema\tinput_row\ttype\tcontig\tstrand\tstart\tend\tphase\ttable\tpartial5\tpartial3\tstratum\tattributes\n";
	private static final String STOP_HEADER="schema\tcontig\tstrand\tterminal\treference_raw_rows\treference_primary_rows\tbaseline_raw_rows\tbaseline_primary_rows\tcandidate_raw_rows\tcandidate_primary_rows\tchange\treference_rows\tbaseline_rows\tcandidate_rows\n";
	private HbmTailChangeAuditTest(){}
}
