package test.prot;

import prot.HbmTailLossAudit;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

/** Literal schema fixture for HbmTailLossAudit. @author Yelan */
public final class HbmTailLossAuditTest{

	public static void main(String[] args) throws Exception{
		final Path dir=Files.createTempDirectory("hbm-tail-loss-audit-test");
		final Path list=dir.resolve("list.tsv"), rows=dir.resolve("rows.tsv"), support=dir.resolve("support.tsv");
		final Path calls=dir.resolve("calls.gff"), ref=dir.resolve("ref.gff"), out=dir.resolve("out.tsv");
		write(rows,
			"query\tquery_source\tfamily_id\tvariant\tlength\tstatus\tscore64\tidentity\tpaired_fraction\tcore_coverage\tpaired\tpaired_core\tterminal_paired\tterminal_identical\tquery_terminal_paired\tquery_terminal_identical\tquery_suffix\tmodel_suffix\tquery_start\tmodel_start\tmodel_length\tpath\tshortlist_rank\tterminal20_score64\tterminal40_score64\tpaired_identity\tcounterfactual_pass\n"+
			"q1\tinput\t7\toriginal\t33\tstop_supported\t123\t45.5\t0.9\t0.8\t31\t30\t9\t2\t20\t8\t0\t1\t2\t3\t40\tmmmm\t4\t-5\t6\t47.7\tfalse\n"+
			"q1\tinput\t7\treplace20\t33\tstop_supported\t123\t45.5\t0.9\t0.8\t31\t30\t9\t2\t20\t8\t0\t1\t2\t3\t40\tmmmm\t4\t-5\t6\t47.7\tfalse\n");
		write(support,
			"query\tgenome\tcontig\tstart\tstop\tstrand\tfamily_id\trepresentative\ttail_length\tscore64\tidentity\tpaired_fraction\tcore_cov\tterminal_paired\tquery_suffix\treference_suffix\tquery_start\treference_start\n"+
			"q1\ttoy\tctg\t5\t12\t+\t7\trep\t8\t123\t45.5\t0.9\t0.8\t9\t0\t1\t2\t3\n");
		write(calls, "##gff-version 3\nctg\tBBTools\tCDS\t5\t12\t.\t+\t0\tID=call\n");
		write(ref, "##gff-version 3\nctg\tRef\tCDS\t5\t12\t.\t+\t0\tID=ref\n");
		write(list, "genome\trows\tsupport\tcalls\tref\ntoy\t"+rows+"\t"+support+"\t"+calls+"\t"+ref+"\n");
		HbmTailLossAudit.main(new String[]{"list="+list, "out="+out});
		final List<String> lines=Files.readAllLines(out, StandardCharsets.UTF_8);
		check(lines.size()==2, "Expected header plus one loss row");
		final String[] header=lines.get(0).split("\t", -1), row=lines.get(1).split("\t", -1);
		check(header.length==30, "Header must describe all 30 emitted fields: "+header.length);
		check(row.length==header.length, "Loss row width must match the header");
		check(row[26].equals("4") && row[27].equals("-5") && row[28].equals("6") && row[29].equals("47.7"),
			"Replay tail columns shifted out of order");
		check(row[7].equals("1") && row[8].equals("1") && row[9].equals("1") && row[10].equals("1") && row[11].equals("1"),
			"Expected exact selected and raw-reference stop matches");
		System.err.println("HBM_TAIL_LOSS_AUDIT_TEST_PASS");
	}

	private static void write(Path path, String s) throws Exception{Files.write(path, s.getBytes(StandardCharsets.UTF_8));}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private HbmTailLossAuditTest(){}
}
