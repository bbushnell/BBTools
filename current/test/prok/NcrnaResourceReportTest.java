package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import consensus.BaseGraph;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.NcrnaResourceReport;
import prok.RrnaResourceIO;
import prok.TrnaConsensusBuilder;
import structures.ByteBuilder;

/** Checks HBM/consensus identity and full accounting even for unassigned members.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaResourceReportTest {
	public static void main(String[] args)throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Supply a fresh fixture directory");}
		final Path root=Paths.get(args[0]);Files.createDirectory(root);
		final String a="ACGTACGTACGTACGTACGT", c="CCCCCCCCCCCCCCCCCCCC";
		write(root.resolve("ref.fa"), ">cobalamin_0 n=2\n"+a+"\n>cobalamin_1 n=1\n"+c+'\n');
		write(root.resolve("input.fa"), ">q0 tid_1\n"+a+"\n>q1 tid_2\n"+c+"\n>q2 tid_3\n"+a+'\n');
		write(root.resolve("assignments.tsv"), "#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus\n"+
			"q0 tid_1\t1\tunknown\tc0\tcobalamin_0 n=2\tkept\nq1 tid_2\t2\tunknown\tc1\tcobalamin_1 n=1\tkept\nq2 tid_3\t3\tunknown\torphan_q2\t\torphan\n");
		final ArrayList<BaseGraph> graphs=new ArrayList<BaseGraph>();
		graphs.add(new BaseGraph("cobalamin_0 n=2", a.getBytes(StandardCharsets.US_ASCII), null, 0, 0));
		graphs.add(new BaseGraph("cobalamin_1 n=1", c.getBytes(StandardCharsets.US_ASCII), null, 1, 0));
		TrnaConsensusBuilder.writeTextModels(graphs, root.resolve("models.hbm").toString());
		run(root, "good", "models.hbm");
		final ByteFile in=RrnaResourceIO.open(root.resolve("good.summary.tsv").toString());final LineParser1 p=new LineParser1('\t');
		try{
			in.nextLine();p.set(in.nextLine());
			check(p.parseInt(0)==3 && p.parseInt(1)==2 && p.parseInt(2)==2 && p.parseInt(4)==1 && p.parseInt(6)==3 && p.parseInt(13)==6,
				"All three members, including the orphan, need exhaustive evidence against both models");
			check(p.parseFloat(8)==1 && p.parseFloat(12)==1, "Exact members must score identity1 independent of builder assignment");
		}finally{check(!in.close(), "Report summary read failed");}
		checkParallel(root, a, c);
		graphs.set(1, new BaseGraph("cobalamin_1 n=1", a.getBytes(StandardCharsets.US_ASCII), null, 1, 0));
		TrnaConsensusBuilder.writeTextModels(graphs, root.resolve("mismatch.hbm").toString());
		boolean failed=false;try{run(root, "bad", "mismatch.hbm");}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Equal model counts/names cannot conceal changed reference bases");
		System.out.println("NCRNA_RESOURCE_REPORT_TEST_PASS");
	}
	/** Crosses several worker batches and checks ordering and dispositions. */
	private static void checkParallel(final Path root, final String a, final String c){
		final ByteBuilder fasta=new ByteBuilder(), assignments=new ByteBuilder();
		assignments.append("#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus\n");
		for(int i=0; i<129; i++){
			fasta.append(">m").append(i).append(" tid_1\n").append(i%2==0?a:c).nl();
			assignments.append('m').append(i).append(" tid_1\t1\tunknown\tc\t");
			if(i%3==0){assignments.append("cobalamin_0 n=2\tkept\n");}
			else{assignments.tab().append(i%3==1?"small":"orphan").nl();}
		}
		write(root.resolve("parallel.fa"), fasta.toString());write(root.resolve("parallel.tsv"), assignments.toString());
		for(int threads:new int[]{1, 16}){
			NcrnaResourceReport.main(new String[]{"in="+root.resolve("parallel.fa"), "ref="+root.resolve("ref.fa"),
				"hbm="+root.resolve("models.hbm"), "assignments="+root.resolve("parallel.tsv"),
				"outprefix="+root.resolve("t"+threads), "expected=129", "reportcut=0.6", "t="+threads});
		}
		for(String suffix:new String[]{".members.tsv.gz", ".summary.tsv", ".models.tsv", ".histogram.tsv"}){
			final ByteFile one=RrnaResourceIO.open(root.resolve("t1"+suffix).toString()), many=RrnaResourceIO.open(root.resolve("t16"+suffix).toString());
			try{
				for(byte[] line; (line=one.nextLine())!=null; ){check(Arrays.equals(line, many.nextLine()), "Thread count changed ordered report: "+suffix);}
				check(many.nextLine()==null, "Parallel report contains extra rows: "+suffix);
			}finally{check(!one.close() && !many.close(), "Parallel parity input failed");}
		}
		// Duplicate input must fail before any output is opened, rather than hiding a missing member.
		write(root.resolve("duplicate.fa"), ">q0 tid_1\n"+a+"\n>q0 tid_1\n"+a+"\n>q2 tid_3\n"+a+'\n');
		boolean failed=false;
		try{
			NcrnaResourceReport.main(new String[]{"in="+root.resolve("duplicate.fa"), "ref="+root.resolve("ref.fa"),
				"hbm="+root.resolve("models.hbm"), "assignments="+root.resolve("assignments.tsv"),
				"outprefix="+root.resolve("duplicate"), "expected=3", "reportcut=0.6", "t=16"});
		}catch(IllegalArgumentException expected){failed=true;}
		check(failed && !Files.exists(root.resolve("duplicate.members.tsv.gz")), "Duplicate members must fail before report publication");
		// Scrabble's 21-bit coordinate limit fails inside a worker, not during input validation.
		final char[] oversized=new char[1<<21];Arrays.fill(oversized, 'A');
		write(root.resolve("oversized.fa"), ">q0 tid_1\n"+new String(oversized)+'\n');
		write(root.resolve("one.tsv"), "#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus\n"+
			"q0 tid_1\t1\tunknown\tc0\tcobalamin_0 n=2\tkept\n");
		failed=false;
		try{
			NcrnaResourceReport.main(new String[]{"in="+root.resolve("oversized.fa"), "ref="+root.resolve("ref.fa"),
				"hbm="+root.resolve("models.hbm"), "assignments="+root.resolve("one.tsv"),
				"outprefix="+root.resolve("worker_failure"), "expected=1", "reportcut=0.6", "t=16"});
		}catch(IllegalStateException expected){failed=expected.getCause() instanceof AssertionError;}
		check(failed && !Files.exists(root.resolve("worker_failure.members.tsv.gz")), "Worker assertions must reach the caller without publishing an incomplete report");
	}
	private static void run(final Path root, final String prefix, final String hbm){
		NcrnaResourceReport.main(new String[]{"in="+root.resolve("input.fa"), "ref="+root.resolve("ref.fa"), "hbm="+root.resolve(hbm),
			"assignments="+root.resolve("assignments.tsv"), "outprefix="+root.resolve(prefix), "expected=3", "reportcut=0.6"});
	}
	private static void write(final Path path, final String text){
		final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());try{out.print(text);}finally{RrnaResourceIO.finish(out);}
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
}
