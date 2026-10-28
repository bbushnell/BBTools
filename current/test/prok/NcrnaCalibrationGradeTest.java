package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.NcrnaCalibrationGrade;
import prok.NcrnaCalibrationPrep;
import prok.RrnaResourceIO;

/** Distinguishes duplicate FPs, unresolved overlaps and genuinely unmatched calls.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCalibrationGradeTest {
	public static void main(String[] args)throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Supply a fresh fixture directory");}
		final Path root=Paths.get(args[0]);Files.createDirectory(root);
		write(root.resolve("families.tsv"), "cobalamin\tcobalamin\n");
		write(root.resolve("source.tsv"), NcrnaCalibrationPrep.TRUTH_HEADER+"\nctg\tcobalamin\t10\t29\t+\nctg\tcobalamin\t50\t69\t+\nctg\tcobalamin\t100\t119\t+\n");
		write(root.resolve("full.tsv"), NcrnaCalibrationPrep.TRUTH_HEADER+"\nctg\tcobalamin\t10\t29\t+\nctg\tcobalamin\t45\t74\t+\n");
		write(root.resolve("unresolved.tsv"), NcrnaCalibrationPrep.TRUTH_HEADER+"\nctg\tcobalamin\t100\t119\t+\n");
		write(root.resolve("calls.gff"), call(10,29)+call(10,29)+call(100,119)+call(200,219)+call(45,74));
		NcrnaCalibrationGrade.main(new String[]{"in="+root.resolve("calls.gff"), "source="+root.resolve("source.tsv"), "full="+root.resolve("full.tsv"),
			"unresolved="+root.resolve("unresolved.tsv"), "familymap="+root.resolve("families.tsv"), "outprefix="+root.resolve("result")});
		checkSummary(root.resolve("result.source.summary.tsv"), 3, 2, 0);
		checkSummary(root.resolve("result.full.summary.tsv"), 2, 2, 0);
		final ByteFile in=RrnaResourceIO.open(root.resolve("result.accounting.tsv").toString());final LineParser1 p=new LineParser1('\t');
		try{
			in.nextLine();p.set(in.nextLine());
			check(p.parseInt(0)==5 && p.parseInt(1)==4 && p.parseInt(2)==1 && p.parseInt(3)==1,
				"Unresolved calls stay unscored; duplicates remain FPs without redundant CM work; the novel physical call needs review");
		}finally{check(!in.close(), "Accounting read failed");}
		write(root.resolve("adjacent.tsv"), NcrnaCalibrationPrep.TRUTH_HEADER+"\nctg\tcobalamin\t10\t29\t+\nctg\tcobalamin\t30\t49\t+\n");
		write(root.resolve("empty.tsv"), NcrnaCalibrationPrep.TRUTH_HEADER+'\n');
		write(root.resolve("duplicate_calls.gff"), call(10,49)+call(10,49));
		NcrnaCalibrationGrade.main(new String[]{"in="+root.resolve("duplicate_calls.gff"), "source="+root.resolve("adjacent.tsv"), "full="+root.resolve("adjacent.tsv"),
			"unresolved="+root.resolve("empty.tsv"), "familymap="+root.resolve("families.tsv"), "outprefix="+root.resolve("adjacent_result")});
		checkSummary(root.resolve("adjacent_result.source.summary.tsv"), 1, 1, 1);
		checkSummary(root.resolve("adjacent_result.full.summary.tsv"), 1, 1, 1);
		System.out.println("NCRNA_CALIBRATION_GRADE_TEST_PASS");
	}
	private static String call(final int start, final int stop){return "ctg\tBBTools\tRNA\t"+start+'\t'+stop+"\t1\t+\t0\tRNA,model:cobalamin_0\n";}
	private static void checkSummary(final Path path, final int tp, final int fp, final int fn){
		final ByteFile in=RrnaResourceIO.open(path.toString());final LineParser1 p=new LineParser1('\t');
		try{
			in.nextLine();p.set(in.nextLine());check(p.termEquals("cobalamin",0) && p.parseInt(1)==tp && p.parseInt(2)==fn && p.parseInt(3)==fp,
				"One-to-one grade must retain duplicate FPs and the explicit view denominator");
		}finally{check(!in.close(), "Grade summary read failed");}
	}
	private static void write(final Path path, final String text){
		final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());try{out.print(text);}finally{RrnaResourceIO.finish(out);}
	}
	private static void check(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
}
