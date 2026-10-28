package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.CobalaminCallTransitions;
import prok.CobalaminPanelPrep;
import prok.GffCallSetDelta;
import prok.NcrnaReviewUnion;
import prok.RrnaResourceIO;

/** Exact survival differs from overlap; assignment recovery differs from relabeling;
 * every removed non-B12 gene row remains visible, including unexplained changes.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminCallTransitionsTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh fixture directory");final Path root=Paths.get(args[0]);Files.createDirectory(root);
		write(root, "contigs.tsv", CobalaminPanelPrep.CONTIG_HEADER+"\ng\tassembly\tseqA\t2000\ng\tassembly\tseqB\t2000\n");
		write(root, "panel.tsv", CobalaminPanelPrep.HEADER+"\nt1\tg\tseqA\tRF00174\tRF00174\t100\t199\t100\t199\t+\n"
			+"t2\tg\tseqA\tRF00174\tRF00174\t400\t499\t400\t499\t+\nt3\tg\tseqB\tRF00174\tRF00174\t100\t199\t100\t199\t+\n");
		final String before=NcrnaReviewUnion.PANEL_CALL_HEADER+'\n'+call("old_tp1", "RF01482", 100, 199, "TP", "t1")
			+call("old_tp2", "RF01689", 400, 499, "TP", "t2")+call("fp_exact", "RF01482", 700, 799, "APPARENT_FP", "NA")
			+call("fp_overlap", "RF01689", 1000, 1099, "APPARENT_FP", "NA")+call("fp_absent", "RF01482", 1500, 1599, "APPARENT_FP", "NA");
		final String after=NcrnaReviewUnion.PANEL_CALL_HEADER+'\n'+call("new_tp1", "RF00174", 100, 199, "TP", "t1")
			+call("new_fp", "RF01689", 700, 799, "APPARENT_FP", "NA")+call("shifted_fp", "RF01689", 1050, 1149, "APPARENT_FP", "NA");
		write(root, "before.tsv", before);write(root, "after.tsv", after);
		write(root, "changes.tsv", GffCallSetDelta.CHANGE_HEADER+'\n'
			+"REMOVED\tseqA\tBBTools\tCDS\t600\t899\t10\t+\t0\tCDS\n"
			+"ADDED\tseqA\tBBTools\tCDS\t801\t899\t9\t+\t0\tCDS\n"
			+"REMOVED\tseqA\tBBTools\tCDS\t150\t170\t2\t-\t0\tCDS\n"
			+"REMOVED\tseqA\tBBTools\ttRNA\t1600\t1669\t2\t+\t.\ttRNA\n"
			+"REMOVED\tseqA\tBBTools\tRNA\t400\t499\t1\t+\t.\tRNA,family:cobalamin_rf01689\n");
		run(root, "after.tsv", "valid", 3);
		check(value(root.resolve("valid.summary.tsv"), "before_FP_exact_present")==1
			&& value(root.resolve("valid.summary.tsv"), "before_FP_overlap_only")==1
			&& value(root.resolve("valid.summary.tsv"), "before_FP_absent")==1,
			"A shifted overlapping call is not an exact survivor and an absent call is retained explicitly");
		check(value(root.resolve("valid.summary.tsv"), "recovered_by_RF00174")==1
			&& value(root.resolve("valid.summary.tsv"), "lost_native_assignment")==1,
			"Recovery follows the same native truth assignment, not a nearby call or the emitted label alone");
		check(contains(root.resolve("valid.fp.tsv"), "fp_exact\tg\tseqA\t700\t799\t+\tRF01482\tEXACT_PRESENT\tnew_fp\tRF01689"),
			"An exact physical survivor keeps its changed emitted subtype visible");
		check(contains(root.resolve("valid.tp.tsv"), "old_tp1\tg\tt1\tRF00174\tRF00174\tRF01482\tRETAINED_CHANGED_LABEL\tnew_tp1\tRF00174\ttrue")
			&& contains(root.resolve("valid.tp.tsv"), "old_tp2\tg\tt2\tRF00174\tRF00174\tRF01689\tLOST_NATIVE_ASSIGNMENT\tNA\tNA\tNA"),
			"Changed emitted labels preserve truth identity; unassigned original TPs remain explicit");
		check(contains(root.resolve("valid.labels.tsv"), "RF01482\t1\t2\t0\t0\t0\t0\t0\t1\t0\t1\t0\t1\t0\t1")
			&& contains(root.resolve("valid.labels.tsv"), "RF01689\t1\t1\t0\t0\t2\t0\t0\t0\t1\t0\t1\t0\t1\t0"),
			"Original emitting labels own losses; a recovered RF00174 call must not be counted as a lost sibling TP");
		check(contains(root.resolve("valid.gene_links.tsv"), "2\tCDS\tg\tnew_fp\tRF01689\tAPPARENT_FP\tSAME\t100\tREMOVED")
			&& contains(root.resolve("valid.gene_links.tsv"), "4\tCDS\tg\tnew_tp1\tRF00174\tTP\tOPPOSITE\t21\tREMOVED")
			&& contains(root.resolve("valid.gene_links.tsv"), "5\ttRNA\tg\tNA\tNA\tNA\tNO_B12_OVERLAP\t0\tREMOVED"),
			"Gene-row links retain TP/FP status, strand relation and removals with no direct B12 overlap");
		write(root, "wrong.tsv", after.replace("\tTP\tt1\t", "\tTP\tt3\t"));
		boolean failed=false;try{run(root, "wrong.tsv", "wrong", 3);}catch(IllegalArgumentException expected){failed=true;}
		check(failed && !Files.exists(root.resolve("wrong.summary.tsv")), "A truth ID on a different contig cannot fabricate recovery");
		write(root, "duplicate.tsv", after+call("duplicate", "RF00174", 700, 799, "APPARENT_FP", "NA"));
		failed=false;try{run(root, "duplicate.tsv", "duplicate", 4);}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Different labels at an exact duplicate physical location require explicit event accounting");
		write(root, "wrong_model.tsv", after.replace("cobalamin_RF00174_fixture", "cobalamin_RF01482_fixture"));
		failed=false;try{run(root, "wrong_model.tsv", "wrong_model", 3);}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "A claimed RF00174 recovery must actually name an RF00174 model");
		System.out.println("COBALAMIN_CALL_TRANSITIONS_TEST_PASS");
	}
	private static String call(String id, String family, int start, int stop, String outcome, String truth){
		return "FULL_CM\tunion\t"+id+"\tg\tseqA\t"+family+'\t'+start+'\t'+stop+"\t+\tcobalamin_"+family+"_fixture\t"+outcome+'\t'+truth+"\tNA\tNA\n";
	}
	private static void run(Path root, String after, String out, int count){
		CobalaminCallTransitions.main(new String[]{"before="+root.resolve("before.tsv"), "after="+root.resolve(after), "panel="+root.resolve("panel.tsv"),
			"contigs="+root.resolve("contigs.tsv"), "changes="+root.resolve("changes.tsv"), "outprefix="+root.resolve(out), "beforecount=5", "aftercount="+count});
	}
	private static long value(Path file, String key){
		final ByteFile in=RrnaResourceIO.open(file.toString());final LineParser1 p=new LineParser1('\t');long result=-1;
		try{in.nextLine();for(byte[] row; (row=in.nextLine())!=null; ){p.set(row);if(p.termEquals(key, 0)){check(result<0, "Summary metrics must be unique");result=p.parseLong(1);}}}
		finally{check(!in.close(), "Fixture summary read failed");}check(result>=0, "Missing summary metric: "+key);return result;
	}
	private static boolean contains(Path file, String prefix){
		final ByteFile in=RrnaResourceIO.open(file.toString());boolean found=false;
		try{for(byte[] row; (row=in.nextLine())!=null; ){if(new String(row, java.nio.charset.StandardCharsets.UTF_8).startsWith(prefix)){found=true;}}}
		finally{check(!in.close(), "Fixture detail read failed");}return found;
	}
	private static void write(Path root, String name, String text){final ByteStreamWriter out=RrnaResourceIO.writer(root.resolve(name).toString());try{out.print(text);}finally{RrnaResourceIO.finish(out);}}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
