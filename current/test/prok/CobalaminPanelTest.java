package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.CobalaminPanelGrade;
import prok.CobalaminPanelPrep;
import prok.NcrnaCalibrationPrep;
import prok.NcrnaClanOverlay;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Tests cross-subtype matches, clan-only truth, unresolved masks and duplicate calls.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminPanelTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh fixture directory");final Path root=Paths.get(args[0]);Files.createDirectory(root);
		fasta(root.resolve("a.fa"), "c1", 1000);fasta(root.resolve("b.fa"), "c2", 300);
		write(root.resolve("roster.tsv"), NcrnaCalibrationPrep.ROSTER_HEADER,
			"0\tg1\ta1\t1\tbacteria\tp\tunused.gff\t"+root.resolve("a.fa")+"\t4\t1000\t1",
			"1\tg2\ta2\t2\tbacteria\tp\tunused.gff\t"+root.resolve("b.fa")+"\t0\t0\t0");
		write(root.resolve("legacy.tsv"), "locus_id\tgenome\tseqid\tstrand\tsource_start1\tsource_stop1\tcm_start1\tcm_stop1\tcm_state",
			"one\tg1\tc1\t+\t100\t199\t100\t199\tCM_FULL_SPAN_SUPPORTED",
			"unknown\tg1\tc1\t+\t600\t699\t0\t0\tNO_REPORTED_HIT");
		write(root.resolve("rf14.tsv"), NcrnaClanOverlay.HEADER,
			overlay("fourteen", "g1", "a1", "c1", "RF01482", 300, 399, 1000, "RF01482"),
			overlay("clan", "g1", "a1", "c1", "RF01482", 400, 499, 1000, "B12clade"),
			overlay("excluded", "g3", "a3", "outside", "RF01482", 10, 109, 1000, "RF01482"));
		write(root.resolve("rf16.tsv"), NcrnaClanOverlay.HEADER,
			overlay("sixteen", "g2", "a2", "c2", "RF01689", 50, 119, 300, "RF01689"),
			overlay("unknown16", "g2", "a2", "c2", "RF01689", 200, 269, 300, "UNRESOLVED"));
		final String[] prep={"roster="+root.resolve("roster.tsv"), "legacy="+root.resolve("legacy.tsv"), "rf14="+root.resolve("rf14.tsv"),
			"rf16="+root.resolve("rf16.tsv"), "outprefix="+root.resolve("panel"), "genomes=2", "contigs=2", "legacyrows=2", "legacyfull=1", "rf14rows=3", "rf16rows=2"};
		CobalaminPanelPrep.main(prep);
		check(lines(root.resolve("panel.loci.tsv"))==7, "Two panel genomes retain six loci; the third genome is excluded, not relabeled");
		write(root.resolve("calls.gff"), "##gff-version 3",
			call("c1", 100, 199, "RF01482"), call("c1", 300, 399, "RF01482"),
			call("c1", 400, 499, "RF00174"), call("c1", 600, 699, "RF01689"),
			call("c2", 50, 119, "RF01689"), call("c2", 50, 119, "RF01689"),
			call("c2", 280, 299, "RF00174"), call("c2", 200, 269, "RF01689"),
			"c1\tCallGenes\tRNA\t800\t899\t.\t+\t.\tmodel:RNaseP_x,family:rnasep");
		CobalaminPanelGrade.main(new String[]{"panel="+root.resolve("panel.loci.tsv"), "contigs="+root.resolve("panel.contigs.tsv"),
			"in="+root.resolve("calls.gff"), "outprefix="+root.resolve("grade")});
		checkSummary(root.resolve("grade.summary.tsv"), "SOURCE", "union", 6, 8, 6, 2, 0, 0, 1);
		checkSummary(root.resolve("grade.summary.tsv"), "FULL_CM", "union", 4, 8, 4, 2, 0, 2, 1);
		checkSummary(root.resolve("grade.summary.tsv"), "FULL_CM", "RF00174", 1, 2, 0, 1, 1, 1, 0);
		checkSummary(root.resolve("grade.summary.tsv"), "FULL_CM", "RF01482", 1, 2, 1, 1, 0, 0, 0);
		checkSummary(root.resolve("grade.summary.tsv"), "FULL_CM", "RF01689", 1, 4, 1, 1, 0, 2, 1);
		check(lines(root.resolve("grade.genomes.tsv"))==17, "Every view/scope reports every genome, including zero-call genomes");
		write(root.resolve("bad.gff"), call("c1", 100, 199, "RF00174").replace("family:cobalamin_rf00174", "family:cobalamin_rf01482"));
		expectFailure(()->CobalaminPanelGrade.main(new String[]{"panel="+root.resolve("panel.loci.tsv"), "contigs="+root.resolve("panel.contigs.tsv"), "in="+root.resolve("bad.gff"), "outprefix="+root.resolve("bad_grade")}));
		write(root.resolve("bad_overlay.tsv"), NcrnaClanOverlay.HEADER, overlay("bad", "g1", "wrong_assembly", "c1", "RF01482", 300, 399, 1000, "RF01482"));
		final String[] badPrep=prep.clone();badPrep[2]="rf14="+root.resolve("bad_overlay.tsv");badPrep[4]="outprefix="+root.resolve("bad_panel");badPrep[9]="rf14rows=1";
		expectFailure(()->CobalaminPanelPrep.main(badPrep));
		write(root.resolve("overlap.panel.tsv"), CobalaminPanelPrep.HEADER,
			"a\tg1\tc1\tRF00174\tRF00174\t100\t199\t100\t199\t+",
			"b\tg1\tc1\tRF00174\tRF00174\t150\t249\t150\t249\t+");
		write(root.resolve("overlap.gff"), call("c1", 140, 209, "RF00174"), call("c1", 140, 209, "RF00174"));
		CobalaminPanelGrade.main(new String[]{"panel="+root.resolve("overlap.panel.tsv"), "contigs="+root.resolve("panel.contigs.tsv"),
			"in="+root.resolve("overlap.gff"), "outprefix="+root.resolve("overlap_grade")});
		checkSummary(root.resolve("overlap_grade.summary.tsv"), "FULL_CM", "union", 2, 2, 1, 1, 1, 0, 1);
		checkSummary(root.resolve("overlap_grade.summary.tsv"), "SOURCE", "RF00174", 2, 2, 1, 1, 1, 0, 1);
		write(root.resolve("collision.panel.tsv"), CobalaminPanelPrep.HEADER,
			"a\tg1\tc1\tRF00174\tRF00174\t100\t199\t100\t199\t+",
			"b\tg1\tc1\tRF00174\tRF00174\t110\t209\t100\t199\t+");
		expectFailure(()->CobalaminPanelGrade.main(new String[]{"panel="+root.resolve("collision.panel.tsv"), "contigs="+root.resolve("panel.contigs.tsv"),
			"in="+root.resolve("overlap.gff"), "outprefix="+root.resolve("collision_grade")}));
		final String[] shiftedFields=overlay("collision", "g1", "a1", "c1", "RF01482", 310, 409, 1000, "RF01482").split("\t", -1);
		shiftedFields[18]="20";shiftedFields[19]="20";shiftedFields[27]="300";shiftedFields[28]="399";shiftedFields[32]="GEOMETRY-CORRECTED";
		final String shifted=String.join("\t", shiftedFields);
		write(root.resolve("collision.overlay.tsv"), NcrnaClanOverlay.HEADER,
			overlay("original", "g1", "a1", "c1", "RF01482", 300, 399, 1000, "RF01482"), shifted);
		final String[] collisionPrep=prep.clone();collisionPrep[2]="rf14="+root.resolve("collision.overlay.tsv");
		collisionPrep[4]="outprefix="+root.resolve("collision_panel");collisionPrep[9]="rf14rows=2";
		expectFailure(()->CobalaminPanelPrep.main(collisionPrep));
		System.out.println("COBALAMIN_PANEL_TEST_PASS");
	}
	private static String overlay(String id, String genome, String assembly, String seqid, String family, int start, int stop, int contigLength, String member){
		final boolean admitted=!member.equals("UNRESOLVED"), clan=member.equals("B12clade");
		return String.join("\t", id, genome, assembly, "1", "bacteria", "p", seqid, "+", start+"-"+stop, ""+(stop-start+1), "1", id,
			"false", "false", "ANNOTATION_FULL_COPY", ""+contigLength, "false", "0", "0", "0", "1", family,
			admitted?clan?"CROSS_LABELED":family:"UNRESOLVED", ""+admitted, admitted?family:"NONE", admitted?family:"NA", "RAW_ADJUDICATED",
			admitted?""+start:"NA", admitted?""+stop:"NA", admitted?""+(stop-start+1):"NA", admitted?"WITHIN_BAND":"NA", member, admitted?"SOURCE_EXACT":"UNRESOLVED");
	}
	private static String call(String seqid, int start, int stop, String family){
		return seqid+"\tCallGenes\tRNA\t"+start+'\t'+stop+"\t.\t+\t.\tmodel:cobalamin_"+family+"_test,family:cobalamin_"+family.toLowerCase(java.util.Locale.ROOT);
	}
	private static void fasta(Path path, String name, int length){
		final ByteBuilder b=new ByteBuilder();b.append('>').append(name).nl();for(int i=0; i<length; i++){b.append('A');}write(path, b.toString());
	}
	private static void write(Path path, String... lines){
		final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());try{for(String s:lines){out.println(s);}}finally{RrnaResourceIO.finish(out);}
	}
	private static int lines(Path path){
		final ByteFile in=RrnaResourceIO.open(path.toString());int n=0;try{while(in.nextLine()!=null){n++;}}finally{check(!in.close(), "Fixture read failed");}return n;
	}
	private static void checkSummary(Path path, String view, String scope, int... counts){
		final ByteFile in=RrnaResourceIO.open(path.toString());final LineParser1 p=new LineParser1('\t');boolean found=false;
		try{for(byte[] row; (row=in.nextLine())!=null; ){
			p.set(row);if(!p.termEquals(view, 0) || !p.termEquals(scope, 1)){continue;}
			check(!found, "Repeated report row");found=true;for(int i=0; i<counts.length; i++){check(p.parseInt(3+i)==counts[i], view+"/"+scope+" count column"+i+" violates the independent fixture expectation");}
		}}finally{check(!in.close(), "Summary read failed");}check(found, "Expected grade scope absent");
	}
	private static void expectFailure(Runnable action){boolean failed=false;try{action.run();}catch(IllegalArgumentException expected){failed=true;}check(failed, "Invalid attribution/assembly must fail before producing a grade");}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
