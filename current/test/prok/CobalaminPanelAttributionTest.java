package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import prok.CobalaminPanelAttribution;
import prok.CobalaminPanelPrep;
import prok.NcrnaReviewUnion;
import prok.RrnaResourceIO;

/** Separates assigned label swaps, physical overlap, shifted boundaries and strand changes.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminPanelAttributionTest {

	public static void main(String[] args)throws Exception{
		require(args.length==1, "Supply a fresh fixture directory");final Path dir=Paths.get(args[0]);Files.createDirectory(dir);
		write(dir.resolve("contigs.tsv"), CobalaminPanelPrep.CONTIG_HEADER, "g\ta\tc\t1000");
		write(dir.resolve("truth.tsv"), CobalaminPanelPrep.HEADER,
			"r1\tg\tc\tRF00174\tRF00174\t10\t109\t10\t109\t+",
			"r16\tg\tc\tRF01689\tRF01689\t500\t579\t500\t579\t+",
			"unknown\tg\tc\tRF01482\tUNRESOLVED\t800\t879\t0\t0\t+");
		for(int k=7; k<=9; k++){
			write(dir.resolve(k+".tsv"), NcrnaReviewUnion.PANEL_CALL_HEADER,
				call("0", 10, 109, "+", k==8?"RF00174":"RF01482", "TP", "r1"),
				call("1", k==8?501:500, k==8?580:579, "+", "RF01689", "TP", "r16"),
				call("2", 300, 399, k==8?"-":"+", "RF01482", "APPARENT_FP", "NA"),
				call("3", 800, 879, "+", "RF01689", "UNSCORED", "NA"));
		}
		write(dir.resolve("runs.tsv"), CobalaminPanelAttribution.RUN_HEADER,
			"k7\t7\t100\t0.64\t"+dir.resolve("7.tsv"), "k8\t8\t100\t0.64\t"+dir.resolve("8.tsv"), "k9\t9\t100\t0.64\t"+dir.resolve("9.tsv"));
		run(dir, "audit", "runs.tsv");
		require(contains(dir.resolve("audit.summary.tsv"), "k7\t7\t100\t0.6400\t4\t2\t1\t1\t4\t1\t2\t2\t1\t1"),
			"A sibling assigned to RF00174 is one direct label swap; unresolved evidence is not a full core");
		require(contains(dir.resolve("audit.overlap.tsv"), "k7\tk8\t4\t4\t2\t1\t3\t1\t1"),
			"An exact label change, shifted same-strand match and opposite-strand nonmatch must stay distinct");
		require(contains(dir.resolve("audit.overlap.tsv"), "k7\tk9\t4\t4\t4\t0\t4\t0\t0"), "Identical predictions must have exact and reciprocal set parity");
		write(dir.resolve("bad.tsv"), NcrnaReviewUnion.PANEL_CALL_HEADER,
			call("0", 300, 399, "+", "RF01482", "TP", "r1"));
		write(dir.resolve("bad_runs.tsv"), CobalaminPanelAttribution.RUN_HEADER,
			"k7\t7\t100\t0.64\t"+dir.resolve("bad.tsv"), "k8\t8\t100\t0.64\t"+dir.resolve("8.tsv"), "k9\t9\t100\t0.64\t"+dir.resolve("9.tsv"));
		boolean failed=false;
		try{run(dir, "bad", "bad_runs.tsv");}catch(IllegalArgumentException e){failed=e.getMessage().contains("Saved union assignment");}
		require(failed, "A saved TP pointing at a nonoverlapping reference must fail rather than manufacture a label swap");
		System.out.println("COBALAMIN_PANEL_ATTRIBUTION_TEST_PASS");
	}
	private static void run(Path dir, String name, String manifest){
		CobalaminPanelAttribution.main(new String[]{"panel="+dir.resolve("truth.tsv"), "contigs="+dir.resolve("contigs.tsv"),
			"runs="+dir.resolve(manifest), "expectedruns=3", "outprefix="+dir.resolve(name)});
	}
	private static String call(String id, int start, int stop, String strand, String family, String outcome, String truth){
		return "FULL_CM\tunion\t"+id+"\tg\tc\t"+family+'\t'+start+'\t'+stop+'\t'+strand+"\tcobalamin_"+family+"_test\t"+outcome+'\t'+truth+"\tNA\tNA";
	}
	private static boolean contains(Path path, String expected){
		final ByteFile in=RrnaResourceIO.open(path.toString());boolean found=false;
		try{for(byte[] row; (row=in.nextLine())!=null; ){if(expected.equals(new String(row, java.nio.charset.StandardCharsets.US_ASCII))){found=true;}}}
		finally{require(!in.close(), "Fixture output read failed");}return found;
	}
	private static void write(Path path, String... lines){
		final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());try{for(String s:lines){out.println(s);}}finally{RrnaResourceIO.finish(out);}
	}
	private static void require(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
