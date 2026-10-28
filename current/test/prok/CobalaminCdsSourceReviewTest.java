package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import prok.CobalaminCdsSourceReview;
import prok.CobalaminPanelPrep;
import prok.GffCallSetDelta;
import prok.NcrnaCalibrationPrep;
import prok.RrnaResourceIO;

/** Source presence and start accuracy use the strand-aware translation stop;
 * multiple source matches stay explicit and shortened genes are not whole losses.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminCdsSourceReviewTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh fixture directory");final Path root=Paths.get(args[0]);Files.createDirectory(root);
		write(root, "contigs.tsv", CobalaminPanelPrep.CONTIG_HEADER+"\ng\tassembly\tseq\t4000\n");
		final String source=feature(160, 399, '+', "plus_best")+feature(130, 399, '+', "plus_alternative")
			+feature(600, 899, '-', "minus")+feature(1200, 1499, '+', "removed")+feature(1800, 2099, '-', "wrong_strand");
		write(root, "source.gff", source);
		write(root, "roster.tsv", NcrnaCalibrationPrep.ROSTER_HEADER+"\n0\tg\tassembly\t1\tbacteria\tfixture\t"+root.resolve("source.gff")+"\tgenome.fa\t0\t4000\t1\n");
		// Synthetic hashes exercise the schema only; the execution producer hashes real inputs.
		write(root, "hashes.tsv", "00000000000000000000\t"+root.resolve("source.gff")+'\n');
		final String changes=GffCallSetDelta.CHANGE_HEADER+'\n'+delta("REMOVED", 100, 399, '+')+delta("ADDED", 160, 399, '+')
			+delta("REMOVED", 600, 899, '-')+delta("ADDED", 600, 839, '-')+delta("REMOVED", 1200, 1499, '+')
			+delta("REMOVED", 1800, 2099, '+')+delta("REMOVED", 2400, 2699, '-')+delta("ADDED", 3000, 3299, '+');
		write(root, "changes.tsv", changes);run(root, "changes.tsv", "result", 5);
		check(contains(root.resolve("result.summary.tsv"), "2\tg\tseq\t+\t100\t399\t160\t399\tSAME_STOP_START_CHANGED\t2\tANNOTATED-START-CHANGE\t1\t0\t1"),
			"Both alternative source starts remain visible; one improves and one has equal absolute error");
		check(contains(root.resolve("result.summary.tsv"), "4\tg\tseq\t-\t600\t899\t600\t839\tSAME_STOP_START_CHANGED\t1\tANNOTATED-START-CHANGE\t0\t1\t0"),
			"Minus-strand translation stop is genomic start; shortening the exact source start worsens its error");
		check(contains(root.resolve("result.summary.tsv"), "6\tg\tseq\t+\t1200\t1499\tNA\tNA\tREMOVED\t1\tREAL-GENE-LOSS"),
			"An absent prediction with a matching annotated translation stop gets the requested source-supported loss label");
		check(contains(root.resolve("result.summary.tsv"), "7\tg\tseq\t+\t1800\t2099\tNA\tNA\tREMOVED\t0\tPREDICTED-ONLY")
			&& contains(root.resolve("result.summary.tsv"), "8\tg\tseq\t-\t2400\t2699\tNA\tNA\tREMOVED\t0\tPREDICTED-ONLY"),
			"Opposite-strand annotation and absent annotation cannot become same-stop source support");
		check(contains(root.resolve("result.added.tsv"), "9\tg\tseq\t3000\t3299\t+\tfalse"), "An unpaired new CDS row must not disappear");
		write(root, "ambiguous.tsv", changes+delta("REMOVED", 103, 399, '+'));
		boolean failed=false;try{run(root, "ambiguous.tsv", "ambiguous", 6);}catch(IllegalArgumentException expected){failed=true;}
		check(failed && !Files.exists(root.resolve("ambiguous.summary.tsv")), "Ambiguous removed CDS anchors must fail before claiming biological gene counts");
		System.out.println("COBALAMIN_CDS_SOURCE_REVIEW_TEST_PASS");
	}
	private static String feature(int start, int stop, char strand, String id){return "seq\tNCBI\tCDS\t"+start+'\t'+stop+"\t.\t"+strand+"\t0\tID="+id+";product=fixture protein;\n";}
	private static String delta(String direction, int start, int stop, char strand){return direction+"\tseq\tBBTools\tCDS\t"+start+'\t'+stop+"\t10\t"+strand+"\t0\tCDS\n";}
	private static void run(Path root, String changes, String out, int removed){
		CobalaminCdsSourceReview.main(new String[]{"changes="+root.resolve(changes), "contigs="+root.resolve("contigs.tsv"), "roster="+root.resolve("roster.tsv"),
			"gffhashes="+root.resolve("hashes.tsv"), "outprefix="+root.resolve(out), "removed="+removed, "added=3"});
	}
	private static boolean contains(Path file, String prefix){
		final ByteFile in=RrnaResourceIO.open(file.toString());boolean found=false;
		try{for(byte[] row; (row=in.nextLine())!=null; ){if(new String(row, java.nio.charset.StandardCharsets.UTF_8).startsWith(prefix)){found=true;}}}
		finally{check(!in.close(), "Fixture detail read failed");}return found;
	}
	private static void write(Path root, String name, String text){final ByteStreamWriter out=RrnaResourceIO.writer(root.resolve(name).toString());try{out.print(text);}finally{RrnaResourceIO.finish(out);}}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
