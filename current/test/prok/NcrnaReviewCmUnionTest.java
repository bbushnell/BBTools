package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import prok.NcrnaReviewCmUnion;
import prok.NcrnaReviewJoin;
import prok.NcrnaReviewUnion;
import prok.RrnaResourceIO;

/** Overlapping model support counts one context; absent evidence never means negative.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaReviewCmUnionTest {

	public static void main(String[] args)throws Exception{
		require(args.length==1, "Supply a fresh fixture directory");final Path dir=Paths.get(args[0]);Files.createDirectory(dir);
		write(dir.resolve("loci.tsv"), NcrnaReviewUnion.LOCUS_HEADER, locus("a"), locus("b"), locus("c"));
		write(dir.resolve("contexts.fa"), ">a lflank=2 rflank=1", "AACGTAA", ">b lflank=2 rflank=1", "GGCCTTA", ">c lflank=2 rflank=1", "TTAACCG");
		write(dir.resolve("RF00174.tsv"), NcrnaReviewJoin.CM_HEADER, evidence("a", true), evidence("b", false), evidence("c", false));
		write(dir.resolve("RF01482.tsv"), NcrnaReviewJoin.CM_HEADER, evidence("a", false), evidence("b", false), evidence("c", false));
		write(dir.resolve("RF01689.tsv"), NcrnaReviewJoin.CM_HEADER, evidence("a", true), evidence("b", false), evidence("c", true));
		manifest(dir, "models.tsv", "RF01482.tsv");run(dir, "models.tsv", "union");
		require(lines(dir.resolve("union.summary.tsv")).get(1).equals("3\t2\t1\t7"), "Two models supporting a count as one context; only b remains unresolved");
		final ArrayList<String> fasta=new ArrayList<String>();
		RrnaResourceIO.stream(dir.resolve("union.unresolved.contexts.fa.gz").toString(), r->{fasta.add(r.id);fasta.add(new String(r.bases, java.nio.charset.StandardCharsets.US_ASCII));});
		require(fasta.size()==2 && fasta.get(0).equals("b lflank=2 rflank=1") && fasta.get(1).equals("GGCCTTA"), "Sensitive subset preserves the original ID, flanks, orientation and bases");
		require(lines(dir.resolve("union.unresolved.loci.tsv")).get(1).equals(locus("b")), "Sensitive subset metadata must be an unchanged original row");
		write(dir.resolve("missing.tsv"), NcrnaReviewJoin.CM_HEADER, evidence("a", false), evidence("c", false));
		manifest(dir, "bad_models.tsv", "missing.tsv");boolean failed=false;
		try{run(dir, "bad_models.tsv", "bad");}catch(IllegalArgumentException e){failed=e.getMessage().contains("Missing model evidence");}
		require(failed, "A missing model row cannot silently mark b unsupported");
		System.out.println("NCRNA_REVIEW_CM_UNION_TEST_PASS");
	}
	private static String locus(String id){return id+"\tg\ta\tgenome.fa\tc\tcobalamin\t-\t10\t13\t100\t2\t1\t7\tPREDICTION_UNRESOLVED";}
	private static String evidence(String id, boolean full){return id+"\t1\t1\t"+(full?"1\tONE_FULL_GA_OVERLAPPING_HIT":"0\tNO_FULL_GA_OVERLAPPING_HIT");}
	private static void manifest(Path dir, String name, String middle){
		write(dir.resolve(name), "model_id\tstage\tloci", "RF00174\tDEFAULT\t"+dir.resolve("RF00174.tsv"),
			"RF01482\tDEFAULT\t"+dir.resolve(middle), "RF01689\tDEFAULT\t"+dir.resolve("RF01689.tsv"));
	}
	private static void run(Path dir, String manifest, String name){
		NcrnaReviewCmUnion.main(new String[]{"loci="+dir.resolve("loci.tsv"), "in="+dir.resolve("contexts.fa"),
			"models="+dir.resolve(manifest), "expected=3", "outprefix="+dir.resolve(name)});
	}
	private static ArrayList<String> lines(Path path){
		final ByteFile in=RrnaResourceIO.open(path.toString());final ArrayList<String> result=new ArrayList<String>();
		try{for(byte[] row; (row=in.nextLine())!=null; ){result.add(new String(row, java.nio.charset.StandardCharsets.US_ASCII));}}
		finally{require(!in.close(), "Fixture output read failed");}return result;
	}
	private static void write(Path path, String... lines){
		final ByteStreamWriter out=RrnaResourceIO.writer(path.toString());try{for(String s:lines){out.println(s);}}finally{RrnaResourceIO.finish(out);}
	}
	private static void require(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
