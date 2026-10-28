package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.NcrnaReviewJoin;
import prok.RrnaResourceIO;

/** Checks subtype/geometry ambiguity and exact evidence joins without truth edits.
 * @author Brian, Raiden
 */
public final class NcrnaReviewJoinTest {
	public static void main(String[] args){
		require(args.length==1, "Supply a fresh fixture directory");final File dir=new File(args[0]);
		require(!dir.exists() && dir.mkdirs(), "Fixture must not overwrite existing evidence");final String base=dir.getAbsolutePath();
		write(base+"/context.tsv", NcrnaReviewJoin.CONTEXT_HEADER, context("r0", 0, 0, 0, 0), context("r1", 0, 0, 1, 0),
			context("r2", 1, 0, 0, 0), context("r3", 0, 1, 0, 0), context("r4", 0, 0, 1, 1), context("r5", 0, 0, 0, 0));
		write(base+"/overlaps.tsv", NcrnaReviewJoin.OVERLAP_HEADER,
			feature("r1", 1, "OTHER_FEATURE", false, "+", "RF01482", 20), feature("r2", 2, "TARGET_RNA_ANNOTATION", true, "+", "RF00174", 8),
			feature("r3", 3, "CDS", false, "+", "none", 20), feature("r4", 4, "OTHER_FEATURE", false, "-", "RF01689", 20));
		write(base+"/hits.tsv", NcrnaReviewJoin.HIT_HEADER, hit("r0", 1, 101, 120, 1, 20, 40, 20, true),
			hit("r1", 2, 101, 120, 5, 20, 35, 20, false), hit("r2", 3, 103, 110, 1, 20, 33, 8, false));
		write(base+"/cm.tsv", NcrnaReviewJoin.CM_HEADER, cm("r0", 1, 1), cm("r1", 1, 0), cm("r2", 1, 0), cm("r3", 0, 0), cm("r4", 0, 0), cm("r5", 0, 0));
		final String[] members={NcrnaReviewJoin.MEMBER_HEADER, member("a", "c0", "r0"), member("a", "c1", "r1"), member("a", "c2", "r2"),
			member("a", "c3", "r3"), member("a", "c4", "r4"), member("a", "c5", "r5"), member("b", "c0", "r0")};
		write(base+"/members.tsv", members);run(base, "members.tsv", "hits.tsv", "overlaps.tsv", "good");
		final ArrayList<String> rows=lines(base+"/good.calls.tsv"), summary=lines(base+"/good.summary.tsv");
		require(rows.size()==8 && summary.size()==7, "Seven memberships remain distinct; only identical run/evidence classes aggregate");
		require(rows.get(1).endsWith("\tFULL_CM_INTERGENIC_SUPPORT\tEVIDENCE_ONLY_NO_TRUTH_EDIT") && rows.get(7).endsWith("\tFULL_CM_INTERGENIC_SUPPORT\tEVIDENCE_ONLY_NO_TRUTH_EDIT"), "Two runs retain the same exact locus without merging their calls");
		require(rows.get(2).contains("\tOTHER_B12_ANNOTATION_SCOPE_UNRESOLVED\t"), "Sibling annotation is not silently assigned to RF00174");
		final LineParser1 p=new LineParser1('\t');p.set(rows.get(3).getBytes(StandardCharsets.UTF_8));
		require(rows.get(3).contains("\tTARGET_ANNOTATION_GEOMETRY_REVIEW\t") && p.parseDouble(43)==33 && p.parseDouble(44)==1.5
			&& p.parseInt(45)==8 && p.termEquals("true", 46) && p.termEquals("false", 47), "A full CM within an oversized prediction preserves the failed reciprocal geometry");
		require(rows.get(4).contains("\tNO_REPORTED_CM_FEATURE_OVERLAP_AMBIGUOUS\t") && rows.get(5).contains("\tNO_REPORTED_CM_FEATURE_OVERLAP_AMBIGUOUS\t"), "CDS and opposite-strand B12 evidence do not establish a TP or FP");
		require(rows.get(6).contains("\tNO_REPORTED_CM_INTERGENIC_AMBIGUOUS\t"), "CM absence in an intergenic call is still unresolved");
		final String[] duplicate=members.clone();duplicate[7]=member("a", "c0", "r0");write(base+"/duplicate.tsv", duplicate);
		fails(base, "duplicate.tsv", "hits.tsv", "overlaps.tsv", "dup", "repeated membership");
		write(base+"/bad_hit.tsv", NcrnaReviewJoin.HIT_HEADER, hit("r0", 1, 101, 120, 1, 20, 40, 20, false));
		fails(base, "members.tsv", "bad_hit.tsv", "overlaps.tsv", "bad_hit", "Stored CM qualification");
		write(base+"/bad_source.tsv", NcrnaReviewJoin.OVERLAP_HEADER, feature("r1", 1, "OTHER_FEATURE", false, "+", "RF01482", 20).replace("\tg\ta\t", "\twrong\ta\t"));
		fails(base, "members.tsv", "hits.tsv", "bad_source.tsv", "bad_source", "Feature source differs");
		System.out.println("NCRNA_REVIEW_JOIN_TEST_PASS");
	}
	private static void run(String b, String members, String hits, String overlaps, String out){
		assert(!out.isEmpty()):"Each invocation writes to its own output prefix";
		NcrnaReviewJoin.main(new String[]{"context="+b+"/context.tsv", "overlaps="+b+'/'+overlaps, "hits="+b+'/'+hits, "cm="+b+"/cm.tsv",
			"members="+b+'/'+members, "outprefix="+b+'/'+out, "expected=6", "events=7", "runs=2", "clen=20", "gabits=31.5"});
	}
	private static void fails(String b, String members, String hits, String overlaps, String out, String message){
		boolean failed=false;try{run(b, members, hits, overlaps, out);}catch(IllegalArgumentException e){failed=e.getMessage().contains(message);}
		require(failed && !new File(b+'/'+out+".calls.tsv").exists(), "Bad input must fail before publishing any joined output: "+message);
	}
	private static String context(String id, int target, int cds, int other, int opposite){
		final int n=target+cds+other;assert(n<=1 && opposite<=n):"Fixture has at most one feature per locus";
		return id+"\tg\ta\tchr\t101\t120\t+\t20\t20\t0\t1.000000\t"+n+"\t0\t0\t"+target+"\t0\t"+cds+'\t'+other+'\t'+(n-opposite)+'\t'+opposite
			+"\t0\t"+(n==0?"INTERGENIC":"FEATURE_OVERLAP")+"\tfixture.gff\t0123456789abcdef0123\tfixture.fa\t0123456789abcdef0123";
	}
	private static String feature(String id, int line, String kind, boolean target, String strand, String model, int overlap){
		assert(line>0 && overlap>0):"Fixture feature has a real line and positive overlap";
		return id+"\tg\ta\tfixture.fa\t0123456789abcdef0123\tfixture.gff\t0123456789abcdef0123\t"+line+'\t'+kind+'\t'+target+'\t'+overlap+'\t'
			+(strand.equals("+")?"SAME":"OPPOSITE")+"\tfalse\tchr\tRefSeq\t"+(kind.equals("CDS")?"CDS":"riboswitch")+"\t101\t"+(100+overlap)+"\t.\t"+strand+"\t.\tID=f"+line+";Dbxref=RFAM:"+model;
	}
	private static String hit(String id, int line, int start, int stop, int from, int to, int bits, int overlap, boolean strict){
		assert(stop>=start):"Fixture CM coordinates preserve inclusive genomic ordering";
		return id+'\t'+line+"\tchr\t"+start+'\t'+stop+"\t+\t"+from+'\t'+to+"\tno\t"+bits+"\t1e-10\ttrue\t"+overlap+'\t'+strict;
	}
	private static String cm(String id, int hits, int strict){return id+'\t'+hits+'\t'+hits+'\t'+strict+'\t'+(strict==0?"NO_FULL_GA_OVERLAPPING_HIT":"ONE_FULL_GA_OVERLAPPING_HIT");}
	private static String member(String run, String call, String id){return run+'\t'+call+'\t'+id+"\tAPPARENT_FP\tAPPARENT_FP";}
	private static void write(String path, String... rows){
		assert(rows.length>0):"Every fixture input has its declared schema";
		final ByteStreamWriter w=RrnaResourceIO.writer(path);try{for(String row:rows){w.println(row);}}finally{RrnaResourceIO.finish(w);}
	}
	private static ArrayList<String> lines(String path){
		final ByteFile in=RrnaResourceIO.open(path);final ArrayList<String> rows=new ArrayList<String>();
		try{for(byte[] row; (row=in.nextLine())!=null; ){rows.add(new String(row, StandardCharsets.UTF_8));}}
		finally{require(!in.close(), "Fixture output read failed");}return rows;
	}
	private static void require(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
