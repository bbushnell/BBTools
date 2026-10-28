package test.prok;

import java.io.File;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import prok.NcrnaCalibrationPrep;
import prok.NcrnaReviewUnion;
import prok.RrnaResourceIO;

/** Exercises exact union, run-local IDs, clipped strand geometry and source ambiguity.
 * @author Brian, Raiden
 */
public final class NcrnaReviewUnionTest {
	public static void main(String[] args){
		require(args.length==1, "Supply a fresh fixture directory");
		final File root=new File(args[0]);require(!root.exists() && root.mkdirs(), "Fixture must not replace existing evidence");
		final String base=root.getAbsolutePath(), fasta=base+"/genome.fa", roster=base+"/roster.tsv", runs=base+"/runs.tsv";
		write(fasta, ">chr", "AAAAAAAAAAAAAAAAAAAACCCCCCCCCCGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTT");
		write(roster, NcrnaCalibrationPrep.ROSTER_HEADER, "0\tg0\ta0\t1\tbacteria\tp\tNA\t"+fasta+"\t0\t70\t1");
		write(base+"/a.tsv", NcrnaReviewUnion.CALL_HEADER,
			row("c0", 2, 8, "+", true), row("c1", 65, 69, "-", true), row("c2", 5, 12, "+", true),
			row("c3", 2, 8, "-", true), row("c4", 30, 40, "+", false));
		write(base+"/b.tsv", NcrnaReviewUnion.CALL_HEADER, row("c0", 2, 8, "+", true));
		write(runs, "run\tcalls\treview_candidates", "a\t"+base+"/a.tsv\t4", "b\t"+base+"/b.tsv\t1");
		run(runs, roster, base+"/union", 1);
		final ArrayList<String> summary=lines(base+"/union.summary.tsv"), fa=new ArrayList<String>();
		RrnaResourceIO.stream(base+"/union.contexts.fa.gz", r->{fa.add('>'+r.id);fa.add(new String(r.bases, java.nio.charset.StandardCharsets.UTF_8));});
		require(summary.size()==2 && summary.get(1).equals("2\t1\t1\t5\t4\t45\t3"), "One exact duplicate merges; overlaps and opposite strands remain, with all five memberships");
		require(fa.size()==8 && fa.get(0).equals(">review_0 lflank=1 rflank=3") && fa.get(1).equals("AAAAAAAAAAA"), "Plus context must clip at the actual left contig edge");
		require(fa.get(2).equals(">review_1 lflank=3 rflank=3") && fa.get(3).equals("AAAAAAAAAAAAAA"), "Overlapping nonidentical call must keep its own complete context");
		require(fa.get(4).equals(">review_2 lflank=3 rflank=1") && fa.get(5).equals("TTTTTTTTTTT"), "Minus-strand context must reverse-complement and swap clipped flanks");
		require(fa.get(6).equals(">review_3 lflank=1 rflank=3") && fa.get(7).equals("AAAAAAAAA"), "Minus context at the right contig edge must clip the sense leader");
		final ArrayList<String> members=lines(base+"/union.members.tsv");
		require(members.size()==6 && members.get(1).startsWith("a\tc0\treview_0\t") && members.get(5).startsWith("b\tc0\treview_0\t"), "Run-local candidate IDs must remain distinct members of the shared physical locus");
		write(base+"/ambiguous.tsv", NcrnaCalibrationPrep.ROSTER_HEADER,
			"0\tg0\ta0\t1\tbacteria\tp\tNA\t"+fasta+"\t0\t70\t1", "1\tg1\ta1\t2\tbacteria\tp\tNA\t"+fasta+"\t0\t70\t1");
		boolean rejected=false;
		try{run(runs, base+"/ambiguous.tsv", base+"/bad", 2);}catch(IllegalArgumentException e){rejected=e.getMessage().contains("ambiguous");}
		require(rejected && !new File(base+"/bad.contexts.fa.gz").exists() && !new File(base+"/bad.summary.tsv").exists(), "Ambiguous source identity must fail before any completed extraction is published");
		write(base+"/panel.tsv", NcrnaReviewUnion.PANEL_CALL_HEADER,
			panel("SOURCE", "union", "0", 2, 8, "APPARENT_FP"), panel("SOURCE", "union", "1", 30, 40, "TP"),
			panel("SOURCE", "union", "2", 50, 60, "APPARENT_FP"),
			panel("FULL_CM", "RF01482", "0", 2, 8, "APPARENT_FP"),
			panel("FULL_CM", "union", "0", 2, 8, "TP"), panel("FULL_CM", "union", "1", 30, 40, "APPARENT_FP"),
			panel("FULL_CM", "union", "2", 50, 60, "UNSCORED"));
		write(base+"/panel_runs.tsv", "run\tcalls\treview_candidates", "panel\t"+base+"/panel.tsv\t1");
		runPanel(base+"/panel_runs.tsv", roster, base+"/panel_union");
		final ArrayList<String> panelMembers=lines(base+"/panel_union.members.tsv");
		require(panelMembers.size()==2 && panelMembers.get(1).equals("panel\t1\treview_0\tTP\tAPPARENT_FP"),
			"Review only FULL_CM union FPs: exclude subtype label-swap FPs and unresolved calls, retain SOURCE outcome");
		write(base+"/mismatch.tsv", NcrnaReviewUnion.PANEL_CALL_HEADER,
			panel("SOURCE", "union", "1", 30, 40, "TP"), panel("FULL_CM", "union", "1", 31, 40, "APPARENT_FP"));
		write(base+"/mismatch_runs.tsv", "run\tcalls\treview_candidates", "mismatch\t"+base+"/mismatch.tsv\t1");
		rejected=false;
		try{runPanel(base+"/mismatch_runs.tsv", roster, base+"/mismatch");}catch(IllegalArgumentException e){rejected=e.getMessage().contains("same call identities");}
		require(rejected, "A call ID cannot silently join different coordinates across reference views");
		System.out.println("NCRNA_REVIEW_UNION_TEST_PASS");
	}
	private static String panel(String view, String scope, String id, int start, int stop, String outcome){
		assert(start>0 && stop>=start):"Native panel fixtures retain inclusive physical coordinates";
		return view+'\t'+scope+'\t'+id+"\tg0\tchr\tRF01482\t"+start+'\t'+stop+"\t+\tcobalamin_RF01482_test\t"+outcome+"\tNA\tNA\tNA";
	}
	private static void runPanel(String runs, String roster, String out){
		NcrnaReviewUnion.main(new String[]{"runs="+runs, "roster="+roster, "outprefix="+out, "family=cobalamin", "panel=t", "genomes=1", "contigs=1", "runcount=1", "flank=3"});
	}
	private static void run(String runs, String roster, String out, int genomes){
		assert(genomes>0):"Fixture declares its complete original source roster";
		NcrnaReviewUnion.main(new String[]{"runs="+runs, "roster="+roster, "outprefix="+out, "family=cobalamin", "genomes="+genomes, "contigs=1", "runcount=2", "flank=3"});
	}
	private static String row(String id, int start, int stop, String strand, boolean review){
		assert(start>0 && stop>=start):"Fixture coordinates are one-based inclusive";
		final String result=review?"APPARENT_FP":"TP";
		return id+"\tchr\tcobalamin\t"+start+'\t'+stop+'\t'+strand+'\t'+result+'\t'+result+"\tfalse\t"+review;
	}
	private static void write(String path, String... rows){
		assert(rows.length>0):"Fixture files always include their declared schema or sequence";
		final ByteStreamWriter out=RrnaResourceIO.writer(path);try{for(String row:rows){out.println(row);}}finally{RrnaResourceIO.finish(out);}
	}
	private static ArrayList<String> lines(String path){
		final ByteFile in=RrnaResourceIO.open(path);final ArrayList<String> rows=new ArrayList<String>();
		try{for(byte[] row; (row=in.nextLine())!=null; ){rows.add(new String(row, java.nio.charset.StandardCharsets.UTF_8));}}
		finally{require(!in.close(), "Fixture result read failed");}return rows;
	}
	private static void require(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
