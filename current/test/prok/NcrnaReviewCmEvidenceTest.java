package test.prok;

import java.io.File;
import java.util.ArrayList;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import prok.NcrnaReviewCmEvidence;
import prok.NcrnaReviewUnion;
import prok.RrnaResourceIO;

/** Checks independent CM attribution, projection and conservation without truth edits.
 * @author Brian, Raiden
 */
public final class NcrnaReviewCmEvidenceTest {
	public static void main(String[] args){
		require(args.length==1, "Supply a fresh fixture directory");final File root=new File(args[0]);
		require(!root.exists() && root.mkdirs(), "Fixture must not overwrite existing evidence");final String base=root.getAbsolutePath();
		write(base+"/loci.tsv", NcrnaReviewUnion.LOCUS_HEADER, locus("r0", 101, "+"), locus("r1", 201, "-"),
			locus("r2", 301, "+"), locus("r3", 401, "+"), locus("r4", 501, "+"));
		final String[] footer={"# Program: cmsearch", "# Version: 1.1.5 (Sep 2023)", "# [ok]"};
		write(base+"/hits.tbl", hit("r0", 1, 6, 25, "+", 40), hit("r0", 1, 6, 25, "+", 30), hit("r0", 1, 1, 4, "+", 40),
			hit("r1", 1, 6, 25, "+", 40), hit("r1", 1, 7, 25, "+", 39), hit("r2", 1, 25, 6, "-", 40), hit("r3", 2, 6, 25, "+", 40),
			footer[0], footer[1], footer[2]);
		run(base, "hits.tbl", "evidence");
		final ArrayList<String> summary=lines(base+"/evidence.summary.tsv"), raw=lines(base+"/evidence.hits.tsv"), loci=lines(base+"/evidence.loci.tsv");
		require(summary.size()==2 && summary.get(1).equals("5\t7\t3\t1\t1\tCM_EVIDENCE_NOT_TRUTH_EDITS_OR_FP_VERDICTS"), "Below-GA, neighbor, opposite-strand and partial-model evidence must not become full qualifying hits");
		require(raw.size()==8 && raw.get(4).startsWith("r1\t4\tchr\t201\t220\t-\t"), "Sense CM coordinates must project back through a minus-strand context");
		require(raw.get(6).startsWith("r2\t6\tchr\t301\t320\t-\t") && raw.get(6).endsWith("\t0\tfalse"), "Opposite-strand hit stays in raw evidence but does not confirm the predicted sense");
		require(loci.get(1).equals("r0\t3\t2\t1\tONE_FULL_GA_OVERLAPPING_HIT") && loci.get(2).endsWith("MULTIPLE_FULL_GA_OVERLAPPING_HITS"), "Every hit and ambiguous multiple-support case must remain explicit");
		write(base+"/foreign.tbl", hit("foreign", 1, 6, 25, "+", 40), footer[0], footer[1], footer[2]);
		boolean rejected=false;try{run(base, "foreign.tbl", "bad");}catch(IllegalArgumentException e){rejected=e.getMessage().contains("frozen review context");}
		require(rejected && !new File(base+"/bad.summary.tsv").exists(), "Foreign context evidence must fail before reporting a complete audit");
		System.out.println("NCRNA_REVIEW_CM_EVIDENCE_TEST_PASS");
	}
	private static void run(String base, String hits, String out){
		assert(!hits.isEmpty()):"Fixture explicitly binds the CM table under test";
		NcrnaReviewCmEvidence.main(new String[]{"loci="+base+"/loci.tsv", "hits="+base+'/'+hits, "outprefix="+base+'/'+out,
			"expected=5", "modelid=RF00174", "family=cobalamin", "clen=20", "gabits=31.5"});
	}
	private static String locus(String id, int start, String strand){
		assert(start>5 && start+24<1000):"Synthetic context has five real flanking bases at both ends";
		return id+"\tg\ta\tfixture.fa\tchr\tcobalamin\t"+strand+'\t'+start+'\t'+(start+19)+"\t1000\t5\t5\t30\tPREDICTION_UNRESOLVED";
	}
	private static String hit(String id, int from, int start, int stop, String strand, int score){
		assert(from>=1 && from<=20):"Fixture model span is within the declared 20-position test CM";
		return id+" - Cobalamin RF00174 cm "+from+" 20 "+start+' '+stop+' '+strand+" no 1 0.5 0.0 "+score+" 1e-12 ! fixture";
	}
	private static void write(String path, String... rows){
		assert(rows.length>0):"Every fixture input has an explicit schema or completion footer";
		final ByteStreamWriter out=RrnaResourceIO.writer(path);try{for(String row:rows){out.println(row);}}finally{RrnaResourceIO.finish(out);}
	}
	private static ArrayList<String> lines(String path){
		final ByteFile in=RrnaResourceIO.open(path);final ArrayList<String> rows=new ArrayList<String>();
		try{for(byte[] row; (row=in.nextLine())!=null; ){rows.add(new String(row, java.nio.charset.StandardCharsets.UTF_8));}}
		finally{require(!in.close(), "Fixture result read failed");}return rows;
	}
	private static void require(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
