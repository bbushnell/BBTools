package prok;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;

import prok.NcrnaGffFamilyGrader.Locus;
import prok.LegacyRrnaGffGrader.GradeResult;

/** The nine minimal synthetic fixture cases from
 * legacy_rrna_comparator_interface_spec_20260902.md section 6, one test method per
 * table row, each isolating exactly one invariant. Mirrors NcrnaGffFamilyGraderTest's
 * pure fixture-testing pattern (grade() takes plain in-memory Locus lists, no file I/O)
 * for eight of the nine; fixture #7 additionally exercises the real file-reading path
 * (loadCallsFromGff), not only the pure extraction function -- Citan's review correction,
 * 2026-09-02: a fixture that only tests extractLegacyFamily never proves the loader
 * actually surfaces this as a fail-loud error when reading a real file. */
public class LegacyRrnaGffGraderTest {

	public static void main(String[] args) throws Exception {
		testSimpleTP();
		testSimpleFN();
		testSimpleFP();
		testDuplicateQueryLocus();
		testCrossFamilySubtypeConfusion();
		test18SScored();
		test18SNegativeSeqid();
		testUnrecognizedSubtypeTokenFailsLoud();
		testReciprocalOverlapBoundary();
		testStrandMismatch();
		testFamilyExtraction();
		testEndToEndFileIO();
		System.out.println("PASS LegacyRrnaGffGraderTest");
	}

	/** Spec fixture #1: one call cleanly overlapping its own truth locus (reciprocal
	 * overlap well above 50% both ways) -> TP for both rows, linked match IDs. */
	private static void testSimpleTP(){
		final ArrayList<Locus> truth=new ArrayList<>();
		truth.add(new Locus("ctg1", "16S", 105, 495, '+'));
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("ctg1", "16S", 100, 500, '+'));

		final GradeResult r=LegacyRrnaGffGrader.grade(truth, calls, new HashSet<String>());
		check("16S tp", 1, r.byFamily.get("16S").tp);
		check("16S fn", 0, r.byFamily.get("16S").fn);
		check("16S fp", 0, r.byFamily.get("16S").fp);
		check("unionTP", 1, r.unionTP);
		check("unionFN", 0, r.unionFN);
		check("unionFP", 0, r.unionFP);
	}

	/** Spec fixture #2: one truth locus with no overlapping call at all -> FN. */
	private static void testSimpleFN(){
		final ArrayList<Locus> truth=new ArrayList<>();
		truth.add(new Locus("ctg1", "23S", 1000, 3000, '+'));
		final ArrayList<Locus> calls=new ArrayList<>();

		final GradeResult r=LegacyRrnaGffGrader.grade(truth, calls, new HashSet<String>());
		check("23S tp", 0, r.byFamily.get("23S").tp);
		check("23S fn", 1, r.byFamily.get("23S").fn);
		check("23S fp", 0, r.byFamily.get("23S").fp);
		check("unionFN", 1, r.unionFN);
	}

	/** Spec fixture #3: one query call with no overlapping truth at all -> FP. Uses a
	 * short (<300bp) 5S-family locus to also exercise the 5S length gate incidentally. */
	private static void testSimpleFP(){
		final ArrayList<Locus> truth=new ArrayList<>();
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("ctg1", "5S", 50, 150, '+'));//length 101, <300, valid 5S

		final GradeResult r=LegacyRrnaGffGrader.grade(truth, calls, new HashSet<String>());
		check("5S tp", 0, r.byFamily.get("5S").tp);
		check("5S fn", 0, r.byFamily.get("5S").fn);
		check("5S fp", 1, r.byFamily.get("5S").fp);
		check("5S dup", 0, r.byFamily.get("5S").dup);
		check("unionFP", 1, r.unionFP);
	}

	/** Spec fixture #4: two calls both overlapping ONE truth locus -- direct regression
	 * test for the CompareGff.overlapMetrics non-consuming defect this whole grader
	 * lineage exists to avoid. Must produce exactly one TP + one DUP, never two TPs and
	 * never a plain FP for the loser. */
	private static void testDuplicateQueryLocus(){
		final ArrayList<Locus> truth=new ArrayList<>();
		truth.add(new Locus("ctg1", "16S", 1000, 1100, '-'));
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("ctg1", "16S", 1000, 1100, '-'));//exact match, highest overlap -> wins TP
		calls.add(new Locus("ctg1", "16S", 1005, 1095, '-'));//also overlaps SAME locus -> DUP

		final GradeResult r=LegacyRrnaGffGrader.grade(truth, calls, new HashSet<String>());
		check("16S tp", 1, r.byFamily.get("16S").tp);
		check("16S fp", 1, r.byFamily.get("16S").fp);
		check("16S dup", 1, r.byFamily.get("16S").dup);
		check("16S fn", 0, r.byFamily.get("16S").fn);
		check("unionTP", 1, r.unionTP);
		check("unionFP", 1, r.unionFP);
	}

	/** Spec fixture #5: a call attributed the WRONG subtype overlapping a truth locus of
	 * a different in-scope family at the same coordinates -- a real caller subtype error.
	 * Must show FN for the true family, FP for the wrong-family call, AND both rows
	 * flagged by the cross-label diagnostic with linked match IDs -- while the
	 * family-blind union pass still correctly counts it as a detection (something WAS
	 * found there, just mislabeled). */
	private static void testCrossFamilySubtypeConfusion(){
		final ArrayList<Locus> truth=new ArrayList<>();
		truth.add(new Locus("ctg1", "23S", 1010, 2490, '+'));
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("ctg1", "16S", 1000, 2500, '+'));//wrong family, overlaps the 23S truth

		final GradeResult r=LegacyRrnaGffGrader.grade(truth, calls, new HashSet<String>());
		check("23S tp", 0, r.byFamily.get("23S").tp);
		check("23S fn", 1, r.byFamily.get("23S").fn);
		check("16S tp", 0, r.byFamily.get("16S").tp);
		check("16S fp", 1, r.byFamily.get("16S").fp);
		check("unionTP", 1, r.unionTP);//family-blind pooled pass: the 16S call DOES match the 23S truth
		check("unionFN", 0, r.unionFN);
		check("crossLabel_23S", 1, r.crossLabelCounts.get("23S"));
		check("crossLabel_16S", 0, r.crossLabelCounts.get("16S"));
	}

	/** An 18S truth locus and matching call are ordinary scored true positives. */
	private static void test18SScored(){
		final ArrayList<Locus> truth=new ArrayList<>();
		truth.add(new Locus("ctg1", "18S", 205, 1795, '+'));
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("ctg1", "18S", 200, 1800, '+'));

		final GradeResult r=LegacyRrnaGffGrader.grade(truth, calls, new HashSet<String>());
		check("18S tp", 1, r.byFamily.get("18S").tp);
		check("18S fn", 0, r.byFamily.get("18S").fn);
		check("18S fp", 0, r.byFamily.get("18S").fp);
		check("unionTP", 1, r.unionTP);
		check("unionFN", 0, r.unionFN);
		check("unionFP", 0, r.unionFP);
	}

	/** 18S participates in the ordinary negative-sequence FP path. */
	private static void test18SNegativeSeqid(){
		final ArrayList<Locus> truth=new ArrayList<>();
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("neg18", "18S", 100, 1700, '+'));
		final HashSet<String> neg=new HashSet<>(); neg.add("neg18");
		final GradeResult r=LegacyRrnaGffGrader.grade(truth, calls, neg);
		check("18S neg fp", 1, r.byFamily.get("18S").fp);
		check("18S neg union fp", 1, r.unionFP);
		check("18S neg seqid fp", 1, r.negFP);
	}

	/** Spec fixture #7: an rRNA-typed GFF row whose attribute text matches none of
	 * {16S,23S,18S,5S} (or matches "5S" text but at length>=300, so the length gate
	 * correctly rejects it as 5S) must fail loud, never silently skip. Checked TWO ways:
	 * directly against extractLegacyFamily (the pure function, cheap and exact), AND
	 * against the real loadCallsFromGff file-reading path via a minimal temporary
	 * malformed rRNA GFF -- the pure-function check alone never proves the loader itself
	 * surfaces this as IllegalStateException when reading an actual file (Citan's review
	 * correction, 2026-09-02). */
	private static void testUnrecognizedSubtypeTokenFailsLoud() throws Exception {
		if(LegacyRrnaGffGrader.extractLegacyFamily("RNA,startScr:0,stopScr:0", 400)!=null){
			throw new AssertionError("Expected null (no recognized subtype token) but got a family");
		}
		//"5S" text present, but length>=300 -> the length gate must still reject it as 5S
		if(LegacyRrnaGffGrader.extractLegacyFamily("5S,startScr:1,len:450", 450)!=null){
			throw new AssertionError("Expected null (5S text present but length>=300 fails the gate) but got a family");
		}

		final File tmp=new File(System.getProperty("java.io.tmpdir"),
			"LegacyRrnaGffGraderTest_malformed_"+System.nanoTime()+".gff");
		try{
			final PrintWriter pw=new PrintWriter(tmp);
			pw.println("##gff-version 3");
			pw.println("ctg1\tBBTools\trRNA\t100\t200\t50.0\t+\t0\tRNA,startScr:0,stopScr:0");
			pw.close();
			boolean threw=false;
			try{
				LegacyRrnaGffGrader.loadCallsFromGff(tmp.getAbsolutePath());
			}catch(IllegalStateException expected){threw=true;}
			if(!threw){
				throw new AssertionError("loadCallsFromGff accepted a malformed rRNA row instead of failing loud");
			}
		}finally{
			tmp.delete();
		}
	}

	/** Spec fixture #8: overlap is 100% of the CALL's length but only 10% of the TRUTH
	 * locus's length -- must NOT match, proving the reciprocal (both-lengths) requirement
	 * is enforced, not just one-directional coverage. */
	private static void testReciprocalOverlapBoundary(){
		final ArrayList<Locus> truth=new ArrayList<>();
		truth.add(new Locus("ctg1", "16S", 100, 1099, '+'));//length 1000
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("ctg1", "16S", 100, 199, '+'));//length 100, fully inside truth: overlap=100 (100% of call, 10% of truth)

		final GradeResult r=LegacyRrnaGffGrader.grade(truth, calls, new HashSet<String>());
		check("16S tp", 0, r.byFamily.get("16S").tp);
		check("16S fn", 1, r.byFamily.get("16S").fn);
		check("16S fp", 1, r.byFamily.get("16S").fp);
	}

	/** Spec fixture #9: identical coordinates, opposite strand -- must NOT match despite
	 * perfect coordinate overlap, proving the same-strand requirement is enforced. */
	private static void testStrandMismatch(){
		final ArrayList<Locus> truth=new ArrayList<>();
		truth.add(new Locus("ctg1", "16S", 100, 500, '+'));
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("ctg1", "16S", 100, 500, '-'));

		final GradeResult r=LegacyRrnaGffGrader.grade(truth, calls, new HashSet<String>());
		check("16S tp", 0, r.byFamily.get("16S").tp);
		check("16S fn", 1, r.byFamily.get("16S").fn);
		check("16S fp", 1, r.byFamily.get("16S").fp);
	}

	/** Direct unit coverage of extractLegacyFamily's check order and the 5S length gate,
	 * beyond the fail-loud case above -- confirms each real subtype resolves correctly. */
	private static void testFamilyExtraction(){
		checkFamily("16S,startScr:1,len:1541", 1541, "16S");
		checkFamily("23S,startScr:1,len:2932", 2932, "23S");
		checkFamily("18S,startScr:1,len:1800", 1800, "18S");
		checkFamily("5S,startScr:1,len:111", 111, "5S");
		//5S text present but length right at the boundary: 299 passes, 300 fails (length<300, strict)
		checkFamily("5S,startScr:1,len:299", 299, "5S");
		checkFamily("5S,startScr:1,len:300", 300, null);
	}

	private static void checkFamily(String attrs, int length, String expected){
		final String observed=LegacyRrnaGffGrader.extractLegacyFamily(attrs, length);
		final boolean ok=(expected==null ? observed==null : expected.equals(observed));
		if(!ok){throw new AssertionError("For attrs='"+attrs+"' len="+length+": expected "+expected+", observed "+observed);}
	}

	/** End-to-end file-I/O reproduction of the hand-run smoke test from this session's
	 * implementation report -- exercises loadTruthTsv + loadCallsFromGff + writeReport
	 * together (not just the pure grade() core the other nine fixtures use), and asserts
	 * on the ACTUAL WRITTEN .summary.tsv/.rawcalls.tsv file contents, not just in-memory
	 * GradeResult fields, so a bug in writeReport's formatting would be caught too.
	 * Citan's durability request, 2026-09-02: fold the hand-run check into the automated
	 * suite so it's reproducible, not a one-off manual verification. The fixture has clean
	 * 16S and 18S TPs, a 23S truth locus with no call (FN),
	 * and an unrelated CDS row (must be silently skipped, never appear in the output). */
	private static void testEndToEndFileIO() throws Exception {
		final File dir=new File(System.getProperty("java.io.tmpdir"),
			"LegacyRrnaGffGraderTest_e2e_"+System.nanoTime());
		dir.mkdirs();
		final File gffFile=new File(dir, "calls.gff");
		final File truthFile=new File(dir, "truth.tsv");
		final String outPrefix=new File(dir, "out").getAbsolutePath();
		final File summaryFile=new File(outPrefix+".summary.tsv");
		final File rawFile=new File(outPrefix+".rawcalls.tsv");
		try{
			PrintWriter pw=new PrintWriter(gffFile);
			pw.println("##gff-version 3");
			pw.println("ctg1\tBBTools\trRNA\t100\t500\t582.34\t+\t0\t16S,startScr:1,len:401");
			pw.println("ctg1\tBBTools\trRNA\t205\t1795\t100.0\t+\t0\t18S,startScr:1,len:1591");
			pw.println("ctg2\tBBTools\tCDS\t10\t100\t50.0\t+\t0\tCDS,fr0,startScr:0,stopScr:0");
			pw.close();

			pw=new PrintWriter(truthFile);
			pw.println("seqid\tfamily\tstart\tstop\tstrand");
			pw.println("ctg1\t16S\t105\t495\t+");
			pw.println("ctg1\t18S\t200\t1800\t+");
			pw.println("ctg3\t23S\t1000\t2000\t+");
			pw.close();

			final ArrayList<Locus> truth=LegacyRrnaGffGrader.loadTruthTsv(truthFile.getAbsolutePath());
			final ArrayList<Locus> calls=LegacyRrnaGffGrader.loadCallsFromGff(gffFile.getAbsolutePath());
			final GradeResult r=LegacyRrnaGffGrader.grade(truth, calls, new HashSet<String>());
			LegacyRrnaGffGrader.writeReport(outPrefix, r);

			//In-memory sanity check first -- cheap, and pinpoints loader-vs-writer if
			//the file-level assertions below ever fail.
			check("16S tp", 1, r.byFamily.get("16S").tp);
			check("23S fn", 1, r.byFamily.get("23S").fn);
			check("18S tp", 1, r.byFamily.get("18S").tp);
			check("unionTP", 2, r.unionTP);

			final String summary=readFile(summaryFile);
			if(!summary.contains("16S\t1\t0\t0\t1.000000\t1.000000")){
				throw new AssertionError("summary.tsv missing expected 16S row:\n"+summary);
			}
			if(!summary.contains("23S\t0\t1\t0")){
				throw new AssertionError("summary.tsv missing expected 23S row:\n"+summary);
			}
			if(!summary.contains("18S\t1\t0\t0\t1.000000\t1.000000")){
				throw new AssertionError("summary.tsv missing expected 18S row:\n"+summary);
			}

			final String raw=readFile(rawFile);
			int scored18SRows=0;
			for(String line : raw.split("\n")){
				if(line.contains("\t18S\t") && line.contains("\tTP\t")){scored18SRows++;}
				if(line.contains("ctg2")){
					throw new AssertionError("rawcalls.tsv unexpectedly contains a row for the "
						+"CDS-only seqid ctg2 -- non-rRNA rows must be silently skipped: "+line);
				}
			}
			if(scored18SRows!=2){//exactly one TRUTH row + one CALL row
				throw new AssertionError("Expected exactly 2 scored 18S rawcalls.tsv rows "
					+"(1 truth + 1 call), found "+scored18SRows+":\n"+raw);
			}
		}finally{
			gffFile.delete(); truthFile.delete(); summaryFile.delete(); rawFile.delete(); dir.delete();
		}
	}

	private static String readFile(File f) throws Exception {
		final StringBuilder sb=new StringBuilder();
		final BufferedReader br=new BufferedReader(new FileReader(f));
		String line;
		while((line=br.readLine())!=null){sb.append(line).append("\n");}
		br.close();
		return sb.toString();
	}

	private static void check(String label, long expected, long observed){
		if(expected!=observed){
			throw new AssertionError("Expected "+label+"="+expected+", observed "+observed);
		}
	}
}
