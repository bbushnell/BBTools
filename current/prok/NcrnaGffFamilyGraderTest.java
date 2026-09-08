package prok;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import prok.NcrnaGffFamilyGrader.GradeResult;
import prok.NcrnaGffFamilyGrader.Locus;

/** Focused tests for NcrnaGffFamilyGrader's matching+grading core and its input
 * validation, mirroring NcrnaCombinedGradingDriverTest's fixture-testing pattern for
 * the pure logic, plus small real-file tests for the loader validation added in the
 * 2026-09-02 corrections round (CRLF tolerance, strand/field-count fail-loud). */
public class NcrnaGffFamilyGraderTest {

	public static void main(String[] args) throws Exception {
		testCleanTpFnFp();
		testDuplicateQueryLocus();
		testCrossFamilyMatch();
		testModelPrefixExtraction();
		testAttributionGapFailsLoud();
		testNegSeqidFP();
		testMaximumCardinalityCounterexample();
		testCrlfTolerance();
		testMalformedStrandFailsLoud();
		testWrongFieldCountFailsLoud();
		testAuditTrailReconstruction();
		testModelTokenBoundary();
		testLeadingTruthComments();
		testEmptyTruth();
		testDuplicateTruthFailsLoud();
		testMalformedGffFieldCountFailsLoud();
		testMissingModelFailsLoud();
		testCrossLabelManyToOneNotOverwritten();
		testEndToEndReconstruction();
		testGenericFamilyMapThreeFamilyFixture();
		testFamilymapWrongFieldCountFailsLoud();
		testFamilymapEmptyFieldFailsLoud();
		testFamilymapDuplicateFamilyFailsLoud();
		testFamilymapDuplicatePrefixFailsLoud();
		testEmptyFamilymapFailsLoud();
		testGenericPrefixAmbiguityFailsLoud();
		testBoundaryOffsetFormulas();
		testBoundaryAggregateBlock();
		System.out.println("PASS NcrnaGffFamilyGraderTest");
	}

	/** One clean TP (call overlaps its own truth >=50% reciprocal), one truth locus
	 * with no matching call (FN), one call with no matching truth (FP) -- on the SAME
	 * family, different seqid, so they cannot interact with each other. */
	private static void testCleanTpFnFp(){
		final ArrayList<Locus> truth=new ArrayList<>();
		truth.add(new Locus("ctg1", "RF00013", 100, 200, '+'));//will be matched -> TP
		truth.add(new Locus("ctg2", "RF00013", 100, 200, '+'));//no call -> FN
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("ctg1", "RF00013", 105, 195, '+'));//overlaps truth on ctg1 -> TP
		calls.add(new Locus("ctg3", "RF00013", 100, 200, '+'));//no truth on ctg3 -> FP

		final GradeResult r=NcrnaGffFamilyGrader.grade(truth, calls, new HashSet<String>());
		check("tpA", 1, r.tpA);
		check("fnA", 1, r.fnA);
		check("fpA", 1, r.fpA);
		check("dupA", 0, r.dupA);
		check("unionTP", 1, r.unionTP);
		check("unionFN", 1, r.unionFN);
		check("unionFP", 1, r.unionFP);
	}

	/** Two calls both overlapping ONE truth locus -- the direct regression test for the
	 * CompareGff.overlapMetrics non-consuming defect this grader is built to avoid.
	 * Must produce exactly one TP + one duplicate, NEVER two TPs. */
	private static void testDuplicateQueryLocus(){
		final ArrayList<Locus> truth=new ArrayList<>();
		truth.add(new Locus("ctg1", "RF01685", 1000, 1100, '-'));
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("ctg1", "RF01685", 1000, 1100, '-'));//exact match, highest overlap -> wins the TP
		calls.add(new Locus("ctg1", "RF01685", 1005, 1095, '-'));//also overlaps the SAME locus -> duplicate, FP

		final GradeResult r=NcrnaGffFamilyGrader.grade(truth, calls, new HashSet<String>());
		check("tpB", 1, r.tpB);
		check("fpB", 1, r.fpB);
		check("dupB", 1, r.dupB);
		check("fnB", 0, r.fnB);
		check("unionTP", 1, r.unionTP);
		check("unionFP", 1, r.unionFP);
	}

	/** A family-B call overlapping a family-A truth locus (wrong subtype label) --
	 * must score as: union TP (something real WAS detected there), subtype-A FN (never
	 * got the correct label), subtype-B FP (fired under the wrong label), AND is
	 * surfaced via the crossAasB diagnostic -- never silently absorbed. */
	private static void testCrossFamilyMatch(){
		final ArrayList<Locus> truth=new ArrayList<>();
		truth.add(new Locus("ctg1", "RF00013", 500, 600, '+'));
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("ctg1", "RF01685", 505, 595, '+'));//wrong family, but geometrically overlaps the RF00013 truth

		final GradeResult r=NcrnaGffFamilyGrader.grade(truth, calls, new HashSet<String>());
		check("tpA", 0, r.tpA);
		check("fnA", 1, r.fnA);
		check("tpB", 0, r.tpB);
		check("fpB", 1, r.fpB);
		check("unionTP", 1, r.unionTP);//family-blind pooled pass: the B call DOES match the A truth
		check("unionFN", 0, r.unionFN);
		check("crossAasB", 1, r.crossAasB);
		check("crossBasA", 0, r.crossBasA);
	}

	private static void testModelPrefixExtraction(){
		checkFamily("startScr:0,stopScr:0,innerScr:0,len:186,model:6S_RF00013_unknown_c1113 n=5",
			NcrnaGffFamilyGrader.FAM_A);
		checkFamily("len:109,model:6S_RF01685_unknown_c19 n=7", NcrnaGffFamilyGrader.FAM_B);
		//A different, non-sixS generic ncRNA family's model -- must resolve to null (silent skip),
		//not one of the two sixS families.
		final String rnasepAttrs="len:380,model:rnasep_consensus_1";
		final int mi=rnasepAttrs.indexOf("model:");
		final String fam=NcrnaGffFamilyGrader.extractFamily(rnasepAttrs, mi);
		if(fam!=null){throw new AssertionError("Expected null (non-sixS family) but got "+fam);}
	}

	private static void checkFamily(String attrs, String expectedFamily){
		final int mi=attrs.indexOf("model:");
		if(mi<0){throw new AssertionError("No model: token found in: "+attrs);}
		final String observed=NcrnaGffFamilyGrader.extractFamily(attrs, mi);
		if(!expectedFamily.equals(observed)){
			throw new AssertionError("Expected "+expectedFamily+", observed "+observed+" for: "+attrs);
		}
	}

	/** A truth locus on a seqid ALSO listed in negSeqids is a data-consistency
	 * violation -- must fail loud, not silently grade around it. */
	private static void testAttributionGapFailsLoud(){
		final ArrayList<Locus> truth=new ArrayList<>();
		truth.add(new Locus("ctg1", "RF00013", 100, 200, '+'));
		final HashSet<String> negSeqids=new HashSet<>();
		negSeqids.add("ctg1");
		boolean threw=false;
		try{
			NcrnaGffFamilyGrader.grade(truth, new ArrayList<Locus>(), negSeqids);
		}catch(IllegalStateException expected){threw=true;}
		if(!threw){throw new AssertionError("Truth locus on a negseqids-listed seqid was accepted, not rejected");}
	}

	/** A call on a negseqids-listed seqid is an unmatchable-by-construction straight FP
	 * (and does not require the grouping/matching path at all). */
	private static void testNegSeqidFP(){
		final ArrayList<Locus> truth=new ArrayList<>();
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("negctg", "RF00013", 100, 200, '+'));
		final HashSet<String> negSeqids=new HashSet<>();
		negSeqids.add("negctg");

		final GradeResult r=NcrnaGffFamilyGrader.grade(truth, calls, negSeqids);
		check("fpA", 1, r.fpA);
		check("negFP", 1, r.negFP);
		check("unionFP", 1, r.unionFP);
		check("tpA", 0, r.tpA);
	}

	/** THE regression test for the greedy-vs-maximum-cardinality correction (Citan's
	 * review, 2026-09-02). Every overlap below is hand-verified against the actual
	 * reciprocal-50%-of-BOTH-lengths rule (not just plausible-looking coordinates):
	 * T1=[100,109] len10, T2=[105,114] len10 (T1/T2 overlap each other by 5bp -- fine
	 * for a synthetic matching-algorithm test, no biological realism required).
	 * C1=[102,111] len10: overlap-with-T1 = min(111,109)-max(102,100)+1 = 8 (80% of
	 * both len-10 loci, passes); overlap-with-T2 = min(111,114)-max(102,105)+1 = 7 (70%
	 * of both, passes, but LOWER than the T1 overlap -- C1's own edges sort
	 * T1(8),T2(7)). C2=[95,104] len10: overlap-with-T1 = min(104,109)-max(95,100)+1 = 5
	 * (exactly 50% of both len-10 loci -- passes at the boundary); overlap-with-T2 =
	 * min(104,114)-max(95,105)+1 = 0 (no overlap) -- C2's ONLY edge is to T1.
	 *
	 * <p>A global greedy highest-overlap-first matcher sorts all candidate pairs
	 * (C1-T1:8, C1-T2:7, C2-T1:5) and assigns C1-T1 first (the single highest), which
	 * consumes both C1 and T1 -- C1-T2 and C2-T1 are both now unusable (C1 and T1 gone),
	 * leaving C2 stranded: only 1 total match. Kuhn's algorithm re-routes correctly:
	 * C1 first claims T1 (its top edge), then when C2 needs T1, the augmenting search
	 * finds that C1 can instead take T2 (its second edge), freeing T1 for C2 -- both
	 * truth loci and both calls end up matched, 2 total, the true maximum cardinality. */
	private static void testMaximumCardinalityCounterexample(){
		final ArrayList<Locus> truth=new ArrayList<>();
		truth.add(new Locus("ctg1", "RF00013", 100, 109, '+'));//T1
		truth.add(new Locus("ctg1", "RF00013", 105, 114, '+'));//T2
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("ctg1", "RF00013", 102, 111, '+'));//C1: edges to T1(overlap8) and T2(overlap7)
		calls.add(new Locus("ctg1", "RF00013", 95, 104, '+'));//C2: edge to T1(overlap5) only

		final GradeResult r=NcrnaGffFamilyGrader.grade(truth, calls, new HashSet<String>());
		//Maximum cardinality is 2 (e.g. C1-T2, C2-T1): both truth loci matched, both
		//calls matched, zero FN, zero FP. A greedy highest-overlap-first matcher would
		//instead grab C1-T1 (overlap 8, the single highest) first and strand C2 -- only
		//1 total match. Assert the FULL maximum here.
		check("tpA", 2, r.tpA);
		check("fnA", 0, r.fnA);
		check("fpA", 0, r.fpA);
		check("dupA", 0, r.dupA);
	}

	/** CRLF line endings are a convention difference, not malformed data -- every
	 * loader must tolerate them. */
	private static void testCrlfTolerance() throws Exception {
		final File f=File.createTempFile("ncrnagraider_truth_crlf", ".tsv");
		f.deleteOnExit();
		try(FileWriter w=new FileWriter(f)){
			w.write("seqid\tfamily\tstart\tstop\tstrand\r\n");
			w.write("ctgX\tRF00013\t10\t20\t+\r\n");
		}
		final ArrayList<Locus> truth=NcrnaGffFamilyGrader.loadTruthTsv(f.getAbsolutePath());
		if(truth.size()!=1){throw new AssertionError("Expected 1 truth locus from CRLF file, got "+truth.size());}
		final Locus l=truth.get(0);
		if(!"ctgX".equals(l.seqid) || !"RF00013".equals(l.family) || l.start!=10 || l.stop!=20 || l.strand!='+'){
			throw new AssertionError("CRLF-tolerant parse produced wrong values: "+l.seqid+" "+l.family
				+" "+l.start+"-"+l.stop+" "+l.strand);
		}
	}

	/** An invalid strand character must fail loud, not default to something silently. */
	private static void testMalformedStrandFailsLoud() throws Exception {
		final File f=File.createTempFile("ncrnagraider_truth_badstrand", ".tsv");
		f.deleteOnExit();
		try(FileWriter w=new FileWriter(f)){
			w.write("seqid\tfamily\tstart\tstop\tstrand\n");
			w.write("ctgX\tRF00013\t10\t20\tX\n");
		}
		boolean threw=false;
		try{
			NcrnaGffFamilyGrader.loadTruthTsv(f.getAbsolutePath());
		}catch(IllegalArgumentException expected){threw=true;}
		if(!threw){throw new AssertionError("Invalid strand character 'X' was accepted, not rejected");}
	}

	/** A truth row with the wrong field count (4, not exactly 5) must fail loud. */
	private static void testWrongFieldCountFailsLoud() throws Exception {
		final File f=File.createTempFile("ncrnagraider_truth_badfields", ".tsv");
		f.deleteOnExit();
		try(FileWriter w=new FileWriter(f)){
			w.write("seqid\tfamily\tstart\tstop\tstrand\n");
			w.write("ctgX\tRF00013\t10\t20\n");//missing strand field
		}
		boolean threw=false;
		try{
			NcrnaGffFamilyGrader.loadTruthTsv(f.getAbsolutePath());
		}catch(IllegalArgumentException expected){threw=true;}
		if(!threw){throw new AssertionError("Truth row with 4 fields (not exactly 5) was accepted, not rejected");}
	}

	/** Round-2 auditability: every aggregate in the summary must be mechanically
	 * reconstructable by scanning rawcalls alone. Reruns testCleanTpFnFp's exact
	 * scenario and (a) reconstructs tpA/fnA/fpA purely by grouping auditRows on
	 * role/family/subtypeOutcome, confirming they match r's own counters, and (b)
	 * checks the matched pair's subtypeMatchId links point at each other's locusId. */
	private static void testAuditTrailReconstruction(){
		final ArrayList<Locus> truth=new ArrayList<>();
		truth.add(new Locus("ctg1", "RF00013", 100, 200, '+'));//matched -> TP
		truth.add(new Locus("ctg2", "RF00013", 100, 200, '+'));//FN
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("ctg1", "RF00013", 105, 195, '+'));//matched -> TP
		calls.add(new Locus("ctg3", "RF00013", 100, 200, '+'));//FP

		final GradeResult r=NcrnaGffFamilyGrader.grade(truth, calls, new HashSet<String>());

		long reconTpA=0, reconFnA=0, reconFpA=0;
		NcrnaGffFamilyGrader.AuditRow matchedTruthRow=null, matchedCallRow=null;
		for(NcrnaGffFamilyGrader.AuditRow ar : r.auditRows){
			if(!"RF00013".equals(ar.locus.family)){continue;}
			if(ar.role.equals("TRUTH")){
				if(ar.subtypeOutcome.equals("TP")){reconTpA++; matchedTruthRow=ar;}
				else if(ar.subtypeOutcome.equals("FN")){reconFnA++;}
			}else{//CALL
				if(ar.subtypeOutcome.equals("TP")){matchedCallRow=ar;}
				else if(ar.subtypeOutcome.equals("FP") || ar.subtypeOutcome.equals("DUP")){reconFpA++;}
			}
		}
		check("reconstructed tpA", r.tpA, reconTpA);
		check("reconstructed fnA", r.fnA, reconFnA);
		check("reconstructed fpA", r.fpA, reconFpA);
		if(matchedTruthRow==null || matchedCallRow==null){
			throw new AssertionError("Expected one matched TRUTH row and one matched CALL row in auditRows");
		}
		if(!matchedTruthRow.subtypeMatchId.equals(matchedCallRow.locusId)){
			throw new AssertionError("Matched TRUTH row's subtypeMatchId ("+matchedTruthRow.subtypeMatchId
				+") does not point at the matched CALL row's own locusId ("+matchedCallRow.locusId+")");
		}
		if(!matchedCallRow.subtypeMatchId.equals(matchedTruthRow.locusId)){
			throw new AssertionError("Matched CALL row's subtypeMatchId ("+matchedCallRow.subtypeMatchId
				+") does not point at the matched TRUTH row's own locusId ("+matchedTruthRow.locusId+")");
		}
	}

	/** Round-2 hardening: a bare prefix match is not enough. "6S_RF000135_x" shares
	 * "6S_RF00013" as a literal prefix but is NOT the RF00013 family -- the character
	 * right after the prefix ('5', a digit) is not a real token boundary ('_' or ',')
	 * and must be rejected. A model name that happens to END exactly at the prefix
	 * (nothing follows) is a valid boundary case and must still match. */
	private static void testModelTokenBoundary(){
		final String badAttrs="len:186,model:6S_RF000135_unknown_c1";
		final int mi1=badAttrs.indexOf("model:");
		final String fam1=NcrnaGffFamilyGrader.extractFamily(badAttrs, mi1);
		if(fam1!=null){
			throw new AssertionError("6S_RF000135_x (digit right after the prefix) must NOT match "
				+"RF00013 -- observed "+fam1);
		}

		final String exactAttrs="len:186,model:6S_RF00013";//prefix runs to end of string
		final int mi2=exactAttrs.indexOf("model:");
		final String fam2=NcrnaGffFamilyGrader.extractFamily(exactAttrs, mi2);
		if(!NcrnaGffFamilyGrader.FAM_A.equals(fam2)){
			throw new AssertionError("model:6S_RF00013 with nothing following is a valid boundary "
				+"case and must match RF00013 -- observed "+fam2);
		}

		final String commaAttrs="model:6S_RF01685,len:109";//boundary is a comma, not underscore
		final int mi3=commaAttrs.indexOf("model:");
		final String fam3=NcrnaGffFamilyGrader.extractFamily(commaAttrs, mi3);
		if(!NcrnaGffFamilyGrader.FAM_B.equals(fam3)){
			throw new AssertionError("A comma right after the prefix is a valid boundary and must "
				+"match RF01685 -- observed "+fam3);
		}
	}

	/** Truth TSV must tolerate leading '#' comment lines (and blank lines) before the
	 * header row. */
	private static void testLeadingTruthComments() throws Exception {
		final File f=File.createTempFile("ncrnagraider_truth_comments", ".tsv");
		f.deleteOnExit();
		try(FileWriter w=new FileWriter(f)){
			w.write("# generated by some tool\n");
			w.write("#another comment\n");
			w.write("\n");
			w.write("seqid\tfamily\tstart\tstop\tstrand\n");
			w.write("ctgX\tRF00013\t10\t20\t+\n");
		}
		final ArrayList<Locus> truth=NcrnaGffFamilyGrader.loadTruthTsv(f.getAbsolutePath());
		if(truth.size()!=1){throw new AssertionError("Expected 1 truth locus past leading comments, got "+truth.size());}
	}

	/** An empty truth set (header only, zero data rows) must not crash -- calls simply
	 * become straight FPs with no possible TP/FN. */
	private static void testEmptyTruth(){
		final ArrayList<Locus> truth=new ArrayList<>();
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("ctg1", "RF00013", 100, 200, '+'));
		final GradeResult r=NcrnaGffFamilyGrader.grade(truth, calls, new HashSet<String>());
		check("tpA", 0, r.tpA);
		check("fnA", 0, r.fnA);
		check("fpA", 1, r.fpA);
		check("unionFP", 1, r.unionFP);
	}

	/** An exact-duplicate truth row (identical seqid+family+start+stop+strand twice)
	 * is virtually certainly a data/extraction bug, not real biology -- must fail loud,
	 * not silently double-count. */
	private static void testDuplicateTruthFailsLoud() throws Exception {
		final File f=File.createTempFile("ncrnagraider_truth_dup", ".tsv");
		f.deleteOnExit();
		try(FileWriter w=new FileWriter(f)){
			w.write("seqid\tfamily\tstart\tstop\tstrand\n");
			w.write("ctgX\tRF00013\t10\t20\t+\n");
			w.write("ctgX\tRF00013\t10\t20\t+\n");//exact duplicate
		}
		boolean threw=false;
		try{
			NcrnaGffFamilyGrader.loadTruthTsv(f.getAbsolutePath());
		}catch(IllegalArgumentException expected){threw=true;}
		if(!threw){throw new AssertionError("Exact-duplicate truth row was accepted, not rejected");}
	}

	/** A GFF row with the wrong field count (8, not exactly 9) must fail loud. */
	private static void testMalformedGffFieldCountFailsLoud() throws Exception {
		final File f=File.createTempFile("ncrnagraider_gff_badfields", ".gff");
		f.deleteOnExit();
		try(FileWriter w=new FileWriter(f)){
			w.write("##gff-version 3\n");
			w.write("ctgX\tBBTools\tRNA\t10\t20\t1.0\t+\t0\n");//missing the 9th (attributes) field
		}
		boolean threw=false;
		try{
			NcrnaGffFamilyGrader.loadCallsFromGff(f.getAbsolutePath());
		}catch(IllegalArgumentException expected){threw=true;}
		if(!threw){throw new AssertionError("GFF row with 8 fields (not exactly 9) was accepted, not rejected");}
	}

	/** An RNA-typed GFF row with no "model:" attribute at all is a genuine attribution
	 * gap (Orf.java always sets trnaModel for a generic-ncRNA hit) -- must fail loud. */
	private static void testMissingModelFailsLoud() throws Exception {
		final File f=File.createTempFile("ncrnagraider_gff_nomodel", ".gff");
		f.deleteOnExit();
		try(FileWriter w=new FileWriter(f)){
			w.write("##gff-version 3\n");
			w.write("ctgX\tBBTools\tRNA\t10\t20\t1.0\t+\t0\tstartScr:0,stopScr:0,innerScr:0,len:10\n");
		}
		boolean threw=false;
		try{
			NcrnaGffFamilyGrader.loadCallsFromGff(f.getAbsolutePath());
		}catch(IllegalStateException expected){threw=true;}
		if(!threw){throw new AssertionError("RNA-typed GFF row with no model: attribute was accepted, not rejected");}
	}

	/** Round-2b regression: a single unmatched truth locus geometrically overlapped by
	 * TWO different unmatched opposite-family calls must retain BOTH links in
	 * crossLabelMatchIds, not just the first (or last) one found. Coordinates
	 * hand-verified: T_multi=[100,140] len41. CallB_1=[100,140] len41: overlap=41
	 * (100% both ways). CallB_2=[105,135] len31: overlap=min(140,135)-max(100,105)+1=31
	 * (100% of CallB_2's own 31, 31/41=75.6% of T_multi's 41 -- both >=50%, passes). No
	 * family-B truth exists in this fixture, so both calls are automatically unmatched
	 * at the subtype-B level and eligible for the cross-label check against T_multi. */
	private static void testCrossLabelManyToOneNotOverwritten(){
		final ArrayList<Locus> truth=new ArrayList<>();
		truth.add(new Locus("ctg5", "RF00013", 100, 140, '+'));//T_multi
		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("ctg5", "RF01685", 100, 140, '+'));//CallB_1
		calls.add(new Locus("ctg5", "RF01685", 105, 135, '+'));//CallB_2

		final GradeResult r=NcrnaGffFamilyGrader.grade(truth, calls, new HashSet<String>());
		check("crossAasB (once per truth locus, not per relationship)", 1, r.crossAasB);

		NcrnaGffFamilyGrader.AuditRow truthRow=null;
		for(NcrnaGffFamilyGrader.AuditRow ar : r.auditRows){
			if(ar.role.equals("TRUTH")){truthRow=ar;}
		}
		if(truthRow==null){throw new AssertionError("Missing truth audit row");}
		if(truthRow.crossLabelMatchIds.size()!=2){
			throw new AssertionError("Expected T_multi to retain BOTH cross-matching call IDs, got "
				+truthRow.crossLabelMatchIds.size()+": "+truthRow.crossLabelMatchIds);
		}
		//Both calls must also carry the reciprocal link back to T_multi -- neither
		//overwritten by the other's processing order.
		int callsWithCrossLink=0;
		for(NcrnaGffFamilyGrader.AuditRow ar : r.auditRows){
			if(ar.role.equals("CALL") && ar.crossLabelOutcome.equals("CROSS_MATCH")){
				callsWithCrossLink++;
				if(ar.crossLabelMatchIds.size()!=1 || !ar.crossLabelMatchIds.get(0).equals(truthRow.locusId)){
					throw new AssertionError("Call's cross-link does not point back at the truth locus: "
						+ar.crossLabelMatchIds);
				}
			}
		}
		check("calls carrying a CROSS_MATCH link back to T_multi", 2, callsWithCrossLink);
	}

	/** End-to-end: writes REAL .rawcalls.tsv/.summary.tsv files via writeReport, then
	 * re-reads them from disk (not the in-memory GradeResult) and reconstructs every
	 * subtype/union/neg/dup/cross aggregate purely from the parsed file content,
	 * checked against hand-computed expected values. Exercises TP/FN/FP/DUP/NEG_FP and
	 * a genuine cross-label many-to-one case together in one fixture.
	 *
	 * <p>Fixture (all overlaps hand-verified, not assumed):
	 * <ul>
	 * <li>T0=[ctg1,100,200]+ <-> C0=[ctg1,105,195]+ (overlap91, both len>=50%) -> TP.
	 * <li>T1=[ctg2,100,200]+ : no call anywhere on ctg2 -> FN, no cross (no RF01685 call
	 * on ctg2 either).
	 * <li>C1=[ctg3,100,200]+ : no truth anywhere on ctg3 -> FP (not dup, zero edges).
	 * <li>C2=[ctg1,106,196]+ : overlaps T0 (overlap91) same as C0, loses the augmenting
	 * contest for T0 (C0 got there first and has no alternate edge) -> DUP.
	 * <li>C_neg=[negctg,50,150]+ family RF00013, negctg is negseqids-listed -> NEG_FP.
	 * <li>T_multi=[ctg5,100,140]+, CallB_1=[ctg5,100,140]+, CallB_2=[ctg5,105,135]+ :
	 * same many-to-one cross-label case as testCrossLabelManyToOneNotOverwritten -- both
	 * calls unmatched at subtype-B (FP each), T_multi unmatched at subtype-A (FN), BOTH
	 * cross-linked (crossAasB counts once).
	 * <li>T_B0=[ctg9,500,600]- <-> C_B0=[ctg9,505,595]- (overlap91) -> TP (family B).
	 * </ul>
	 * Hand-computed expected aggregates: tpA=1, fnA=2 (T1,T_multi), fpA=3 (C1,C2,C_neg),
	 * dupA=1 (C2); tpB=1, fnB=0, fpB=2 (CallB_1,CallB_2), dupB=0; negFP=1;
	 * unionTP=3 (T0,T_multi,T_B0 all union-matched -- T_multi via the family-blind
	 * pooled pass finding CallB_1), unionFN=1 (T1 only), unionFP=4 (C1,C2,CallB_2, plus
	 * C_neg counted separately in the neg-handling path); crossAasB=1, crossBasA=0. */
	private static void testEndToEndReconstruction() throws Exception {
		final ArrayList<Locus> truth=new ArrayList<>();
		truth.add(new Locus("ctg1", "RF00013", 100, 200, '+'));//T0
		truth.add(new Locus("ctg2", "RF00013", 100, 200, '+'));//T1
		truth.add(new Locus("ctg5", "RF00013", 100, 140, '+'));//T_multi
		truth.add(new Locus("ctg9", "RF01685", 500, 600, '-'));//T_B0

		final ArrayList<Locus> calls=new ArrayList<>();
		calls.add(new Locus("ctg1", "RF00013", 105, 195, '+'));//C0 -> TP
		calls.add(new Locus("ctg3", "RF00013", 100, 200, '+'));//C1 -> FP
		calls.add(new Locus("ctg1", "RF00013", 106, 196, '+'));//C2 -> DUP
		calls.add(new Locus("negctg", "RF00013", 50, 150, '+'));//C_neg -> NEG_FP
		calls.add(new Locus("ctg5", "RF01685", 100, 140, '+'));//CallB_1
		calls.add(new Locus("ctg5", "RF01685", 105, 135, '+'));//CallB_2
		calls.add(new Locus("ctg9", "RF01685", 505, 595, '-'));//C_B0 -> TP

		final HashSet<String> negSeqids=new HashSet<>();
		negSeqids.add("negctg");

		final GradeResult r=NcrnaGffFamilyGrader.grade(truth, calls, negSeqids);

		final File tmpDir=File.createTempFile("ncrnagraider_e2e", "");
		if(!tmpDir.delete() || !tmpDir.mkdir()){throw new AssertionError("Could not create temp dir "+tmpDir);}
		tmpDir.deleteOnExit();
		final String prefix=tmpDir.getAbsolutePath()+"/out";
		NcrnaGffFamilyGrader.writeReport(prefix, r);

		final File rawFile=new File(prefix+".rawcalls.tsv");
		final File sumFile=new File(prefix+".summary.tsv");
		rawFile.deleteOnExit(); sumFile.deleteOnExit();
		if(!rawFile.exists() || !sumFile.exists()){
			throw new AssertionError("writeReport did not create both expected files");
		}

		//Reconstruct purely from the WRITTEN rawcalls.tsv content -- not from r directly.
		long reconTpA=0, reconFnA=0, reconFpA=0, reconDupA=0;
		long reconTpB=0, reconFnB=0, reconFpB=0, reconDupB=0;
		long reconNegFP=0;
		long reconUnionTP=0, reconUnionFN=0, reconUnionFP=0;
		long reconCrossAasB=0, reconCrossBasA=0;

		final ArrayList<String[]> rows=readTsv(rawFile, 13);
		for(String[] f : rows){
			final String family=f[2], role=f[6], subtypeOutcome=f[7], unionOutcome=f[9], crossOutcome=f[11];
			final boolean isA=family.equals("RF00013");
			if(role.equals("TRUTH")){
				if(subtypeOutcome.equals("TP")){if(isA){reconTpA++;}else{reconTpB++;}}
				else if(subtypeOutcome.equals("FN")){if(isA){reconFnA++;}else{reconFnB++;}}
				if(unionOutcome.equals("TP")){reconUnionTP++;}else if(unionOutcome.equals("FN")){reconUnionFN++;}
				if(crossOutcome.equals("CROSS_DETECTED")){if(isA){reconCrossAasB++;}else{reconCrossBasA++;}}
			}else{//CALL
				if(subtypeOutcome.equals("TP")){/* counted on the TRUTH side already */}
				else if(subtypeOutcome.equals("FP")){if(isA){reconFpA++;}else{reconFpB++;}}
				else if(subtypeOutcome.equals("DUP")){if(isA){reconFpA++; reconDupA++;}else{reconFpB++; reconDupB++;}}
				else if(subtypeOutcome.equals("NEG_FP")){reconNegFP++; if(isA){reconFpA++;}else{reconFpB++;}}
				if(unionOutcome.equals("FP")){reconUnionFP++;}
			}
		}

		check("e2e tpA", 1, reconTpA);
		check("e2e fnA", 2, reconFnA);
		check("e2e fpA", 3, reconFpA);
		check("e2e dupA", 1, reconDupA);
		check("e2e tpB", 1, reconTpB);
		check("e2e fnB", 0, reconFnB);
		check("e2e fpB", 2, reconFpB);
		check("e2e dupB", 0, reconDupB);
		check("e2e negFP", 1, reconNegFP);
		check("e2e unionTP", 3, reconUnionTP);
		check("e2e unionFN", 1, reconUnionFN);
		check("e2e unionFP", 4, reconUnionFP);
		check("e2e crossAasB", 1, reconCrossAasB);
		check("e2e crossBasA", 0, reconCrossBasA);

		//Cross-check against the WRITTEN summary.tsv's own printed numbers too, so both
		//output files are verified, not just rawcalls.
		final HashMap<String, String[]> summaryRows=readSummaryKeyed(sumFile);
		final String[] aRow=summaryRows.get("RF00013"), bRow=summaryRows.get("RF01685"), uRow=summaryRows.get("union");
		if(aRow==null || bRow==null || uRow==null){throw new AssertionError("summary.tsv missing an expected scope row");}
		check("summary.tsv RF00013 TP", 1, Long.parseLong(aRow[0]));
		check("summary.tsv RF00013 FN", 2, Long.parseLong(aRow[1]));
		check("summary.tsv RF00013 FP", 3, Long.parseLong(aRow[2]));
		check("summary.tsv union TP", 3, Long.parseLong(uRow[0]));
		check("summary.tsv union FN", 1, Long.parseLong(uRow[1]));
		check("summary.tsv union FP", 4, Long.parseLong(uRow[2]));
	}

	/** The development family-map path is intentionally separate from the legacy 6S
	 * API.  Exercise three families: one correct R58 call, one correct LSU call, and
	 * one R58-labelled call at an OTHER truth locus.  The latter is a union TP but a
	 * subtype OTHER FN plus R58 FP, proving the generic path keeps both scopes. */
	private static void testGenericFamilyMapThreeFamilyFixture() throws Exception {
		final File map=File.createTempFile("ncrna_familymap", ".tsv");
		final File truth=File.createTempFile("ncrna_truth", ".tsv");
		final File gff=File.createTempFile("ncrna_calls", ".gff");
		map.deleteOnExit(); truth.deleteOnExit(); gff.deleteOnExit();
		try(FileWriter w=new FileWriter(map)){
			w.write("R58\tr58_consensus\nLSU\tlsu_consensus\nOTHER\tother_consensus\n");
		}
		try(FileWriter w=new FileWriter(truth)){
			w.write("seqid\tfamily\tstart\tstop\tstrand\nctg1\tR58\t100\t200\t+\n"
				+"ctg2\tLSU\t300\t500\t-\nctg3\tOTHER\t600\t700\t+\n");
		}
		try(FileWriter w=new FileWriter(gff)){
			w.write("ctg1\tCallGenes\tRNA\t105\t195\t.\t+\t.\tmodel:r58_consensus_0\n"
				+"ctg2\tCallGenes\tRNA\t305\t495\t.\t-\t.\tmodel:lsu_consensus_0\n"
				+"ctg3\tCallGenes\tRNA\t605\t695\t.\t+\t.\tmodel:r58_consensus_1\n");
		}
		final NcrnaGffFamilyGrader.FamilyMap fm=NcrnaGffFamilyGrader.loadFamilyMap(map.getAbsolutePath());
		final NcrnaGffFamilyGrader.GenericGradeResult r=NcrnaGffFamilyGrader.gradeGeneric(
			NcrnaGffFamilyGrader.loadTruthTsv(truth.getAbsolutePath(), fm),
			NcrnaGffFamilyGrader.loadCallsFromGff(gff.getAbsolutePath(), fm), new HashSet<String>(), fm);
		check("generic R58 TP", 1, r.counts.get("R58").tp);
		check("generic R58 FP", 1, r.counts.get("R58").fp);
		check("generic LSU TP", 1, r.counts.get("LSU").tp);
		check("generic OTHER FN", 1, r.counts.get("OTHER").fn);
		check("generic union TP", 3, r.unionTP);
		check("generic union FN", 0, r.unionFN);
		check("generic union FP", 0, r.unionFP);
	}

	/** A familymap row with the wrong field count (not exactly 2 tab-delimited fields)
	 * must fail loud. */
	private static void testFamilymapWrongFieldCountFailsLoud() throws Exception {
		final File f=File.createTempFile("ncrnagraider_familymap_badfields", ".tsv");
		f.deleteOnExit();
		try(FileWriter w=new FileWriter(f)){
			w.write("R58\tr58_consensus\textra\n");//3 fields, not 2
		}
		boolean threw=false;
		try{
			NcrnaGffFamilyGrader.loadFamilyMap(f.getAbsolutePath());
		}catch(IllegalArgumentException expected){threw=true;}
		if(!threw){throw new AssertionError("Familymap row with 3 fields (not exactly 2) was accepted, not rejected");}
	}

	/** A familymap row with an empty family or an empty prefix field must fail loud --
	 * either half is worth checking independently since they're distinct sub-conditions
	 * in the same guard. */
	private static void testFamilymapEmptyFieldFailsLoud() throws Exception {
		final File emptyFamily=File.createTempFile("ncrnagraider_familymap_emptyfamily", ".tsv");
		emptyFamily.deleteOnExit();
		try(FileWriter w=new FileWriter(emptyFamily)){
			w.write("\tr58_consensus\n");
		}
		boolean threwEmptyFamily=false;
		try{
			NcrnaGffFamilyGrader.loadFamilyMap(emptyFamily.getAbsolutePath());
		}catch(IllegalArgumentException expected){threwEmptyFamily=true;}
		if(!threwEmptyFamily){throw new AssertionError("Familymap row with an empty family field was accepted, not rejected");}

		final File emptyPrefix=File.createTempFile("ncrnagraider_familymap_emptyprefix", ".tsv");
		emptyPrefix.deleteOnExit();
		try(FileWriter w=new FileWriter(emptyPrefix)){
			w.write("R58\t\n");
		}
		boolean threwEmptyPrefix=false;
		try{
			NcrnaGffFamilyGrader.loadFamilyMap(emptyPrefix.getAbsolutePath());
		}catch(IllegalArgumentException expected){threwEmptyPrefix=true;}
		if(!threwEmptyPrefix){throw new AssertionError("Familymap row with an empty prefix field was accepted, not rejected");}
	}

	/** Two rows sharing the same family key (even with different prefixes) is a
	 * data/config bug, not a legitimate familymap -- must fail loud, not silently keep
	 * the last-written mapping. */
	private static void testFamilymapDuplicateFamilyFailsLoud() throws Exception {
		final File f=File.createTempFile("ncrnagraider_familymap_dupfamily", ".tsv");
		f.deleteOnExit();
		try(FileWriter w=new FileWriter(f)){
			w.write("R58\tr58_consensus\nR58\tr58_alt_consensus\n");
		}
		boolean threw=false;
		try{
			NcrnaGffFamilyGrader.loadFamilyMap(f.getAbsolutePath());
		}catch(IllegalArgumentException expected){threw=true;}
		if(!threw){throw new AssertionError("Familymap with a duplicate family key was accepted, not rejected");}
	}

	/** Two DIFFERENT family keys sharing the same GFF model-prefix value would make
	 * every GFF row carrying that prefix ambiguous by construction -- must be rejected
	 * at map-load time, not deferred to a per-row ambiguity throw later. */
	private static void testFamilymapDuplicatePrefixFailsLoud() throws Exception {
		final File f=File.createTempFile("ncrnagraider_familymap_dupprefix", ".tsv");
		f.deleteOnExit();
		try(FileWriter w=new FileWriter(f)){
			w.write("R58\tshared_consensus\nLSU\tshared_consensus\n");
		}
		boolean threw=false;
		try{
			NcrnaGffFamilyGrader.loadFamilyMap(f.getAbsolutePath());
		}catch(IllegalArgumentException expected){threw=true;}
		if(!threw){throw new AssertionError("Familymap with two family keys sharing one prefix was accepted, not rejected");}
	}

	/** A familymap with zero valid rows (comments/blank lines only) must fail loud, not
	 * silently produce a FamilyMap that matches nothing. */
	private static void testEmptyFamilymapFailsLoud() throws Exception {
		final File f=File.createTempFile("ncrnagraider_familymap_empty", ".tsv");
		f.deleteOnExit();
		try(FileWriter w=new FileWriter(f)){
			w.write("# no real rows here\n\n");
		}
		boolean threw=false;
		try{
			NcrnaGffFamilyGrader.loadFamilyMap(f.getAbsolutePath());
		}catch(IllegalArgumentException expected){threw=true;}
		if(!threw){throw new AssertionError("Empty familymap (comments/blank lines only) was accepted, not rejected");}
	}

	/** Generic-path prefix ambiguity: two configured prefixes ("r58" and "r58_v2") where
	 * one is itself a token-boundary-respecting prefix of a real model string that also
	 * boundary-matches the other -- "model:r58_v2_x" boundary-matches BOTH "r58" (next
	 * char '_') and "r58_v2" (next char '_'). loadCallsFromGff must refuse to guess and
	 * fail loud rather than silently picking whichever family happened to be inserted
	 * first. */
	private static void testGenericPrefixAmbiguityFailsLoud() throws Exception {
		final File map=File.createTempFile("ncrnagraider_familymap_ambiguous", ".tsv");
		final File gff=File.createTempFile("ncrnagraider_ambiguous_calls", ".gff");
		map.deleteOnExit(); gff.deleteOnExit();
		try(FileWriter w=new FileWriter(map)){
			w.write("A\tr58\nB\tr58_v2\n");
		}
		try(FileWriter w=new FileWriter(gff)){
			w.write("ctg1\tCallGenes\tRNA\t100\t200\t.\t+\t.\tmodel:r58_v2_x\n");
		}
		final NcrnaGffFamilyGrader.FamilyMap fm=NcrnaGffFamilyGrader.loadFamilyMap(map.getAbsolutePath());
		boolean threw=false;
		try{
			NcrnaGffFamilyGrader.loadCallsFromGff(gff.getAbsolutePath(), fm);
		}catch(IllegalStateException expected){threw=true;}
		if(!threw){throw new AssertionError("A model string boundary-matching two configured familymap "
			+"prefixes was accepted, not rejected as ambiguous");}
	}

	/** Direct unit test for the MatchedPair d5/d3 formulas (closes the coverage gap named
	 * in Citan's boundary_metrics_review_20260906.md -- until now these formulas were only
	 * verified by external re-derivation over real run output, never by a test). Every
	 * expected value below is hand-derived from the documented transcript-orientation
	 * convention (NcrnaGffFamilyGrader.MatchedPair javadoc): d5/d3 positive = call end
	 * DOWNSTREAM (toward 3') of the corresponding truth end.
	 * <ul>
	 * <li>Plus strand, call extends both ends: truth [100,200]+, call [95,206] --
	 * call 5' starts 5bp upstream (d5=95-100=-5), call 3' runs 6bp past (d3=206-200=+6).
	 * <li>Plus strand, call truncated both ends: truth [100,200]+, call [103,198] --
	 * d5=+3 (starts late), d3=-2 (ends early).
	 * <li>Minus strand: the transcript runs stop->start, so the 5' end is the STOP
	 * coordinate. Truth [300,500]-, call [305,495]: call 5' at 495 is 5bp downstream of
	 * truth 5' at 500 (d5=500-495=+5, truncation); call 3' at 305 stops 5bp short of
	 * truth 3' at 300 (d3=300-305=-5, truncation).
	 * <li>Minus strand, extension: truth [300,500]-, call [290,510]: d5=500-510=-10,
	 * d3=300-290=+10.
	 * <li>Exact match, both strands: d5=d3=0.
	 * <li>Orientation invariance: a plus-strand pair and its coordinate MIRROR on the
	 * minus strand must produce the IDENTICAL (d5,d3) -- that is the entire point of the
	 * transcript-orientation convention. Mirror of (truth[100,200], call[98,203]) about
	 * the interval is (truth[100,200], call[97,202]) on '-': plus gives (-2,+3); minus
	 * gives (200-202, 100-97)=(-2,+3).
	 * </ul> */
	private static void testBoundaryOffsetFormulas(){
		checkPair("plus/extension", new Locus("c","R58",100,200,'+'), new Locus("c","R58",95,206,'+'), -5, 6);
		checkPair("plus/truncation", new Locus("c","R58",100,200,'+'), new Locus("c","R58",103,198,'+'), 3, -2);
		checkPair("minus/truncation", new Locus("c","R58",300,500,'-'), new Locus("c","R58",305,495,'-'), 5, -5);
		checkPair("minus/extension", new Locus("c","R58",300,500,'-'), new Locus("c","R58",290,510,'-'), -10, 10);
		checkPair("plus/exact", new Locus("c","R58",100,200,'+'), new Locus("c","R58",100,200,'+'), 0, 0);
		checkPair("minus/exact", new Locus("c","R58",300,500,'-'), new Locus("c","R58",300,500,'-'), 0, 0);
		//Orientation invariance: mirrored geometry, identical transcript-space offsets.
		checkPair("plus/mirror-ref", new Locus("c","R58",100,200,'+'), new Locus("c","R58",98,203,'+'), -2, 3);
		checkPair("minus/mirror", new Locus("c","R58",100,200,'-'), new Locus("c","R58",97,202,'-'), -2, 3);
	}

	private static void checkPair(String label, Locus truth, Locus call, int expD5, int expD3){
		final NcrnaGffFamilyGrader.MatchedPair mp=new NcrnaGffFamilyGrader.MatchedPair("R58", truth, call);
		if(mp.d5!=expD5 || mp.d3!=expD3){
			throw new AssertionError(label+": expected (d5,d3)=("+expD5+","+expD3+"), observed ("
				+mp.d5+","+mp.d3+") for truth ["+truth.start+","+truth.stop+"]"+truth.strand
				+" call ["+call.start+","+call.stop+"]"+call.strand);
		}
	}

	/** End-to-end test of writeGenericReport's boundary-accuracy block AND .matches.tsv
	 * (the aggregate half of the same coverage gap): grades a fixture with hand-computed
	 * offsets, then re-reads BOTH written files and checks every boundary number against
	 * hand-derived values -- including a zero-TP family, whose mean/max columns must be
	 * the "." placeholder, never a division by zero or a fabricated 0.
	 * <p>Fixture offsets (formulas verified by testBoundaryOffsetFormulas):
	 * R58: (0,0) exact + (-5,+6) plus-extension + (+5,-5) minus-truncation
	 *   -> nTP=3 exact5=1 exact3=1 both=1 meanAbs=(10/3, 11/3) max=(5,6);
	 * LSU: (-10,+10) minus-extension -> nTP=1, means/max 10;
	 * OTHER: truth with no call -> nTP=0, "." placeholders.
	 * All reciprocal overlaps >=50% of both lengths (hand-checked: 101/112, 91/101,
	 * 201/221 -- each pair's smaller fraction). Expected means are rendered through the
	 * same String.format("%.2f") the grader uses, so the comparison is exact-by-format,
	 * not float-fuzzy. */
	private static void testBoundaryAggregateBlock() throws Exception {
		final File map=File.createTempFile("ncrna_boundary_familymap", ".tsv");
		final File truth=File.createTempFile("ncrna_boundary_truth", ".tsv");
		final File gff=File.createTempFile("ncrna_boundary_calls", ".gff");
		map.deleteOnExit(); truth.deleteOnExit(); gff.deleteOnExit();
		try(FileWriter w=new FileWriter(map)){
			w.write("R58\tr58_consensus\nLSU\tlsu_consensus\nOTHER\tother_consensus\n");
		}
		try(FileWriter w=new FileWriter(truth)){
			w.write("seqid\tfamily\tstart\tstop\tstrand\n"
				+"ctg1\tR58\t100\t200\t+\n"//exact call -> (0,0)
				+"ctg2\tR58\t100\t200\t+\n"//call [95,206] -> (-5,+6)
				+"ctg4\tR58\t1000\t1100\t-\n"//call [1005,1095] -> (+5,-5)
				+"ctg3\tLSU\t300\t500\t-\n"//call [290,510] -> (-10,+10)
				+"ctg5\tOTHER\t600\t700\t+\n");//no call -> FN; OTHER has zero TPs
		}
		try(FileWriter w=new FileWriter(gff)){
			w.write("ctg1\tCallGenes\tRNA\t100\t200\t.\t+\t.\tmodel:r58_consensus_0\n"
				+"ctg2\tCallGenes\tRNA\t95\t206\t.\t+\t.\tmodel:r58_consensus_0\n"
				+"ctg4\tCallGenes\tRNA\t1005\t1095\t.\t-\t.\tmodel:r58_consensus_1\n"
				+"ctg3\tCallGenes\tRNA\t290\t510\t.\t-\t.\tmodel:lsu_consensus_0\n");
		}
		final NcrnaGffFamilyGrader.FamilyMap fm=NcrnaGffFamilyGrader.loadFamilyMap(map.getAbsolutePath());
		final NcrnaGffFamilyGrader.GenericGradeResult r=NcrnaGffFamilyGrader.gradeGeneric(
			NcrnaGffFamilyGrader.loadTruthTsv(truth.getAbsolutePath(), fm),
			NcrnaGffFamilyGrader.loadCallsFromGff(gff.getAbsolutePath(), fm), new HashSet<String>(), fm);
		check("boundary fixture R58 TP", 3, r.counts.get("R58").tp);
		check("boundary fixture LSU TP", 1, r.counts.get("LSU").tp);
		check("boundary fixture OTHER FN", 1, r.counts.get("OTHER").fn);

		final File tmpDir=File.createTempFile("ncrna_boundary_out", "");
		if(!tmpDir.delete() || !tmpDir.mkdir()){throw new AssertionError("Could not create temp dir "+tmpDir);}
		tmpDir.deleteOnExit();
		final String prefix=tmpDir.getAbsolutePath()+"/out";
		NcrnaGffFamilyGrader.writeGenericReport(prefix, r);
		final File sumFile=new File(prefix+".summary.tsv"), mFile=new File(prefix+".matches.tsv");
		sumFile.deleteOnExit(); mFile.deleteOnExit();

		//Written .matches.tsv: exactly the 4 TP pairs, family-map order then truth order,
		//with the hand-computed d5/d3 in the last two columns.
		final String[][] expMatches={
			{"R58","ctg1","+","100","200","100","200","0","0"},
			{"R58","ctg2","+","100","200","95","206","-5","6"},
			{"R58","ctg4","-","1000","1100","1005","1095","5","-5"},
			{"LSU","ctg3","-","300","500","290","510","-10","10"}};
		final ArrayList<String[]> matchRows=readTsv(mFile, 9);
		check("matches.tsv row count", expMatches.length, matchRows.size());
		for(int i=0; i<expMatches.length; i++){
			for(int j=0; j<9; j++){
				if(!expMatches[i][j].equals(matchRows.get(i)[j])){
					throw new AssertionError("matches.tsv row "+i+" col "+j+": expected "
						+expMatches[i][j]+", observed "+matchRows.get(i)[j]);
				}
			}
		}

		//Written summary boundary block: keyed rows past the "# boundary accuracy" marker.
		final HashMap<String,String[]> bRows=new HashMap<>();
		try(BufferedReader br=new BufferedReader(new FileReader(sumFile))){
			boolean inBlock=false; String line;
			while((line=br.readLine())!=null){
				if(line.startsWith("# boundary accuracy")){inBlock=true; continue;}
				if(!inBlock || line.isEmpty() || line.startsWith("family\t")){continue;}
				final String[] f=line.split("\t", -1);
				if(f.length!=9){throw new AssertionError("Boundary block row needs 9 fields: "+line);}
				bRows.put(f[0], f);
			}
		}
		final String m5r58=String.format("%.2f", 10/3.0), m3r58=String.format("%.2f", 11/3.0);
		checkRow("R58 boundary", bRows.get("R58"), "3","1","1","1",m5r58,m3r58,"5","6");
		final String ten=String.format("%.2f", 10.0);
		checkRow("LSU boundary", bRows.get("LSU"), "1","0","0","0",ten,ten,"10","10");
		checkRow("OTHER boundary (zero-TP placeholders)", bRows.get("OTHER"), "0","0","0","0",".",".",".",".");
	}

	private static void checkRow(String label, String[] row, String... expected){
		if(row==null){throw new AssertionError(label+": row missing from summary boundary block");}
		for(int i=0; i<expected.length; i++){
			if(!expected[i].equals(row[i+1])){
				throw new AssertionError(label+" col "+(i+1)+": expected "+expected[i]+", observed "+row[i+1]);
			}
		}
	}

	/** Reads a TSV file, skipping the header row, requiring every data row to have
	 * exactly `expectedFields` columns (fails loud if not -- keeps the test itself
	 * honest about the file format it's parsing). */
	private static ArrayList<String[]> readTsv(File f, int expectedFields) throws Exception {
		final ArrayList<String[]> rows=new ArrayList<>();
		try(BufferedReader r=new BufferedReader(new FileReader(f))){
			String line=r.readLine();//header
			while((line=r.readLine())!=null){
				if(line.isEmpty()){continue;}
				final String[] parts=line.split("\t", -1);
				if(parts.length!=expectedFields){
					throw new AssertionError("Expected "+expectedFields+" fields, got "+parts.length+": "+line);
				}
				rows.add(parts);
			}
		}
		return rows;
	}

	/** summary.tsv-specific reader: keys the three "scope TP FN FP precision recall"
	 * rows by their scope name, ignoring the header/blank/comment/diagnostic lines. */
	private static HashMap<String, String[]> readSummaryKeyed(File f) throws Exception {
		final HashMap<String, String[]> out=new HashMap<>();
		try(BufferedReader r=new BufferedReader(new FileReader(f))){
			String line;
			while((line=r.readLine())!=null){
				if(line.isEmpty() || line.charAt(0)=='#'){continue;}
				final String[] parts=line.split("\t", -1);
				if(parts.length==6 && (parts[0].equals("RF00013") || parts[0].equals("RF01685") || parts[0].equals("union"))){
					out.put(parts[0], new String[]{parts[1], parts[2], parts[3], parts[4], parts[5]});
				}
			}
		}
		return out;
	}

	private static void check(String label, long expected, long observed){
		if(expected!=observed){
			throw new AssertionError("Expected "+label+"="+expected+", observed "+observed);
		}
	}
}
