package prot;

/**
 * Synthetic self-test for {@link MagQC}, the CheckM1-style completeness/
 * contamination oracle. Constructs {@link MarkerVector}s with known copy counts
 * (no file I/O), runs the oracle, and checks the two percentages, the raw counts,
 * and the frozen context fields against hand-computed truth.
 *
 * <h3>Cases (4-family marker set, effective denominator 4)</h3>
 * <ul>
 * <li><b>[1,1,1,1]</b> -&gt; complete 100%, contamination 0% (excess 0).</li>
 * <li><b>[1,1,0,0]</b> -&gt; complete 50% (2/4), contamination 0%.</li>
 * <li><b>[1,2,2,1]</b> -&gt; complete 100%, excess = 0+1+1+0 = 2, contamination
 * 50% (2/4); multi-copy families 2, so both formulas agree here.</li>
 * <li><b>[1,3,1,1]</b> -&gt; complete 100%, excess = 0+2+0+0 = 2, contamination
 * 50% (excess-copy) but multi-copy contamination only 25% (1/4) — showing the
 * headline excess-copy formula counts 3 copies more than 2.</li>
 * <li><b>[] (empty)</b> -&gt; insufficient evidence: NaN %, sufficientEvidence
 * false, ood_status unknown.</li>
 * </ul>
 *
 * <p>Fails loud (throws) on any mismatch; prints each case's report and PASS.</p>
 *
 * <p>Run: {@code java -cp current prot.MagQCTest}</p>
 *
 * @author Eru
 */
public final class MagQCTest {

	/**
	 * Runs the synthetic test.
	 * @param args Ignored.
	 */
	public static void main(String[] args){
		final MagQC qc=new MagQC();

		//Case 1: [1,1,1,1] -> 100% complete, 0% contamination.
		final MagQCResult r1=qc.estimate(vec(new int[]{1,1,1,1}));
		print("[1,1,1,1]", r1);
		checkEq(r1.completeness, 100.0, "case1 completeness");
		checkEq(r1.contamination, 0.0, "case1 contamination");
		check(r1.detectedMarkers==4, "case1 detected != 4: "+r1.detectedMarkers);
		check(r1.excessCopies==0, "case1 excess != 0: "+r1.excessCopies);
		check(r1.effectiveDenominator==4, "case1 denom != 4: "+r1.effectiveDenominator);
		check(r1.sufficientEvidence, "case1 should have evidence");
		check(r1.oodStatus.equals("unknown"), "case1 ood_status != unknown: "+r1.oodStatus);
		check(Double.isNaN(r1.assignmentConfidence), "case1 assignment_confidence not NaN");
		check(r1.domainAssignment.equals("Bacteria"), "case1 domain: "+r1.domainAssignment);

		//Case 2: [1,1,0,0] -> 50% complete, 0% contamination.
		final MagQCResult r2=qc.estimate(vec(new int[]{1,1,0,0}));
		print("[1,1,0,0]", r2);
		checkEq(r2.completeness, 50.0, "case2 completeness");
		checkEq(r2.contamination, 0.0, "case2 contamination");
		check(r2.detectedMarkers==2, "case2 detected != 2: "+r2.detectedMarkers);
		check(r2.excessCopies==0, "case2 excess != 0: "+r2.excessCopies);

		//Case 3: [1,2,2,1] -> 100% complete, excess 2 -> 50% contamination.
		final MagQCResult r3=qc.estimate(vec(new int[]{1,2,2,1}));
		print("[1,2,2,1]", r3);
		checkEq(r3.completeness, 100.0, "case3 completeness");
		check(r3.excessCopies==2, "case3 excess != 2: "+r3.excessCopies);
		checkEq(r3.contamination, 50.0, "case3 contamination (excess-copy)");
		check(r3.multiCopyMarkers==2, "case3 multiCopy != 2: "+r3.multiCopyMarkers);
		checkEq(r3.contaminationMultiCopy, 50.0, "case3 contamination (multi-copy)");

		//Case 4: [1,3,1,1] -> excess 2 -> 50% excess-copy, but 25% multi-copy.
		final MagQCResult r4=qc.estimate(vec(new int[]{1,3,1,1}));
		print("[1,3,1,1]", r4);
		checkEq(r4.completeness, 100.0, "case4 completeness");
		check(r4.excessCopies==2, "case4 excess != 2: "+r4.excessCopies);
		checkEq(r4.contamination, 50.0, "case4 contamination (excess-copy)");
		check(r4.multiCopyMarkers==1, "case4 multiCopy != 1: "+r4.multiCopyMarkers);
		checkEq(r4.contaminationMultiCopy, 25.0, "case4 contamination (multi-copy)");

		//Case 5: empty vector -> insufficient evidence.
		final MagQCResult r5=qc.estimate(vec(new int[]{}));
		print("[] (empty)", r5);
		check(!r5.sufficientEvidence, "case5 should be insufficient evidence");
		check(Double.isNaN(r5.completeness), "case5 completeness should be NaN");
		check(Double.isNaN(r5.contamination), "case5 contamination should be NaN");
		check(r5.effectiveDenominator==0, "case5 denom != 0: "+r5.effectiveDenominator);
		check(r5.oodStatus.equals("unknown"), "case5 ood_status != unknown: "+r5.oodStatus);

		System.out.println("\nPASS: all 5 cases match hand-computed truth; "+
			"all frozen fields populated (raw counts, domain, ood_status=unknown).");
	}

	/**
	 * Builds a marker vector with the given per-family copy counts. Family ids and
	 * representative names are synthetic (fam0..famN-1), domain Bacteria; proteins
	 * matched/unmatched are computed from the counts.
	 * @param counts Per-family copy counts.
	 * @return The synthetic vector.
	 */
	private static MarkerVector vec(final int[] counts){
		final int[] fid=new int[counts.length];
		final String[] rep=new String[counts.length];
		int matched=0;
		for(int i=0; i<counts.length; i++){
			fid[i]=i;
			rep[i]="fam"+i;
			matched+=counts[i];
		}
		return new MarkerVector(counts, fid, rep, "Bacteria", matched, 0);
	}

	/** Prints a case's report block. */
	private static void print(final String label, final MagQCResult r){
		System.out.println("\nvector "+label+":");
		System.out.println("  completeness_pct="+MagQCCLI.fmt(r.completeness)+
			"  contamination_pct="+MagQCCLI.fmt(r.contamination)+
			" (multicopy="+MagQCCLI.fmt(r.contaminationMultiCopy)+")");
		System.out.println("  expected="+r.expectedMarkers+" detected="+r.detectedMarkers+
			" multicopy="+r.multiCopyMarkers+" excess="+r.excessCopies+
			" denom="+r.effectiveDenominator);
		System.out.println("  domain="+r.domainAssignment+" marker_set_id="+r.markerSetId+
			" lineage_taxid="+r.lineageTaxID+" rank="+r.rank);
		System.out.println("  assignment_confidence="+MagQCCLI.fmt(r.assignmentConfidence)+
			" ("+r.assignmentConfidenceModel+") ood_status="+r.oodStatus+
			" sufficient_evidence="+r.sufficientEvidence);
	}

	/** Throws if the condition is false. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new RuntimeException("TEST FAILED: "+message);}
	}

	/** Throws if two doubles differ by more than 1e-9. */
	private static void checkEq(final double actual, final double expected, final String what){
		if(Math.abs(actual-expected)>1e-9){
			throw new RuntimeException("TEST FAILED: "+what+" expected "+expected+
				" got "+actual);
		}
	}
}
