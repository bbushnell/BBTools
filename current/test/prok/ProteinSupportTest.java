package test.prok;

import prok.GeneCaller;
import prok.Orf;
import prok.ProkObject;
import prok.ProteinSupport;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import dna.GeneticCode;

/** Small arithmetic/translation regressions; real-reference calling is a separate pilot.
 * @author Keqing
 */
public final class ProteinSupportTest {

	public static void main(String[] args){
		dpHistograms();
		tailCaps();
		searchPolicy();
		modelBonus();
		pruningProtection();
		System.setProperty("bbtools.diagnostic.minorStartMultiplier", "0.8");
		final byte[] bases=new byte[600]; Arrays.fill(bases, (byte)'A');
		bases[0]='C'; bases[1]='T'; bases[2]='G'; bases[597]='T';
		final Orf orf=new Orf("fixture", 0, 599, 0, 0, bases, false, ProkObject.CDS, false, false);
		orf.startScore=0.7f; orf.stopScore=0.5f; orf.kmerScore=200;
		{
			orf.orfScore=731.25f; // Distinct cached score exercises the no-overlap NN path.
			final float cached=GeneCaller.proteinAwareScore(orf, 0);
			final float overlap=GeneCaller.proteinAwareScore(orf, 12);
			require(overlap==orf.calcOrfScore(12)*0.8f, "CTG overlap baseline must actually exercise minor-start attenuation");
			ProteinSupport.addBonuses(orf, 12.5f, 25f);
			require(orf.stopScore==0.5f && orf.startScore==0.7f && orf.kmerScore==200,
				"Homology contributions must not shift pretrained NN inputs");
			require(GeneCaller.proteinAwareScore(orf, 0)==cached+37.5f, "Cached score must receive the evidence bonus exactly once");
			require(GeneCaller.proteinAwareScore(orf, 12)==overlap+37.5f, "Overlap rescoring must retain the full separate bonus after minor-start scaling");
			boolean duplicate=false;
			try{ProteinSupport.addBonuses(orf, 12.5f, 25f);}catch(IllegalArgumentException expected){duplicate=true;}
			require(duplicate, "A nonzero evidence contribution cannot be applied twice");
			orf.proteinBonus=0; orf.orfScore=cached;
			require(Float.floatToIntBits(GeneCaller.proteinAwareScore(orf, 12))==Float.floatToIntBits(overlap), "Removing support must restore exact overlap-score bits");
			final Orf rna=new Orf("rna", 0, 599, 0, 0, bases, false, ProkObject.RNA);
			rna.orfScore=cached;
			require(GeneCaller.proteinAwareScore(rna, 12)==cached, "Generic RNA must preserve its cached scoring path");
			orf.orfScore=-0.0f;
			ProteinSupport.addBonuses(orf, 0f, 0f);
			require(Float.floatToRawIntBits(orf.orfScore)==Float.floatToRawIntBits(-0.0f), "Zero evidence must preserve signed-zero native score bytes");
		}
		require(ProteinSupport.translationReason("MAKR".getBytes(StandardCharsets.US_ASCII))==null, "Canonical peptide must remain searchable");
		require("ambiguous_codon".equals(ProteinSupport.translationReason(new byte[]{'M','X'})), "Ambiguous codons must not become conserved-gene evidence");
		require("internal_stop".equals(ProteinSupport.translationReason(new byte[]{'M','*','A'})), "Internal stops must not reach the family aligner");
		require("empty_translation".equals(ProteinSupport.translationReason(new byte[0])), "Empty translation must not reach the family aligner");
		final GeneticCode standard=GeneticCode.forTable(11);
		require(Arrays.equals(tail("GTGAAATAA", standard, 500), bytes("VK")), "A tail begins with elongation V, not artificial initiation M");
		require(Arrays.equals(tail("GTGAAATAA", standard, 1), bytes("K")), "Tail cap retains the residue immediately before the stop");
		require(Arrays.equals(tail("ATGTGATAA", GeneticCode.forTable(4), 500), bytes("MW")), "Table4 TGA sense residue must survive internal translation");
		require(Arrays.equals(tail("ATGTGATAA", GeneticCode.forTable(25), 500), bytes("MG")), "Table25 TGA sense residue must survive internal translation");
		require("internal_stop".equals(ProteinSupport.translationReason(tail("ATGTGATAA", standard, 500))), "A wrong-code internal stop cannot become tail evidence");
		final byte[] longOrf=new byte[1803]; Arrays.fill(longOrf, (byte)'A');
		longOrf[1800]='T';
		final byte[] capped=ProteinSupport.translateTail(longOrf, 0, 1802, standard, 500);
		require(capped.length==500 && capped[0]=='K' && capped[499]=='K', "Long ORF translation must be bounded at five hundred sense codons");
		boolean badStop=false;
		try{tail("GTGAAAAAA", standard, 500);}catch(IllegalArgumentException expected){badStop=true;}
		require(badStop, "A synthetic or non-stop endpoint cannot be translated as a real stop");
		System.err.println("PROTEIN_SUPPORT_SCORE_PASS cached_once=true overlap_once=true zero_bits=true native_nn_inputs=true rna_unchanged=true tail_translation=true");
	}
	private static byte[] bytes(String s){return s.getBytes(StandardCharsets.US_ASCII);}
	/** Literal histogram oracle, independent of aligner output or summary calculations. */
	private static void dpHistograms(){
		require(ProteinSupport.DP_HISTOGRAM_BINS==101, "Top50 scoring plus at most50 trace reruns needs bins0..100");
		final long[] matched=new long[101], unmatched=new long[101];
		ProteinSupport.recordDpCounts(matched, unmatched, 0, false);
		ProteinSupport.recordDpCounts(matched, unmatched, 4, false);
		ProteinSupport.recordDpCounts(matched, unmatched, 100, true);
		ProteinSupport.recordDpCounts(matched, unmatched, 50, true);
		final long[] expectedMatched=new long[101], expectedUnmatched=new long[101];
		expectedMatched[50]=1; expectedMatched[100]=1; expectedUnmatched[0]=1; expectedUnmatched[4]=1;
		require(Arrays.equals(matched, expectedMatched) && Arrays.equals(unmatched, expectedUnmatched),
			"Each searched query increments exactly one final-assignment bin, including zero-DP no-hits");
		final long[] left=new long[101], right=new long[101], merged=new long[101];
		left[0]=2; left[5]=5; left[50]=1; left[100]=1;
		right[0]=2; right[1]=6; right[5]=3;
		final long[] savedLeft=left.clone(), savedRight=right.clone();
		ProteinSupport.mergeDpCounts(merged, left); ProteinSupport.mergeDpCounts(merged, right);
		final long[] expected=new long[101]; expected[0]=4; expected[1]=6; expected[5]=8; expected[50]=1; expected[100]=1;
		require(Arrays.equals(merged, expected) && Arrays.equals(left, savedLeft) && Arrays.equals(right, savedRight),
			"Merge adds literal bins without altering donor arrays");
		require(ProteinSupport.dpSummary(merged).equals("queries=20 dp_fills=196 zero_dp=4 mean=9.8 median=3.0 p95=50.0 max=100 aligned_queries=16 aligned_mean=12.25"),
			"Twenty observations have196 fills, middle pair1/5 and nearest-rank95th observation50");
		final long[] odd=new long[101]; odd[0]=1; odd[4]=1; odd[100]=1;
		require(ProteinSupport.dpSummary(odd).equals("queries=3 dp_fills=104 zero_dp=1 mean=34.666666666666664 median=4.0 p95=100.0 max=100 aligned_queries=2 aligned_mean=52.0"),
			"Odd median and95th percentile must come from literal ranks, not interpolated bins");
		require(ProteinSupport.dpSummary(new long[101]).equals("queries=0 dp_fills=0 zero_dp=0 mean=NaN median=NaN p95=NaN max=0 aligned_queries=0 aligned_mean=NaN"),
			"Empty cohorts must not report invented zero means or percentiles");
		final long[] zeros=new long[101]; zeros[0]=5;
		require(ProteinSupport.dpSummary(zeros).equals("queries=5 dp_fills=0 zero_dp=5 mean=0.0 median=0.0 p95=0.0 max=0 aligned_queries=0 aligned_mean=NaN"),
			"All-zero work is a real cohort but has no aligned-query mean");
		zeros[0]=Long.MAX_VALUE;
		require(ProteinSupport.dpSummary(zeros).contains("median=0.0 p95=0.0 max=0"),
			"Nearest ranks must not overflow when multiplied percentages would");
		for(final long[] bad : new long[][]{null, new long[100], new long[102]}){
			rejectDp(()->ProteinSupport.dpSummary(bad));
			rejectDp(()->ProteinSupport.mergeDpCounts(merged, bad));
			rejectDp(()->ProteinSupport.mergeDpCounts(bad, merged));
			rejectDp(()->ProteinSupport.recordDpCounts(bad, unmatched, 1, true));
			rejectDp(()->ProteinSupport.recordDpCounts(matched, bad, 1, false));
		}
		rejectDp(()->ProteinSupport.recordDpCounts(matched, unmatched, -1, false));
		rejectDp(()->ProteinSupport.recordDpCounts(matched, unmatched, 101, true));
		rejectDp(()->ProteinSupport.recordDpCounts(matched, matched, 1, true));
		rejectDp(()->ProteinSupport.mergeDpCounts(merged, merged));
		final long[] bad=new long[101]; bad[100]=-1;
		rejectDp(()->ProteinSupport.dpSummary(bad));
		rejectDp(()->ProteinSupport.mergeDpCounts(merged, bad));
		rejectDp(()->ProteinSupport.mergeDpCounts(bad, merged));
		rejectDp(()->ProteinSupport.recordDpCounts(bad, unmatched, 100, true));
		final long[] overflowing=new long[101]; overflowing[100]=Long.MAX_VALUE;
		rejectDp(()->ProteinSupport.mergeDpCounts(merged, overflowing));
		require(Arrays.equals(merged, expected), "A late invalid/overflowed bin must leave the entire merge destination untouched");
		rejectDp(()->ProteinSupport.recordDpCounts(overflowing, unmatched, 100, true));
		rejectDp(()->ProteinSupport.dpSummary(overflowing)); // Weighted DP total overflow.
		final long[] sumOverflow=new long[101]; sumOverflow[50]=Long.MAX_VALUE/50; sumOverflow[100]=1;
		rejectDp(()->ProteinSupport.dpSummary(sumOverflow)); // Individual products fit, but their sum does not.
		zeros[1]=1;
		rejectDp(()->ProteinSupport.dpSummary(zeros)); // Query total overflow.
		require(Arrays.equals(matched, expectedMatched) && Arrays.equals(unmatched, expectedUnmatched),
			"Rejected observations must not mutate earlier matched/unmatched counts");
		System.err.println("PROTEIN_DP_HISTOGRAM_TEST_PASS literal_bins=true zero_nohit=true merge=true ranks=true overflow=true");
	}
	private static void rejectDp(Runnable action){
		boolean rejected=false;
		try{action.run();}catch(IllegalArgumentException | ArithmeticException expected){rejected=true;}
		require(rejected, "Invalid or overflowing histogram arithmetic must fail loudly");
	}
	/** Nondefault windows retain the stop-facing codons and declare their GFF provenance. */
	private static void tailCaps(){
		final byte[] bases=new byte[1806]; Arrays.fill(bases, (byte)'A');
		bases[1800]='G'; bases[1801]='T'; bases[1802]='G'; bases[1803]='T';
		for(int cap : new int[]{200, 300, 400, 500, 600}){
			final byte[] peptide=ProteinSupport.translateTail(bases, 0, 1805, GeneticCode.forTable(11), cap);
			require(peptide.length==cap && peptide[0]=='K' && peptide[cap-1]=='V',
				"Every requested cap retains the last sense codon but omits the real stop");
			final structures.ByteBuilder row=new structures.ByteBuilder(); row.append("existing");
			ProteinSupport.appendTailCap(row, cap);
			require(row.toString().equals(cap==500 ? "existing" : "existing;bbgap_tail_cap="+cap),
				"Default500 must preserve bytes; every other cap must be explicit without replacing existing attributes");
		}
		System.err.println("PROTEIN_TAIL_CAP_TEST_PASS translation_200_600=true default_gff_bytes=true nondefault_provenance=true");
	}
	/** CLI auto settings must leave the exhaustive comparator unfiltered and reject ambiguous contracts. */
	private static void searchPolicy(){
		final prot.ProteinTailAssigner.SearchPolicy fast=ProteinSupport.searchPolicy("FAST", -1, 0, -1);
		final prot.ProteinTailAssigner.SearchPolicy exhaustive=ProteinSupport.searchPolicy("exhaustive", -1, 0, -1);
		require(fast.mode==prot.ProteinTailAssigner.Mode.FAST && fast.minSeeds==1 && fast.lookahead==4,
			"Fast CLI auto selects the explicit experimental seed floor and bounded-four search");
		require(exhaustive.mode==prot.ProteinTailAssigner.Mode.EXHAUSTIVE && exhaustive.minSeeds==0 && exhaustive.lookahead==0,
			"Exhaustive CLI auto must reproduce original candidate extent without seed or stopping filters");
		final prot.ProteinTailAssigner.SearchPolicy first=ProteinSupport.searchPolicy("fast", 0, 0, 0);
		require(first.minSeeds==0 && first.lookahead==0, "Explicit zero selects a no-seed, first-accepted-hit control");
		for(String mode : new String[]{null, "unknown"}){
			boolean failed=false;
			try{ProteinSupport.searchPolicy(mode, -1, 0, -1);}catch(IllegalArgumentException expected){failed=true;}
			require(failed, "Unknown search modes must fail before resource loading or output creation");
		}
		for(double density : new double[]{-0.01, 1.01, Double.NaN, Double.POSITIVE_INFINITY}){
			boolean failed=false;
			try{ProteinSupport.searchPolicy("fast", 1, density, 4);}catch(IllegalArgumentException expected){failed=true;}
			require(failed, "Invalid seed density cannot silently discard or admit every family");
		}
		boolean filteredExhaustive=false;
		try{ProteinSupport.searchPolicy("exhaustive", 1, 0, 0);}catch(IllegalArgumentException expected){filteredExhaustive=true;}
		require(filteredExhaustive, "Filtered searches must not be mislabeled as the frozen exhaustive comparator");
		System.err.println("PROTEIN_SEARCH_POLICY_TEST_PASS explicit_auto=true legacy_unfiltered=true invalid=true");
	}
	/** Literal orthogonal arms: neither, reward plus retention, retention only, reward only. */
	private static void pruningProtection(){
		final ProteinSupport.Protection auto=ProteinSupport.Protection.forName("auto");
		final ProteinSupport.Protection on=ProteinSupport.Protection.forName("ON"), off=ProteinSupport.Protection.forName("off");
		require(!auto.protects(0, true) && auto.protects(50, true), "Default zero/positive behavior must remain native/current respectively");
		require(on.protects(0, true) && !off.protects(50, true), "Retention-only and reward-only arms must differ solely in the pruning exemption");
		final ProteinSupport.StopEvidence noHit=new ProteinSupport.StopEvidence(null, null, 11, 50);
		final ProteinSupport.StopEvidence unsupported=new ProteinSupport.StopEvidence(null, "synthetic_stop", 11, 0);
		require(!on.protects(0, noHit.supported()) && !on.protects(0, unsupported.supported()),
			"Annotation metadata alone must never protect a no-hit or unsupported stop");
		final byte[] bases=bytes("ATGAAATAA");
		final Orf candidate=new Orf("protection", 0, 8, 0, 0, bases, false, ProkObject.CDS, false, false);
		candidate.orfScore=-0.0f;
		ProteinSupport.addBonuses(candidate, 0, 0);
		require(!ProteinSupport.protectsFromPruning(candidate), "Unannotated zero candidate retains the original pruning path");
		require(on.protects(candidate.proteinBonus, true)
			&& Float.floatToRawIntBits(candidate.orfScore)==Float.floatToRawIntBits(-0.0f),
			"Retention-only protection must not alter cached native score bits");
		ProteinSupport.addBonuses(candidate, 25, 25);
		require(ProteinSupport.protectsFromPruning(candidate) && candidate.orfScore==50,
			"Legacy positive-bonus callers without metadata preserve the old auto predicate");
		require(!off.protects(candidate.proteinBonus, true) && candidate.orfScore==50,
			"Turning protection off must not remove the independent additive reward");
		boolean rejected=false;
		try{ProteinSupport.Protection.forName("sometimes");}catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("hbmprotect");}
		require(rejected, "Unknown protection policy must fail before output creation");
		System.err.println("PROTEIN_PRUNING_PROTECTION_TEST_PASS four_arms=true zero_bits=true nohit_unsupported=true invalid=true");
	}
	/** Literal score oracle for the proposed 20+30*L/mean formula and unchanged flat/zero modes. */
	private static void modelBonus(){
		require(!ProteinSupport.modelWeighted("flat") && ProteinSupport.modelWeighted("MODEL"), "Bonus modes must retain flat default semantics and case-insensitive CLI values");
		final int[] lengths={150, 300, 600}; final float[] totals={35, 50, 80};
		for(int i=0; i<lengths.length; i++){
			require(20+ProteinSupport.middleBonus(30, lengths[i], 300, true)==totals[i], "Model-length reward differs from literal 35/50/80 oracle");
			require(Float.floatToIntBits(ProteinSupport.middleBonus(25, lengths[i], 300, false))==Float.floatToIntBits(25),
				"Flat bonus must preserve its original bits for every model length");
			require(ProteinSupport.middleBonus(0, lengths[i], 300, true)==0, "Model weighting cannot turn zero support into an active rescue");
		}
		boolean badMode=false;
		try{ProteinSupport.modelWeighted("length");}catch(IllegalArgumentException expected){badMode=expected.getMessage().contains("flat or model");}
		require(badMode, "Unknown bonus mode must fail before loading resources or opening output");
		for(double mean : new double[]{0, -1, Double.NaN, Double.POSITIVE_INFINITY}){
			boolean failed=false;
			try{ProteinSupport.middleBonus(30, 300, mean, true);}catch(IllegalArgumentException expected){failed=true;}
			require(failed, "Undefined reference mean must not become a numeric protein bonus");
		}
		boolean overflow=false;
		try{ProteinSupport.middleBonus(Float.MAX_VALUE, 600, 300, true);}catch(IllegalArgumentException expected){overflow=true;}
		require(overflow, "Unrepresentable model reward must fail rather than enter DP as infinity");
		System.err.println("PROTEIN_MODEL_BONUS_TEST_PASS literal_35_50_80=true flat_bits=true zero=true invalid=true");
	}
	private static byte[] tail(String s, GeneticCode code, int cap){final byte[] b=bytes(s); return ProteinSupport.translateTail(b, 0, b.length-1, code, cap);}
	private static void require(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private ProteinSupportTest(){}
}
