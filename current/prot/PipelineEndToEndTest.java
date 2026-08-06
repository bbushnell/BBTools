package prot;

import java.util.ArrayList;
import java.util.List;

/**
 * End-to-end INTEGRATION self-test that composes the ENTIRE basic MAG-QC pipeline
 * in memory, with no disk I/O, to prove the five pieces actually compose — not just
 * that each works alone. The output object of every stage is fed as the input of
 * the next, so a type/format mismatch between any two stages fails this test.
 *
 * <p>Chain exercised (each arrow is a real hand-off):</p>
 * <ol>
 * <li>{@link GenomeProteins} (synthetic Bacteria genomes) -&gt;
 * <li>{@link MarkerFactory#build} -&gt; per-domain {@link MarkerSet} -&gt;
 * <li>{@link MarkerVectorizer#vectorize} (bin + marker set) -&gt;
 * {@link MarkerVector} -&gt;
 * <li>{@link MagQC#estimate} -&gt; {@link MagQCResult}.
 * </ol>
 *
 * <h3>Ground truth built into the fixtures</h3>
 * <p>Six Bacteria genomes G1..G6, each carrying exactly one copy of five
 * mutually-dissimilar proteins A,B,C,D,E. Every family is therefore universal
 * single-copy (present exactly once in 6/6 genomes, fraction 1.0 &ge; the 0.97
 * selection threshold), so the factory MUST select all 5 -&gt; marker-set
 * dimension N = 5.</p>
 *
 * <p>The synthetic bin is built with a KNOWN marker profile against those 5
 * markers:</p>
 * <ul>
 * <li>A -&gt; 1 copy (present exactly once)
 * <li>B -&gt; 2 copies (DUPLICATED — the one contamination event)
 * <li>C -&gt; 1 copy (present exactly once)
 * <li>D -&gt; 1 copy (present exactly once)
 * <li>E -&gt; 0 copies (MISSING — the one incompleteness event)
 * <li>F -&gt; 1 copy of a sixth, dissimilar NON-marker protein (matches nothing)
 * </ul>
 *
 * <h3>Hand-computed expected end result (order-independent)</h3>
 * <p>Vector count multiset = {one 2, three 1s, one 0}; dimension N = 5.</p>
 * <ul>
 * <li>familiesPresent = 4 (A,B,C,D; E absent)
 * <li>familiesExactlyOnce = 3 (A,C,D)
 * <li>familiesMultiCopy = 1 (B)
 * <li>excessCopies = (2-1) = 1
 * <li>proteinsMatched = 5 (A,B,B,C,D), proteinsUnmatched = 1 (F)
 * <li><b>completeness = 100 * detected/denom = 100 * 4/5 = 80.0 %</b>
 *   = (N-1)/N * 100, the incompleteness from the one missing family.
 * <li><b>contamination (headline, excess-copy) = 100 * 1/5 = 20.0 %</b>
 *   from the one duplicated family.
 * <li>contaminationMultiCopy (secondary) = 100 * 1/5 = 20.0 %
 * </ul>
 *
 * <p>Assertions use {@code assert}; run with {@code -ea}. Prints the marker set,
 * the vector, the raw counts, and the final percentages, then PASS.</p>
 *
 * <p>Run: {@code java -ea --add-modules jdk.incubator.vector -cp current
 * prot.PipelineEndToEndTest}</p>
 *
 * @author Eru
 */
public final class PipelineEndToEndTest {

	//Five mutually-dissimilar ~60-aa base sequences (shared design with the
	//per-stage tests): A..E are the universal single-copy marker families.
	private static final String A=
		"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVK";
	private static final String B=
		"MSDNKGLAEFTQGVMRHGDVLLKQGIVDESLYNRMLDEACGKHFPWLLSQFDNGYWETIN";
	private static final String C=
		"MTEQKLISEEDLNSAVDHHFLKNMGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGST";
	private static final String D=
		"MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITK";
	private static final String E=
		"MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIW";
	//A sixth, dissimilar NON-marker protein present only in the bin (matches no
	//family): distinct residue composition from A..E, standard residues only.
	private static final String F=
		"WWFFYYCCPPHHGGMMNNQQKKRRDDEESSTTAAIILLVVWWFFYYCCPPHHGGMMNNQQ";

	/**
	 * Runs the end-to-end integration test.
	 * @param args Ignored.
	 */
	public static void main(String[] args){
		System.out.println("=== MAG-QC pipeline end-to-end integration test ===");

		//---- Stage 1: synthetic genomes (GenomeProteins) ----
		final String[] markerSeqs={A, B, C, D, E};
		final String[] markerNames={"markerA", "markerB", "markerC", "markerD", "markerE"};
		final ArrayList<GenomeProteins> genomes=new ArrayList<GenomeProteins>();
		for(int gi=1; gi<=6; gi++){
			final String gid="G"+gi;
			final ArrayList<ProteinSequence> prots=new ArrayList<ProteinSequence>();
			for(int m=0; m<markerSeqs.length; m++){
				//One copy of each universal single-copy family in every genome.
				prots.add(new ProteinSequence(gid+"_"+markerNames[m], markerSeqs[m]));
			}
			genomes.add(new GenomeProteins(gid, "Bacteria", "synthetic", prots));
		}
		System.out.println("Stage 1: built "+genomes.size()+" Bacteria genomes, "+
			markerSeqs.length+" proteins each.");

		//---- Stage 2: MarkerFactory.build -> MarkerSet ----
		final MarkerFactory factory=new MarkerFactory();
		final List<MarkerSet> sets=factory.build(genomes, "e2e-v1",
			"2026-07-31T00:00:00Z", "NA");
		assert(sets.size()==1) : "expected exactly 1 domain marker set, got "+sets.size();
		final MarkerSet ms=sets.get(0);
		System.out.println("Stage 2: marker set domain="+ms.domain+" version="+ms.version+
			" genomes="+ms.genomeCount+" families="+ms.families.size()+
			" selected="+ms.selectedCount());
		System.out.println("  "+MarkerFactoryCLI.header());
		for(final MarkerFamily f : ms.families){
			System.out.println("  "+MarkerFactoryCLI.row(ms, f));
		}
		//The factory must have selected all 5 universal single-copy families.
		final int N=ms.selectedCount();
		assert(ms.domain.equals("Bacteria")) : "wrong domain: "+ms.domain;
		assert(ms.genomeCount==6) : "expected 6 genomes, got "+ms.genomeCount;
		assert(N==5) : "expected 5 selected markers (N), got "+N;

		//---- Stage 3: build a synthetic bin with a KNOWN marker profile ----
		final ArrayList<ProteinSequence> bin=new ArrayList<ProteinSequence>();
		bin.add(new ProteinSequence("bin_A_1", A));      //family A: 1 copy
		bin.add(new ProteinSequence("bin_B_1", B));      //family B: 2 copies (duplicate)
		bin.add(new ProteinSequence("bin_B_2", B));
		bin.add(new ProteinSequence("bin_C_1", C));      //family C: 1 copy
		bin.add(new ProteinSequence("bin_D_1", D));      //family D: 1 copy
		//family E intentionally absent from the bin (incompleteness)
		bin.add(new ProteinSequence("bin_nonMarker_F", F)); //matches nothing
		System.out.println("Stage 3: bin has "+bin.size()+" proteins "+
			"(A x1, B x2, C x1, D x1, E x0, non-marker F x1).");

		//---- Stage 4: MarkerVectorizer.vectorize -> MarkerVector ----
		final MarkerVectorizer vec=new MarkerVectorizer();
		final MarkerVector mv=vec.vectorize(bin, ms);
		System.out.println("Stage 4: MarkerVector dimension="+mv.dimension());
		System.out.println("  index\tfamily_id\trepresentative\tcount");
		for(int i=0; i<mv.dimension(); i++){
			System.out.println("  "+i+"\t"+mv.familyIds[i]+"\t"+
				mv.representativeIds[i]+"\t"+mv.counts[i]);
		}
		System.out.println("  familiesPresent="+mv.familiesPresent()+
			" familiesExactlyOnce="+mv.familiesExactlyOnce()+
			" familiesMultiCopy="+mv.familiesMultiCopy());
		System.out.println("  proteinsMatched="+mv.proteinsMatched+
			" proteinsUnmatched="+mv.proteinsUnmatched);

		//Vector-level ground truth (order-independent: assert on the count multiset).
		assert(mv.dimension()==N) : "vector dimension "+mv.dimension()+" != N "+N;
		int ones=0, twos=0, zeros=0, sum=0;
		for(final int c : mv.counts){
			if(c==0){zeros++;}else if(c==1){ones++;}else if(c==2){twos++;}
			else{throw new RuntimeException("unexpected count "+c);}
			sum+=c;
		}
		assert(twos==1) : "expected exactly one duplicated family (count 2), got "+twos;
		assert(ones==3) : "expected exactly three single-copy families (count 1), got "+ones;
		assert(zeros==1) : "expected exactly one missing family (count 0), got "+zeros;
		assert(sum==5) : "sum of copy counts != 5: "+sum;
		assert(mv.familiesPresent()==4) : "familiesPresent != 4: "+mv.familiesPresent();
		assert(mv.familiesExactlyOnce()==3) :
			"familiesExactlyOnce != 3: "+mv.familiesExactlyOnce();
		assert(mv.familiesMultiCopy()==1) : "familiesMultiCopy != 1: "+mv.familiesMultiCopy();
		assert(mv.proteinsMatched==5) : "proteinsMatched != 5: "+mv.proteinsMatched;
		assert(mv.proteinsUnmatched==1) :
			"proteinsUnmatched != 1 (non-marker F leaked into a family): "+mv.proteinsUnmatched;
		assert(mv.domain.equals("Bacteria")) : "vector domain: "+mv.domain;

		//---- Stage 5: MagQC.estimate -> MagQCResult ----
		final MagQC qc=new MagQC();
		final MagQCResult r=qc.estimate(mv, ms);
		System.out.println("Stage 5: MagQC result");
		System.out.println("  completeness_pct   = "+MagQCCLI.fmt(r.completeness));
		System.out.println("  contamination_pct  = "+MagQCCLI.fmt(r.contamination)+
			" (multicopy="+MagQCCLI.fmt(r.contaminationMultiCopy)+")");
		System.out.println("  expected="+r.expectedMarkers+" detected="+r.detectedMarkers+
			" multicopy="+r.multiCopyMarkers+" excess="+r.excessCopies+
			" denom="+r.effectiveDenominator);
		System.out.println("  domain="+r.domainAssignment+" marker_set_id="+r.markerSetId+
			" ood_status="+r.oodStatus+" sufficient_evidence="+r.sufficientEvidence);

		//Final ground truth: N=5, N-1=4 present with one duplicated.
		// completeness = (N-1)/N * 100 = 4/5 * 100 = 80.0
		// contamination = excess/N * 100 = 1/5 * 100 = 20.0
		assert(r.effectiveDenominator==5) : "denom != 5: "+r.effectiveDenominator;
		assert(r.detectedMarkers==4) : "detected != 4: "+r.detectedMarkers;
		assert(r.excessCopies==1) : "excess != 1: "+r.excessCopies;
		assert(r.multiCopyMarkers==1) : "multiCopy != 1: "+r.multiCopyMarkers;
		assert(near(r.completeness, 80.0)) : "completeness != 80.0: "+r.completeness;
		assert(near(r.contamination, 20.0)) : "contamination != 20.0: "+r.contamination;
		assert(near(r.contaminationMultiCopy, 20.0)) :
			"contaminationMultiCopy != 20.0: "+r.contaminationMultiCopy;
		assert(r.sufficientEvidence) : "should have sufficient evidence";
		assert(r.domainAssignment.equals("Bacteria")) : "result domain: "+r.domainAssignment;
		assert(r.markerSetId.equals("e2e-v1")) : "marker_set_id: "+r.markerSetId;

		System.out.println("\nPASS: the 5 pipeline pieces compose end-to-end. "+
			"Hand-computed truth confirmed: completeness=80.0% (4/5 families, one missing) "+
			"and contamination=20.0% (1 excess copy, one duplicated) out of N=5 markers.");
	}

	/** True if two doubles are within 1e-9. */
	private static boolean near(final double a, final double b){return Math.abs(a-b)<1e-9;}
}
