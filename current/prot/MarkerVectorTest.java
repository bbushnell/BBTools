package prot;

import java.util.ArrayList;
import java.util.List;

/**
 * Synthetic self-test for {@link MarkerVectorizer}. Constructs a marker set with
 * four known single-copy marker families directly (no file I/O), then builds a
 * synthetic bin whose proteins have a known ground-truth copy count per family,
 * and verifies the produced {@link MarkerVector} and its derived scalars.
 *
 * <p>Ground truth — four selected families A, B, C, D (family ids 0..3); the bin
 * contains:</p>
 * <ul>
 * <li><b>A</b> — one copy -&gt; count 1 (present exactly once).</li>
 * <li><b>B</b> — two copies -&gt; count 2 (multi-copy).</li>
 * <li><b>C</b> — no copy -&gt; count 0 (absent).</li>
 * <li><b>D</b> — one copy -&gt; count 1 (present exactly once).</li>
 * <li><b>nonMarker</b> — a fifth, dissimilar protein that matches no family.</li>
 * </ul>
 *
 * <p>Expected vector [1,2,0,1]; present=3, exactly-once=2, multi-copy=1;
 * matched=4, unmatched=1. Fails loud (throws) on any mismatch; prints the vector
 * and PASS on success.</p>
 *
 * <p>Run: {@code java -cp current prot.MarkerVectorTest}</p>
 *
 * @author Eru
 */
public final class MarkerVectorTest {

	//Five mutually-dissimilar ~60-aa base sequences (shared with MarkerFactoryTest's
	//design): A..D are the marker family representatives; E is a non-marker protein.
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

	/**
	 * Runs the synthetic test.
	 * @param args Ignored.
	 */
	public static void main(String[] args){
		//Build a marker set with 4 selected single-copy families directly.
		final MarkerSet ms=buildMarkerSet();
		System.out.println("Marker set: domain="+ms.domain+" families="+ms.families.size()+
			" selected="+ms.selectedCount());

		//Build a synthetic bin with the known copy counts.
		final ArrayList<ProteinSequence> bin=new ArrayList<ProteinSequence>();
		bin.add(new ProteinSequence("bin_A_1", A));    //family A: count 1
		bin.add(new ProteinSequence("bin_B_1", B));    //family B: count 2
		bin.add(new ProteinSequence("bin_B_2", B));
		bin.add(new ProteinSequence("bin_D_1", D));    //family D: count 1
		bin.add(new ProteinSequence("bin_nonMarker", E)); //matches nothing
		//(family C intentionally absent from the bin)

		final MarkerVectorizer vec=new MarkerVectorizer();
		final MarkerVector mv=vec.vectorize(bin, ms);

		//Print the produced vector and scalars.
		System.out.println("\nBin proteins: "+bin.size());
		System.out.println("index\tfamily_id\trepresentative\tcount");
		for(int i=0; i<mv.dimension(); i++){
			System.out.println(i+"\t"+mv.familyIds[i]+"\t"+mv.representativeIds[i]+"\t"+mv.counts[i]);
		}
		System.out.println("dimension="+mv.dimension());
		System.out.println("familiesPresent="+mv.familiesPresent()+
			" familiesExactlyOnce="+mv.familiesExactlyOnce()+
			" familiesMultiCopy="+mv.familiesMultiCopy());
		System.out.println("proteinsMatched="+mv.proteinsMatched+
			" proteinsUnmatched="+mv.proteinsUnmatched);

		//Verify the vector against ground truth.
		check(mv.dimension()==4, "expected dimension 4, got "+mv.dimension());
		check(mv.counts[0]==1, "family A (0) count != 1: "+mv.counts[0]);
		check(mv.counts[1]==2, "family B (1) count != 2: "+mv.counts[1]);
		check(mv.counts[2]==0, "family C (2) count != 0: "+mv.counts[2]);
		check(mv.counts[3]==1, "family D (3) count != 1: "+mv.counts[3]);
		check(mv.familiesPresent()==3, "familiesPresent != 3: "+mv.familiesPresent());
		check(mv.familiesExactlyOnce()==2, "familiesExactlyOnce != 2: "+mv.familiesExactlyOnce());
		check(mv.familiesMultiCopy()==1, "familiesMultiCopy != 1: "+mv.familiesMultiCopy());
		check(mv.proteinsMatched==4, "proteinsMatched != 4: "+mv.proteinsMatched);
		check(mv.proteinsUnmatched==1, "proteinsUnmatched != 1 (non-marker leaked): "+
			mv.proteinsUnmatched);
		check(mv.domain.equals("Bacteria"), "wrong domain: "+mv.domain);

		System.out.println("\nPASS: vector=[1,2,0,1]; present=3, exactly-once=2, multi-copy=1; "+
			"non-marker protein matched nothing (matched=4, unmatched=1).");
	}

	/**
	 * Builds a marker set with four selected single-copy families A,B,C,D whose
	 * representatives are the base sequences, provenance thresholds 90% / 0.8.
	 * @return The synthetic marker set.
	 */
	private static MarkerSet buildMarkerSet(){
		final String[] seqs={A, B, C, D};
		final String[] names={"famA", "famB", "famC", "famD"};
		final ArrayList<MarkerFamily> families=new ArrayList<MarkerFamily>();
		for(int i=0; i<seqs.length; i++){
			final CopyNumberDistribution dist=new CopyNumberDistribution();
			dist.add(1);//nominal: carried once in the (single) reference genome
			final ProteinSequence rep=new ProteinSequence(names[i], seqs[i]);
			families.add(new MarkerFamily(i, rep, dist, true));//selected single-copy
		}
		final List<String> ids=new ArrayList<String>();
		ids.add("refGenome1");
		final MarkerSetProvenance prov=new MarkerSetProvenance("2026-07-31T00:00:00Z",
			0.97, 90.0, 0.8, 5, ids, "NA");
		return new MarkerSet("Bacteria", "test-v1", prov, 1, families);
	}

	/** Throws with the message if the condition is false. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new RuntimeException("TEST FAILED: "+message);}
	}
}
