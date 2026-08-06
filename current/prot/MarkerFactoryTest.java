package prot;

import java.util.ArrayList;
import java.util.List;

/**
 * Synthetic self-test for {@link MarkerFactory}. Constructs a small set of
 * Bacteria "genomes" with a known ground truth and verifies the factory selects
 * exactly the universal single-copy families and computes correct prevalence and
 * copy-number distributions. Fails loud (throws) on any mismatch; prints the
 * marker-set table and PASS on success.
 *
 * <p>Ground truth over 6 genomes (G1..G6):</p>
 * <ul>
 * <li><b>markerA, markerB, markerC</b> — present exactly once in all 6 -&gt;
 * MUST be selected.</li>
 * <li><b>rareD</b> — present exactly once in only G1, G2 (prevalence 2/6) -&gt;
 * NOT selected (low prevalence).</li>
 * <li><b>multiE</b> — present twice in all 6 -&gt; NOT selected (multi-copy;
 * exactly-once fraction 0).</li>
 * </ul>
 *
 * <p>Run: {@code java -cp current prot.MarkerFactoryTest}</p>
 *
 * @author Eru
 */
public final class MarkerFactoryTest {

	//Five mutually-dissimilar ~60-aa base sequences; copies are identical within a
	//family (so clustering groups them) and distinct between families.
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
		final ArrayList<GenomeProteins> genomes=new ArrayList<GenomeProteins>();
		for(int gi=1; gi<=6; gi++){
			final String gid="G"+gi;
			final ArrayList<ProteinSequence> prots=new ArrayList<ProteinSequence>();
			//Universal single-copy families: one copy each in every genome.
			prots.add(new ProteinSequence(gid+"_markerA", A));
			prots.add(new ProteinSequence(gid+"_markerB", B));
			prots.add(new ProteinSequence(gid+"_markerC", C));
			//Multi-copy family: two copies in every genome.
			prots.add(new ProteinSequence(gid+"_multiE_1", E));
			prots.add(new ProteinSequence(gid+"_multiE_2", E));
			//Low-prevalence family: only in G1 and G2, one copy.
			if(gi<=2){prots.add(new ProteinSequence(gid+"_rareD", D));}
			genomes.add(new GenomeProteins(gid, "Bacteria", "synthetic", prots));
		}

		final MarkerFactory factory=new MarkerFactory();
		final List<MarkerSet> sets=factory.build(genomes, "test-v1", "2026-07-31T00:00:00Z", "NA");

		System.out.println("Marker sets built: "+sets.size());
		check(sets.size()==1, "expected 1 domain marker set, got "+sets.size());
		final MarkerSet ms=sets.get(0);
		System.out.println("Domain="+ms.domain+" version="+ms.version+
			" genomes="+ms.genomeCount+" families="+ms.families.size()+
			" selected="+ms.selectedCount());
		System.out.println("Build timestamp (provenance)="+ms.provenance.buildTimestamp+
			" selectionThreshold="+ms.provenance.selectionThreshold);
		System.out.println(MarkerFactoryCLI.header());
		for(final MarkerFamily f : ms.families){System.out.println(MarkerFactoryCLI.row(ms, f));}

		check(ms.domain.equals("Bacteria"), "wrong domain: "+ms.domain);
		check(ms.genomeCount==6, "expected 6 genomes, got "+ms.genomeCount);
		//5 distinct families should be present in this domain.
		check(ms.families.size()==5, "expected 5 families, got "+ms.families.size());

		//Classify families by their ground-truth signature and verify each.
		int universal=0, rare=0, multi=0, selected=0;
		for(final MarkerFamily f : ms.families){
			final int[] b=f.dist.bins;
			if(f.selectedSingleCopy){selected++;}
			if(b[1]==6 && f.dist.present()==6){
				//Universal single-copy: all 6 genomes exactly once.
				universal++;
				check(f.selectedSingleCopy, "universal family "+f.familyId+" not selected");
				check(near(f.prevalence(), 1.0), "universal prevalence != 1.0: "+f.prevalence());
				check(near(f.fractionExactlyOnce(), 1.0),
					"universal exactly-once != 1.0: "+f.fractionExactlyOnce());
			}else if(b[0]==4 && b[1]==2){
				//rareD: absent in 4, once in 2.
				rare++;
				check(!f.selectedSingleCopy, "low-prevalence family "+f.familyId+" wrongly selected");
				check(near(f.prevalence(), 2.0/6.0), "rare prevalence wrong: "+f.prevalence());
				check(near(f.fractionExactlyOnce(), 2.0/6.0),
					"rare exactly-once wrong: "+f.fractionExactlyOnce());
			}else if(b[2]==6){
				//multiE: two copies in all 6.
				multi++;
				check(!f.selectedSingleCopy, "multi-copy family "+f.familyId+" wrongly selected");
				check(near(f.prevalence(), 1.0), "multi prevalence != 1.0: "+f.prevalence());
				check(near(f.fractionExactlyOnce(), 0.0),
					"multi exactly-once != 0.0: "+f.fractionExactlyOnce());
			}else{
				throw new RuntimeException("Unexpected family shape: family "+f.familyId+
					" bins=["+b[0]+","+b[1]+","+b[2]+","+b[3]+","+b[4]+"]");
			}
		}

		check(universal==3, "expected 3 universal single-copy families, got "+universal);
		check(rare==1, "expected 1 low-prevalence family, got "+rare);
		check(multi==1, "expected 1 multi-copy family, got "+multi);
		check(selected==3, "expected exactly 3 selected markers, got "+selected);
		check(ms.selectedCount()==3, "selectedCount()!=3: "+ms.selectedCount());

		System.out.println("\nResult: universal(selected)="+universal+" rare(excluded)="+rare+
			" multi(excluded)="+multi+" total_selected="+selected);
		System.out.println("PASS: exactly the 3 universal single-copy families were selected; "+
			"low-prevalence and multi-copy families were correctly excluded.");
	}

	/** True if two doubles are within 1e-9. */
	private static boolean near(final double a, final double b){return Math.abs(a-b)<1e-9;}

	/** Throws with the message if the condition is false. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new RuntimeException("TEST FAILED: "+message);}
	}
}
