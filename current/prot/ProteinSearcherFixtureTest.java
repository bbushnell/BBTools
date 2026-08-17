package prot;

import java.util.ArrayList;
import java.util.List;

/**
 * Disposable differential fixture for the ProteinSearcher de-boxing rewrite (RULE #12
 * baseline). Not a permanent test -- exercises seeding (shared/disjoint k-mers), the
 * short-query (less-than-k) fallback path, reducedSeed, maxTargetSeqs culling, and
 * minSeedHits>1, then dumps every hit's TSV in the frozen deterministic order.
 */
public class ProteinSearcherFixtureTest {

	public static void main(String[] args){
		List<ProteinSequence> targets=new ArrayList<ProteinSequence>();
		targets.add(new ProteinSequence("t1", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR"));
		targets.add(new ProteinSequence("t2", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKG"));
		targets.add(new ProteinSequence("t3", "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWA"));
		targets.add(new ProteinSequence("t4", "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGAT"));
		targets.add(new ProteinSequence("t5", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQPDGKG"));
		final char[] fillSuffix={'A','C','D','E','F','G','H'};
		for(int i=0; i<50; i++){
			targets.add(new ProteinSequence("filler"+i, "MDLKVHNQTAWRSCFYGPIEDLKVHNQTAWRSCFYGPIEDMDLKVHNQTAWRSCFYGPIE"+fillSuffix[i%7]));
		}

		List<ProteinSequence> queries=new ArrayList<ProteinSequence>();
		queries.add(new ProteinSequence("q1_full_match", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR"));
		queries.add(new ProteinSequence("q2_partial", "MKTAYIAKQRQISFVKSHFSRQLEERL"));
		queries.add(new ProteinSequence("q3_short_below_k", "MKT"));
		queries.add(new ProteinSequence("q4_disjoint", "WWWWWWWWWWWWWWWWWWWWWWWWWWWW"));
		queries.add(new ProteinSequence("q5_after_short", "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWA"));
		queries.add(new ProteinSequence("q6_manyhits", "MDLKVHNQTAWRSCFYGPIEDLKVHNQTAWRSCFYGPIEDMDLKVHNQTAWRSCFYGPIEA"));

		double searchSpace0;

		//Scenario A: default params (k=5, minSeedHits=1)
		ProteinSearcher a=new ProteinSearcher();
		List<ProteinHit> ra=a.search(queries, targets);
		System.out.println("===SCENARIO_A_default===");
		for(ProteinHit h : ra){System.out.println(h.toTsv());}

		//Scenario B: reducedSeed=true
		ProteinSearcher b=new ProteinSearcher();
		b.reducedSeed=true;
		List<ProteinHit> rb=b.search(queries, targets);
		System.out.println("===SCENARIO_B_reducedSeed===");
		for(ProteinHit h : rb){System.out.println(h.toTsv());}

		//Scenario C: minSeedHits=3, maxTargetSeqs=2
		ProteinSearcher c=new ProteinSearcher();
		c.minSeedHits=3;
		c.maxTargetSeqs=2;
		List<ProteinHit> rc=c.search(queries, targets);
		System.out.println("===SCENARIO_C_minseed3_maxtargets2===");
		for(ProteinHit h : rc){System.out.println(h.toTsv());}

		//Scenario D: short query FIRST (exercises the touched-latch fallback path),
		//then a normal query right after, to confirm the seed filter still works post-fallback.
		List<ProteinSequence> qOrderTest=new ArrayList<ProteinSequence>();
		qOrderTest.add(new ProteinSequence("d1_short", "MK"));
		qOrderTest.add(new ProteinSequence("d2_normal", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR"));
		ProteinSearcher d=new ProteinSearcher();
		List<ProteinHit> rd=d.search(qOrderTest, targets);
		System.out.println("===SCENARIO_D_short_then_normal===");
		for(ProteinHit h : rd){System.out.println(h.toTsv());}

		System.out.println("===COUNTS===");
		System.out.println("A="+ra.size()+" B="+rb.size()+" C="+rc.size()+" D="+rd.size());
	}
}
