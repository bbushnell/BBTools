package prot;

import java.util.ArrayList;
import java.util.List;

/**
 * Disposable differential fixture for the ProteinFamilyClusterer HashMap-to-array-index
 * de-boxing (RULE #12 baseline). 24 sequences across several divergent families plus a
 * few singletons/orphans, enough width to stress the groups-map iteration order (dense
 * int keys 0..n) that the rewrite changes from HashMap.entrySet() to array-index order.
 */
public class ProteinFamilyClustererFixtureTest {

	public static void main(String[] args){
		List<ProteinSequence> seqs=new ArrayList<ProteinSequence>();
		String[] famA={
			"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR",
			"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKG",
			"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKA",
			"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKC"};
		String[] famB={
			"MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWA",
			"MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWG",
			"MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWC"};
		String[] famC={
			"GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGAT",
			"GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGAA",
			"GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGAC"};
		String[] famD={
			"ACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC",
			"ACACACACACACACACACACACACACACACACACACACACACACACACACACACACAG",
			"ACACACACACACACACACACACACACACACACACACACACACACACACACACACACAA"};
		String[] famE={
			"CAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGA",
			"CAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGG",
			"CAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGC"};
		String[] famF={
			"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC",
			"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGA",
			"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGG"};
		String[][] allFams={famA, famB, famC, famD, famE, famF};
		char[] famNames={'A','B','C','D','E','F'};
		for(int f=0; f<allFams.length; f++){
			for(int i=0; i<allFams[f].length; i++){
				seqs.add(new ProteinSequence("f"+famNames[f]+"_"+i, allFams[f][i]));
			}
		}
		//A few orphan singletons, distinct from everything else.
		seqs.add(new ProteinSequence("orphan1", "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW"));
		seqs.add(new ProteinSequence("orphan2", "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY"));

		final ProteinFamilyClusterer pfc=new ProteinFamilyClusterer();
		final List<ProteinFamily> families=pfc.cluster(seqs);

		System.out.println("total_input="+seqs.size()+" total_families="+families.size());
		for(ProteinFamily fam : families){
			final StringBuilder sb=new StringBuilder();
			sb.append("family id=").append(fam.id).append(" consensus_id=").append(fam.consensus.id);
			sb.append(" members=[");
			for(int i=0; i<fam.members.size(); i++){
				if(i>0){sb.append(',');}
				sb.append(fam.members.get(i).id);
			}
			sb.append(']');
			System.out.println(sb.toString());
		}
	}
}
