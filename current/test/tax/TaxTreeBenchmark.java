package test.tax;

import java.util.ArrayList;
import java.util.Random;

import parse.Parse;
import tax.TaxNode;
import tax.TaxTree;

/** Repeated common-ancestor workload with deterministic pairs and a topology oracle.
 * Input loading and oracle construction are outside the measured loop.
 * @author Brian Bushnell, Shinobu */
public final class TaxTreeBenchmark{

	/** Creates an instance with no retained state. */
	public TaxTreeBenchmark(){}

	/** Measures lookup on a supplied tree; each JVM reports its observed checksum.
	 * @param args in=tree, pairs=65536, rounds=64, related=f, expect=baseline|fixed,
	 * and optional diagnostics=fallback */
	public static void main(String[] args){
		String input="auto";
		String diagnostics=null;
		int pairs=65536, rounds=64;
		boolean related=false, baseline=false;
		for(String arg : args){
			String[] split=arg.split("=", 2);
			if(split.length!=2){throw new IllegalArgumentException("Expected key=value: "+arg);}
			String key=split[0], value=split[1];
			if(key.equals("in")){input=value;}
			else if(key.equals("pairs")){pairs=Parse.parseIntKMG(value);}
			else if(key.equals("rounds")){rounds=Parse.parseIntKMG(value);}
			else if(key.equals("related")){related=Parse.parseBoolean(value);}
			else if(key.equals("diagnostics")){
				if(!value.equals("fallback")){throw new IllegalArgumentException(arg);}
				diagnostics=value;
			}
			else if(key.equals("expect")){
				if(!value.equals("baseline") && !value.equals("fixed")){throw new IllegalArgumentException(arg);}
				baseline=value.equals("baseline");
			}else{throw new IllegalArgumentException("Unknown parameter "+arg);}
		}
		if(pairs<1 || rounds<1){throw new IllegalArgumentException("pairs and rounds must be positive");}
		TaxTree.ALLOW_SHARED_TREE=false;
		TaxTree.LOAD_FAKE_TREE=false;
		TaxTree tree=TaxTree.loadTaxTree(input, null, false, false);
		TaxNode[] a=new TaxNode[pairs], b=new TaxNode[pairs];
		Random random=new Random(812376L);
		ArrayList<TaxNode> ancestors=new ArrayList<TaxNode>(64);
		long expected=0;
		int mismatches=0;
		for(int i=0; i<pairs; i++){
			a[i]=choose(tree, random);
			b[i]=related ? ancestor(tree, a[i], 1+random.nextInt(8)) : choose(tree, random);
			int wanted=oracle(tree, a[i], b[i], ancestors);
			int actual=tree.commonAncestor(a[i], b[i]).id;
			if(actual!=wanted){mismatches++;}
			expected+=wanted;
		}
		if(!baseline && mismatches!=0){throw new AssertionError("Topology oracle mismatches: "+mismatches);}
		if(diagnostics!=null){
			final int fallbackPairs=countFallbackPairs(tree, a, b);
			System.out.println("diagnostic=fallback type=structural fallback_pairs="+fallbackPairs+
					" pairs="+pairs+" fallback_rate="+(fallbackPairs/(double)pairs)+
					" related="+related+" oracle_mismatches="+mismatches);
			return;
		}
		long warm=run(tree, a, b, 4);
		long start=System.nanoTime();
		long sum=run(tree, a, b, rounds);
		long elapsed=System.nanoTime()-start;
		if(!baseline && sum!=expected*rounds){throw new AssertionError("Measured checksum disagrees with topology oracle");}
		if(sum*4!=warm*rounds){throw new AssertionError("Measured and warmup rounds disagree");}
		System.out.println("checksum="+sum+" expected="+(expected*rounds)+" oracle_mismatches="+mismatches+
				" pairs="+pairs+" rounds="+rounds+" related="+related+" compute_ns="+elapsed);
	}

	/** Selects uniformly among occupied native array slots, retrying holes.
	 * @param tree Nonempty native tree
	 * @param random Deterministic generator
	 * @return Borrowed node */
	private static TaxNode choose(TaxTree tree, Random random){
		TaxNode node;
		do{node=tree.nodes[random.nextInt(tree.nodes.length)];}while(node==null);
		return node;
	}

	/** Walks at most hops parent links, stopping at a self-parent root.
	 * @param tree Tree containing the complete parent chain
	 * @param node Starting node
	 * @param hops Maximum number of parent steps
	 * @return Borrowed ancestor, possibly the starting node or root */
	private static TaxNode ancestor(TaxTree tree, TaxNode node, int hops){
		for(int i=0; i<hops && node.pid!=node.id; i++){node=tree.getNode(node.pid);}
		return node;
	}

	/** Finds first shared lineage membership without using rank or production LCA.
	 * @param tree Tree containing both parent chains
	 * @param a First node
	 * @param b Second node
	 * @param path Reusable scratch list, cleared and filled with the first lineage
	 * @return Closest shared ancestor's TaxID
	 * @throws AssertionError If a chain cycles or the two chains have no shared node */
	private static int oracle(TaxTree tree, TaxNode a, TaxNode b, ArrayList<TaxNode> path){
		path.clear();
		for(int steps=0;; steps++){
			if(steps>=tree.nodes.length){throw new AssertionError("Cycle in first oracle lineage");}
			path.add(a);
			if(a.pid==a.id){break;}
			a=tree.getNode(a.pid);
		}
		for(int steps=0; !path.contains(b); steps++){
			if(steps>=tree.nodes.length || b.pid==b.id){throw new AssertionError("No common ancestor in oracle");}
			b=tree.getNode(b.pid);
		}
		return b.id;
	}

	/** Counts pairs that would enter TaxTree's rank-guard fallback.
	 * This structural diagnostic mirrors the guard before timing begins; it is not
	 * hot-loop profiling or a speed test.
	 * @param tree Tree containing all parents for both arrays
	 * @param a First node of each pair
	 * @param b Second nodes, with the same array length as a
	 * @return Number of pairs whose rank walk reaches a non-increasing parent edge */
	private static int countFallbackPairs(TaxTree tree, TaxNode[] a, TaxNode[] b){
		int count=0;
		for(int i=0; i<a.length; i++){
			if(hasRankGuardFallback(tree, a[i], b[i])){count++;}
		}
		return count;
	}

	/** Mirrors TaxTree.commonAncestor's fallback guard without calling the fallback.
	 * @param tree Tree containing both parent chains
	 * @param a First node
	 * @param b Second node
	 * @return True if the rank walk would switch to the depth-aligned fallback */
	private static boolean hasRankGuardFallback(TaxTree tree, TaxNode a, TaxNode b){
		while(a!=b){
			if(a.levelExtended<b.levelExtended){
				TaxNode parent=tree.getNode(a.pid);
				if(parent.levelExtended<=a.levelExtended){return true;}
				a=parent;
			}else{
				TaxNode parent=tree.getNode(b.pid);
				if(parent.levelExtended<=b.levelExtended){return true;}
				b=parent;
			}
		}
		return false;
	}

	/** Computes repeated lookups without allocating per pair or per round.
	 * @param tree Tree queried for common ancestors
	 * @param a First node of each pair
	 * @param b Second nodes, with the same array length as a
	 * @param rounds Number of complete passes over the pairs
	 * @return Long checksum of returned ancestor IDs */
	private static long run(TaxTree tree, TaxNode[] a, TaxNode[] b, int rounds){
		long sum=0;
		for(int round=0; round<rounds; round++){
			for(int i=0; i<a.length; i++){sum+=tree.commonAncestor(a[i], b[i]).id;}
		}
		return sum;
	}
}
