package test.tax;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

import structures.ByteBuilder;
import tax.TaxNode;
import tax.TaxTree;
import tax.TaxTreeText;

/** Finite parent-chain and count-gathering regressions, independent of rank ordering.
 * @author Brian Bushnell, Shinobu */
public final class TaxTreeContractTest{

	/** Creates an instance with no retained state. */
	public TaxTreeContractTest(){}

	/** Runs known-answer cases; expect=baseline records the original defects.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If temporary tree creation, loading or deletion fails */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter "+arg);}
		}
		int mismatches=0;
		for(int sharedRank : new int[]{TaxTree.NO_RANK_E, TaxTree.SPECIES_E, TaxTree.GENUS_E}){
			for(int firstRank : new int[]{TaxTree.NO_RANK_E, TaxTree.SPECIES_E, TaxTree.FAMILY_E}){
				for(int secondRank : new int[]{TaxTree.NO_RANK_E, TaxTree.SPECIES_E, TaxTree.FAMILY_E}){
					TaxTree tree=tree(sharedRank, firstRank, secondRank);
					for(int a=1; a<tree.nodes.length; a++){
						for(int b=1; b<tree.nodes.length; b++){
							int expected=ancestorOracle(tree, a, b);
							int observed=tree.commonAncestor(a, b);
							if(observed!=expected){
								mismatches++;
								if(!baseline){throw new AssertionError("LCA of "+a+","+b+" expected "+expected+" but got "+observed+" at shared rank "+sharedRank);}
							}
						}
					}
					check(tree.commonAncestor(4, 4)==4, "A node is its own common ancestor");
					check(tree.commonAncestor(1, 5)==1, "Root is an ancestor of every node");
				}
			}
		}
		if(baseline){check(mismatches>0, "Baseline must expose the rank-order assumption");}
		else{check(mismatches==0, "Every pair must match its actual parent-chain intersection");}
		TaxTree counts=tree(TaxTree.GENUS_E);
		counts.incrementRaw(4, 7);
		counts.percolateUp();
		ArrayList<TaxNode> control=counts.gatherNodesAtLeastLimit(5, TaxTree.NO_RANK, TaxTree.LIFE);
		check(control.size()==1 && control.get(0).id==4, "Explicit full-rank gathering must select leaf4 only");
		counts=tree(TaxTree.GENUS_E);
		counts.incrementRaw(4, 7);
		counts.percolateUp();
		boolean failed=false;
		try{
			ArrayList<TaxNode> result=counts.gatherNodesAtLeastLimit(5);
			check(result.size()==1 && result.get(0).id==4, "Default gathering must equal explicit full-rank gathering");
		}catch(ArrayIndexOutOfBoundsException expected){
			failed=true;
			if(!baseline){throw expected;}
		}
		check(failed==baseline, "Default gather failure must occur only in the captured baseline");
		System.out.println("PASS lca_mismatches="+mismatches+" default_gather_failed="+failed);
	}

	/** Builds a five-node tree whose topology stays fixed while its shared rank varies.
	 * @param sharedRank Extended rank of node3, the parent shared by nodes4 and5
	 * @return Fresh tree with independent zeroed counters
	 * @throws Exception If temporary tree creation, loading or deletion fails */
	public static TaxTree tree(int sharedRank) throws Exception{
		return tree(sharedRank, TaxTree.SPECIES_E, TaxTree.SPECIES_E);
	}

	/** Builds the same parent topology with independently variable child ranks.
	 * @param sharedRank Extended rank of shared parent node3
	 * @param firstRank Extended rank of child node4
	 * @param secondRank Extended rank of child node5
	 * @return Fresh connected tree with independent zeroed counters
	 * @throws Exception If temporary tree creation, loading or deletion fails */
	private static TaxTree tree(int sharedRank, int firstRank, int secondRank) throws Exception{
		ByteBuilder bb=new ByteBuilder();
		bb.append("##TaxTree2\t1\n##Encoding\tA48\n##Tree\t6\t5\n");
		bb.append("##Config\t0\t0\t0\t0\t0\n##MergedMap\t0\t0\n");
		bb.append("##LevelExtended\tTID\tPID\tNumChildren\tMinParentLevelExtended\tMaxChildLevelExtended\tFlag\tName\n");
		row(bb, 1, 1, TaxTree.LIFE_E, 1, "Life");
		row(bb, 2, 1, TaxTree.DOMAIN_E, 1, "Bacteria");
		row(bb, 3, 2, sharedRank, 2, "Shared");
		row(bb, 4, 3, firstRank, 0, "First");
		row(bb, 5, 3, secondRank, 0, "Second");
		Path file=Files.createTempFile("tax-parent-contract-", ".tsv");
		try{
			Files.write(file, bb.toBytes());
			return TaxTreeText.load(file.toString(), 1);
		}finally{Files.delete(file);}
	}

	/** Appends an A48 tree row with topology supplied explicitly by the fixture.
	 * @param bb Destination builder
	 * @param id Native TaxID
	 * @param parent Parent TaxID; the root points to itself
	 * @param rank Extended taxonomic rank
	 * @param children Number of immediate children
	 * @param name Literal node name */
	private static void row(ByteBuilder bb, int id, int parent, int rank, int children, String name){
		bb.appendA48(rank).tab().appendA48(id).tab().appendA48(parent).tab().appendA48(children)
				.tab().appendA48(TaxTree.LIFE_E).tab().appendA48(0).tab().appendA48(rank).tab().append(name).nl();
	}

	/** Finds the first shared ancestor by literal path membership, without rank comparisons.
	 * @param tree Connected fixture tree with native IDs indexing its node array
	 * @param a First native TaxID
	 * @param b Second native TaxID
	 * @return Closest shared ancestor's native TaxID */
	private static int ancestorOracle(TaxTree tree, int a, int b){
		boolean[] ancestors=new boolean[tree.nodes.length];
		while(true){
			ancestors[a]=true;
			int parent=tree.getNode(a).pid;
			if(parent==a){break;}
			a=parent;
		}
		while(!ancestors[b]){b=tree.getNode(b).pid;}
		return b;
	}

	/** Keeps test checks enabled even when the tested JVM runs with -da.
	 * @param condition Required fixture invariant
	 * @param message Diagnostic for a violated invariant */
	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
