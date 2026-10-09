package tax;

import java.io.File;
import java.io.PrintWriter;

import structures.IntHashMap;

/** Tests auxiliary fake-TaxID overlays. */
public final class TaxTreeAuxTest {

	public static void main(String[] args) throws Exception{
		testLookupAndLineage();
		testDuplicate();
		testMissingParent();
		testLowFakeID();
	}

	private static void testLookupAndLineage() throws Exception{
		TaxTree tree=load("fake_taxid\tparent_taxid\trank\tname\n"+
				"2000000000\t2\tphylum\tFakeota\n"+
				"2000000001\t2000000000\tclass\tFakeia\n"+
				"2000000002\t2000000001\tgenus\tFakeus\n"+
				"2000000003\t2000000002\tstrain\tGenomeA\n"+
				"2000000004\t2000000002\tstrain\tGenomeB\n");
		check(tree.getNode(2000000003).name.equals("GenomeA"), "durable fake leaf did not resolve");
		check(tree.parseNodeFromHeader("2000000003", true).name.equals("GenomeA"), "10-digit fake TaxID did not parse");
		expectParseFailure(tree, "2147483648");
		check(tree.getNode(562)==null, "private packed slot leaked as a public TaxID");
		check(tree.getParentID(2000000003)==2000000002, "parent was not returned as a durable fake TaxID");
		check(tree.commonAncestor(2000000003, 2000000004)==2000000002, "fake leaves did not coalesce at their fake genus");
		check(tree.descendsFrom(2000000003, 2000000000), "leaf does not descend from fake phylum");
		check(tree.descendsFrom(2000000003, TaxTree.BACTERIA_ID), "leaf does not descend from real graft point");
		check(tree.getNode(TaxTree.BACTERIA_ID).numChildren==2, "real parent child count was not incremented");
		tree.incrementRaw(2000000003, 7);
		tree.percolateUp();
		check(tree.getNode(2000000003).countSum==7, "leaf count did not percolate to itself");
		check(tree.getNode(TaxTree.BACTERIA_ID).countSum==7, "leaf count did not percolate through fake parents");
	}

	private static void testDuplicate() throws Exception{
		expectFailure("fake_taxid\tparent_taxid\trank\tname\n"+
				"2000000000\t2\tphylum\tFakeota\n"+
				"2000000000\t2\tphylum\tOtherota\n");
	}

	private static void testMissingParent() throws Exception{
		expectFailure("fake_taxid\tparent_taxid\trank\tname\n"+
				"2000000000\t2000000999\tphylum\tFakeota\n");
	}

	private static void testLowFakeID() throws Exception{
		expectFailure("fake_taxid\tparent_taxid\trank\tname\n"+
				"6000000\t2\tphylum\tFakeota\n");
	}

	private static void expectFailure(String body) throws Exception{
		boolean ok=false;
		try{
			load(body);
		}catch(RuntimeException e){
			ok=true;
		}
		check(ok, "auxiliary tree load should have failed");
	}

	private static void expectParseFailure(TaxTree tree, String s){
		boolean ok=false;
		try{
			tree.parseNodeFromHeader(s, true);
		}catch(RuntimeException e){
			ok=true;
		}
		check(ok, s+" should have failed during TaxID parsing");
	}

	private static TaxTree load(String body) throws Exception{
		File f=File.createTempFile("taxtree_aux_test_", ".tsv");
		f.deleteOnExit();
		PrintWriter pw=new PrintWriter(f);
		pw.print(body);
		pw.close();
		return TaxTreeAux.apply(baseTree(), f.getAbsolutePath());
	}

	private static TaxTree baseTree(){
		TaxNode[] nodes=new TaxNode[562];
		nodes[TaxTree.LIFE_ID]=new TaxNode(TaxTree.LIFE_ID, TaxTree.LIFE_ID,
				TaxTree.LIFE, TaxTree.LIFE_E, "Life");
		nodes[TaxTree.BACTERIA_ID]=new TaxNode(TaxTree.BACTERIA_ID, TaxTree.LIFE_ID,
				TaxTree.DOMAIN, TaxTree.DOMAIN_E, "Bacteria");
		nodes[561]=new TaxNode(561, TaxTree.BACTERIA_ID,
				TaxTree.GENUS, TaxTree.GENUS_E, "Escherichia");
		nodes[TaxTree.LIFE_ID].numChildren=1;
		nodes[TaxTree.BACTERIA_ID].numChildren=1;
		return new TaxTree(nodes, new IntHashMap(), 0, false, false, false, 0);
	}

	private static void check(boolean b, String message){
		if(!b){throw new AssertionError(message);}
	}

}
