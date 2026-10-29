package test.tax;

import java.io.File;
import java.io.PrintWriter;

import structures.IntHashMap;
import tax.TaxNode;
import tax.TaxTree;
import tax.TaxTreeAux;

/** Tests auxiliary fake-TaxID overlays. */
public final class TaxTreeAuxTest{

	/** Runs finite auxiliary-ID lookup, ancestry and invalid-overlay checks.
	 * @param args Optional expectation flag, expect=baseline or expect=fixed
	 * @throws Exception If a temporary fixture cannot be created */
	public static void main(String[] args) throws Exception{
		String expectation=parseExpectation(args);
		testLookupAndLineage();
		testSuccessiveOverlayNativeParent(expectation);
		testSuccessiveOverlayAuxParent(expectation);
		testSuccessiveOverlayPrivateSlotParent(expectation);
		testSuccessiveOverlayNoRankPropagation(expectation);
		testDuplicate();
		testMissingParent();
		testLowFakeID();
	}

	/** Prevents construction of this standalone fixture. */
	private TaxTreeAuxTest(){}

	/** First overlay used by successive-apply checks. */
	private static final String FIRST_OVERLAY=
			"fake_taxid\tparent_taxid\trank\tname\n"+
			"2000000000\t2\tphylum\tFakeota\n"+
			"2000000001\t2000000000\tclass\tFakeia\n"+
			"2000000002\t2000000001\tgenus\tFakeus\n"+
			"2000000003\t2000000002\tstrain\tGenomeA\n"+
			"2000000004\t2000000002\tstrain\tGenomeB\n";

	/** Second overlay whose parent is native Bacteria. */
	private static final String SECOND_NATIVE_PARENT=
			"fake_taxid\tparent_taxid\trank\tname\n"+
			"2000000010\t2\tstrain\tGenomeC\n";

	/** Second overlay whose parent is a durable ID from the first overlay. */
	private static final String SECOND_AUX_PARENT=
			"fake_taxid\tparent_taxid\trank\tname\n"+
			"2000000010\t2000000000\tstrain\tGenomeC\n";

	/** Second overlay whose parent is an old private packed slot, not a public TaxID. */
	private static final String SECOND_PRIVATE_SLOT_PARENT=
			"fake_taxid\tparent_taxid\trank\tname\n"+
			"2000000010\t562\tstrain\tGenomeC\n";

	/** First overlay with no-rank ancestors used by propagation checks. */
	private static final String FIRST_NORANK_OVERLAY=
			"fake_taxid\tparent_taxid\trank\tname\n"+
			"2000000020\t2\tno rank\tNoRankParent\n"+
			"2000000021\t2000000020\tno rank\tNoRankChild\n";

	/** Second overlay whose ranked node hangs below old no-rank ancestors. */
	private static final String SECOND_NORANK_DESCENDANT=
			"fake_taxid\tparent_taxid\trank\tname\n"+
			"2000000022\t2000000021\tstrain\tGenomeD\n";

	/** Checks durable auxiliary IDs, lineage, child counts and count propagation.
	 * @throws Exception If the fixture cannot be written */
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

	/** Checks a second overlay attached to a native parent.
	 * @param expectation Expected behavior: baseline or fixed
	 * @throws Exception If a fixture cannot be written */
	private static void testSuccessiveOverlayNativeParent(String expectation) throws Exception{
		if("baseline".equals(expectation)){
			String message=expectNativeParentExposedSlotFailure();
			System.out.println("successive-native-parent expect=baseline failure="+message);
		}else{
			TaxTree tree=apply(apply(baseTree(), FIRST_OVERLAY), SECOND_NATIVE_PARENT);
			check(tree.getNode(2000000003).name.equals("GenomeA"),
					"first-overlay durable fake leaf did not survive a second overlay");
			check(tree.parseNodeFromHeader("2000000003", true).name.equals("GenomeA"),
					"first-overlay header TaxID did not survive a second overlay");
			check(tree.getNode(562)==null, "first old private slot leaked after a second overlay");
			check(tree.getNode(563)==null, "second old private slot leaked after a second overlay");
			check(tree.getParentID(2000000003)==2000000002,
					"old fake parent was not returned as a durable fake TaxID");
			check(tree.commonAncestor(2000000003, 2000000004)==2000000002,
					"old fake leaves did not coalesce after a second overlay");
			check(tree.descendsFrom(2000000003, 2000000000),
					"old fake leaf lost ancestry to the first fake phylum");
			check(tree.descendsFrom(2000000003, TaxTree.BACTERIA_ID),
					"old fake leaf lost ancestry to the real graft point");
			check(tree.getNode(2000000010).name.equals("GenomeC"),
					"second overlay attached to a native parent did not resolve");
			check(tree.getParentID(2000000010)==TaxTree.BACTERIA_ID,
					"new native-parent overlay returned the wrong parent");
			System.out.println("successive-native-parent expect=fixed oldID=GenomeA "+
					"slots=hidden newID=GenomeC parent="+TaxTree.BACTERIA_ID);
		}
	}

	/** Checks a second overlay attached to a prior auxiliary parent.
	 * @param expectation Expected behavior: baseline or fixed
	 * @throws Exception If a fixture cannot be written */
	private static void testSuccessiveOverlayAuxParent(String expectation) throws Exception{
		if("baseline".equals(expectation)){
			String message=expectUnknownPriorAuxParent();
			System.out.println("successive-aux-parent expect=baseline failure="+message);
		}else{
			TaxTree tree=apply(apply(baseTree(), FIRST_OVERLAY), SECOND_AUX_PARENT);
			check(tree.getNode(2000000010).name.equals("GenomeC"),
					"second overlay attached to an old fake parent did not resolve");
			check(tree.getParentID(2000000010)==2000000000,
					"new fake child did not keep its old fake parent");
			check(tree.descendsFrom(2000000010, 2000000000),
					"new fake child does not descend from its old fake parent");
			check(tree.descendsFrom(2000000010, TaxTree.BACTERIA_ID),
					"new fake child does not descend from the native graft point");
			System.out.println("successive-aux-parent expect=fixed newID=GenomeC "+
					"parent=2000000000 nativeAncestor="+TaxTree.BACTERIA_ID);
		}
	}

	/** Checks fixed rejection of an old private slot as a parent TaxID.
	 * Baseline never reaches this path because second overlays fail first.
	 * @param expectation Expected behavior: baseline or fixed
	 * @throws Exception If a fixture cannot be written */
	private static void testSuccessiveOverlayPrivateSlotParent(String expectation) throws Exception{
		if("baseline".equals(expectation)){return;}
		TaxTree first=apply(baseTree(), FIRST_OVERLAY);
		check(first.getNode(2000000003).name.equals("GenomeA"),
				"first overlay changed before invalid private-slot parent test");
		check(first.getNode(562)==null,
				"private packed slot was externally visible before the invalid parent test");
		check(first.externalToInternalMap.size()==5,
				"first overlay should have exactly five durable ID mappings");
		String message=expectPrivateSlotParentFailure(first);
		check(first.externalToInternalMap.size()==5 && !first.externalToInternalMap.contains(2000000010),
				"rejected overlay inserted a mapping into the borrowed base tree");
		check(first.getNode(2000000003).name.equals("GenomeA"),
				"first overlay lookup changed after rejected private-slot parent");
		check(first.getNode(562)==null,
				"private packed slot leaked after rejected private-slot parent");
		check(first.getParentID(2000000003)==2000000002,
				"first overlay parent map changed after rejected private-slot parent");
		System.out.println("successive-private-slot-parent expect=fixed failure="+message+
				" base=unchanged");
	}

	/** Checks fixed propagation over old auxiliary no-rank edges.
	 * Baseline never reaches this path because old auxiliary parents do not resolve.
	 * @param expectation Expected behavior: baseline or fixed
	 * @throws Exception If a fixture cannot be written */
	private static void testSuccessiveOverlayNoRankPropagation(String expectation) throws Exception{
		if("baseline".equals(expectation)){return;}
		TaxTree tree=apply(apply(baseTree(), FIRST_NORANK_OVERLAY), SECOND_NORANK_DESCENDANT);
		TaxNode directParent=tree.getNode(2000000021);
		TaxNode grandparent=tree.getNode(2000000020);
		check(directParent.maxChildLevelExtended==TaxTree.STRAIN_E,
				"direct old no-rank parent did not see new ranked child");
		check(grandparent.maxChildLevelExtended==TaxTree.STRAIN_E,
				"old no-rank grandparent did not see new ranked descendant");
		System.out.println("successive-norank-propagation expect=fixed directMax="+
				directParent.maxChildLevelExtended+" grandparentMax="+grandparent.maxChildLevelExtended);
	}

	/** Checks the baseline failure when a second overlay exposes old packed slots.
	 * @return Exact exception message from the baseline failure
	 * @throws Exception If a fixture cannot be written */
	private static String expectNativeParentExposedSlotFailure() throws Exception{
		IllegalArgumentException failure=null;
		try{
			apply(apply(baseTree(), FIRST_OVERLAY), SECOND_NATIVE_PARENT);
		}catch(IllegalArgumentException e){
			failure=e;
		}
		check(failure!=null,
				"baseline should expose old private slots before populating level data");
		check("TaxNode 2000000000 is stored at index 562".equals(failure.getMessage()),
				"baseline exposed-slot failure had the wrong message");
		check(hasStackFrame(failure, "tax.TaxTree", "populateLevelData"),
				"baseline exposed-slot failure did not originate in TaxTree.populateLevelData");
		return failure.getMessage();
	}

	/** Checks the baseline failure when a second overlay parents to a prior fake ID.
	 * @return Exact exception message from the baseline failure
	 * @throws Exception If a fixture cannot be written */
	private static String expectUnknownPriorAuxParent() throws Exception{
		String message=null;
		try{
			apply(apply(baseTree(), FIRST_OVERLAY), SECOND_AUX_PARENT);
		}catch(IllegalArgumentException e){
			message=e.getMessage();
		}
		check("Auxiliary TaxID 2000000010 has unknown parent 2000000000".equals(message),
				"baseline should reject a second overlay parented to a prior fake ID");
		return message;
	}

	/** Checks fixed rejection of a second overlay parented to an old private slot.
	 * @param base Previously loaded overlay tree
	 * @return Exact exception message from the fixed failure
	 * @throws Exception If a fixture cannot be written */
	private static String expectPrivateSlotParentFailure(TaxTree base) throws Exception{
		String message=null;
		try{
			apply(base, SECOND_PRIVATE_SLOT_PARENT);
		}catch(IllegalArgumentException e){
			message=e.getMessage();
		}
		check("Auxiliary TaxID 2000000010 has unknown parent 562".equals(message),
				"fixed path should reject an old private slot as a parent ID");
		return message;
	}

	/** Rejects duplicate durable IDs in one overlay.
	 * @throws Exception If the fixture cannot be written */
	private static void testDuplicate() throws Exception{
		expectFailure("fake_taxid\tparent_taxid\trank\tname\n"+
				"2000000000\t2\tphylum\tFakeota\n"+
				"2000000000\t2\tphylum\tOtherota\n");
	}

	/** Rejects an auxiliary parent absent from both base and overlay.
	 * @throws Exception If the fixture cannot be written */
	private static void testMissingParent() throws Exception{
		expectFailure("fake_taxid\tparent_taxid\trank\tname\n"+
				"2000000000\t2000000999\tphylum\tFakeota\n");
	}

	/** Rejects IDs outside the durable auxiliary-ID range.
	 * @throws Exception If the fixture cannot be written */
	private static void testLowFakeID() throws Exception{
		expectFailure("fake_taxid\tparent_taxid\trank\tname\n"+
				"6000000\t2\tphylum\tFakeota\n");
	}

	/** Checks that loading an invalid overlay fails rather than accepting it.
	 * @param body Literal invalid TSV content
	 * @throws Exception If the fixture cannot be written */
	private static void expectFailure(String body) throws Exception{
		boolean ok=false;
		try{
			load(body);
		}catch(RuntimeException e){
			ok=true;
		}
		check(ok, "auxiliary tree load should have failed");
	}

	/** Checks rejection of an out-of-range integer header.
	 * @param tree Loaded fixture tree
	 * @param s Invalid numeric header */
	private static void expectParseFailure(TaxTree tree, String s){
		boolean ok=false;
		try{
			tree.parseNodeFromHeader(s, true);
		}catch(RuntimeException e){
			ok=true;
		}
		check(ok, s+" should have failed during TaxID parsing");
	}

	/** Checks whether an exception stack contains a named method origin.
	 * @param t Throwable to inspect
	 * @param className Expected class name
	 * @param methodName Expected method name
	 * @return True if the stack trace contains that frame */
	private static boolean hasStackFrame(Throwable t, String className, String methodName){
		for(StackTraceElement ste : t.getStackTrace()){
			if(ste.getClassName().equals(className) && ste.getMethodName().equals(methodName)){
				return true;
			}
		}
		return false;
	}

	/** Parses the BBTools-style expectation flag.
	 * @param args Command-line arguments
	 * @return baseline or fixed */
	private static String parseExpectation(String[] args){
		if(args.length<1){return "fixed";}
		check(args.length==1, "unknown extra argument "+(args.length>1 ? args[1] : ""));
		String arg=args[0];
		String prefix="expect=";
		check(arg.regionMatches(true, 0, prefix, 0, prefix.length()), "unknown argument "+arg);
		String expectation=arg.substring(prefix.length());
		check("baseline".equals(expectation) || "fixed".equals(expectation),
				"expectation must be baseline or fixed, not "+expectation);
		return expectation;
	}

	/** Writes a process-lifetime temporary overlay and applies it to a fresh base.
	 * @param body Literal TSV content
	 * @return Independently constructed augmented tree
	 * @throws Exception If the fixture cannot be written or loaded */
	private static TaxTree load(String body) throws Exception{
		return apply(baseTree(), body);
	}

	/** Writes a process-lifetime temporary overlay and applies it to an existing tree.
	 * @param tree Base tree
	 * @param body Literal TSV content
	 * @return Augmented tree
	 * @throws Exception If the fixture cannot be written or loaded */
	private static TaxTree apply(TaxTree tree, String body) throws Exception{
		File f=File.createTempFile("taxtree_aux_test_", ".tsv");
		f.deleteOnExit();
		PrintWriter pw=new PrintWriter(f);
		pw.print(body);
		pw.close();
		return TaxTreeAux.apply(tree, f.getAbsolutePath());
	}

	/** Creates Life, Bacteria and Escherichia with literal native IDs and child counts.
	 * @return Fresh base tree with room for native IDs through561 */
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

	/** Keeps fixture checks active under both assertion modes.
	 * @param b Required condition
	 * @param message Failure explanation */
	private static void check(boolean b, String message){
		if(!b){throw new AssertionError(message);}
	}

}
