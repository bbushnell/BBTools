package test.tax;

import tax.IDNode;

/** Finite baseline for IDNode all-zero self-merge behavior.
 * @author Jean */
public final class IDNodeAllZeroSimilarityTest{

	/** Runs baseline or fixed expectations for zero-only similarity rows.
	 * @param args Optional expect=baseline or expect=fixed */
	public static void main(String[] args){
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		checkPositiveSimilarityControl();
		checkAllZeroPair(baseline);
		checkAllZeroTriplet(baseline);
		checkDisconnectedZeroComponents(baseline);
		System.out.println("ID_NODE_ALL_ZERO_SIMILARITY_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private IDNodeAllZeroSimilarityTest(){}

	/** Verifies an ordinary positive-identity two-leaf cluster is intact. */
	private static void checkPositiveSimilarityControl(){
		final IDNode root=IDNode.makeTree(new IDNode[]{
				new IDNode(new double[0], 0, "positive_a"),
				new IDNode(new double[]{0.75}, 1, "positive_b")});
		checkTwoLeafTree(root, "positive_a", "positive_b");
		checkPositiveNewick(root.toNewick().toString());
	}

	/** Verifies two all-zero rows hit the baseline self-merge or form a fixed tree.
	 * @param baseline Whether to assert the baseline or fixed expectation */
	private static void checkAllZeroPair(final boolean baseline){
		IDNode root=null;
		Throwable failure=null;
		try{
			root=IDNode.makeTree(new IDNode[]{
					new IDNode(new double[0], 0, "zero_a"),
					new IDNode(new double[]{0}, 1, "zero_b")});
		}catch(Throwable t){
			failure=t;
		}

		if(baseline){
			checkBaselineSelfMerge(failure, "two all-zero leaves");
		}else{
			check(failure==null, "fixed two-zero rows must form a tree: "+failure);
			checkTwoLeafTree(root, "zero_a", "zero_b");
		}
	}

	/** Verifies three all-zero rows cannot strand leaf zero as the chosen partner.
	 * @param baseline Whether to assert the baseline or fixed expectation */
	private static void checkAllZeroTriplet(final boolean baseline){
		IDNode root=null;
		Throwable failure=null;
		try{
			root=IDNode.makeTree(new IDNode[]{
					new IDNode(new double[0], 0, "zero_a"),
					new IDNode(new double[]{0}, 1, "zero_b"),
					new IDNode(new double[]{0, 0}, 2, "zero_c")});
		}catch(Throwable t){
			failure=t;
		}

		if(baseline){
			checkBaselineSelfMerge(failure, "three all-zero leaves");
		}else{
			check(failure==null, "fixed three-zero rows must form a tree: "+failure);
			checkTree(root, "zero_a", "zero_b", "zero_c");
		}
	}

	/** Verifies disconnected zero-similarity components merge once each.
	 * @param baseline Whether to assert the baseline or fixed expectation */
	private static void checkDisconnectedZeroComponents(final boolean baseline){
		IDNode root=null;
		Throwable failure=null;
		try{
			root=IDNode.makeTree(new IDNode[]{
					new IDNode(new double[0], 0, "component_a"),
					new IDNode(new double[]{0.90}, 1, "component_b"),
					new IDNode(new double[]{0, 0}, 2, "component_c"),
					new IDNode(new double[]{0, 0, 0.80}, 3, "component_d")});
		}catch(Throwable t){
			failure=t;
		}

		if(baseline){
			if(assertionsEnabled()){
				check(failure instanceof AssertionError,
						"disconnected -ea baseline must trip a self-merge assertion: "+failure);
				check(failure.getStackTrace().length>0 &&
						failure.getStackTrace()[0].getClassName().equals("tax.IDNode") &&
						failure.getStackTrace()[0].getMethodName().equals("<init>"),
						"disconnected -ea baseline failed outside IDNode merge construction: "+
						failure);
			}else{
				check(failure==null,
						"disconnected -da baseline should duplicate A/B without throwing: "+
						failure);
				check(countLabel(root, "component_a")==2,
						"disconnected -da baseline should duplicate component_a twice: "+
						root.toNewick());
				check(countLabel(root, "component_b")==2,
						"disconnected -da baseline should duplicate component_b twice: "+
						root.toNewick());
				check(countLabel(root, "component_c")==1,
						"disconnected -da baseline should retain component_c once: "+
						root.toNewick());
				check(countLabel(root, "component_d")==1,
						"disconnected -da baseline should retain component_d once: "+
						root.toNewick());
				check(countLeaves(root)==6,
						"disconnected -da baseline should contain six leaves: "+
						root.toNewick());
			}
		}else{
			check(failure==null,
					"fixed disconnected zero components must form a tree: "+failure);
			checkTree(root, "component_a", "component_b", "component_c", "component_d");
		}
	}

	/** Verifies a baseline self-merge fails in the merge constructor.
	 * @param failure Captured failure
	 * @param scenario Diagnostic scenario name */
	private static void checkBaselineSelfMerge(final Throwable failure,
			final String scenario){
		final boolean assertions=assertionsEnabled();
		check(failure!=null, scenario+" baseline must fail on a self-merge");
		if(assertions){
			check(failure instanceof AssertionError,
					scenario+" -ea baseline must trip self-merge: "+failure);
		}else{
			check(failure instanceof ArrayIndexOutOfBoundsException,
					scenario+" -da baseline must overrun the self-merged row: "+
					failure);
		}
		check(failure.getStackTrace().length>0 &&
				failure.getStackTrace()[0].getClassName().equals("tax.IDNode") &&
				failure.getStackTrace()[0].getMethodName().equals("<init>"),
				scenario+" baseline failed outside IDNode merge construction: "+
				failure);
	}

	/** Verifies a cluster tree retains each label exactly once.
	 * @param root Cluster root
	 * @param labels Expected leaf labels in row order */
	private static void checkTree(final IDNode root, final String... labels){
		check(root!=null, "Root is null");
		check(root.parent==null, "Root should not have a parent");
		check(root.bs.cardinality()==labels.length,
				"Root membership should contain exactly "+labels.length+
				" leaves: "+root.bs);
		for(int i=0; i<labels.length; i++){
			check(root.bs.get(i), "Root membership should contain leaf "+i+
					": "+root.bs);
			check(countLabel(root, labels[i])==1,
					"Tree must contain exactly one "+labels[i]+
					" leaf: "+root.toNewick());
		}
		check(countLeaves(root)==labels.length,
				"Tree should contain exactly "+labels.length+
				" leaves: "+root.toNewick());
		checkSubtree(root, null);
	}

	/** Verifies parent links and distinct children recursively.
	 * @param node Cluster node
	 * @param parent Expected parent */
	private static void checkSubtree(final IDNode node, final IDNode parent){
		check(node!=null, "Tree contains a null node");
		check(node.parent==parent, "Wrong parent for "+node);
		check((node.left==null)==(node.right==null),
				"Node should have zero or two children: "+node);
		if(node.left!=null){
			check(node.left!=node.right, "Sibling links should be distinct: "+node);
			checkSubtree(node.left, node);
			checkSubtree(node.right, node);
		}
	}

	/** Verifies a two-leaf tree retains each label exactly once.
	 * Does not constrain left/right topology.
	 * @param root Cluster root
	 * @param a First expected leaf label
	 * @param b Second expected leaf label */
	private static void checkTwoLeafTree(final IDNode root, final String a, final String b){
		check(root!=null, "Root is null");
		check(root.parent==null, "Root should not have a parent");
		check(root.left!=null && root.right!=null,
				"Root should have two children: "+root);
		check(root.left!=root.right, "Root children should be distinct nodes");
		check(root.left.left==null && root.left.right==null &&
				root.right.left==null && root.right.right==null,
				"Root children should both be leaves: "+root.toNewick());
		checkTree(root, a, b);
	}

	/** Verifies the positive control's exact branch lengths in either topology.
	 * @param newick Newick rendering without final semicolon */
	private static void checkPositiveNewick(final String newick){
		check(newick.equals("(positive_a:0.2500,positive_b:0.2500)") ||
				newick.equals("(positive_b:0.2500,positive_a:0.2500)"),
				"Positive control Newick should contain literal 0.2500 branches: "+
				newick);
	}

	/** Counts matching labels in a cluster tree.
	 * @param node Cluster node, or null
	 * @param label Label to count
	 * @return Number of nodes with the requested label */
	private static int countLabel(final IDNode node, final String label){
		if(node==null){return 0;}
		int count=label.equals(node.name) ? 1 : 0;
		count+=countLabel(node.left, label);
		count+=countLabel(node.right, label);
		return count;
	}

	/** Counts leaves in a cluster tree.
	 * @param node Cluster node, or null
	 * @return Number of leaves */
	private static int countLeaves(final IDNode node){
		if(node==null){return 0;}
		if(node.left==null){return 1;}
		return countLeaves(node.left)+countLeaves(node.right);
	}

	/** Returns true when assertions are enabled for this class.
	 * @return Assertion-enabled status */
	private static boolean assertionsEnabled(){
		boolean assertions=false;
		assert(assertions=true);
		return assertions;
	}

	/** Keeps fixture checks active under -da.
	 * @param condition Required condition
	 * @param message Failure message */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

}
