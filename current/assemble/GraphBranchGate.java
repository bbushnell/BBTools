package assemble;

/**
 * Experimental final-graph veto for implicit paths at outward branch ends.
 * A single completed contig edge does not rule out other unresolved kmer paths.
 * This policy neither alters exact-overlap decisions nor consults another table.
 * @author Fischl
 */
final class GraphBranchGate {

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks the consumed ends before destination flipping; rejection never mutates graph state. */
	boolean accept(final Contig left, final Contig right, final Edge edge, final int k){
		assert(left!=null && right!=null && edge!=null && k>0 && edge.sourceRight()) :
				"BubblePopper.expandRightSimple supplies a rightward edge before destination flipping.";
		if(edge.overlap>0){return true;}
		// A path shorter than K still chooses a repeat exit. Citricoccus supplied
		// a false inversion through 19 steps at K32, despite the resulting 13bp overlap.
		assert(edge.length>=0 && (edge.bases==null ? edge.length==0 : edge.bases.length==edge.length)) :
				"Implicit path length must match the bytes consumed by BubblePopper.merge(), including an empty path.";
		evaluations++;
		final int leftCode=left.rightCode, rightCode=edge.destRight() ? right.rightCode : right.leftCode;
		if(!outwardBranch(leftCode) && !outwardBranch(rightCode)){return true;}
		rejected++;
		if(BubblePopper.verbose){
			System.err.println("Graph branch rejected: left="+left.id+", right="+right.id+
					", k="+k+", path="+edge.length+", overlap="+(k-edge.length)+
					", codes="+leftCode+","+rightCode+".");
		}
		return false;
	}

	/** Matches the outward-branch states produced by Tadpole.classifyGraphEnd. */
	private static boolean outwardBranch(final int code){
		assert(code>=0 && code<Tadpole.MAX_CODE) : "Endpoint codes index ShaveObject.codeStrings: "+code;
		return code==Tadpole.F_BRANCH || code==Tadpole.D_BRANCH;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	long evaluations, rejected;
}
