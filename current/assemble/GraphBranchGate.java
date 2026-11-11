package assemble;

/**
 * Vetoes ordinary direct implicit joins at unresolved outward branch ends.
 * BubblePopper bypasses this gate for cross-K merges; explicit-overlap edges
 * are accepted without counting. An implicit join is rejected if either
 * consumed endpoint has F_BRANCH or D_BRANCH, even when its path is shorter
 * than K. TadpoleMulti can install it for final graph simplification;
 * MultiKBranchRescue can also install it for initial cleanup.
 * @author Fischl
 */
public final class GraphBranchGate{

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Checks consumed ends before the caller flips or merges the destination.
	 * Only this gate's counters can change; contigs and reciprocal graph edges
	 * remain untouched by this check.
	 */
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

	/** Returns whether {@code code} names forward branching or branching in both directions. */
	private static boolean outwardBranch(final int code){
		assert(code>=0 && code<Tadpole.MAX_CODE) : "Endpoint codes index ShaveObject.codeStrings: "+code;
		return code==Tadpole.F_BRANCH || code==Tadpole.D_BRANCH;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Candidate implicit edges that reached this gate, excluding explicit overlaps. */
	public long evaluations;

	/** Candidate implicit edges vetoed for consuming an unresolved outward branch. */
	public long rejected;
}
