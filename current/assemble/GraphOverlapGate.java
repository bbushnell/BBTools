package assemble;

import dna.AminoAcid;

/**
 * Experimental veto for graph walks conflicting with an existing exact overlap.
 * An overlap is evidence of ambiguity, not permission to replace the graph walk.
 * The serial final-graph phase owns this policy and its counters.
 * @author Fischl
 */
final class GraphOverlapGate {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Bounds the search by the original assemble-K minus one. */
	GraphOverlapGate(final int maxOverlap_){
		if(maxOverlap_<1){throw new IllegalArgumentException("Overlap search requires a positive ceiling.");}
		maxOverlap=maxOverlap_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Checks merge-time sequences before BubblePopper flips the destination.
	 * Explicit overlap edges retain their own fusion policy. Ordinary walks are
	 * vetoed only when their spacing differs from an exact graph-K-or-longer overlap.
	 * No input, edge, or graph state is changed, even on rejection.
	 */
	boolean accept(final Contig left, final Contig right, final Edge edge, final int k){
		assert(left!=null && right!=null && edge!=null && k>0 && edge.sourceRight()) :
				"BubblePopper.expandRightSimple supplies a rightward edge before flipping its destination.";
		if(edge.overlap>0 || k>maxOverlap){return true;}
		assert(edge.sourceTrim>=0 && edge.destTrim>=0 && edge.sourceTrim<left.length()
				&& edge.destTrim<right.length()) : "Trim bounds must preserve both parents for merge() spelling.";
		final int leftEnd=left.length()-edge.sourceTrim;
		final int limit=Math.min(maxOverlap, Math.min(leftEnd-1, right.length()-edge.destTrim-1));
		if(limit<k){return true;}
		evaluations++;
		final int pathLength=edge.bases==null ? 0 : edge.bases.length;
		assert(pathLength==edge.length) : "Graph edge length must match the bytes consumed by BubblePopper.merge().";
		for(int overlap=limit; overlap>=k; overlap--){
			if((long)pathLength-k==-(long)overlap){continue;}// Same spacing is not a conflict.
			if(exact(left.bases, leftEnd-overlap, right.bases, edge.destRight(), edge.destTrim, overlap)){
				rejected++;
				if(BubblePopper.verbose){
					System.err.println("Graph overlap rejected: left="+left.id+", right="+right.id+
							", k="+k+", overlap="+overlap+", edge="+pathLength+".");
				}
				return false;
			}
		}
		return true;
	}

	/** Compares an oriented suffix/prefix without allocating a reverse complement. */
	private static boolean exact(final byte[] left, final int start, final byte[] right,
			final boolean reverse, final int rightStart, final int length){
		assert(start>=0 && start+length<=left.length && rightStart>=0 && rightStart+length<=right.length) :
				"Overlap comparison must stay inside both retained parental sequences.";
		for(int i=0; i<length; i++){
			final byte a=left[start+i];
			final byte raw=right[reverse ? right.length-1-rightStart-i : rightStart+i];
			if(!AminoAcid.isFullyDefined(a) || !AminoAcid.isFullyDefined(raw)){return false;}
			final byte b=reverse ? AminoAcid.baseToComplementExtended[raw] : raw;
			if(a!=b){return false;}
		}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	final int maxOverlap;
	long evaluations, rejected;
}
