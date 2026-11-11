package assemble;

import dna.AminoAcid;

/**
 * Experimental veto for implicit graph walks that conflict with an exact overlap.
 * Explicit-overlap and cross-K joins keep their own policies. Ordinary implicit
 * paths are rejected when a bounded exact parent overlap of at least K implies
 * different spacing. An overlap is evidence of ambiguity, not permission to
 * replace the graph walk. The final-graph caller supplies assemble-K minus one
 * as the configured search ceiling.
 * @author Fischl
 */
public final class GraphOverlapGate{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Stores a positive caller-supplied overlap search ceiling. */
	public GraphOverlapGate(final int maxOverlap_){
		if(maxOverlap_<1){throw new IllegalArgumentException("Overlap search requires a positive ceiling.");}
		maxOverlap=maxOverlap_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Checks consumed ends before the caller flips or merges the destination.
	 * The search ceiling is also bounded below both retained parent lengths.
	 * Explicit-overlap edges and empty search ranges are accepted without counting.
	 * No contig or edge is mutated; only this gate's counters can change.
	 */
	public boolean accept(final Contig left, final Contig right, final Edge edge, final int k){
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

	/** Compares an oriented left suffix and right prefix without allocating. */
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

	/** Largest exact parent overlap allowed by the caller. */
	final int maxOverlap;

	/** Eligible implicit walks entering the search. */
	long evaluations;

	/** Eligible implicit walks vetoed by a conflicting exact overlap. */
	public long rejected;
}
