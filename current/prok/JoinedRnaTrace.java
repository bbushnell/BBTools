package prok;

/** Structural checks for an insertion hypothesis, independent of benchmark tools.
 * Consensus-query reference deletions represent the proposed genomic insertion.
 * @author Brian Bushnell, Raiden */
final class JoinedRnaTrace {
	static int longestDeletion(byte[] trace){
		if(trace==null){throw new IllegalArgumentException("Joined validity requires the retained expanded MSA trace");}
		int longest=0,run=0;for(byte op:trace){if(op=='D'){run++;longest=Math.max(longest,run);}else{run=0;}}return longest;
	}
	static boolean valid(Euk18sPacBioAligner.Result r,int windowLength){
		if(r==null || !r.validBounds || r.leftOverhang!=0 || r.rightOverhang!=0){return false;}
		final int span=r.rStop-r.rStart+1,longest=longestDeletion(r.matchString);
		assert(longest<=span):"A real reference-deletion run cannot exceed the physical alignment span";
		return longest>=256 && span>0 && span<=windowLength && span-longest>=916;
	}
}
