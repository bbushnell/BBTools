package shared;

import java.util.BitSet;

/**
 * Marks tightly bounded low-depth k-mer runs without attempting correction.
 * The caller supplies the tool-specific relative-depth error predicate.
 */
public final class ErrorRunMarker {

	private ErrorRunMarker() {}

	/** Tool-specific implementation of the existing high-vs-low error test. */
	public interface ErrorPredicate {
		boolean isError(int high, int low);
	}

	/** Per-read marking and skip diagnostics. */
	public static final class Result {
		public int marked;
		public int exactKRuns;
		public int exactKMinus1Runs;
		public int edgeRunsMarked;
		public int skippedEndRuns;
		public int skippedMergedRuns;
	}

	/**
	 * Marks one base for isolated runs of exactly k or k-1 low-count kmers.
	 * A valid run has two consecutive high-depth kmers on each side.  The low
	 * threshold is max(1, flankDepth/ratio), and both flank depths must pass
	 * the caller's existing error predicate.  Runs longer than k are treated as
	 * merged and skipped; runs lacking a two-kmer flank are treated as end runs.
	 */
	public static Result mark(final byte[] bases, final byte[] quals, final int[] counts,
			final int countSize, final int k, final float ratio, final ErrorPredicate predicate) {
		final Result result=new Result();
		if(bases==null || counts==null || countSize<1 || k<2 || bases.length<k){return result;}
		final float safeRatio=Tools.max(1f, ratio);
		final BitSet marked=new BitSet(bases.length);

		// First identify longer isolated runs so they are explicitly reported as merged.
		for(int s=2; s+k+1<countSize; s++){
			final int left=Tools.min(counts[s-2], counts[s-1]);
			final int threshold=Tools.max(1, (int)(left/safeRatio));
			int e=s;
			while(e<countSize && counts[e]<=threshold && predicate.isError(left, counts[e])){e++;}
			final int end=e-1;
			if(end-s+1>k && end+2<countSize && qualifies(counts, countSize, s, end, safeRatio, predicate, true)){
				result.skippedMergedRuns++;
				s=end;
			}
		}

		// Mark only exact-length runs.  K is considered before K-1 so a valid
		// substitution/insert signature wins over an overlapping deletion signature.
		for(final int runLen : new int[] {k, k-1}){
			for(int s=0, end=s+runLen-1; end<countSize; s++, end++){
				if(!qualifies(counts, countSize, s, end, safeRatio, predicate, true)){continue;}
				final int pos=(runLen==k ? s+k-1 : s+k-2);
				if(pos<k-1 || pos>bases.length-k){continue;}
				if(marked.get(pos)){continue;}
				if(bases[pos]!='N'){
					bases[pos]='N';
					if(quals!=null && pos<quals.length){quals[pos]=0;}
					marked.set(pos);
					result.marked++;
					if(runLen==k){result.exactKRuns++;}
					else{result.exactKMinus1Runs++;}
				}
			}
		}

		// At a read edge, an eligible run need not be exactly K long.  Mark the
		// inferred error base for any run of length 1..K.  For a left-edge run
		// starting at k-mer 0, that base is the rightmost low k-mer's right edge.
		boolean leftEdgeHandled=false;
		for(int runLen=1; !leftEdgeHandled && runLen<=k && runLen<countSize; runLen++){
			final int end=runLen-1;
			if(qualifies(counts, countSize, 0, end, safeRatio, predicate, false)){
				final int pos=end;
				if(pos>=0 && pos<bases.length && bases[pos]!='N'){
					bases[pos]='N';
					if(quals!=null && pos<quals.length){quals[pos]=0;}
					result.marked++;
					result.edgeRunsMarked++;
				}
				leftEdgeHandled=true;
			}
		}
		boolean rightEdgeHandled=false;
		for(int runLen=1; !rightEdgeHandled && runLen<=k && runLen<countSize; runLen++){
			final int s=countSize-runLen;
			final int end=countSize-1;
			if(qualifies(counts, countSize, s, end, safeRatio, predicate, false)){
				final int pos=s+k-1;
				if(pos>=0 && pos<bases.length && bases[pos]!='N'){
					bases[pos]='N';
					if(quals!=null && pos<quals.length){quals[pos]=0;}
					result.marked++;
					result.edgeRunsMarked++;
				}
				rightEdgeHandled=true;
			}
		}

		// Report exact-looking edge runs that exceed the permitted edge length or
		// otherwise lack a valid one-sided high-depth flank.
		final BitSet endSeen=new BitSet(countSize);
		for(final int runLen : new int[] {k, k-1}){
			for(int s=0, end=s+runLen-1; end<countSize; s++, end++){
				if(s>=2 && end+2<countSize){continue;}
				if((s==0 && leftEdgeHandled) || (end==countSize-1 && rightEdgeHandled)){continue;}
				if(endSeen.get(s) || !qualifies(counts, countSize, s, end, safeRatio, predicate, false)){continue;}
				result.skippedEndRuns++;
				endSeen.set(s);
			}
		}
		return result;
	}

	private static boolean qualifies(final int[] counts, final int countSize, final int s, final int e,
			final float ratio, final ErrorPredicate predicate, final boolean requireBoth){
		if(s<0 || e>=countSize || e-s+1<1){return false;}
		final boolean leftPresent=s>=2, rightPresent=e+2<countSize;
		if(requireBoth && (!leftPresent || !rightPresent)){return false;}
		if(!leftPresent && !rightPresent){return false;}

		int left=Integer.MAX_VALUE, right=Integer.MAX_VALUE, lowMax=0;
		if(leftPresent){left=Tools.min(counts[s-2], counts[s-1]);}
		if(rightPresent){right=Tools.min(counts[e+1], counts[e+2]);}
		for(int i=s; i<=e; i++){lowMax=Tools.max(lowMax, counts[i]);}
		final int flank=Tools.min(left, right);
		final int threshold=Tools.max(1, (int)(flank/ratio));
		final int highRequired=Tools.max(1, (int)Math.ceil(lowMax*ratio));
		if(leftPresent && left<=highRequired){return false;}
		if(rightPresent && right<=highRequired){return false;}
		for(int i=s; i<=e; i++){
			final int low=counts[i];
			if(low>threshold){return false;}
			if(leftPresent && !predicate.isError(left, low)){return false;}
			if(rightPresent && !predicate.isError(right, low)){return false;}
		}
		return true;
	}
}
