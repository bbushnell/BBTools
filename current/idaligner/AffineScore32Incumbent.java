package idaligner;

/** Allocation-free construction of achievable incumbents for affine scoring. */
public final class AffineScore32Incumbent {

	/** Returns the score of one contiguous MSA-like query placement.
	 * Rejects lengths whose worst-case positive or negative sum could overflow int.
	 */
	public static int msaLikeUngapped(final byte[] query, final byte[] ref,
			final int refStart) {
		if(query==null) {throw new NullPointerException("A query is required.");}
		if(ref==null) {throw new NullPointerException("A reference is required.");}
		if(refStart<0 || refStart>ref.length-query.length) {
			throw new IndexOutOfBoundsException("Invalid ungapped placement "+
					refStart+".."+(refStart+query.length-1)+" for length "+ref.length);
		}
		validateLength(query.length);
		int score=0;
		for(int i=0; i<query.length; i++) {
			final int q=CrossCutScore32.baseCode(query[i]);
			final int r=CrossCutScore32.baseCode(ref[refStart+i]);
			score+=q==r && q<=8 ? MATCH : (q>8 || r>8 ? AMBIGUITY : SUBSTITUTION);
		}
		return score;
	}

	static void validateLength(final int queryLength) {
		if(queryLength<0 || (long)queryLength*MATCH>Integer.MAX_VALUE ||
				(long)queryLength*SUBSTITUTION<Integer.MIN_VALUE) {
			throw new IllegalArgumentException("Ungapped query length exceeds safe int score sums: "+queryLength);
		}
	}

	private AffineScore32Incumbent() {}
	private static final int MATCH=100;
	private static final int SUBSTITUTION=-127;
	private static final int AMBIGUITY=0;
}
