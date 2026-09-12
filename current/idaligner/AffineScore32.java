package idaligner;

import java.util.Arrays;

/**
 * Configurable score-only affine aligner using six reusable int arrays.
 * The full query is aligned while reference prefix and suffix are free.
 * One instance is intended per worker thread and is not thread-safe.
 *
 * <p>Scoring rejects dimensions/parameters whose conservative path bounds
 * cannot fit strictly above the unreachable sentinel and at most
 * {@link Integer#MAX_VALUE}. For a nonempty query/reference window, the lower
 * bound is {@code queryLength*min(0,substitution,ambiguity,insertionOpen,
 * insertionExtend) + windowLength*min(deletionOpen,deletionExtend)};
 * the upper bound is {@code queryLength*max(match,ambiguity)}. The lower bound
 * must exceed {@code Integer.MIN_VALUE/4}. Empty-reference calls instead check
 * their exact all-insertion score. The query/window length sum must be at most
 * {@code Integer.MAX_VALUE-64} to protect shared diagonal/word indexing.
 * These bounds can reject inputs whose optimum alone would fit.</p>
 */
public final class AffineScore32 {

	public AffineScore32(final int match_, final int substitution_, final int ambiguity_,
			final int insertionOpen_, final int insertionExtend_,
			final int deletionOpen_, final int deletionExtend_) {
		if(match_<=0 || substitution_>0 || insertionOpen_>0 || insertionExtend_>0 ||
				deletionOpen_>0 || deletionExtend_>0) {
			throw new IllegalArgumentException("Scores must use positive matches and nonpositive penalties.");
		}
		match=match_;
		substitution=substitution_;
		ambiguity=ambiguity_;
		insertionOpen=insertionOpen_;
		insertionExtend=insertionExtend_;
		deletionOpen=deletionOpen_;
		deletionExtend=deletionExtend_;
	}

	/**
	 * Simplified affine approximation using MSA's steady-state match reward,
	 * substitution cost, and first insertion/deletion open and extension tiers.
	 */
	public static AffineScore32 msaLike() {
		return new AffineScore32(100, -127, 0, -395, -39, -472, -33);
	}

	/** Returns the best score over the final query row. */
	public int score(final byte[] query, final byte[] ref) {
		if(query==null) {throw new NullPointerException("A query is required.");}
		if(ref==null) {throw new NullPointerException("A reference is required.");}
		return ref.length==0 ? scoreEmptyReference(query.length) :
			score(query, ref, 0, ref.length-1);
	}

	/** Returns the best score in an inclusive, nonempty reference window. */
	public int score(final byte[] query, final byte[] ref, final int refStart,
			final int refEnd) {
		if(query==null) {throw new NullPointerException("A query is required.");}
		validateRange(ref, refStart, refEnd);
		final int qLen=query.length;
		final int rLen=refEnd-refStart+1;
		validateScoreRange(qLen, rLen, match, substitution, ambiguity,
				insertionOpen, insertionExtend, deletionOpen, deletionExtend);
		if(qLen==0) {return 0;}
		ensureCapacity(rLen+1);

		int[] prevM=matchA;
		int[] prevI=insertionA;
		int[] prevD=deletionA;
		int[] currM=matchB;
		int[] currI=insertionB;
		int[] currD=deletionB;
		Arrays.fill(prevM, 0, rLen+1, 0);
		Arrays.fill(prevI, 0, rLen+1, BAD);
		Arrays.fill(prevD, 0, rLen+1, 0);

		for(int row=1; row<=qLen; row++) {
			currM[0]=BAD;
			currD[0]=BAD;
			currI[0]=(row==1 ? insertionOpen : add(prevI[0], insertionExtend));
			for(int col=1; col<=rLen; col++) {
				final int baseScore=baseScore(query[row-1], ref[refStart+col-1]);
				currM[col]=add(max3(prevM[col-1], prevI[col-1], prevD[col-1]), baseScore);
				currI[col]=max3(add(prevM[col], insertionOpen),
						add(prevI[col], insertionExtend), add(prevD[col], insertionOpen));
				currD[col]=max3(add(currM[col-1], deletionOpen),
						add(currD[col-1], deletionExtend), add(currI[col-1], deletionOpen));
			}

			int[] swap=prevM; prevM=currM; currM=swap;
			swap=prevI; prevI=currI; currI=swap;
			swap=prevD; prevD=currD; currD=swap;
		}

		int best=BAD;
		for(int col=1; col<=rLen; col++) {
			best=Math.max(best, max3(prevM[col], prevI[col], prevD[col]));
		}
		return best;
	}

	/** Returns the all-insertion score used when the reference is empty. */
	int scoreEmptyReference(final int queryLength) {
		validateScoreRange(queryLength, 0, match, substitution, ambiguity,
				insertionOpen, insertionExtend, deletionOpen, deletionExtend);
		if(queryLength<=0) {return 0;}
		return insertionOpen+(queryLength-1)*insertionExtend;
	}

	int baseScore(final byte queryBase, final byte refBase) {
		final int q=CrossCutScore32.baseCode(queryBase);
		final int r=CrossCutScore32.baseCode(refBase);
		return q==r && q<=8 ? match : (q>8 || r>8 ? ambiguity : substitution);
	}

	private void ensureCapacity(final int length) {
		if(matchA.length>=length) {return;}
		matchA=new int[length];
		matchB=new int[length];
		insertionA=new int[length];
		insertionB=new int[length];
		deletionA=new int[length];
		deletionB=new int[length];
	}

	private static int add(final int score, final int delta) {
		return score<=BAD ? BAD : score+delta;
	}

	private static int max3(final int a, final int b, final int c) {
		return Math.max(a, Math.max(b, c));
	}

	private static void validateRange(final byte[] ref, final int refStart, final int refEnd) {
		if(ref==null) {throw new NullPointerException("A reference is required.");}
		if(refStart<0 || refEnd<refStart || refEnd>=ref.length) {
			throw new IndexOutOfBoundsException("Invalid inclusive reference range "+
					refStart+".."+refEnd+" for length "+ref.length);
		}
	}

	/** Checks the shared numeric contract once, outside all DP loops. */
	static void validateScoreRange(final int qLen, final int rLen,
			final int match, final int substitution, final int ambiguity,
			final int insertionOpen, final int insertionExtend,
			final int deletionOpen, final int deletionExtend) {
		if(qLen<0 || rLen<0 || (long)qLen+rLen>Integer.MAX_VALUE-64L) {
			throw new IllegalArgumentException("Affine dimensions exceed safe diagonal/word indexing: query="+
					qLen+", window="+rLen);
		}
		if(qLen==0) {return;}
		final long lower, upper;
		if(rLen==0) {
			lower=upper=insertionOpen+(qLen-1L)*insertionExtend;
		}else {
			final int queryMin=Math.min(0, Math.min(Math.min(substitution, ambiguity),
					Math.min(insertionOpen, insertionExtend)));
			// Each transition consumes a query base or is a reference deletion.
			// Thus these bounds cover every prefix and every candidate addition,
			// including paths that will lose a max/pruning comparison.
			lower=(long)qLen*queryMin+(long)rLen*Math.min(deletionOpen, deletionExtend);
			upper=(long)qLen*Math.max(match, ambiguity);
		}
		if(lower<=BAD || upper>Integer.MAX_VALUE) {
			throw new IllegalArgumentException("Affine score bounds outside safe int domain: query="+
					qLen+", window="+rLen+", lower="+lower+", upper="+upper+
					"; require lower>"+BAD+" and upper<="+Integer.MAX_VALUE);
		}
		assert(lower>BAD && upper<=Integer.MAX_VALUE):
				"Reachable DP states must stay above BAD and int additions must not wrap.";
	}

	private int[] matchA=new int[0];
	private int[] matchB=new int[0];
	private int[] insertionA=new int[0];
	private int[] insertionB=new int[0];
	private int[] deletionA=new int[0];
	private int[] deletionB=new int[0];

	final int match;
	final int substitution;
	final int ambiguity;
	final int insertionOpen;
	final int insertionExtend;
	final int deletionOpen;
	final int deletionExtend;
	// validateScoreRange keeps reachable states distinct from this sentinel.
	static final int BAD=Integer.MIN_VALUE/4;
}
