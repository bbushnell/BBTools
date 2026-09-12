package idaligner;

import java.util.Arrays;

/** Exact incumbent-bounded affine score DP over reusable packed-word frontiers.
 * Scoring enforces the conservative numeric/dimension contract of {@link AffineScore32}.
 */
public final class AffineScore32ActiveWords {

	public AffineScore32ActiveWords(final int match_, final int substitution_,
			final int ambiguity_, final int insertionOpen_,
			final int insertionExtend_, final int deletionOpen_,
			final int deletionExtend_) {
		if(match_<=0 || substitution_>0 || ambiguity_>match_ ||
				insertionOpen_>0 || insertionExtend_>0 ||
				deletionOpen_>0 || deletionExtend_>0) {
			throw new IllegalArgumentException("Invalid affine scores for the upper bound.");
		}
		match=match_; substitution=substitution_; ambiguity=ambiguity_;
		insertionOpen=insertionOpen_; insertionExtend=insertionExtend_;
		deletionOpen=deletionOpen_; deletionExtend=deletionExtend_;
	}

	public static AffineScore32ActiveWords msaLike() {
		return new AffineScore32ActiveWords(100, -127, 0, -395, -39, -472, -33);
	}

	/**
	 * Returns the exact full score when {@code incumbent} is an achievable path
	 * score for these exact bytes. The intended caller supplies a same-input
	 * banded score.
	 */
	public int score(final byte[] query, final byte[] ref, final int incumbent) {
		if(query==null) {throw new NullPointerException("A query is required.");}
		if(ref==null) {throw new NullPointerException("A reference is required.");}
		if(ref.length>0) {return score(query, ref, 0, ref.length-1, incumbent);}
		AffineScore32.validateScoreRange(query.length, 0, match, substitution, ambiguity,
				insertionOpen, insertionExtend, deletionOpen, deletionExtend);
		resetStats(query.length, 0);
		final long ceiling=(long)query.length*match;
		if(incumbent>ceiling) {
			throw new IllegalArgumentException("Incumbent exceeds the score ceiling: "+incumbent);
		}
		if(query.length==0) {return 0;}
		final int exact=insertionOpen+(query.length-1)*insertionExtend;
		if(incumbent>exact) {throw new IllegalArgumentException("Incumbent is not a valid path score.");}
		return exact;
	}

	/** Returns the exact score in an inclusive, nonempty reference window. */
	public int score(final byte[] query, final byte[] ref, final int refStart,
			final int refEnd, final int incumbent) {
		if(query==null) {throw new NullPointerException("A query is required.");}
		validateRange(ref, refStart, refEnd);
		final int rLen=refEnd-refStart+1;
		AffineScore32.validateScoreRange(query.length, rLen, match, substitution, ambiguity,
				insertionOpen, insertionExtend, deletionOpen, deletionExtend);
		resetStats(query.length, rLen);
		final long ceiling=(long)query.length*match;
		if(incumbent>ceiling) {
			throw new IllegalArgumentException("Incumbent exceeds the score ceiling: "+incumbent);
		}
		if(query.length==0) {return 0;}
		if(incumbent==ceiling) {return incumbent;}
		return scoreNonempty(query, ref, refStart, rLen, incumbent);
	}

	private int scoreNonempty(final byte[] query, final byte[] ref,
			final int refStart, final int rLen, final int incumbent) {
		final int qLen=query.length, wordCount=(rLen+64)>>>6;
		ensureCapacity(rLen+1, wordCount);
		int[] prevM=matchA, prevI=insertionA, prevD=deletionA;
		int[] currM=matchB, currI=insertionB, currD=deletionB;
		Arrays.fill(prevM, 0, rLen+1, 0);
		Arrays.fill(prevI, 0, rLen+1, BAD);
		Arrays.fill(prevD, 0, rLen+1, 0);
		long[] prevActive=activeA, currActive=activeB;
		Arrays.fill(prevActive, 0, wordCount, -1L);
		maskTail(prevActive, wordCount, rLen+1);
		int best=incumbent;

		for(int row=1; row<=qLen; row++) {
			final int remaining=qLen-row;
			Arrays.fill(currActive, 0, wordCount, 0L);
			currM[0]=currD[0]=BAD;
			final int rawLeftI=get(prevActive, 0) ?
					(row==1 ? insertionOpen : add(prevI[0], insertionExtend)) : BAD;
			currI[0]=prune(rawLeftI, remaining, best);
			countState(rawLeftI, currI[0]);
			if(currI[0]!=BAD) {set(currActive, 0);}

			buildCandidates(prevActive, wordCount, rLen);
			int col=nextSetBit(candidateWords, 1, wordCount, rLen);
			while(col>=1) {
				final int rawM=get(prevActive, col-1) ?
						add(max3(prevM[col-1], prevI[col-1], prevD[col-1]),
								baseScore(query[row-1], ref[refStart+col-1])) : BAD;
				final int rawI=get(prevActive, col) ?
						max3(add(prevM[col], insertionOpen),
								add(prevI[col], insertionExtend),
								add(prevD[col], insertionOpen)) : BAD;
				currM[col]=prune(rawM, remaining, best);
				currI[col]=prune(rawI, remaining, best);
				final int rawD=get(currActive, col-1) ?
						max3(add(currM[col-1], deletionOpen),
								add(currD[col-1], deletionExtend),
								add(currI[col-1], deletionOpen)) : BAD;
				currD[col]=prune(rawD, remaining, best);
				lastProcessedCells++;
				countState(rawM, currM[col]);
				countState(rawI, currI[col]);
				countState(rawD, currD[col]);
				if(currM[col]!=BAD || currI[col]!=BAD || currD[col]!=BAD) {
					set(currActive, col);
				}
				if(row==qLen) {best=Math.max(best, max3(currM[col], currI[col], currD[col]));}
				final int nextCandidate=nextSetBit(candidateWords, col+1, wordCount, rLen);
				if(get(currActive, col) && col<rLen &&
						(nextCandidate<0 || col+1<nextCandidate)) {
					col++;
				}else {col=nextCandidate;}
			}

			int[] swap=prevM; prevM=currM; currM=swap;
			swap=prevI; prevI=currI; currI=swap;
			swap=prevD; prevD=currD; currD=swap;
			final long[] wordSwap=prevActive; prevActive=currActive; currActive=wordSwap;
		}
		return best;
	}

	private void buildCandidates(final long[] previous, final int wordCount,
			final int lastColumn) {
		long carry=0;
		for(int word=0; word<wordCount; word++) {
			final long bits=previous[word];
			lastFrontierBits+=Long.bitCount(bits);
			candidateWords[word]=bits | (bits<<1) | carry;
			carry=bits>>>63;
		}
		candidateWords[0]&=~1L;
		maskTail(candidateWords, wordCount, lastColumn+1);
	}

	private static int nextSetBit(final long[] words, final int from,
			final int wordCount, final int lastColumn) {
		if(from>lastColumn) {return -1;}
		int word=from>>>6;
		long bits=words[word] & (-1L<<(from&63));
		while(true) {
			if(bits!=0) {
				final int bit=(word<<6)+Long.numberOfTrailingZeros(bits);
				return bit<=lastColumn ? bit : -1;
			}
			if(++word>=wordCount) {return -1;}
			bits=words[word];
		}
	}

	private static void maskTail(final long[] words, final int wordCount,
			final int bitCount) {
		final int remainder=bitCount&63;
		if(remainder!=0) {words[wordCount-1]&=(1L<<remainder)-1L;}
	}

	private static boolean get(final long[] words, final int bit) {
		return (words[bit>>>6]&(1L<<(bit&63)))!=0;
	}

	private static void set(final long[] words, final int bit) {
		words[bit>>>6]|=1L<<(bit&63);
	}

	private int prune(final int score, final int remaining, final int best) {
		if(score==BAD) {return BAD;}
		return (long)score+(long)remaining*match>best ? score : BAD;
	}

	private void countState(final int raw, final int kept) {
		if(raw==BAD) {return;}
		if(kept==BAD) {lastPrunedStates++;}
		else {lastKeptStates++;}
	}

	private int baseScore(final byte queryBase, final byte refBase) {
		final int q=CrossCutScore32.baseCode(queryBase);
		final int r=CrossCutScore32.baseCode(refBase);
		return q==r && q<=8 ? match : (q>8 || r>8 ? ambiguity : substitution);
	}

	private void ensureCapacity(final int length, final int words) {
		if(matchA.length<length) {
			matchA=new int[length]; matchB=new int[length];
			insertionA=new int[length]; insertionB=new int[length];
			deletionA=new int[length]; deletionB=new int[length];
		}
		if(activeA.length<words) {
			activeA=new long[words]; activeB=new long[words]; candidateWords=new long[words];
		}
	}

	private void resetStats(final int qLen, final int rLen) {
		lastDenseCells=(long)qLen*rLen;
		lastProcessedCells=lastFrontierBits=lastKeptStates=lastPrunedStates=0;
	}

	private static void validateRange(final byte[] ref, final int refStart,
			final int refEnd) {
		if(ref==null) {throw new NullPointerException("A reference is required.");}
		if(refStart<0 || refEnd<refStart || refEnd>=ref.length) {
			throw new IndexOutOfBoundsException("Invalid inclusive reference range "+
					refStart+".."+refEnd+" for length "+ref.length);
		}
	}

	private static int add(final int score, final int delta) {
		return score==BAD ? BAD : score+delta;
	}

	private static int max3(final int a, final int b, final int c) {
		return Math.max(a, Math.max(b, c));
	}

	long lastDenseCells() {return lastDenseCells;}
	long lastProcessedCells() {return lastProcessedCells;}
	long lastFrontierBits() {return lastFrontierBits;}
	long lastKeptStates() {return lastKeptStates;}
	long lastPrunedStates() {return lastPrunedStates;}

	private int[] matchA=new int[0], matchB=new int[0];
	private int[] insertionA=new int[0], insertionB=new int[0];
	private int[] deletionA=new int[0], deletionB=new int[0];
	private long[] activeA=new long[0], activeB=new long[0], candidateWords=new long[0];
	private long lastDenseCells, lastProcessedCells, lastFrontierBits;
	private long lastKeptStates, lastPrunedStates;
	private final int match, substitution, ambiguity;
	private final int insertionOpen, insertionExtend, deletionOpen, deletionExtend;
	private static final int BAD=AffineScore32.BAD;
}
