package idaligner;

import java.util.Arrays;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/** Direct-word candidate for exact active anti-diagonal affine SIMD scoring.
 * Scoring enforces the conservative numeric/dimension contract of {@link AffineScore32}.
 */
@SuppressWarnings("restriction")
public final class AffineScore32DiagonalActiveSIMDWords {

	public AffineScore32DiagonalActiveSIMDWords(final int match_, final int substitution_,
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
		matchVector=IntVector.broadcast(SPECIES, match);
		substitutionVector=IntVector.broadcast(SPECIES, substitution);
		ambiguityVector=IntVector.broadcast(SPECIES, ambiguity);
		insertionOpenVector=IntVector.broadcast(SPECIES, insertionOpen);
		insertionExtendVector=IntVector.broadcast(SPECIES, insertionExtend);
		deletionOpenVector=IntVector.broadcast(SPECIES, deletionOpen);
		deletionExtendVector=IntVector.broadcast(SPECIES, deletionExtend);
	}

	public static AffineScore32DiagonalActiveSIMDWords msaLike() {
		return new AffineScore32DiagonalActiveSIMDWords(100, -127, 0, -395, -39, -472, -33);
	}

	/**
	 * Returns the exact full score when {@code incumbent} is an achievable path
	 * score for these exact query/reference bytes. A same-input banded score is
	 * the intended source; the numerical ceiling alone cannot prove achievability.
	 */
	public int score(final byte[] query, final byte[] ref, final int incumbent) {
		if(query==null) {throw new NullPointerException("A query is required.");}
		if(ref==null) {throw new NullPointerException("A reference is required.");}
		if(ref.length>0) {return score(query, ref, 0, ref.length-1, incumbent);}
		AffineScore32.validateScoreRange(query.length, 0, match, substitution, ambiguity,
				insertionOpen, insertionExtend, deletionOpen, deletionExtend);
		resetStats(query.length, ref.length);
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
		final int qLen=query.length;
		ensureCapacity(qLen, rLen);
		encode(query, ref, refStart, rLen);
		int[] mKm2=matchA, mKm1=matchB, mK=matchC;
		int[] iKm2=insertionA, iKm1=insertionB, iK=insertionC;
		int[] dKm2=deletionA, dKm1=deletionB, dK=deletionC;
		long[] wordsKm2=activeWordsA, wordsKm1=activeWordsB, wordsK=activeWordsC;
		final int wordCount=(rLen>>>6)+2;
		Arrays.fill(wordsKm2, 0, wordCount, 0L);
		Arrays.fill(wordsKm1, 0, wordCount, 0L);
		Arrays.fill(wordsK, 0, wordCount, 0L);
		int leftInsertion=BAD;
		int best=incumbent;

		for(int k=0; k<=qLen+rLen; k++) {
			Arrays.fill(wordsK, 0, wordCount, 0L);
			if(k<=rLen) {
				mK[k]=0; iK[k]=BAD; dK[k]=0;
				set(wordsK, k);
			}
			if(k>=1 && k<=qLen) {
				final int rawLeft=(k==1 ? insertionOpen : add(leftInsertion, insertionExtend));
				leftInsertion=prune(rawLeft, qLen-k, best);
				mK[0]=dK[0]=BAD; iK[0]=leftInsertion;
				if(leftInsertion!=BAD) {set(wordsK, 0);}
			}

			buildCandidates(wordsKm2, wordsKm1, candidateWords, wordCount, rLen);
			final int minCol=Math.max(1, k-qLen);
			final int maxCol=Math.min(rLen, k-1);
			for(int runStart=nextSetBit(candidateWords, minCol, maxCol+1);
					runStart>=0;) {
				final int runEnd=nextClearBit(candidateWords, runStart, maxCol+1);
				int col=runStart;
				for(; col+WIDTH<=runEnd; col+=WIDTH) {
					scoreVectorChunk(qLen, k, col, best, mKm2, iKm2, dKm2,
							mKm1, iKm1, dKm1, mK, iK, dK,
							wordsKm2, wordsKm1, wordsK);
					lastVectorChunks++;
					lastVectorCells+=WIDTH;
					lastProcessedCells+=WIDTH;
					final int bottomCol=k-qLen;
					if(bottomCol>=col && bottomCol<col+WIDTH) {
						best=Math.max(best, max3(mK[bottomCol], iK[bottomCol], dK[bottomCol]));
					}
				}
				for(; col<runEnd; col++) {
					final int row=k-col;
					final int rawM=get(wordsKm2, col-1) ?
							add(max3(mKm2[col-1], iKm2[col-1], dKm2[col-1]),
									baseScore(query[row-1], ref[refStart+col-1])) : BAD;
					final int rawI=get(wordsKm1, col) ?
							max3(add(mKm1[col], insertionOpen),
									add(iKm1[col], insertionExtend),
									add(dKm1[col], insertionOpen)) : BAD;
					final int rawD=get(wordsKm1, col-1) ?
							max3(add(mKm1[col-1], deletionOpen),
									add(dKm1[col-1], deletionExtend),
									add(iKm1[col-1], deletionOpen)) : BAD;
					final int remaining=qLen-row;
					mK[col]=prune(rawM, remaining, best);
					iK[col]=prune(rawI, remaining, best);
					dK[col]=prune(rawD, remaining, best);
					lastProcessedCells++;
					if(mK[col]!=BAD || iK[col]!=BAD || dK[col]!=BAD) {
						set(wordsK, col);
					}
					if(row==qLen) {best=Math.max(best, max3(mK[col], iK[col], dK[col]));}
				}
				runStart=nextSetBit(candidateWords, runEnd, maxCol+1);
			}

			int[] swap=mKm2; mKm2=mKm1; mKm1=mK; mK=swap;
			swap=iKm2; iKm2=iKm1; iKm1=iK; iK=swap;
			swap=dKm2; dKm2=dKm1; dKm1=dK; dK=swap;
			final long[] wordSwap=wordsKm2;
			wordsKm2=wordsKm1; wordsKm1=wordsK; wordsK=wordSwap;
		}
		return best;
	}

	private void scoreVectorChunk(final int qLen, final int k, final int col,
			final int best, final int[] mKm2, final int[] iKm2, final int[] dKm2,
			final int[] mKm1, final int[] iKm1, final int[] dKm1,
			final int[] mK, final int[] iK, final int[] dK,
			final long[] wordsKm2, final long[] wordsKm1, final long[] wordsK) {
		final IntVector q=IntVector.fromArray(SPECIES, reversedQuery, qLen+col-k);
		final IntVector r=IntVector.fromArray(SPECIES, encodedRef, col-1);
		final VectorMask<Integer> ambiguous=q.compare(VectorOperators.GT, EIGHT)
				.or(r.compare(VectorOperators.GT, EIGHT));
		final VectorMask<Integer> matching=q.compare(VectorOperators.EQ, r)
				.and(ambiguous.not());
		final IntVector base=substitutionVector.blend(matchVector, matching)
				.blend(ambiguityVector, ambiguous);

		final VectorMask<Integer> matchSource=mask(wordsKm2, col-1);
		final IntVector diagonal=IntVector.fromArray(SPECIES, mKm2, col-1)
				.max(IntVector.fromArray(SPECIES, iKm2, col-1))
				.max(IntVector.fromArray(SPECIES, dKm2, col-1));
		final IntVector rawM=BAD_VECTOR.blend(addReachable(diagonal, base), matchSource);

		final VectorMask<Integer> insertionSource=mask(wordsKm1, col);
		final IntVector rawI=BAD_VECTOR.blend(
				addReachable(IntVector.fromArray(SPECIES, mKm1, col), insertionOpenVector)
				.max(addReachable(IntVector.fromArray(SPECIES, iKm1, col), insertionExtendVector))
				.max(addReachable(IntVector.fromArray(SPECIES, dKm1, col), insertionOpenVector)),
				insertionSource);

		final VectorMask<Integer> deletionSource=mask(wordsKm1, col-1);
		final IntVector rawD=BAD_VECTOR.blend(
				addReachable(IntVector.fromArray(SPECIES, mKm1, col-1), deletionOpenVector)
				.max(addReachable(IntVector.fromArray(SPECIES, iKm1, col-1), deletionOpenVector))
				.max(addReachable(IntVector.fromArray(SPECIES, dKm1, col-1), deletionExtendVector)),
				deletionSource);

		for(int lane=0; lane<WIDTH; lane++) {
			final int remaining=qLen+col+lane-k;
			final long threshold=(long)best-(long)remaining*match;
			thresholdScratch[lane]=threshold<Integer.MIN_VALUE ? Integer.MIN_VALUE :
					(threshold>Integer.MAX_VALUE ? Integer.MAX_VALUE : (int)threshold);
		}
		final IntVector threshold=IntVector.fromArray(SPECIES, thresholdScratch, 0);
		final VectorMask<Integer> keepM=rawM.compare(VectorOperators.NE, BAD_VECTOR)
				.and(rawM.compare(VectorOperators.GT, threshold));
		final VectorMask<Integer> keepI=rawI.compare(VectorOperators.NE, BAD_VECTOR)
				.and(rawI.compare(VectorOperators.GT, threshold));
		final VectorMask<Integer> keepD=rawD.compare(VectorOperators.NE, BAD_VECTOR)
				.and(rawD.compare(VectorOperators.GT, threshold));
		BAD_VECTOR.blend(rawM, keepM).intoArray(mK, col);
		BAD_VECTOR.blend(rawI, keepI).intoArray(iK, col);
		BAD_VECTOR.blend(rawD, keepD).intoArray(dK, col);
		final long live=keepM.or(keepI).or(keepD).toLong()&LANE_MASK;
		orBits(wordsK, col, live);
	}

	private static VectorMask<Integer> mask(final long[] words, final int start) {
		final int word=start>>>6;
		final int shift=start&63;
		long lanes=words[word]>>>shift;
		if(shift>64-WIDTH) {lanes|=words[word+1]<<(64-shift);}
		return VectorMask.fromLong(SPECIES, lanes&LANE_MASK);
	}

	private static void set(final long[] words, final int bit) {
		words[bit>>>6]|=1L<<(bit&63);
	}

	private static boolean get(final long[] words, final int bit) {
		return (words[bit>>>6]&(1L<<(bit&63)))!=0;
	}

	private static void orBits(final long[] words, final int start, final long bits) {
		final int word=start>>>6, shift=start&63;
		words[word]|=bits<<shift;
		if(shift>64-WIDTH) {words[word+1]|=bits>>>(64-shift);}
	}

	private static IntVector addReachable(final IntVector prior, final IntVector delta) {
		return BAD_VECTOR.blend(prior.add(delta),
				prior.compare(VectorOperators.NE, BAD_VECTOR));
	}

	private static void buildCandidates(final long[] km2, final long[] km1,
			final long[] out, final int wordCount, final int maxColumn) {
		long carry1=0, carry2=0;
		for(int word=0; word<wordCount; word++) {
			final long a=km1[word], b=km2[word];
			out[word]=a|(a<<1)|carry1|(b<<1)|carry2;
			carry1=a>>>63; carry2=b>>>63;
		}
		out[0]&=~1L;
		final int lastWord=maxColumn>>>6, used=(maxColumn&63)+1;
		if(used<64) {out[lastWord]&=(1L<<used)-1L;}
		Arrays.fill(out, lastWord+1, wordCount, 0L);
	}

	private static int nextSetBit(final long[] words, final int from,
			final int limit) {
		if(from>=limit) {return -1;}
		int word=from>>>6;
		long bits=words[word]&(-1L<<(from&63));
		while(true) {
			if(bits!=0) {
				final int bit=(word<<6)+Long.numberOfTrailingZeros(bits);
				return bit<limit ? bit : -1;
			}
			word++;
			if((word<<6)>=limit) {return -1;}
			bits=words[word];
		}
	}

	private static int nextClearBit(final long[] words, final int from,
			final int limit) {
		int word=from>>>6;
		long clear=~words[word]&(-1L<<(from&63));
		while(true) {
			if(clear!=0) {
				return Math.min(limit, (word<<6)+Long.numberOfTrailingZeros(clear));
			}
			word++;
			if((word<<6)>=limit) {return limit;}
			clear=~words[word];
		}
	}

	private int prune(final int score, final int remaining, final int best) {
		if(score==BAD) {return BAD;}
		return (long)score+(long)remaining*match>best ? score : BAD;
	}

	private int baseScore(final byte queryBase, final byte refBase) {
		final int q=CrossCutScore32.baseCode(queryBase);
		final int r=CrossCutScore32.baseCode(refBase);
		return q==r && q<=8 ? match : (q>8 || r>8 ? ambiguity : substitution);
	}

	private void ensureCapacity(final int qLen, final int rLen) {
		if(reversedQuery.length<qLen) {reversedQuery=new int[qLen];}
		if(encodedRef.length<rLen) {encodedRef=new int[rLen];}
		final int wordLength=(rLen>>>6)+2;
		if(activeWordsA.length<wordLength) {
			activeWordsA=new long[wordLength]; activeWordsB=new long[wordLength];
			activeWordsC=new long[wordLength];
			candidateWords=new long[wordLength];
		}
		final int length=rLen+1;
		if(matchA.length>=length) {return;}
		matchA=new int[length]; matchB=new int[length]; matchC=new int[length];
		insertionA=new int[length]; insertionB=new int[length]; insertionC=new int[length];
		deletionA=new int[length]; deletionB=new int[length]; deletionC=new int[length];
	}

	private void encode(final byte[] query, final byte[] ref, final int refStart,
			final int rLen) {
		for(int i=0; i<query.length; i++) {
			reversedQuery[query.length-1-i]=CrossCutScore32.baseCode(query[i]);
		}
		for(int i=0; i<rLen; i++) {
			encodedRef[i]=CrossCutScore32.baseCode(ref[refStart+i]);
		}
	}

	private void resetStats(final int qLen, final int rLen) {
		lastDenseCells=(long)qLen*rLen;
		lastProcessedCells=lastVectorChunks=lastVectorCells=0;
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

	static int vectorWidth() {return WIDTH;}
	long lastDenseCells() {return lastDenseCells;}
	long lastProcessedCells() {return lastProcessedCells;}
	long lastVectorChunks() {return lastVectorChunks;}
	long lastVectorCells() {return lastVectorCells;}

	private int[] reversedQuery=new int[0], encodedRef=new int[0];
	private int[] matchA=new int[0], matchB=new int[0], matchC=new int[0];
	private int[] insertionA=new int[0], insertionB=new int[0], insertionC=new int[0];
	private int[] deletionA=new int[0], deletionB=new int[0], deletionC=new int[0];
	private final int[] thresholdScratch=new int[WIDTH];
	private long[] activeWordsA=new long[0], activeWordsB=new long[0], activeWordsC=new long[0];
	private long[] candidateWords=new long[0];
	private long lastDenseCells, lastProcessedCells, lastVectorChunks, lastVectorCells;
	private final int match, substitution, ambiguity;
	private final int insertionOpen, insertionExtend, deletionOpen, deletionExtend;
	private final IntVector matchVector, substitutionVector, ambiguityVector;
	private final IntVector insertionOpenVector, insertionExtendVector;
	private final IntVector deletionOpenVector, deletionExtendVector;
	private static final VectorSpecies<Integer> SPECIES=IntVector.SPECIES_256;
	private static final int WIDTH=SPECIES.length();
	private static final long LANE_MASK=(1L<<WIDTH)-1L;
	private static final int BAD=AffineScore32.BAD;
	private static final IntVector BAD_VECTOR=IntVector.broadcast(SPECIES, BAD);
	private static final IntVector EIGHT=IntVector.broadcast(SPECIES, 8);
}
