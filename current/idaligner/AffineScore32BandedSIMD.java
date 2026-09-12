package idaligner;

import java.util.Arrays;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * Fixed-corridor form of {@link AffineScore32SIMD}. Only cells satisfying
 * {@code |(col-row)-centerDiagonal| <= halfWidth} are evaluated. The full
 * query is aligned while reference prefix and suffix remain free.
 *
 * <p>One instance owns nine reusable score arrays and is not thread-safe.</p>
 * <p>Scoring enforces the conservative numeric/dimension contract of
 * {@link AffineScore32}, using the full selected window even for a narrow band.</p>
 */
@SuppressWarnings("restriction")
public final class AffineScore32BandedSIMD {

	public AffineScore32BandedSIMD(final int match_, final int substitution_,
			final int ambiguity_, final int insertionOpen_,
			final int insertionExtend_, final int deletionOpen_,
			final int deletionExtend_) {
		if(match_<=0 || substitution_>0 || insertionOpen_>0 ||
				insertionExtend_>0 || deletionOpen_>0 || deletionExtend_>0) {
			throw new IllegalArgumentException(
					"Scores must use positive matches and nonpositive penalties.");
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

	public static AffineScore32BandedSIMD msaLike() {
		return new AffineScore32BandedSIMD(100, -127, 0, -395, -39, -472, -33);
	}

	/** Writes {@code {score, inclusive reference stop, terminal state}}. */
	public int scoreWithEndpoint(final byte[] query, final byte[] ref,
			final int centerDiagonal, final int halfWidth, final int[] result) {
		validateResult(result);
		final int score=score(query, ref, centerDiagonal, halfWidth);
		writeEndpoint(result, score);
		return score;
	}

	/** Writes {@code {score, absolute inclusive reference stop, terminal state}}. */
	public int scoreWithEndpoint(final byte[] query, final byte[] ref,
			final int refStart, final int refEnd, final int centerDiagonal,
			final int halfWidth, final int[] result) {
		validateResult(result);
		final int score=score(query, ref, refStart, refEnd, centerDiagonal, halfWidth);
		writeEndpoint(result, score);
		return score;
	}

	public int score(final byte[] query, final byte[] ref,
			final int centerDiagonal, final int halfWidth) {
		if(query==null) {throw new NullPointerException("A query is required.");}
		if(ref==null) {throw new NullPointerException("A reference is required.");}
		validateWidth(halfWidth);
		resetStats(query.length, ref.length);
		if(ref.length==0) {
			AffineScore32.validateScoreRange(query.length, 0, match, substitution, ambiguity,
					insertionOpen, insertionExtend, deletionOpen, deletionExtend);
			final int score=emptyReferenceScore(query.length, centerDiagonal, halfWidth);
			if(score!=BAD) {lastBestState=query.length==0 ? MATCH_STATE : INSERTION_STATE;}
			return score;
		}
		return score(query, ref, 0, ref.length-1, centerDiagonal, halfWidth);
	}

	public int score(final byte[] query, final byte[] ref, final int refStart,
			final int refEnd, final int centerDiagonal, final int halfWidth) {
		if(query==null) {throw new NullPointerException("A query is required.");}
		validateRange(ref, refStart, refEnd);
		validateWidth(halfWidth);
		final int qLen=query.length;
		final int rLen=refEnd-refStart+1;
		AffineScore32.validateScoreRange(qLen, rLen, match, substitution, ambiguity,
				insertionOpen, insertionExtend, deletionOpen, deletionExtend);
		resetStats(qLen, rLen);
		if(qLen==0) {
			for(int col=1; col<=rLen; col++) {
				if(inCorridor(0, col, centerDiagonal, halfWidth)) {
					lastBestColumn=col;
					lastBestRefStop=refStart+col-1;
					lastBestState=MATCH_STATE;
					lastBestOnCorridorEdge=isOnCorridorEdge(0, col,
							centerDiagonal, halfWidth);
					return 0;
				}
			}
			return BAD;
		}
		ensureCapacity(qLen, rLen);
		encode(query, ref, refStart, rLen);
		clearAll(rLen+1);

		int[] mKm2=matchA, mKm1=matchB, mK=matchC;
		int[] iKm2=insertionA, iKm1=insertionB, iK=insertionC;
		int[] dKm2=deletionA, dKm1=deletionB, dK=deletionC;
		int minKm2=1, maxKm2=0, topKm2=-1; boolean leftKm2=false;
		int minKm1=1, maxKm1=0, topKm1=-1; boolean leftKm1=false;
		int minK=1, maxK=0, topK=-1; boolean leftK=false;
		int leftInsertion=BAD;
		int best=BAD, bestCol=-1, bestState=-1;

		for(int k=0; k<=qLen+rLen; k++) {
			clearTouched(mK, iK, dK, minK, maxK, topK, leftK);
			final int newTop=k<=rLen && inCorridor(0, k, centerDiagonal, halfWidth) ? k : -1;
			if(newTop>=0) {
				mK[newTop]=0;
				iK[newTop]=BAD;
				dK[newTop]=0;
			}
			final boolean newLeft=k>=1 && k<=qLen &&
					inCorridor(k, 0, centerDiagonal, halfWidth);
			if(newLeft) {
				if(k==1) {
					leftInsertion=inCorridor(0, 0, centerDiagonal, halfWidth) ?
							insertionOpen : BAD;
				}else {
					leftInsertion=add(leftInsertion, insertionExtend);
				}
				mK[0]=BAD;
				iK[0]=leftInsertion;
				dK[0]=BAD;
			}else {
				leftInsertion=BAD;
			}

			final int fullMin=Math.max(1, k-qLen);
			final int fullMax=Math.min(rLen, k-1);
			final int corridorMin=clampToInt(ceilHalf(
					(long)k+centerDiagonal-halfWidth));
			final int corridorMax=clampToInt(Math.floorDiv(
					(long)k+centerDiagonal+halfWidth, 2));
			final int activeMin=Math.max(fullMin, corridorMin);
			final int activeMax=Math.min(fullMax, corridorMax);
			if(activeMin<=activeMax) {
				lastVisitedCells+=activeMax-activeMin+1L;
				int col=activeMin;
				for(; col+WIDTH-1<=activeMax; col+=WIDTH) {
					scoreVectorCells(qLen, k, col, mKm2, iKm2, dKm2,
							mKm1, iKm1, dKm1, mK, iK, dK);
					lastVectorChunks++;
					lastVectorCells+=WIDTH;
				}
				for(; col<=activeMax; col++) {
					scoreScalarCell(qLen, k, col, mKm2, iKm2, dKm2,
							mKm1, iKm1, dKm1, mK, iK, dK);
				}
			}

			final int bottomCol=k-qLen;
			if(bottomCol>=activeMin && bottomCol<=activeMax) {
				final int m=mK[bottomCol], i=iK[bottomCol], d=dK[bottomCol];
				final int score=max3(m, i, d);
				if(score>best) {
					best=score;
					bestCol=bottomCol;
					bestState=m>=i && m>=d ? MATCH_STATE :
							(i>=d ? INSERTION_STATE : DELETION_STATE);
				}
			}

			int[] swap=mKm2; mKm2=mKm1; mKm1=mK; mK=swap;
			swap=iKm2; iKm2=iKm1; iKm1=iK; iK=swap;
			swap=dKm2; dKm2=dKm1; dKm1=dK; dK=swap;
			int swapInt=minKm2; minKm2=minKm1; minKm1=activeMin; minK=swapInt;
			swapInt=maxKm2; maxKm2=maxKm1; maxKm1=activeMax; maxK=swapInt;
			swapInt=topKm2; topKm2=topKm1; topKm1=newTop; topK=swapInt;
			final boolean swapBoolean=leftKm2;
			leftKm2=leftKm1; leftKm1=newLeft; leftK=swapBoolean;
		}

		if(bestCol>=1) {
			lastBestColumn=bestCol;
			lastBestRefStop=refStart+bestCol-1;
			lastBestState=bestState;
			lastBestOnCorridorEdge=isOnCorridorEdge(qLen, bestCol,
					centerDiagonal, halfWidth);
		}
		return best;
	}

	private void writeEndpoint(final int[] result, final int score) {
		result[0]=score;
		result[1]=lastBestRefStop;
		result[2]=lastBestState;
	}

	private void scoreVectorCells(final int qLen, final int k, final int col,
			final int[] mKm2, final int[] iKm2, final int[] dKm2,
			final int[] mKm1, final int[] iKm1, final int[] dKm1,
			final int[] mK, final int[] iK, final int[] dK) {
		final int queryIndex=qLen+col-k;
		final IntVector q=IntVector.fromArray(SPECIES, reversedQuery, queryIndex);
		final IntVector r=IntVector.fromArray(SPECIES, encodedRef, col-1);
		final VectorMask<Integer> ambiguous=q.compare(VectorOperators.GT, EIGHT)
				.or(r.compare(VectorOperators.GT, EIGHT));
		final VectorMask<Integer> matching=q.compare(VectorOperators.EQ, r)
				.and(ambiguous.not());
		final IntVector base=substitutionVector.blend(matchVector, matching)
				.blend(ambiguityVector, ambiguous);

		final IntVector diagonalM=IntVector.fromArray(SPECIES, mKm2, col-1);
		final IntVector diagonalI=IntVector.fromArray(SPECIES, iKm2, col-1);
		final IntVector diagonalD=IntVector.fromArray(SPECIES, dKm2, col-1);
		final IntVector diagonalBest=diagonalM.max(diagonalI).max(diagonalD);
		BAD_VECTOR.blend(diagonalBest.add(base),
				diagonalBest.compare(VectorOperators.NE, BAD_VECTOR)).intoArray(mK, col);

		final IntVector upM=IntVector.fromArray(SPECIES, mKm1, col);
		final IntVector upI=IntVector.fromArray(SPECIES, iKm1, col);
		final IntVector upD=IntVector.fromArray(SPECIES, dKm1, col);
		addReachable(upM, insertionOpenVector).max(
				addReachable(upI, insertionExtendVector)).max(
				addReachable(upD, insertionOpenVector)).intoArray(iK, col);

		final IntVector leftM=IntVector.fromArray(SPECIES, mKm1, col-1);
		final IntVector leftI=IntVector.fromArray(SPECIES, iKm1, col-1);
		final IntVector leftD=IntVector.fromArray(SPECIES, dKm1, col-1);
		addReachable(leftM, deletionOpenVector).max(
				addReachable(leftI, deletionOpenVector)).max(
				addReachable(leftD, deletionExtendVector)).intoArray(dK, col);
	}

	private static IntVector addReachable(final IntVector prior,
			final IntVector delta) {
		return BAD_VECTOR.blend(prior.add(delta),
				prior.compare(VectorOperators.NE, BAD_VECTOR));
	}

	private void scoreScalarCell(final int qLen, final int k, final int col,
			final int[] mKm2, final int[] iKm2, final int[] dKm2,
			final int[] mKm1, final int[] iKm1, final int[] dKm1,
			final int[] mK, final int[] iK, final int[] dK) {
		final int queryIndex=qLen+col-k;
		final int base=baseScore(reversedQuery[queryIndex], encodedRef[col-1]);
		mK[col]=add(max3(mKm2[col-1], iKm2[col-1], dKm2[col-1]), base);
		iK[col]=max3(add(mKm1[col], insertionOpen),
				add(iKm1[col], insertionExtend), add(dKm1[col], insertionOpen));
		dK[col]=max3(add(mKm1[col-1], deletionOpen),
				add(dKm1[col-1], deletionExtend), add(iKm1[col-1], deletionOpen));
	}

	private void ensureCapacity(final int qLen, final int rLen) {
		if(reversedQuery.length<qLen) {reversedQuery=new int[qLen];}
		if(encodedRef.length<rLen) {encodedRef=new int[rLen];}
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

	private static void clearTouched(final int[] m, final int[] i, final int[] d,
			final int min, final int max, final int top, final boolean left) {
		if(min<=max) {
			Arrays.fill(m, min, max+1, BAD);
			Arrays.fill(i, min, max+1, BAD);
			Arrays.fill(d, min, max+1, BAD);
		}
		if(top>=0) {m[top]=BAD; i[top]=BAD; d[top]=BAD;}
		if(left) {m[0]=BAD; i[0]=BAD; d[0]=BAD;}
	}

	private void clearAll(final int length) {
		Arrays.fill(matchA, 0, length, BAD); Arrays.fill(matchB, 0, length, BAD);
		Arrays.fill(matchC, 0, length, BAD); Arrays.fill(insertionA, 0, length, BAD);
		Arrays.fill(insertionB, 0, length, BAD); Arrays.fill(insertionC, 0, length, BAD);
		Arrays.fill(deletionA, 0, length, BAD); Arrays.fill(deletionB, 0, length, BAD);
		Arrays.fill(deletionC, 0, length, BAD);
	}

	private void resetStats(final int qLen, final int rLen) {
		lastDenseCells=(long)qLen*rLen;
		lastVisitedCells=lastVectorChunks=lastVectorCells=0;
		lastBestColumn=lastBestRefStop=lastBestState=-1;
		lastBestOnCorridorEdge=false;
	}

	private int emptyReferenceScore(final int qLen, final int center,
			final int halfWidth) {
		for(int row=0; row<=qLen; row++) {
			if(!inCorridor(row, 0, center, halfWidth)) {return BAD;}
		}
		return qLen==0 ? 0 : insertionOpen+(qLen-1)*insertionExtend;
	}

	private int baseScore(final int q, final int r) {
		return q==r && q<=8 ? match : (q>8 || r>8 ? ambiguity : substitution);
	}

	private static int add(final int prior, final int delta) {
		return prior==BAD ? BAD : prior+delta;
	}

	private static int max3(final int a, final int b, final int c) {
		return Math.max(a, Math.max(b, c));
	}

	private static boolean inCorridor(final int row, final int col,
			final int center, final int halfWidth) {
		final long difference=(long)col-row-center;
		return difference>=-halfWidth && difference<=halfWidth;
	}

	private static boolean isOnCorridorEdge(final int row, final int col,
			final int center, final int halfWidth) {
		final long difference=(long)col-row-center;
		return difference==-halfWidth || difference==halfWidth;
	}

	private static long ceilHalf(final long value) {return -Math.floorDiv(-value, 2);}
	private static int clampToInt(final long value) {
		return value<Integer.MIN_VALUE ? Integer.MIN_VALUE :
				(value>Integer.MAX_VALUE ? Integer.MAX_VALUE : (int)value);
	}

	private static void validateWidth(final int halfWidth) {
		if(halfWidth<0) {throw new IllegalArgumentException("Negative half-width: "+halfWidth);}
	}

	private static void validateResult(final int[] result) {
		if(result==null || result.length<3) {
			throw new IllegalArgumentException("A three-int endpoint result buffer is required.");
		}
	}

	private static void validateRange(final byte[] ref, final int refStart,
			final int refEnd) {
		if(ref==null) {throw new NullPointerException("A reference is required.");}
		if(refStart<0 || refEnd<refStart || refEnd>=ref.length) {
			throw new IndexOutOfBoundsException("Invalid inclusive reference range "+
					refStart+".."+refEnd+" for length "+ref.length);
		}
	}

	static int vectorWidth() {return WIDTH;}
	static int scoreArrayCount() {return 9;}
	long lastDenseCells() {return lastDenseCells;}
	long lastVisitedCells() {return lastVisitedCells;}
	long lastVectorChunks() {return lastVectorChunks;}
	long lastVectorCells() {return lastVectorCells;}
	int lastBestColumn() {return lastBestColumn;}
	int lastBestRefStop() {return lastBestRefStop;}
	int lastBestState() {return lastBestState;}
	/** Returns whether the most recent winning terminal cell lay on either corridor edge. */
	public boolean lastBestOnCorridorEdge() {return lastBestOnCorridorEdge;}

	private int[] reversedQuery=new int[0], encodedRef=new int[0];
	private int[] matchA=new int[0], matchB=new int[0], matchC=new int[0];
	private int[] insertionA=new int[0], insertionB=new int[0], insertionC=new int[0];
	private int[] deletionA=new int[0], deletionB=new int[0], deletionC=new int[0];
	private long lastDenseCells, lastVisitedCells, lastVectorChunks, lastVectorCells;
	private int lastBestColumn=-1, lastBestRefStop=-1, lastBestState=-1;
	private boolean lastBestOnCorridorEdge;

	private final int match, substitution, ambiguity;
	private final int insertionOpen, insertionExtend, deletionOpen, deletionExtend;
	private final IntVector matchVector, substitutionVector, ambiguityVector;
	private final IntVector insertionOpenVector, insertionExtendVector;
	private final IntVector deletionOpenVector, deletionExtendVector;

	public static final int MATCH_STATE=0;
	public static final int DELETION_STATE=1;
	public static final int INSERTION_STATE=2;
	private static final VectorSpecies<Integer> SPECIES=IntVector.SPECIES_256;
	private static final int WIDTH=SPECIES.length();
	private static final int BAD=AffineScore32.BAD;
	private static final IntVector BAD_VECTOR=IntVector.broadcast(SPECIES, BAD);
	private static final IntVector EIGHT=IntVector.broadcast(SPECIES, 8);
}
