package simd;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorSpecies;
import static jdk.incubator.vector.VectorOperators.*;
import static align2.MultiStateAligner11ts.*;

/**
 * Vectorized horizontal fill for the live interior of a limited MSA row.
 * <p>The caller advances one query row at a time and has already selected a
 * live inclusive {@code [from, to]} span from the previous row, the half-band
 * and the minimum-score frontier. This helper fills complete eight-column
 * blocks of the match/substitution and insertion states inside that span,
 * records compact ancestry flags, and returns the first column left for the
 * scalar remainder. Deletion states use a left dependency in the current row,
 * so the caller computes them after this method returns.</p>
 * <p>All score arrays use the packed {@code score|time} representation from
 * {@code MultiStateAligner11ts}. The method reads prior row states at
 * {@code col-1} and {@code col}, reference codes at {@code col-1} and
 * {@code col}, and writes current row match/insertion states plus trace flags at
 * {@code col}. The three output arrays must be distinct. Input arrays must not
 * overlap with outputs in a way that lets a write to {@code [from, to]} change
 * a later vector read.</p>
 *
 * @author Collei
 */
public final class SIMDMSALimited{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs this stateless helper. */
	public SIMDMSALimited(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Fills complete vector blocks in one live row span.
	 *
	 * @param pm Prior row match-state scores
	 * @param pd Prior row deletion-state scores
	 * @param pi Prior row insertion-state scores
	 * @param cm Current row match-state output
	 * @param ci Current row insertion-state output
	 * @param refs Reference bases indexed by MSA column
	 * @param flags Compact predecessor flags, indexed by MSA column
	 * @param q Current query base
	 * @param q0 Previous query base, or a nonmatching sentinel for the first row
	 * @param row Current query row, 1-based
	 * @param rows Number of query rows
	 * @param columns Number of reference columns
	 * @param floor Packed dead-cell score written into blocked cells
	 * @param insertionBarrier Row band that blocks insertions near query edges
	 * @param gapSymbol Reference sentinel that marks a gap column
	 * @param from First live score column to fill, inclusive
	 * @param to Last live score column to fill, inclusive
	 * @param target Packed minimum score target used to prune unreachable cells
	 * @return First column not filled by the vector loop
	 */
	public static int fillRow(final int[] pm, final int[] pd, final int[] pi,
			final int[] cm, final int[] ci, final int[] refs, final int[] flags,
			final int q, final int q0, final int row, final int rows, final int columns,
			final int floor, final int insertionBarrier, final int gapSymbol,
			final int from, final int to, final int target){
		assert(from>=1 && to<=columns && row>=1 && row<=rows) :
			"Limited vector span must lie within the retained row buffers";
		final int limit=to-7;
		int col=from;
		if(col>limit){return col;}
		final IntVector query=C(q), previousQuery=C(q0), gapCode=C(gapSymbol),
				dead=C(floor), lastColumn=C(columns-1);
		final IntVector targetScore=C(target), remainingQuery=C(rows-row),
				referenceColumns=C(columns);
		IntVector positions=LANES.add(col);
		for(; col<=limit; col+=8){
			// Bound remaining reference matches by query length before shifted-score multiplication.
			final IntVector cutoff=targetScore.sub(remainingQuery
					.min(referenceColumns.sub(positions)).mul(MATCH2));
			final IntVector ref=IntVector.fromArray(S, refs, col),
					prevRef=IntVector.fromArray(S, refs, col-1);
			final VectorMask<Integer> gap=ref.eq(gapCode),
					match=ref.eq(query).and(ref.compare(NE, N)),
					prevMatch=prevRef.eq(previousQuery).and(prevRef.compare(NE, N));
			final IntVector diag=IntVector.fromArray(S, pm, col-1),
					a=diag.and(SCORE_MASK),
					b=IntVector.fromArray(S, pd, col-1).and(SCORE_MASK),
					c=IntVector.fromArray(S, pi, col-1).and(SCORE_MASK),
					streak=diag.and(TIME_MASK);
			IntVector sub=SUB3.blend(SUB2, streak.compare(LT, COST3))
					.blend(SUB, streak.eq(ZERO));
			sub=sub.blend(SUB.blend(SUBR, streak.compare(LE, ONE)), prevMatch);
			if(q=='N'){sub=NOCALL;}else{sub=sub.blend(NOCALL, ref.eq(N));}
			final IntVector gain=sub.blend(MATCH.blend(MATCH2, prevMatch), match);
			final IntVector fromMS=a.add(gain), otherGain=SUB.blend(MATCH, match),
					fromD=b.add(otherGain), fromI=c.add(otherGain);
			final VectorMask<Integer> chooseMS=fromMS.compare(GE, fromD)
					.and(fromMS.compare(GE, fromI));
			IntVector time=ONE.blend(streak.add(ONE), match.eq(prevMatch).and(chooseMS));
			time=time.blend(WRAPPED_TIME, time.compare(GT, TIME_MAX)).blend(ZERO, gap);
			IntVector ms=fromMS.max(fromD).max(fromI).or(time).blend(dead, gap);
			ms.blend(dead, ms.compare(LT, cutoff)).intoArray(cm, col);
			IntVector prevMS=TWO.blend(ONE, b.compare(GE, c))
					.blend(ZERO, a.compare(GE, b).and(a.compare(GE, c)))
					.blend(ZERO, time.compare(GT, ONE));
			final IntVector up=IntVector.fromArray(S, pm, col).and(SCORE_MASK),
					upPacked=IntVector.fromArray(S, pi, col),
					upI=upPacked.and(SCORE_MASK), is=upPacked.and(TIME_MASK);
			final IntVector icost=INS4.blend(INS3, is.compare(LT, COST4))
					.blend(INS2, is.compare(LT, COST3)).blend(INS, is.eq(ZERO));
			final IntVector x=up.add(INS), y=upI.add(icost);
			IntVector itime=ONE.blend(is.add(ONE), y.compare(GT, x));
			itime=itime.blend(WRAPPED_TIME, itime.compare(GT, TIME_MAX));
			VectorMask<Integer> blocked=gap;
			if(row<insertionBarrier){blocked=blocked.or(positions.compare(GT, ONE));}
			if(row>rows-insertionBarrier){blocked=blocked.or(positions.compare(LT, lastColumn));}
			itime=itime.blend(ZERO, blocked);
			IntVector ins=x.max(y).or(itime).blend(dead, blocked);
			ins.blend(dead, ins.compare(LT, cutoff)).intoArray(ci, col);
			prevMS.or(ZERO.blend(EIGHT, itime.compare(GT, ONE)
					.or(upI.compare(GT, up)))).intoArray(flags, col);
			positions=positions.add(EIGHT);
		}
		return col;
	}

	/**
	 * Broadcasts a scalar value into every integer lane.
	 *
	 * @param value Value to broadcast
	 * @return Vector with {@code value} in every lane
	 */
	private static IntVector C(final int value){return IntVector.broadcast(S, value);}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** 256-bit integer vector species; all limited row loops assume its eight lanes. */
	private static final VectorSpecies<Integer> S=IntVector.SPECIES_256;
	/** Lane indices added to the first score column of the active vector block. */
	private static final IntVector LANES=
			IntVector.fromArray(S, new int[]{0, 1, 2, 3, 4, 5, 6, 7}, 0);
	/** Zero in every lane. */
	private static final IntVector ZERO=C(0);
	/** One in every lane. */
	private static final IntVector ONE=C(1);
	/** Two in every lane. */
	private static final IntVector TWO=C(2);
	/** Trace flag for insertion transitions that extend a prior insertion. */
	private static final IntVector EIGHT=C(8);
	/** No-call reference/query sentinel. */
	private static final IntVector N=C('N');
	/** Packed-score bit mask. */
	private static final IntVector SCORE_MASK=C(SCOREMASK);
	/** Packed time/streak bit mask. */
	private static final IntVector TIME_MASK=C(TIMEMASK);
	/** Maximum packed time before wrap. */
	private static final IntVector TIME_MAX=C(MAX_TIME);
	/** Wrapped packed time value. */
	private static final IntVector WRAPPED_TIME=C(MAX_TIME-MASK5);
	/** First-substitution penalty. */
	private static final IntVector SUB=C(POINTSoff_SUB);
	/** Second-substitution penalty. */
	private static final IntVector SUB2=C(POINTSoff_SUB2);
	/** Later-substitution penalty. */
	private static final IntVector SUB3=C(POINTSoff_SUB3);
	/** Repeat-substitution penalty after a previous same base. */
	private static final IntVector SUBR=C(POINTSoff_SUBR);
	/** No-call penalty. */
	private static final IntVector NOCALL=C(POINTSoff_NOCALL);
	/** Match reward. */
	private static final IntVector MATCH=C(POINTSoff_MATCH);
	/** Match reward after a previous match. */
	private static final IntVector MATCH2=C(POINTSoff_MATCH2);
	/** First-insertion penalty. */
	private static final IntVector INS=C(POINTSoff_INS);
	/** Second-insertion penalty. */
	private static final IntVector INS2=C(POINTSoff_INS2);
	/** Third-insertion penalty. */
	private static final IntVector INS3=C(POINTSoff_INS3);
	/** Fourth-or-later insertion penalty. */
	private static final IntVector INS4=C(POINTSoff_INS4);
	/** Repeat threshold for the third substitution/insertion cost. */
	private static final IntVector COST3=C(LIMIT_FOR_COST_3);
	/** Repeat threshold for the fourth insertion cost. */
	private static final IntVector COST4=C(LIMIT_FOR_COST_4);
}
