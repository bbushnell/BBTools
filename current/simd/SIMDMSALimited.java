package simd;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorSpecies;
import static jdk.incubator.vector.VectorOperators.*;
import static align2.MultiStateAligner11ts.*;

/** Lazily linked8-lane MS/INS row fill; DEL's left dependency is handled by the caller.
 * @author Collei */
public final class SIMDMSALimited {
	/** Fills complete vectors and returns the first remaining scalar column. */
	public static int fillRow(int[] pm, int[] pd, int[] pi, int[] cm, int[] ci, int[] refs, int[] flags, int q, int q0, int row, int rows, int columns, int floor, int insertionBarrier, int gapSymbol, int from, int to, int target){
		assert(from>=1 && to<=columns && row>=1 && row<=rows) : "Limited vector span must lie within the retained row buffers";
		final int limit=to-7;int col=from;
		if(col>limit){return col;}
		final IntVector query=C(q), previousQuery=C(q0), gapCode=C(gapSymbol), dead=C(floor), lastColumn=C(columns-1);
		final IntVector targetScore=C(target), remainingQuery=C(rows-row), referenceColumns=C(columns);
		IntVector positions=LANES.add(col);
		for(; col<=limit; col+=8){
			// Bound remaining reference matches by query length before shifted-score multiplication.
			final IntVector cutoff=targetScore.sub(remainingQuery.min(referenceColumns.sub(positions)).mul(MATCH2));
			final IntVector ref=IntVector.fromArray(S, refs, col), prevRef=IntVector.fromArray(S, refs, col-1);
			final VectorMask<Integer> gap=ref.eq(gapCode), match=ref.eq(query).and(ref.compare(NE, N)), prevMatch=prevRef.eq(previousQuery).and(prevRef.compare(NE, N));
			final IntVector diag=IntVector.fromArray(S, pm, col-1), a=diag.and(SCORE_MASK), b=IntVector.fromArray(S, pd, col-1).and(SCORE_MASK), c=IntVector.fromArray(S, pi, col-1).and(SCORE_MASK), streak=diag.and(TIME_MASK);
			IntVector sub=SUB3.blend(SUB2, streak.compare(LT, COST3)).blend(SUB, streak.eq(ZERO));
			sub=sub.blend(SUB.blend(SUBR, streak.compare(LE, ONE)), prevMatch);
			if(q=='N'){sub=NOCALL;}else{sub=sub.blend(NOCALL, ref.eq(N));}
			final IntVector gain=sub.blend(MATCH.blend(MATCH2, prevMatch), match);
			final IntVector fromMS=a.add(gain), otherGain=SUB.blend(MATCH, match), fromD=b.add(otherGain), fromI=c.add(otherGain);
			final VectorMask<Integer> chooseMS=fromMS.compare(GE, fromD).and(fromMS.compare(GE, fromI));
			IntVector time=ONE.blend(streak.add(ONE), match.eq(prevMatch).and(chooseMS));
			time=time.blend(WRAPPED_TIME, time.compare(GT, TIME_MAX)).blend(ZERO, gap);
			IntVector ms=fromMS.max(fromD).max(fromI).or(time).blend(dead, gap);
			ms.blend(dead, ms.compare(LT, cutoff)).intoArray(cm, col);
			IntVector prevMS=TWO.blend(ONE, b.compare(GE, c)).blend(ZERO, a.compare(GE, b).and(a.compare(GE, c))).blend(ZERO, time.compare(GT, ONE));
			final IntVector up=IntVector.fromArray(S, pm, col).and(SCORE_MASK), upPacked=IntVector.fromArray(S, pi, col), upI=upPacked.and(SCORE_MASK), is=upPacked.and(TIME_MASK);
			final IntVector icost=INS4.blend(INS3, is.compare(LT, COST4)).blend(INS2, is.compare(LT, COST3)).blend(INS, is.eq(ZERO));
			final IntVector x=up.add(INS), y=upI.add(icost);
			IntVector itime=ONE.blend(is.add(ONE), y.compare(GT, x));itime=itime.blend(WRAPPED_TIME, itime.compare(GT, TIME_MAX));
			VectorMask<Integer> blocked=gap;
			if(row<insertionBarrier){blocked=blocked.or(positions.compare(GT, ONE));}
			if(row>rows-insertionBarrier){blocked=blocked.or(positions.compare(LT, lastColumn));}
			itime=itime.blend(ZERO, blocked);IntVector ins=x.max(y).or(itime).blend(dead, blocked);
			ins.blend(dead, ins.compare(LT, cutoff)).intoArray(ci, col);
			prevMS.or(ZERO.blend(EIGHT, itime.compare(GT, ONE).or(upI.compare(GT, up)))).intoArray(flags, col);
			positions=positions.add(EIGHT);
		}
		return col;
	}
	private static IntVector C(int value){return IntVector.broadcast(S, value);}
	private static final VectorSpecies<Integer> S=IntVector.SPECIES_256;
	private static final IntVector LANES=IntVector.fromArray(S, new int[]{0, 1, 2, 3, 4, 5, 6, 7}, 0);
	/** Immutable species-specific vectors remain inside this lazily linked SIMD class. */
	private static final IntVector ZERO=C(0), ONE=C(1), TWO=C(2), EIGHT=C(8), N=C('N');
	private static final IntVector SCORE_MASK=C(SCOREMASK), TIME_MASK=C(TIMEMASK), TIME_MAX=C(MAX_TIME), WRAPPED_TIME=C(MAX_TIME-MASK5);
	private static final IntVector SUB=C(POINTSoff_SUB), SUB2=C(POINTSoff_SUB2), SUB3=C(POINTSoff_SUB3), SUBR=C(POINTSoff_SUBR), NOCALL=C(POINTSoff_NOCALL);
	private static final IntVector MATCH=C(POINTSoff_MATCH), MATCH2=C(POINTSoff_MATCH2), INS=C(POINTSoff_INS), INS2=C(POINTSoff_INS2), INS3=C(POINTSoff_INS3), INS4=C(POINTSoff_INS4);
	private static final IntVector COST3=C(LIMIT_FOR_COST_3), COST4=C(LIMIT_FOR_COST_4);
}
