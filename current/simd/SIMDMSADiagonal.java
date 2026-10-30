package simd;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorSpecies;
import static jdk.incubator.vector.VectorOperators.*;
import static align2.MultiStateAligner11ts.*;

/**
 * Vectorized antidiagonal fill for the three affine MSA states.
 * <p>Each eight-lane block advances adjacent columns on one antidiagonal.
 * The caller owns diagonal scheduling and the state buffers; this class computes
 * match/substitution, deletion and insertion cells, plus one compact trace byte
 * per column. Any columns left after the last full vector are intended to use the
 * same transition rules through the scalar remainder path.</p>
 *
 * @author Collei
 */
public final class SIMDMSADiagonal{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Constructs this stateless helper. */
	public SIMDMSADiagonal(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Fills the inclusive column range of one antidiagonal.
	 * <p>All score arrays use the packed {@code score|time} representation from
	 * {@code MultiStateAligner11ts}. The method reads old match/deletion/insertion
	 * states at {@code col-1}, prior match/deletion/insertion states at
	 * {@code col-1} or {@code col}, reference bases at {@code col-1} and
	 * {@code col}, and query bases at {@code rows-diagonal+col} and the following
	 * query slot. It writes current match, deletion and insertion states at
	 * {@code col}, and writes trace flags to {@code trace[offset+col]}.</p>
	 * <p>The three score output arrays must be distinct. Input and output arrays must not overlap in
	 * any way that lets an earlier write in {@code [lo, hi]} change a later read.
	 * Required live ranges are {@code refs[lo-1..hi]},
	 * {@code query[rows-diagonal+lo..rows-diagonal+hi+1]},
	 * old states {@code [lo-1..hi-1]}, prior states {@code [lo-1..hi]}, current
	 * states {@code [lo..hi]} and {@code trace[offset+lo..offset+hi]}.</p>
	 *
	 * @param om Source match-state antidiagonal for match/substitution transitions
	 * @param od Source deletion-state antidiagonal for match/substitution transitions
	 * @param oi Source insertion-state antidiagonal for match/substitution transitions
	 * @param pm Prior match-state antidiagonal for deletion/insertion transitions
	 * @param pd Prior deletion-state antidiagonal for deletion extensions
	 * @param pi Prior insertion-state antidiagonal for insertion extensions
	 * @param cm Current match-state output
	 * @param cd Current deletion-state output
	 * @param ci Current insertion-state output
	 * @param refs Reference bases indexed by MSA column
	 * @param query Query bases indexed by MSA row in reverse row order
	 * @param trace Compact traceback flags, addressed as {@code offset+col}
	 * @param offset Offset from score columns to the first trace byte
	 * @param diagonal Antidiagonal number; each row is {@code diagonal-col}
	 * @param rows Number of query rows
	 * @param columns Number of reference columns
	 * @param lo First score column to fill, inclusive
	 * @param hi Last score column to fill, inclusive
	 * @param floor Packed dead-cell score written into blocked cells
	 * @param insertionBarrier Row band that blocks insertions near query edges
	 * @param deletionBarrier Row band that blocks deletions near query edges
	 * @param gapSymbol Reference sentinel that marks a gap column
	 */
	public static void fill(final int[] om, final int[] od, final int[] oi,
			final int[] pm, final int[] pd, final int[] pi, final int[] cm,
			final int[] cd, final int[] ci, final int[] refs, final int[] query,
			final byte[] trace, final int offset, final int diagonal, final int rows,
			final int columns, final int lo, final int hi, final int floor,
			final int insertionBarrier, final int deletionBarrier, final int gapSymbol){
		assert(lo>=1 && hi<cm.length && diagonal-lo<=rows && diagonal-hi>=1) :
			"Diagonal bounds must address valid query rows and score columns";
		int col=lo;
		for(; col+7<=hi; col+=8){
			final IntVector pos=LANES.add(col), row=C(diagonal).sub(pos);
			final IntVector ref=IntVector.fromArray(S, refs, col),
					r0=IntVector.fromArray(S, refs, col-1),
					q=IntVector.fromArray(S, query, rows-diagonal+col),
					q0=IntVector.fromArray(S, query, rows-diagonal+col+1);
			final VectorMask<Integer> gap=ref.eq(gapSymbol),
					match=q.eq(ref).and(ref.compare(NE, 'N')),
					prior=q0.eq(r0).and(r0.compare(NE, 'N'));
			final IntVector packed=IntVector.fromArray(S, om, col-1),
					a=packed.and(SCOREMASK),
					b=IntVector.fromArray(S, od, col-1).and(SCOREMASK),
					c=IntVector.fromArray(S, oi, col-1).and(SCOREMASK),
					streak=packed.and(TIMEMASK);
			IntVector sub=C(POINTSoff_SUB3).blend(POINTSoff_SUB2,
					streak.compare(LT, LIMIT_FOR_COST_3)).blend(POINTSoff_SUB, streak.eq(0));
			sub=sub.blend(C(POINTSoff_SUB).blend(POINTSoff_SUBR, streak.compare(LE, 1)),
					prior).blend(POINTSoff_NOCALL, q.eq('N').or(ref.eq('N')));
			final IntVector x=a.add(sub.blend(C(POINTSoff_MATCH).blend(POINTSoff_MATCH2,
					prior), match)),
					gain=C(POINTSoff_SUB).blend(POINTSoff_MATCH, match),
					y=b.add(gain), z=c.add(gain);
			IntVector time=C(1).blend(streak.add(1),
					match.eq(prior).and(x.compare(GE, y)).and(x.compare(GE, z)));
			time=time.blend(MAX_TIME-MASK5, time.compare(GT, MAX_TIME)).blend(0, gap);
			x.max(y).max(z).or(time).blend(floor, gap).intoArray(cm, col);
			IntVector flags=C(2).blend(1, b.compare(GE, c))
					.blend(0, a.compare(GE, b).and(a.compare(GE, c)))
					.blend(0, time.compare(GT, 1));
			final IntVector lm=IntVector.fromArray(S, pm, col-1).and(SCOREMASK),
					ld=IntVector.fromArray(S, pd, col-1),
					ds=ld.and(TIMEMASK), dscore=ld.and(SCOREMASK);
			final IntVector extra=C(0).blend(POINTSoff_DEL_REF_N, ref.eq('N'))
					.blend(POINTSoff_GAP, gap);
			final IntVector dcost=C(0).blend(POINTSoff_DEL5, ds.and(MASK5).eq(0))
					.blend(POINTSoff_DEL4, ds.compare(LT, LIMIT_FOR_COST_5))
					.blend(POINTSoff_DEL3, ds.compare(LT, LIMIT_FOR_COST_4))
					.blend(POINTSoff_DEL2, ds.compare(LT, LIMIT_FOR_COST_3))
					.blend(POINTSoff_DEL, ds.eq(0));
			final IntVector dx=lm.add(POINTSoff_DEL).add(extra), dy=dscore.add(dcost).add(extra);
			VectorMask<Integer> blocked=row.compare(LT, deletionBarrier)
					.or(row.compare(GT, rows-deletionBarrier));
			IntVector dt=C(1).blend(ds.add(1), dy.compare(GT, dx));
			dt=dt.blend(MAX_TIME-MASK5, dt.compare(GT, MAX_TIME)).blend(0, blocked);
			dx.max(dy).or(dt).blend(floor, blocked).intoArray(cd, col);
			flags=flags.or(C(0).blend(4, dt.compare(GT, 1).or(dscore.compare(GT, lm))));
			final IntVector um=IntVector.fromArray(S, pm, col).and(SCOREMASK),
					ui=IntVector.fromArray(S, pi, col),
					is=ui.and(TIMEMASK), iscore=ui.and(SCOREMASK);
			final IntVector icost=C(POINTSoff_INS4)
					.blend(POINTSoff_INS3, is.compare(LT, LIMIT_FOR_COST_4))
					.blend(POINTSoff_INS2, is.compare(LT, LIMIT_FOR_COST_3))
					.blend(POINTSoff_INS, is.eq(0));
			final IntVector ix=um.add(POINTSoff_INS), iy=iscore.add(icost);
			blocked=gap.or(row.compare(LT, insertionBarrier).and(pos.compare(GT, 1)))
					.or(row.compare(GT, rows-insertionBarrier).and(pos.compare(LT, columns-1)));
			IntVector it=C(1).blend(is.add(1), iy.compare(GT, ix));
			it=it.blend(MAX_TIME-MASK5, it.compare(GT, MAX_TIME)).blend(0, blocked);
			ix.max(iy).or(it).blend(floor, blocked).intoArray(ci, col);
			flags=flags.or(C(0).blend(8, it.compare(GT, 1).or(iscore.compare(GT, um))));
			((ByteVector)flags.convertShape(I2B, ByteVector.SPECIES_64, 0))
					.intoArray(trace, offset+col);
		}
		for(; col<=hi; col++){
			cell(om, od, oi, pm, pd, pi, cm, cd, ci, refs, query, trace, offset,
					diagonal, rows, columns, floor, col, insertionBarrier, deletionBarrier, gapSymbol);
		}
	}

	/**
	 * Fills one column in the scalar remainder after the final full vector block.
	 *
	 * @param om Source match-state antidiagonal for match/substitution transitions
	 * @param od Source deletion-state antidiagonal for match/substitution transitions
	 * @param oi Source insertion-state antidiagonal for match/substitution transitions
	 * @param pm Prior match-state antidiagonal for deletion/insertion transitions
	 * @param pd Prior deletion-state antidiagonal for deletion extensions
	 * @param pi Prior insertion-state antidiagonal for insertion extensions
	 * @param cm Current match-state output
	 * @param cd Current deletion-state output
	 * @param ci Current insertion-state output
	 * @param refs Reference bases indexed by MSA column
	 * @param query Query bases indexed by MSA row in reverse row order
	 * @param trace Compact traceback flags, addressed as {@code offset+col}
	 * @param offset Offset from score columns to the first trace byte
	 * @param diagonal Antidiagonal number; each row is {@code diagonal-col}
	 * @param rows Number of query rows
	 * @param columns Number of reference columns
	 * @param floor Packed dead-cell score written into blocked cells
	 * @param col Score column to fill
	 * @param insertionBarrier Row band that blocks insertions near query edges
	 * @param deletionBarrier Row band that blocks deletions near query edges
	 * @param gapSymbol Reference sentinel that marks a gap column
	 */
	private static void cell(final int[] om, final int[] od, final int[] oi,
			final int[] pm, final int[] pd, final int[] pi, final int[] cm,
			final int[] cd, final int[] ci, final int[] refs, final int[] query,
			final byte[] trace, final int offset, final int diagonal, final int rows,
			final int columns, final int floor, final int col, final int insertionBarrier,
			final int deletionBarrier, final int gapSymbol){
		final int row=diagonal-col, q=query[rows-row], q0=query[rows-row+1],
				r=refs[col], r0=refs[col-1];
		final boolean gap=r==gapSymbol, match=q==r && r!='N', prior=q0==r0 && r0!='N';
		final int a=om[col-1]&SCOREMASK, b=od[col-1]&SCOREMASK,
				c=oi[col-1]&SCOREMASK, s=om[col-1]&TIMEMASK;
		int value=floor;
		if(!gap){
			final int gain=match ? (prior ? POINTSoff_MATCH2 : POINTSoff_MATCH) :
					q=='N' || r=='N' ? POINTSoff_NOCALL :
					prior ? (s<=1 ? POINTSoff_SUBR : POINTSoff_SUB) : POINTSoff_SUB_ARRAY[s+1];
			final int x=a+gain, y=b+(match ? POINTSoff_MATCH : POINTSoff_SUB),
					z=c+(match ? POINTSoff_MATCH : POINTSoff_SUB);
			int t=x>=y && x>=z && match==prior ? s+1 : 1;
			if(t>MAX_TIME){t=MAX_TIME-MASK5;}
			value=Math.max(x, Math.max(y, z))|t;
		}
		cm[col]=value;
		int flags=(value&TIMEMASK)>1 ? 0 : a>=b && a>=c ? 0 : b>=c ? 1 : 2;
		final int lm=pm[col-1]&SCOREMASK, ld=pd[col-1]&SCOREMASK, ds=pd[col-1]&TIMEMASK;
		value=floor;
		if(row>=deletionBarrier && row<=rows-deletionBarrier){
			final int extra=r=='N' ? POINTSoff_DEL_REF_N : gap ? POINTSoff_GAP : 0;
			final int cost=ds==0 ? POINTSoff_DEL : ds<LIMIT_FOR_COST_3 ? POINTSoff_DEL2 :
					ds<LIMIT_FOR_COST_4 ? POINTSoff_DEL3 :
					ds<LIMIT_FOR_COST_5 ? POINTSoff_DEL4 :
					(ds&MASK5)==0 ? POINTSoff_DEL5 : 0;
			final int x=lm+POINTSoff_DEL+extra, y=ld+cost+extra;
			int t=x>=y ? 1 : ds+1;
			if(t>MAX_TIME){t=MAX_TIME-MASK5;}
			value=Math.max(x, y)|t;
		}
		cd[col]=value;
		flags|=(value&TIMEMASK)>1 || ld>lm ? 4 : 0;
		final int um=pm[col]&SCOREMASK, ui=pi[col]&SCOREMASK, is=pi[col]&TIMEMASK;
		value=floor;
		if(!gap && !(row<insertionBarrier && col>1) &&
				!(row>rows-insertionBarrier && col<columns-1)){
			final int x=um+POINTSoff_INS, y=ui+POINTSoff_INS_ARRAY[is+1];
			int t=x>=y ? 1 : is+1;
			if(t>MAX_TIME){t=MAX_TIME-MASK5;}
			value=Math.max(x, y)|t;
		}
		ci[col]=value;
		flags|=(value&TIMEMASK)>1 || ui>um ? 8 : 0;
		trace[offset+col]=(byte)flags;
	}

	/**
	 * Broadcasts a scalar value into every integer lane.
	 *
	 * @param x Value to broadcast
	 * @return Vector with {@code x} in every lane
	 */
	private static IntVector C(final int x){return IntVector.broadcast(S, x);}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** 256-bit integer vector species; all loop bounds assume its eight lanes. */
	private static final VectorSpecies<Integer> S=IntVector.SPECIES_256;
	/** Lane indices added to the first score column of the active vector block. */
	private static final IntVector LANES=IntVector.fromArray(S, new int[]{0, 1, 2, 3, 4, 5, 6, 7}, 0);
}
