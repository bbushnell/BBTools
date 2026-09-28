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
	public static int fillRow(int[] pm,int[] pd,int[] pi,int[] cm,int[] ci,int[] refs,int[] flags,int q,int q0,int row,int rows,int columns,int floor,int insertionBarrier,int gapSymbol,int from,int to,int target){
		assert(from>=1 && to<=columns && row>=1 && row<=rows) : "Limited vector span must lie within the retained row buffers";
		final int limit=to-7;int col=from;
		for(;col<=limit;col+=8){
			final IntVector positions=LANES.add(col);
			// Every remaining query/reference base can contribute at most MATCH2.
			// Bound by the query before multiplying:wide references can otherwise overflow shifted scores.
			final IntVector cutoff=C(target).sub(C(rows-row).min(C(columns).sub(positions)).mul(POINTSoff_MATCH2));
			final IntVector ref=IntVector.fromArray(S,refs,col),prevRef=IntVector.fromArray(S,refs,col-1);
			final VectorMask<Integer> gap=ref.eq(gapSymbol),match=ref.eq(q).and(ref.compare(NE,'N')),prevMatch=prevRef.eq(q0).and(prevRef.compare(NE,'N'));
			final IntVector diag=IntVector.fromArray(S,pm,col-1),a=diag.and(SCOREMASK),b=IntVector.fromArray(S,pd,col-1).and(SCOREMASK),c=IntVector.fromArray(S,pi,col-1).and(SCOREMASK),streak=diag.and(TIMEMASK);
			IntVector sub=C(POINTSoff_SUB3).blend(POINTSoff_SUB2,streak.compare(LT,LIMIT_FOR_COST_3)).blend(POINTSoff_SUB,streak.eq(0));
			sub=sub.blend(C(POINTSoff_SUB).blend(POINTSoff_SUBR,streak.compare(LE,1)),prevMatch);
			if(q=='N'){sub=C(POINTSoff_NOCALL);}else{sub=sub.blend(POINTSoff_NOCALL,ref.eq('N'));}
			final IntVector gain=sub.blend(C(POINTSoff_MATCH).blend(POINTSoff_MATCH2,prevMatch),match);
			final IntVector fromMS=a.add(gain),otherGain=C(POINTSoff_SUB).blend(POINTSoff_MATCH,match),fromD=b.add(otherGain),fromI=c.add(otherGain);
			final VectorMask<Integer> chooseMS=fromMS.compare(GE,fromD).and(fromMS.compare(GE,fromI));
			IntVector time=C(1).blend(streak.add(1),match.eq(prevMatch).and(chooseMS));
			time=time.blend(MAX_TIME-MASK5,time.compare(GT,MAX_TIME)).blend(0,gap);
			IntVector ms=fromMS.max(fromD).max(fromI).or(time).blend(floor,gap);
			ms.blend(floor,ms.compare(LT,cutoff)).intoArray(cm,col);
			IntVector prevMS=C(2).blend(1,b.compare(GE,c)).blend(0,a.compare(GE,b).and(a.compare(GE,c))).blend(0,time.compare(GT,1));
			final IntVector up=IntVector.fromArray(S,pm,col).and(SCOREMASK),upPacked=IntVector.fromArray(S,pi,col),upI=upPacked.and(SCOREMASK),is=upPacked.and(TIMEMASK);
			final IntVector icost=C(POINTSoff_INS4).blend(POINTSoff_INS3,is.compare(LT,LIMIT_FOR_COST_4)).blend(POINTSoff_INS2,is.compare(LT,LIMIT_FOR_COST_3)).blend(POINTSoff_INS,is.eq(0));
			final IntVector x=up.add(POINTSoff_INS),y=upI.add(icost);
			IntVector itime=C(1).blend(is.add(1),y.compare(GT,x));itime=itime.blend(MAX_TIME-MASK5,itime.compare(GT,MAX_TIME));
			VectorMask<Integer> blocked=gap;
			if(row<insertionBarrier){blocked=blocked.or(positions.compare(GT,1));}
			if(row>rows-insertionBarrier){blocked=blocked.or(positions.compare(LT,columns-1));}
			itime=itime.blend(0,blocked);IntVector ins=x.max(y).or(itime).blend(floor,blocked);
			ins.blend(floor,ins.compare(LT,cutoff)).intoArray(ci,col);
			prevMS.or(C(0).blend(8,itime.compare(GT,1).or(upI.compare(GT,up)))).intoArray(flags,col);
		}
		return col;
	}
	private static IntVector C(int value){return IntVector.broadcast(S,value);}
	private static final VectorSpecies<Integer> S=IntVector.SPECIES_256;
	private static final IntVector LANES=IntVector.fromArray(S,new int[]{0,1,2,3,4,5,6,7},0);
}
