package simd;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.ByteVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorSpecies;
import static jdk.incubator.vector.VectorOperators.*;
import static align2.MultiStateAligner11ts.*;

/** All three affine states vectorized across an antidiagonal;contiguous trace stores.
 * @author Collei */
public final class SIMDMSADiagonal {
	public static void fill(int[] om,int[] od,int[] oi,int[] pm,int[] pd,int[] pi,int[] cm,int[] cd,int[] ci,int[] refs,int[] query,byte[] trace,int offset,int diagonal,int rows,int columns,int lo,int hi,int floor,int insertionBarrier,int deletionBarrier,int gapSymbol){
		assert(lo>=1 && hi<cm.length && diagonal-lo<=rows && diagonal-hi>=1) : "Diagonal bounds must address valid query rows and score columns";
		int col=lo;
		for(;col+7<=hi;col+=8){
			final IntVector pos=LANES.add(col),row=C(diagonal).sub(pos);
			final IntVector ref=IntVector.fromArray(S,refs,col),r0=IntVector.fromArray(S,refs,col-1),q=IntVector.fromArray(S,query,rows-diagonal+col),q0=IntVector.fromArray(S,query,rows-diagonal+col+1);
			final VectorMask<Integer> gap=ref.eq(gapSymbol),match=q.eq(ref).and(ref.compare(NE,'N')),prior=q0.eq(r0).and(r0.compare(NE,'N'));
			final IntVector packed=IntVector.fromArray(S,om,col-1),a=packed.and(SCOREMASK),b=IntVector.fromArray(S,od,col-1).and(SCOREMASK),c=IntVector.fromArray(S,oi,col-1).and(SCOREMASK),streak=packed.and(TIMEMASK);
			IntVector sub=C(POINTSoff_SUB3).blend(POINTSoff_SUB2,streak.compare(LT,LIMIT_FOR_COST_3)).blend(POINTSoff_SUB,streak.eq(0));
			sub=sub.blend(C(POINTSoff_SUB).blend(POINTSoff_SUBR,streak.compare(LE,1)),prior).blend(POINTSoff_NOCALL,q.eq('N').or(ref.eq('N')));
			final IntVector x=a.add(sub.blend(C(POINTSoff_MATCH).blend(POINTSoff_MATCH2,prior),match)),gain=C(POINTSoff_SUB).blend(POINTSoff_MATCH,match),y=b.add(gain),z=c.add(gain);
			IntVector time=C(1).blend(streak.add(1),match.eq(prior).and(x.compare(GE,y)).and(x.compare(GE,z)));
			time=time.blend(MAX_TIME-MASK5,time.compare(GT,MAX_TIME)).blend(0,gap);
			x.max(y).max(z).or(time).blend(floor,gap).intoArray(cm,col);
			IntVector flags=C(2).blend(1,b.compare(GE,c)).blend(0,a.compare(GE,b).and(a.compare(GE,c))).blend(0,time.compare(GT,1));
			final IntVector lm=IntVector.fromArray(S,pm,col-1).and(SCOREMASK),ld=IntVector.fromArray(S,pd,col-1),ds=ld.and(TIMEMASK),dscore=ld.and(SCOREMASK);
			final IntVector extra=C(0).blend(POINTSoff_DEL_REF_N,ref.eq('N')).blend(POINTSoff_GAP,gap);
			final IntVector dcost=C(0).blend(POINTSoff_DEL5,ds.and(MASK5).eq(0)).blend(POINTSoff_DEL4,ds.compare(LT,LIMIT_FOR_COST_5)).blend(POINTSoff_DEL3,ds.compare(LT,LIMIT_FOR_COST_4)).blend(POINTSoff_DEL2,ds.compare(LT,LIMIT_FOR_COST_3)).blend(POINTSoff_DEL,ds.eq(0));
			final IntVector dx=lm.add(POINTSoff_DEL).add(extra),dy=dscore.add(dcost).add(extra);
			VectorMask<Integer> blocked=row.compare(LT,deletionBarrier).or(row.compare(GT,rows-deletionBarrier));
			IntVector dt=C(1).blend(ds.add(1),dy.compare(GT,dx));dt=dt.blend(MAX_TIME-MASK5,dt.compare(GT,MAX_TIME)).blend(0,blocked);
			dx.max(dy).or(dt).blend(floor,blocked).intoArray(cd,col);
			flags=flags.or(C(0).blend(4,dt.compare(GT,1).or(dscore.compare(GT,lm))));
			final IntVector um=IntVector.fromArray(S,pm,col).and(SCOREMASK),ui=IntVector.fromArray(S,pi,col),is=ui.and(TIMEMASK),iscore=ui.and(SCOREMASK);
			final IntVector icost=C(POINTSoff_INS4).blend(POINTSoff_INS3,is.compare(LT,LIMIT_FOR_COST_4)).blend(POINTSoff_INS2,is.compare(LT,LIMIT_FOR_COST_3)).blend(POINTSoff_INS,is.eq(0));
			final IntVector ix=um.add(POINTSoff_INS),iy=iscore.add(icost);
			blocked=gap.or(row.compare(LT,insertionBarrier).and(pos.compare(GT,1))).or(row.compare(GT,rows-insertionBarrier).and(pos.compare(LT,columns-1)));
			IntVector it=C(1).blend(is.add(1),iy.compare(GT,ix));it=it.blend(MAX_TIME-MASK5,it.compare(GT,MAX_TIME)).blend(0,blocked);
			ix.max(iy).or(it).blend(floor,blocked).intoArray(ci,col);
			flags=flags.or(C(0).blend(8,it.compare(GT,1).or(iscore.compare(GT,um))));
			((ByteVector)flags.convertShape(I2B,ByteVector.SPECIES_64,0)).intoArray(trace,offset+col);
		}
		for(;col<=hi;col++){cell(om,od,oi,pm,pd,pi,cm,cd,ci,refs,query,trace,offset,diagonal,rows,columns,floor,col,insertionBarrier,deletionBarrier,gapSymbol);}
	}
	private static void cell(int[] om,int[] od,int[] oi,int[] pm,int[] pd,int[] pi,int[] cm,int[] cd,int[] ci,int[] refs,int[] query,byte[] trace,int offset,int diagonal,int rows,int columns,int floor,int col,int insertionBarrier,int deletionBarrier,int gapSymbol){
		final int row=diagonal-col,q=query[rows-row],q0=query[rows-row+1],r=refs[col],r0=refs[col-1];
		final boolean gap=r==gapSymbol,match=q==r && r!='N',prior=q0==r0 && r0!='N';
		final int a=om[col-1]&SCOREMASK,b=od[col-1]&SCOREMASK,c=oi[col-1]&SCOREMASK,s=om[col-1]&TIMEMASK;
		int value=floor;
		if(!gap){final int gain=match?(prior?POINTSoff_MATCH2:POINTSoff_MATCH):q=='N'||r=='N'?POINTSoff_NOCALL:prior?(s<=1?POINTSoff_SUBR:POINTSoff_SUB):POINTSoff_SUB_ARRAY[s+1];
			int x=a+gain,y=b+(match?POINTSoff_MATCH:POINTSoff_SUB),z=c+(match?POINTSoff_MATCH:POINTSoff_SUB),t=x>=y && x>=z && match==prior?s+1:1;
			if(t>MAX_TIME){t=MAX_TIME-MASK5;}value=Math.max(x,Math.max(y,z))|t;}
		cm[col]=value;int flags=(value&TIMEMASK)>1?0:a>=b && a>=c?0:b>=c?1:2;
		final int lm=pm[col-1]&SCOREMASK,ld=pd[col-1]&SCOREMASK,ds=pd[col-1]&TIMEMASK;
		value=floor;
		if(row>=deletionBarrier && row<=rows-deletionBarrier){int extra=r=='N'?POINTSoff_DEL_REF_N:gap?POINTSoff_GAP:0;
			int cost=ds==0?POINTSoff_DEL:ds<LIMIT_FOR_COST_3?POINTSoff_DEL2:ds<LIMIT_FOR_COST_4?POINTSoff_DEL3:ds<LIMIT_FOR_COST_5?POINTSoff_DEL4:(ds&MASK5)==0?POINTSoff_DEL5:0;
			int x=lm+POINTSoff_DEL+extra,y=ld+cost+extra,t=x>=y?1:ds+1;if(t>MAX_TIME){t=MAX_TIME-MASK5;}value=Math.max(x,y)|t;}
		cd[col]=value;flags|=(value&TIMEMASK)>1 || ld>lm?4:0;
		final int um=pm[col]&SCOREMASK,ui=pi[col]&SCOREMASK,is=pi[col]&TIMEMASK;
		value=floor;
		if(!gap && !(row<insertionBarrier && col>1) && !(row>rows-insertionBarrier && col<columns-1)){
			int x=um+POINTSoff_INS,y=ui+POINTSoff_INS_ARRAY[is+1],t=x>=y?1:is+1;if(t>MAX_TIME){t=MAX_TIME-MASK5;}value=Math.max(x,y)|t;}
		ci[col]=value;flags|=(value&TIMEMASK)>1 || ui>um?8:0;trace[offset+col]=(byte)flags;
	}
	private static IntVector C(int x){return IntVector.broadcast(S,x);}
	private static final VectorSpecies<Integer> S=IntVector.SPECIES_256;
	private static final IntVector LANES=IntVector.fromArray(S,new int[]{0,1,2,3,4,5,6,7},0);
}
