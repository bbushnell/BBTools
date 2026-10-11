package prok;

import java.util.Arrays;
import java.util.Random;
import stream.Read;

/** Real refinement checks: final inclusive span, no observer dependency,
 * successful fallback, rejection and non16S isolation on both aligner routes.
 * @author Raiden */
public final class GeneCallerSsuScoreTest {
	public static void main(String[] args){
		final boolean oldIda=GeneCaller.useIDAligner,oldFallback=ProkObject.rrnaFallback;
		final Read[] old16=ProkObject.r16SSequence,old18=ProkObject.r18SSequence;
		final float oldMin=ProkObject.min16SIdentity;final int oldStart=ProkObject.ssuStartSlop,oldStop=ProkObject.ssuStopSlop;
		final byte[] gene=new byte[300],bases=new byte[400],alphabet={'A','C','G','T'};final Random random=new Random(193619);
		for(int i=0;i<bases.length;i++){bases[i]=alphabet[random.nextInt(4)];}System.arraycopy(bases,50,gene,0,gene.length);
		final byte[] bad=new byte[300];Arrays.fill(bad,(byte)'A');
		try{
			ProkObject.min16SIdentity=.9f;ProkObject.ssuStartSlop=ProkObject.ssuStopSlop=0;
			for(boolean ida:new boolean[]{false,true}){
				GeneCaller.useIDAligner=ida;ProkObject.rrnaFallback=false;
				ProkObject.r16SSequence=new Read[]{new Read(gene,null,"exact",0)};
				final GeneCaller caller=new GeneCaller(60,1,1,0f,0f,0f,0f,0f,new GeneModel(false));
				Orf orf=orf(bases,ProkObject.r16S);
				check(caller.refineByAlignment(orf,bases,0,new StatsContainer(ProkObject.r16S)),"Perfect model must pass real refinement");
				check(orf.start==50 && orf.stop==349 && orf.orfScore==300f,"Score must use actual refined inclusive span and measured perfect identity with no observer");
				check(orf.ssuFullModel==gene && "exact".equals(orf.ssuFullModelName),"SSU comparison must retain the actual accepted full model without an observer");
				final float[] identity={Float.NaN};
				caller.setAttemptSink((index,total,preStart,preStop,preScore,accepted,postStart,postStop,strand,scaffold,label,length,id,reason,candidate)->{if(accepted){identity[0]=id;}});
				orf=orf(bases,ProkObject.r16S);check(caller.refineByAlignment(orf,bases,0,new StatsContainer(ProkObject.r16S)) && orf.orfScore==300f && identity[0]==1f,"Adding an observer cannot change accepted scoring");
				caller.setAttemptSink(null);ProkObject.rrnaFallback=true;ProkObject.r16SSequence=new Read[]{new Read(bad,null,"reject",0),new Read(gene,null,"exact",1)};
				orf=orf(bases,ProkObject.r16S);check(caller.refineByAlignment(orf,bases,0,new StatsContainer(ProkObject.r16S)) && orf.orfScore==300f,"Accepted fallback identity must replace failed-attempt state");
				check(orf.ssuFullModel==gene && "exact".equals(orf.ssuFullModelName),"Fallback provenance must identify the accepted model, never the rejected universal attempt");
				ProkObject.r16SSequence=new Read[]{new Read(bad,null,"reject",0)};orf=orf(bases,ProkObject.r16S);
				check(!caller.refineByAlignment(orf,bases,0,new StatsContainer(ProkObject.r16S)) && orf.orfScore<=-999f,"Rejected alignment must never receive a positive identity score");
				ProkObject.r18SSequence=new Read[]{new Read(gene,null,"exact",0)};orf=orf(bases,ProkObject.r18S);
				check(caller.refineByAlignment(orf,bases,0,new StatsContainer(ProkObject.r18S)) && orf.orfScore==17f,"Legacy18S is outside the16S scoring change");
				ProkObject.r16SSequence=null;orf=orf(bases,ProkObject.r16S);
				check(caller.refineByAlignment(orf,bases,0,new StatsContainer(ProkObject.r16S)) && orf.orfScore==17f,"No-consensus compatibility path has no identity and retains its existing score");
			}
		}finally{GeneCaller.useIDAligner=oldIda;ProkObject.rrnaFallback=oldFallback;ProkObject.r16SSequence=old16;ProkObject.r18SSequence=old18;ProkObject.min16SIdentity=oldMin;ProkObject.ssuStartSlop=oldStart;ProkObject.ssuStopSlop=oldStop;}
		System.out.println("SSU_SCORE_TEST_PASS both_aligners refined_inclusive_span observer_parity fallback rejection non16S no_consensus");
	}
	static Orf orf(byte[] bases,int type){final Orf r=new Orf("fixture",80,319,0,0,bases,false,type);r.orfScore=17f;return r;}
	static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
