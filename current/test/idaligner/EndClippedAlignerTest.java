package test.idaligner;

import idaligner.*;

import java.util.Arrays;
import java.util.Random;

/** Coordinate, score, operation and scratch-reuse fixtures for free model ends.
 * @author Raiden
 */
public final class EndClippedAlignerTest {
	public static void main(String[] args){
		final EndClippedAligner aligner=new EndClippedAligner();
		final byte[] query=random(1800, 419);
		for(int[] bounds:new int[][]{{0,1799},{0,1199},{400,1799},{400,1315},{400,1314}}){
			final int from=bounds[0], to=bounds[1], len=to-from+1;
			final byte[] ref=padded(Arrays.copyOfRange(query, from, to+1), 23, 41);
			final EndClippedAligner.Result r=aligner.align(query, ref);
			check(r!=null && r.qStart==from && r.qStop==to && r.rStart==23 && r.rStop==22+len,
				"Clipped query and target coordinates must be independently correct");
			check(r.identity==1f && r.score==len && r.matchString.length==len && r.qLen==1800,
				"Clipped bases must not enter identity or traceback; original model length is retained");
			replay(query, ref, r);
		}
		final byte[] exact=random(200, 31), insertion=new byte[201];
		System.arraycopy(exact,0,insertion,0,80);insertion[80]='N';System.arraycopy(exact,80,insertion,81,120);
		EndClippedAligner.Result r=aligner.align(exact,insertion);
		check(r.score==199 && r.dels==1 && r.matches==200 && r.identity==200f/201,
			"Internal target insertion is penalized, not clipped");replay(exact,insertion,r);
		r=aligner.align(insertion,exact);
		check(r.score==199 && r.ins==1 && r.matches==200,"Internal model insertion must consume query only");replay(insertion,exact,r);
		final byte[] unknown=exact.clone();unknown[90]='N';r=aligner.align(exact,unknown);
		check(r.score==199 && r.ns==1 && r.identity==199.5f/200,"Internal N is neutral with half identity credit");replay(exact,unknown,r);
		check(aligner.align(new byte[0],exact)==null && aligner.align(new byte[]{'A'},new byte[]{'C'})==null
			&& aligner.align(new byte[]{'N'},new byte[]{'N'})==null,"Empty/no-positive paths must not expose false endpoints");
		final byte[] twice=new byte[400];System.arraycopy(exact,0,twice,0,200);System.arraycopy(exact,0,twice,200,200);
		r=aligner.align(exact,twice);check(r.rStart==0 && r.rStop==199,"Equal local maxima retain the first row-major endpoint");
		// A better short local hit can mask a weaker long one. This is a known
		// one-winner limitation, not proof that every >=916 alignment failed.
		final byte[] competing=new byte[4000];Arrays.fill(competing,(byte)'N');
		System.arraycopy(query,0,competing,0,700);System.arraycopy(query,800,competing,3000,1000);
		for(int p=3004;p<4000;p+=5){competing[p]=competing[p]=='A'?(byte)'C':(byte)'A';}
		r=aligner.align(query,competing);
		check(r!=null && r.rStop-r.rStart+1<916,"Best short hit is returned even with a weaker known 1000-base homolog elsewhere");
		weightedEdges(aligner, exact);
		localSupport();
		ContigEdgeAlignmentTest.main(new String[0]);
		System.out.println("PASS EndClippedAlignerTest full/prefix/suffix/both/minlength_edges/indels/N/nohit/ties/scratch_reuse");
	}
	static void localSupport(){
		final String good=repeat('m',80), gap=repeat('N',100), weak=repeat('S',100);
		final EndClippedAligner.Result left=traced(good+weak,1), right=traced(weak+good,2);
		final EndClippedAligner.Result a=EndClippedAligner.trimEdgeSupport(left,true,false,20,.6f);
		final EndClippedAligner.Result b=EndClippedAligner.trimEdgeSupport(right,false,true,20,.6f);
		check(a!=null && b!=null && a.rStop-a.rStart==b.rStop-b.rStart
			&& a.rStart==0 && b.rStop==right.rLen-1 && a.physicalClips==1 && b.physicalClips==2,
			"The same local-support rule must retain mirrored physical ends and never invent an interior truncation flag");
		check(a.rStop<100 && b.rStart>80 && a.identity>.8f && b.identity>.8f,
			"A low-support extension must stop before the distant weak region can inflate the accepted span");
		final EndClippedAligner.Result exact=traced(good,1), shortExact=traced("mmmm",2);
		check(EndClippedAligner.trimEdgeSupport(exact,true,false,20,.6f)==exact
			&& EndClippedAligner.trimEdgeSupport(shortExact,false,true,64,.6f)==shortExact,
			"Supported exact fragments, including shorter-than-window controls, are unchanged");
		check(EndClippedAligner.trimEdgeSupport(traced(weak,1),true,false,20,.6f)==null,
			"No supported physical anchor cannot become an edge rescue");
		check(EndClippedAligner.trimEdgeSupport(traced(gap+good,1),true,false,20,.54f)==null,
			"Strong interior sequence cannot excuse an all-N physical anchor or create missing local evidence");
		final EndClippedAligner.Result border=traced(repeat('m',12)+repeat('S',8),1);
		check(EndClippedAligner.trimEdgeSupport(border,true,false,20,.6f)==border
			&& EndClippedAligner.trimEdgeSupport(border,true,false,20,.6001f)==null,
			"A window exactly at the declared identity floor passes without multiplication-rounding bias");
		final EndClippedAligner.Result interior=traced(weak+good,0);
		check(EndClippedAligner.trimEdgeSupport(interior,true,true,20,.6f)==interior,
			"Physical contact alone cannot authorize trimming a full model or an interior rescue");
		final EndClippedAligner.Result split=EndClippedAligner.trimEdgeSupport(traced(good+gap+repeat('m',120),3),true,true,20,.6f);
		check(split!=null && split.physicalClips==2 && split.rStart>80 && split.rStop==299,
			"Neutral-N gaps break local support; disconnected edge fragments must not be joined across the gap");
		check(EndClippedAligner.trimEdgeSupport(traced(good+gap+good,3),true,true,20,.6f)==null,
			"Equal supported disconnected ends are ambiguous, never resolved by strand-dependent left preference");
		final EndClippedAligner.Result mirror=EndClippedAligner.trimEdgeSupport(traced(repeat('m',120)+gap+good,3),true,true,20,.6f);
		check(mirror!=null && mirror.physicalClips==1 && mirror.rStop==299-split.rStart,
			"Mirroring unequal disconnected anchors must mirror the selected component");
		final EndClippedAligner.Result indel=traced("SSSID"+repeat('S',50)+"mmID"+good,2);
		final EndClippedAligner.Result cropped=EndClippedAligner.trimEdgeSupport(indel,false,true,20,.6f);
		check(cropped!=null && cropped.qStop==indel.qStop && cropped.rStop==indel.rStop
			&& cropped.qStop-cropped.qStart+1==cropped.matches+cropped.subs+cropped.ins+cropped.ns
			&& cropped.rStop-cropped.rStart+1==cropped.matches+cropped.subs+cropped.dels+cropped.ns,
			"Asymmetric I/D consumption before and inside the retained trace must preserve both coordinate axes");
		System.out.println("PASS boundary-local support mirrored ends, exact/short controls, N breaks, disconnected ends and gap coordinates");
	}
	static String repeat(char c,int n){final char[] a=new char[n];Arrays.fill(a,c);return new String(a);}
	static EndClippedAligner.Result traced(String ops,int physical){
		final EndClippedAligner.Result r=new EndClippedAligner.Result();
		r.setFromMatchString(ops.getBytes(java.nio.charset.StandardCharsets.US_ASCII));
		r.qStart=(physical&1)!=0?40:0;r.qStop=r.qStart+r.matches+r.subs+r.ins+r.ns-1;
		r.qLen=r.qStop+1+((physical&2)!=0?40:0);r.rStart=0;r.rStop=r.matches+r.subs+r.dels+r.ns-1;r.rLen=r.rStop+1;
		r.identity=(r.matches+.5f*r.ns)/ops.length();r.matchReward=46;r.nonmatchPenalty=54;
		r.optimizationScore=r.matches*46-(r.subs+r.ins+r.dels)*54;r.physicalClips=physical;return r;
	}
	static void weightedEdges(EndClippedAligner aligner,byte[] model){
		for(int[] profile:new int[][]{{1,1},{46,54},{2,3}}){
			for(int[] bounds:new int[][]{{0,199},{0,149},{40,199},{40,149}}){
				final byte[] fragment=Arrays.copyOfRange(model,bounds[0],bounds[1]+1);
				final EndClippedAligner.Result r=aligner.alignContigEdges(model,fragment,bounds[0]>0,bounds[1]<199,profile[0],profile[1]);
				check(r!=null && r.qStart==bounds[0] && r.qStop==bounds[1] && r.rStart==0 && r.rStop==fragment.length-1
					&& r.identity==1f && r.score==fragment.length && r.optimizationScore==fragment.length*profile[0],
					"Weighted edge objectives retain exact full/left/right/both fragments and separate elementary from weighted scores");
				replay(model,fragment,r);
			}
			final byte[] insertion=new byte[model.length+1];System.arraycopy(model,0,insertion,0,80);
			insertion[80]='N';System.arraycopy(model,80,insertion,81,model.length-80);
			final EndClippedAligner.Result r=aligner.alignContigEdges(model,insertion,false,false,profile[0],profile[1]);
			check(r!=null && r.matches==model.length && r.dels==1 && r.score==model.length-1
				&& r.optimizationScore==model.length*profile[0]-profile[1],"Weighted internal-gap penalty must conserve original traceback score and operations");
			replay(model,insertion,r);
		}
		final EndClippedAligner.Result old=aligner.alignContigEdges(model,model,true,true), explicit=aligner.alignContigEdges(model,model,true,true,1,1);
		check(old.score==explicit.score && old.optimizationScore==explicit.optimizationScore && Arrays.equals(old.matchString,explicit.matchString),
			"Legacy physical-edge overload remains exactly the elementary objective");
		final byte[] longModel=random(3405,719), noisy=random(1159,1219);
		System.arraycopy(longModel,0,noisy,1055,104);
		final EndClippedAligner.Result supported=aligner.alignContigEdges(longModel,noisy,true,true,46,54);
		check(supported!=null && supported.qStart==0 && supported.qStop==103 && supported.rStart==1055 && supported.rStop==1158,
			"A real model-prefix fragment at the contig end must beat extending through unrelated upstream sequence");
		replay(longModel,noisy,supported);
	}
	/** Replay actual bases, not the DP recurrence, to catch swapped axes and clipping counted as gaps. */
	static void replay(byte[] query,byte[] ref,EndClippedAligner.Result r){
		int q=r.qStart,t=r.rStart,score=0;
		for(byte op:r.matchString){
			if(op=='I'){q++;score--;}
			else if(op=='D'){t++;score--;}
			else{
				final byte a=query[q++],b=ref[t++];
				check(op=='N' ? a=='N'||b=='N' : op=='m' ? a==b : a!=b,"Trace symbol must describe original bases");
				score+=op=='m'?1:op=='S'?-1:0;
			}
		}
		check(q==r.qStop+1 && t==r.rStop+1 && score==r.score,"Trace replay must end at both reported endpoints with the same score");
	}
	static byte[] random(int n,long seed){final Random r=new Random(seed);final byte[] b=new byte[n],a={'A','C','G','T'};for(int i=0;i<n;i++){b[i]=a[r.nextInt(4)];}return b;}
	static byte[] padded(byte[] b,int left,int right){final byte[] out=new byte[left+b.length+right];Arrays.fill(out,(byte)'N');System.arraycopy(b,0,out,left,b.length);return out;}
	static void check(boolean value,String why){if(!value){throw new AssertionError(why);}}
}
