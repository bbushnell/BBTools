package test.idaligner;

import idaligner.*;

import java.util.Arrays;

/** Virtual padding removal must conserve real bases and both coordinate axes.
 * @author Raiden */
public final class ContigEdgeAlignmentTest {
	public static void main(String[] args){
		final byte[] model=EndClippedAlignerTest.random(200,931);
		for(boolean scrabble:new boolean[]{false,true}){
			for(int[] span:new int[][]{{0,199},{0,149},{40,199},{40,149}}){
				final byte[] ref=Arrays.copyOfRange(model,span[0],span[1]+1);
				final EndClippedAligner.Result r=ContigEdgeAlignment.padded(model,ref,span[0]>0,span[1]<199,1,scrabble);
				check(r!=null && r.qStart==span[0] && r.qStop==span[1] && r.rStart==0 && r.rStop==ref.length-1 && r.identity==1f,
					"Normal padded alignment must recover exact physical fragments without artificial-N confidence: scrabble="+scrabble);
				EndClippedAlignerTest.replay(model,ref,r);
			}
		}
		// Three virtual Ns, then one genuine genomic N adjacent to them.
		AlignmentStats raw=trace("NNNNmm",0,6,6);
		EndClippedAligner.Result r=ContigEdgeAlignment.observed(raw,6,3,3,true,false);
		check(r.qStart==3 && r.qStop==5 && r.ns==1 && r.matches==2 && r.identity==2.5f/3 && r.physicalClips==1,
			"Real genomic Ns must survive coordinate cropping even when adjacent to indistinguishable virtual-N operations");
		raw=trace("NImDmN",0,5,5);r=ContigEdgeAlignment.observed(raw,5,1,3,true,true);
		check(r.qStart==2 && r.qStop==3 && r.rStart==0 && r.rStop==2 && r.dels==1 && r.ins==0 && r.physicalClips==3,
			"Query insertions at the virtual junction are excluded; internal reference deletions preserve their distinct axis consumption");
		raw=trace("NmImDN",0,5,5);r=ContigEdgeAlignment.observed(raw,5,1,3,true,true);
		check(r!=null && r.ins==1 && r.dels==1 && r.qStart==1 && r.qStop==3 && r.rStop==2,
			"Internal query insertions remain observed after removing virtual padding at both sides");
		check(ContigEdgeAlignment.observed(trace("NN",0,2,5),2,2,3,true,false)==null,
			"A positive or neutral path entirely in artificial padding cannot become a genomic call");
		r=ContigEdgeAlignment.observed(trace("Imm",2,3,6),3,0,4,false,true);
		check(r!=null && r.qStart==0 && r.qStop==2 && r.ins==1 && r.physicalClips==0,
			"A query deletion anchored inside real sequence must remain observed, not become fictitious physical clipping");
		r=ContigEdgeAlignment.observed(trace("mmI",0,3,6),3,0,4,false,true);
		check(r!=null && r.qStop==2 && r.ins==1 && r.rStop==1,
			"Trailing query consumption inside real genomic context must survive padding removal");
		final EndClippedAligner.Result accepted=new EndClippedAligner().alignContigEdges(model,Arrays.copyOfRange(model,40,150),true,true);
		final ContigEdgeAlignment.Bounds bounds=ContigEdgeAlignment.refineInterior(model,Arrays.copyOfRange(model,40,150),accepted,10);
		check(bounds.start==accepted.rStart && bounds.stop==accepted.rStop && !bounds.changed,
			"A fragment missing both model ends has no interior endpoint to refine");
		interiorWithLongQuery(model);
		interiorBeyondCandidate(model);
		System.out.println("PASS padded normal alignment, observed-N identity, junction I/D coordinates, padding-only rejection and both-edge refinement exemption");
	}
	static void interiorWithLongQuery(byte[] model){
		final byte[] ref=EndClippedAlignerTest.random(150,1931);
		System.arraycopy(model,40,ref,0,60);System.arraycopy(model,140,ref,60,60);
		final String middle=EndClippedAlignerTest.repeat('m',60)+EndClippedAlignerTest.repeat('I',40)+EndClippedAlignerTest.repeat('m',60);
		final EndClippedAligner.Result left=EndClippedAlignerTest.traced(middle+"DDDDDDDD",1);left.rLen=ref.length;
		final ContigEdgeAlignment.Bounds a=ContigEdgeAlignment.refineInterior(model,ref,left,10);
		check(left.qStop-left.qStart+1>Math.min(ref.length,left.rStop+11),
			"This fixture must exercise a consumed model slice longer than the available refinement reference");
		System.out.println("Contiguous40-base-gap Scrabble bounds="+a.start+"-"+a.stop+" changed="+a.changed);
		check(a.start==0 && a.stop>=a.start && a.stop<ref.length,
			"A bounded-band reanchor must preserve the physical side and valid coordinates even when a long internal gap defeats exact alignment");
		final EndClippedAligner.Result right=EndClippedAlignerTest.traced("DDDDDDDD"+middle,2);
		right.rLen=ref.length;right.rStart=ref.length-1-left.rStop;right.rStop=ref.length-1;
		final ContigEdgeAlignment.Bounds b=ContigEdgeAlignment.refineInterior(reverse(model),reverse(ref),right,10);
		check(b.start>=0 && b.start<=b.stop && b.stop==ref.length-1,
			"The opposite physical endpoint must remain fixed even when heuristic alignment differs by direction");
		// Distributed single-base deletions exercise a long query without requiring
		// Scrabble's bounded drift to cross one40-base gap as a single event.
		final StringBuilder ops=new StringBuilder();int position=0;
		for(int k=0;k<160;k++){
			if(k<120 && k%3==2){ops.append('I');}else{ref[position++]=model[40+k];ops.append('m');}
		}
		check(position==120,"The explicit deletion fixture must remove40 model bases");
		final EndClippedAligner.Result distributedLeft=EndClippedAlignerTest.traced(ops+"DDDDDDDD",1);distributedLeft.rLen=ref.length;
		final ContigEdgeAlignment.Bounds da=ContigEdgeAlignment.refineInterior(model,ref,distributedLeft,10);
		final EndClippedAligner.Result distributedRight=EndClippedAlignerTest.traced("DDDDDDDD"+new StringBuilder(ops).reverse(),2);
		distributedRight.rLen=ref.length;distributedRight.rStart=22;distributedRight.rStop=149;
		final ContigEdgeAlignment.Bounds db=ContigEdgeAlignment.refineInterior(reverse(model),reverse(ref),distributedRight,10);
		check(da.start==0 && da.stop==119 && da.changed && db.start==30 && db.stop==149 && db.changed,
			"Distributed-indel model longer than reference must reanchor the interior end on both sides: "+da.start+"-"+da.stop+" / "+db.start+"-"+db.stop);
		final byte[] lower=ref.clone();for(int i=0;i<lower.length;i++){lower[i]=(byte)Character.toLowerCase((char)lower[i]);}
		final ContigEdgeAlignment.Bounds low=ContigEdgeAlignment.refineInterior(model,lower,distributedLeft,10);
		check(low.start==da.start && low.stop==da.stop,"Soft-masked genomic bases must use the same uppercase private refinement context");
		final byte[] full=EndClippedAlignerTest.random(190,2931);System.arraycopy(reverse(ref),0,full,40,ref.length);
		final ContigEdgeAlignment.Bounds translated=ContigEdgeAlignment.refineInterior(reverse(model),full,distributedRight,10,40);
		check(translated.start==40+db.start && translated.stop==189,
			"Window-relative alignment must refine in the original full-contig frame without moving its physical side");
	}
	static void interiorBeyondCandidate(byte[] model){
		final byte[] full=EndClippedAlignerTest.random(190,3931);
		System.arraycopy(model,0,full,32,158);
		// The original window begins at40, eight bases into the real model span.
		final EndClippedAligner.Result right=EndClippedAlignerTest.traced("IIIIIIII"+EndClippedAlignerTest.repeat('m',150),2);
		right.qLen=model.length;
		final ContigEdgeAlignment.Bounds a=ContigEdgeAlignment.refineInterior(model,full,right,10,40);
		check(a.start==32 && a.start<40 && a.stop==189 && a.changed,
			"Full-contig refinement must recover the interior endpoint outside its original candidate window while fixing the physical right end: "+a.start+"-"+a.stop);
		final EndClippedAligner.Result left=EndClippedAlignerTest.traced(EndClippedAlignerTest.repeat('m',150)+"IIIIIIII",1);
		left.qLen=model.length;left.qStart=42;left.qStop=199;
		final ContigEdgeAlignment.Bounds b=ContigEdgeAlignment.refineInterior(reverse(model),reverse(full),left,10,0);
		check(b.start==0 && b.stop==157 && b.stop>=left.rLen && b.changed,
			"Mirrored full-contig refinement must extend beyond the original window without moving the physical left end: "+b.start+"-"+b.stop);
	}
	static byte[] reverse(byte[] source){final byte[] r=source.clone();for(int i=0;i<source.length;i++){r[i]=source[source.length-1-i];}return r;}
	static AlignmentStats trace(String ops,int start,int qLength,int rLength){
		final AlignmentStats r=new AlignmentStats(true);r.qLen=qLength;r.rLen=rLength;r.rStart=start;
		r.setFromMatchString(ops.getBytes(java.nio.charset.StandardCharsets.US_ASCII));r.rStop=start+r.matches+r.subs+r.dels+r.ns-1;return r;
	}
	static void check(boolean ok,String message){if(!ok){throw new AssertionError(message);}}
}
