package test.prok;

import prok.*;

import java.util.ArrayList;
import java.util.Arrays;
import consensus.BaseGraph;
import dna.AminoAcid;
import idaligner.AlignmentStats;
import idaligner.EndClippedAligner;
import idaligner.ScrabbleAligner;
import map.LongHashSet;
import stream.Read;

/** Physical-edge rescue must preserve model origin, ordinary verification and endpoints.
 * These are mechanism fixtures, not threshold calibration.
 * @author Raiden */
public final class NcrnaPhysicalEdgeTest {
	public static void main(String[] args){
		final byte[] model=Prok16sPartialRescueTest.random(1000);
		modelCoordinates(model);policy(model);edge18Defaults(model);maskedLegacyNets(model);
		for(boolean reused:new boolean[]{false,true}){
			for(int mode=1; mode<=3; mode++){
				final int from=mode==2?0:200, to=mode==1?1000:850;
				final byte[] bases=Prok16sPartialRescueTest.embed(Arrays.copyOfRange(model,from,to),mode==2?80:0,mode==1?80:0);
				for(int strand=0; strand<2; strand++){
					final NcrnaScavenger off=caller(model,reused,false,null,1f,1f), on=caller(model,reused,true,null,1f,1f);
					check(off.scavenge("edge",bases,strand,new ArrayList<int[]>()).isEmpty(),"Unchanged primary rejects the deliberately stringent clipped fixture");
					final Prok16sPartialRescueTest.Recorder observer=reused?Prok16sPartialRescueTest.observe(on,model):null;
					final ArrayList<Orf> result=on.scavenge("edge",bases,strand,new ArrayList<int[]>());
					check(result.size()==1,"One physical fragment must survive both primary caller shapes; mode="+mode+" reused="+reused);
					final Orf orf=result.get(0);check(orf.ncrnaPartial==mode,"Sense clipping bits must retain both model ends");
					check(orf.start==(mode==2?80:0) && orf.stop==bases.length-1-(mode==1?80:0),"Boundary trim/votes/NN must not move accepted physical clipped ends");
					check(observer==null || observer.count==(mode==3?0:1),"Single-edge partials must process their interior endpoint; both missing ends have no refinement target");
					if(strand==1){orf.flip();}
					final int mask=strand==0?mode:((mode&1)<<1)|((mode&2)>>>1);
					check(orf.toGff().contains(";partial="+(mask==1?"10":mask==2?"01":"11")),"GFF partials must follow genomic orientation");
				}
			}
			final byte[] full=Prok16sPartialRescueTest.embed(model,100,100);
			final NcrnaScavenger off=caller(model,reused,false,null,1f,1f), on=caller(model,reused,true,null,1f,1f);
			final ArrayList<Orf> a=off.scavenge("full",full,0,new ArrayList<int[]>()), b=on.scavenge("full",full,0,new ArrayList<int[]>());
			check(a.size()==1 && Prok16sPartialRescueTest.calls(a).equals(Prok16sPartialRescueTest.calls(b)) && off.alignmentCount()==on.alignmentCount(),"Successful primary output and alignment work must stay identical");
			for(boolean gap:new boolean[]{false,true}){
				final byte[] internal=Prok16sPartialRescueTest.embed(Arrays.copyOfRange(model,200,850),4000,4000);
				if(gap){Arrays.fill(internal,0,4000,(byte)'N');Arrays.fill(internal,4650,internal.length,(byte)'N');}
				final NcrnaScavenger x=caller(model,reused,false,null,1f,1f), y=caller(model,reused,true,null,1f,1f);
				check(x.scavenge("internal",internal,0,new ArrayList<int[]>()).isEmpty() && y.scavenge("internal",internal,0,new ArrayList<int[]>()).isEmpty()
					&& x.alignmentCount()==y.alignmentCount(),"Internal/N-gap candidate ends must not authorize extra clipping");
			}
			final byte[] claimed=Prok16sPartialRescueTest.embed(Arrays.copyOfRange(model,200,1000),0,80);
			final ArrayList<int[]> claims=new ArrayList<int[]>();claims.add(new int[]{0,9});
			check(caller(model,reused,true,null,1f,1f).scavenge("claimed",claimed,0,claims).isEmpty(),"Claim subtraction cannot create a new physical left boundary");
		}
		final BaseGraph graph=graph(model);
		final byte[] internal18=Prok16sPartialRescueTest.embed(Arrays.copyOfRange(model,200,850),4000,4000);
		// Random flanks extend the optimal local path and lower identity below .99.
		// Neutral N flanks isolate free model-end clipping from incidental matches.
		Arrays.fill(internal18,0,4000,(byte)'N');Arrays.fill(internal18,4650,internal18.length,(byte)'N');
		final NcrnaScavenger legacy18=caller(model,true,false,null,1f,1f), edge18=caller(model,true,true,null,1f,1f);
		for(NcrnaScavenger s:new NcrnaScavenger[]{legacy18,edge18}){s.family=Euk18sRuntimeConfig.FAMILY;s.outputType=ProkObject.r18S;s.modelClipRescue=true;s.modelClipRescueId=.99f;}
		final int legacyCount=legacy18.scavenge("internal18",internal18,0,new ArrayList<int[]>()).size();
		final int edgeCount=edge18.scavenge("internal18",internal18,0,new ArrayList<int[]>()).size();
		check(legacyCount==1 && edgeCount==0,"TARGET must suppress the actual18S unrestricted interior rescue; legacy="+legacyCount+" target="+edgeCount);
		final byte[] mutated=Arrays.copyOfRange(model,200,850);for(int i=40;i<mutated.length-10;i+=50){mutated[i]=(byte)(mutated[i]=='A'?'C':'A');}
		final NcrnaScavenger hbm=caller(model,false,true,graph,.999f,.9f);
		final ArrayList<Orf> rescued=hbm.scavenge("hbm",mutated,0,new ArrayList<int[]>());
		check(rescued.size()==1 && rescued.get(0).ncrnaPartial==3 && hbm.hbmScoreCalls()>0,"Physical-edge HBM rescue must work without changing legacy primary verification");
		its();
		System.out.println("PASS physical left/right/both edges and strands, internal/N-gap/claim negatives, unchanged primary output/work, model-origin HBM rescue, endpoint bypass, GFF and ITS coordinates");
	}
	static void maskedLegacyNets(byte[] model){
		final int[] offsets={-1,0,1};
		for(int mask=1;mask<=3;mask++){
			final LengthNet active=new LengthNet();
			final NcrnaBoundaryScorer.RefinementResult r=NcrnaBoundaryScorer.refineBoundariesWithCounts(
				(mask&1)!=0?null:active,(mask&2)!=0?null:active,model,100,300,model,null,null,null,
				0,0,0,0,.5f,200f,offsets,offsets,0,0,false,false,NcrnaBoundaryScorer.FEATURES_V3,.8f,
				null,null,null,null,mask);
			check(((mask&1)==0 || r.startOffset==0 && r.startSitesScored==0)
				&& ((mask&2)==0 || r.stopOffset==0 && r.stopSitesScored==0),
				"A missing physical endpoint must never score a net or influence the other endpoint's sequential sweep");
			check(mask==3?active.calls==0:active.calls==3,
				"Only the interior end should score the complete three-position sweep");
			if(mask==1){check(r.stopOffset==1,"Interior stop retains ordinary length-sensitive NN selection");}
			if(mask==2){check(r.startOffset==-1,"Interior start retains ordinary length-sensitive NN selection");}
		}
	}
	static final class LengthNet extends ml.CellNet{
		float[] values;int calls;
		@Override public ml.CellNet applyInput(float[] input){values=input;return this;}
		@Override public float feedForward(){calls++;return values[4];}
	}
	static void edge18Defaults(byte[] model){
		final LongHashSet seeds=new LongHashSet();
		for(int i=0; i+17<=model.length; i++){
			long key=0;for(int j=0; j<17; j++){key=(key<<2)|AminoAcid.baseToNumber[model[i+j]];}seeds.add(key);
		}
		final byte[] fragment=Arrays.copyOfRange(model, 200, 1000);
		for(int i=40; i<fragment.length-10; i+=50){fragment[i]=(byte)(fragment[i]=='A'?'C':'A');}
		String expected=null;
		for(boolean legacyTarget:new boolean[]{true, false}){
			final NcrnaFamily f=new NcrnaFamily("euk18S", new byte[][]{model}, null,
				new String[]{"euk18S_mitochondria"}, seeds, 17, 600, 1500, 9, 1, false, 1f, 0f, 0f, 1, 0f, 1f, 1f, 1f);
			f.outputType=ProkObject.r18S;f.reuseConsensusAlignment=true;f.trimAlignmentExtent=false;f.scavengePass2=false;
			f.modelClipRescue=legacyTarget;f.modelClipRescueId=legacyTarget ? .66f : Float.NaN;
			NcrnaFamilyRegistry.enablePhysicalEdges(f);
			check(f.modelClipRescue==legacyTarget && f.modelClipRescueId==.66f,
				"Physical-edge registration must retain the measured18S cutoff without enabling interior rescue");
			final ArrayList<Orf> result=Prok16sPartialRescueTest.scavenger(f).scavenge("edge18", fragment, 0, new ArrayList<int[]>());
			check(result.size()==1 && result.get(0).ssuClippedMitoRescue && result.get(0).ssuFullModel==model,
				"Default edge rescue must preserve mitochondrial reconciliation provenance as well as acceptance below the primary cutoff");
			final String calls=Prok16sPartialRescueTest.calls(result);
			if(expected==null){expected=calls;}else{check(expected.equals(calls), "Disabling interior rescue must not alter the saved TARGET physical-edge GFF");}
		}
	}
	static void modelCoordinates(byte[] model){
		final EndClippedAligner aligner=new EndClippedAligner();final BaseGraph graph=graph(model);
		for(int mode=1;mode<=3;mode++){
			final int from=mode==2?0:200, to=mode==1?1000:850;final byte[] fragment=Arrays.copyOfRange(model,from,to);
			final EndClippedAligner.Result stats=aligner.alignContigEdges(model,fragment,mode!=2,mode!=1);
			check(stats!=null && stats.qStart==from && stats.qStop==to-1,"Expected true query/model clipping");
			final Read r=TrnaConsensusBuilder.invertToModelFrame(fragment,stats,stats.qStart);
			check(r.start==from && r.stop==to-1 && Arrays.equals(r.bases,fragment),"HBM conversion must retain original model coordinates on both sides");
			check(TrnaConsensusBuilder.scoreAlignedAgainstModel(fragment,stats,graph,stats.qStart)>.99f,"Correctly positioned exact fragment must score against its supported HBM span");
		}
		final byte[] window=Prok16sPartialRescueTest.embed(model,30,40);final AlignmentStats full=new AlignmentStats(true);
		ScrabbleAligner.alignAndTraceStatic(model,window,full);
		final Read a=TrnaConsensusBuilder.invertToModelFrame(window,full), b=TrnaConsensusBuilder.invertToModelFrame(window,full,0);
		check(a.start==0 && a.stop==model.length-1 && a.start==b.start && a.stop==b.stop && Arrays.equals(a.bases,b.bases) && Arrays.equals(a.match,b.match),"Legacy zero-origin conversion must retain exact bases/match/coordinates");
		check(Float.floatToIntBits(TrnaConsensusBuilder.scoreAlignedAgainstModel(window,full,graph))==Float.floatToIntBits(TrnaConsensusBuilder.scoreAlignedAgainstModel(window,full,graph,0)),"Zero-origin HBM score must remain bit-identical");
		final AlignmentStats deletions=new AlignmentStats(true);deletions.qLen=20;deletions.rLen=5;deletions.rStart=0;deletions.rStop=4;deletions.setFromMatchString("IImmmmm".getBytes(java.nio.charset.StandardCharsets.US_ASCII));
		final Read offset=TrnaConsensusBuilder.invertToModelFrame("ACGTA".getBytes(java.nio.charset.StandardCharsets.US_ASCII),deletions,7);
		check(offset.start==9 && offset.stop==13,"Leading inverted model deletions advance from qStart rather than zero");
	}
	static void policy(byte[] model){
		final NcrnaScavenger unnamed=new NcrnaScavenger(new byte[][]{model},null,null,null,17,600,1500);unnamed.physicalEdgeMinLen=600;
		boolean failed=false;try{unnamed.scavenge("unnamed",model,0,new ArrayList<int[]>());}catch(IllegalStateException e){failed=e.getMessage().contains("requires named models");}
		check(failed,"Unannotated models must not let a clipped end enter the ordinary voting fallback");
		for(String name:new String[]{"16S","23S","prok5S","euk5S","euk18S","r58","lsu","rnasep","tmrna","srp_small","srp_large","sixs_rf00013","sixs_rf01685"}){
			final NcrnaFamily f=new NcrnaFamily(name,new byte[][]{model},null,new String[]{"model"},null,17,600,1500,9,1,false,1,0,0,1);
			if(name.equals("16S") || name.equals("23S")){f.prok16sPartialMinLen=500;}
			NcrnaFamilyRegistry.enablePhysicalEdges(f);final NcrnaScavenger s=Prok16sPartialRescueTest.scavenger(f);
			check(s.physicalEdgeMinLen==(f.prok16sPartialMinLen>0?500:600),"Registered policy must reach actual worker without tuning length");
		}
	}
	static NcrnaScavenger caller(byte[] model,boolean reused,boolean edge,BaseGraph graph,float pass,float borderline){
		final LongHashSet seeds=new LongHashSet();for(int i=0;i+17<=model.length;i++){long key=0;for(int j=0;j<17;j++){key=(key<<2)|AminoAcid.baseToNumber[model[i+j]];}seeds.add(key);}
		final NcrnaScavenger s=new NcrnaScavenger(new byte[][]{model},graph==null?null:new BaseGraph[]{graph},new String[]{"bacteria"},seeds,17,600,1500,
			9,1,false,1,0,0,1,0,1,pass,borderline);
		s.family=reused?"16S":"r58";s.outputType=reused?ProkObject.r16S:ProkObject.RNA;s.reuseConsensusAlignment=reused;s.trimAlignmentExtent=!reused;
		s.scavengePass2=false;s.physicalEdgeMinLen=edge?600:0;s.hbmPass=.9f;return s;
	}
	static void its(){
		final byte[] bases=Prok16sPartialRescueTest.random(400);
		for(int strand=0;strand<2;strand++){
			final Orf a=new Orf("its",0,99,strand,0,bases,false,ProkObject.r18S), b=new Orf("its",150,249,strand,0,bases,false,ProkObject.RNA), c=new Orf("its",300,399,strand,0,bases,false,ProkObject.RNA);
			a.ncrnaFamily="euk18S";b.ncrnaFamily="r58";c.ncrnaFamily="lsu";a.ncrnaPartial=1;c.ncrnaPartial=2;
			if(strand==1){a.flip();b.flip();c.flip();}
			final ArrayList<Orf> calls=new ArrayList<Orf>();calls.add(a);calls.add(b);calls.add(c);
			final ArrayList<RrnaItsDeriver.Derived> results=CallGenes.deriveItsCore(new Read(bases,null,"its",0),calls);
			check(results.size()==3 && a.ncrnaPartial==1 && c.ncrnaPartial==2,"ITS derivation must preserve called anchors and their partial provenance");
			for(RrnaItsDeriver.Derived d:results){check(d.start>=100 && d.stop<=299,"Physical outer clipping must not extend an interior spacer beyond its anchors");}
		}
	}
	static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
	static BaseGraph graph(byte[] model){
		assert(model!=null && model.length>0):"Synthetic HBM must mirror production's initialized probability tables";
		final BaseGraph graph=new BaseGraph("model",model,null,0,0);consensus.BaseGraphHelper.initForScoring(graph);return graph;
	}
}
