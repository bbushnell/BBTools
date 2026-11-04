package test.prok;


import java.util.Arrays;
import map.LongHashSet;
import ml.CellNet;
import prok.GeneCaller;
import prok.NcrnaFamily;
import prok.NcrnaScavenger;
import prok.Orf;
import prok.ProkObject;
import prok.RrnaEndpointCallerFeatures;
import prok.RrnaPositionalKmerTable;

/** Exercises the actual accepted-candidate caller seam, including disabled and
 * original-model paths, worker copies and length fallback.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden */
public final class MitoEndpointCallerTest {
	public static void main(String[] args){
		final byte[] bases=new byte[2000];Arrays.fill(bases,(byte)'A');final byte[] original=bases.clone();
		final String[] names=new String[158];final String[] first={"universal","bacteria","archaea","plastid","mitochondria"};
		final byte[][] refs=new byte[158][];
		for(int i=0;i<158;i++){names[i]=i<5?first[i]:"mitochondria_coverage_fixture_"+i;refs[i]=new byte[100];Arrays.fill(refs[i],(byte)'A');}
		final MitoEndpointTables tables=new MitoEndpointTables(158,9,9,100f);tables.add(bases,600,699,5);
		final int[] use=new int[2];final MitoEndpointCaller.Resources resources=new MitoEndpointCaller.Resources(names,refs,new LengthNet(use),tables,192);
		final NcrnaFamily family=family(names,refs);
		family.boundaryRefinerFactory=resources;
		final int[] prokCalls={0};installProkEndpoints(family,bases,prokCalls);
		final NcrnaScavenger caller=caller(family),second=caller(family);
		check(use[1]==2,"Each GeneCaller worker must copy the net instead of sharing mutable inference state");
		for(int model=0;model<5;model++){
			final Orf o=candidate(bases,names[model],600,699,0);finish(caller,o,bases,model,true);
			check(o.start==600 && o.stop==699 && o.orfScore==42 && use[0]==0,"Every original16Smodel must bypass the mitochondrial net completely");
		}
		check(prokCalls[0]==8,"The four original refined models retain both existing position nets; the original mitochondrial model remains raw");
		for(int strand=0;strand<2;strand++){
			final Orf o=candidate(bases,names[5],600,699,strand);finish(strand==0?caller:second,o,bases,5,true);
			check(o.start==408 && o.stop==891 && o.orfScore==42,"Added mitochondrial models must refine sense coordinates in the real caller without changing detection score");
		}
		check(use[0]==1328 && prokCalls[0]==8,"Each short-anchor refinement scores277legal start sites,385stop sites after expansion,and two initial confidences without invoking prok nets");
		final Orf nonQuantum=candidate(bases,names[5],600,699,0);finish(caller,nonQuantum,bases,5,false);
		check(nonQuantum.start==600 && nonQuantum.stop==699 && use[0]==1328,"A clipped/non-Quantum identity must not be fed into the full-model V3 net");
		final Orf partial=candidate(bases,names[5],0,699,0);partial.ncrnaPartial=1;
		new MitoEndpointCaller(resources).refine(partial,bases,5,0,.8f,true,1,60,1000);
		check(partial.start==0 && partial.stop==699 && use[0]==1328,"A source-truncated call retains its raw physical endpoint and does not use the complete-source net");
		family.maxLen=200;final NcrnaScavenger bounded=caller(family);final Orf invalid=candidate(bases,names[5],600,699,0);finish(bounded,invalid,bases,5,true);
		check(invalid.start==600 && invalid.stop==699 && invalid.orfScore==42,"A jointly overlong NN prediction must retain the raw accepted call");
		family.boundaryRefinerFactory=null;
		final NcrnaScavenger off=caller(family);
		final int before=use[0];
		final Orf unchanged=candidate(bases,names[5],600,699,0);finish(off,unchanged,bases,5,true);
		check(unchanged.start==600 && unchanged.stop==699 && use[0]==before && Arrays.equals(original,bases),"NN-off leaves mitochondrial calls and genomic bases unchanged");
		family.trimAlignmentExtent=true;
		boolean rejected=false;
		try{resources.create(family);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "The factory must reject competing endpoint policies before worker inference");
		System.out.println("PASS MitoEndpointCallerTest original_models_disabled_path_both_strands_worker_copies_identity_partial_length_fallback");
	}
	private static void installProkEndpoints(NcrnaFamily f,byte[] bases,int[] calls){
		final RrnaPositionalKmerTable.Table[] five=new RrnaPositionalKmerTable.Table[158],three=new RrnaPositionalKmerTable.Table[158];
		final CellNet[] starts=new CellNet[158],stops=new CellNet[158];final boolean[] raw=new boolean[158];Arrays.fill(raw,true);
		for(int i=0;i<4;i++){
			five[i]=new RrnaPositionalKmerTable.Table(f.modelNames[i],"5prime",6,0,.5,12);
			three[i]=new RrnaPositionalKmerTable.Table(f.modelNames[i],"3prime",6,-5,.5,12);
			five[i].add(bases,600);three[i].add(bases,699);five[i].finish();three[i].finish();
			starts[i]=new PositionNet(calls);stops[i]=new PositionNet(calls);raw[i]=false;
		}
		f.setRrnaEndpointInference(new RrnaEndpointCallerFeatures.Resources("16S",f.modelNames,f.library,five,three,starts,stops,raw,true),null);
	}
	private static NcrnaFamily family(String[] names,byte[][] refs){
		final LongHashSet seeds=new LongHashSet();seeds.add(0L);
		final NcrnaFamily f=new NcrnaFamily("16S",refs,null,names,seeds,17,60,256,7,1,false,1,0,0,1,0,1,.66f,.66f);
		f.reuseConsensusAlignment=true;f.trimAlignmentExtent=false;f.maxLen=500;f.outputType=ProkObject.r16S;return f;
	}
	private static NcrnaScavenger caller(NcrnaFamily f){
		final NcrnaScavenger s=new NcrnaScavenger(f.library,f.models,f.modelNames,f.kmerSet,f.kLong,f.minLen,f.windowPad,
			f.indexK,f.indexTopN,f.adaptive,f.adaptFloor,f.adaptTopFrac,f.adaptQFrac,f.fixedMinHits,f.scoreA,f.scoreB,f.idPass,f.idBorderline);
		GeneCaller.applyNcrnaFamilyControls(s,f);return s;
	}
	private static Orf candidate(byte[] bases,String model,int start,int stop,int strand){
		final Orf o=new Orf("contig",start,stop,strand,0,bases,false,ProkObject.r16S);o.ncrnaFamily="16S";o.trnaModel=model;o.orfScore=42;return o;
	}
	private static void finish(NcrnaScavenger caller,Orf o,byte[] bases,int model,boolean quantum){
		caller.finishAcceptedBoundary(o,bases,model,0,bases.length-1,0,.8f,o.stop-o.start+1,.8f,quantum,null,null);
	}
	private static final class LengthNet extends CellNet {
		LengthNet(int[] use_){use=use_;}
		@Override public int numInputs(){return 10;}
		@Override public int numOutputs(){return 1;}
		@Override public synchronized CellNet copy(boolean copyDelta){use[1]++;return new LengthNet(use);}
		@Override public CellNet applyInput(float[] values){MitoEndpointCallerTest.check(values.length==10 && values[6]==.8f,"Actual V3 identity must reach the network");input=values;return this;}
		@Override public float feedForward(){use[0]++;return input[4];}
		final int[] use;float[] input;
	}
	private static final class PositionNet extends CellNet {
		PositionNet(int[] calls_){calls=calls_;output[12]=1;}
		@Override public int numInputs(){return 28;}
		@Override public int numOutputs(){return 25;}
		@Override public synchronized CellNet copy(boolean copyDelta){return new PositionNet(calls);}
		@Override public CellNet applyInput(float[] input){MitoEndpointCallerTest.check(input.length==28,"Existing prok endpoint geometry must remain radius12");return this;}
		@Override public float feedForward(){calls[0]++;return 0;}
		@Override public float[] getOutput(){return output;}
		final int[] calls;final float[] output=new float[25];
	}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
