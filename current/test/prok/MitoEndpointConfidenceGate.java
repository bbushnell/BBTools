package test.prok;

import prok.NcrnaBoundaryScorer;

import ml.CellNet;

/** Exact cached absolute-score gating for the existing two-step V3 scorer.
 * The first-end decision has two possible states: moved or raw. Therefore only
 * two second-end sweeps are needed for any number of score thresholds.
 * Scores are regression outputs, not calibrated probabilities.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden */
public final class MitoEndpointConfidenceGate {
	public MitoEndpointConfidenceGate(CellNet net_,byte[] bases_,int start,int stop,int model_,byte[] consensus_,
			float gc_,float identity_,MitoEndpointTables tables_,int[] offsets_,int edges_){
		require(net_!=null && tables_!=null && bases_!=null && start>=0 && stop<bases_.length && stop-start>=15
			&& model_>=0 && model_<tables_.models && (edges_&~3)==0,"Gate cache must retain a valid complete raw anchor, bound model and physical context");
		net=net_;bases=bases_;model=model_;consensus=consensus_;gc=gc_;identity=identity_;tables=tables_;offsets=offsets_;edges=edges_;
		startConfidence=score(start,stop,false);stopConfidence=score(start,stop,true);
		require(Float.isFinite(startConfidence) && Float.isFinite(stopConfidence),"Initial endpoint scores must be finite before ordering refinement");
		startFirst=startConfidence<=stopConfidence;
		first=sweep(start,stop,!startFirst);
		final int move=offset(first,Float.NEGATIVE_INFINITY);
		secondMoved=sweep(start+(startFirst?move:0),stop+(startFirst?0:move),startFirst);
		secondRaw=move==0?secondMoved:sweep(start,stop,startFirst);
		physicalScores=2+(int)first[3]+(int)secondMoved[3]+(secondRaw==secondMoved?0:(int)secondRaw[3]);
	}
	public NcrnaBoundaryScorer.RefinementResult at(float threshold){
		require(!Float.isNaN(threshold),"NaN cannot define a confidence gate; infinities are explicit ungated/raw controls");
		final int firstOffset=offset(first,threshold);final float[] second=firstOffset==0?secondRaw:secondMoved;
		final int secondOffset=offset(second,threshold);
		return startFirst?new NcrnaBoundaryScorer.RefinementResult(firstOffset,secondOffset,(int)first[3],(int)second[3],chosen(first,threshold),chosen(second,threshold))
			:new NcrnaBoundaryScorer.RefinementResult(secondOffset,firstOffset,(int)second[3],(int)first[3],chosen(second,threshold),chosen(first,threshold));
	}
	private float score(int start,int stop,boolean end){
		final int e=end?1:0,inside=tables.insideCount(e);
		return NcrnaBoundaryScorer.score(net,bases,start,stop,end,consensus,null,tables.table(model,e),inside,9-inside,gc,tables.meanLength,0,
			NcrnaBoundaryScorer.FEATURES_V3,identity,null,null,tables.table(-1,e));
	}
	private float[] sweep(int start,int stop,boolean end){
		final int e=end?1:0,inside=tables.insideCount(e);
		final float[] result=NcrnaBoundaryScorer.bestOffset(net,bases,start,stop,end,consensus,null,tables.table(model,e),inside,9-inside,gc,tables.meanLength,
			offsets,(edges&(end?2:1))!=0,0,NcrnaBoundaryScorer.FEATURES_V3,identity,null,null,tables.table(-1,e));
		require(result.length==4 && Float.isFinite(result[1]) && Float.isFinite(result[2]),"Cached sweep must retain finite best and raw scores");return result;
	}
	private static int offset(float[] sweep,float threshold){return sweep[1]>=threshold && sweep[1]-sweep[2]>0?(int)sweep[0]:0;}
	private static float chosen(float[] sweep,float threshold){return offset(sweep,threshold)==0?sweep[2]:sweep[1];}
	public float firstBest(){return first[1];}public float firstRaw(){return first[2];}
	public float secondMovedBest(){return secondMoved[1];}public float secondMovedRaw(){return secondMoved[2];}
	public float secondUnmovedBest(){return secondRaw[1];}public float secondUnmovedRaw(){return secondRaw[2];}
	public int firstBestOffset(){return (int)first[0];}public int secondMovedBestOffset(){return (int)secondMoved[0];}public int secondUnmovedBestOffset(){return (int)secondRaw[0];}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	public final boolean startFirst;public final float startConfidence,stopConfidence;public final int physicalScores;
	private final float[] first,secondMoved,secondRaw;private final CellNet net;private final byte[] bases,consensus;private final int model,edges;
	private final float gc,identity;private final MitoEndpointTables tables;private final int[] offsets;
}
