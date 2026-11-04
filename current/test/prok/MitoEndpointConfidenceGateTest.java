package test.prok;



import java.util.Arrays;
import ml.CellNet;
import prok.NcrnaBoundaryScorer;

/** Gating the first end changes the second end's length feature.
 * Preserved validation of mitochondrial endpoint refinement. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. This fixture remains
 * useful if that investigation resumes; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden
 */
public final class MitoEndpointConfidenceGateTest {
	public static void main(String[] args){
		for(boolean stopFirst:new boolean[]{false,true}){for(boolean negative:new boolean[]{false,true}){test(stopFirst,negative);}}
		System.out.println("PASS MitoEndpointConfidenceGateTest sequential_raw_or_moved_signed_scores_infinite_controls_no_rescoring");
	}
	private static void test(boolean stopFirst,boolean negative){
		final byte[] bases=new byte[2000];Arrays.fill(bases,(byte)'A');final int[] offsets=new int[129];for(int i=0;i<129;i++){offsets[i]=i-64;}
		final MitoEndpointTables tables=new MitoEndpointTables(1,9,9,100);final LengthNet net=new LengthNet(stopFirst,negative);
		final MitoEndpointConfidenceGate cache=new MitoEndpointConfidenceGate(net,bases,600,699,0,new byte[100],.5f,.8f,tables,offsets,0);
		final NcrnaBoundaryScorer.RefinementResult nativeResult=tables.refine(new LengthNet(stopFirst,negative),bases,600,699,0,new byte[100],.5f,.8f,offsets,offsets,0);
		final NcrnaBoundaryScorer.RefinementResult ungated=cache.at(Float.NEGATIVE_INFINITY);
		check(ungated.startOffset==nativeResult.startOffset && ungated.stopOffset==nativeResult.stopOffset,"Ungated cache must exactly reproduce actual production sequential refinement");
		check(cache.startFirst!=stopFirst && (stopFirst?ungated.stopOffset==-20 && ungated.startOffset==-40:ungated.startOffset==20 && ungated.stopOffset==40),"Synthetic weak-first/strong-second targets must exercise both refinement orders");
		final NcrnaBoundaryScorer.RefinementResult gated=cache.at(negative?-.3f:.7f);
		check(stopFirst?gated.stopOffset==0 && gated.startOffset==-20:gated.startOffset==0 && gated.stopOffset==20,"Rejecting the weak first-end move must select the second sweep computed at the RAW first endpoint");
		final NcrnaBoundaryScorer.RefinementResult raw=cache.at(Float.POSITIVE_INFINITY);
		check(raw.startOffset==0 && raw.stopOffset==0,"The all-reject control must preserve both raw anchors exactly");
		final NcrnaBoundaryScorer.RefinementResult equality=cache.at(cache.firstBest());
		check(stopFirst?equality.stopOffset==-20:equality.startOffset==20,"A winning score exactly equal to the cutoff must be accepted by the declared >= gate");
		final NcrnaBoundaryScorer.RefinementResult above=cache.at(Math.nextUp(cache.firstBest()));
		check(stopFirst?above.stopOffset==0:above.startOffset==0,"The next representable cutoff above the first winner must reject that move");
		for(float cut:new float[]{Float.NEGATIVE_INFINITY,cache.firstBest(),Math.nextUp(cache.firstBest()),negative?-.3f:.7f,Float.POSITIVE_INFINITY}){
			final NcrnaBoundaryScorer.RefinementResult direct=tables.refine(new ThresholdNet(stopFirst,negative,cut),bases,600,699,0,new byte[100],.5f,.8f,offsets,offsets,0),cached=cache.at(cut);
			check(direct.startOffset==cached.startOffset && direct.stopOffset==cached.stopOffset && direct.startSitesScored==cached.startSitesScored && direct.stopSitesScored==cached.stopSitesScored,
				"Cached finite gates must match actual sequential rescoring with fixed initial order, including equality and first-end rejection");
		}
		final CellNet flat=new CellNet(){@Override public CellNet applyInput(float[] in){return this;}@Override public float feedForward(){return 0;}};
		final MitoEndpointConfidenceGate tied=new MitoEndpointConfidenceGate(flat,bases,600,699,0,new byte[100],.5f,.8f,tables,offsets,0);
		check(tied.at(0).startOffset==0 && tied.at(0).stopOffset==0,"Equal candidate scores must retain raw anchors through the existing strict-improvement margin");
		for(final float invalid:new float[]{Float.NaN,Float.POSITIVE_INFINITY,Float.NEGATIVE_INFINITY}){
			final CellNet broken=new CellNet(){int n;@Override public CellNet applyInput(float[] in){return this;}@Override public float feedForward(){return ++n==3?invalid:0;}};
			boolean failed=false;try{new MitoEndpointConfidenceGate(broken,bases,600,699,0,new byte[100],.5f,.8f,tables,offsets,0);}catch(IllegalArgumentException error){failed=error.getMessage().contains("Nonfinite boundary score");}
			check(failed,"Every nonfinite candidate score, even after valid ordering scores, must fail loudly");
		}
		final int before=net.calls;for(int i=-10;i<=20;i++){cache.at(i/10f);}check(net.calls==before && before==cache.physicalScores,"Threshold count must not multiply network scoring work");
	}
	/** Native direct scorer oracle: preserve its first two ordering scores, then
	 * collapse subthreshold sites to one equal floor so its margin keeps raw. */
	private static final class ThresholdNet extends CellNet{
		ThresholdNet(boolean stop,boolean negative,float threshold_){delegate=new LengthNet(stop,negative);threshold=threshold_;}
		@Override public CellNet applyInput(float[] in){delegate.applyInput(in);return this;}
		@Override public float feedForward(){final float score=delegate.feedForward();return ++calls<=2 || score>=threshold?score:-Float.MAX_VALUE;}
		final LengthNet delegate;final float threshold;int calls;
	}
	private static final class LengthNet extends CellNet{
		LengthNet(boolean stop,boolean negative_){stopFirst=stop;negative=negative_;}
		@Override public CellNet applyInput(float[] in){MitoEndpointConfidenceGateTest.check(in.length==10,"Fixture uses the actual V3 feature layout");input=in;return this;}
		@Override public float feedForward(){calls++;final boolean weak=(input[3]==1)==stopFirst;return (weak?.6f:.9f)-(negative?1:0)-Math.abs(input[4]*100-(weak?80:120))/100;}
		float[] input;int calls;final boolean stopFirst,negative;
	}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
