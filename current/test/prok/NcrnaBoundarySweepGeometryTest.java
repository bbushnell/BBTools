package test.prok;

import java.util.Arrays;
import ml.CellNet;
import prok.NcrnaBoundaryScorer;
import prok.TrnaBoundaryFeatures;

/** Wide sweeps may cross the opposite end without losing real genomic context.
 * Underpadding an interior candidate must still fail loudly.
 * @author Brian Bushnell, Raiden */
public final class NcrnaBoundarySweepGeometryTest {
	public static void main(String[] args){
		boolean enabled=false;assert(enabled=true);check(enabled,"Run this invariant regression with assertions enabled");
		final NcrnaBoundaryScorer.RefinementResult longSpan=score(new ConstantNet(),600,1000,2000,64,false,false);
		check(longSpan.startSitesScored==129 && longSpan.stopSitesScored==129 && longSpan.startOffset==0 && longSpan.stopOffset==0,
			"Long interior spans must retain the unchanged complete sweep and tie policy");
		final NcrnaBoundaryScorer.RefinementResult shortSpan=score(new ConstantNet(),600,630,2000,64,false,false);
		check(shortSpan.startSitesScored==80 && shortSpan.stopSitesScored==80,"A31base span admits80legal positions per end at radius64");
		final NcrnaBoundaryScorer.RefinementResult minimum=score(new ConstantNet(),600,615,2000,64,false,false);
		check(minimum.startSitesScored==65 && minimum.stopSitesScored==65,"A16base span keeps offset0 and only outward moves");
		for(boolean stopFirst:new boolean[]{false,true}){
			final NcrnaBoundaryScorer.RefinementResult moved=score(new MovingNet(stopFirst),600,700,2000,64,false,false);
			check(stopFirst ? moved.stopOffset==-64 && moved.startOffset==0 && moved.stopSitesScored==129 && moved.startSitesScored==86
				: moved.startOffset==64 && moved.stopOffset==0 && moved.startSitesScored==129 && moved.stopSitesScored==86,
				"The second sweep must count legal positions against the first end's inward move, in either refinement order");
		}
		final NcrnaBoundaryScorer.RefinementResult edge=score(new ConstantNet(),5,500,1000,32,true,false);
		check(edge.startSitesScored==38 && edge.stopSitesScored==65,"Only a real physical left edge may clip the outward start sweep");
		boolean failed=false;try{score(new ConstantNet(),5,500,1000,32,false,false);}catch(AssertionError e){
			failed=e.getMessage()!=null && e.getMessage().contains("missing context for a valid candidate");}
		check(failed,"An underpadded interior window must still fail, even though candidate spans themselves are valid");
		failed=false;try{score(new ConstantNet(),500,994,1000,32,false,false);}catch(AssertionError e){
			failed=e.getMessage()!=null && e.getMessage().contains("missing context for a valid candidate");}
		check(failed,"Missing right context must remain detectable too");
		final NcrnaBoundaryScorer.RefinementResult right=score(new ConstantNet(),500,994,1000,32,false,true);
		check(right.startSitesScored==65 && right.stopSitesScored==38,"A real physical right edge may clip only the outward stop sweep");
		System.out.println("PASS NcrnaBoundarySweepGeometryTest short_span_sequential_both_orders_physical_edges_underpadding");
	}
	private static NcrnaBoundaryScorer.RefinementResult score(CellNet net,int start,int stop,int length,int radius,boolean left,boolean right){
		check(start>=0 && stop<length && stop-start>=15,"Fixture spans must begin valid");
		final byte[] bases=new byte[length];Arrays.fill(bases,(byte)'A');final int[] offsets=new int[2*radius+1];
		for(int i=0;i<offsets.length;i++){offsets[i]=i-radius;}
		final TrnaBoundaryFeatures.NinemerTable table=(seq,pos)->0;
		return NcrnaBoundaryScorer.refineBoundariesWithCounts(net,net,bases,start,stop,new byte[100],null,
			table,table,9,0,9,0,.5f,1000f,offsets,offsets,0,0,left,right,NcrnaBoundaryScorer.FEATURES_V3,.8f,null,null);
	}
	private static class ConstantNet extends CellNet {
		@Override public CellNet applyInput(float[] values){NcrnaBoundarySweepGeometryTest.check(values.length==10,"Fixture uses the real V3 layout");input=values;return this;}
		@Override public float feedForward(){return 0;}
		float[] input;
	}
	private static final class MovingNet extends ConstantNet {
		MovingNet(boolean stopFirst_){stopFirst=stopFirst_;}
		@Override public float feedForward(){return (input[3]==1)==stopFirst ? 1-input[4] : 100;}
		final boolean stopFirst;
	}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
