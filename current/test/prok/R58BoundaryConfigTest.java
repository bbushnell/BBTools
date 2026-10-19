package test.prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import ml.CellNet;
import prok.CallGenes;
import prok.GeneCaller;
import prok.NcrnaFamily;
import prok.NcrnaScavenger;
import prok.Orf;
import prok.ProkObject;
import prok.R58BoundaryConfig;
import prok.TrnaBoundaryFeatures;

/** Range routing and actual NN invocation after a rejected degrading trim.
 * @author Raiden */
public final class R58BoundaryConfigTest {
	public static void main(String[] args){
		final R58BoundaryConfig c=CallGenes.R58_BOUNDARY_CONFIG;c.clear();
		try{
			final int[] baseline={-3,-2,-1,0,1,2};
			final int[] start=CallGenes.boundaryStartOffsets("r58"),stop=CallGenes.boundaryStopOffsets("r58");
			check(start.length==68 && start[0]==-11 && start[67]==56 && stop.length==84 && stop[0]==-46 && stop[83]==37, "Selected correction coverage and margin must reach production factories");
			start[0]=99;check(CallGenes.boundaryStartOffsets("r58")[0]==-11, "Factories must not expose mutable shared configuration");
			check(c.parse("r58boundarystartrange", "-5,8") && c.parse("r58boundarystoprange", "-7,9"), "Both CLI ranges must be recognized");
			check(CallGenes.boundaryStartOffsets("r58").length==14 && CallGenes.boundaryStopOffsets("r58").length==17, "Explicit ranges override the selected default geometry");
			check(Arrays.equals(CallGenes.boundaryStartOffsets("lsu"), baseline), "5.8S overrides cannot change LSU");
			fails(()->c.validate(false));c.validate(true);
			for(String bad:new String[]{"1,4","-5,-1","-300,2","-200,200","0","0,1,2","2,-2"}){fails(()->R58BoundaryConfig.range(bad));}
			check(R58BoundaryConfig.range("0,0").length==1, "A fixed endpoint remains a valid control");
			c.clear();invocation(false);invocation(true);
		}finally{c.clear();}
		System.out.println("PASS R58 selected ranges, overrides, family isolation and real retained-anchor NN invocation");
	}
	private static void invocation(boolean raw){
		final byte[] model=bytes("ACGTTGCAAGTCGATCGTACGATGC"),bases=bytes("TTTTTTTTTTACGTTGCAAGTCGATCGTACGATGCGGGGGGGGGG");bases[20]='A';
		final CountingNet net=new CountingNet();final int[] offsets={-1,0,1};
		final TrnaBoundaryFeatures.NinemerTable table=(seq,pos)->1;
		final NcrnaFamily f=new NcrnaFamily("r58",new byte[][]{model},null,null,null,17,10,10,7,1,false,0f,0f,0f,1,0f,1f,.5f,.5f,.5f,.85f,
			net,net,table,table,7,0,7,0,model.length,offsets,offsets);
		f.trimAlignmentExtent=true;
		final R58BoundaryConfig c=new R58BoundaryConfig();c.parse("r58boundaryraw", Boolean.toString(raw));c.apply(f,true);
		final NcrnaScavenger s=new NcrnaScavenger(f.library,null,null,null,17,10,10,7,1,false,0f,0f,0f,1,0f,1f,.5f,.5f,
			net,net,table,table,7,0,7,0,model.length,offsets,offsets);
		GeneCaller.applyNcrnaFamilyControls(s,f);
		final Orf o=new Orf("fixture",8,36,0,0,bases,false,ProkObject.RNA);
		final boolean trimmed=s.finishAcceptedBoundary(o,bases,0,8,36,0,1f,29,1f,true,null,null);
		check(!trimmed && o.start==8 && o.stop==36, "Rejected trim and tied NN scores must retain the verified anchor");
		check(raw?net.calls>0:net.calls==0, "Configured raw-anchor eligibility must reach actual NN scoring");
	}
	private static final class CountingNet extends CellNet{
		@Override public CellNet copy(boolean copyWeights){return this;}
		@Override public CellNet applyInput(float[] input){return this;}
		@Override public float feedForward(){calls++;return 0f;}
		int calls;
	}
	private static byte[] bytes(String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static void fails(Runnable action){try{action.run();}catch(IllegalArgumentException expected){return;}throw new AssertionError("Invalid boundary configuration was accepted");}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
