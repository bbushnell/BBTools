package test.prok;


import java.util.Arrays;
import ml.CellNet;
import prok.NcrnaBoundaryScorer;

/** Count pooling and feature parity against the real V3 inference path.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden
 */
public final class MitoEndpointTablesTest {
	public static void main(String[] args){
		final byte[] seq=new byte[200];Arrays.fill(seq,(byte)'A');
		for(int i=40;i<49;i++){seq[i]=(byte)'C';}for(int i=132;i<=140;i++){seq[i]=(byte)'G';}
		final float mean=MitoEndpointTables.medianLength(new byte[][]{new byte[101],new byte[173]},0);
		check(mean==137,"Unequal model lengths must still use one family median");
		final MitoEndpointTables t=new MitoEndpointTables(2,9,9,mean);
		t.add(seq,40,140,0);t.add(seq,40,140,1);t.add(seq,40,140,1);
		check(t.table(0,0).rawCount(seq,40)==1 && t.table(1,0).rawCount(seq,40)==2 && t.table(-1,0).rawCount(seq,40)==3,"Model/family observations must preserve actual member multiplicity");
		check(t.table(0,1).rawCount(seq,132)==1 && t.table(-1,1).rawCount(seq,132)==3,"Stop k9 must cover nine inside bases ending at the true stop");
		check(t.table(0,0).rawCount(seq,-1)==0 && t.table(0,0).rawCount(seq,199)==0,"Clipped lookup windows return no observation");
		final byte[] ambiguous=seq.clone();ambiguous[44]=(byte)'N';check(t.table(0,0).rawCount(ambiguous,40)==0,"Ambiguous bases cannot be counted as an exact kmer");
		final CaptureNet start=new CaptureNet(),stop=new CaptureNet();
		NcrnaBoundaryScorer.refineBoundariesWithCounts(start,stop,seq,41,139,new byte[101],null,t.table(0,0),t.table(0,1),9,0,9,0,
			.42f,mean,new int[]{-1,0,1},new int[]{-1,0,1},0,0,false,false,NcrnaBoundaryScorer.FEATURES_V3,.77f,null,null,t.table(-1,0),t.table(-1,1));
		check(Arrays.equals(start.first,t.features(seq,41,139,0,0,.42f,.77f,-1)),"Packed-table start vectors must equal actual V3 scorer inputs bit for bit");
		check(Arrays.equals(stop.first,t.features(seq,41,139,0,1,.42f,.77f,-1)),"Packed-table stop vectors must equal actual V3 scorer inputs bit for bit");
		check(t.features(seq,41,139,0,0,.42f,.77f,-1)[4]==99/mean
			&& t.features(seq,41,139,1,0,.42f,.77f,-1)[4]==99/mean,"Model choice must never replace the shared length denominator");
		final float[] training=t.features(seq,40,140,0,0,.42f,.77f,40),inference=t.features(seq,40,140,0,0,.42f,.77f,-1);
		check(training[0]<inference[0],"Declared training leave-one-out must remove the member from both profile observations");
		check(!NcrnaBoundaryScorer.usesPerSiteScrabble(NcrnaBoundaryScorer.FEATURES_V3),"Packed profiles must retain the alignment-free candidate-site path");
		final SequentialNet sequential=new SequentialNet();
		final NcrnaBoundaryScorer.RefinementResult refined=t.refine(sequential,seq,41,139,0,new byte[101],.42f,.77f,new int[]{-1,0,1},new int[]{-1,0,1},0);
		check(refined.startOffset==-1 && sequential.changedStopSpan,"Second endpoint must receive the first endpoint's changed span through the actual sequential scorer");
		final int edgeMask=MitoEndpointTables.windowEdgeMask(new MitoEndpointGeometry("41..141:+","+",false),500,40,140,200);
		check(edgeMask==1,"A window may touch a physical contig end while the raw anchor is inside it");
		final int[] edgeOffsets=new int[101];for(int i=0;i<edgeOffsets.length;i++){edgeOffsets[i]=i-50;}
		final NcrnaBoundaryScorer.RefinementResult edge=t.refine(new CaptureNet(),seq,41,139,0,new byte[101],.42f,.77f,edgeOffsets,new int[]{-1,0,1},edgeMask);
		check(edge.startSitesScored==92,"Physical-edge provenance must allow exactly the92in-bounds start candidates");
		check(MitoEndpointTables.windowEdgeMask(new MitoEndpointGeometry("360..460:-","-",false),500,40,140,200)==1,"Minus-strand physical edge must be converted into sense coordinates");
		check(args.length==1,"Provide one fresh table round-trip path in the run directory");
		final String pin="0123456789abcdef0123";t.write(args[0],pin);final MitoEndpointTables copy=MitoEndpointTables.load(args[0],2,pin);
		check(copy.meanLength==mean && Arrays.equals(copy.features(seq,41,139,0,0,.42f,.77f,-1),t.features(seq,41,139,0,0,.42f,.77f,-1)),"Serialization must preserve normalization and every table-driven feature");
		System.out.println("PASS MitoEndpointTablesTest counts_pooling_k9_V3_parity_leave_one_out_family_mean_sequential_roundtrip");
	}
	private static final class SequentialNet extends CellNet{
		@Override public CellNet applyInput(float[] values){input=values;return this;}
		@Override public float feedForward(){
			if(input[3]==1){if(input[4]==100/137f){changedStopSpan=true;}return 2;}
			return input[1]-input[2];
		}
		float[] input;boolean changedStopSpan;
	}
	private static final class CaptureNet extends CellNet{
		@Override public CellNet applyInput(float[] values){if(first==null){first=values.clone();}return this;}
		@Override public float feedForward(){return 0;}
		float[] first;
	}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
