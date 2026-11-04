package test.prok;


import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import ml.CellNet;
import prok.NcrnaBoundaryScorer;
import prok.TrnaBoundaryFeatures;

/** Source geometry and correction signs are checked against actual scorer sweeps.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden
 */
public final class MitoEndpointGeometryTest {
	public static void main(String[] args){
		final byte[] genome="AACCGGTTAACCGGTT".getBytes(StandardCharsets.US_ASCII);
		final MitoEndpointGeometry plus=new MitoEndpointGeometry("3..6:+,9..10:+", "+", false);
		plus.verify(genome, bytes("CCGGAA"), "joined-plus");
		final MitoEndpointGeometry.Window p=plus.window(genome,2,100);
		check(p!=null && p.start==2 && p.stop==9 && Arrays.equals(p.bases,bytes("AACCGGTTAACC")), "Joined genomic flank must retain its intron");
		final MitoEndpointGeometry minus=new MitoEndpointGeometry("9..10:-,3..6:-", "-", false);
		minus.verify(genome,bytes("TTCCGG"),"joined-minus");
		final MitoEndpointGeometry.Window m=minus.window(genome,2,100);
		check(m!=null && m.start==2 && m.stop==9 && Arrays.equals(m.bases,bytes("GGTTAACCGGTT")), "Minus window must be transcript-oriented with genuine genomic flanks");
		final MitoEndpointGeometry origin=new MitoEndpointGeometry("14..16:+,1..3:+", "+", true);
		origin.verify(genome,bytes("GTTAAC"),"origin-plus");
		final MitoEndpointGeometry.Window o=origin.window(genome,2,100);
		check(o!=null && o.start==2 && o.stop==7 && Arrays.equals(o.bases,bytes("CGGTTAACCG")), "Origin crossing must unwrap exactly one genome turn");
		final MitoEndpointGeometry reverseOrigin=new MitoEndpointGeometry("1..3:-,14..16:-", "-", true);
		reverseOrigin.verify(genome,bytes("GTTAAC"),"origin-minus");
		check(reverseOrigin.window(genome,2,100)!=null, "Reverse origin crossing is an ordered sense span");
		check(origin.window(genome,100,100).bases.length==genome.length, "Circular flanks must not repeat the same genome");
		check(plus.window(genome,2,7)==null, "Budget cap must flag a long genomic span, never splice it to fit");
		check(new MitoEndpointGeometry("9..10:+,3..6:+", "+", false).window(genome,2,100)==null, "Undeclared origin wrapping cannot be invented");
		check(new MitoEndpointGeometry("3..6:+,9..10:-", "MIXED", false).window(genome,2,100)==null, "Mixed strands cannot become one linear caller window");
		final MitoEndpointGeometry.Window edge=new MitoEndpointGeometry("1..4:+", "+", false).window(genome,10,100);
		check(edge.start==0 && edge.stop==3, "Physical edge must retain the actual clipped flank");
		boolean failed=false;try{plus.verify(genome,bytes("CCGGAT"),"mismatch");}catch(IllegalArgumentException e){failed=true;}
		check(failed,"Wrong original genomic bases must fail before training");
		final MitoEndpointGeometry.Window extendedPlus=plus.extendWindow(genome,p,20);
		check(Arrays.equals(extendedPlus.bases,genome) && extendedPlus.start==p.start,"Physical left clipping must not invent padding or move the truth");
		final MitoEndpointGeometry.Window extendedMinus=minus.extendWindow(genome,m,3);
		check(Arrays.equals(extendedMinus.bases,bytes("ACCGGTTAACCGGTT")) && extendedMinus.start==m.start+3,"Minus context must extend the sense sequence and translate all coordinates equally");
		final MitoEndpointGeometry.Window circular=origin.extendWindow(genome,o,20);
		check(circular.bases.length==o.bases.length+40 && circular.start==o.start+20
			&& Arrays.equals(Arrays.copyOfRange(circular.bases,20,20+o.bases.length),o.bases),"Explicit circular source may supply repeated extra context while preserving the one original window");
		check(circular.bases[0]==genome[7] && circular.bases[16]==genome[7],"Circular flank bases must follow the actual source period, not zero padding");
		final MitoEndpointGeometry.Window reverseOld=reverseOrigin.window(genome,2,100),reverseWide=reverseOrigin.extendWindow(genome,reverseOld,20);
		check(reverseWide.start==reverseOld.start+20 && reverseWide.stop==reverseOld.stop+20
			&& Arrays.equals(Arrays.copyOfRange(reverseWide.bases,20,20+reverseOld.bases.length),reverseOld.bases),"Reverse-origin expansion must preserve the original sense bases and translate both true ends identically");
		extendedSweep(false);extendedSweep(true);sign(false);sign(true);
		System.out.println("PASS MitoEndpointGeometryTest plus_minus_join_origin_edge_base_parity_actual_scorer_sign");
	}
	private static void extendedSweep(boolean reverse){
		final byte[] genome=new byte[2000];Arrays.fill(genome,(byte)'A');
		final MitoEndpointGeometry geometry=new MitoEndpointGeometry("701..900:"+(reverse?"-":"+"),reverse?"-":"+",false);
		final MitoEndpointGeometry.Window old=geometry.window(genome,100,1000),wide=geometry.extendWindow(genome,old,203);
		final int shift=wide.start-old.start,s=0,e=old.bases.length-1;final int[] offsets=new int[385];for(int i=0;i<385;i++){offsets[i]=i-192;}
		final TrnaBoundaryFeatures.NinemerTable zero=(seq,pos)->0;
		boolean failed=false;try{NcrnaBoundaryScorer.refineBoundariesWithCounts(new PeakNet(),new PeakNet(),old.bases,s,e,new byte[200],null,
			zero,zero,9,0,9,0,.5f,1000f,offsets,offsets,0,0,false,false,NcrnaBoundaryScorer.FEATURES_V3,.8f,null,null);
		}catch(AssertionError error){failed=error.getMessage().contains("missing context");}
		check(failed,"A raw stop at an artificial crop edge must fail before real genomic context is supplied");
		final TrnaBoundaryFeatures.NinemerTable strict=(seq,pos)->{check(pos>=0 && pos+9<=seq.length,"Every profile's +/-2 neighbor must retain its whole9mer, including extreme candidates");return 0;};
		for(int inside:new int[]{0,9}){
			final NcrnaBoundaryScorer.RefinementResult result=NcrnaBoundaryScorer.refineBoundariesWithCounts(new PeakNet(),new PeakNet(),wide.bases,s+shift,e+shift,new byte[200],null,
				strict,strict,inside,9-inside,inside,9-inside,.5f,1000f,offsets,offsets,0,0,false,false,NcrnaBoundaryScorer.FEATURES_V3,.8f,null,null);
			check(result.startSitesScored==385 && result.stopSitesScored==385 && result.startOffset==0 && result.stopOffset==0,
				"Real expanded context must restore every candidate and profile without changing fixed anchors or claiming physical edges, on either strand");
		}
	}
	private static void sign(boolean reverse){
		final int len=256,truthLo=50,truthHi=149,callLo=53,callHi=145;
		final int ts=reverse?len-1-truthHi:truthLo,te=reverse?len-1-truthLo:truthHi;
		final int s=reverse?len-1-callHi:callLo,e=reverse?len-1-callLo:callHi;
		final byte[] genome=new byte[len];Arrays.fill(genome,(byte)'A');
		final TrnaBoundaryFeatures.NinemerTable start=(seq,pos)->pos==ts?100:0;
		final TrnaBoundaryFeatures.NinemerTable stop=(seq,pos)->pos==te-8?100:0;
		final int[] offsets={-4,-3,0,3,4};
		final NcrnaBoundaryScorer.RefinementResult result=NcrnaBoundaryScorer.refineBoundariesWithCounts(new PeakNet(),new PeakNet(),genome,s,e,
			new byte[100],null,start,stop,9,0,9,0,.5f,100f,offsets,offsets,0,0,false,false,NcrnaBoundaryScorer.FEATURES_V3,.8f,null,null);
		check(result.startOffset==(reverse?-4:-3) && result.stopOffset==(reverse?3:4), "Actual caller sweep must add the corrections on "+(reverse?"minus":"plus"));
		check(result.startOffset==MitoEndpointGeometry.correction(ts,s) && result.stopOffset==MitoEndpointGeometry.correction(te,e), "Measured reachability sign must match actual caller scoring");
		check(s+result.startOffset==ts && e+result.stopOffset==te,"Refined coordinates must equal the independently specified truth");
		check(result.startSitesScored==offsets.length && result.stopSitesScored==offsets.length,"Both actual sweeps must score every candidate");
	}
	private static final class PeakNet extends CellNet{
		@Override public CellNet applyInput(float[] in){MitoEndpointGeometryTest.check(in.length==10 && in[6]==.8f && in[7]==0 && in[8]==0 && in[9]==0,"Actual V3 scorer feature contract changed");score=in[1]-in[2];return this;}
		@Override public float feedForward(){return score;}
		private float score;
	}
	private static byte[] bytes(String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
