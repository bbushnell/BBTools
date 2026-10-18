package prok;

import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.Random;
import idaligner.EndClippedAligner;

/** Physical-end clipping and actual caller plumbing; synthetic fixtures are not calibration.
 * @author Raiden */
public final class Euk18sContigEdgeClipTest {
	public static void main(String[] args)throws Exception{
		final byte[] model=random(1800,419);
		final EndClippedAligner aligner=new EndClippedAligner();
		final byte[] fragment=Arrays.copyOfRange(model,400,1600);
		final byte[] internal=embed(fragment,400,400);
		Arrays.fill(internal,0,400,(byte)'N');Arrays.fill(internal,1600,internal.length,(byte)'N');
		final EndClippedAligner.Result legacy=aligner.align(model,internal);
		check(legacy.qStart==400 && legacy.qStop==1599 && legacy.rStart==400 && legacy.rStop==1599,
			"Neutral flanks must expose the exact internal local trim: "+legacy.qStart+".."+legacy.qStop+" / "+legacy.rStart+".."+legacy.rStop);
		for(boolean left:new boolean[]{false,true}){for(boolean right:new boolean[]{false,true}){
			final EndClippedAligner.Result r=aligner.alignContigEdges(model,internal,left,right);
			check(r!=null && (r.qStart==0 || (left && r.rStart==0))
				&& (r.qStop==1799 || (right && r.rStop==internal.length-1)),"Even a boundary-containing window cannot restart inside its reference");
		}}
		for(int mode=1; mode<=3; mode++){
			final int from=mode==2 ? 0 : 400, to=mode==1 ? 1800 : 1600;
			final byte[] part=Arrays.copyOfRange(model,from,to);
			final byte[] ref=embed(part,mode==2 ? 300 : 0,mode==1 ? 300 : 0);
			final EndClippedAligner.Result r=aligner.alignContigEdges(model,ref,mode!=2,mode!=1);
			check(r!=null && r.qStart==from && r.qStop==to-1 && r.identity==1f,"Physical absent model ends must be clipped without trimming present sequence");
			for(int strand=0; strand<2; strand++){
				final Orf call=call(model,ref,strand,0,ref.length-1);
				check(call!=null && call.ncrnaPartial==mode,"Actual rejected-primary rescue must preserve correct sense partial flags on both strands");
				if(strand==1){call.flip();check(call.ncrnaPartial==mode,"Coordinate flips must not exchange biological model ends");}
				final int genomic=strand==0 ? mode : ((mode&1)<<1)|((mode&2)>>>1);
				final String annotation=genomic==1 ? "10" : genomic==2 ? "01" : "11";
				check(call.toGff().contains(";partial="+annotation),"GFF must expose missing physical left/right ends in genomic orientation");
			}
		}
		final byte[] full=embed(model,200,300);
		final EndClippedAligner.Result complete=aligner.alignContigEdges(model,full,false,false);
		check(complete.qStart==0 && complete.qStop==1799 && complete.rStart==200 && complete.rStop==199+model.length
			&& complete.identity==1f,"Full-model glocal alignment must keep free reference flanks");
		// A cropped candidate window with both apparent ends is still INTERNAL to the real contig.
		final byte[] containing=embed(fragment,400,400);
		check(call(model,containing,0,400,1599)==null,"Window edges cannot masquerade as physical contig ends in the actual scavenger");
		check(call(model,Arrays.copyOfRange(model,400,1315),0,0,914)==null,"Physical clipping must retain the existing916-base minimum");
		final byte[] primaryPartial=Arrays.copyOfRange(model,0,1750);
		for(int strand=0; strand<2; strand++){
			final Orf primary=call(model,primaryPartial,strand,0,1749,.8f,true);
			final Orf baseline=call(model,primaryPartial,strand,0,1749,.8f,true,false);
			check(primary!=null && baseline!=null && primary.start==0 && primary.stop==1749
				&& primary.toGff().equals(baseline.toGff()) && primary.ncrnaPartial==baseline.ncrnaPartial,
				"Physical-edge retry must leave an accepted primary's complete GFF and partial metadata unchanged");
		}
		System.out.println("PASS Euk18sContigEdgeClipTest physical_ends/no_internal_restart/window_guard/strands/partial_GFF/minimum_length");
	}
	private static Orf call(byte[] model,byte[] bases,int strand,int start,int stop)throws Exception{
		return call(model,bases,strand,start,stop,1f,false);
	}
	private static Orf call(byte[] model,byte[] bases,int strand,int start,int stop,float cutoff,boolean pacBio)throws Exception{
		return call(model,bases,strand,start,stop,cutoff,pacBio,true);
	}
	private static Orf call(byte[] model,byte[] bases,int strand,int start,int stop,float cutoff,boolean pacBio,boolean physical)throws Exception{
		final NcrnaScavenger s=new NcrnaScavenger(new byte[][]{model},null,new String[]{"edge_fixture"},null,17,916,3109,
			9,1,false,1,0,0,1,0,1,cutoff,cutoff);
		s.family=Euk18sRuntimeConfig.FAMILY;s.outputType=ProkObject.r18S;s.reuseConsensusAlignment=true;
		s.modelClipRescue=true;s.modelClipRescueId=.99f;s.trimAlignmentExtent=false;
		s.physicalEdgeMinLen=physical ? 916 : 0;
		s.pacBioConsensusAlignment=pacBio;s.pacBioRolling=pacBio;
		final Method method=NcrnaScavenger.class.getDeclaredMethod("alignWindow",String.class,byte[].class,int.class,int.class,
			int.class,int.class,int.class,int[].class,long[].class);
		method.setAccessible(true);
		return (Orf)method.invoke(s,"edge_fixture",bases,strand,1,start,stop,0,new int[0],new long[0]);
	}
	private static byte[] random(int length,int seed){
		assert(length>0):"Alignment fixtures need nonempty DNA";
		final byte[] out=new byte[length], alphabet={'A','C','G','T'};final Random rng=new Random(seed);
		for(int i=0; i<length; i++){out[i]=alphabet[rng.nextInt(4)];}return out;
	}
	private static byte[] embed(byte[] sequence,int left,int right){
		assert(left>=0 && right>=0):"Fixture flanks must preserve physical coordinates";
		final byte[] out=random(left+sequence.length+right,9281);
		System.arraycopy(sequence,0,out,left,sequence.length);return out;
	}
	private static void check(boolean condition,String why){if(!condition){throw new AssertionError(why);}}
}
