package test.prok;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import prok.GeneCaller;
import prok.GeneModel;
import prok.Orf;
import prok.ProkObject;

/** RNA has no coding frame; retain all ordinary predecessor geometry limits.
 * Exercises both real DP predecessor paths, not just the edge predicate.
 * @author Brian, Raiden */
public final class NcrnaPathOverlapTest {
	public static void main(String[] args)throws Exception{
		check(args.length==0,"This finite fixture takes no external data or configuration");
		final byte[] bases=new byte[500];Arrays.fill(bases,(byte)'A');
		for(int strand=0;strand<2;strand++){
			final Orf cds=orf(bases,100,199,strand,0,ProkObject.CDS);
			for(int type:new int[]{ProkObject.RNA,ProkObject.tRNA,ProkObject.r5S,ProkObject.r16S,ProkObject.r23S,ProkObject.r18S}){
				final Orf rna=orf(bases,198,293,strand,0,type);
				check(rna.isValidPrev(cds,2),"Two-base CDS/RNA overlap fits the explicit limit; placeholder frame0 is not coding evidence");
				check(!rna.isValidPrev(cds,1),"RNA must still obey the ordinary overlap limit");
				check(!rna.isValidPrev(orf(bases,198,250,strand,1,ProkObject.CDS),10),"Equal starts cannot become a predecessor edge");
				check(!rna.isValidPrev(orf(bases,100,293,strand,1,ProkObject.CDS),100),"Equal stops cannot become a predecessor edge");
				check(orf(bases,292,400,strand,0,ProkObject.CDS).isValidPrev(rna,2),"The RNA exemption must also apply when RNA is the predecessor");
			}
			final Orf same=orf(bases,198,293,strand,0,ProkObject.CDS);
			check(!same.isValidPrev(cds,10),"Overlapping same-frame CDS/CDS remains forbidden");
			check(orf(bases,200,293,strand,0,ProkObject.CDS).isValidPrev(cds,0),"Nonoverlapping same-frame CDS remains allowed");
			check(orf(bases,198,293,strand,1,ProkObject.CDS).isValidPrev(cds,2),"Different-frame CDS overlap retains its prior behavior");
			check(orf(bases,198,293,1-strand,0,ProkObject.CDS).isValidPrev(cds,2),"Opposite-strand CDS is not a same-frame edge");
			path(bases,strand,ProkObject.RNA,true);
			path(bases,strand,ProkObject.CDS,false);
		}
		System.out.println("NCRNA_PATH_OVERLAP_TEST_PASS both_strands actual_DP_RNA_retained same_frame_CDS_rejected both_orders all_RNA_types overlap_and_order_limits");
	}
	private static void path(byte[] bases,int strand,int nextType,boolean both)throws Exception{
		assert(strand==0 || strand==1) : "The fixture exercises the actual plus and minus predecessor tables";
		final GeneCaller caller=new GeneCaller(40,2,2,0,0,0,0,0,new GeneModel(false));
		@SuppressWarnings("unchecked")
		final ArrayList<Orf>[] frames=new ArrayList[6];
		for(int i=0;i<frames.length;i++){frames[i]=new ArrayList<Orf>();}
		final Orf cds=orf(bases,100,199,strand,0,ProkObject.CDS),next=orf(bases,198,293,strand,0,nextType);
		cds.orfScore=1000;next.orfScore=100;frames[3*strand].add(cds);frames[3*strand].add(next);
		final Method method=GeneCaller.class.getDeclaredMethod("findPath",ArrayList[].class,byte[].class);method.setAccessible(true);
		@SuppressWarnings("unchecked")
		final ArrayList<Orf> result=(ArrayList<Orf>)method.invoke(caller,frames,bases);
		check(result.size()==(both?2:1) && result.contains(cds) && result.contains(next)==both,
			"Actual DP must retain the positive RNA/CDS pair but reject the competing same-frame CDS");
	}
	private static Orf orf(byte[] bases,int start,int stop,int strand,int frame,int type){
		assert(start>=0 && stop<bases.length && stop-start>=2) : "Fixture coordinates must support the constructor's codon reads";
		return new Orf("fixture",start,stop,strand,frame,bases,false,type);
	}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
