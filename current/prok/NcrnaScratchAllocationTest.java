package prok;

import java.util.ArrayList;
import java.util.Arrays;
import map.LongHashSet;

/** Run in a fresh -Xmx32m JVM. The former12M-float scratch alone exceeds that heap.
 * Exercises both modern tRNA and generic conserved-seed paths on a real long buffer.
 * @author Raiden
 */
public final class NcrnaScratchAllocationTest {
	public static void main(String[] args){
		check(Runtime.getRuntime().maxMemory()<=40L*1024*1024,"Run with -Xmx32m so the fixture actually excludes a whole-contig float scratch");
		ProkObject.call16S=ProkObject.call23S=ProkObject.call5S=ProkObject.call18S=false;
		ProkObject.calltRNA=true;TrnaCaller.SCAVENGE_ONLY=true;
		final byte[] model=new byte[80];Arrays.fill(model,(byte)'A');
		GeneCaller.trnaLibrary=new byte[][]{model};ProkObject.trnaKmers=new LongHashSet(8);ProkObject.trnaKmers.add(0);
		check(GeneCaller.ncrnaFamilies.isEmpty(),"Fixture requires fresh process family state");
		final LongHashSet seeds=new LongHashSet(8);seeds.add(0);
		GeneCaller.ncrnaFamilies.add(new NcrnaFamily("scratch_fixture",new byte[][]{model},null,null,seeds,17,20,50,7,1,false,0f,0f,0f,1,0f,1f,.9f,.8f));
		GeneCaller.initializeConservedRnaSeedIndex();
		final byte[] bases=new byte[12000000];Arrays.fill(bases,(byte)'N');
		final GeneModel pgm=new GeneModel(false);pgm.statstRNA.setInner(1,1);
		final GeneCaller caller=new GeneCaller(20,0,0,0f,0f,0f,0f,0f,pgm);
		final ArrayList<Orf>[] calls=caller.makeRnas("long_nohit",bases);
		check(calls[0].isEmpty() && calls[1].isEmpty() && caller.sharedSweepPasses()==2,"Both actual strand scans must complete without invented calls");
		for(byte b:bases){check(b=='N',"The two strand passes must restore the input");}
		System.out.println("PASS NcrnaScratchAllocationTest 12M_bases 32m_heap modern_trna generic_rna both_strands");
	}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
