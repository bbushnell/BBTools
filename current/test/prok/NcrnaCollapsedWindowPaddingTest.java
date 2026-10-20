package test.prok;

import java.util.ArrayList;
import java.util.Arrays;
import prok.CallGenes;
import prok.GeneCaller;
import prok.NcrnaFamily;
import prok.NcrnaFamilyRegistry;
import prok.NcrnaScavenger;

/** Tests retained-window invariance, short-cluster median geometry and real CLI wiring.
 * @author Raiden */
public final class NcrnaCollapsedWindowPaddingTest {
	public static void main(String[] args){
		final int[] chain={100,110,120,130,140,150,160,170,180,190};
		check(scavenger(0).prepareCandidateWindows(chain,null,1000).isEmpty(),"Explicit zero must retain the historical short-intersection drop");
		for(String id:new String[]{"prok5S","23S","16S","euk5S","tRNA","euk18S","rnasep","tmrna","r58","lsu","sixs_rf00013","sixs_rf01685","srp_small","srp_large"}){
			check(NcrnaFamilyRegistry.get(id).collapsedWindowPad==(id.equals("srp_small")?75:0),
				"Only the measured SRP-small family may default to median-seed padding: "+id);
		}
		final ArrayList<int[]> padded=scavenger(75).prepareCandidateWindows(chain,null,1000);
		check(padded.size()==1 && padded.get(0)[0]==65 && padded.get(0)[1]==232,
			"Lower median seed140 with pad75/K17 must restore exactly[65,232],168inclusive bases");
		final int[] retained={100,110};
		final ArrayList<int[]> old=scavenger(0).prepareCandidateWindows(retained,null,1000);
		final ArrayList<int[]> same=scavenger(75).prepareCandidateWindows(retained,null,1000);
		check(old.size()==1 && same.size()==1 && Arrays.equals(old.get(0),same.get(0)),
			"An already admitted intersection must not be widened or recentered");
		final int[] two={100,110,120,130,140,150,160,170,180,190,500};
		final ArrayList<int[]> split=scavenger(75).prepareCandidateWindows(two,null,1000);
		check(split.size()==2 && split.get(0)[0]==65 && split.get(1)[0]==425 && split.get(1)[1]==592,
			"Fallback must use its own cluster median and preserve the following independent seed window");
		final int[] nearby={100,110,120,130,140,150,160,170,180,190,215};
		final ArrayList<int[]> redundant=scavenger(75).prepareCandidateWindows(nearby,null,1000);
		check(redundant.size()==2 && redundant.get(0)[0]==65 && redundant.get(0)[1]==232
			&& redundant.get(1)[0]==140 && redundant.get(1)[1]==307,
			"Padding can overlap a neighboring retained window; preserve both and measure redundant alignment/duplicate-call effects");
		final NcrnaScavenger edge=scavenger(75);
		final ArrayList<int[]> clipped=edge.prepareCandidateWindows(new int[]{8,18,28,38,48,58,68,78,88,98},null,110);
		check(clipped.size()==1 && clipped.get(0)[0]==0 && clipped.get(0)[1]==109,
			"Short internal intersection must pad around its median then clip at the actual110base contig bounds");
		fails(()->scavenger(-1).prepareCandidateWindows(chain,null,1000));
		fails(()->scavenger(10).prepareCandidateWindows(chain,null,1000));
		final NcrnaScavenger voting=scavenger(75);voting.voteWindows=true;
		fails(()->voting.prepareCandidateWindows(chain,null,1000));
		check(scavenger(75).prepareCandidateWindows(new int[0],null,1000).isEmpty(),"No seeds must remain no windows");
		if(args.length>0){
			check(args[0].startsWith("expectpad="),"Optional CLI fixture begins with expectpad=N");
			final int expected=Integer.parseInt(args[0].substring(10));
			new CallGenes(Arrays.copyOfRange(args,1,args.length));
			int seen=0;
			for(NcrnaFamily f:GeneCaller.ncrnaFamilies){
				if(f.name.equals("srp_small")){
					seen++;check(f.collapsedWindowPad==expected && f.minLen==80 && f.windowPad==75,
						"CLI must change only the requested fallback padding, never the gene floor or original seed pad");
					final NcrnaScavenger actual=scavenger(0);GeneCaller.applyNcrnaFamilyControls(actual,f);
					check(actual.collapsedWindowPad==expected,"GeneCaller must transfer the loaded family's exact control");
					check(actual.prepareCandidateWindows(chain,null,1000).size()==(expected==0?0:1),
						"CLI-selected control must reach the actual production preparation seam");
				}else{check(f.collapsedWindowPad==0,"Other loaded families must retain disabled padding");}
			}
			check(seen==1,"Combined caller fixture must retain exactly one SRP-small family");
			System.out.println("NCRNA_COLLAPSED_PADDING_CLI_PASS pad="+expected+" gene_min=80 other_families_unchanged=true");
		}
		System.out.println("NCRNA_COLLAPSED_PADDING_TEST_PASS default_exact median168 retained_exact separate_clusters physical_edges invalid_controls");
	}
	private static NcrnaScavenger scavenger(int pad){
		final NcrnaScavenger s=new NcrnaScavenger(null,null,null,null,17,80,75);
		s.collapseFrac=.85f;s.physicalEdgeMinLen=80;s.collapsedWindowPad=pad;return s;
	}
	private static void fails(Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid padding policy was accepted");
	}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
