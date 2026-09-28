package prok;

import java.util.ArrayList;
import java.util.Arrays;

import map.LongHashSet;

/** Checks work counters against callbacks and actual shared-scanner control flow.
 * @author Raiden
 */
public class NcrnaWorkStatsTest {

	public static void main(String[] args){
		testWindowCounters(false);
		testWindowCounters(true);
		testSharedScanCounters();
		System.out.println("PASS NcrnaWorkStatsTest");
	}

	private static void testWindowCounters(boolean refresh){
		final LongHashSet seeds=new LongHashSet(2);seeds.add(0L);//Exactly A^17.
		final NcrnaScavenger scavenger=new NcrnaScavenger(new byte[][]{repeat('A',24)},
			null,null,seeds,17,20,20,7,1,false,0f,0f,0f,1,0f,1f,.9f,.8f);
		scavenger.refreshClaimedWindows=refresh;
		scavenger.minKmerHits=1000;
		final byte[] bases=repeat('N',120);
		Arrays.fill(bases,22,39,(byte)'A');
		Arrays.fill(bases,82,99,(byte)'A');
		final long[] observed=new long[4];
		scavenger.setWorkloadSink(new NcrnaWorkloadInstrumentSink(){
			@Override public void seedHits(String contig,int strand,int[] hits){observed[0]+=hits.length;}
			@Override public void scheduledWindow(String contig,int strand,int pass,int start,int stop){observed[pass]++;}
			@Override public void postSeedWindow(String contig,int strand,int pass,int start,int stop,int hits){observed[3]++;}
		});
		final ArrayList<int[]> claimed=new ArrayList<int[]>();claimed.add(new int[]{40,59});
		scavenger.scavenge("split",bases,0,claimed);
		check(observed[0]==2 && observed[1]==2 && observed[2]==2 && observed[3]==0,
			"Both passes must schedule the two surviving windows before the seed gate: "+Arrays.toString(observed));
		check(scavenger.kmerHitCount()==observed[0] && scavenger.windowCount()==observed[1]+observed[2],
			"Counters disagree with independent callbacks; refresh="+refresh);
		check(scavenger.alignmentCount()==0,"Seed-rejected windows must count without alignments");
		scavenger.setWorkloadSink(null);
		scavenger.scavenge("sinkOff",bases,0,claimed);
		check(scavenger.kmerHitCount()==4 && scavenger.windowCount()==8,
			"Always-on counters must accumulate while the optional sink is off");
		claimed.clear();claimed.add(new int[]{0,119});
		scavenger.scavenge("claimed",bases,0,claimed);
		check(scavenger.kmerHitCount()==6 && scavenger.windowCount()==8,
			"Fully claimed input has seed hits but schedules no windows");
		scavenger.scavenge("short",repeat('A',19),0,new ArrayList<int[]>());
		scavenger.scavenge("noHits",repeat('N',120),0,new ArrayList<int[]>());
		check(scavenger.kmerHitCount()==6 && scavenger.windowCount()==8,
			"Minimum-length and no-hit returns must not invent work");
	}

	private static void testSharedScanCounters(){
		final ArrayList<NcrnaFamily> saved=new ArrayList<NcrnaFamily>(GeneCaller.ncrnaFamilies);
		final boolean trna=ProkObject.calltRNA, s16=ProkObject.call16S, s23=ProkObject.call23S;
		final boolean s5=ProkObject.call5S, s18=ProkObject.call18S;
		try{
			ProkObject.calltRNA=ProkObject.call16S=ProkObject.call23S=ProkObject.call5S=ProkObject.call18S=false;
			GeneCaller.ncrnaFamilies.clear();
			GeneCaller.ncrnaFamilies.add(family("first",17,20));
			GeneCaller.ncrnaFamilies.add(family("second",17,20));
			GeneCaller.initializeConservedRnaSeedIndex();
			final GeneCaller caller=caller(), untouched=caller();
			caller.makeRnas("short",repeat('A',16));
			check(caller.sharedSweepPasses()==0,"Fewer than17 bases return before the shared scan loop");
			caller.makeRnas("belowFamilyMinLen",repeat('A',17));
			check(caller.sharedSweepPasses()==2,"A17 is scanned on both strands even below family minLen");
			check(Arrays.equals(caller.ncrnaKmerHitCounts(),new long[]{0,0}),
				"Family-local hit counters must retain their own minimum-length guard");
			caller.makeRnas("noHits",repeat('N',30));
			check(caller.sharedSweepPasses()==4,"No-hit strands still perform a shared scan");
			caller.makeRnas("sharedHits",repeat('A',30));
			check(caller.sharedSweepPasses()==6,"Two families share two strand scans, not four");
			check(Arrays.equals(caller.ncrnaKmerHitCounts(),new long[]{14,14}),
				"A30 contains14 occurrences of A17, each routed to both families");
			check(Arrays.equals(caller.ncrnaWindowCounts(),new long[]{1,1}),
				"Each family schedules one collapsed window, even when the seed-count filter rejects it");
			check(untouched.sharedSweepPasses()==0 && Arrays.equals(untouched.ncrnaKmerHitCounts(),new long[]{0,0}),
				"A shared index must not share mutable worker counters");
			GeneCaller.ncrnaFamilies.clear();
			GeneCaller.ncrnaFamilies.add(family("legacy16",16,20));
			GeneCaller.initializeConservedRnaSeedIndex();
			final GeneCaller legacy=caller();legacy.makeRnas("legacy",repeat('A',30));
			check(legacy.sharedSweepPasses()==0 && Arrays.equals(legacy.ncrnaKmerHitCounts(),new long[]{15}),
				"A non17-mer family retains its local scanner and must not count a shared pass");
			GeneCaller.ncrnaFamilies.clear();GeneCaller.initializeConservedRnaSeedIndex();
			final GeneCaller disabled=caller();disabled.makeRnas("disabled",repeat('A',30));
			check(disabled.sharedSweepPasses()==0 && disabled.ncrnaWindowCounts().length==0,
				"Disabled families must leave no shared scan or family counters");
		}finally{
			GeneCaller.ncrnaFamilies.clear();GeneCaller.ncrnaFamilies.addAll(saved);
			ProkObject.calltRNA=trna;ProkObject.call16S=s16;ProkObject.call23S=s23;
			ProkObject.call5S=s5;ProkObject.call18S=s18;
			GeneCaller.initializeConservedRnaSeedIndex();
		}
	}

	private static NcrnaFamily family(String name,int k,int minLen){
		final LongHashSet seeds=new LongHashSet(2);seeds.add(0L);
		final NcrnaFamily family=new NcrnaFamily(name,new byte[][]{repeat('A',24)},null,null,
			seeds,k,minLen,50,7,1,false,0f,0f,0f,1,0f,1f,.9f,.8f);
		family.seedMinHits=1000;
		return family;
	}

	private static GeneCaller caller(){return new GeneCaller(20,0,0,0f,0f,0f,0f,0f,new GeneModel(false));}
	private static byte[] repeat(char base,int length){final byte[] bases=new byte[length];Arrays.fill(bases,(byte)base);return bases;}
	private static void check(boolean ok,String reason){if(!ok){throw new AssertionError(reason);}}
}
