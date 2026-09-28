package prok;

import java.lang.reflect.Constructor;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import dna.AminoAcid;
import map.LongHashSet;
import map.LongObjectMap;

/** Recovery B55-B58: guard bounds, real work counts, length admission and family isolation.
 * Runs alone in a fresh JVM because CallGenes configuration is process-global.
 * @author Raiden
 */
public class NcrnaRecoveryControlsTest {
	private static final byte[] CONS="ACGTTGCAAGTCGATCGTACGATGC".getBytes(StandardCharsets.US_ASCII);

	public static void main(String[] args) throws Exception{
		testValidation();testCosts();testMaxLen();testFamilyWiring();testOverrides();
		System.out.println("PASS NcrnaRecoveryControlsTest");
	}

	private static void testValidation(){
		for(int k:new int[]{-1,0,16,32,Integer.MAX_VALUE}){
			reject(()->new TrnaKmerIndex(new byte[][]{CONS},k,false,0f,0f,0f,1));
		}
		for(float bad:new float[]{-1f,Float.NaN,Float.POSITIVE_INFINITY,Float.NEGATIVE_INFINITY}){
			reject(()->new TrnaKmerIndex(new byte[][]{CONS},7,false,bad,0f,0f,1));
			reject(()->new TrnaKmerIndex(new byte[][]{CONS},7,false,0f,bad,0f,1));
			reject(()->new TrnaKmerIndex(new byte[][]{CONS},7,false,0f,0f,bad,1));
		}
		reject(()->new TrnaKmerIndex(new byte[][]{CONS},7,false,0f,0f,0f,-1));
		TrnaKmerIndex.validateConfiguration(15,0f,2f,2f,0);
		final TrnaKmerIndex index=new TrnaKmerIndex(new byte[][]{CONS},7,false,0f,0f,0f,1);
		reject(()->index.shortlist(CONS,0));reject(()->index.shortlist(CONS,-1));
		check(index.queriesProcessed()==0,"Invalid topN must not mutate query statistics");
		check(Arrays.equals(index.shortlist(CONS,1),new int[]{0}),"Valid one-model default rejected exact sequence");
		reject(()->new NcrnaScavenger(null,null,null,null,17,20,50,7,0,false,0f,0f,0f,1));
		reject(()->new NcrnaScavenger(null,null,null,null,17,20,50,16,1,false,0f,0f,0f,1));
		reject(()->CallGenes.parseSweepInt("topn","4294967297",1,Integer.MAX_VALUE));
		for(String bad:new String[]{"0","-1","1001","NaN","Infinity"}){
			reject(()->CallGenes.parseSweepInt("r58minlen",bad,1,1000));
		}
	}

	private static void testCosts(){
		for(boolean reuse:new boolean[]{false,true}){
			for(boolean quantum:new boolean[]{false,true}){
				final NcrnaScavenger s=scavenger(new byte[][]{CONS},false,reuse,.9f,.8f);
				s.quantumThresh=quantum ? 0 : 120;
				check(run(s,CONS).size()==1,"Exact sequence must be accepted");
				check(s.alignmentCount()==1 && s.alignedBases()==CONS.length,"Detection must count one actual input span");
				check(s.hbmScoreCalls()==0 && s.hbmBasesScored()==0,"No model means no HBM work");
			}
			final byte[] mutant=CONS.clone();mutant[mutant.length-1]='A';
			final NcrnaScavenger hbm=scavenger(new byte[][]{CONS},true,reuse,.99f,.5f);hbm.hbmPass=.5f;
			final ArrayList<Orf> rescued=run(hbm,mutant);
			check(rescued.size()==1,"Borderline locus must exercise HBM rescue");
			final int scoredSpan=rescued.get(0).length();
			check(hbm.hbmScoreCalls()==1 && hbm.hbmBasesScored()==scoredSpan,
				"One scored ORF must count once: reuse="+reuse+" calls="+hbm.hbmScoreCalls()+" bases="+hbm.hbmBasesScored()+" span="+scoredSpan);
			check(hbm.alignmentCount()==(reuse ? 1 : 2) && hbm.alignedBases()==CONS.length+(reuse ? 0 : scoredSpan),
				"Reused traceback must not add an alignment; legacy rescue must count its real alignment");
		}
		final NcrnaScavenger trim=scavenger(new byte[][]{CONS},false,false,.9f,.8f);
		trim.trimAlignmentExtent=true;
		check(run(trim,CONS).size()==1 && trim.alignmentCount()==2 && trim.alignedBases()==2*CONS.length,
			"Endpoint alignment must count its input once");
	}

	private static void testMaxLen() throws Exception{
		final byte[] longer="ACGTTGCAAGTCGATCGTACGATGCTTAGGCATCGAT".getBytes(StandardCharsets.US_ASCII);
		final NcrnaScavenger s=scavenger(new byte[][]{longer,CONS},false,true,.9f,.8f);
		s.maxLen=CONS.length;s.rankedModelFallback=true;
		final ArrayList<Orf> result=run(s,longer);
		check(result.size()==1 && "model1".equals(result.get(0).trnaModel),
			"Overlong rank-one model must not prevent later short-model acceptance");
		check(s.alignmentCount()==2 && s.alignedBases()==2*longer.length,"Must try both models once");
		final NcrnaScavenger rejected=scavenger(new byte[][]{longer},false,true,.9f,.8f);rejected.maxLen=CONS.length;
		check(run(rejected,longer).isEmpty(),"An overlong sole model must be rejected");
		final NcrnaScavenger voted=scavenger(new byte[][]{CONS},false,true,.9f,.8f);
		voted.maxLen=CONS.length;voted.voteEnds=true;
		final LongObjectMap<SeedOffsetTable.KmerInfo> map=new LongObjectMap<SeedOffsetTable.KmerInfo>(2,SeedOffsetTable.KmerInfo.class);
		map.put(seed(),new SeedOffsetTable.KmerInfo(0,19,0,0,1));
		final Constructor<SeedOffsetTable> ctor=SeedOffsetTable.class.getDeclaredConstructor(int.class,long.class,String.class,LongObjectMap.class);
		ctor.setAccessible(true);voted.voteTable=ctor.newInstance(17,1L,null,map);
		final byte[] flanked=Arrays.copyOf(CONS,CONS.length+20);Arrays.fill(flanked,CONS.length,flanked.length,(byte)'N');
		check(run(voted,flanked).isEmpty() && voted.votedStops==1,"Length cap must also reject expansion from endpoint voting");
	}

	private static void testFamilyWiring(){
		ProkObject.calltRNA=ProkObject.call16S=ProkObject.call23S=ProkObject.call5S=ProkObject.call18S=false;
		GeneCaller.ncrnaFamilies.clear();
		final NcrnaFamily family=new NcrnaFamily("fixture",new byte[][]{CONS},null,new String[]{"model0"},
			seeds(),17,20,100,7,1,false,0f,0f,0f,1,0f,1f,.9f,.8f);
		family.reuseConsensusAlignment=true;family.trimAlignmentExtent=false;family.scavengePass2=false;
		family.maxLen=CONS.length;
		GeneCaller.ncrnaFamilies.add(family);GeneCaller.initializeConservedRnaSeedIndex();
		final GeneCaller caller=new GeneCaller(20,0,0,0f,0f,0f,0f,0f,new GeneModel(false));
		check(caller.makeRnas("fixture",CONS.clone())[0].size()==1,"Family controls must reach the lazy scavenger");
		check(Arrays.equals(caller.ncrnaAlignmentCounts(),new long[]{1}) && Arrays.equals(caller.ncrnaAlignedBases(),new long[]{CONS.length}),
			"Worker aggregation must match actual family alignment work");
		check(caller.ncrnaHbmScoreCalls()[0]==0 && caller.ncrnaHbmBasesScored()[0]==0,"Worker must not invent HBM work");
		GeneCaller.ncrnaFamilies.clear();GeneCaller.initializeConservedRnaSeedIndex();
	}

	private static void testOverrides(){
		CallGenes.NCRNA_FAMILIES_ENABLED=true;CallGenes.R58LSU_ENABLED=true;
		CallGenes.NCRNA_BOUNDARY_NN_ENABLED=false;
		CallGenes.loadNcrnaResources();
		final NcrnaFamily originalR58=find("r58"),originalLsu=find("lsu"),originalP=find("rnasep");
		check(originalR58.minLen==140 && originalR58.indexK==7 && originalLsu.indexK==9,"Packaged defaults changed");
		CallGenes.NCRNA_FAMILY_FILTER="r58";CallGenes.NCRNA_INDEX_K_OVERRIDE=6;
		CallGenes.NCRNA_INDEX_TOP_N_OVERRIDE=3;CallGenes.NCRNA_INDEX_MIN_HITS_OVERRIDE=4;
		CallGenes.NCRNA_ADAPTIVE_MINHITS_OVERRIDE=true;CallGenes.NCRNA_ADAPT_FLOOR_OVERRIDE=2f;
		CallGenes.NCRNA_ADAPT_TOPFRAC_OVERRIDE=.3f;CallGenes.NCRNA_ADAPT_QFRAC_OVERRIDE=.2f;
		CallGenes.R58_MINLEN_OVERRIDE=130;
		check(CallGenes.parseFamilyIndexOverride("R58NCRNAMAPK","8"),"Historical case-insensitive mapk alias lost");
		check(CallGenes.parseFamilyIndexOverride("lsuncrnafixedminhits","449"),"Historical minhits alias lost");
		CallGenes.validateNcrnaSweepOverrides();GeneCaller.ncrnaFamilies.clear();CallGenes.loadNcrnaResources();
		final NcrnaFamily r=find("r58"),l=find("lsu"),p=find("rnasep");
		check(r.minLen==130 && r.indexK==8 && r.indexTopN==3 && r.fixedMinHits==4 && r.adaptive
			&& r.adaptFloor==2f && r.adaptTopFrac==.3f && r.adaptQFrac==.2f,"Specific override must beat targeted generic configuration");
		check(l.indexK==originalLsu.indexK && l.indexTopN==originalLsu.indexTopN && l.fixedMinHits==449
			&& l.minLen==originalLsu.minLen && !l.adaptive,"R58 controls leaked to LSU");
		check(p.indexK==originalP.indexK && p.fixedMinHits==originalP.fixedMinHits && p.minLen==originalP.minLen,
			"Overrides leaked to RNase P");
		CallGenes.NCRNA_FAMILY_FILTER=null;reject(()->CallGenes.validateNcrnaSweepOverrides());
		CallGenes.NCRNA_FAMILY_FILTER="r58";CallGenes.R58LSU_ENABLED=false;reject(()->CallGenes.validateNcrnaSweepOverrides());
	}

	private static NcrnaFamily find(String name){
		for(NcrnaFamily f:GeneCaller.ncrnaFamilies){if(f.name.equals(name)){return f;}}
		throw new AssertionError("Missing family "+name);
	}
	private static NcrnaScavenger scavenger(byte[][] library,boolean hbm,boolean reuse,float pass,float borderline){
		final String[] names=new String[library.length];final BaseGraph[] models=hbm ? new BaseGraph[library.length] : null;
		for(int i=0;i<library.length;i++){
			names[i]="model"+i;
			if(hbm){models[i]=new BaseGraph(names[i],library[i],null,i,0);BaseGraphHelper.initForScoring(models[i]);}
		}
		final NcrnaScavenger s=new NcrnaScavenger(library,models,names,seeds(),17,20,100,7,library.length,false,0f,0f,0f,1,0f,1f,pass,borderline);
		s.family="fixture";s.reuseConsensusAlignment=reuse;s.trimAlignmentExtent=false;s.scavengePass2=false;
		return s;
	}
	private static long seed(){long seed=0;for(int i=0;i<17;i++){seed=(seed<<2)|AminoAcid.baseToNumber[CONS[i]];}return seed;}
	private static LongHashSet seeds(){final LongHashSet set=new LongHashSet(2);set.add(seed());return set;}
	private static ArrayList<Orf> run(NcrnaScavenger s,byte[] bases){return s.scavenge("fixture",bases,0,new ArrayList<int[]>());}
	private static void reject(Runnable r){try{r.run();}catch(IllegalArgumentException expected){return;}throw new AssertionError("Invalid configuration was accepted");}
	private static void check(boolean ok,String message){if(!ok){throw new AssertionError(message);}}
}
