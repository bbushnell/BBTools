package test.prok;

import prok.*;

import java.util.ArrayList;
import java.util.Arrays;
import map.LongHashSet;

/** Packaged resource/default binding, explicit overrides and disabled-family isolation.
 * Biological validation remains the whole-genome caller comparison.
 * @author Brian, Raiden, Ganyu */
public final class Prok16sReleaseTest {
	public static void main(String[] args){
		if(ml.Function.TYPE_RATES_CUM==null){ml.Function.normalizeTypeRates();}
		final ArrayList<NcrnaFamily> empty=new ArrayList<NcrnaFamily>();
		new Prok16sEndpointConfig().apply(empty);
		final Prok16sEndpointConfig disabled=new Prok16sEndpointConfig();disabled.parse("rrna16sendpoint","f");disabled.apply(empty);
		final Prok16sEndpointConfig absent=new Prok16sEndpointConfig();absent.parse("rrna16sendpointtables","/missing16s");reject(()->absent.apply(empty));
		CallGenes.setRrna17Profile("16s");CallGenes.loadRrna17Resources("16s");
		check(GeneCaller.ncrnaFamilies.size()==1,"The16S-only profile must not register other seeded rRNA families");
		final NcrnaFamily f=GeneCaller.ncrnaFamilies.get(0);
		checkPackaged(f);
		check(f.name.equals("16S") && f.kLong==17 && f.seedMinHits==2
			&& f.preserveContainingWindows && f.indexK==9 && f.indexTopN==1 && f.minLen==1400
			&& f.idPass==.66f && f.idBorderline==.66f && f.prok16sPartialMinLen==1200,
			"Packaged caller must bind the measured seeds, windows, identity and physical-edge floor");
		final Prok16sWindowConfig raw=new Prok16sWindowConfig();raw.parse("rrna16spartialminlen","0");raw.apply(GeneCaller.ncrnaFamilies);
		check(f.prok16sPartialMinLen==0 && f.preserveContainingWindows,"Partial override must preserve the omitted selected window default");
		final Prok16sWindowConfig legacy=new Prok16sWindowConfig();legacy.parse("rrna16swindow","legacy");legacy.apply(GeneCaller.ncrnaFamilies);
		check(!f.preserveContainingWindows && f.prok16sPartialMinLen==0,"Window override must preserve the omitted floor");
		final Prok16sWindowConfig selected=new Prok16sWindowConfig();selected.parse("rrna16swindow","containment");selected.parse("rrna16spartialminlen","1200");selected.apply(GeneCaller.ncrnaFamilies);
		for(String floor:new String[]{"-1","1401","NaN"}){reject(()->new Prok16sWindowConfig().parse("rrna16spartialminlen",floor));}
		reject(()->new Prok16sWindowConfig().parse("rrna16swindow","vote"));
		disabled.apply(GeneCaller.ncrnaFamilies);check(f.rrnaEndpointFeatures==null,"Explicitly disabled endpoint processing must load no resources");
		new Prok16sEndpointConfig().apply(GeneCaller.ncrnaFamilies);
		final RrnaEndpointCallerFeatures.Resources r=f.rrnaEndpointFeatures;
		check(r!=null && r.acceptedVerifierIdentity && r.inputs==28 && r.sites==25,"Packaged nets must use the trained accepted-identity contract");
		int nets=0,rawModels=0;for(int i=0; i<r.names.length; i++){
			check(r.names[i].equals(RrnaResourceIO.first(f.modelNames[i])) && r.raw[i]==(i>=4),"Endpoint indices must preserve the original four refined models and all raw mitochondrial models");
			if(r.raw[i]){check((i==4 ? r.names[i].equals("mitochondria") : r.names[i].startsWith("mitochondria_coverage_")) && r.fiveNets[i]==null && r.threeNets[i]==null,"Only original or coverage-selected mitochondrial16S retains raw endpoints without nets");rawModels++;}
			else{check(r.five[i].counts==null && r.three[i].counts==null && r.fiveNets[i]!=null && r.threeNets[i]!=null,"Production must load both scores-only tables and nets for each refined model");nets+=2;}
		}
		check(r.names.length==f.library.length && nets==8 && rawModels==f.library.length-4,"Every packaged model binds either its original endpoint pair or explicit raw policy");
		final Prok16sEndpointConfig missing=new Prok16sEndpointConfig();missing.parse("rrna16sendpointnets","/missing16s.bbnets");reject(()->missing.apply(GeneCaller.ncrnaFamilies));
		partialCds();
		System.out.println("PASS Prok16sReleaseTest packaged_defaults overrides disabled_missing_resources accepted_identity raw_mito scores_only partial_CDS");
	}
	static void checkPackaged(NcrnaFamily f){
		check(f.library.length==5 || f.library.length==158,"Only the original five-model or selected153-addition package is supported");
		// The loader checks the registry; this independent literal also prevents a
		// changed registry and changed resource from silently agreeing on a new prefix.
		for(int i=0;i<5;i++){check(RrnaResourceIO.first(f.modelNames[i]).equals(ORIGINAL_NAMES[i]),"Added models must not replace or reorder an original model");}
		final LongHashSet old=originalSeeds();
		for(long seed:old.toArray()){check(f.kmerSet.contains(seed),"Every original forward17mer must remain in the packaged seed set: "+seed);}
		check(f.kmerSet.size()==(f.library.length==5?285:680),"Model count and exact selected seed count must belong to the same resource generation");
		if(f.library.length==158){
			check(MitoSsuModelPolicy.coverage(f.name,f.modelNames) && f.modelThresholds!=null,"All153 additions require the explicit coverage policy and length bounds");
			for(int i=0;i<5;i++){check(f.modelThresholds.minimum(i)==0 && f.modelThresholds.maximum(i)==Integer.MAX_VALUE,"Coverage length limits must not alter original model floors");}
			for(int i=5;i<158;i++){check(f.modelThresholds.minimum(i)>0 && f.modelThresholds.maximum(i)>=f.modelThresholds.minimum(i),"Every added model needs a valid measured length interval");}
		}
	}
	static LongHashSet originalSeeds(){
		final LongHashSet result=new LongHashSet(512);final long mask=(1L<<34)-1;
		for(String sequence:ORIGINAL_SEED_CONTIGS){long key=0;int length=0;
			for(int i=0;i<sequence.length();i++){final int base=RrnaResourceIO.base((byte)sequence.charAt(i));
				if(base<0){key=0;length=0;}else{key=((key<<2)|base)&mask;if(++length>=17){result.add(key);}}
			}
		}
		check(result.size()==285,"Literal pre-expansion seed fixture must reconstruct exactly285 distinct forward17mers");return result;
	}
	private static final String[] ORIGINAL_NAMES={"universal","bacteria","archaea","plastid","mitochondria"};
	// Exact four contigs from the pre-expansion packaged seeds.fa.gz, sha80
	// 1ba6638d6cf944e0cbf2. N separators forbid synthetic junction seed keys.
	private static final String[] ORIGINAL_SEED_CONTIGS={
		"TGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTNNGCGG"
		+"TGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCANNAGCAAACAGGATTAGATACCCTGGTA"
		+"GTCCACGCCGTAAACGNNACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGNNGAATACGTC"
		+"CCTGCTCCTTGCACACACCGCCCGTCANNCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGA",
		"TGAAACTTAAAGGAATTGGCGGGGGAGCACNNCGAACCGGATTAGATACCCGGGTAGTCCNNTAAAACTC"
		+"AAAGGAATAGACGGGGACNNGAATTGACGGGGGCCCGCACAAGCGGNNGGAGGAAGGTGGGGATGACGTC"
		+"AANNCGTGCCAGCAGCCGCGGTAATACGNNGGGCTACACACGTGCTACAATGGNNAAACTCAAAGGAATT"
		+"GACGGGGNNTGGGACTGAGACACGGCCCANNGCCAGCCGCCGCGGTAATAC",
		"GCAACGCGAAGAACCTTACCNNGTGCCAGCAGCAGCGGTAANNGAGAGGAGGTGCATGGCCGNNAAAAAC"
		+"CTCACCAGGGTTTNNGTGTAGCGGTGAAATGCGNNGCTGGCGGCGTGGATAAGNNGATAACCTCGGGAAA"
		+"CTGNNCAGGATGAACGCTGGCGGNNCAATGCGCGAAAGCGTGANNAACGCTGGCGGCGTGGATNNAACCT"
		+"CGGGAAACTGAGGNNGTAGGGGTGAAATCCGTNNGGAATCGCTAGTAATCGNNGAAGTCGTAACAAGGTA",
		"CCGCCTGGGGAGTACGGNNACACACCGCCCGTCAAANNAAACTACTGTGGTAGGGNNAAAAACCCTTTCTCAGT"
	};
	static void partialCds(){
		final byte[] bases=new byte[900];Arrays.fill(bases,(byte)'A');
		for(int frame=0; frame<3; frame++){for(int strand=0; strand<2; strand++){
			final ArrayList<Orf> calls=GeneCaller.makeOrfsForFrame("edge",bases,frame,strand,60,null);
			check(calls.size()==1 && calls.get(0).toGff().endsWith(";partial=11"),"Actual physical-edge CDS generators must emit Prodigal partial bits on either strand/frame");
			final Orf internal=new Orf("internal",100,799,strand,frame,bases,true,ProkObject.CDS,true,true);
			check(!internal.toGff().contains("partial="),"Internal missing-start/stop evidence is not physical contig truncation");
		}}
	}
	static void reject(Runnable action){boolean failed=false;try{action.run();}catch(RuntimeException e){failed=true;}check(failed,"Invalid explicit resource/configuration must fail loudly");}
	static void check(boolean value,String why){if(!value){throw new AssertionError(why);}}
}
