package prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import map.LongHashSet;

/** Model-specific rejection must not conceal a different passing model.
 * @author Raiden
 */
public final class NcrnaModelThresholdsTest {
	public static void main(String[] args)throws Exception{
		final Path dir=Files.createTempDirectory("ncrna_cutoffs_");final Path table=dir.resolve("cutoffs.tsv");
		final byte[] gene=new byte[120];final Random random=new Random(74063);final byte[] alphabet={'A','C','G','T'};
		for(int i=0;i<gene.length;i++){gene[i]=alphabet[random.nextInt(4)];}
		final byte[] high=gene.clone(),low=gene.clone();for(int i=40;i<120;i+=10){high[i]=other(high[i]);}for(int i=40;i<120;i+=5){low[i]=other(low[i]);}
		final byte[] bases=new byte[200];Arrays.fill(bases,(byte)'N');System.arraycopy(gene,0,bases,40,120);
		final NcrnaFamily family=family(high,low,gene),unaffected=family(high,low,gene);final ArrayList<NcrnaFamily> families=new ArrayList<NcrnaFamily>();families.add(family);
		final ArrayList<Orf> baseline=call(family,bases);check(baseline.size()==1 && baseline.get(0).trnaModel.equals("high"),"Ordinary scalar cutoff must choose the higher-identity model in this fixture");
		write(table,"euk5S\tlow\t.6\t.6\neuk5S\thigh\t1\t1\n");NcrnaModelThresholds.load(table.toString(),families);
		final ArrayList<Orf> specific=call(family,bases);check(specific.size()==1 && specific.get(0).trnaModel.equals("low"),"Higher-identity model failing its own cutoff cannot hide a different passing model; bind by name, not table row order");
		write(table,"euk5S\thigh\t.68\t.68\neuk5S\tlow\t.68\t.68\n");NcrnaModelThresholds.load(table.toString(),families);
		check(RrnaEndpointCallerFeaturesTest.calls(baseline).equals(RrnaEndpointCallerFeaturesTest.calls(call(family,bases))),"Common .68 must still admit both fixture models and preserve the selected result");
		check(unaffected.modelThresholds==null && RrnaEndpointCallerFeaturesTest.calls(baseline).equals(RrnaEndpointCallerFeaturesTest.calls(call(unaffected,bases))),"Other family instances retain their scalar policy");
		final NcrnaModelThresholds good=family.modelThresholds;
		for(String invalid:new String[]{"euk5S\thigh\t.7\t.7\n","euk5S\thigh\t.6\t.7\neuk5S\tlow\t.6\t.6\n","euk5S\thigh\tNaN\t.6\neuk5S\tlow\t.6\t.6\n","euk5S\thigh\t.6\t.6\neuk5S\thigh\t.6\t.6\n","unknown\thigh\t.6\t.6\n"}){
			write(table,invalid);reject(()->NcrnaModelThresholds.load(table.toString(),families));check(family.modelThresholds==good,"Malformed cutoff files must not partially mutate the live family");
		}
		final NcrnaScavenger nativeCaller=RrnaEndpointCallerFeaturesTest.scavenger(unaffected);
		reject(()->nativeCaller.setModelThresholds(new NcrnaModelThresholds(new String[]{"low","high"},new float[]{.6f,.6f},new float[]{.6f,.6f})));
		family.reuseConsensusAlignment=false;write(table,"euk5S\thigh\t.6\t.6\neuk5S\tlow\t.6\t.6\n");reject(()->NcrnaModelThresholds.load(table.toString(),families));
		System.out.println("PASS NcrnaModelThresholdsTest actual_alignment_model_rejection name_binding common68 scalar_unchanged atomic_invalid_rows unsupported_path");
	}
	static NcrnaFamily family(byte[] high,byte[] low,byte[] gene){
		final LongHashSet seeds=new LongHashSet(32);for(int start=0;start<20;start++){long word=0;for(int j=0;j<17;j++){word=(word<<2)|dna.AminoAcid.baseToNumber[gene[start+j]];}seeds.add(word);}
		final NcrnaFamily f=new NcrnaFamily("euk5S",new byte[][]{high,low},null,new String[]{"high","low"},seeds,17,90,150,7,2,false,1,0,0,1,0,1,.6f,.6f);
		f.reuseConsensusAlignment=true;f.trimAlignmentExtent=false;f.scavengePass2=false;return f;
	}
	static ArrayList<Orf> call(NcrnaFamily family,byte[] bases){final NcrnaScavenger s=RrnaEndpointCallerFeaturesTest.scavenger(family);s.setMappingOracleExhaustive(true);return s.scavenge("contig",bases,0,new ArrayList<int[]>());}
	static byte other(byte b){return b=='A'?(byte)'C':(byte)'A';}
	static void write(Path path,String rows)throws Exception{Files.write(path,("family\tmodel\tidpass\tidborderline\n"+rows).getBytes(StandardCharsets.US_ASCII));}
	static void reject(Runnable r){try{r.run();}catch(IllegalArgumentException e){return;}throw new AssertionError("Expected invalid model cutoff rejection");}
	static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
