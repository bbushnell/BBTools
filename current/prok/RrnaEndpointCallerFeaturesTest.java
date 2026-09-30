package prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import map.LongHashSet;
import parse.LineParser1;
import structures.ByteBuilder;

/** Actual family/scavenger plumbing plus byte-for-byte materializer feature parity.
 * @author Raiden
 */
public final class RrnaEndpointCallerFeaturesTest {
	public static void main(String[] args) throws Exception{
		final Path dir=Files.createTempDirectory("rrna_caller_features_");final byte[] bases=new byte[180];Arrays.fill(bases,(byte)'A');
		System.arraycopy(bytes("CGTACG"),0,bases,47,6);System.arraycopy(bytes("TGCATG"),0,bases,125,6);final byte[] original=bases.clone(),consensus=Arrays.copyOfRange(bases,50,130);
		final RrnaPositionalKmerTable.Table five=new RrnaPositionalKmerTable.Table(MODEL,"5prime",6,-3,.5),three=new RrnaPositionalKmerTable.Table(MODEL,"3prime",6,-4,.5);
		for(int i=0;i<100;i++){five.add(bases,50);three.add(bases,129);}five.finish();three.finish();
		final Recorder recorder=new Recorder();final NcrnaFamily family=family(consensus);final NcrnaScavenger off=scavenger(family);
		off.captureRrnaEndpointFeatures(null,null,-1);check(off.alignmentCount()==0,"Default-off hook must not inspect inputs or align");
		family.setRrnaEndpointFeatures(new RrnaPositionalKmerTable.Table[]{five},new RrnaPositionalKmerTable.Table[]{three},recorder);
		final NcrnaScavenger on=scavenger(family);final Orf candidate=candidate(bases,53,125,0);on.captureRrnaEndpointFeatures(candidate,bases,0);
		check(recorder.count==1 && recorder.fiveUsable && recorder.threeUsable && on.alignmentCount()==1 && on.alignedBases()==73,"One candidate produces two features using one counted candidate-span alignment");
		check(candidate.start==53 && candidate.stop==125 && Arrays.equals(bases,original),"Feature observation must not mutate raw ends or caller bases");
		final byte[] expectedFive=recorder.five.clone(),expectedThree=recorder.three.clone();
		final float gc=shared.Tools.calcGC(bases);final Path input=dir.resolve("loci.fa");
		Files.write(input,bytes(">a lflank=50 rflank=50 gccontig="+Float.toString(gc)+"\n"+new String(bases,StandardCharsets.US_ASCII)+"\n"));
		final RrnaEndpointVector.Distribution[] hist={new RrnaEndpointVector.Distribution(),new RrnaEndpointVector.Distribution()};hist[0].add(3,1);hist[1].add(-4,1);
		final String prefix=dir.resolve("training").toString();new RrnaEndpointVector(MODEL,consensus,five,three,hist,74003,1).run(input.toString(),prefix,1);
		check(Arrays.equals(expectedFive,features(prefix+".5prime.tsv.gz")) && Arrays.equals(expectedThree,features(prefix+".3prime.tsv.gz")),"All28 serialized caller feature columns must be byte-identical to materializer rows at nonzero anchors");
		final RrnaPositionalKmerTable.Table mixedThree=new RrnaPositionalKmerTable.Table(MODEL,"3prime",9,-8,.5);
		for(int i=0;i<100;i++){mixedThree.add(bases,129);}mixedThree.finish();
		final Recorder mixedRecorder=new Recorder();final NcrnaFamily mixedFamily=family(consensus);
		mixedFamily.setRrnaEndpointFeatures(new RrnaPositionalKmerTable.Table[]{five},new RrnaPositionalKmerTable.Table[]{mixedThree},mixedRecorder);
		final NcrnaScavenger mixedCaller=scavenger(mixedFamily);mixedCaller.captureRrnaEndpointFeatures(candidate(bases,53,125,0),bases,0);
		final Path fivePath=dir.resolve("mixed.5prime.tsv.gz"),threePath=dir.resolve("mixed.3prime.tsv.gz"),refPath=dir.resolve("refs.fa"),histPath=dir.resolve("hist.tsv");
		five.write(fivePath.toString());mixedThree.write(threePath.toString());
		Files.write(refPath,bytes('>'+MODEL+'\n'+new String(consensus,StandardCharsets.US_ASCII)+'\n'));
		Files.write(histPath,bytes("model\tend\tsigned_error\tcount\n"+MODEL+"\t5prime\t3\t1\n"+MODEL+"\t3prime\t-4\t1\n"));
		final String mixedPrefix=dir.resolve("mixed_training").toString();
		RrnaEndpointVector.main(new String[]{"in="+input,"model="+MODEL,"consensus="+refPath,"starttable="+fivePath,"stoptable="+threePath,"histogram="+histPath,"outprefix="+mixedPrefix,"expected=1","samples=1","seed=74003"});
		check(mixedRecorder.fiveUsable && mixedRecorder.threeUsable && mixedCaller.alignmentCount()==1
			&& Arrays.equals(mixedRecorder.five,features(mixedPrefix+".5prime.tsv.gz")) && Arrays.equals(mixedRecorder.three,features(mixedPrefix+".3prime.tsv.gz")),"Mixed k6/k9 and independent anchors must preserve all28 training/caller feature values at each end");
		check(RrnaEndpointVector.endpointFromClass(53,9)==50 && RrnaEndpointVector.endpointFromClass(125,16)==129,"Class decoding must recover gene ends rather than k-mer anchor starts");
		final byte[] lower=bytes(new String(bases,StandardCharsets.US_ASCII).toLowerCase(java.util.Locale.ROOT)),lowerOriginal=lower.clone();
		on.captureRrnaEndpointFeatures(candidate(lower,53,125,1),lower,0);
		check(recorder.strand==1 && Arrays.equals(expectedFive,recorder.five) && Arrays.equals(expectedThree,recorder.three) && Arrays.equals(lower,lowerOriginal),"Reverse-pass sense coordinates and soft masking must preserve features without mutating input");
		on.captureRrnaEndpointFeatures(candidate(bases,0,70,0),bases,0);check(!recorder.fiveUsable,"Clipped context must be marked unscorable, never silently padded or zero-filled");
		final long before=on.alignmentCount(),beforeRows=recorder.count;final ArrayList<Orf> expected=off.scavenge("contig",bases,0,new ArrayList<int[]>()),actual=on.scavenge("contig",bases,0,new ArrayList<int[]>());
		check(!actual.isEmpty() && calls(expected).equals(calls(actual)) && recorder.count>beforeRows,"Real accepted-candidate hook must fire without changing calls");
		check(on.alignmentCount()-before==off.alignmentCount()+recorder.count-beforeRows,"Feature alignments must remain visible in actual caller workload");
		boolean failed=false;family.trimAlignmentExtent=true;try{scavenger(family);}catch(IllegalArgumentException e){failed=true;}check(failed,"Raw-endpoint feature contract must reject trimming-enabled configuration");
		System.out.println("PASS RrnaEndpointCallerFeaturesTest actual_family_hook byte_identical_28cols mixed_k6_k9 nonzero_anchors strands lowercase clipping default_off calls_unchanged alignment_accounting");
	}
	static NcrnaFamily family(byte[] consensus){
		final LongHashSet seeds=new LongHashSet(32);for(int i=0;i+17<=consensus.length;i++){long word=0;for(int j=0;j<17;j++){word=(word<<2)|dna.AminoAcid.baseToNumber[consensus[i+j]];}seeds.add(word);}
		final NcrnaFamily f=new NcrnaFamily("euk5S",new byte[][]{consensus},null,new String[]{MODEL},seeds,17,20,50,7,1,false,1f,0f,0f,1,0f,1f,.6f,.6f);
		f.reuseConsensusAlignment=true;f.trimAlignmentExtent=false;f.scavengePass2=false;return f;
	}
	static NcrnaScavenger scavenger(NcrnaFamily f){
		final NcrnaScavenger s=new NcrnaScavenger(f.library,f.models,f.modelNames,f.kmerSet,f.kLong,f.minLen,f.windowPad,f.indexK,f.indexTopN,f.adaptive,f.adaptFloor,f.adaptTopFrac,f.adaptQFrac,f.fixedMinHits,f.scoreA,f.scoreB,f.idPass,f.idBorderline);GeneCaller.applyNcrnaFamilyControls(s,f);return s;
	}
	static Orf candidate(byte[] bases,int start,int stop,int strand){final Orf o=new Orf("contig",start,stop,strand,0,bases,false,ProkObject.RNA);o.ncrnaFamily="euk5S";o.trnaModel=MODEL;return o;}
	static String calls(ArrayList<Orf> list){final StringBuilder b=new StringBuilder();for(Orf o:list){b.append(o.start).append(':').append(o.stop).append(':').append(o.strand).append(':').append(o.orfScore).append(':').append(o.trnaModel).append('\n');}return b.toString();}
	static byte[] features(String path){final fileIO.ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');byte[] found=null;
		try{for(byte[] row;(row=in.nextLine())!=null;){if(row[0]=='#'){continue;}check(found==null,"Fixture has one materialized row per end");p.set(row);check(p.terms()==53,"Materializer schema retains28inputs/25labels");p.setBounds(27);found=Arrays.copyOf(row,p.b());}}finally{check(!in.close(),"Fixture read must finish");}check(found!=null,"Missing materialized feature row");return found;}
	static final class Recorder implements RrnaEndpointCallerFeatures.Sink{
		@Override public void capture(String contig,int strand_,int model,int start,int stop,float[] f,boolean fu,float[] t,boolean tu){count++;strand=strand_;five=serialize(f);three=serialize(t);fiveUsable=fu;threeUsable=tu;}
		static byte[] serialize(float[] values){final ByteBuilder b=new ByteBuilder();for(int i=0;i<values.length;i++){if(i>0){b.tab();}b.appendSlow(values[i]);}return b.toBytes();}
		long count;int strand;byte[] five,three;boolean fiveUsable,threeUsable;
	}
	static byte[] bytes(String s){return s.getBytes(StandardCharsets.US_ASCII);}
	static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
	static final String MODEL="euk5S_fixture";
}
