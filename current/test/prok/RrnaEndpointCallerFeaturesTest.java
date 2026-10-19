package test.prok;

import prok.*;

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
		fixedEnds(bases,consensus,five,three,gc);
		System.out.println("PASS RrnaEndpointCallerFeaturesTest actual_family_hook byte_identical_28cols mixed_k6_k9 nonzero_anchors strands lowercase clipping default_off calls_unchanged alignment_accounting");
	}
	static void fixedEnds(byte[] bases,byte[] consensus,RrnaPositionalKmerTable.Table five,RrnaPositionalKmerTable.Table three,float gc){
		final ml.CellNet[] fiveNets={constantNet(five.inputs,five.sites,9)},threeNets={constantNet(three.inputs,three.sites,16)};
		final RrnaEndpointCallerFeatures.Resources resources=new RrnaEndpointCallerFeatures.Resources("prok5S",new String[]{MODEL},new byte[][]{consensus},
			new RrnaPositionalKmerTable.Table[]{five},new RrnaPositionalKmerTable.Table[]{three},fiveNets,threeNets,null,true);
		final Recorder oldRecorder=new Recorder(),zeroRecorder=new Recorder();
		final RrnaEndpointCallerFeatures old=new RrnaEndpointCallerFeatures(resources,oldRecorder),zero=new RrnaEndpointCallerFeatures(resources,zeroRecorder);
		final Orf a=prokCandidate(bases,53,125,0),b=prokCandidate(bases,53,125,0);a.orfScore=b.orfScore=123f;
		old.capture(a,bases,0,gc,70,200,.73f);zero.capture(b,bases,0,gc,70,200,.73f,0);
		check(a.start==50 && a.stop==129 && b.start==a.start && b.stop==a.stop
			&& Arrays.equals(oldRecorder.five,zeroRecorder.five) && Arrays.equals(oldRecorder.three,zeroRecorder.three),
			"Legacy overload and explicit mask0 must retain exact feature bytes and real inference results");
		for(int strand=0;strand<2;strand++){
			for(int mask=1;mask<=2;mask++){
				final Recorder recorder=new Recorder();final RrnaEndpointCallerFeatures inference=new RrnaEndpointCallerFeatures(resources,recorder);
				// A disabled network must never be dereferenced, even when its table would score.
				if(mask==1){inference.fiveNets[0]=null;}else{inference.threeNets[0]=null;}
				final Orf candidate=prokCandidate(bases,53,125,strand);candidate.orfScore=123f;
				inference.capture(candidate,bases,0,gc,70,200,.73f,mask);
				check(candidate.start==(mask==1?53:50) && candidate.stop==(mask==2?125:129) && candidate.orfScore==123f,
					"Only unfixed oriented endpoint moves; accepted detection score remains unchanged on both strands");
				check(recorder.fiveUsable==(mask==2) && recorder.threeUsable==(mask==1),"Observation cannot claim a fixed endpoint was scored");
				final float[] active=mask==1?inference.three:inference.five,disabled=mask==1?inference.five:inference.three;
				check(active[resources.sites]==73f/consensus.length && active[resources.sites+2]==.73f,
					"Interior endpoint retains full trained consensus length and supplied accepted identity, not cropped-model substitutes");
				for(float value:disabled){check(Float.isNaN(value),"Disabled end must not leak preceding-candidate feature values");}
			}
		}
		final RrnaEndpointCallerFeatures scoped=new RrnaEndpointCallerFeatures(resources,new Recorder());
		final Orf tooShort=prokCandidate(bases,53,125,0);scoped.capture(tooShort,bases,0,gc,80,200,.73f,1);
		check(tooShort.start==53 && tooShort.stop==125,"Ordinary full-model minimum rejects77bp proposal, while explicit partial70 floor above permits it");
		// Nulling the disabled table proves feature extraction is skipped, not just its output discarded.
		resources.five[0]=null;
		final Orf edge=prokCandidate(bases,0,125,0);scoped.capture(edge,bases,0,gc,70,200,.73f,1);
		check(edge.start==0 && edge.stop==129,"Physical5prime coordinate stays fixed while normal3prime inference runs");
		resources.five[0]=five;
		final Recorder silent=new Recorder();final RrnaEndpointCallerFeatures both=new RrnaEndpointCallerFeatures(resources,silent);
		final Orf fixed=prokCandidate(bases,0,179,1);both.capture(fixed,null,0,Float.NaN,70,200,Float.NaN,3);
		check(silent.count==0 && fixed.start==0 && fixed.stop==179,"Both fixed means no features, identity computation, NN inference or sink observation");
		boolean rejected=false;try{both.capture(fixed,bases,0,gc,70,200,.73f,4);}catch(IllegalArgumentException e){rejected=true;}
		check(rejected,"Unknown mask bits cannot silently suppress endpoint work");
		final String rawName="mitochondria";
		final RrnaEndpointCallerFeatures.Resources mitoResources=new RrnaEndpointCallerFeatures.Resources("16S",new String[]{MODEL,rawName},new byte[][]{consensus,consensus},
			new RrnaPositionalKmerTable.Table[]{five,null},new RrnaPositionalKmerTable.Table[]{three,null},new ml.CellNet[]{fiveNets[0],null},new ml.CellNet[]{threeNets[0],null},new boolean[]{false,true},true);
		final Recorder mitoRecorder=new Recorder();final RrnaEndpointCallerFeatures mito=new RrnaEndpointCallerFeatures(mitoResources,mitoRecorder);
		for(int mask=0;mask<=3;mask++){final Orf candidate=candidate(bases,53,125,0);candidate.ncrnaFamily="16S";candidate.trnaModel=rawName;
			mito.capture(candidate,bases,1,gc,70,200,.73f,mask);check(candidate.start==53 && candidate.stop==125,"Raw mitochondrial model must retain both endpoints for every mask");}
		check(mitoRecorder.count==0,"Raw mitochondrial exemption must not produce endpoint feature rows");
		System.out.println("PASS fixed-end masks, mask0 byte parity, both strands, partial minimum, full-reference/accepted identity, no disabled work and raw mitochondrial exemption");
	}
	static ml.CellNet constantNet(int inputs,int sites,int winner){
		final ml.CellNet net=new ml.CellNet(new int[]{inputs,sites},931,1f,0f,1,null);
		// CellNet.randomize requires normalized activation-selection rates even
		// without hidden layers. Preserve shared settings around this test setup.
		final float[] rates=ml.Function.TYPE_RATES.clone(),cumulative=ml.Function.TYPE_RATES_CUM;
		try{
			if(cumulative==null){ml.Function.normalizeTypeRates();}
			net.randomize();
		}finally{
			System.arraycopy(rates,0,ml.Function.TYPE_RATES,0,rates.length);
			ml.Function.TYPE_RATES_CUM=cumulative;
		}
		for(int i=0;i<net.finalLayer.length;i++){final ml.Cell cell=net.finalLayer[i];Arrays.fill(cell.weights,0f);cell.setBias(i==winner?1f:0f,true);cell.function=ml.Function.getFunction(ml.Function.LINEAR);}
		return net;
	}
	static Orf prokCandidate(byte[] bases,int start,int stop,int strand){final Orf o=candidate(bases,start,stop,strand);o.ncrnaFamily="prok5S";return o;}
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
