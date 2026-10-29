package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import dna.AminoAcid;
import map.LongHashSet;
import prok.CallGenes;
import prok.GeneCaller;
import prok.NcrnaFamily;
import prok.NcrnaScavenger;
import prok.Orf;
import structures.ByteBuilder;

/** Real HBM and identity-rescue paths, including a later shortlisted model.
 * @author Ganyu */
public final class NcrnaTmrnaRescueMarginTest {
	private static final byte[] MODEL=bytes("ACGTCAGTGCATCGACTGACGATCGTAGCTAGCATGCTACGATCGATGCTAGCTAGCTACGTAGC");
	private static final byte[] QUERY=bytes("TTGAGCGATCGTACGTTAGCATCGGATGCTACGTTCAGCGATGCTACGATCGTAGCATGCGATCA");
	// Deterministic synthetic fixture: primary max=.76209676; rescued identities
	// are .7815126 then .79565215, in K1 shortlist order. No HBM can mask this path.
	private static final byte[] WINDOW=bytes("TGTACCGCCAACCACGACTTCACGGTACTTCAACACAAAGGAAATCCCGTCTTGGGCCTCACCCCGGGGGGGTTGGGGGACTTCACAACTCATAATTGCCTTAAGTTGTGCCGTAATGTATAAGCAAGATTCCTCGTAGCCATATTTCCGCAGAACTATCCTACCATGCTGTGAGATCCCTTTGTGCAATGATGAAGCCACGGCGGGCCCGTCCTCCGATCCCCATGTGTCTTGCCAGCCCACCCGGCAGCGACACCTCTAACGTTAACCATATGAAATCCTACGAACTAGGGCAGTATACACTCTGAGCGACGCCTACAGTGAATTGCGCTGGTCTGTTTCCCTCAATTCGTCGGAGAT");
	private static final byte[] FIRST=bytes("AGCTGCGATCCCGAGGTTGGCCGGCGCAACGTGCAGAGGTCTTCAAAGGCCGTATCCTGCTTAAGTTGTGCCGTAATGTATATTCAAGATTCCTCGTAGCCCTATTTCGGCAGCACTATCCTACCATGCTGTGGGATCGCTTTGTCCAATGATGAAGCCACGGCGGGCCCGTACTTCGATCCCCATGTGTCTTGTCATCCCACCTCGCAGCGACGCCTCTCTTGGGAGCGATTATTCACTTGTGTTTATGCTATCTTTTAAGCGCTGTACGCTACTAGGC");
	private static final byte[] LATER=bytes("TCAGTTATCTGCATGAAAGGTGTTCTACGTACGAACGCTGTGAAGCTGTCCCGTAATGTCTAAGCAAGATCCCTCGTCCCCATATTTCGGCAGAACTATCCTACCATGCTGTGAGATCCATGTGTGCAATGATGAAGCCACGGCGGGCCCGTTCTCCAATCCCCATGTGTCTTGCCAGCCCACCCGGCAGCGACACCTGTATAGCTGCATAACGTGAAGGGTAATCACCGTGCCTTATTA");

	public static void main(String[] args)throws Exception{
		check(CallGenes.TMRNA_RESCUE_ID_MARGIN==.005f && CallGenes.TMRNA_RESCUE_HBM_MARGIN==.005f,
			"The accepted production curve selected 0.005 for both rescued-acceptance margins");
		final boolean oldEnabled=CallGenes.NCRNA_FAMILIES_ENABLED,oldTmrna=CallGenes.TMRNA_ENABLED;
		try{
			CallGenes.NCRNA_FAMILIES_ENABLED=false;CallGenes.TMRNA_ENABLED=false;
			CallGenes.validateNcrnaSweepOverrides();
			CallGenes.TMRNA_RESCUE_ID_MARGIN_EXPLICIT=true;
			CallGenes.TMRNA_RESCUE_ID_MARGIN=0;
			CallGenes.validateNcrnaSweepOverrides();
			CallGenes.TMRNA_RESCUE_ID_MARGIN_EXPLICIT=false;
			CallGenes.TMRNA_RESCUE_ID_MARGIN=.005f;
			CallGenes.TMRNA_RESCUE_HBM_MARGIN_EXPLICIT=true;
			CallGenes.TMRNA_RESCUE_HBM_MARGIN=0;
			CallGenes.validateNcrnaSweepOverrides();
			CallGenes.TMRNA_RESCUE_HBM_MARGIN=.005f;
			boolean rejected=false;
			try{CallGenes.validateNcrnaSweepOverrides();}catch(IllegalArgumentException e){rejected=e.getMessage().contains("tmRNA-specific overrides");}
			check(rejected,"Default margins must allow disabled tmRNA, while explicit nonzero overrides require it enabled");
		}finally{
			CallGenes.NCRNA_FAMILIES_ENABLED=oldEnabled;CallGenes.TMRNA_ENABLED=oldTmrna;
			CallGenes.TMRNA_RESCUE_ID_MARGIN_EXPLICIT=false;CallGenes.TMRNA_RESCUE_HBM_MARGIN_EXPLICIT=false;
			CallGenes.TMRNA_RESCUE_ID_MARGIN=.005f;CallGenes.TMRNA_RESCUE_HBM_MARGIN=.005f;
		}
		testHbmAndPrimary();testIdentityContinuation();testInvalidMargins();testInvalidCli();
		System.out.println("PASS zero-margin byte/work parity, HBM OR fallback, identity rescue continuation/rejection, unchanged primary/other-family behavior, production plumbing and early CLI validation");
	}
	private static void testHbmAndPrimary(){
		final NcrnaScavenger baseline=scavenger("tmrna",QUERY);
		final ArrayList<Orf> base=run(baseline,QUERY);
		check(base.size()==1 && baseline.hbmScoreCalls()==1,"Fixture must reach and pass HBM rescue, not primary identity");
		final NcrnaScavenger zero=scavenger("tmrna",QUERY);zero.setTmrnaRescueMargins(0,0);
		check(gff(base).equals(gff(run(zero,QUERY))) && zero.alignmentCount()==baseline.alignmentCount()
			&& zero.hbmScoreCalls()==baseline.hbmScoreCalls(),"Explicit zero margins preserve emitted bytes and verifier work");
		final NcrnaScavenger idOnly=scavenger("tmrna",QUERY);idOnly.setTmrnaRescueMargins(.02f,0);
		check(gff(base).equals(gff(run(idOnly,QUERY))),"Raising rescued identity must preserve independently passing HBM fallback");
		final NcrnaScavenger raised=scavenger("tmrna",QUERY);raised.setTmrnaRescueMargins(0,.02f);
		check(run(raised,QUERY).isEmpty() && raised.hbmScoreCalls()==1,"Raised HBM margin must reject borderline HBM rescue");
		final NcrnaScavenger exactBase=scavenger("tmrna",MODEL);
		final String exact=gff(run(exactBase,MODEL));
		final NcrnaScavenger primary=scavenger("tmrna",MODEL);primary.setTmrnaRescueMargins(.2f,.2f);
		check(!exact.isEmpty() && exact.equals(gff(run(primary,MODEL))) && primary.hbmScoreCalls()==0
			&& primary.alignmentCount()==exactBase.alignmentCount(),"Primary acceptance bytes and work must ignore rescue margins");
		// Exercise the actual GeneCaller field-to-instance path for both family scopes.
		final float oldId=CallGenes.TMRNA_RESCUE_ID_MARGIN,oldHbm=CallGenes.TMRNA_RESCUE_HBM_MARGIN;
		try{
			CallGenes.TMRNA_RESCUE_ID_MARGIN=.02f;CallGenes.TMRNA_RESCUE_HBM_MARGIN=.02f;
			final NcrnaScavenger routed=scavenger("tmrna",QUERY);applyControls(routed,"tmrna");
			check(run(routed,QUERY).isEmpty(),"Production family controls must deliver the HBM margin to the verifier");
			final NcrnaScavenger other=scavenger("rnasep",QUERY);applyControls(other,"rnasep");
			check(gff(base).equals(gff(run(other,QUERY)).replace("family:RNase P","family:tmRNA")),"Global tmRNA margins must not affect another family");
		}finally{CallGenes.TMRNA_RESCUE_ID_MARGIN=oldId;CallGenes.TMRNA_RESCUE_HBM_MARGIN=oldHbm;}
	}
	private static void testIdentityContinuation(){
		final NcrnaScavenger baseline=identityScavenger();
		final ArrayList<Orf> first=run(baseline,WINDOW);
		check(first.size()==1 && "first".equals(first.get(0).trnaModel) && baseline.alignmentCount()==3,
			"Two failed primary alignments must lead to first-model identity rescue with early exit");
		final NcrnaScavenger zero=identityScavenger();zero.setTmrnaRescueMargins(0,0);
		check(gff(first).equals(gff(run(zero,WINDOW))) && zero.alignmentCount()==3,"Explicit zero must preserve identity-rescue early exit");
		final NcrnaScavenger raised=identityScavenger();raised.setTmrnaRescueMargins(.0215f,0);
		final ArrayList<Orf> later=run(raised,WINDOW);
		check(later.size()==1 && "later".equals(later.get(0).trnaModel) && raised.alignmentCount()==4
			&& raised.hbmScoreCalls()==0,"Raised identity cutoff must continue past the first model and accept the later model");
		final NcrnaScavenger reject=identityScavenger();reject.setTmrnaRescueMargins(.05f,0);
		check(run(reject,WINDOW).isEmpty() && reject.alignmentCount()==4,"Both rescues below the raised identity cutoff must be rejected");
	}
	private static NcrnaScavenger identityScavenger(){
		final NcrnaScavenger s=new NcrnaScavenger(new byte[][]{FIRST,LATER},null,new String[]{"first","later"},
			seeds(WINDOW),17,40,400,1,2,false,0,0,0,0,0,1,.76709676f,.70f);
		s.family="tmrna";s.scavengePass2=false;s.trimAlignmentExtent=false;return s;
	}
	private static void testInvalidMargins(){
		for(float bad:new float[]{-.01f,Float.NaN,Float.POSITIVE_INFINITY,.9f}){
			boolean rejected=false;
			try{scavenger("tmrna",MODEL).setTmrnaRescueMargins(bad,0);}catch(IllegalArgumentException e){rejected=true;}
			check(rejected,"Invalid identity margin must fail before an otherwise passing primary: "+bad);
			rejected=false;
			try{scavenger("tmrna",MODEL).setTmrnaRescueMargins(0,bad);}catch(IllegalArgumentException e){rejected=true;}
			check(rejected,"Invalid HBM margin must fail before an otherwise passing primary: "+bad);
		}
	}
	private static void testInvalidCli()throws Exception{
		final Path dir=Files.createTempDirectory("tmrna-rescue-cli-");
		final Path in=dir.resolve("input.fa"),out=dir.resolve("output.gff");
		Files.write(in,bytes(">exact\n"+new String(MODEL,StandardCharsets.US_ASCII)+"\n"));
		final float oldId=CallGenes.TMRNA_RESCUE_ID_MARGIN,oldHbm=CallGenes.TMRNA_RESCUE_HBM_MARGIN;
		try{
			for(boolean identity:new boolean[]{true,false}){
				CallGenes.TMRNA_RESCUE_ID_MARGIN=0;CallGenes.TMRNA_RESCUE_HBM_MARGIN=0;
				boolean rejected=false;
				try{new CallGenes(new String[]{"in="+in,"out="+out,"ncrna=t","tmrna=t","taxonomy=f","t=1",
					identity?"tmrnarescueidmargin=.9":"tmrnarescuehbmmargin=.9"});}
				catch(IllegalArgumentException e){rejected=e.getMessage().contains("tmRNA effective rescue");}
				check(rejected && !Files.exists(out),"Effective cutoff >1 must fail during CLI setup before scanning/output");
				check((identity?CallGenes.TMRNA_RESCUE_ID_MARGIN:CallGenes.TMRNA_RESCUE_HBM_MARGIN)==.9f,
					"Normal CLI spelling must reach its owned margin field");
			}
		}finally{
			CallGenes.TMRNA_RESCUE_ID_MARGIN=oldId;CallGenes.TMRNA_RESCUE_HBM_MARGIN=oldHbm;
			Files.deleteIfExists(out);Files.delete(in);Files.delete(dir);
		}
	}
	private static void applyControls(NcrnaScavenger s,String name){
		final NcrnaFamily f=new NcrnaFamily(name,new byte[][]{MODEL},null,null,seeds(QUERY),17,60,100,
			7,1,false,0,0,0,0,0,1,.65f,.60f,.55f,.9f,null,null,null,null,-1,-1,-1,-1,0);
		f.scavengePass2=false;GeneCaller.applyNcrnaFamilyControls(s,f);
	}
	private static NcrnaScavenger scavenger(String family,byte[] query){
		final BaseGraph graph=new BaseGraph("fixture",MODEL,null,0,0);
		for(int i=0;i<MODEL.length;i++){
			BaseGraphHelper.setCounts(graph,i,MODEL[i]=='A'?10:1,MODEL[i]=='C'?10:1,MODEL[i]=='G'?10:1,MODEL[i]=='T'?10:1,0);
		}
		BaseGraphHelper.initForScoring(graph);
		final NcrnaScavenger s=new NcrnaScavenger(new byte[][]{MODEL},new BaseGraph[]{graph},null,seeds(query),17,60,100,
			7,1,false,0f,0f,0f,0,0f,1f,.65f,.60f);
		s.family=family;s.hbmPass=.55f;s.scavengePass2=false;return s;
	}
	private static LongHashSet seeds(byte[] query){
		assert(query.length>=17):"Fixture requires a complete seed inside its supplied query";
		long key=0;
		for(int i=0;i<17;i++){final int b=AminoAcid.baseToNumber[query[i]];check(b>=0,"Fixture seed must be unambiguous");key=(key<<2)|b;}
		final LongHashSet seeds=new LongHashSet(2);seeds.add(key);return seeds;
	}
	private static ArrayList<Orf> run(NcrnaScavenger s,byte[] query){return s.scavenge("fixture",query.clone(),0,new ArrayList<int[]>());}
	private static String gff(ArrayList<Orf> calls){final ByteBuilder b=new ByteBuilder();for(Orf o:calls){o.appendGff(b);b.nl();}return b.toString();}
	private static byte[] bytes(String s){return s.getBytes(StandardCharsets.US_ASCII);}
	private static void check(boolean ok,String message){if(!ok){throw new AssertionError(message);}}
}
