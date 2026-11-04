package test.prok;

import prok.*;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import idaligner.ScrabbleAligner;

/** Proves the G28 v3 width-10 layout and zero-alignment per-site path. */
public final class NcrnaBoundaryFeatureV3Test {

	public static void main(String[] args) throws Exception{
		final float[] direct=NcrnaBoundaryScorer.buildV3Input(new float[]{.1f,.2f,.3f},true,.4f,.5f,.6f);
		check(Arrays.equals(direct,new float[]{.1f,.2f,.3f,1f,.4f,.5f,.6f,0f,0f,0f}),"v3 direct feature order changed: "+Arrays.toString(direct));
		final float[] consensusProfile={.2f,.4f,.6f};NcrnaBoundaryVectorGen.blendProfiles(consensusProfile,new float[]{.6f,.4f,.2f});
		check(Arrays.equals(consensusProfile,new float[]{.4f,.4f,.4f}),"50/50 family/consensus profile blend changed: "+Arrays.toString(consensusProfile));
		final byte[] window=new byte[100];Arrays.fill(window,(byte)'A');
		final NcrnaBoundaryVectorGen.VectorStats stats=new NcrnaBoundaryVectorGen.VectorStats();
		stats.records=1;stats.processedLoci=1;
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		final ScrabbleAligner counter=new ScrabbleAligner();final long loops0=counter.loops();
		try(PrintStream out=new PrintStream(bytes,false,"UTF-8")){
			final int emitted=NcrnaBoundaryVectorGen.emitVectors(out,window,20,80,true,null,null,0,0,null,.4f,61f,
				new int[]{0,1},false,stats,true,NcrnaBoundaryScorer.FEATURES_V3,.625f,null);
			check(emitted==2,"v3 emitted count="+emitted);
		}
		check(counter.loops()==loops0,"v3 candidate scoring invoked ScrabbleAligner");
		check(!NcrnaBoundaryScorer.usesPerSiteScrabble(NcrnaBoundaryScorer.FEATURES_V3),
			"G28 V3 endpoint inference must structurally bypass per-site Scrabble alignment");
		check(NcrnaBoundaryScorer.usesPerSiteScrabble(NcrnaBoundaryScorer.FEATURES_V1),
			"G16 V1 alternate unexpectedly lost its explicit per-site alignment classification");
		final String[] lines=new String(bytes.toByteArray(),StandardCharsets.UTF_8).trim().split("\\R");
		check(lines.length==2,"v3 line count="+lines.length);
		final String[] a=lines[0].split("\\t"),b=lines[1].split("\\t");
		check(a.length==11&&b.length==11,"v3 vector width");
		check(a[6].equals("0.625000")&&b[6].equals("0.625000"),"Quantum identity must occupy feature 7 and remain locus-constant");
		for(int i=7;i<=9;i++){check(a[i].equals("0.000000")&&b[i].equals("0.000000"),"v3 spare column must be zero: "+i);}
		check(!a[4].equals(b[4]),"length ratio must remain candidate-specific at feature 5");
		testModelColumn(window);
		testMismatchedKFailsLoudly();
		System.out.println("PASS NcrnaBoundaryFeatureV3Test");
	}

	private static void testMismatchedKFailsLoudly(){
		final TrnaBoundaryFeatures.NinemerTable zero=(seq,pos)->0;
		final TrnaNinemerTableBuilder.LoadedTable consensus=new TrnaNinemerTableBuilder.LoadedTable(9,
			TrnaBoundaryFeatures.BoundaryType.START,4,5,zero);
		final TrnaNinemerTableBuilder.LoadedTable family=new TrnaNinemerTableBuilder.LoadedTable(7,
			TrnaBoundaryFeatures.BoundaryType.START,4,5,zero);
		boolean failed=false;
		try{
			NcrnaBoundaryVectorGen.validateFamilyConsensusTable(consensus,family,
				TrnaBoundaryFeatures.BoundaryType.START,"start");
		}catch(IllegalArgumentException e){
			failed=e.getMessage()!=null&&e.getMessage().contains("k mismatch");
		}
		check(failed,"mismatched family/consensus k did not fail loudly");
	}

	private static void testModelColumn(byte[] window) throws Exception{
		final NcrnaBoundaryVectorGen.VectorStats stats=new NcrnaBoundaryVectorGen.VectorStats();stats.records=1;stats.processedLoci=1;
		final ByteArrayOutputStream taggedBytes=new ByteArrayOutputStream(),trainerBytes=new ByteArrayOutputStream();
		try(PrintStream tagged=new PrintStream(taggedBytes,false,"UTF-8");PrintStream trainer=new PrintStream(trainerBytes,false,"UTF-8")){
			NcrnaBoundaryVectorGen.emitVectors(tagged,window,20,80,false,null,null,0,0,null,.4f,61f,new int[]{0},false,stats,true,
				NcrnaBoundaryScorer.FEATURES_V3,.625f,null,null,trainer,"archaea",true);
		}
		final String[] tagged=new String(taggedBytes.toByteArray(),StandardCharsets.UTF_8).trim().split("\\t");
		final String[] trainer=new String(trainerBytes.toByteArray(),StandardCharsets.UTF_8).trim().split("\\t");
		check(tagged.length==12&&tagged[10].equals("archaea")&&tagged[11].equals("1"),"model-tagged vector schema changed");
		check(trainer.length==11&&trainer[10].equals("1"),"numeric trainer vector schema changed");
	}

	private static void check(boolean condition,String message){if(!condition){throw new AssertionError(message);}}
}
