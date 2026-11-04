package test.prok;

import prok.*;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import idaligner.ScrabbleAligner;

/** Proves v2 vectors reuse locus/model constants and emit no candidate alignment dependency. */
public final class NcrnaBoundaryFeatureV2Test {

	public static void main(String[] args) throws Exception{
		final byte[] window=new byte[100];Arrays.fill(window,(byte)'A');
		final float[] fuzz={.25f,.5f,.75f};
		final NcrnaBoundaryVectorGen.VectorStats stats=new NcrnaBoundaryVectorGen.VectorStats();
		stats.records=1;stats.processedLoci=1;
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		final ScrabbleAligner counter=new ScrabbleAligner();
		final long loops0=counter.loops();
		try(PrintStream out=new PrintStream(bytes,false,"UTF-8")){
			final int emitted=NcrnaBoundaryVectorGen.emitVectors(out,window,20,80,true,null,null,0,0,null,.4f,61f,
				new int[]{0,1},false,stats,true,NcrnaBoundaryScorer.FEATURES_V2,.625f,fuzz);
			check(emitted==2,"v2 emitted count="+emitted);
		}
		check(counter.loops()==loops0,"v2 candidate scoring invoked ScrabbleAligner");
		final String[] lines=new String(bytes.toByteArray(),StandardCharsets.UTF_8).trim().split("\\R");
		check(lines.length==2,"v2 line count="+lines.length);
		final String[] a=lines[0].split("\\t"),b=lines[1].split("\\t");
		check(a.length==11 && b.length==11,"v2 vector width");
		check(a[0].equals("0.625000") && b[0].equals("0.625000"),"ANI must be locus-constant");
		for(int i=5; i<=7; i++){check(a[i].equals(b[i]),"HBM fuzz must be model/end-constant at column "+i);}
		check(!a[8].equals(b[8]),"length ratio must remain candidate-specific");
		System.out.println("PASS NcrnaBoundaryFeatureV2Test");
	}

	private static void check(boolean condition,String message){if(!condition){throw new AssertionError(message);}}
}
