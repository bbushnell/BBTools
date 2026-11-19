package test.prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

import ml.CellNet;
import ml.CellNetParser;

/** Exercises stored-label shapes and every exported output against analytic predictions.
 * @author Yoimiya
 */
public final class SubnetGpuRoundTripFixtureTest {
	/** Uses a linear six-head network so expected values are independent of inference code. */
	public static void main(String[] args) throws Exception{
		if(args.length!=0){throw new IllegalArgumentException("Fixture accepts no arguments");}
		final Path root=Files.createTempDirectory("subnet-roundtrip-");
		final Path vectors=root.resolve("vectors.tsv"),predictions=root.resolve("predictions.tsv");
		try{
			final ArrayList<byte[]> lines=new ArrayList<byte[]>();
			final String header="##bbnet\n#version 1\n#concise\n#dense\n#seed 1\n#layers 2\n#dims 2 6\n#coding decimal\n";
			for(String line : header.split("\n")){lines.add(line.getBytes(StandardCharsets.UTF_8));}
			for(int i=0; i<6; i++){
				lines.add(("C"+(i+3)+" LINEAR "+i+" 1 2").getBytes(StandardCharsets.UTF_8));
			}
			final CellNet net=CellNetParser.loadFromLines(lines);
			Files.write(predictions,"0\t5\t6\t7\t8\t9\t10\n".getBytes(StandardCharsets.UTF_8));
			Files.write(vectors,"#dims\t2\t2\t0\n1\t2\t0.7\t0.1\n".getBytes(StandardCharsets.UTF_8));
			final double error=SubnetGpuRoundTripTest.compare(net,vectors.toString(),predictions.toString(),1,0.00002,null);
			if(error!=0){throw new AssertionError("Analytic linear predictions must match all six heads exactly");}
			final CellNet reference=CellNetParser.loadFromLines(lines);
			if(SubnetGpuRoundTripTest.compare(net,vectors.toString(),null,1,0,reference)!=0){
				throw new AssertionError("Network-to-network mode must compare all six outputs exactly without a prediction file");
			}
			final ArrayList<byte[]> changed=new ArrayList<byte[]>(lines);
			changed.set(changed.size()-1,"C8 LINEAR 6 1 2".getBytes(StandardCharsets.UTF_8));
			boolean changedRejected=false;
			try{SubnetGpuRoundTripTest.compare(net,vectors.toString(),null,1,0,CellNetParser.loadFromLines(changed));}
			catch(IllegalStateException expected){changedRejected=expected.getMessage().contains("output=5");}
			if(!changedRejected){throw new AssertionError("Reference-network mode must reject a changed sixth-head bias");}
			final ArrayList<byte[]> sparseLines=new ArrayList<byte[]>();
			for(String line : header.replace("#dense", "#sparse").split("\n")){
				sparseLines.add(line.getBytes(StandardCharsets.UTF_8));
			}
			for(int i=0; i<6; i++){
				sparseLines.add(("C"+(i+3)+" LINEAR "+i+" 1 2").getBytes(StandardCharsets.UTF_8));
				sparseLines.add(("I"+(i+3)+" 0 1").getBytes(StandardCharsets.UTF_8));
			}
			final CellNet localDense=CellNetParser.loadInferenceFromLines(lines);
			final CellNet localSparse=CellNetParser.loadInferenceFromLines(sparseLines);
			for(boolean unrelatedDense : new boolean[]{true, false}){
				final CellNet unrelated=CellNetParser.loadFromLines(unrelatedDense ? lines : sparseLines);
				if(unrelated.dense!=unrelatedDense || !localDense.dense || localSparse.dense){
					throw new AssertionError("Each loaded network must retain its own dense/sparse representation");
				}
				compareMixed(localDense, localSparse, vectors);
				compareMixed(localSparse, localDense, vectors);
			}
			Files.write(predictions,"0\t5\t6\t7\t8\t9\t11\n".getBytes(StandardCharsets.UTF_8));
			boolean rejected=false;
			try{SubnetGpuRoundTripTest.compare(net,vectors.toString(),predictions.toString(),1,0.00002,null);}
			catch(IllegalStateException expected){rejected=expected.getMessage().contains("output=5");}
			if(!rejected){throw new AssertionError("A mismatch in the sixth head must fail comparison");}
			Files.write(vectors,"#dims\t2\t3\t0\n1\t2\t0.7\t0.1\t0\n".getBytes(StandardCharsets.UTF_8));
			rejected=false;
			try{SubnetGpuRoundTripTest.compare(net,vectors.toString(),predictions.toString(),1,0.00002,null);}
			catch(IllegalArgumentException expected){rejected=expected.getMessage().startsWith("Stored targets");}
			if(!rejected){throw new AssertionError("Unsupported three-label/six-head vectors must fail");}
			System.out.println("SIX_OUTPUT_ROUNDTRIP_FIXTURE_PASS");
		}finally{
			Files.deleteIfExists(vectors);
			Files.deleteIfExists(predictions);
			Files.delete(root);
		}
	}

	/** Reference comparison must dispatch from each network's own representation. */
	private static void compareMixed(CellNet net, CellNet reference, Path vectors){
		if(SubnetGpuRoundTripTest.compare(net, vectors.toString(), null, 1, 0, reference)!=0){
			throw new AssertionError("Dense/sparse native comparison must preserve all six output bits after unrelated model loads");
		}
	}
}
