package test.ml;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

import ml.CellNet;
import ml.CellNetParser;
import parse.PreParser;
import shared.Shared;

/** Exercises stored-label shapes and every exported output against analytic predictions.
 * @author Yoimiya, Raiden
 */
public final class SubnetGpuRoundTripFixtureTest {
	/** Uses a linear six-head network so expected values are independent of inference code. */
	public static void main(String[] args) throws Exception{
		// This no-argument fixture must run rather than take the normal tool help exit.
		final PreParser pre=new PreParser(args, System.err, SubnetGpuRoundTripFixtureTest.class, false, true, false);
		try{
			if(pre.args.length!=0){throw new IllegalArgumentException("Fixture accepts no arguments");}
			run();
		}finally{Shared.closeStream(pre.outstream);}
	}

	/** Runs the retained small analytic inputs after standard argument handling. */
	private static void run() throws Exception{
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
			checkExportCopyContract(net, vectors, predictions);
			final Path image=root.resolve("ordinary.bbnet");
			Files.write(image, net.toBytes().toBytes());
			try{
				SubnetGpuRoundTripTest.main(new String[]{"net="+image, "vectors="+vectors,
					"predictions="+predictions, "rows=1", "ordinary=t", "argmax=t"});
				checkPredictionBits(net,image,vectors,predictions,root);
			}finally{Files.deleteIfExists(image);}
			final ArrayList<byte[]> tiedLines=new ArrayList<byte[]>(lines);
			tiedLines.set(tiedLines.size()-2, "C7 LINEAR 5 1 2".getBytes(StandardCharsets.UTF_8));
			Files.write(predictions, "0\t5\t6\t7\t8\t9.999999\t10\n".getBytes(StandardCharsets.UTF_8));
			final long[] argmaxChanges=new long[1];
			final double tiedError=SubnetGpuRoundTripTest.compare(CellNetParser.loadFromLines(tiedLines),
				vectors.toString(), predictions.toString(), 1, 0.00002, null, argmaxChanges);
			if(tiedError<=0 || argmaxChanges[0]!=1){
				throw new AssertionError("A sub-tolerance difference can change first argmax and must be reported separately");
			}
			final CellNet tied=CellNetParser.loadFromLines(tiedLines),near=tied.copy(false);
			near.net[1][4].bias=Math.nextDown(Math.nextDown(near.net[1][4].bias));
			final SubnetGpuRoundTripTest.Comparison selection=SubnetGpuRoundTripTest.compareExportPaths(
				tied,vectors.toString(),predictions.toString(),1,.00002,near);
			if(selection.nativeArgmaxChanges!=1 || selection.oracleArgmaxChanges!=1 ||
					selection.referenceOracleArgmaxChanges!=0 || selection.nativeMaximum<=0){
				throw new AssertionError("Every native and oracle first-maximum difference must be measured separately");
			}
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

	/** Captures actual checked bits, preserves failed prefixes, and never replaces output. */
	private static void checkPredictionBits(CellNet net,Path image,Path vectors,Path predictions,Path root) throws Exception{
		final Path accepted=root.resolve("checked.bits.tsv"),failed=root.resolve("failed.bits.tsv"),
				late=root.resolve("late.bits.tsv");
		final Path acceptedPartial=root.resolve("checked.bits.tsv.partial"),
				failedPartial=root.resolve("failed.bits.tsv.partial"),latePartial=root.resolve("late.bits.tsv.partial");
		final byte[] originalVectors=Files.readAllBytes(vectors),originalPredictions=Files.readAllBytes(predictions);
		final String expected="#prediction_bits\t1\t6\n0\t1084227584\t1086324736\t1088421888\t1090519040\t1091567616\t1092616192\n";
		try{
			// The saved oracle differs within tolerance. Output must still contain
			// Java's analytic5..10, not rounded values or a copy of the oracle.
			Files.write(predictions,"0\t5.000001\t6\t7\t8\t9\t10\n".getBytes(StandardCharsets.UTF_8));
			SubnetGpuRoundTripTest.main(new String[]{"net="+image,"vectors="+vectors,
				"predictions="+predictions,"rows=1","ordinary=t","argmax=t","outbits="+accepted});
			if(!expected.equals(new String(Files.readAllBytes(accepted),StandardCharsets.UTF_8)) || Files.exists(acceptedPartial)){
				throw new AssertionError("Published capture must contain the actual checked IEEE754 bits and retire its partial name");
			}
			boolean rejected=false;
			try{SubnetGpuRoundTripTest.main(new String[]{"net="+image,"vectors="+vectors,
				"predictions="+predictions,"rows=1","outbits="+accepted});}
			catch(IllegalArgumentException e){rejected=e.getMessage().contains("Refusing to replace");}
			if(!rejected || !expected.equals(new String(Files.readAllBytes(accepted),StandardCharsets.UTF_8))){
				throw new AssertionError("An existing checked output must survive a repeated invocation unchanged");
			}
			Files.write(vectors,"#dims\t2\t2\t0\n1\t2\t0.7\t0.1\n1\t2\t0.7\t0.1\n".getBytes(StandardCharsets.UTF_8));
			Files.write(predictions,"0\t5\t6\t7\t8\t9\t10\n1\t5\t6\t7\t8\t9\t11\n".getBytes(StandardCharsets.UTF_8));
			rejected=false;
			try{SubnetGpuRoundTripTest.main(new String[]{"net="+image,"vectors="+vectors,
				"predictions="+predictions,"rows=2","outbits="+failed});}
			catch(IllegalStateException e){rejected=e.getMessage().contains("row=1 output=5");}
			if(!rejected || Files.exists(failed) || !expected.equals(new String(Files.readAllBytes(failedPartial),StandardCharsets.UTF_8))){
				throw new AssertionError("A later-row mismatch must retain the checked prefix only as unaccepted partial evidence");
			}
			// The main entry publishes only after its later ordinary/copy checks.
			// Model a later gate throwing after a complete captured first pass.
			try(SubnetGpuRoundTripTest.PredictionBits capture=new SubnetGpuRoundTripTest.PredictionBits(late.toString(),6)){
				capture.add(0,net);
				throw new IllegalStateException("later gate failed");
			}catch(IllegalStateException e){
				if(!e.getMessage().equals("later gate failed")){throw e;}
			}
			if(Files.exists(late) || !Files.exists(latePartial)){
				throw new AssertionError("Closing a complete first-pass capture cannot publish it after a later gate failed");
			}
			System.out.println("CHECKED_PREDICTION_BITS_FIXTURE_PASS");
		}finally{
			Files.write(vectors,originalVectors);
			Files.write(predictions,originalPredictions);
			Files.deleteIfExists(accepted); Files.deleteIfExists(acceptedPartial);
			Files.deleteIfExists(failed); Files.deleteIfExists(failedPartial);
			Files.deleteIfExists(late); Files.deleteIfExists(latePartial);
		}
	}

	/** Bounds cannot conceal changed copy parameters or weaken the exact reference gate. */
	private static void checkExportCopyContract(CellNet net,Path vectors,Path predictions){
		final CellNet exact=net.copy(false);
		SubnetGpuRoundTripTest.requireForwardCopy(net,exact);
		final CellNet tiny=net.copy(false);
		tiny.net[tiny.net.length-1][5].bias=Math.nextUp(Math.nextUp(tiny.net[tiny.net.length-1][5].bias));
		// A known sub-tolerance output difference exercises reporting independent
		// of whether this JVM happens to change its SIMD reduction order in a test.
		final SubnetGpuRoundTripTest.Comparison result=SubnetGpuRoundTripTest.compareExportPaths(
			tiny,vectors.toString(),predictions.toString(),1,.00002,net);
		if(result.maximum<=0 || result.nativeMaximum<=0 || result.maximum>.00002 ||
				result.referenceMaximum!=0 || result.oracleArgmaxChanges!=0 || result.nativeArgmaxChanges!=0){
			throw new AssertionError("Bounded export comparison must retain the measured tiny difference and both oracle checks");
		}
		boolean rejected=false;
		try{SubnetGpuRoundTripTest.requireForwardCopy(net,tiny);}
		catch(IllegalStateException expected){rejected=expected.getMessage().contains("forward parameters");}
		if(!rejected){throw new AssertionError("Even a sub-tolerance copy bias change must fail the exact parameter gate");}
		rejected=false;
		try{SubnetGpuRoundTripTest.compare(tiny,vectors.toString(),null,1,.00002,net);}
		catch(IllegalStateException expected){rejected=expected.getMessage().contains("output bits differ");}
		if(!rejected){throw new AssertionError("Existing reference comparison must remain bit-exact even below tolerance");}
		tiny.net[tiny.net.length-1][5].bias+=1;
		rejected=false;
		try{SubnetGpuRoundTripTest.compareExportPaths(tiny,vectors.toString(),predictions.toString(),1,.00002,net);}
		catch(IllegalStateException expected){rejected=expected.getMessage().contains("output=5");}
		if(!rejected){throw new AssertionError("Bounded mode must reject a beyond-bound sixth-head output");}
		// The reference path must also pass the oracle on the very evaluation used
		// by the native pair; an earlier successful evaluation cannot cover it.
		rejected=false;
		try{SubnetGpuRoundTripTest.compareExportPaths(net,vectors.toString(),predictions.toString(),1,.00002,tiny);}
		catch(IllegalStateException expected){rejected=expected.getMessage().contains("output=5");}
		if(!rejected){throw new AssertionError("A good primary output cannot hide a bad reference-path output");}
	}

	/** Reference comparison must dispatch from each network's own representation. */
	private static void compareMixed(CellNet net, CellNet reference, Path vectors){
		if(SubnetGpuRoundTripTest.compare(net, vectors.toString(), null, 1, 0, reference)!=0){
			throw new AssertionError("Dense/sparse native comparison must preserve all six output bits after unrelated model loads");
		}
	}
}
