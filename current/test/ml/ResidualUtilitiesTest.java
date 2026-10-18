package test.ml;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import ml.CellNet;
import ml.CellNetParser;
import ml.CellNetWeightBits;
import ml.RegCalEval;
import ml.SequenceToVector;
import structures.ByteBuilder;

/** Regression controls for deprecated scoring, isolated writer policy and missing models.
 * @author Shinobu */
public final class ResidualUtilitiesTest{

	/** Runs one case against a baseline or repaired runtime in a fresh output directory. */
	public static void main(String[] args) throws Exception{
		if(args.length!=3){throw new IllegalArgumentException("baseline|fixed case fresh-directory");}
		final boolean fixed=args[0].equals("fixed");
		check(fixed || args[0].equals("baseline"), "Unknown phase");
		final Path root=Paths.get(args[2]);
		Files.createDirectory(root);
		final String name=args[1];
		if(name.equals("score-forward") || name.equals("score-reverse")){
			score(fixed, name.equals("score-reverse"));
		}else if(name.equals("models-absent") || name.equals("models-present")){
			models(root, fixed, name.equals("models-absent"));
		}else if(name.startsWith("writer-")){
			writer(root, fixed, name);
		}else{throw new IllegalArgumentException(name);}
		System.out.println("ResidualUtilitiesTest PASS "+name+" "+args[0]);
	}

	/** The rejected orientation must not touch bases; forward inference stays numerically exact. */
	@SuppressWarnings("deprecation")
	private static void score(boolean fixed, boolean reverse){
		final byte[] bases="ACG".getBytes(StandardCharsets.US_ASCII), before=bases.clone();
		final CellNet net=parse(model(true, false));
		boolean rejected=false;
		try{
			final float result=SequenceToVector.score(bases, new float[]{2, 1, 4}, net, reverse);
			check(result==0.25f, "Literal model prediction changed: "+result);
		}catch(AssertionError e){
			if(e.getMessage()==null || !e.getMessage().contains("does not rebuild vec")){throw e;}
			rejected=true;
		}
		check(rejected==(fixed && reverse), "Wrong reverse-complement guard behavior");
		check(Arrays.equals(bases, before), "Scoring changed caller bases");
	}

	/** A dangling listed model deterministically exercises absent-at-load handling without a race. */
	private static void models(Path root, boolean fixed, boolean absent) throws Exception{
		final Path nets=Files.createDirectory(root.resolve("nets")), val=root.resolve("val.tsv");
		Files.write(val, "2\t1\t4\t0.25\n".getBytes(StandardCharsets.US_ASCII));
		if(absent){Files.createSymbolicLink(nets.resolve("net_one.bbnet"), root.resolve("not-created").toAbsolutePath());}
		else{Files.write(nets.resolve("net_one.bbnet"), model(true, false).getBytes(StandardCharsets.US_ASCII));}
		RuntimeException failure=null;
		try{RegCalEval.main(new String[]{"nets="+nets, "val="+val});}
		catch(RuntimeException e){failure=e;}
		if(!absent){check(failure==null, "Ordinary calibration failed: "+failure);}
		else if(fixed){
			check(failure instanceof IllegalStateException && failure.getMessage().contains("No models successfully loaded"),
				"Missing explicit empty-model diagnostic: "+failure);
		}else{check(failure instanceof IndexOutOfBoundsException, "Baseline empty-result failure changed: "+failure);}
	}

	/** Conflicting inherited flags must neither affect conversion nor leak on success/failure. */
	private static void writer(Path root, boolean fixed, String name) throws Exception{
		final boolean sparse=name.contains("sparse"), invalid=name.equals("writer-failure");
		final int bits=name.endsWith("32") ? 32 : name.endsWith("24") ? 24 : 18;
		final Path in=root.resolve("in.bbnet"), out=root.resolve("out.bbnet");
		Files.write(in, model(!sparse, invalid).getBytes(StandardCharsets.US_ASCII));
		final boolean[] original=flags();
		try{
			setFlags(new boolean[]{false, false, true, true, true, true});
			final boolean[] inherited=flags();
			Throwable failure=null;
			try{CellNetWeightBits.main(new String[]{"in="+in, "out="+out, "bits="+bits});}
			catch(AssertionError | IllegalArgumentException | IllegalStateException e){failure=e;}
			if(!fixed){
				check(failure!=null && !Arrays.equals(inherited, flags()), "Baseline writer defect was not exposed");
			}else{
				check(Arrays.equals(inherited, flags()), "Conversion leaked global writer policy");
				if(invalid){
					check(failure instanceof IllegalArgumentException && failure.getMessage().contains("nonfinite"),
						"Expected nonfinite serialization rejection: "+failure);
					check(!Files.exists(out), "Failed serialization published output");
				}else{
					check(failure==null, "Inherited settings broke conversion: "+failure);
					final CellNet decoded=CellNetParser.load(out.toString(), false);
					check(decoded.dense==!sparse && decoded.weightBits()==bits, "Representation/precision changed");
					check(decoded.applyInput(new float[]{2, 1, 4}).feedForward()==0.25f, "Converted prediction changed");
					final String text=new String(Files.readAllBytes(out), StandardCharsets.US_ASCII);
					check(text.contains("#concise\n") && text.contains("#coding A48\n") && !text.contains("#indexencoding"),
						"Conversion did not use its declared writer policy");
				}
			}
		}finally{setFlags(original);}
	}

	/** Builds an exact three-input linear model; the failure case has an infinite bias. */
	private static String model(boolean dense, boolean invalid){
		final String header="##bbnet\n#version 1\n#concise\n#seed 23\n#layers 2\n#blocksize 1\n#dims 3 1\n"+
			(dense ? "#dense\n" : "#sparse\n");
		if(invalid){
			return new ByteBuilder().append(header).append("#coding A48\nC4 LINEAR ")
				.appendFloatA48(Float.POSITIVE_INFINITY).space().appendFloatA48(0.5f)
				.space().appendFloatA48(0).space().appendFloatA48(-0.25f).nl().toString();
		}
		return header+"#coding decimal\n"+(dense ? "C4 LINEAR 0.25 0.5 0 -0.25\n" :
			"I4 0 2\nW4 LINEAR 0.25 0.5 -0.25\n");
	}

	/** Parses literal fixture lines through the production parser. */
	private static CellNet parse(String text){
		final ArrayList<byte[]> lines=new ArrayList<byte[]>();
		for(String line:text.split("\n")){lines.add(line.getBytes(StandardCharsets.US_ASCII));}
		return CellNetParser.loadFromLines(lines);
	}

	/** Captures all six settings changed by precision conversion. */
	private static boolean[] flags(){
		return new boolean[]{CellNet.codingA48Out, CellNet.CONCISE, CellNet.OUT_DENSE,
			CellNet.OUT_SPARSE, CellNet.OUT_HEX, CellNet.OUT_DELTA_A48};
	}
	/** Restores the saved writer settings in their documented order. */
	private static void setFlags(boolean[] flags){
		CellNet.codingA48Out=flags[0]; CellNet.CONCISE=flags[1]; CellNet.OUT_DENSE=flags[2];
		CellNet.OUT_SPARSE=flags[3]; CellNet.OUT_HEX=flags[4]; CellNet.OUT_DELTA_A48=flags[5];
	}
	/** Throws independently of assertion enablement when a test oracle fails. */
	private static void check(boolean condition, String message){if(!condition){throw new AssertionError(message);}}
}
