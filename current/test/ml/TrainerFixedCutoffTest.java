package test.ml;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;

import ml.Cell;
import ml.CellNet;
import ml.CellNetParser;
import ml.Trainer;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Fixed evaluation cutoffs must initialize fresh networks as well as loaded ones.
 * Finite one-batch runs verify saved metadata without changing backpropagation.
 * @author Brian, Ganyu
 */
public final class TrainerFixedCutoffTest{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, TrainerFixedCutoffTest.class, false);
		final boolean oldExplicit=Trainer.setCutoffForEvaluation;
		final float oldCutoff=Trainer.cutoffForEvaluation;
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String installation=null;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				if(a.equals("installation")){installation=b;}
				else{check(parser.parse(arg, a, b), "Unknown cutoff fixture argument");}
			}
			check(parser.out1!=null && installation!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Supply installation= out=fresh_directory and native runtime options only");
			final File root=new File(parser.out1);
			check(root.mkdir(), "Fixture output directory must be fresh");
			run(root, new File(installation));
			pp.outstream.println("TRAINER_FIXED_CUTOFF_TEST_PASS");
		}finally{
			Trainer.setCutoffForEvaluation=oldExplicit;
			Trainer.cutoffForEvaluation=oldCutoff;
			Shared.closeStream(pp.outstream);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Finite Runs         ----------------*/
	/*--------------------------------------------------------------*/

	private static void run(final File root, final File installation){
		final String data=write(root, "train.tsv", "#dims\t1\t1\t0\n0\t0\n0.25\t0\n0.75\t1\n1\t1\n");
		final String monitor=write(root, "monitor.tsv", "#dims\t1\t1\t0\n0.1\t0\n0.9\t1\n");
		final float trainingCutoff=Cell.cutoffForTraining;
		Trainer.setCutoffForEvaluation=false;
		Trainer.cutoffForEvaluation=.5f;
		final Trainer legacy=new Trainer(new String[]{"in="+data, "dims=1,2,1", "batches=1", "t=1", "simd=f"});
		final CellNet old=legacy.randomNetwork(58);
		check(old.cutoff==-1, "Without an explicit evaluation cutoff, legacy automatic-cutoff initialization remains unchanged");
		// Activation-rate initialization is intentionally once per JVM. Reuse the
		// initializer here; the real CLI cases below each get their own JVM.
		Trainer.setCutoffForEvaluation=true;
		Trainer.cutoffForEvaluation=.25f;
		final CellNet fresh=legacy.randomNetwork(58);
		check(fresh.cutoff==.25f, "A fresh network must carry the requested cutoff before TrainerThread retains it");
		comparePredictions(old, fresh);
		check(Cell.cutoffForTraining==trainingCutoff, "Evaluation-cutoff initialization must not change the training cutoff");

		final String a=path(root, "quarter.bbnet");
		final String b=path(root, "three_quarters.bbnet");
		train(installation, data, monitor, a, .25f);
		train(installation, data, monitor, b, .75f);
		final CellNet first=CellNetParser.load(a, false);
		final CellNet second=CellNetParser.load(b, false);
		check(first.cutoff==.25f && second.cutoff==.75f,
			"The actual one-batch training/statistics/export path preserves both explicit evaluation cutoffs");
		comparePredictions(first, second);
		check(Cell.cutoffForTraining==trainingCutoff, "The reporting fix cannot change backpropagation's cutoff");
	}

	private static void train(final File installation, final String data, final String monitor, final String output, final float cutoff){
		final ArrayList<String> command=new ArrayList<String>();
		command.add("bash");
		command.add(new File(installation, "train.sh").getAbsolutePath());
		Collections.addAll(command, "in="+data, "validate="+monitor, "out="+output, "dims=1,2,1", "batches=1",
			"seed=58", "nets=1", "cycles=1", "scan=0", "alpha=0.001", "alphamult=1", "anneal=0", "balance=0",
			"ptriage=0", "ntriage=0", "sig=1", "tanh=0", "rslog=0", "msig=0", "final=sig", "simd=f",
			"t=1", "ordered=t", "cutoffeval="+cutoff, "ow=f", "-Xmx256m", "-eoom");
		try{
			final Process process=new ProcessBuilder(command).directory(installation).inheritIO().start();
			check(process.waitFor()==0, "Isolated one-batch trainer must finish successfully");
		}catch(InterruptedException failure){
			Thread.currentThread().interrupt();
			throw new IllegalStateException("Interrupted isolated trainer fixture", failure);
		}catch(IOException failure){throw new IllegalStateException("Could not launch isolated trainer fixture", failure);}
	}

	/** Identical seeded weights and one-batch updates must give identical raw values. */
	private static void comparePredictions(final CellNet a, final CellNet b){
		check(a.numInputs()==1 && b.numInputs()==1 && a.numOutputs()==1 && b.numOutputs()==1, "Fixture layout is one input and one output");
		final float[] input=new float[1];
		for(final float value:new float[]{-.5f, 0, .5f, 1, 2}){
			input[0]=value;
			a.applyInput(input);
			b.applyInput(input);
			a.feedForward();
			b.feedForward();
			check(Float.floatToRawIntBits(a.getOutput(0))==Float.floatToRawIntBits(b.getOutput(0)),
				"Changing only the evaluation cutoff must preserve initialized and trained raw predictions");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          File Helpers        ----------------*/
	/*--------------------------------------------------------------*/

	private static String path(final File root, final String name){return new File(root, name).getPath();}
	private static String write(final File root, final String name, final String text){
		final String path=path(root, name);
		try(Output output=new Output(path)){output.writer.print(text);}
		return path;
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
