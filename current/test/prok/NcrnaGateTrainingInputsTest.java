package test.prok;

import java.io.File;
import java.util.Locale;

import ml.CellNet;
import ml.CellNetParser;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaGateTrainingInputs;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Literal training-only scaling, unchanged concatenation and serialized parity.
 * @author Brian, Ganyu
 */
public final class NcrnaGateTrainingInputsTest{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaGateTrainingInputsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				check(parser.parse(arg, a, b), "Unknown gate-input fixture argument");
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use out=fresh_directory and native runtime flags only");
			final File root=new File(parser.out1);
			check(root.mkdir(), "Fresh fixture directory required");
			run(root);
			pp.outstream.println("NCRNA_GATE_TRAINING_INPUTS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Finite Cases        ----------------*/
	/*--------------------------------------------------------------*/

	private static void run(final File root){
		final String first=row(0, 0);
		final String second=row(2, 1);
		final String a=write(root, "first.tsv", HEADER+first);
		final String b=write(root, "second.tsv", HEADER+second);
		final String monitor=write(root, "monitor.tsv", HEADER+row(100, 1));
		final String prepared=path(root, "prepared");
		final String[] args={"in="+a+','+b, "monitor="+monitor, "expected=1,1", "monitorexpected=1", "seed=58", "out="+prepared};
		new NcrnaGateTrainingInputs(args).process();
		final String[] duplicate=args.clone();
		duplicate[0]="in="+a+','+a;
		duplicate[5]="out="+path(root, "duplicate_training");
		failsDuplicate(()->new NcrnaGateTrainingInputs(duplicate));
		duplicate[0]=args[0];
		duplicate[1]="monitor="+a;
		duplicate[5]="out="+path(root, "monitor_leak");
		failsDuplicate(()->new NcrnaGateTrainingInputs(duplicate));
		check(read(prepared+"/train.raw.tsv.gz").equals(HEADER+first+second), "Raw concatenation preserves every input byte, target and order");
		check(read(prepared+"/training_feature_ranges.tsv").contains("TRAINING_FOLDS_ONLY\tfirst_stage_scores\tscore0\t0\t2\t0.0\t0.0\t0.0\t2.0\t2.0\t0\n"),
			"Training range quantiles use the two teaching values and nearest ranks only");
		check(read(prepared+"/monitor_feature_ranges.tsv").contains("MONITOR_NOT_USED_TO_FIT_SCALING\tfirst_stage_scores\tscore0\t0\t1\t100.0\t100.0\t100.0\t100.0\t100.0\t1\n"),
			"Monitor range reports its actual outlier without clipping");
		check(read(prepared+"/all_feature_ranges.tsv").contains("ALL_FOLDS_DIAGNOSTIC_ONLY\tfirst_stage_scores\tscore0\t0\t3\t0.0\t0.0\t2.0\t100.0\t100.0\t1\n"),
			"Combined diagnostics include the monitor but cannot alter the earlier frozen transform");
		check(read(prepared+"/summary.tsv").contains("\n2\t1\t1\t1\t32\t8\t1\tTRAINING_MINMAX_FLOAT32_EIGHT_DISTINCT_UNTRAINED_SEEDS_NO_CLIPPING\n"),
			"The monitor outlier is counted once and retained in the full population");
		for(int seed=58; seed<66; seed++){
			final String file=prepared+"/seed"+seed+".bbnet";
			final CellNet net=CellNetParser.load(file, false);
			check(net.numInputs()==32 && net.numOutputs()==1 && net.hasInputNormalization(), "Each fresh start carries the explicit32-input transform");
			final float[] mean=net.inputMeanCopy();
			final float[] inverse=net.inputInverseStdCopy();
			check(mean[0]==0 && inverse[0]==0.5f && mean[1]==5 && inverse[1]==1,
				"Training0/2 maps to0/1; constant5 maps to0 with scale1; monitor100 must not fit either statistic");
			for(int i=2; i<32; i++){check(mean[i]==0 && inverse[i]==1, "Constant-zero columns normalize to zero without enormous inverse deviations");}
			final String serialized=read(file);
			check(serialized.contains("\n#seed "+seed+"\n") && !serialized.contains("\n#epochs ")
				&& serialized.contains("##initialization fresh-random-untrained"), "Eight explicit seeds are untrained starts, not repeated trained checkpoints");
			final CellNet manual=net.copy(false);
			manual.setInputNormalization(null, null);
			final float[] raw=new float[32];
			raw[0]=100;
			raw[1]=5;
			final float[] normalized=raw.clone();
			for(int i=0; i<32; i++){normalized[i]=(normalized[i]-mean[i])*inverse[i];}
			check(normalized[0]==50 && normalized[1]==0, "Out-of-range monitor values remain extrapolated, never silently clipped");
			net.applyInput(raw);
			manual.applyInput(normalized);
			net.feedForward();
			manual.feedForward();
			check(Float.floatToRawIntBits(net.getOutput(0))==Float.floatToRawIntBits(manual.getOutput(0))
				&& raw[0]==100 && raw[1]==5, "Serialized normalization equals manual float32 preprocessing and preserves borrowed raw inputs");
		}
		args[2]="expected=2,1";
		args[5]="out="+path(root, "wrong_count");
		fails(()->new NcrnaGateTrainingInputs(args).process());
		check(!new File(path(root, "wrong_count/summary.tsv")).exists(), "A truncated training population cannot publish a completed preparation");
		args[2]="expected=1,1";
		args[1]="monitor="+write(root, "bad_monitor.tsv", HEADER+row(100, 2));
		args[5]="out="+path(root, "bad_label");
		fails(()->new NcrnaGateTrainingInputs(args).process());
	}

	/*--------------------------------------------------------------*/
	/*----------------          File Helpers        ----------------*/
	/*--------------------------------------------------------------*/

	private static String row(final int first, final int label){
		final ByteBuilder text=new ByteBuilder().append(first).append("\t5\t");
		for(int i=2; i<32; i++){text.append("0\t");}
		return text.append(label).nl().toString();
	}
	private static String path(final File root, final String name){return new File(root, name).getPath();}
	private static String write(final File root, final String name, final String text){
		final String file=path(root, name);
		try(Output output=new Output(file)){output.writer.print(text);}
		return file;
	}
	private static String read(final String path){
		final ByteBuilder text=new ByteBuilder();
		try(Input input=new Input(path)){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){text.append(line).nl();}
		}
		return text.toString();
	}
	private static void fails(final Runnable run){
		try{run.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid gate training inputs must fail loudly");
	}
	private static void failsDuplicate(final Runnable run){
		try{run.run();}
		catch(RuntimeException expected){
			check(expected.getMessage().contains("Duplicate file"), "The native input guard must reject duplicate training/monitor names");
			return;
		}
		throw new AssertionError("Duplicate gate populations must fail before output creation");
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private static final String HEADER="#dims\t32\t1\t0\n";
}
