package test.prok;

import java.io.PrintStream;
import java.util.Arrays;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Streams two aligned NcrnaMlCompatibilityCapture files without evaluating nets.
 * Every nonprediction record must match exactly. Prediction differences retain
 * network, mode, row and output identity; first-maximum changes are independent
 * of absolute error. This reports evidence, never a tolerance or noise waiver.
 * Adapted from template.A_SampleByteFile for two lockstep byte readers.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaMlCaptureCompare{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns redirected diagnostics even when parsing or comparison fails. */
	public static void main(final String[] args){
		final PreParser pre=new PreParser(args, null, false);
		try{
			new NcrnaMlCaptureCompare(pre.args, pre.outstream).process();
		}finally{
			Shared.closeStream(pre.outstream);
		}
	}

	/** Uses standard input/output flags; both inputs and the report must differ. */
	private NcrnaMlCaptureCompare(final String[] args, final PrintStream diagnostics){
		outstream=diagnostics;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase();
			final String b=split.length>1 ? split[1] : null;
			if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter "+arg);}
		}
		leftPath=parser.in1;
		rightPath=parser.in2;
		report=parser.out1;
		if(leftPath==null || rightPath==null || report==null || parser.out2!=null
				|| parser.overwrite || parser.append){
			throw new IllegalArgumentException("Require in=baseline in2=candidate out=fresh_report ow=f");
		}
		if(!Tools.testInputFiles(false, true, leftPath, rightPath)
				|| !Tools.testOutputFiles(false, false, false, report)
				|| !Tools.testForDuplicateFiles(true, leftPath, rightPath, report)){
			throw new IllegalArgumentException("Capture inputs and fresh report must be distinct accessible files");
		}
		assert(!leftPath.equals(rightPath)) : "A self-file comparison cannot establish cross-capture compatibility";
	}

	/*--------------------------------------------------------------*/
	/*----------------          Comparison          ----------------*/
	/*--------------------------------------------------------------*/

	/** Preserves all differences and refuses malformed, truncated or misaligned files. */
	private void process(){
		assert(leftPath!=null && rightPath!=null) : "Both capture readers depend on constructor path validation";
		final LineParser1 left=new LineParser1('\t');
		final LineParser1 right=new LineParser1('\t');
		long lines=0;
		int networks=0;
		boolean complete=false;
		try(NcrnaTextIO.Input a=new NcrnaTextIO.Input(leftPath);
				NcrnaTextIO.Input b=new NcrnaTextIO.Input(rightPath);
				NcrnaTextIO.Output out=new NcrnaTextIO.Output(report)){
			write(out, buffer.append("kind\tnetwork\tmode\trow\toutput\tleft\tright\tabs_difference\n"));
			for(byte[] x=a.reader.nextLine(), y=b.reader.nextLine(); x!=null || y!=null;
					x=a.reader.nextLine(), y=b.reader.nextLine()){
				lines++;
				if(x==null || y==null || complete){throw new IllegalArgumentException("Capture length/end mismatch at line "+lines);}
				left.set(x);
				right.set(y);
				final boolean prediction=left.termEquals("PREDICTIONS", 0);
				final boolean training=left.termEquals("TRAIN_PREDICTION", 0);
				if(prediction || training){
					compare(out, left, right, training);
				}else{
					if(!Arrays.equals(x, y)){
						stateDifferences++;
						write(out, buffer.append("STATE\t").append(network).tab().append(mode)
							.tab().append(row).tab().append("line=").append(lines).tab()
							.append("field=").append(left.parseString(0)).tab().append("NA\tNA\n"));
						write(out, buffer.append("STATE_LEFT\t").append(lines).tab().append(x).nl());
						write(out, buffer.append("STATE_RIGHT\t").append(lines).tab().append(y).nl());
					}
					if(left.termEquals("NETWORK", 0)){
						flush(out);
						network=left.parseString(1);
						networks++;
					}else if(left.termEquals("LEGACY_TANH", 0)){
						flush(out);
						network="legacy_tanh_"+left.parseString(1);
					}else if(left.termEquals("MODE", 0)){
						flush(out);
						mode=left.parseString(1);
						row=-1;
						trainingRow=0;
					}else if(left.termEquals("DIMS", 0)){
						outputs=left.parseInt(left.terms()-1);
					}else if(left.termEquals("ROW", 0)){
						final int next=left.parseInt(1);
						if(pendingPrediction || next!=row+1){throw new IllegalArgumentException("Missing/repeated/out-of-order captured row for "+network+" "+mode);}
						row=next;
						pendingPrediction=true;
						rawSeen=false;
					}else if(left.termEquals("RAW", 0)){
						if(!pendingPrediction || rawSeen){throw new IllegalArgumentException("Input row framing mismatch for "+network+" "+mode);}
						rawSeen=true;
					}else if(left.termEquals("COMPLETE", 0)){
						if(left.parseInt(1)!=networks || networks<1){throw new IllegalArgumentException("Capture network count does not conserve");}
						complete=true;
					}
				}
			}
			flush(out);
			if(!complete){throw new IllegalArgumentException("Capture lacks its terminal COMPLETE record");}
			write(out, buffer.append("TOTAL\tvalues=").append(totalValues).tab()
				.append("different_bits=").append(totalDifferences).tab()
				.append("max_abs=").append(Double.toString(totalMaximum)).tab()
				.append("first_maximum_changes=").append(totalSelections).tab()
				.append("state_record_differences=").append(stateDifferences).nl());
		}
		if(stateDifferences>0){throw new IllegalStateException("Nonprediction capture records differ: "+stateDifferences+"; see "+report);}
		outstream.println("CAPTURE_COMPARISON_COMPLETE_NOT_ACCEPTANCE networks="+networks);
	}

	/** Decodes every raw FP32 output and retains the first maximum on ties. */
	private void compare(final NcrnaTextIO.Output out, final LineParser1 left,
			final LineParser1 right, final boolean training){
		final int first=training ? 2 : 1;
		final int count=left.terms()-first;
		if(count<1 || count!=outputs || right.terms()!=left.terms()
				|| !right.termEquals(training ? "TRAIN_PREDICTION" : "PREDICTIONS", 0)
				|| (training && (left.parseInt(1)!=count || right.parseInt(1)!=count))
				|| (!training && (count!=outputs || row<0 || !pendingPrediction || !rawSeen))){
			throw new IllegalArgumentException("Prediction shape/framing mismatch for "+network+" "+mode+" row "+row);
		}
		assert(count>0) : "A first-maximum comparison requires at least one recorded output";
		final int position=training ? trainingRow++ : row;
		int leftBest=0;
		int rightBest=0;
		float leftMax=Float.NEGATIVE_INFINITY;
		float rightMax=Float.NEGATIVE_INFINITY;
		for(int i=0; i<count; i++){
			final int x=left.parseInt(first+i);
			final int y=right.parseInt(first+i);
			final float a=Float.intBitsToFloat(x);
			final float b=Float.intBitsToFloat(y);
			if(!Float.isFinite(a) || !Float.isFinite(b)){throw new IllegalArgumentException("Nonfinite captured prediction");}
			if(a>leftMax){leftMax=a; leftBest=i;}
			if(b>rightMax){rightMax=b; rightBest=i;}
			final double delta=Math.abs((double)a-b);
			if(training){
				trainingMaximum=Math.max(trainingMaximum, delta);
				trainingValues++;
			}else{
				maximum=Math.max(maximum, delta);
				values++;
			}
			if(x!=y){
				if(training){trainingDifferences++;}
				else{differences++;}
				write(out, buffer.append(training ? "TRAIN_OUTPUT\t" : "OUTPUT\t")
					.append(network).tab().append(mode).tab().append(position).tab().append(i)
					.tab().append(Float.toString(a)).tab().append(Float.toString(b))
					.tab().append(Double.toString(delta)).nl());
			}
		}
		if(leftBest!=rightBest){
			if(training){trainingSelections++;}
			else{selections++;}
			write(out, buffer.append(training ? "TRAIN_SELECTION\t" : "SELECTION\t").append(network).tab().append(mode)
				.tab().append(position).tab().append("NA").tab().append(leftBest).tab().append(rightBest).tab().append("NA\n"));
		}
		if(!training){pendingPrediction=false;}
	}

	/** Emits one bounded per-network/mode summary, including unchanged groups. */
	private void flush(final NcrnaTextIO.Output out){
		assert(differences<=values && totalDifferences<=totalValues) : "Per-mode summaries must conserve differing versus inspected output counts";
		if(pendingPrediction){throw new IllegalArgumentException("Capture ended a mode without its pending prediction");}
		if(values>0){summary(out, "SUMMARY", values, differences, maximum, selections);}
		if(trainingValues>0){summary(out, "TRAIN_SUMMARY", trainingValues, trainingDifferences, trainingMaximum, trainingSelections);}
		totalValues+=values+trainingValues;
		totalDifferences+=differences+trainingDifferences;
		totalSelections+=selections+trainingSelections;
		totalMaximum=Math.max(totalMaximum, Math.max(maximum, trainingMaximum));
		values=0;
		differences=0;
		selections=0;
		maximum=0;
		trainingValues=0;
		trainingDifferences=0;
		trainingSelections=0;
		trainingMaximum=0;
	}

	/** Keeps inference and seeded training-control denominators distinct. */
	private void summary(final NcrnaTextIO.Output out, final String kind, final long count,
			final long different, final double max, final long selected){
		assert(different<=count && count>0) : "Summary consumers require changed outputs to be a subset of inspected outputs";
		write(out, buffer.append(kind).tab().append(network).tab().append(mode)
			.tab().append("values=").append(count).tab().append("different_bits=").append(different)
			.tab().append("max_abs=").append(Double.toString(max))
			.tab().append("first_maximum_changes=").append(selected).nl());
	}

	/** The native writer copies the reused byte builder before return. */
	private void write(final NcrnaTextIO.Output out, final ByteBuilder text){
		assert(text.length()>0) : "A comparison report record must be framed before writing";
		out.writer.print(text);
		text.clear();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final String leftPath;
	private final String rightPath;
	private final String report;
	private final PrintStream outstream;
	private final ByteBuilder buffer=new ByteBuilder();
	private String network="none";
	private String mode="none";
	private int row=-1;
	private int outputs=-1;
	private int trainingRow;
	private boolean pendingPrediction, rawSeen;
	private long values, differences, selections, totalValues, totalDifferences, totalSelections, stateDifferences;
	private double maximum, totalMaximum;
	private long trainingValues, trainingDifferences, trainingSelections;
	private double trainingMaximum;
}
