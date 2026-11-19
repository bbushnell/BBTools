package test.prok;

import java.io.File;
import java.util.Locale;

import ml.CellNet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.ScoreNcrnaCenteredGate;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Literal scored/held counts and raw-bit cache validation with an injected
 * constant network. No fit or scientific first-stage inference.
 * @author Brian, Ganyu
 */
public final class ScoreNcrnaCenteredGateTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, ScoreNcrnaCenteredGateTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, pair[0].toLowerCase(Locale.ROOT), pair.length>1 ? pair[1] : null), "Unknown finite score-test argument");
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null, "Finite controls require out=fresh_directory only");
			final File root=new File(parser.out1);
			check(root.mkdir(), "Fresh finite score-test directory required");
			run(root);
			pp.outstream.println("NCRNA_CENTERED_GATE_SCORES_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final File root){
		final File inputs=new File(root, "input_fixture");
		NcrnaCenteredGateInputsTest.main(new String[]{"out="+inputs.getPath(), "t=1"});
		for(int fold=0; fold<5; fold++){
			final ByteBuilder scores=new ByteBuilder();
			try(Input audit=new Input(new File(inputs, "fold"+fold+".rows.tsv.gz").getPath())){
				for(byte[] row; (row=audit.reader.nextLine())!=null && row[0]=='#';){scores.append(row).nl();}
			}
			scores.append("row1\tstatus");
			for(int i=0; i<22; i++){scores.append("\tscore_bits_").append(i);}
			scores.nl();
			for(int trial=1; trial<=2; trial++){
				scores.append(trial).append("\tPROPOSAL");
				for(int i=0; i<22; i++){scores.append("\t0");}
				scores.nl();
			}
			scores.append("3\tFIXED\n#COMPLETE\t3\t2\n");
			write(new File(inputs, "fold"+fold+".scores.tsv.gz").getPath(), scores);
		}
		final String net=new File(root, "injected.net").getPath();
		write(net, new ByteBuilder("Finite injected network; never parsed.\n"));
		final String[] args={"in="+new File(inputs, "manifest.tsv").getPath(), "out="+new File(root, "result").getPath(),
			"net="+net, "member=RF00174", "end=5prime", "seed=59", "partition=internal", "t=1", "unpigz=f", "pigz=f"};
		final FixedNet fixed=new FixedNet(.75f);
		new ScoreNcrnaCenteredGate(args).process(fixed);
		check(fixed.calls==10, "Exactly ten proposals invoke the gate; the five fixed holds never do");
		int rows=0, proposals=0;
		boolean complete=false;
		try(Input cache=new Input(new File(root, "result/scores.tsv.gz").getPath())){
			final LineParser1 p=new LineParser1('\t');
			for(byte[] row; (row=cache.reader.nextLine())!=null;){
				p.set(row);
				if(p.termEquals("#COMPLETE", 0)){
					check(p.terms()==3 && p.parseInt(1)==15 && p.parseInt(2)==10 && cache.reader.nextLine()==null,
						"Cache completion preserves all-call and inference counts separately");
					complete=true;
					break;
				}
				if(row[0]=='#' || p.termEquals("dataset", 0)){continue;}
				check(p.terms()==53 && p.parseInt(1)==rows/3 && p.parseInt(2)==rows%3+1, "Cache retains original dataset/fold/row order and complete53-field layout");
				check(p.length(17)==0, "A new scored cache must omit the unbounded maximum rather than copy an extreme control");
				if(rows%3<2){
					check(p.parseInt(30)==Float.floatToRawIntBits(.75f), "Each cached gate score retains exact FP32 bits");
					for(int i=31; i<53; i++){check(p.parseInt(i)==0, "Raw first-stage evidence is retained without reevaluation");}
					proposals++;
				}else{
					for(int i=30; i<53; i++){check(p.length(i)==0, "A fixed hold must never fabricate gate or first-stage scores");}
				}
				rows++;
			}
		}
		check(complete && rows==15 && proposals==10, "Every matched call survives the gate-score cache");
		args[1]="out="+new File(root, "bad_nonfinite").getPath();
		boolean rejected=false;
		try{new ScoreNcrnaCenteredGate(args).process(new FixedNet(Float.NaN));}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("finite gate evidence");}
		check(rejected, "Nonfinite gate output must fail loudly, never become an implicit hold");
		args[1]="out="+new File(root, "bad_partition").getPath();
		args[6]="partition=outer";
		rejected=false;
		try{new ScoreNcrnaCenteredGate(args).process(new FixedNet(.75f));}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("cannot mix internal and outer");}
		check(rejected, "Internal folds cannot be relabelled as an outer comparison");
	}

	private static final class FixedNet extends CellNet{
		FixedNet(final float value_){value=value_;}
		@Override public int numInputs(){return 35;}
		@Override public int numOutputs(){return 1;}
		@Override public CellNet applyInput(final float[] values){
			ScoreNcrnaCenteredGateTest.check(values.length==35 && values[32]==2 && values[33]==2 && values[34]==0,
				"The gate receives the distinct centered35layout");
			return this;
		}
		@Override public float feedForward(){calls++;return value;}
		@Override public float getOutput(final int index){return value;}
		final float value;
		int calls;
	}
	private static void write(final String path, final ByteBuilder rows){try(Output out=new Output(path)){out.writer.print(rows);}}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
