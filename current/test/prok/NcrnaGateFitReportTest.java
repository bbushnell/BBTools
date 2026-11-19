package test.prok;

import java.io.File;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaGateFitReport;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Independent two-fit collection, byte preservation and corrupted-receipt controls.
 * @author Brian Bushnell, Ganyu, Raiden
 */
public final class NcrnaGateFitReportTest{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaGateFitReportTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown collection fixture option");
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && parser.in1==null && parser.in2==null
				&& parser.out2==null && parser.maxReads<0 && parser.samplerate==1, "Use out=fresh_directory and native runtime flags");
			final File root=new File(parser.out1);
			check(root.mkdir(), "Fixture directory must be fresh");
			run(root);
			pp.outstream.println("NCRNA_GATE_FIT_REPORT_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	/*--------------------------------------------------------------*/
	/*----------------          Finite Cases        ----------------*/
	/*--------------------------------------------------------------*/

	/** Two complete fake receipts test collection without creating or fitting nets. */
	private static void run(final File root){
		final String input=file(root, "cells.tsv");
		write(input, "cell\tarm\tend\tseed\toutput\trate\trole\tjob\nr\tR\t5prime\t58\trslog\t0.003\tCANARY\tunused\nt\tT\t5prime\t58\trslog\t0.001\tREMAINING\tunused\n");
		for(final String cell:new String[]{"r", "t"}){
			final File run=new File(root, "runs/"+cell);
			check(new File(run, "replay").mkdirs(), "Fresh run/replay fixture directory");
			check(new File(root, "submitted/"+cell).mkdirs(), "Fresh submitted receipt directory");
			write(file(run, "exit.txt"), "0\n");
			write(file(run, "complete.txt"), "GATE_CONVERGENCE_FIT_AND_SAVED_REPLAY_COMPLETE_OUTER_UNPERFORMED\n");
			write(file(run, "source.txt"), SOURCE+"\n");
			write(file(run, "build.seal"), SEAL+"\n");
			write(file(run, "job.txt"), cell.equals("r") ? "100\n" : "101\n");
			write(file(root, "submitted/"+cell+"/job.txt"), cell.equals("r") ? "100\n" : "101\n");
			write(file(run, "gate.bbnet"), "##bbnet\n#seed "+(cell.equals("r") ? 58 : 65)+"\n#dims 32 32 32 16 1\n");
			write(file(run, "train.stderr"), log(cell.equals("r")));
			for(final String stage:new String[]{"train", "score", "replay"}){write(file(run, stage+".time.tsv"), "0.01\t0.02\t0.03\n");}
			write(file(run, "replay/summary.tsv"), curves());
			write(file(run, "replay/matched_rates.tsv"), "matched_to\ttarget_moves\ttarget_move_pct\tmatch_status\tcutoff_operator\ttied_at_cutoff\t"+curves().split("\n", 2)[0]+"\n"
				+matched("margin0")+matched("margin0.1")+matched("all_proposals"));
			write(file(run, "replay/ranking.tsv"), "end\tproposals\timproving\tnon_improving\timproving_pct\tauc_higher_score_predicts_improvement\tqualification\n5prime\t50183\t37230\t12953\t74.188470\t0.8203283472\tTIES_COUNT_HALF_NO_NEW_INFERENCE\n");
		}
		collect(input, file(root, "good"), 2);
		check(read(file(root, "good/summary.tsv")).contains("\n2\t52\t6\t2\t"), "Two cells contribute all52curves,6matched rows and2AUC rows");
		check(read(file(root, "good/ranking.tsv")).contains("t\tT\t5prime\t58\trslog\t0.001\t5prime\t50183\t37230\t12953\t74.188470\t0.8203283472\t"),
			"Scientific values and precision are copied literally, not rounded or recomputed");
		check(read(file(root, "good/index.tsv")).contains("t\tT\t5prime\t58\trslog\t0.001\t101\t65\tfpr\t"),
			"The winning seed65 is distinct from the manifest's seed-set anchor58");
		check(read(file(root, "good/monitor_curves.tsv")).contains("\t150000\t8\t\t\t0.1\t0.8\t"),
			"Eight unlabelled starts produce one honest min/max envelope, not invented per-seed curves");
		try(Input monitor=new Input(file(root, "good/monitor_curves.tsv"))){
			monitor.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line; (line=monitor.reader.nextLine())!=null;){
				p.set(line);
				check(p.terms()==21, "All declared monitor columns survive optional metric omission");
				if(p.termEquals("R", 1)){
					check(p.parseDouble(8)==.12 && p.parseDouble(9)==.12 && p.parseDouble(10)==.14 && p.parseDouble(11)==.14
						&& p.parseDouble(18)==0 && p.parseDouble(19)==0, "The first observation supplies real MSE bounds and preserves measured zero learning rate");
					for(int i=12; i<=17; i++){check(p.length(i)==0, "Regression supplies no WER/FPR/FNR bounds");}
				}else{
					check(p.length(8)==0 && p.length(9)==0 && p.parseDouble(10)==.1 && p.parseDouble(11)==.8,
						"Classification supplies no training MSE and retains all eight native-error extrema");
				}
			}
		}
		for(final String invalid:new String[]{"NaN", "Infinity", "-Infinity"}){
			write(file(root, "runs/r/train.stderr"), log(true).replace("trainMSE=0.12", "trainMSE="+invalid));
			final String name="nonfinite_monitor_"+invalid;
			fails(()->collect(input, file(root, name), 2));
			check(!new File(root, name+"/summary.tsv").exists(), "A computed nonfinite monitor value cannot become a missing measurement");
		}
		write(file(root, "runs/r/train.stderr"), log(true));
		fails(()->collect(input, file(root, "missing_cell"), 3));
		write(file(root, "runs/t/source.txt"), "1111111111111111111111111111111111111111\n");
		fails(()->collect(input, file(root, "wrong_source"), 2));
		write(file(root, "runs/t/source.txt"), SOURCE+"\n");
		write(file(root, "runs/t/replay/summary.tsv"), curves().replaceFirst("0.123456789", "0.987654321"));
		fails(()->collect(input, file(root, "changed_control"), 2));
		write(file(root, "runs/t/replay/summary.tsv"), curves());
		write(file(root, "runs/t/exit.txt"), "1\n");
		fails(()->collect(input, file(root, "failed_job"), 2));
		write(file(root, "runs/t/exit.txt"), "0\n");
		write(file(root, "runs/t/train.stderr"), log(false).replaceFirst("Batch 150k:[^\\n]*\\n", ""));
		fails(()->collect(input, file(root, "missing_monitor_start"), 2));
		for(final String name:new String[]{"missing_cell", "wrong_source", "changed_control", "failed_job", "missing_monitor_start"}){
			check(!new File(root, name+"/summary.tsv").exists(), "A broken population or receipt cannot publish a completion summary");
		}
		write(file(root, "runs/t/train.stderr"), log(false).replace("150000 batches", "300000 batches").replace("Batch 150k:", "Batch 300k:"));
		fails(()->collect(input, file(root, "long_wrong_budget"), 2));
		new NcrnaGateFitReport(new String[]{"in="+input, "out="+file(root, "long_good"), "expected=2",
			"source="+SOURCE, "seal="+SEAL, "trainerselection=wer", "trainerbatches=300k"}).process();
		check(read(file(root, "long_good/monitor_curves.tsv")).contains("\t300000\t8\t\t\t0.1\t0.8\t"),
			"Longer training retains all eight actual terminal observations at its declared budget");
		check(read(file(root, "long_good/index.tsv")).contains("\t300000\tINTERNAL_MONITOR_SELECTION_OUTER_UNPERFORMED"),
			"The index must state the longer training budget explicitly");
		write(file(root, "runs/t/train.stderr"), log(false).replace("Batch 150k:", "Batch 300k:"));
		fails(()->new NcrnaGateFitReport(new String[]{"in="+input, "out="+file(root, "header_mismatch"), "expected=2",
			"source="+SOURCE, "seal="+SEAL, "trainerselection=wer", "trainerbatches=300k"}).process());
		fails(()->new NcrnaGateFitReport(new String[]{"in="+input, "out="+file(root, "zero_budget"), "expected=2",
			"source="+SOURCE, "seal="+SEAL, "trainerselection=wer", "trainerbatches=0"}));
		for(final String name:new String[]{"long_wrong_budget", "header_mismatch", "zero_budget"}){
			check(!new File(root, name+"/summary.tsv").exists(), "A mismatched training budget cannot publish completion");
		}
	}
	/** Literal terminal logs distinguish one MSE observation from eight ERR rows. */
	private static String log(final boolean regression){
		if(regression){return "loaded 196595 training samples, 50183 external validation samples, dims=[32, 32, 32, 16, 1]\nepoch 800 lr=0.00000 trainMSE=0.12 valMSE=0.14 best=0.13\nround-trip check vs in-memory model: maxDiff=0.0 over 1000 samples (OK)\n";}
		final ByteBuilder text=new ByteBuilder("Loaded samples: training=196595 validation=50183\nUsing 64 worker threads, 8 trainer threads, 22 jobs/trainer/batch, 150000 batches, 1 cycle to train 8 networks total.\n");
		for(int i=1; i<=8; i++){text.append("Batch 150k:\tERR= 0.").append(i).append("\tWER= 1.2\tFPR= 0.2\tFNR= 0.1\tCTF= 0.5\tALP= 0.00003\t0.05s\n");}
		return text.toString();
	}
	private static void collect(final String input, final String output, final int expected){
		new NcrnaGateFitReport(new String[]{"in="+input, "out="+output, "expected="+expected, "source="+SOURCE, "seal="+SEAL, "trainerselection=fpr"}).process();
	}
	/** Synthetic forty-column framing retains extra decimal digits as a byte oracle. */
	private static String curves(){
		final ByteBuilder text=new ByteBuilder("policy\tthreshold\tend\ttrials\tproposals");
		for(int i=5; i<40; i++){text.tab().append("field").append(i);}
		text.nl();
		for(int i=0; i<26; i++){text.append(row(i<4 ? new String[]{"do_nothing", "margin0", "margin0.1", "all_proposals"}[i] : "gate"));}
		return text.toString();
	}
	private static String matched(final String name){return name+"\t1\t0.001\tEXACT\t>=\t1\t"+row("gate_matched_rate");}
	private static String row(final String policy){
		final ByteBuilder text=new ByteBuilder(policy).append("\tNA\t5prime\t88966\t50183\t0.123456789");
		for(int i=6; i<40; i++){text.append("\t0");}
		return text.nl().toString();
	}
	/*--------------------------------------------------------------*/
	/*----------------          File Helpers        ----------------*/
	/*--------------------------------------------------------------*/

	private static String file(final File root, final String name){return new File(root, name).getPath();}
	/** Deliberately replace synthetic receipts inside this fixture's fresh root. */
	private static void write(final String path, final String text){
		try(Output out=new Output(new ByteStreamWriter(path, true, false, false))){out.writer.print(text);}
	}
	private static String read(final String path){
		final ByteBuilder text=new ByteBuilder();
		try(Input in=new Input(path)){for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){text.append(line).nl();}}
		return text.toString();
	}
	private static void fails(final Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Broken completed-fit evidence must fail loudly");
	}
	private static void check(final boolean good, final String why){if(!good){throw new AssertionError(why);}}
	/*--------------------------------------------------------------*/
	/*----------------       Fixture Constants      ----------------*/
	/*--------------------------------------------------------------*/

	/** Deliberately synthetic Git/source and build identities, never real evidence. */
	private static final String SOURCE="0000000000000000000000000000000000000000";
	private static final String SEAL="00000000000000000000";
}
