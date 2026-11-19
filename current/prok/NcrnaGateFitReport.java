package prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Lossless collection of completed gate-fit tables. No inference, threshold
 * selection, endpoint-metric recomputation or pooling of models/populations.
 * Native monitor logs additionally yield explicitly labelled min/max envelopes.
 * Each row retains the literal fit identity from the submitted-cell manifest.
 * Failed attempts remain in their native directories, outside that manifest.
 * Monitor metrics absent from the declared trainer arm have empty bound cells;
 * the first actual observation initializes each supplied metric, including zero.
 * @author Brian Bushnell, Ganyu, Raiden
 */
public final class NcrnaGateFitReport{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaGateFitReport.class, false);
		try{
			new NcrnaGateFitReport(pp.args).process();
			pp.outstream.println("NCRNA_GATE_FIT_REPORT_COMPLETE_SAVED_TABLES_ONLY_NO_SELECTION");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Bind one complete cell manifest and the immutable fit provenance. */
	public NcrnaGateFitReport(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("expected")){expected=Integer.parseInt(b);}
			else if(a.equals("source")){source=b;}
			else if(a.equals("seal")){seal=b;}
			else if(a.equals("trainerselection")){trainerSelection=b;}
			else if(a.equals("trainerbatches")){trainerBatches=Parse.parseIntKMG(b);}
			else{require(parser.parse(arg, a, b), "Unknown gate-fit report argument: "+arg);}
		}
		input=parser.in1;
		output=parser.out1;
		require(input!=null && output!=null && expected>0 && source!=null && source.matches("[0-9a-f]{40}")
			&& seal!=null && seal.matches("[0-9a-f]{20}") && !new File(output).exists()
			&& ("fpr".equals(trainerSelection) || "wer".equals(trainerSelection))
			&& trainerBatches>0 && trainerBatches<=1000000
			&& parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
			&& parser.maxReads<0 && parser.samplerate==1,
			"Use complete in=cells.tsv out=fresh_directory expected=count source=fit_commit seal=fit_build_sha80 trainerselection=fpr|wer trainerbatches=1..1000000");
		require(Tools.testInputFiles(false, true, input), "The literal submitted-cell manifest must exist");
	}

	/*--------------------------------------------------------------*/
	/*----------------       Table Collection       ----------------*/
	/*--------------------------------------------------------------*/

	/** Publish the completion summary only after every declared cell has closed. */
	public void process(){
		final File root=new File(input).getAbsoluteFile().getParentFile();
		require(new File(output).mkdir(), "Report output must be fresh");
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		final ByteBuilder prefix=new ByteBuilder();
		try(Input manifest=new Input(input);
			Output curves=new Output(path("curves.tsv"));
			Output matched=new Output(path("matched_rates.tsv"));
			Output ranking=new Output(path("ranking.tsv"));
			Output monitor=new Output(path("monitor_curves.tsv"));
			Output costs=new Output(path("costs.tsv"));
			Output index=new Output(path("index.tsv"))){
			RrnaResourceIO.header(manifest.reader, "cell\tarm\tend\tseed\toutput\trate\trole\tjob", input);
			costs.writer.println(PREFIX+"stage\tuser_cpu_seconds\tsystem_cpu_seconds\twall_seconds");
			monitor.writer.println(PREFIX+"step\tobservations\ttrain_mse_min\ttrain_mse_max\tmonitor_native_error_min\tmonitor_native_error_max\tmonitor_wer_min\tmonitor_wer_max\tfpr_min\tfpr_max\tfnr_min\tfnr_max\tlr_min\tlr_max\tinterpretation");
			index.writer.println(PREFIX+"job_id\tselected_model_seed\tdeclared_checkpoint_selection\tfit_source\tfit_build_seal\tmodel_bytes\tnative_directory\tdeclared_training_steps\tqualification");
			for(byte[] line=manifest.reader.nextLine(); line!=null; line=manifest.reader.nextLine()){
				p.set(line);
				require(p.terms()==8, "Every fit manifest row must retain its eight submitted fields");
				final String cell=p.parseString(0);
				final String arm=p.parseString(1);
				final String end=p.parseString(2);
				final int seed=p.parseInt(3);
				final String activation=p.parseString(4);
				final String rate=p.parseString(5);
				require(cell.matches("[A-Za-z0-9_.]+") && names.add(cell) && (arm.equals("R") || arm.equals("T"))
					&& (end.equals("5prime") || end.equals("3prime")) && seed>=0
					&& (activation.equals("sigmoid") || activation.equals("rslog"))
					&& Double.isFinite(Double.parseDouble(rate)) && Double.parseDouble(rate)>0,
					"Fit identities must be unique and retain their declared trainer, end, seed, activation and positive rate");
				final File run=new File(root, "runs/"+cell);
				require(one(new File(run, "exit.txt")).equals("0")
					&& one(new File(run, "complete.txt")).equals("GATE_CONVERGENCE_FIT_AND_SAVED_REPLAY_COMPLETE_OUTER_UNPERFORMED"),
					"Only closed fit/score/replay cells may enter the report: "+cell);
				require(one(new File(run, "source.txt")).equals(source) && one(new File(run, "build.seal")).equals(seal),
					"Every compared fit must use the declared source and qualified build: "+cell);
				final String job=one(new File(run, "job.txt"));
				require(job.matches("[0-9]+") && job.equals(one(new File(root, "submitted/"+cell+"/job.txt"))),
					"Executed and submitted job receipts must agree: "+cell);
				prefix.clear().append(cell).tab().append(arm).tab().append(end).tab().append(seed).tab().append(activation).tab().append(rate).tab();
				final long trials=end.equals("5prime") ? 88966 : 88979;
				final long proposals=end.equals("5prime") ? 50183 : 61986;
				merge(new File(run, "replay/summary.tsv"), curves, prefix, 0, end, trials, proposals);
				merge(new File(run, "replay/matched_rates.tsv"), matched, prefix, 1, end, trials, proposals);
				merge(new File(run, "replay/ranking.tsv"), ranking, prefix, 2, end, trials, proposals);
				final int trainingSteps=arm.equals("R") ? 800 : trainerBatches;
				trainingLog(new File(run, "train.stderr"), monitor, prefix, arm, end, trainingSteps);
				for(final String stage:new String[]{"train", "score", "replay"}){
					final String time=one(new File(run, stage+".time.tsv"));
					final String[] values=time.split("\t", -1);
					require(values.length==3, "Require separate native user/system/wall values for "+cell+"/"+stage);
					for(final String value:values){
						final double number=Double.parseDouble(value);
						require(Double.isFinite(number) && number>=0, "Native times must be finite and nonnegative");
					}
					costs.writer.print(new ByteBuilder().append(prefix).append(stage).tab().append(time).nl());
				}
				final File model=new File(run, "gate.bbnet");
				require(model.isFile() && model.length()>0, "Each completed cell must retain its selected model");
				final long selectedSeed=modelSeed(model);
				require(arm.equals("R") ? selectedSeed==seed : selectedSeed>=58 && selectedSeed<=65,
					"The selected network must come from this cell's declared fresh starts");
				index.writer.print(new ByteBuilder().append(prefix).append(job).tab().append(selectedSeed)
					.tab().append(arm.equals("R") ? "monitor_mse" : trainerSelection).tab().append(source).tab().append(seal)
					.tab().append(model.length()).tab().append(run.getPath()).tab().append(trainingSteps)
					.append("\tINTERNAL_MONITOR_SELECTION_OUTER_UNPERFORMED\n"));
				cells++;
			}
		}
		require(cells==expected, "The completed report must contain every declared cell: expected="+expected+" actual="+cells);
		try(Output summary=new Output(path("summary.tsv"))){
			summary.writer.println("cells\tcurve_rows\tmatched_rows\tranking_rows\tqualification");
			summary.writer.print(new ByteBuilder().append(cells).tab().append(rows[0]).tab().append(rows[1]).tab().append(rows[2])
				.append("\tLOSSLESS_ENDPOINT_TABLES_MONITOR_ENVELOPES_NO_SELECTION\n"));
		}
	}

	/** Validate table framing/populations while copying original data bytes exactly. */
	private void merge(final File file, final Output out, final ByteBuilder prefix, final int kind,
			final String end, final long trials, final long proposals){
		final LineParser1 p=new LineParser1('\t');
		int count=0;
		try(Input input=new Input(file.getPath())){
			final byte[] header=input.reader.nextLine();
			require(header!=null, "Completed table needs its schema: "+file);
			if(headers[kind]==null){
				headers[kind]=header.clone();
				out.writer.print(new ByteBuilder().append(PREFIX).append(header).nl());
			}else{require(Arrays.equals(headers[kind], header), "Compared cells must use identical report schemas: "+file);}
			p.set(header);
			final int width=p.terms();
			require(width==(kind==0 ? 40 : kind==1 ? 46 : 7), "Require the complete qualified gate report schema: "+file);
			int endColumn=-1, trialColumn=-1, proposalColumn=-1;
			for(int i=0; i<width; i++){
				if(p.termEquals("end", i)){endColumn=i;}
				if(p.termEquals("trials", i)){trialColumn=i;}
				if(p.termEquals("proposals", i)){proposalColumn=i;}
			}
			require(endColumn>=0 && proposalColumn>=0 && (kind==2 || trialColumn>=0), "Missing population fields in "+file);
			final ByteBuilder text=new ByteBuilder();
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				require(p.terms()==width && p.termEquals(end, endColumn) && p.parseLong(proposalColumn)==proposals
					&& (trialColumn<0 || p.parseLong(trialColumn)==trials), "Incomplete or mismatched population in "+file);
				if(kind==0 && count<4){
					final int endIndex=end.equals("5prime") ? 0 : 1;
					require(p.termEquals(CONTROLS[count], 0), "Fixed-rule controls must retain their declared order");
					if(baselines[endIndex][count]==null){baselines[endIndex][count]=line.clone();}
					else{require(Arrays.equals(baselines[endIndex][count], line), "Same-end fixed controls must match exactly across every fit: "+file);}
				}
				out.writer.print(text.clear().append(prefix).append(line).nl());
				count++;
			}
		}
		require(kind==0 ? count==26 : kind==1 ? count>=3 && count<=5 : count==1,
			"Require all coarse policies, every attainable matched-count bound and one AUC row: "+file+" rows="+count);
		rows[kind]+=count;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Monitor Evidence       ----------------*/
	/*--------------------------------------------------------------*/

	/** Read the entire trainer log. Trainer does not label individual seed lines:
	 * its eight observations are a per-setting envelope, never a claimed seed curve.
	 * Native Trainer ERR and RegressionTrainer MSE are deliberately not equated.
	 */
	private static void trainingLog(final File file, final Output out, final ByteBuilder prefix, final String arm, final String end, final int finalStep){
		final boolean regression=arm.equals("R");
		final LogPoint[] points=new LogPoint[finalStep+1];
		final long teaching=end.equals("5prime") ? 196595 : 245638;
		final long validation=end.equals("5prime") ? 50183 : 61986;
		boolean loaded=false, workers=regression, roundtrip=!regression;
		try(Input input=new Input(file.getPath())){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				final String text=new String(line, java.nio.charset.StandardCharsets.UTF_8);
				require(!text.contains("Exception in thread") && !text.contains("AssertionError") && !text.startsWith("Unknown parameter"),
					"A completed training log contains an error: "+file);
				if(regression){
					loaded|=text.equals("loaded "+teaching+" training samples, "+validation+" external validation samples, dims=[32, 32, 32, 16, 1]");
					roundtrip|=text.startsWith("round-trip check vs in-memory model:") && text.endsWith("(OK)");
				}else{
					loaded|=text.equals("Loaded samples: training="+teaching+" validation="+validation);
					workers|=text.startsWith("Using 64 worker threads, 8 trainer threads,") && text.endsWith("to train 8 networks total.");
					if(text.startsWith("Using 64 worker threads,")){
						require(text.contains(", "+finalStep+" batches, "), "Native Trainer budget must equal declared trainerbatches: "+file);
					}
				}
				if(!text.startsWith(regression ? "epoch " : "Batch ")){continue;}
				final Matcher matcher=(regression ? REGRESSION_LINE : TRAINER_LINE).matcher(text);
				require(matcher.matches(), "Unrecognized monitor line; do not silently drop an observation: "+text);
				final int step=Parse.parseIntKMG(matcher.group(1));
				require(step>0 && step<=finalStep, "Monitor step outside the declared training budget");
				LogPoint point=points[step];
				if(point==null){points[step]=point=new LogPoint();}
				point.count++;
				if(regression){
					point.add(5, matcher.group(2));
					point.add(0, matcher.group(3));
					point.add(1, matcher.group(4));
				}else{
					point.add(1, matcher.group(2));
					point.add(2, matcher.group(3));
					point.add(3, matcher.group(4));
					point.add(4, matcher.group(5));
					point.add(5, matcher.group(7));
				}
			}
		}
		require(loaded && workers && roundtrip && points[finalStep]!=null, "Require actual populations, worker contract, roundtrip and final training step: "+file);
		final ByteBuilder text=new ByteBuilder();
		for(int step=1; step<points.length; step++){
			final LogPoint point=points[step];
			if(point==null){continue;}
			require(point.count==(regression ? 1 : 8), "Each printed step must contain every trained network once: "+file+" step="+step+" count="+point.count);
			text.clear().append(prefix).append(step).tab().append(point.count);
			for(int metric=0; metric<6; metric++){
				text.tab();
				final boolean supplied=regression ? metric==0 || metric==1 || metric==5 : metric!=0;
				if(supplied){text.appendSlow(point.min[metric]).tab().appendSlow(point.max[metric]);}
				else{text.tab();}
			}
			out.writer.print(text.tab().append(regression ? "SINGLE_FRESH_MODEL_MONITOR_MSE" : "EIGHT_START_SETTING_ENVELOPE_NATIVE_ERR_NOT_SEED_IDENTIFIED").nl());
		}
	}

	/** Primitive bounds for one printed step; observations have no seed identity. */
	private static final class LogPoint{
		void add(final int index, final String text){
			assert(count>0 && index>=0 && index<6) : "trainingLog counts each parsed observation before adding its trainer-defined metrics";
			final double value=Double.parseDouble(text);
			require(Double.isFinite(value), "Nonfinite monitor metric: "+text);
			min[index]=count==1 ? value : Math.min(min[index], value);
			max[index]=count==1 ? value : Math.max(max[index], value);
		}
		final double[] min=new double[6], max=new double[6];
		int count;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Metadata Helpers       ----------------*/
	/*--------------------------------------------------------------*/

	/** Read only the serialized identity/layout; do not instantiate or infer a net. */
	private static long modelSeed(final File file){
		long seed=-1;
		try(Input input=new Input(file.getPath())){
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				final String text=new String(line, java.nio.charset.StandardCharsets.UTF_8);
				if(text.startsWith("#seed ")){seed=Long.parseLong(text.substring(6));}
				if(text.startsWith("#dims ")){
					require(text.equals("#dims 32 32 32 16 1") && seed>=0, "Selected model must retain the declared gate layout and seed");
					return seed;
				}
			}
		}
		throw new IllegalArgumentException("Missing selected network identity/layout: "+file);
	}

	/** Reject missing, empty or multi-row completion metadata instead of guessing. */
	private static String one(final File file){
		try(Input input=new Input(file.getPath())){
			final byte[] line=input.reader.nextLine();
			require(line!=null && input.reader.nextLine()==null, "Expected one complete metadata row: "+file);
			return new String(line, java.nio.charset.StandardCharsets.UTF_8);
		}
	}
	private String path(final String name){return new File(output, name).getPath();}
	private static void require(final boolean good, final String why){if(!good){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private static final String PREFIX="cell\ttrainer\tdeclared_end\tinitial_seed\toutput_activation\tlearning_rate\t";
	private static final String[] CONTROLS={"do_nothing", "margin0", "margin0.1", "all_proposals"};
	private static final Pattern REGRESSION_LINE=Pattern.compile("epoch ([0-9]+) lr=(\\S+) trainMSE=(\\S+) valMSE=(\\S+) best=(\\S+)");
	private static final Pattern TRAINER_LINE=Pattern.compile("Batch ([0-9]+k?):\\s+ERR=\\s*(\\S+)\\s+WER=\\s*(\\S+)\\s+FPR=\\s*(\\S+)\\s+FNR=\\s*(\\S+)\\s+CTF=\\s*(\\S+)\\s+ALP=\\s*(\\S+)\\s+\\S+s");
	/** First encountered fixed-control rows, independently retained for each end. */
	private final byte[][][] baselines=new byte[2][4][];
	/** First encountered schemas and copied row counts for the three table kinds. */
	private final byte[][] headers=new byte[3][];
	private final long[] rows=new long[3];
	private String input, output, source, seal;
	/** Explicit caller declaration; a fixed cutoff alone does not choose WER ranking. */
	private String trainerSelection;
	/** Explicit uniform budget for Trainer cells in this report; never inferred from an incomplete log. */
	private int trainerBatches=150000;
	private int expected=-1, cells;
}
