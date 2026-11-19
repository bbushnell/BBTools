package prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ReadWrite;
import ml.CellNet;
import ml.CellNetParser;
import ml.Trainer;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.FloatList;

/** Combines unchanged gate rows and binds training-only input scaling to fresh
 * random networks. CellNet applies the same serialized transform in training
 * and inference; no normalized matrix or folded-weight approximation is needed.
 * This prepares eight distinct untrained starts, not eight retrains of one net.
 * @author Brian, Ganyu
 */
public final class NcrnaGateTrainingInputs{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaGateTrainingInputs.class, false);
		try{
			new NcrnaGateTrainingInputs(pp.args).process();
			pp.outstream.println("NCRNA_GATE_TRAINING_INPUTS_COMPLETE_EIGHT_FRESH_STARTS_NO_FIT");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public NcrnaGateTrainingInputs(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("monitor")){monitor=b;}
			else if(a.equals("expected")){
				final String[] parts=b.split(",");
				expected=new long[parts.length];
				for(int i=0; i<parts.length; i++){expected[i]=Long.parseLong(parts[i]);}
			}else if(a.equals("monitorexpected")){monitorExpected=Long.parseLong(b);}
			else if(a.equals("seed")){seed=Long.parseLong(b);}
			else{require(parser.parse(arg, a, b), "Unknown gate-input argument: "+arg);}
		}
		training=parser.in1==null ? null : parser.in1.split(",");
		output=parser.out1;
		require(training!=null && training.length>0 && expected!=null && expected.length==training.length
			&& monitor!=null && monitorExpected>0 && seed>=0 && seed<=Long.MAX_VALUE-7 && output!=null,
			"Bind in=training_matrices expected=per_file_counts monitor= monitorexpected= seed= out=fresh_directory");
		for(final long count:expected){require(count>0, "Each accepted training matrix must be nonempty");}
		require(!parser.overwrite && !parser.append && parser.in2==null && parser.out2==null
			&& parser.maxReads<0 && parser.samplerate==1 && !new File(output).exists(), "Prepare complete populations into a fresh directory");
		final String[] inputs=Arrays.copyOf(training, training.length+1);
		inputs[training.length]=monitor;
		require(Tools.testInputFiles(false, true, inputs), "Training/monitor matrices must be readable and distinct");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Data Preparation      ----------------*/
	/*--------------------------------------------------------------*/

	public void process(){
		require(new File(output).mkdir(), "Cannot create fresh gate training-input directory");
		Arrays.fill(minimums, Float.POSITIVE_INFINITY);
		Arrays.fill(maximums, Float.NEGATIVE_INFINITY);
		for(int i=0; i<32; i++){
			trainingColumns[i]=new FloatList();
			monitorColumns[i]=new FloatList();
		}
		final boolean oldBgzf=ReadWrite.USE_BGZF;
		final boolean oldForce=ReadWrite.FORCE_BGZIP;
		ReadWrite.USE_BGZF=ReadWrite.FORCE_BGZIP=true;
		final String combined=path("train.raw.tsv.gz");
		try(Output raw=new Output(combined)){
			raw.writer.println(HEADER);
			final LineParser1 p=new LineParser1('\t');
			final float[] values=new float[32];
			final ByteBuilder text=new ByteBuilder();
			for(int file=0; file<training.length; file++){
				long count=0;
				try(Input input=new Input(training[file])){
					RrnaResourceIO.header(input.reader, HEADER, training[file]);
					for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
						final int label=parse(line, p, values);
						count++;
						rows++;
						if(rows==1){System.arraycopy(values, 0, probe, 0, 32);}
						positives+=label;
						for(int i=0; i<32; i++){
							minimums[i]=Math.min(minimums[i], values[i]);
							maximums[i]=Math.max(maximums[i], values[i]);
							trainingColumns[i].add(values[i]);
						}
						raw.writer.print(text.clear().append(line).nl());
					}
				}
				require(count==expected[file], "Training input count differs: "+training[file]+" expected="+expected[file]+" actual="+count);
			}
		}finally{ReadWrite.USE_BGZF=oldBgzf;ReadWrite.FORCE_BGZIP=oldForce;}
		for(int i=0; i<32; i++){
			final double range=(double)maximums[i]-minimums[i];
			means[i]=minimums[i];
			inverse[i]=(float)(range==0 ? 1 : 1/range);
			require(Float.isFinite(means[i]) && Float.isFinite(inverse[i]) && inverse[i]>0, "Input scaling must be finite");
			// Positive scaling is monotone: finite transformed extrema cover every
			// training value already validated and emitted in the first pass.
			require(Float.isFinite((minimums[i]-means[i])*inverse[i]) && Float.isFinite((maximums[i]-means[i])*inverse[i]),
				"Normalized training extrema overflow at column "+i);
		}
		writeRanges(trainingColumns, "training_feature_ranges.tsv", "TRAINING_FOLDS_ONLY", null);
		validateMonitor();
		writeRanges(monitorColumns, "monitor_feature_ranges.tsv", "MONITOR_NOT_USED_TO_FIT_SCALING", outsideMonitorColumns);
		for(int i=0; i<32; i++){trainingColumns[i].addAll(monitorColumns[i]);}
		writeRanges(trainingColumns, "all_feature_ranges.tsv", "ALL_FOLDS_DIAGNOSTIC_ONLY", outsideMonitorColumns);
		writeScaling();
		writeFreshNetworks(combined);
		try(Output report=new Output(path("summary.tsv"))){
			report.writer.println("training_rows\ttraining_improving\tmonitor_rows\tmonitor_improving\tinputs\trandom_starts\tmonitor_rows_outside_training_range\tqualification");
			report.writer.print(new ByteBuilder().append(rows).tab().append(positives).tab().append(monitorExpected)
				.tab().append(monitorPositives).append("\t32\t8\t").append(outsideMonitor)
				.append("\tTRAINING_MINMAX_FLOAT32_EIGHT_DISTINCT_UNTRAINED_SEEDS_NO_CLIPPING\n"));
		}
	}

	/** The monitor is checked only after scaling is frozen; its values cannot fit it. */
	private void validateMonitor(){
		final LineParser1 p=new LineParser1('\t');
		final float[] values=new float[32];
		long count=0;
		long positive=0;
		try(Input input=new Input(monitor)){
			RrnaResourceIO.header(input.reader, HEADER, monitor);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				positive+=parse(line, p, values);
				boolean outside=false;
				for(int i=0; i<32; i++){
					final float normalized=(values[i]-means[i])*inverse[i];
					require(Float.isFinite(normalized), "Normalized monitor input overflow at column "+i);
					monitorColumns[i].add(values[i]);
					if(values[i]<minimums[i] || values[i]>maximums[i]){
						outsideMonitorColumns[i]++;
						outside=true;
					}
				}
				if(outside){outsideMonitor++;}
				count++;
			}
		}
		require(count==monitorExpected, "Complete monitor count differs for "+monitor);
		monitorPositives=positive;
	}

	private void writeScaling(){
		try(Output report=new Output(path("normalization.tsv"))){
			report.writer.println("column\ttraining_rows\toffset_float_bits\tinverse_range_float_bits\toffset\tinverse_range\traw_min\traw_max\tnormalized_min\tnormalized_max");
			final ByteBuilder text=new ByteBuilder();
			for(int i=0; i<32; i++){
				report.writer.print(text.clear().append(i).tab().append(rows).tab().append(Float.floatToRawIntBits(means[i]))
					.tab().append(Float.floatToRawIntBits(inverse[i])).tab().appendSlow(means[i]).tab().appendSlow(inverse[i])
					.tab().appendSlow(minimums[i]).tab().appendSlow(maximums[i])
					.tab().appendSlow((minimums[i]-means[i])*inverse[i]).tab().appendSlow((maximums[i]-means[i])*inverse[i]).nl());
			}
		}
	}

	/** Exact per-column order statistics; diagnostic monitor values never alter scaling. */
	private void writeRanges(final FloatList[] columns, final String file, final String role, final long[] outside){
		try(Output report=new Output(path(file))){
			report.writer.println("role\tgroup\tfeature\tcolumn\trows\tmin\tp01\tp50\tp99\tmax\toutside_training_range");
			final ByteBuilder text=new ByteBuilder();
			for(int i=0; i<32; i++){
				final FloatList values=columns[i];
				require(values.size()>0, "Input range diagnostics require a nonempty declared population");
				values.sort();
				final String group=i<22 ? "first_stage_scores" : i<24 ? "length_ratios" : i<28 ? "gc" : i==28 ? "identity" : "fuzziness";
				final String feature=i<22 ? "score"+i : FEATURES[i-22];
				text.clear().append(role).tab().append(group).tab().append(feature).tab().append(i).tab().append(values.size());
				for(final double fraction:new double[]{0, .01, .5, .99, 1}){
					final int index=Math.max(0, (int)Math.ceil(fraction*values.size())-1);
					text.tab().appendSlow(values.get(index));
				}
				report.writer.print(text.tab().append(outside==null ? 0 : outside[i]).nl());
			}
		}
	}

	/** Constructs eight independent seeds without invoking Trainer.process. */
	private void writeFreshNetworks(final String combined){
		final Trainer initializer=new Trainer(new String[]{"in="+combined, "dims=32,32,32,16,1", "batches=1",
			"final=rslog", "sig=0.6", "tanh=0.4", "rslog=0.02", "msig=0.02", "density=1", "cutoffeval=0.5", "t=1"});
		try(Output manifest=new Output(path("initial_networks.tsv"))){
			manifest.writer.println("seed\tnetwork\tstate");
			for(int i=0; i<8; i++){
				final CellNet net=initializer.randomNetwork(seed+i);
				net.setInputNormalization(means, inverse);
				net.setTag("initialization", "fresh-random-untrained");
				net.setTag("normalization-source", path("normalization.tsv"));
				net.setTag("normalization-method", "training-minmax-no-clipping");
				final String file=path("seed"+(seed+i)+".bbnet");
				try(Output writer=new Output(file)){writer.writer.print(net.toBytes());}
				final CellNet loaded=CellNetParser.load(file, false);
				require(loaded.hasInputNormalization() && Arrays.equals(loaded.inputMeanCopy(), means)
					&& Arrays.equals(loaded.inputInverseStdCopy(), inverse), "Fresh network normalization must survive exact serialization");
				checkRoundTrip(net, loaded);
				manifest.writer.print(new ByteBuilder().append(seed+i).tab().append(file).append("\tUNTRAINED_RANDOM_START\n"));
			}
		}
	}

	private void checkRoundTrip(final CellNet before, final CellNet after){
		final float[] values=means.clone();
		for(int i=0; i<32; i++){
			values[i]=probe[i];
			before.applyInput(values);
			after.applyInput(values);
			before.feedForward();
			after.feedForward();
			require(Float.isFinite(before.getOutput(0)) && Float.floatToRawIntBits(before.getOutput(0))==Float.floatToRawIntBits(after.getOutput(0)),
				"Fresh normalized network must preserve raw-input predictions through export");
			values[i]=means[i];
		}
	}

	private static int parse(final byte[] line, final LineParser1 p, final float[] values){
		p.set(line);
		require(p.terms()==33, "Gate matrices require32inputs and one binary label, without weights");
		for(int i=0; i<32; i++){
			require(p.length(i)>0, "Empty gate input field");
			p.setBounds(i);
			boolean exponent=false;
			for(int j=p.a(); j<p.b(); j++){exponent|=p.line()[j]=='e' || p.line()[j]=='E';}
			values[i]=exponent ? (float)Parse.parseDoubleSlow(p.line(), p.a(), p.b()) : p.parseFloat(i);
			require(Float.isFinite(values[i]), "Nonfinite gate input");
		}
		require(p.termEquals("0", 32) || p.termEquals("1", 32), "Gate target must remain a strict binary label");
		return p.termEquals("1", 32) ? 1 : 0;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private String path(final String name){return new File(output, name).getPath();}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private static final String HEADER="#dims\t32\t1\t0";
	private String[] training;
	private String monitor, output;
	private long[] expected;
	private long monitorExpected=-1, seed=-1, rows, positives, monitorPositives, outsideMonitor;
	private final float[] means=new float[32], inverse=new float[32];
	private final float[] probe=new float[32];
	private final float[] minimums=new float[32], maximums=new float[32];
	private final FloatList[] trainingColumns=new FloatList[32], monitorColumns=new FloatList[32];
	private final long[] outsideMonitorColumns=new long[32];
	private static final String[] FEATURES={"rawLengthRatio", "proposedLengthRatio", "contigGC", "rawSpanGC", "referenceGC",
		"spanMinusReferenceGC", "rescaledRawIdentity", "fuzzMinus", "fuzzAt", "fuzzPlus"};
}
