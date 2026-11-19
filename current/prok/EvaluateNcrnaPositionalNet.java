package prok;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import fileIO.FileFormat;
import fileIO.ByteStreamWriter;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Scores one exported positional network on held-out jittered endpoint trials.
 * Reads the numeric matrix and its identity sidecar together. This conditional
 * vector assay measures neither unique-genome recall nor joint-end caller harm.
 * Adapted from template.A_SampleByteFile; all streams have explicit ownership.
 * @author Brian Bushnell, Raiden
 */
public final class EvaluateNcrnaPositionalNet{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses native arguments, scores the complete input and closes diagnostics. */
	public static void main(final String[] args){
		final Timer timer=new Timer();
		final EvaluateNcrnaPositionalNet tool=new EvaluateNcrnaPositionalNet(args);
		try{
			if(tool.replayFull){tool.processReplay(timer);}
			else if(tool.fullErrors){tool.processFull(timer);}
			else{tool.process(timer);}
		}finally{
			Shared.closeStream(tool.outstream);
		}
	}

	/** Accepts one matrix, paired rows, network, endpoint and summary destination. */
	public EvaluateNcrnaPositionalNet(String[] args){
		final PreParser pre=new PreParser(args, null, false);
		args=pre.args;
		outstream=pre.outstream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase();
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("rows")){
				rows=b;
			}else if(a.equals("raw")){
				rawPath=b;
			}else if(a.equals("comparison")){
				comparisonPath=b;
			}else if(a.equals("net")){
				netPath=b;
			}else if(a.equals("fullerrors")){
				fullErrors=Parse.parseBoolean(b);
			}else if(a.equals("featurevariant")){
				featureVariant=b;
				fullInputs=fullInputs(b);
			}else if(a.equals("extended")){
				extended=Parse.parseBoolean(b);
			}else if(a.equals("replayfull")){
				replayFull=Parse.parseBoolean(b);
			}else if(a.equals("end")){
				if("5prime".equals(b)){end=0;}
				else if("3prime".equals(b)){end=1;}
				else{throw new IllegalArgumentException("end= must be 5prime or 3prime");}
			}else if(parser.parse(arg, a, b)){
				// Native input/output, overwrite and shared options.
			}else{
				throw new IllegalArgumentException("Unknown parameter "+arg);
			}
		}
		in=Tools.fixExtension(parser.in1);
		out=parser.out1;
		if(!replayFull && rawPath==null && out!=null){rawPath=out+".scores.tsv";}
		if(!fullErrors && (replayFull || extended || !featureVariant.equals("baseline"))){
			throw new IllegalArgumentException("Feature variants, extended metrics and saved22 replay require fullerrors=t");
		}
		if(replayFull && (netPath!=null || rows!=null || rawPath!=null || !featureVariant.equals("baseline"))){
			throw new IllegalArgumentException("Replay reads one complete saved22score file at in=; no net/rows/raw/feature override is consumed");
		}
		if(fullErrors && comparisonPath!=null){throw new IllegalArgumentException("fullerrors=t writes its complete policy table at out=; comparison= belongs to the 21-output path");}
		if(!fullErrors && comparisonPath==null && out!=null){comparisonPath=out+".comparison.tsv";}
		overwrite=parser.overwrite;
		if(in==null || (!replayFull && (rows==null || netPath==null)) || out==null || end<0 || parser.append
				|| parser.in2!=null || parser.out2!=null || parser.maxReads>=0 || parser.samplerate!=1){
			throw new IllegalArgumentException("Require in=matrix rows=sidecar net=export end=5prime|3prime out=summary; append is unsupported");
		}
		if(!Tools.testInputFiles(false, true, in, rows, netPath)
				|| !Tools.testOutputFiles(overwrite, false, false, out, rawPath, comparisonPath)
				|| !Tools.testForDuplicateFiles(true, in, rows, netPath, out, rawPath, comparisonPath)){
			throw new IllegalArgumentException("Input/output paths must be readable, writable and distinct");
		}
		assert(end==0 || end==1) : "Endpoint controls which saved sense-coordinate pair supplies the label";
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes a summary only after both complete input streams close successfully. */
	private void process(final Timer timer){
		final CellNet net=CellNetParser.load(netPath, false);
		final Comparison comparison;
		final FileFormat rawFormat=FileFormat.testOutput(rawPath, FileFormat.TXT, null, true, overwrite, false, false);
		try(final NcrnaTextIO.Output raw=new NcrnaTextIO.Output(rawFormat)){
			final ByteBuilder provenance=new ByteBuilder("#NCRNA_POSITIONAL_SCORES_V1\n#coding\tFP32_RAW_SIGNED_INT\n");
			provenance.append("#matrix\t").append(in).nl().append("#rows\t").append(rows).nl()
				.append("#network\t").append(netPath).nl().append("#end\t").append(end==0 ? "5prime" : "3prime").nl();
			raw.writer.print(provenance);
			comparison=compare(net, in, rows, end, raw.writer);
		}
		final Result result=comparison.ungated;
		final FileFormat format=FileFormat.testOutput(out, FileFormat.TXT, null, true, overwrite, false, false);
		try(final NcrnaTextIO.Output output=new NcrnaTextIO.Output(format)){
			output.writer.print(result.summary(end).toBytes());
		}
		final FileFormat compared=FileFormat.testOutput(comparisonPath, FileFormat.TXT, null, true, overwrite, false, false);
		try(final NcrnaTextIO.Output output=new NcrnaTextIO.Output(compared)){
			output.writer.print(comparison.summary(end));
		}
		timer.stop();
		outstream.println(Tools.thingsProcessed(timer.elapsed, result.trials, 8, "Trials"));
		outstream.println("Usable jittered trials; caller detection and joint-end validity are not measured.");
	}

	/** Scores one declared end. The borrowed net is mutated and must be thread confined.
	 * Sidecar labels/counts are checked; numeric-only rows carry no identities with
	 * which to detect a same-label sidecar permutation. Their pairing comes from
	 * the accepted producer and input manifest, not from this reader.
	 */
	public static Result evaluate(final CellNet net, final String matrix, final String sidecar, final int end){
		return compare(net, matrix, sidecar, end, null).ungated;
	}

	/** Evaluates each trial once, retaining all scores before applying either margin.
	 * The optional raw writer is borrowed. A COMPLETE footer is emitted only after
	 * both readers close successfully; callers must also check writer closure.
	 * Input/network fingerprints and execution flags belong in the run manifest.
	 */
	public static Comparison compare(final CellNet net, final String matrix, final String sidecar,
			final int end, final ByteStreamWriter raw){
		if(net==null || (net.numInputs()!=24 && net.numInputs()!=69) || net.numOutputs()!=CLASSES || end<0 || end>1){
			throw new IllegalArgumentException("Require a 24- or 69-input/21-output net and one declared endpoint");
		}
		final Comparison comparison=new Comparison();
		final Result result=comparison.ungated;
		final LineParser1 values=new LineParser1('\t'), metadata=new LineParser1('\t');
		final int width=net.numInputs();
		final float[] input=new float[width], scores=new float[CLASSES];
		final ByteBuilder buffer=new ByteBuilder();
		try(final NcrnaTextIO.Input data=new NcrnaTextIO.Input(matrix);
				final NcrnaTextIO.Input identities=new NcrnaTextIO.Input(sidecar)){
			header(data.reader.nextLine(), matrixHeader(width), matrix);
			header(identities.reader.nextLine(), ROW_HEADER, sidecar);
			if(raw!=null){
				buffer.append(ROW_HEADER).append("\ttrue_class\tmargin0_class\tmargin01_class");
				for(int i=0; i<CLASSES; i++){buffer.append("\tscore_bits_").append(i);}
				raw.print(buffer.nl());
			}
			for(byte[] line=data.reader.nextLine(); line!=null; line=data.reader.nextLine()){
				final byte[] row=identities.reader.nextLine();
				if(row==null){throw new IllegalArgumentException("Sidecar ended before matrix trial "+(result.trials+1));}
				values.set(line);
				final int truth=parseVector(values, input, result.trials+1);
				metadata.set(row);
				validateIdentity(metadata, truth, end, result.trials+1);
				net.applyInput(input);// Serialized normalization applies exactly once here.
				net.feedForward();
				for(int i=0; i<CLASSES; i++){scores[i]=net.getOutput(i);}
				comparison.accept(scores, truth);
				if(raw!=null){
					buffer.clear().append(row).tab().append(truth).tab().append(comparison.predicted)
						.tab().append(comparison.gatedPrediction);
					for(final float score:scores){buffer.tab().append(Float.floatToRawIntBits(score));}
					raw.print(buffer.nl());
				}
			}
			if(identities.reader.nextLine()!=null){throw new IllegalArgumentException("Sidecar has more rows than the matrix");}
		}
		if(result.trials==0){throw new IllegalArgumentException("No usable held-out trials");}
		assert(result.exact<=result.within1 && result.within1<=result.within3 && result.within3<=result.trials)
			: "Nested endpoint tolerances must conserve the same complete trial population";
		if(raw!=null){raw.print(buffer.clear().append("#COMPLETE\t").append(result.trials).nl());}
		return comparison;
	}

	/** The input width is part of the saved matrix schema, never inferred per row. */
	static byte[] matrixHeader(final int width){
		return ("#dims\t"+width+"\t21\t0").getBytes(StandardCharsets.US_ASCII);
	}

	/** Validates the exact producer schema, without allocating strings per row. */
	static void header(final byte[] actual, final byte[] expected, final String path){
		if(!Arrays.equals(actual, expected)){throw new IllegalArgumentException("Unexpected header: "+path);}
	}

	/** Parses finite raw features and exactly one unit label at the declared width. */
	static int parseVector(final LineParser1 parser, final float[] input, final long row){
		final int width=input.length;
		if(parser.terms()!=width+CLASSES){throw new IllegalArgumentException("Expected "+(width+CLASSES)+" numeric fields at trial "+row);}
		int truth=-1;
		for(int i=0; i<width+CLASSES; i++){
			if(parser.length(i)==0){throw new IllegalArgumentException("Empty numeric field at trial "+row);}
			final float value=numericFloat(parser, i);
			if(!Float.isFinite(value)){throw new IllegalArgumentException("Nonfinite numeric field at trial "+row);}
			if(i<width){
				input[i]=value;
			}else if(value==1){
				if(truth>=0){throw new IllegalArgumentException("Multiple unit targets at trial "+row);}
				truth=i-width;
			}else if(value!=0){
				throw new IllegalArgumentException("Target must be one-hot at trial "+row);
			}
		}
		if(truth<0){throw new IllegalArgumentException("Missing unit target at trial "+row);}
		assert(truth<CLASSES) : "The output index comes only from the declared 21 target fields";
		return truth;
	}

	/** Float.toString may use scientific notation for valid small table scores.
	 * Keep the existing native path for ordinary tokens; only exponent tokens use
	 * the native slow fallback, as RrnaPositionalKmerTable.savedDouble already does.
	 */
	private static float numericFloat(final LineParser1 parser, final int term){
		parser.setBounds(term);
		final byte[] line=parser.line();
		for(int i=parser.a(); i<parser.b(); i++){
			if(line[i]=='E' || line[i]=='e'){return (float)Parse.parseDoubleSlow(line, parser.a(), parser.b());}
		}
		return parser.parseFloat(term);
	}

	/** Checks evaluation provenance and zero-based inclusive sense-coordinate labels. */
	static void validateIdentity(final LineParser1 parser, final int truth, final int end, final long row){
		if(parser.terms()!=11 || !parser.termEquals("evaluation", 3)){
			throw new IllegalArgumentException("Require an eleven-field evaluation sidecar at trial "+row);
		}
		for(int i=0; i<11; i++){
			if(parser.length(i)==0){throw new IllegalArgumentException("Empty sidecar field at trial "+row);}
		}
		final int model=parser.parseInt(4), sample=parser.parseInt(6);
		final int rawStart=parser.parseInt(7), rawStop=parser.parseInt(8);
		final int trueStart=parser.parseInt(9), trueStop=parser.parseInt(10);
		final long label=(end==0 ? (long)trueStart-rawStart : (long)trueStop-rawStop)+RADIUS;
		if(model<0 || sample<0 || rawStart<0 || rawStop<rawStart || trueStart<0 || trueStop<trueStart || label!=truth){
			throw new IllegalArgumentException("Sidecar coordinates/label disagree at trial "+row);
		}
	}

	/** The opt-in 22-answer route leaves legacy parsing, selection and files unchanged. */
	private void processFull(final Timer timer){
		final CellNet net=CellNetParser.load(netPath, false);
		final FullComparison result;
		final FileFormat format=FileFormat.testOutput(rawPath, FileFormat.TXT, null, true, overwrite, false, false);
		try(final NcrnaTextIO.Output raw=new NcrnaTextIO.Output(format)){
			raw.writer.print(new ByteBuilder("#NCRNA_FULL_ERROR_SCORES_V1\n#coding\tFP32_RAW_SIGNED_INT\n")
				.append("#matrix\t").append(in).nl().append("#rows\t").append(rows).nl()
				.append("#network\t").append(netPath).nl().append("#end\t").append(end==0 ? "5prime" : "3prime").nl());
			if(!featureVariant.equals("baseline")){
				raw.writer.print(new ByteBuilder("#featurevariant\t").append(featureVariant).nl());
			}
			result=compareFull(net, in, rows, end, raw.writer, fullInputs);
		}
		try(final NcrnaTextIO.Output output=new NcrnaTextIO.Output(FileFormat.testOutput(out,
				FileFormat.TXT, null, true, overwrite, false, false))){
			output.writer.print(result.summary(end, extended));
		}
		timer.stop();
		outstream.println(Tools.thingsProcessed(timer.elapsed, result.margin0.trials, 8, "Trials"));
		outstream.println("Full signed-error trials; all-row RMS includes unreachable truth. This is not whole-caller accuracy.");
	}

	/** Replays already saved outputs; no network is opened or evaluated again. */
	private void processReplay(final Timer timer){
		final FullComparison result=replayFull(in, end);
		try(final NcrnaTextIO.Output output=new NcrnaTextIO.Output(FileFormat.testOutput(out,
				FileFormat.TXT, null, true, overwrite, false, false))){
			output.writer.print(result.summary(end, extended));
		}
		timer.stop();
		outstream.println(Tools.thingsProcessed(timer.elapsed, result.margin0.trials, 8, "Saved trials"));
		outstream.println("NCRNA_FULL_ERROR_REPLAY_COMPLETE_NO_INFERENCE");
	}

	/** Validates complete provenance, every raw score and stored decision, and the
	 * terminal row count before publishing replay metrics. Truncated streams fail.
	 */
	public static FullComparison replayFull(final String path, final int end){
		if(path==null || end<0 || end>1){throw new IllegalArgumentException("Replay needs a saved22score path and an endpoint");}
		final FullComparison result=new FullComparison();
		final LineParser1 p=new LineParser1('\t');
		final LineParser1 metadata=new LineParser1('\t');
		final float[] scores=new float[22];
		try(final NcrnaTextIO.Input input=new NcrnaTextIO.Input(path)){
			RrnaResourceIO.header(input.reader, "#NCRNA_FULL_ERROR_SCORES_V1", path);
			RrnaResourceIO.header(input.reader, "#coding\tFP32_RAW_SIGNED_INT", path);
			for(final String field:new String[]{"#matrix", "#rows", "#network", "#end"}){
				final byte[] row=input.reader.nextLine();
				if(row==null){throw new IllegalArgumentException("Missing saved-score provenance: "+field);}
				p.set(row);
				if(p.terms()!=2 || !p.termEquals(field, 0) || p.length(1)==0
						|| field.equals("#end") && !p.termEquals(end==0 ? "5prime" : "3prime", 1)){
					throw new IllegalArgumentException("Saved-score provenance or endpoint mismatch: "+field);
				}
			}
			byte[] line=input.reader.nextLine();
			if(line!=null){
				p.set(line);
				if(p.termEquals("#featurevariant", 0)){
					if(p.terms()!=2){throw new IllegalArgumentException("Malformed saved feature variant");}
					fullInputs(p.parseString(1));
					line=input.reader.nextLine();
				}
			}
			header(line, fullRawHeader().toBytes(), path);
			boolean complete=false;
			for(line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				if(p.termEquals("#COMPLETE", 0)){
					if(p.terms()!=2 || p.parseLong(1)!=result.margin0.trials || result.margin0.trials==0
							|| input.reader.nextLine()!=null){throw new IllegalArgumentException("Saved-score completion count or terminal position is invalid");}
					complete=true;
					break;
				}
				if(p.terms()!=38){throw new IllegalArgumentException("Saved22score rows require 38 fields");}
				for(int i=0; i<38; i++){if(p.length(i)==0){throw new IllegalArgumentException("Empty saved-score field");}}
				final int truth=p.parseInt(11);
				if(truth<0 || truth>21){throw new IllegalArgumentException("Invalid saved22target");}
				metadata.set(line, 10);
				final long correction=fullCorrection(metadata, truth, end, result.margin0.trials+1);
				if(correction!=p.parseLong(12)){throw new IllegalArgumentException("Saved signed correction differs from coordinates");}
				for(int i=0; i<22; i++){scores[i]=Float.intBitsToFloat(p.parseInt(16+i));}
				result.accept(scores, correction);
				if(result.predicted0!=p.parseInt(13) || result.predicted01!=p.parseInt(14) || result.argmax22!=p.parseInt(15)){
					throw new IllegalArgumentException("Saved strict-policy choices disagree with the preserved scores");
				}
			}
			if(!complete){throw new IllegalArgumentException("Saved22score stream lacks a terminal completion marker");}
		}
		return result;
	}

	private static ByteBuilder fullRawHeader(){
		final ByteBuilder out=new ByteBuilder().append(ROW_HEADER)
			.append("\ttrue_class\ttrue_minus_raw_bp\tmargin0_class\tmargin01_class\targmax22_class");
		for(int i=0; i<22; i++){out.append("\tscore_bits_").append(i);}
		return out;
	}

	/** Explicit single-variable layouts; no inference from an arbitrary net width. */
	private static int fullInputs(final String variant){
		if("baseline".equals(variant) || "identity".equals(variant)){return 69;}
		if("length".equals(variant)){return 70;}
		if("gc".equals(variant)){return 71;}
		throw new IllegalArgumentException("featurevariant= must be baseline, length, gc or identity");
	}

	/** One inference per trial supplies both strict policies, the baseline and all raw scores. */
	public static FullComparison compareFull(final CellNet net, final String matrix, final String sidecar,
			final int end, final ByteStreamWriter raw){
		return compareFull(net, matrix, sidecar, end, raw, 69);
	}

	/** Opt-in feature widths preserve the same22targets, cohorts and strict rules. */
	public static FullComparison compareFull(final CellNet net, final String matrix, final String sidecar,
			final int end, final ByteStreamWriter raw, final int inputs){
		if(inputs<69 || inputs>71 || net==null || net.numInputs()!=inputs || net.numOutputs()!=22 || end<0 || end>1){
			throw new IllegalArgumentException("Full-error evaluation requires the explicit69/70/71-input layout,22outputs and one endpoint");
		}
		final FullComparison result=new FullComparison();
		final float[] input=new float[inputs], scores=new float[22];
		final LineParser1 values=new LineParser1('\t'), metadata=new LineParser1('\t');
		final ByteBuilder buffer=new ByteBuilder();
		try(final NcrnaTextIO.Input data=new NcrnaTextIO.Input(matrix);
				final NcrnaTextIO.Input identities=new NcrnaTextIO.Input(sidecar)){
			header(data.reader.nextLine(), ("#dims\t"+inputs+"\t22\t0").getBytes(StandardCharsets.US_ASCII), matrix);
			header(identities.reader.nextLine(), ROW_HEADER, sidecar);
			if(raw!=null){
				buffer.append(ROW_HEADER).append("\ttrue_class\ttrue_minus_raw_bp\tmargin0_class\tmargin01_class\targmax22_class");
				for(int i=0; i<22; i++){buffer.append("\tscore_bits_").append(i);}
				raw.print(buffer.nl());
			}
			for(byte[] line=data.reader.nextLine(); line!=null; line=data.reader.nextLine()){
				final long trial=result.margin0.trials+1;
				final byte[] row=identities.reader.nextLine();
				if(row==null){throw new IllegalArgumentException("Sidecar ended before full-error trial "+trial);}
				values.set(line);
				final int truth=parseFullVector(values, input, trial);
				metadata.set(row);
				final long correction=fullCorrection(metadata, truth, end, trial);
				net.applyInput(input);
				net.feedForward();
				for(int i=0; i<22; i++){scores[i]=net.getOutput(i);}
				result.accept(scores, correction);
				if(raw!=null){
					buffer.clear().append(row).tab().append(truth).tab().append(correction)
						.tab().append(result.predicted0).tab().append(result.predicted01).tab().append(result.argmax22);
					for(final float score:scores){buffer.tab().append(Float.floatToRawIntBits(score));}
					raw.print(buffer.nl());
				}
			}
			if(identities.reader.nextLine()!=null){throw new IllegalArgumentException("Sidecar has more rows than the full-error matrix");}
		}
		if(result.margin0.trials==0){throw new IllegalArgumentException("No usable full-error held-out trials");}
		assert(result.margin0.trials==result.margin01.trials && result.margin0.trials==result.doNothing.trials)
			: "A COMPLETE raw stream must account for the same rows under both policies and the baseline";
		if(raw!=null){raw.print(buffer.clear().append("#COMPLETE\t").append(result.margin0.trials).nl());}
		return result;
	}

	/** Full schema validation rejects weights, unavailable features and non-unit targets. */
	private static int parseFullVector(final LineParser1 parser, final float[] input, final long row){
		final int fields=input.length+22;
		if(parser.terms()!=fields){throw new IllegalArgumentException("Expected "+fields+" numeric fields at full-error trial "+row);}
		int truth=-1;
		for(int i=0; i<fields; i++){
			if(parser.length(i)==0){throw new IllegalArgumentException("Empty full-error numeric field at trial "+row);}
			final float value=numericFloat(parser, i);
			if(!Float.isFinite(value)){throw new IllegalArgumentException("Nonfinite full-error numeric field at trial "+row);}
			if(i<input.length){input[i]=value;}
			else if(value==1){
				if(truth>=0){throw new IllegalArgumentException("Multiple full-error unit targets at trial "+row);}
				truth=i-input.length;
			}else if(value!=0){throw new IllegalArgumentException("Full-error target must be one-hot at trial "+row);}
		}
		if(truth<0){throw new IllegalArgumentException("Missing full-error unit target at trial "+row);}
		assert(truth<22) : "The explicit schema has 21 reachable offsets and one non-positional answer";
		return truth;
	}

	/** Returns the real signed correction; class21 never substitutes a +11-base distance. */
	private static long fullCorrection(final LineParser1 parser, final int truth, final int end, final long row){
		assert((end==0 || end==1) && truth>=0 && truth<22) : "Validated endpoint/schema determines the one coordinate pair and target meaning";
		if(parser.terms()!=11 || !parser.termEquals("evaluation", 3)){
			throw new IllegalArgumentException("Require an eleven-field evaluation sidecar at full-error trial "+row);
		}
		for(int i=0; i<11; i++){
			if(parser.length(i)==0){throw new IllegalArgumentException("Empty full-error sidecar field at trial "+row);}
		}
		final int model=parser.parseInt(4), sample=parser.parseInt(6);
		final int rawStart=parser.parseInt(7), rawStop=parser.parseInt(8);
		final int trueStart=parser.parseInt(9), trueStop=parser.parseInt(10);
		final long correction=end==0 ? (long)trueStart-rawStart : (long)trueStop-rawStop;
		final int expected=correction<-RADIUS || correction>RADIUS ? CLASSES : (int)correction+RADIUS;
		if(model<0 || sample<0 || rawStart<0 || rawStop<rawStart || trueStart<0 || trueStop<trueStart || truth!=expected){
			throw new IllegalArgumentException("Full-error sidecar coordinates/target disagree at trial "+row);
		}
		return correction;
	}

	/** Three small counter sets over the same streamed trials; no per-row allocation. */
	public static final class FullComparison{
		/** Validates one reused score buffer before deriving either policy. */
		private void accept(final float[] scores, final long correction){
			assert(scores!=null && scores.length==22) : "Full-error inference must supply all21offsets and its non-positional head";
			int best=0;
			for(int i=0; i<22; i++){
				if(!Float.isFinite(scores[i])){throw new IllegalArgumentException("Full-error inference produced a nonfinite score");}
				if(i<CLASSES && scores[i]>scores[best]){best=i;}
			}
			float runner=-Float.MAX_VALUE;
			for(int i=0; i<CLASSES; i++){
				if(i!=best){runner=Math.max(runner, scores[i]);}
			}
			final boolean headWin=scores[CLASSES]>scores[best];
			argmax22=headWin ? CLASSES : best;
			predicted0=choose(scores, best, runner, 0, correction, headWin, margin0);
			predicted01=choose(scores, best, runner, .1f, correction, headWin, margin01);
			doNothing.accept(RADIUS, correction, false, false);
		}

		/** Strict offset/current/runner/OOR comparisons apply even at zero margin. */
		private static int choose(final float[] scores, final int best, final float runner, final float margin,
				final long correction, final boolean headWin, final FullResult result){
			final double top=scores[best];
			final boolean otherwise=best!=RADIUS && top>(double)scores[RADIUS]+margin && top>(double)runner+margin;
			final boolean passes=top>(double)scores[CLASSES]+margin;
			final int selected=NcrnaEndpointConfidence.selectFull(scores, RADIUS, margin);
			result.accept(selected, correction, headWin, otherwise && !passes && !headWin);
			return selected;
		}

		/** A do-nothing row and both policies use exactly the same complete denominator. */
		public ByteBuilder summary(final int end){
			return summary(end, false);
		}

		/** Optional paired changes leave all existing summary columns and defaults intact. */
		public ByteBuilder summary(final int end, final boolean extended){
			assert(margin0.trials>0 && margin0.trials==margin01.trials && margin0.trials==doNothing.trials)
				: "All policy rates must describe the same fully read trial stream";
			final ByteBuilder text=new ByteBuilder(FULL_HEADER);
			if(extended){text.append(CHANGE_HEADER);}
			text.nl();
			doNothing.row(text, "do_nothing", end, extended, doNothing.exact);
			margin0.row(text, "margin0", end, extended, doNothing.exact);
			margin01.row(text, "margin0.1", end, extended, doNothing.exact);
			return text;
		}

		public final FullResult margin0=new FullResult(), margin01=new FullResult(), doNothing=new FullResult();
		private int predicted0, predicted01, argmax22;
	}

	/** Cohorts come from actual coordinates, not from an out-of-range target index. */
	public static final class FullResult{
		/** Counts a real selected offset in its disjoint actual-coordinate cohort. */
		private void accept(final int selected, final long correction, final boolean headWin, final boolean vetoOnly){
			assert(selected>=0 && selected<CLASSES && !(headWin && vetoOnly))
				: "Only real offsets move an end; explicit OOR wins and additional OOR vetoes are disjoint";
			final long shift=selected-RADIUS, error=shift-correction, absolute=Math.abs(error);
			final boolean moved=shift!=0, reachable=correction>=-RADIUS && correction<=RADIUS;
			trials++;
			final long before=Math.abs(correction);
			if(absolute<before){better++;}
			else if(absolute>before){worse++;}
			else{tied++;}
			if(absolute==0){exact++;}
			if(absolute<=1){within1++;}
			if(absolute<=3){within3++;}
			if(correction==0){
				alreadyCorrect++;
				if(moved){movedCorrect++;}
			}else if(reachable){
				reachableWrong++;
				if(error==0){fixedReachableWrong++;}
			}else{
				outOfReach++;
				if(absolute<before){improvedOutOfReach++;}
				else if(absolute>before){worsenedOutOfReach++;}
				else{tiedOutOfReach++;}
				if(moved){movedOutOfReach++;}
				else{heldOutOfReach++;}
			}
			if(reachable){
				if(headWin){reachableOorWins++;}
				if(vetoOnly){reachableOorVetoOnly++;}
			}
			squaredBefore+=(double)correction*correction;
			squaredAfter+=(double)error*error;
			assert(alreadyCorrect+reachableWrong+outOfReach==trials && heldOutOfReach+movedOutOfReach==outOfReach
				&& better+worse+tied==trials && improvedOutOfReach+worsenedOutOfReach+tiedOutOfReach==outOfReach
				&& exact<=within1 && within1<=within3 && within3<=trials)
				: "Disjoint coordinate cohorts and nested endpoint tolerances must conserve every full-tail trial";
		}

		private void row(final ByteBuilder out, final String policy, final int end, final boolean extended, final long baselineExact){
			assert(trials>0 && (end==0 || end==1)) : "A complete declared endpoint stream precedes full-error percentages";
			out.append(policy).tab().append(end==0 ? "5prime" : "3prime").tab().append(trials)
				.tab().append(exact).tab().append(within1).tab().append(within3);
			Result.percentage(out.tab(), exact, trials);
			Result.percentage(out.tab(), within1, trials);
			Result.percentage(out.tab(), within3, trials);
			out.tab().append(alreadyCorrect).tab().append(movedCorrect);
			Result.percentage(out.tab(), movedCorrect, alreadyCorrect);
			out.tab().append(reachableWrong).tab().append(fixedReachableWrong);
			Result.percentage(out.tab(), fixedReachableWrong, reachableWrong);
			out.tab().append(outOfReach).tab().append(heldOutOfReach).tab().append(movedOutOfReach);
			Result.percentage(out.tab(), heldOutOfReach, outOfReach);
			Result.percentage(out.tab(), movedOutOfReach, outOfReach);
			final long reachable=alreadyCorrect+reachableWrong;
			out.tab().append(reachable).tab().append(reachableOorWins);
			Result.percentage(out.tab(), reachableOorWins, reachable);
			out.tab().append(reachableOorVetoOnly);
			Result.percentage(out.tab(), reachableOorVetoOnly, reachable);
			final double rmsBefore=Math.sqrt(squaredBefore/trials);
			final double rmsAfter=Math.sqrt(squaredAfter/trials);
			out.tab().append(rmsBefore, 6).tab().append(rmsAfter, 6);
			if(extended){
				out.tab().append(rmsAfter-rmsBefore, 6).tab().append(better).tab().append(worse).tab().append(tied).tab();
				if(better+worse==0){out.append("NA");}
				else{out.append(better/(double)(better+worse), 6);}
				out.tab().append(improvedOutOfReach).tab().append(worsenedOutOfReach).tab().append(tiedOutOfReach)
					.tab().append(exact-baselineExact).tab().append(100.0*(exact-baselineExact)/trials, 6);
			}
			out.nl();
		}

		public long trials, exact, within1, within3, alreadyCorrect, movedCorrect;
		public long reachableWrong, fixedReachableWrong, outOfReach, heldOutOfReach, movedOutOfReach;
		public long reachableOorWins, reachableOorVetoOnly;
		public long better, worse, tied, improvedOutOfReach, worsenedOutOfReach, tiedOutOfReach;
		public double squaredBefore, squaredAfter;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Trial Counts         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primitive counts over usable jittered trials, never unique loci or genomes. */
	public static final class Result{
		/** Adds one caller-compatible first-maximum class prediction. */
		private void accept(final int predicted, final int truth){
			assert(predicted>=0 && predicted<CLASSES && truth>=0 && truth<CLASSES)
				: "Endpoint errors require class indices in the same centered output window";
			final int error=predicted-truth, absolute=Math.abs(error);
			trials++;
			if(absolute==0){exact++;}
			if(absolute<=1){within1++;}
			if(absolute<=3){within3++;}
			if(truth==RADIUS){
				alreadyCorrect++;
				if(predicted!=RADIUS){movedAlreadyCorrect++;}
			}
			errorSum+=error;
			squaredError+=(long)error*error;
		}

		/** One header/row; movement percentage uses only initially correct trials. */
		public ByteBuilder summary(final int end){
			assert(trials>0 && (end==0 || end==1)) : "A complete nonempty evaluation and its declared endpoint precede summary output";
			final ByteBuilder buffer=new ByteBuilder();
			buffer.append(RESULT_HEADER).nl();
			appendRow(buffer, end);
			return buffer.nl();
		}

		/** Appends one policy's numeric values without imposing a surrounding schema. */
		private void appendRow(final ByteBuilder buffer, final int end){
			assert(trials>0) : "All policy rows summarize the same nonempty validated trial stream";
			buffer.append(end==0 ? "5prime" : "3prime").tab().append(trials).tab().append(exact)
				.tab().append(within1).tab().append(within3);
			percentage(buffer.tab(), exact, trials);
			percentage(buffer.tab(), within1, trials);
			percentage(buffer.tab(), within3, trials);
			buffer.tab().append(alreadyCorrect).tab().append(movedAlreadyCorrect);
			percentage(buffer.tab(), movedAlreadyCorrect, alreadyCorrect);
			buffer.tab().append(errorSum/(double)trials, 6).tab().append(Math.sqrt(squaredError/(double)trials), 6);
		}

		/** An absent denominator is unavailable rather than a zero-percent success. */
		private static void percentage(final ByteBuilder buffer, final long numerator, final long denominator){
			if(denominator==0){buffer.append("NA");}
			else{buffer.append(100.0*numerator/denominator, 6);}
		}

		public long trials, exact, within1, within3, alreadyCorrect, movedAlreadyCorrect;
		public long errorSum, squaredError;
	}

	/** Three policies applied to the same scores; none performs another inference. */
	public static final class Comparison{
		private void accept(final float[] scores, final int truth){
			predicted=RrnaResourceIO.argmax(scores);
			gatedPrediction=NcrnaEndpointConfidence.select(scores, RADIUS, predicted, .1f);
			ungated.accept(predicted, truth);
			gated.accept(gatedPrediction, truth);
			doNothing.accept(RADIUS, truth);
			if(truth!=RADIUS && predicted==truth){
				if(gatedPrediction==truth){fixesKept++;}
				else{fixesLost++;}
			}
		}

		/** Each policy retains the legacy numeric columns and their denominators. */
		public ByteBuilder summary(final int end){
			final ByteBuilder text=new ByteBuilder("policy\t").append(RESULT_HEADER)
				.append("\tungated_fixes_kept\tungated_fixes_lost\n");
			row(text, "margin0", ungated, end, fixesKept+fixesLost, 0);
			row(text, "margin0.1", gated, end, fixesKept, fixesLost);
			row(text, "do_nothing", doNothing, end, 0, fixesKept+fixesLost);
			return text;
		}

		private static void row(final ByteBuilder text, final String policy, final Result result,
				final int end, final long kept, final long lost){
			text.append(policy).tab();
			result.appendRow(text, end);
			text.tab().append(kept).tab().append(lost).nl();
		}

		public final Result ungated=new Result(), gated=new Result(), doNothing=new Result();
		public long fixesKept, fixesLost;
		private int predicted, gatedPrediction;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in, rows, netPath, out, rawPath, comparisonPath;
	private String featureVariant="baseline";
	private int end=-1;
	private int fullInputs=69;
	private boolean overwrite, fullErrors, replayFull, extended;
	private PrintStream outstream=System.err;
	private static final int CLASSES=21, RADIUS=10;
	private static final String RESULT_HEADER="end\tusable_jitter_trials\texact\twithin1\twithin3\texact_pct\twithin1_pct\twithin3_pct"
		+"\talready_correct\tmoved_already_correct\tmoved_already_correct_pct\tmean_error_bp\trms_error_bp";
	private static final String FULL_HEADER="policy\tend\tusable_jitter_trials\texact\twithin1\twithin3\texact_pct\twithin1_pct\twithin3_pct"
		+"\talready_correct\tmoved_already_correct\tmoved_already_correct_pct\treachable_wrong\tfixed_reachable_wrong\tfixed_reachable_wrong_pct"
		+"\tout_of_reach\theld_out_of_reach\tmoved_out_of_reach\theld_out_of_reach_pct\tmoved_out_of_reach_pct"
		+"\treachable\treachable_oor_head_wins\treachable_oor_head_wins_pct\treachable_oor_veto_only\treachable_oor_veto_only_pct"
		+"\trms_before_bp\trms_after_bp";
	private static final String CHANGE_HEADER="\trms_change_bp\tbetter\tworse\ttied\tbetter_fraction"
		+"\timproved_out_of_reach\tworsened_out_of_reach\ttied_out_of_reach\texact_change\texact_change_pp";
	static final byte[] ROW_HEADER=("locus_id\tsource\tassembly\tsplit\tmodel_index\tmodel\tsample"
		+"\traw_start0\traw_stop0\ttrue_start0\ttrue_stop0").getBytes(StandardCharsets.US_ASCII);
}
