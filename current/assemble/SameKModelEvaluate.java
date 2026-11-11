package assemble;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.LongList;

/**
 * Scores a frozen same-K corpus fold and emits exact candidate/control curves.
 * Unknown labels were excluded upstream; neither population is resampled here.
 * Validates the full provenance sidecar while consuming feature rows only for the
 * selected fold. Inference is serial and reuses one feature vector and model state.
 * Scores are retained in memory, separately by kind; curve generation sorts them
 * in place. Raw prediction rows retain input order, while tied scores enter a curve
 * together. Outputs may be partial after failure; complete.txt is written last.
 * @author Fischl
 */
public final class SameKModelEvaluate{

	/** Runs inference without training or modifying the model file; closes the diagnostic stream. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SameKModelEvaluate.class, false);
		try{new SameKModelEvaluate(pp.args, pp.outstream).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires an explicit fold, model and matching full corpus sidecar. */
	private SameKModelEvaluate(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		for(String arg : args){
			final String[] parts=arg.split("=", 2);
			String name=parts[0].toLowerCase(Locale.ROOT);
			while(name.startsWith("-")){name=name.substring(1);}
			final String value=parts.length>1 ? parts[1] : null;
			if(name.equals("rows")){
				rows=value;
			}else if(name.equals("net")){
				network=value;
			}else if(name.equals("fold")){
				fold=Arrays.asList(FOLDS).indexOf(value);
			}else if(name.equals("simdff")){
				simdFF=Parse.parseBoolean(value);
			}else if(!parser.parse(arg, name, value)){throw new IllegalArgumentException("Unknown evaluation option: "+arg);}
		}
		input=parser.in1; output=parser.out1;
		if(input==null || rows==null || network==null || output==null || fold<0 || parser.in2!=null || parser.out2!=null ||
				!Tools.testInputFiles(false, true, input, rows, network) || new File(output).exists()){
			throw new IllegalArgumentException("Required in=fold.tsv rows=rows.tsv net=model.bbnet fold=validation out=NEW_directory.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Joins selected-fold rows in sidecar order, requiring exact headers, matching targets and
	 * finite nonnegative features. Validates every sidecar row's fold, kind and ordinal even
	 * outside the selected fold. Saves normalized-zero score bits and decimal scores before
	 * generating curves; a successful prediction file alone is not a complete evaluation.
	 * Failed feature-reader acquisition closes metadata and preserves the original unchecked
	 * failure, with any metadata-cleanup failure suppressed.
	 * Failed prediction-output acquisition attempts both input closes and rethrows the
	 * original failure with cleanup failures suppressed. Partial output ownership is separate.
	 */
	private void run(){
		final CellNet net=CellNetParser.load(network, false);
		if(net==null || net.numInputs()!=27 || net.numOutputs()!=1){throw new IllegalArgumentException("Expected a27-input binary model.");}
		net.simdFF=simdFF;
		if(!new File(output).mkdir()){throw new IllegalArgumentException("Cannot create new evaluation directory.");}
		final LongList[] scores={new LongList(), new LongList()};
		final long[] ordinals=new long[3];
		final float[] vector=new float[27];
		final LineParser1 mp=new LineParser1('\t'), dp=new LineParser1('\t');
		final ByteBuilder buffer=new ByteBuilder();
		final ByteFile meta=ByteFile.makeByteFile(rows, true);
		final ByteFile data;
		try{data=ByteFile.makeByteFile(input, true);}
		catch(RuntimeException|Error failure){
			try{
				if(meta.close()){failure.addSuppressed(new IllegalStateException("Evaluation metadata cleanup failed."));}
			}catch(RuntimeException|Error cleanup){
				if(cleanup!=failure){failure.addSuppressed(cleanup);}
			}
			throw failure;
		}
		final ByteStreamWriter predictions;
		try{predictions=open(output+"/scores.tsv");}
		catch(RuntimeException|Error failure){
			// ASM024: output acquisition can fail while both inputs are already owned.
			try{
				try{
					if(meta.close()){failure.addSuppressed(new IllegalStateException("Evaluation metadata cleanup failed."));}
				}catch(RuntimeException|Error cleanup){
					if(cleanup!=failure){failure.addSuppressed(cleanup);}
				}
			}finally{
				try{
					try{
						if(data.close()){failure.addSuppressed(new IllegalStateException("Evaluation feature input cleanup failed."));}
					}catch(RuntimeException|Error cleanup){
						if(cleanup!=failure){failure.addSuppressed(cleanup);}
					}
				}finally{throw failure;}// Preserve the primary even if recording a cleanup failure fails.
			}
		}
		boolean error=false;
		try{
			require(meta.nextLine(), "fold\trow\tgenome\tkind\tevent\ttarget");
			require(data.nextLine(), "#dims\t27\t1");
			require(data.nextLine(), "#schema=samek27-v1; unknown labels excluded; no class balancing");
			predictions.print("fold\trow\tgenome\tkind\tevent\ttarget\tscore_bits\tscore\n");
			for(byte[] line=meta.nextLine(); line!=null; line=meta.nextLine()){
				mp.set(line);
				if(mp.terms()!=6){throw new IllegalArgumentException("Wrong corpus provenance width.");}
				int rowFold=-1;
				for(int i=0; i<3; i++){if(mp.termEquals(FOLDS[i], 0)){rowFold=i;}}
				final int kind=mp.termEquals("candidate", 3) ? 0 : mp.termEquals("control", 3) ? 1 : -1;
				final int label=mp.parseInt(5);
				if(rowFold<0 || kind<0 || mp.parseLong(1)!=++ordinals[rowFold] || mp.parseLong(4)<1 || (label!=0 && label!=1)){
					throw new IllegalArgumentException("Invalid fold/kind/ordinal/event/target in sidecar.");
				}
				if(rowFold!=fold){continue;}
				final byte[] values=data.nextLine();
				if(values==null){throw new IllegalArgumentException("Missing paired feature row.");}
				dp.set(values);
				if(dp.terms()!=28 || dp.parseInt(27)!=label){throw new IllegalArgumentException("Feature width/target disagrees with provenance.");}
				for(int i=0; i<27; i++){
					vector[i]=dp.parseFloat(i);
					if(!Float.isFinite(vector[i]) || vector[i]<0){throw new IllegalArgumentException("Invalid normalized feature.");}
				}
				final float score=net.applyInput(vector).feedForward();
				final long packed=pack(score, label);
				scores[kind].add(packed);
				predictions.print(buffer.clear().append(line).tab().append(packed>>>1).tab().appendSlow(score).nl());
			}
			if(data.nextLine()!=null || ordinals[fold]==0){throw new IllegalArgumentException("Unpaired extra rows or empty selected fold.");}
			//TODO: Probable bug ASM024 - an unchecked close failure can skip later cleanup; separate from second-reader acquisition.
		}finally{error=meta.close() | data.close() | predictions.poisonAndWait();}
		if(error){throw new IllegalStateException("Incomplete input/prediction output.");}
		final ByteStreamWriter summary=open(output+"/summary.tsv");
		try{
			summary.print("kind\trows\tpositive\tnegative\tdistinct_scores\tsimd_ff\n");
			for(int kind=0; kind<2; kind++){
				final String name=kind==0 ? "candidate" : "control";
				final ByteStreamWriter curve=open(output+"/"+name+".tsv");
				final long[] totals;
				try{totals=writeCurve(scores[kind], curve);}
				finally{if(curve.poisonAndWait()){throw new IllegalStateException("Curve write failed.");}}
				final ByteStreamWriter points=open(output+"/"+name+"_points.tsv");
				try{writePoints(scores[kind], totals[0], totals[1], points);}
				finally{if(points.poisonAndWait()){throw new IllegalStateException("Operating-point write failed.");}}
				summary.print(buffer.clear().append(name).tab().append(scores[kind].size).tab().append(totals[0])
						.tab().append(totals[1]).tab().append(totals[2]).tab().append(simdFF).nl());
			}
		}finally{if(summary.poisonAndWait()){throw new IllegalStateException("Summary write failed.");}}
		final ByteStreamWriter complete=open(output+"/complete.txt");
		final String message="SAME_K_EVALUATE_PASS fold="+FOLDS[fold]+" candidates="+scores[0].size+" controls="+scores[1].size+" simdff="+simdFF;
		complete.print(message+"\n");
		if(complete.poisonAndWait()){throw new IllegalStateException("Completion marker write failed.");}
		outstream.println(message);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Encodes a finite probability in [0,1] and a binary label for numeric sorting.
	 * Normalizes signed zero, shifts the float bits left once and stores the label
	 * in the low bit; the result is the input representation required by writeCurve. */
	public static long pack(final float score, final int label){
		if(!Float.isFinite(score) || score<0 || score>1 || (label!=0 && label!=1)){
			throw new IllegalArgumentException("Expected finite score in[0,1] and binary label: "+score+", "+label);
		}
		return ((long)Float.floatToRawIntBits(score==0 ? 0 : score)<<1)|label;
	}

	/**
	 * Sorts the supplied packed population in place and emits every attainable decision set
	 * for score>=cutoff, including reject-all and whole tied groups. Returns positive count,
	 * negative count and distinct-score count, in that order; an empty population has only
	 * the reject-all row. Inputs must have the encoding and range enforced by pack.
	 * The caller owns scores and the already-started writer. This method changes the
	 * score order and appends output; it neither starts nor closes the writer.
	 */
	public static long[] writeCurve(final LongList scores, final ByteStreamWriter writer){
		assert(scores!=null && writer!=null) : "Exact threshold curves require a complete population and output writer.";
		scores.sort();
		long positive=0;
		for(int i=0; i<scores.size; i++){positive+=scores.array[i]&1L;}
		final long negative=scores.size-positive;
		long tp=0, fp=0, distinct=0;
		final ByteBuilder buffer=new ByteBuilder();
		writer.print("threshold_bits\tthreshold\tgroup_positive\tgroup_negative\ttp\tfp\ttn\tfn\tprecision\trecall\tfpr\tfnr\n");
		final float maximum=scores.size==0 ? 1 : Float.intBitsToFloat((int)(scores.array[scores.size-1]>>>1));
		append(buffer, Math.nextUp(maximum), 0, 0, tp, fp, positive, negative);
		for(int stop=scores.size; stop>0; ){
			final int bits=(int)(scores.array[stop-1]>>>1);
			int start=stop-1;
			long gp=0, gn=0;
			while(start>=0 && (int)(scores.array[start]>>>1)==bits){
				if((scores.array[start]&1L)==1){gp++;}else{gn++;}
				start--;
			}
			tp+=gp; fp+=gn; distinct++;
			append(buffer, Float.intBitsToFloat(bits), gp, gn, tp, fp, positive, negative);
			if(buffer.length()>65536){writer.print(buffer); buffer.clear();}
			stop=start+1;
		}
		assert(tp==positive && fp==negative) : "The final attainable cutoff must accept every score in this population.";
		writer.print(buffer);
		return new long[]{positive, negative, distinct};
	}

	/**
	 * Emits the fixed cutoff grid from the already-sorted population and its complete totals.
	 * Group counts include only scores equal to each cutoff; cumulative counts include all
	 * scores at or above it. Full exact curves remain the authority for threshold selection.
	 */
	private static void writePoints(final LongList scores, final long positive, final long negative, final ByteStreamWriter writer){
		assert(positive+negative==scores.size) : "Operating-point denominators must match the population used by the full curve.";
		writer.print("threshold_bits\tthreshold\tgroup_positive\tgroup_negative\ttp\tfp\ttn\tfn\tprecision\trecall\tfpr\tfnr\n");
		final ByteBuilder buffer=new ByteBuilder();
		for(float cutoff : CUTOFFS){
			long tp=0, fp=0, gp=0, gn=0;
			for(int i=scores.size-1; i>=0; i--){
				final long packed=scores.array[i];
				final float score=Float.intBitsToFloat((int)(packed>>>1));
				if(score<cutoff){break;}
				if((packed&1L)==1){tp++; if(score==cutoff){gp++;}}else{fp++; if(score==cutoff){gn++;}}
			}
			append(buffer, cutoff, gp, gn, tp, fp, positive, negative);
		}
		writer.print(buffer);
	}

	/** Computes class-conditional rates; undefined ratios are explicitly NA, not zero. */
	private static void append(final ByteBuilder out, final float cutoff, final long gp, final long gn,
			final long tp, final long fp, final long positive, final long negative){
		assert(tp<=positive && fp<=negative) : "Cumulative groups cannot exceed their complete population counts.";
		out.append(Float.floatToRawIntBits(cutoff)).tab().appendSlow(cutoff).tab().append(gp).tab().append(gn)
				.tab().append(tp).tab().append(fp).tab().append(negative-fp).tab().append(positive-tp).tab();
		ratio(out, tp, tp+fp); ratio(out, tp, positive); ratio(out, fp, negative); ratio(out, positive-tp, positive);
		out.setLength(out.length()-1); out.nl();
	}

	/** Appends one exact-denominator rate and a delimiter, preserving undefined values. */
	private static void ratio(final ByteBuilder out, final long numerator, final long denominator){
		if(denominator==0){out.append("NA");}else{out.append(numerator/(double)denominator, 12);}
		out.tab();
	}

	/** Opens only new output files under the newly created evaluation directory. */
	private static ByteStreamWriter open(final String path){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start();
		return writer;
	}

	/** Requires the exact upstream schema before joining rows from separate files. */
	private static void require(final byte[] actual, final String expected){
		if(!Arrays.equals(actual, expected.getBytes(StandardCharsets.UTF_8))){throw new IllegalArgumentException("Expected header: "+expected);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Diagnostic stream owned and closed by main. */
	private final PrintStream outstream;
	/** Selected feature table, full provenance sidecar, saved model and new output directory. */
	private String input, rows, network, output;
	/** Required index into FOLDS; -1 means no supported fold was selected. */
	private int fold=-1;
	/** Per-model feed-forward switch; this does not set the independent global SIMD/FMA flags. */
	private boolean simdFF=false;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	/** Exact accepted sidecar fold names, also used to index independent ordinal counters. */
	private static final String[] FOLDS={"train", "validation", "test"};
	/** Descending reporting grid, beginning above the permitted score range to reject all rows. */
	private static final float[] CUTOFFS={Math.nextUp(1.0f), 1, 0.9999f, 0.9995f, 0.999f, 0.995f, 0.99f,
			0.975f, 0.95f, 0.925f, 0.9f, 0.85f, 0.8f, 0.75f, 0.7f, 0.65f, 0.6f, 0.5f, 0.25f, 0};
}
