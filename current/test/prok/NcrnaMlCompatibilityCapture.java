package test.prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import map.ObjectSet;
import ml.Cell;
import ml.CellNet;
import ml.CellNetParser;
import ml.Function;
import ml.Sample;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO;
import prok.RrnaEndpointVector;
import prok.RrnaPositionalKmerTable;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Captures loaded network state and exact prediction bits for a reviewed roster.
 * This diagnostic never selects, trains or re-exports a supplied network. Each
 * roster row binds one standalone network or explicit zero-based bundle frame
 * to saved raw feature rows or explicitly SYNTHETIC seeded feature bounds.
 * Synthetic edges supplement those rows. Normal, ordinary worker-copy, inference-only and dense optimized
 * modes are recorded independently, without requiring their arithmetic to agree.
 * A separate seeded legacy tanh fixture records initialization and one update
 * entirely in memory; it is not an additional biological fit.
 * Adapted from template.A_SampleByteFile for a small manifest-driven diagnostic.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaMlCompatibilityCapture{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs the capture and closes redirected diagnostics on every exit. */
	public static void main(final String[] args){
		final PreParser pre=new PreParser(args, null, false);
		try{
			new NcrnaMlCompatibilityCapture(pre.args, pre.outstream).process();
		}finally{
			Shared.closeStream(pre.outstream);
		}
	}

	/** Accepts a frozen manifest, fresh destination and bounded real-row count. */
	private NcrnaMlCompatibilityCapture(final String[] args, final PrintStream diagnostics){
		outstream=diagnostics;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase();
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("rows")){
				rowLimit=Integer.parseInt(b);
			}else if(a.equals("legacy")){
				legacy=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				// Standard in/out, thread and SIMD controls.
			}else{
				throw new IllegalArgumentException("Unknown parameter "+arg);
			}
		}
		manifest=parser.in1;
		out=parser.out1;
		if(manifest==null || out==null || rowLimit<1 || rowLimit>1024 || parser.overwrite
				|| parser.append || parser.in2!=null || parser.out2!=null){
			throw new IllegalArgumentException("Require in=manifest out=fresh_capture rows=1..1024 ow=f");
		}
		if(!Tools.testInputFiles(false, true, manifest)
				|| !Tools.testOutputFiles(false, false, false, out)
				|| !Tools.testForDuplicateFiles(true, manifest, out)){
			throw new IllegalArgumentException("Manifest and fresh capture paths must be distinct and accessible");
		}
		assert(rowLimit>0) : "Real rows are mandatory; synthetic edges cannot replace observed feature inputs";
	}

	/*--------------------------------------------------------------*/
	/*----------------           Capture            ----------------*/
	/*--------------------------------------------------------------*/

	/** Captures every manifest row; a partial file is never a completion receipt. */
	private void process(){
		final LineParser1 parser=new LineParser1('\t');
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		int count=0;
		float[][] trainingRows=null;
		try(NcrnaTextIO.Input input=new NcrnaTextIO.Input(manifest);
				NcrnaTextIO.Output output=new NcrnaTextIO.Output(out)){
			final byte[] header=input.reader.nextLine();
			if(header==null || !new String(header, java.nio.charset.StandardCharsets.US_ASCII).equals(HEADER)){
				throw new IllegalArgumentException("Compatibility manifest header must be "+HEADER);
			}
			write(output, buffer.append("SIMD\t").append(Shared.SIMD).nl());
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				if(parser.terms()!=8){throw new IllegalArgumentException("Require eight manifest fields at network "+count);}
				final String id=parser.parseString(0);
				final String path=parser.parseString(1);
				final String rows=parser.parseString(3);
				final String provenance=parser.parseString(6);
				final int frame=parser.parseInt(2);
				final int width=parser.parseInt(4);
				final int outputs=parser.parseInt(5);
				final long seed=parser.parseLong(7);
				if(id.isEmpty() || !ids.add(id) || frame< -1 || width<1 || outputs<1){
					throw new IllegalArgumentException("Invalid or repeated network descriptor: "+id);
				}
				final boolean synthetic="SYNTHETIC".equals(provenance);
				if((synthetic && seed<0) || (!synthetic && (!"REAL".equals(provenance) || seed!= -1))){
					throw new IllegalArgumentException("Require REAL/-1 or SYNTHETIC/nonnegative seed for "+id);
				}
				final ArrayList<byte[]> image=networkImage(path, frame);
				final CellNet ordinary=frame<0 ? CellNetParser.load(path, false) : CellNetParser.loadFromLines(image);
				if(ordinary.numInputs()!=width || ordinary.numOutputs()!=outputs){
					throw new IllegalArgumentException("Manifest shape disagrees with loaded network "+id);
				}
				final float[][] inputs=synthetic ? syntheticRows(rows, width, seed) : realRows(rows, width);
				if(!synthetic && width==10 && trainingRows==null){trainingRows=inputs;}
				write(output, buffer.append("NETWORK\t").append(id).tab().append(path).tab().append(frame).nl());
				write(output, buffer.append("FEATURE_SOURCE\t").append(provenance).tab().append(rows).tab().append(seed).nl());
				capture(output, "ordinary", ordinary, inputs);
				capture(output, "worker_copy", ordinary.copy(false), inputs);
				capture(output, "inference", CellNetParser.loadInferenceFromLines(image), inputs);
				capture(output, "dense_optimized", ordinary.copyDenseForInference(), inputs);
				count++;
			}
			if(count==0){throw new IllegalArgumentException("A compatibility capture cannot have an empty network roster");}
			if(legacy){
				if(trainingRows==null){throw new IllegalArgumentException("The legacy tanh fixture requires saved ten-input rows; use legacy=f for a network-only diagnostic");}
				legacyUpdate(output, trainingRows);
			}
			write(output, buffer.append("COMPLETE\t").append(count).nl());
		}
		outstream.println("Captured "+count+" networks; state and predictions are raw FP32 bits.");
	}

	/** Records native seeded initialization and one two-sample tanh update.
	 * Reflection exposes only the existing package-private transpose operation,
	 * required by CellNet's documented training contract; no ML API is modified.
	 */
	private void legacyUpdate(final NcrnaTextIO.Output output, final float[][] rows){
		assert(rows.length>=1 && rows[0].length==10) : "The legacy fixture uses actual ten-input feature rows";
		Arrays.fill(Function.TYPE_RATES, 0f);
		Function.TYPE_RATES_CUM=null;
		Cell.defaultActivationType=Function.TANH;
		Cell.finalLayerType=Function.RSLOG;
		Cell.randomTypeRate=0;
		for(final boolean dense:new boolean[]{true, false}){
			final CellNet net=new CellNet(new int[]{10, 16, 1}, 58, 1f, 1f, 1, new ArrayList<String>(), dense);
			net.randomize();
			write(output, buffer.append("LEGACY_TANH\t").append(dense).nl());
			capture(output, "initialized", net, rows);
			try{
				final java.lang.reflect.Method transpose=CellNet.class.getDeclaredMethod("transpose");
				transpose.setAccessible(true);
				transpose.invoke(net);
			}catch(final ReflectiveOperationException failure){
				throw new IllegalStateException("Cannot prepare native legacy transpose for backpropagation", failure);
			}
			for(int i=0; i<2; i++){
				final Sample sample=new Sample(rows[i%rows.length], new float[]{i==0 ? 0.25f : 0.75f}, 1f, i);
				net.processSample(sample, true, 1f);
				bits(output, "TRAIN_PREDICTION", sample.result);
			}
			net.applyChanges(2, 0.001f);
			capture(output, "updated", net, rows);
		}
	}

	/** Reads one normal image or checks the complete framing of a requested bundle. */
	private static ArrayList<byte[]> networkImage(final String path, final int frame){
		assert(frame>= -1) : "The manifest uses -1 for a standalone image and nonnegative bundle indices";
		final ArrayList<byte[]> image=new ArrayList<byte[]>();
		try(NcrnaTextIO.Input input=new NcrnaTextIO.Input(path)){
			if(frame<0){
				for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){image.add(line);}
			}else{
				final LineParser1 parser=new LineParser1('\t');
				int declared=-1;
				int seen=0;
				boolean ended=false;
				for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
					parser.set(line);
					if(parser.termEquals("networks", 0)){
						if(declared>=0 || parser.terms()!=2){throw new IllegalArgumentException("Repeated/invalid bundle count "+path);}
						declared=parser.parseInt(1);
					}else if(parser.termEquals("network", 0)){
						if(ended || declared<1 || parser.terms()!=6){throw new IllegalArgumentException("Invalid bundle frame "+path);}
						final int lines=parser.parseInt(5);
						if(lines<1 || lines>1000000){throw new IllegalArgumentException("Invalid frame size "+lines);}
						for(int i=0; i<lines; i++){
							final byte[] payload=input.reader.nextLine();
							if(payload==null){throw new IllegalArgumentException("Truncated bundle "+path);}
							if(seen==frame){image.add(payload);}
						}
						seen++;
					}else if(parser.termEquals("#END_RRNA_ENDPOINT_BBNETS", 0)){
						if(ended || seen!=declared){throw new IllegalArgumentException("Bundle end/count mismatch "+path);}
						ended=true;
					}else if(ended){throw new IllegalArgumentException("Trailing bundle content "+path);}
				}
				if(!ended || seen!=declared || frame>=seen){throw new IllegalArgumentException("Missing bundle frame/end "+path);}
			}
		}
		if(image.isEmpty()){throw new IllegalArgumentException("Missing network image "+path+" frame "+frame);}
		return image;
	}

	/** Retains exactly the first requested raw rows, with no padding or width conversion. */
	private float[][] realRows(final String path, final int width){
		assert(width>0) : "CellNet.applyInput requires an exact positive input width";
		final float[][] rows=new float[rowLimit][width];
		final LineParser1 parser=new LineParser1('\t');
		int labels=-1;
		int count=0;
		try(NcrnaTextIO.Input input=new NcrnaTextIO.Input(path)){
			for(byte[] line=input.reader.nextLine(); line!=null && count<rowLimit; line=input.reader.nextLine()){
				if(line.length==0){continue;}
				parser.set(line);
				if(parser.termEquals("#dims", 0)){
					if(labels>=0 || (parser.terms()!=3 && parser.terms()!=4) || parser.parseInt(1)!=width
							|| (parser.terms()==4 && parser.parseInt(3)!=0)){
						throw new IllegalArgumentException("Wrong raw matrix shape/header "+path);
					}
					labels=parser.parseInt(2);
					if(labels<1){throw new IllegalArgumentException("Invalid label width "+path);}
				}else if(line[0]!='#'){
					if(labels<0 || parser.terms()!=width+labels){throw new IllegalArgumentException("Wrong raw row width "+path);}
					for(int i=0; i<width; i++){
						final float value=parser.parseFloat(i);
						if(!Float.isFinite(value)){throw new IllegalArgumentException("Nonfinite raw feature "+path);}
						rows[count][i]=value;
					}
					count++;
				}
			}
		}
		if(count!=rowLimit){throw new IllegalArgumentException("Not enough real rows in "+path+": "+count);}
		return rows;
	}

	/** Generates labeled synthetic compatibility rows, never biological training data.
	 * Recipes preserve coupled tRNA count profiles or reuse the planted endpoint
	 * table from RrnaEndpointVectorTest, generalized to radius16. These are fixture
	 * sampling intervals, not biological bounds. Every generated input bit enters
	 * the capture so a changed recipe or seed cannot pass replay silently.
	 */
	private float[][] syntheticRows(final String path, final int width, final long seed){
		assert(width>0 && seed>=0) : "Synthetic rows require explicit dimensionality and a repeatable manifest seed";
		final shared.Random randy=Shared.threadLocalRandom(seed);
		if(width==14 && path.equals("synthetic:trna_profiles_v1")){return trnaRows(randy);}
		if((width==28 && path.equals("synthetic:rrna_planted_r12_v1"))
				|| (width==36 && path.equals("synthetic:rrna_planted_r16_v1"))){
			return endpointRows(randy, width);
		}
		throw new IllegalArgumentException("Synthetic recipe/width mismatch: "+path+" / "+width);
	}

	/** Samples coupled count profiles, binary endpoint type and both fuzziness sentinels.
	 * TrnaBoundaryScorer.buildInput supplies this14-column order; the profile
	 * formula is TrnaBoundaryFeatures.enrichmentProfileCore. Count0..1,000,000
	 * and40..150base length are fixture intervals, not asserted biological limits.
	 */
	private float[][] trnaRows(final shared.Random randy){
		assert(randy!=null) : "A manifest seed owns every synthetic count and feature draw";
		final float[][] rows=new float[rowLimit+4][14];
		for(int row=0; row<rows.length; row++){
			final float[] values=rows[row];
			values[0]=randy.nextFloat();
			values[1]=randy.nextFloat();
			profile(values, 2, row==0 ? 0 : row==1 ? 1000000 : randy.nextInt(1000001),
				row<2 ? 0 : randy.nextInt(1000001));
			values[5]=row%2;
			for(int i=6; i<9; i++){values[i]=row%6==0 ? -1f : row%6==1 ? -2f : randy.nextFloat();}
			values[9]=(40+randy.nextInt(111))/76f;
			values[10]=randy.nextFloat();
			profile(values, 11, row<2 ? 0 : randy.nextInt(1000001),
				row==1 ? 1000000 : row==0 ? 0 : randy.nextInt(1000001));
		}
		return rows;
	}

	/** Preserves the shared denominator and a unit ratio in every profile. */
	private static void profile(final float[] row, final int offset, final int center, final int other){
		assert(center>=0 && other>=0) : "Synthetic profiles are generated from nonnegative observation counts";
		final float denominator=1f+Math.max(center, other);
		row[offset]=(float)(Math.log(denominator)/Math.log(2)/8);
		row[offset+1]=(center+1)/denominator;
		row[offset+2]=(other+1)/denominator;
		assert(row[offset+1]==1f || row[offset+2]==1f) : "TrnaBoundaryFeatures profile ratios share their maximum-count denominator";
	}

	/** Reuses the existing planted six-mer feature fixture at the requested radius.
	 * Table score extrema include every key, including absent keys. Extras include
	 * GC/identity0and1 and length-ratio.25and4 as explicit fixture extremes only.
	 */
	private float[][] endpointRows(final shared.Random randy, final int width){
		assert(width==28 || width==36) : "Endpoint fixture width is2R+4, for the shipped R12/R16 networks";
		final int radius=(width-4)/2;
		final int sites=width-3;
		final byte[] bases=new byte[180];
		Arrays.fill(bases, (byte)'A');
		System.arraycopy("CGTACG".getBytes(java.nio.charset.StandardCharsets.US_ASCII), 0, bases, 47, 6);
		System.arraycopy("TGCATG".getBytes(java.nio.charset.StandardCharsets.US_ASCII), 0, bases, 125, 6);
		final RrnaPositionalKmerTable.Table five=new RrnaPositionalKmerTable.Table("test", "5prime", 6, -3, 0.5, radius);
		final RrnaPositionalKmerTable.Table three=new RrnaPositionalKmerTable.Table("test", "3prime", 6, -4, 0.5, radius);
		for(int i=0; i<100; i++){
			five.add(bases, 50);
			three.add(bases, 129);
		}
		five.finish();
		three.finish();
		final float[][] rows=new float[rowLimit+4][width];
		final float[] scratch=new float[sites];
		for(int row=0; row<rows.length; row++){
			final RrnaPositionalKmerTable.Table table=row%2==0 ? five : three;
			final int start=50+randy.nextInt(2*radius+1)-radius;
			final int stop=129+randy.nextInt(2*radius+1)-radius;
			if(!RrnaEndpointVector.fill(table, bases, start, stop, 80, randy.nextFloat(), randy.nextFloat(), scratch, rows[row])){
				throw new IllegalStateException("Complete planted sequence must supply every synthetic endpoint feature");
			}
			if(row>=rowLimit){
				float low=Float.POSITIVE_INFINITY;
				float high=Float.NEGATIVE_INFINITY;
				for(final double score:table.scores){
					low=Math.min(low, (float)score);
					high=Math.max(high, (float)score);
				}
				for(int i=0; i<sites; i++){rows[row][i]=row==rowLimit ? low : row==rowLimit+1 ? high : i%2==0 ? low : high;}
				rows[row][sites]=row%2==0 ? 0.25f : 4f;
				rows[row][sites+1]=row%2;
				rows[row][sites+2]=row%2==0 ? 1f : 0f;
			}
		}
		return rows;
	}

	/** Records state, observed inputs and small deterministic edge vectors in one mode. */
	private void capture(final NcrnaTextIO.Output out, final String mode, final CellNet net, final float[][] rows){
		assert(net!=null && rows.length>0) : "Every manifest network and mode requires real feature rows";
		write(out, buffer.append("MODE\t").append(mode).tab().append(net.dense).nl());
		buffer.append("DIMS");
		for(final int value:net.dims){buffer.tab().append(value);}
		write(out, buffer.nl());
		bits(out, "MEAN", net.inputMeanCopy());
		bits(out, "SCALE", net.inputInverseStdCopy());
		write(out, buffer.append("CLIP\t").append(Float.floatToRawIntBits(net.inputClip())).nl());
		for(int layer=0; layer<net.net.length; layer++){
			for(int index=0; index<net.net[layer].length; index++){
				final Cell cell=net.net[layer][index];
				write(out, buffer.append("CELL\t").append(layer).tab().append(index).tab()
					.append(cell.function==null ? "null" : cell.function.name()).tab()
					.append(Float.floatToRawIntBits(cell.bias)).nl());
				buffer.append("INPUTS\t").append(cell.inputs==null ? -1 : cell.inputs.length);
				if(cell.inputs!=null){for(final int value:cell.inputs){buffer.tab().append(value);}}
				write(out, buffer.nl());
				bits(out, "WEIGHTS", cell.weights);
			}
		}
		for(int row=0; row<rows.length; row++){predict(out, net, row, rows[row]);}
		final float[] edge=new float[net.numInputs()];
		for(int row=0; row<8; row++){
			for(int i=0; i<edge.length; i++){
				edge[i]=row==0 ? 0f : row==1 ? -0f : row==2 ? 1f : row==3 ? -1f :
					row==4 ? (i%2==0 ? 4f : -4f) : row==5 ? Float.MIN_NORMAL :
					row==6 ? (i==edge.length/2 ? 1f : 0f) : (i%2==0 ? 0.000001f : -0.000001f);
			}
			predict(out, net, rows.length+row, edge);
		}
	}

	/** Saves raw inputs and every raw output bit, including signed zero. */
	private void predict(final NcrnaTextIO.Output out, final CellNet net, final int row, final float[] input){
		assert(input.length==net.numInputs()) : "Width binding must precede CellNet inference";
		write(out, buffer.append("ROW\t").append(row).nl());
		bits(out, "RAW", input);
		net.applyInput(input);
		net.feedForward();
		buffer.append("PREDICTIONS");
		for(int i=0; i<net.numOutputs(); i++){
			final float value=net.getOutput(i);
			if(!Float.isFinite(value)){throw new IllegalStateException("Nonfinite capture output at row "+row+" output "+i);}
			buffer.tab().append(Float.floatToRawIntBits(value));
		}
		write(out, buffer.nl());
	}

	/** Writes a nullable float array as its length and exact signed bit words. */
	private void bits(final NcrnaTextIO.Output out, final String name, final float[] values){
		assert(name!=null) : "Capture array records need stable field names for byte comparison";
		buffer.append(name).tab().append(values==null ? -1 : values.length);
		if(values!=null){for(final float value:values){buffer.tab().append(Float.floatToRawIntBits(value));}}
		write(out, buffer.nl());
	}

	/** The native writer copies this buffer before it is cleared for reuse. */
	private static void write(final NcrnaTextIO.Output out, final ByteBuilder buffer){
		assert(buffer.length()>0) : "Only framed capture records may enter the output stream";
		out.writer.print(buffer);
		buffer.clear();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final PrintStream outstream;
	private final String manifest;
	private final String out;
	private int rowLimit=32;
	/** Default preserves the seeded training control; false permits a one-network capture. */
	private boolean legacy=true;
	private final ByteBuilder buffer=new ByteBuilder();
	private static final String HEADER="id\tnet\tframe\tvectors\tinputs\toutputs\tprovenance\tseed";
}
