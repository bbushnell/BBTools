package ml;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Trains a feed-forward network for CONTINUOUS outputs and saves it in BBNet format
 * (via CellNet's own serializer, so the file is byte-compatible with everything that
 * consumes BBNets).
 * <p>
 * Differences from ml.Trainer: pure MSE regression (no balancing, cutoffs, FPR/FNR
 * machinery), Adam optimizer with cosine learning-rate decay, mini-batches, and input
 * standardization learned from the data and FOLDED INTO the first layer's weights at
 * export, so the saved net takes raw inputs.  Hidden layers are tanh; the output
 * activation is rslog (default), linear, or sigmoid, selected with final=.  Sigmoid
 * requires targets in [0,1]; normalize upstream as '#dims N 1' data always is.
 * After writing, the net is reloaded through CellNetParser and checked against the
 * in-memory model (round-trip verification).
 * <p>
 * Input is streamed with ByteFile and parsed with LineParser1, so memory scales with
 * the sample count rather than the file text. Each sample owns one input float array,
 * avoiding Java's single-array length limit for large training matrices.
 *
 * Usage: java ml.RegressionTrainer in=&lt;data.tsv&gt; out=&lt;net.bbnet&gt; dims=16,32,1
 *          [epochs=60] [batch=8192] [lr=0.003] [wd=1e-4] [seed=1] [vfraction=0.1]
 *          [valin=&lt;heldout.tsv&gt;] [simd=t] [threads=1]
 *          [final=rslog|linear|sigmoid] [netin=&lt;start.bbnet&gt;] [hidden=tanh,swish,...]
 *
 * hidden= sets the hidden-layer activations: one name is uniform, several are drawn per cell
 * at random from that set (seeded, so reproducible). Names: tanh sig rslog msig swish esig
 * emsig bell linear. Default 'tanh' keeps the hard-coded fast path and is bit-identical to
 * results recorded before this option existed. Derivatives use derivativeXFX(x,fx) because
 * derivativeFX alone is unimplemented for some functions (Swish), so pre-activations are
 * retained during the forward pass -- but only when a non-default activation is in use.
 * Output activation is separate and stays under final=.
 *
 * netin= continues training from an existing net instead of random initialization. In that
 * mode dims= may be omitted (the net supplies them) and training runs in RAW input space:
 * the saved net already has its original standardization folded into layer-1, so
 * re-standardizing would apply a second transform over one already baked in. mean=0/sd=1
 * makes both the forward transform and the export fold exact no-ops. The extraction is
 * verified against CellNet's own feedForward before training starts.
 *
 * Input format: same as train.sh — '#dims &lt;in&gt; &lt;out&gt;' header, then tab-delimited
 * floats, inputs first, then &lt;out&gt; target columns.  MULTIPLE OUTPUTS are supported: the
 * last dims entry sets the output width, loss is summed squared error across outputs, and
 * every output is checked on the round trip.  One output remains the common case and its
 * behaviour is unchanged.
 *
 * @author Amber (for Brian), with final= modes by Neptune
 * @date July 2026
 */
public class RegressionTrainer {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Program entry point.
	 * @param args Command line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		RegressionTrainer x=new RegressionTrainer(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructor; parses arguments, validates them, and prepares file formats.
	 * @param args Command line arguments
	 */
	public RegressionTrainer(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		parse(args);
		validateParams();
		checkFileExistence();

		ffin=FileFormat.testInput(in, FileFormat.TXT, null, true, true);
		ffval=(valIn==null ? null : FileFormat.testInput(valIn, FileFormat.TXT, null, true, true));
		ffout=FileFormat.testOutput(out, FileFormat.TXT, null, true, overwrite, false, false);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Parses command line arguments.
	 * @param args Command line arguments
	 */
	private void parse(String[] args){
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=", 2);
			final String a=split[0].toLowerCase();
			final String b=(split.length>1) ? split[1] : null;

			if(a.equals("in") || a.equals("input")){
				in=b;
			}else if(a.equals("valin") || a.equals("validationin")){
				valIn=b;
			}else if(a.equals("out") || a.equals("output")){
				out=b;
			}else if(a.equals("dims")){
				dims=Parse.parseIntArray(b, ",", "x");
			}else if(a.equals("epochs")){
				epochs=Integer.parseInt(b);
			}else if(a.equals("batch")){
				batch=Integer.parseInt(b);
			}else if(a.equals("lr") || a.equals("alpha")){
				lr=Double.parseDouble(b);
			}else if(a.equals("wd")){
				wd=Double.parseDouble(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("vfraction")){
				vfraction=Double.parseDouble(b);
				vfractionSet=true;
			}else if(a.equals("final")){
				finalType=parseFinalType(b);
			}else if(a.equals("hidden") || a.equals("hiddenfunctions")){
				hiddenSpec=b;
			}else if(a.equals("netin") || a.equals("innet")){
				netIn=b;
			}else if(a.equals("sparse")){
				useSparse=(b==null || Parse.parseBoolean(b));
			}else if(a.equals("normalize") || a.equals("normalizeweights")){
				normalize=(b==null || Parse.parseBoolean(b));
			}else if(a.equals("normfactor")){
				normFactor=Float.parseFloat(b);
			}else if(a.equals("normshrink")){
				normShrink=Float.parseFloat(b);
			}else if(a.equals("pad8") || a.equals("padlayers")){
				padLayers=(b==null || Parse.parseBoolean(b));
			}else if(a.equals("simd")){
				useSimd=(b==null || Parse.parseBoolean(b));
			}else if(a.equals("threads") || a.equals("t")){
				threads=Integer.parseInt(b);
			}else if(a.equals("sort")){
				sortSamples=(b==null || Parse.parseBoolean(b));
			}else if(a.equals("setsize") || a.equals("subsetsize")){
				setSize=Integer.parseInt(b);
			}else if(a.equals("density") || a.equals("maxdensity")){
				density=Float.parseFloat(b);
			}else if(a.equals("density1")){
				density1=Float.parseFloat(b);
			}else if(a.equals("edgeblock") || a.equals("edgeblocksize")){
				edgeBlockSize=Integer.parseInt(b);
			}else if(a.equals("checkexponents") || a.equals("checkexp")){
				checkExponents=(b==null || Parse.parseBoolean(b));
			}else if(a.equals("forcejavaparsedouble") || a.equals("fjpd")){
				Tools.FORCE_JAVA_PARSE_DOUBLE=Parse.parseBoolean(b);
			}else if(a.equals("ow") || a.equals("overwrite")){
				overwrite=Parse.parseBoolean(b);
			}else{
				outstream.println("Unknown parameter "+arg);
				assert(false) : "Unknown parameter "+arg;
				throw new IllegalArgumentException("Unknown parameter "+arg);
			}
		}
	}

	/** Translates a final= term into a Function constant selector.
	 * @param b Argument value
	 * @return FINAL_LINEAR, FINAL_RSLOG, or FINAL_SIGMOID */
	private static int parseFinalType(String b){
		if(b==null){throw new IllegalArgumentException("final requires a value");}
		if(b.equalsIgnoreCase("linear")){return FINAL_LINEAR;}
		if(b.equalsIgnoreCase("rslog")){return FINAL_RSLOG;}
		if(b.equalsIgnoreCase("sigmoid") || b.equalsIgnoreCase("sig")){return FINAL_SIGMOID;}
		throw new IllegalArgumentException("final must be linear, rslog, or sigmoid: "+b);
	}

	/** Validates parameter ranges and required arguments. */
	private void validateParams(){
		if(in==null || out==null){
			throw new IllegalArgumentException("required: in=, out=");
		}
		if(dims==null && netIn==null){
			throw new IllegalArgumentException("required: dims= (or netin=, which supplies them)");
		}
		if(vfraction<0 || vfraction>=1){
			throw new IllegalArgumentException("vfraction must be in [0,1): "+vfraction);
		}
		if(valIn!=null && vfractionSet && vfraction>0){
			throw new IllegalArgumentException("Use either valin= or vfraction>0, not both.");
		}
		if(threads<1){throw new IllegalArgumentException("threads must be positive: "+threads);}
		if(threads>1 && !useSimd){
			throw new IllegalArgumentException("threads>1 currently requires simd=t.");
		}
		if(dims==null){return;}//remaining checks run after the net supplies dims
		if(dims.length<2){
			throw new IllegalArgumentException("dims needs at least an input and output layer");
		}
		if(epochs<1){throw new IllegalArgumentException("epochs must be positive: "+epochs);}
		if(batch<1){throw new IllegalArgumentException("batch must be positive: "+batch);}
		if(padLayers){
			//Round hidden layers up to whole SIMD lanes so the vector kernels have no scalar
			//tail.  Input and output widths are fixed by the data and are left alone; the
			//added neurons are real and train normally, so this changes the model, not just
			//its layout, which is why it is opt-in.
			final StringBuilder sb=new StringBuilder();
			for(int i=1; i<dims.length-1; i++){
				final int padded=((dims[i]+LANE-1)/LANE)*LANE;
				if(padded!=dims[i]){sb.append(' ').append(dims[i]).append("->").append(padded);}
				dims[i]=padded;
			}
			if(sb.length()>0){outstream.println("pad8: hidden layers"+sb);}
		}
		if(density<=0 || density>1){
			throw new IllegalArgumentException("density must be in (0,1]: "+density);
		}
		if(density1>1){
			throw new IllegalArgumentException("density1 must be <=1: "+density1);
		}
		if(edgeBlockSize<1){
			throw new IllegalArgumentException("edgeblock must be positive: "+edgeBlockSize);
		}
	}

	/** Ensures input is readable and output is writable. */
	private void checkFileExistence(){
		if(!Tools.testOutputFiles(overwrite, false, false, out)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out+"\n");
		}
		if(!Tools.testInputFiles(false, true, in, valIn)){
			throw new RuntimeException("\nCan't read an input file.\n");
		}
		if(!Tools.testForDuplicateFiles(true, in, valIn, out)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Loads data, trains the network, exports it, and verifies the round trip.
	 * @param t Timer for tracking execution time
	 */
	void process(Timer t){
		if(netIn!=null){loadInitialNet();}
		loadData();
		standardize();
		shuffleAndSplit();
		allocateModel();
		if(netIn!=null){verifyLoadedNet();}
		train();
		exportNet();
		roundTripCheck();

		t.stop();
		outstream.println("Time: \t"+t);

		if(errorState){
			throw new RuntimeException(getClass().getName()
				+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads training data and, when supplied, a separate external validation set. */
	private void loadData(){
		numInputs=dims[0];
		numOutputs=dims[dims.length-1];
		trainingData=loadData(ffin, in);
		numSamples=trainingData.size();
		validationData=(ffval==null ? null : loadData(ffval, valIn));

		outstream.println("loaded "+numSamples+" training samples"
			+(validationData==null ? "" : ", "+validationData.size()+" external validation samples")
			+", dims="+Arrays.toString(dims));
		if(numSamples<1){throw new RuntimeException("No samples were loaded from "+in);}
		if(validationData!=null && validationData.size()<1){
			throw new RuntimeException("No samples were loaded from "+valIn);
		}
	}

	/**
	 * Loads one data file, retrying with exact float parsing if requested exponent
	 * detection finds scientific notation.
	 * @param ff Input format
	 * @param name Input path, used in diagnostics
	 * @return Parsed regression data
	 */
	private RegressionData loadData(final FileFormat ff, final String name){
		RegressionData data=readVectors(ff, name);
		if(data==null){
			//Parse's fast path does not support exponents: outside of assertions it reads
			//'e' as a digit and silently returns a wrong number, so re-read exactly.
			Tools.FORCE_JAVA_PARSE_DOUBLE=true;
			outstream.println("Note: exponent notation detected in "+name
				+"; re-reading with exact float parsing.");
			data=readVectors(ff, name);
		}
		return data;
	}

	/**
	 * Streams one input file into row-oriented input and target arrays.
	 * @param ff Input format
	 * @param name Input path, used in diagnostics
	 * @return Parsed data, or null if exact parsing must be enabled and retried
	 */
	private RegressionData readVectors(final FileFormat ff, final String name){
		final ArrayList<float[]> vectorList=new ArrayList<float[]>();
		final ArrayList<float[]> targetList=new ArrayList<float[]>();
		final LineParser1 lp=new LineParser1('\t');

		//Hoisted out of the loop so the common case costs one test, not one per line
		final boolean check=checkExponents && !Tools.FORCE_JAVA_PARSE_DOUBLE;

		final ByteFile bf=ByteFile.makeByteFile(ff);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length==0){continue;}
			if(line[0]=='#'){
				if(Tools.startsWith(line, "#dims")){checkHeader(line, lp, name);}
				continue;
			}
			if(check && hasExponent(line)){
				errorState|=bf.close();
				return null;
			}
			lp.set(line);
			if(lp.terms()<numInputs+numOutputs){continue;}
			final float[] vector=new float[numInputs];
			final float[] target=new float[numOutputs];
			for(int i=0; i<numInputs; i++){vector[i]=lp.parseFloat(i);}
			for(int k=0; k<numOutputs; k++){target[k]=lp.parseFloat(numInputs+k);}
			vectorList.add(vector);
			targetList.add(target);
		}
		errorState|=bf.close();
		return new RegressionData(vectorList, targetList);
	}

	/** @param line A data line
	 * @return True if the line contains exponent notation */
	private static boolean hasExponent(final byte[] line){
		for(int i=0; i<line.length; i++){
			final byte c=line[i];
			if(c=='e' || c=='E'){return true;}
		}
		return false;
	}

	/**
	 * Warns if the file's own '#dims' header disagrees with dims=.  The header is
	 * advisory here: dims= remains authoritative, because it also carries the hidden
	 * layer sizes.  A mismatch in the input count means the columns are being read
	 * with the wrong stride, which produces a silently wrong net.
	 * @param line The '#dims' header line
	 * @param lp Parser to use
	 * @param source Input path, used in warnings
	 */
	private void checkHeader(byte[] line, LineParser1 lp, String source){
		lp.set(line);
		if(lp.terms()<2){return;}
		if(lp.terms()>=3){
			final int headerOutputs=lp.parseInt(2);
			if(headerOutputs!=numOutputs){
				outstream.println("WARNING: "+source+" header declares "+headerOutputs
					+" outputs but dims= declares "+numOutputs+".");
			}
		}
		final int headerInputs=lp.parseInt(1);
		if(headerInputs!=numInputs){
			outstream.println("WARNING: "+source+" header declares "+headerInputs
				+" inputs but dims= declares "+numInputs
				+"; columns will be read with the dims= stride.");
		}
	}

	/** Computes per-column mean and standard deviation of the inputs. */
	private void standardize(){
		mean=new double[numInputs];
		sd=new double[numInputs];
		if(netIn!=null){
			//Continuing from an existing net, so train in RAW input space.
			//The saved net already has the ORIGINAL standardization folded into its layer-1
			//weights, and those folded weights are what we just loaded.  Re-standardizing here
			//would apply a second transform on top of one already baked in, and the export fold
			//would apply a third -- silently, since the round-trip check compares the written net
			//against the in-memory model and both would be wrong identically.
			//mean=0, sd=1 makes both the forward transform and the export fold exact no-ops.
			Arrays.fill(sd, 1);
			return;
		}
		for(float[] vector : trainingData.inputs){
			for(int i=0; i<numInputs; i++){mean[i]+=vector[i];}
		}
		for(int i=0; i<numInputs; i++){mean[i]/=numSamples;}
		for(float[] vector : trainingData.inputs){
			for(int i=0; i<numInputs; i++){
				final double d=vector[i]-mean[i];
				sd[i]+=d*d;
			}
		}
		for(int i=0; i<numInputs; i++){sd[i]=Math.max(1e-9, Math.sqrt(sd[i]/numSamples));}
	}

	/** Shuffles training indices and configures internal or external validation. */
	private void shuffleAndSplit(){
		randy=new Random(seed);
		perm=new int[numSamples];
		for(int i=0; i<numSamples; i++){perm[i]=i;}
		for(int i=numSamples-1; i>0; i--){
			final int j=randy.nextInt(i+1);
			final int t=perm[i]; perm[i]=perm[j]; perm[j]=t;
		}
		if(validationData==null){
			numValid=(int)(numSamples*vfraction);
			numTrain=numSamples-numValid;
		}else{
			numTrain=numSamples;
			numValid=validationData.size();
		}
		if(numTrain<1){throw new RuntimeException("vfraction leaves no training samples.");}
	}

	/**
	 * Allocates weights, biases, Adam state, and activation scratch.
	 * Each layer's weights are one flat row-major array of size rows*cols, indexed
	 * [i*cols+j], which keeps a neuron's inputs contiguous.
	 */
	private void allocateModel(){
		layers=dims.length-1;
		weights=new double[layers][];
		bias=new double[layers][];
		mW=new double[layers][]; vW=new double[layers][];
		mB=new double[layers][]; vB=new double[layers][];
		gW=new double[layers][]; gB=new double[layers][];

		for(int l=0; l<layers; l++){
			final int rows=dims[l+1], cols=dims[l];
			weights[l]=new double[rows*cols];
			bias[l]=new double[rows];
			mW[l]=new double[rows*cols]; vW[l]=new double[rows*cols];
			mB[l]=new double[rows]; vB[l]=new double[rows];
			gW[l]=new double[rows*cols]; gB[l]=new double[rows];

			if(netIn!=null){//continue from the loaded net rather than random init
				final Cell[] layer=net.net[l+1];
				for(int i=0, p=0; i<rows; i++, p+=cols){
					bias[l][i]=layer[i].bias();
					final float[] w=layer[i].weights;
					for(int j=0; j<cols; j++){weights[l][p+j]=w[j];}
				}
			}else{
				final double scale=1.0/Math.sqrt(cols);
				for(int i=0, p=0; i<rows; i++){
					for(int j=0; j<cols; j++, p++){weights[l][p]=randy.nextGaussian()*scale;}
				}
			}
		}

		act=new double[dims.length][];
		preAct=new double[dims.length][];
		delta=new double[dims.length][];
		for(int i=0; i<dims.length; i++){
			act[i]=new double[dims[i]];
			preAct[i]=new double[dims[i]];
			delta[i]=new double[dims[i]];
		}
		scratchA=new double[maxDim()];
		scratchB=new double[maxDim()];

		assignHiddenFunctions();
		buildNetAndMask();
	}

	/** Maps an activation name to its Function type constant.
	 * @param name Activation name
	 * @return One of the Function type constants */
	private static int functionType(String name){
		final String n=name.trim().toLowerCase();
		if(n.equals("tanh")){return Function.TANH;}
		if(n.equals("sig") || n.equals("sigmoid")){return Function.SIG;}
		if(n.equals("rslog")){return Function.RSLOG;}
		if(n.equals("msig")){return Function.MSIG;}
		if(n.equals("swish")){return Function.SWISH;}
		if(n.equals("esig")){return Function.ESIG;}
		if(n.equals("emsig")){return Function.EMSIG;}
		if(n.equals("bell")){return Function.BELL;}
		if(n.equals("linear")){return Function.LINEAR;}
		throw new IllegalArgumentException("Unknown hidden activation: "+name);
	}

	/**
	 * Assigns hidden-layer activations.  A single name leaves hiddenFunc null so the
	 * hard-coded tanh path runs — that keeps the default both fast and bit-identical to
	 * every result recorded before this option existed.  Several names mean each hidden cell
	 * draws one at random from the set, seeded, so a mixed net is reproducible.
	 * <p>
	 * Output-layer activation is NOT touched here; it stays under final=.
	 */
	private void assignHiddenFunctions(){
		final String[] names=hiddenSpec.split(",");
		final int[] types=new int[names.length];
		for(int i=0; i<names.length; i++){types[i]=functionType(names[i]);}

		if(types.length==1 && types[0]==Function.TANH){return;}//fast path, bit-exact

		hiddenFunc=new Function[dims.length][];
		hiddenType=new int[dims.length][];
		final int[] counts=new int[Function.LINEAR+1];
		for(int l=1; l<layers; l++){//hidden activation layers only
			hiddenFunc[l]=new Function[dims[l]];
			hiddenType[l]=new int[dims[l]];
			for(int i=0; i<dims[l]; i++){
				final int t=(types.length==1) ? types[0] : types[randy.nextInt(types.length)];
				hiddenType[l][i]=t;
				hiddenFunc[l][i]=Function.getFunction(t);
				counts[t]++;
			}
		}
		final StringBuilder sb=new StringBuilder("hidden activations:");
		for(int t=0; t<counts.length; t++){
			if(counts[t]>0){sb.append(' ').append(Function.getFunction(t).name())
				.append('=').append(counts[t]);}
		}
		outstream.println(sb);
	}

	/**
	 * Loads an existing net to continue training from, instead of random initialization.
	 * Supplies dims if they were not given and validates them if they were.
	 */
	private void loadInitialNet(){
		net=CellNetParser.load(netIn);
		if(net==null){throw new RuntimeException("Could not load a net from "+netIn);}
		if(padLayers){
			throw new IllegalArgumentException("pad8 cannot be combined with netin: padding "
				+"would change the architecture away from the loaded net");
		}
		final int[] netDims=new int[net.net.length];
		for(int i=0; i<netDims.length; i++){netDims[i]=net.net[i].length;}
		if(dims==null){
			dims=netDims;
		}else if(!Arrays.equals(dims, netDims)){
			throw new IllegalArgumentException("dims="+Arrays.toString(dims)
				+" does not match the net in "+netIn+" "+Arrays.toString(netDims));
		}
		outstream.println("netin: continuing from "+netIn+" "+Arrays.toString(dims)
			+"; training in RAW input space, no re-standardization");
	}

	/**
	 * Confirms the weights were extracted correctly, by checking that this trainer's own
	 * forward pass reproduces the loaded net's output on real samples.
	 * <p>
	 * This is the check that matters for netin, and it is NOT the same as the export
	 * round-trip: that one compares the written net against the in-memory model, so if both
	 * are wrong in the same way it passes anyway. This one compares against an INDEPENDENT
	 * implementation — CellNet's own feedForward — so a bad extraction cannot hide.
	 */
	private void verifyLoadedNet(){
		final int n=Math.min(1000, numSamples);
		double maxDiff=0;
		for(int s=0; s<n; s++){
			final float[] vector=trainingData.inputs[s];
			net.applyInput(vector);
			net.feedForward();
			final double[] b=predict(vector, weights, bias);
			for(int k=0; k<numOutputs; k++){
				maxDiff=Math.max(maxDiff, Math.abs(net.getOutput(k)-b[k]));
			}
		}
		final boolean ok=(maxDiff<1e-3);
		outstream.println(String.format(
			"netin load check vs source net: maxDiff=%.2e over %d samples %s",
			maxDiff, n, (ok ? "(OK)" : "(FAILED - weights were not extracted correctly)")));
		if(!ok){
			throw new RuntimeException("netin verification failed; refusing to train from a "
				+"model that does not reproduce the net it was loaded from.");
		}
	}

	/**
	 * Creates the CellNet that will later be exported and, when density is below 1,
	 * derives the edge mask from the topology CellNet generated for itself rather than
	 * reimplementing its edge selection.  Taking the mask from the net is what keeps the
	 * trained model and the written file agreeing about which edges exist.
	 * Masked-out weights are zeroed here and never updated afterward.
	 */
	private void buildNetAndMask(){
		if(netIn==null){
			//randomize() consults the activation-rate tables that ml.Trainer normally
			//initializes; every cell's function is assigned explicitly at export, so zero
			//the rates (skips random activation selection and satisfies the assertion).
			Arrays.fill(Function.TYPE_RATES, 0f);
			net=new CellNet(dims, Math.max(0, seed), density, density1, edgeBlockSize,
				new ArrayList<String>());
			net.randomize();//allocates edge topology + weight arrays
		}//else the net was already loaded by loadInitialNet(), and IS the export target

		//With netin, any weight the source net left at exactly 0 is a pruned edge and must
		//stay pruned, exactly as it would from density=.
		if(netIn==null && density>=1 && density1<=0){return;}

		edgeMask=new boolean[layers][];
		long kept=0, total=0;
		for(int l=0; l<layers; l++){
			final int rows=dims[l+1], cols=dims[l];
			final boolean[] mask=new boolean[rows*cols];
			final Cell[] layer=net.net[l+1];
			final double[] w=weights[l];
			for(int i=0, p=0; i<rows; i++, p+=cols){
				final float[] netWeights=layer[i].weights;
				for(int j=0; j<cols; j++){
					//makeEdgesDense leaves unselected edges at exactly 0 and randomizes the rest
					final boolean on=(netWeights[j]!=0);
					mask[p+j]=on;
					if(on){kept++;}else{w[p+j]=0;}
					total++;
				}
			}
			edgeMask[l]=mask;
		}
		outstream.println("density="+density+": kept "+kept+" of "+total+" edges ("
			+String.format("%.1f%%", 100.0*kept/Math.max(1, total))+")");
	}

	/** @return The largest layer width, for scratch allocation */
	private int maxDim(){
		int max=0;
		for(int d : dims){max=Math.max(max, d);}
		return max;
	}

	/** Runs the training loop, retaining the best-validating weights. */
	private void train(){
		final int[] order=new int[numTrain];
		for(int i=0; i<numTrain; i++){order[i]=perm[i];}

		if(sortSamples){
			//Unseen samples start at the same default error magnitude ml.Sample uses, so
			//nothing has been "learned" until it has actually been visited once.
			sampleError=new float[numSamples];
			sampleEpoch=new int[numSamples];
			Arrays.fill(sampleError, 1f);
			if(setSize<1 || setSize>numTrain){setSize=numTrain;}
			outstream.println("sort=t: "+setSize+" of "+numTrain+" training samples per epoch");
		}

		if(useSimd){
			setupFloat();
			setupGradientWorkers();
			outstream.println("simd=t: training in float (Shared.SIMD="+shared.Shared.SIMD
				+"); results will be close to but not identical to the scalar path");
		}

		long step=0;
		bestValid=Double.MAX_VALUE;

		try{
			for(int ep=1; ep<=epochs; ep++){
				final int active=(sortSamples ? prioritize(order, ep) : numTrain);
				if(!sortSamples){
					for(int i=numTrain-1; i>0; i--){
						final int j=randy.nextInt(i+1);
						final int t=order[i]; order[i]=order[j]; order[j]=t;
					}
				}
				final double lrNow=lr*0.5*(1+Math.cos(Math.PI*(ep-1)/epochs));
				double trainMse=0;

				for(int start=0; start<active; start+=batch){
					final int end=Math.min(active, start+batch);
					final int bs=end-start;
					if(useSimd){
						if(gradientPool==null){//Exact legacy order for the default one-thread path
							zeroGradientsF();
							for(int s=start; s<end; s++){
								trainMse+=accumulateGradientF(order[s], actF, preActF,
									deltaF, gwF, gbF);
							}
						}else{
							trainMse+=accumulateParallelBatchF(order, start, end);
						}
					}else{
						zeroGradients();
						for(int s=start; s<end; s++){
							trainMse+=accumulateGradient(order[s]);
						}
					}

					step++;
					if(useSimd){applyAdamF(bs, step, lrNow);}else{applyAdam(bs, step, lrNow);}
					if(normalize && normFactor>1e-4f){normalizeWeights();}
				}

				if(useSimd){syncFloatToDouble();}
				final double validMse=validate();
				if(validMse<bestValid){
					bestValid=validMse;
					bestWeights=deepCopy(weights);
					bestBias=deepCopy(bias);
					//Checkpoint: persist the improved net immediately, so a job killed mid-run
					//(node failure, preemption) doesn't lose all epochs of progress -- only
					//whatever improvement happened since the last checkpoint (Brian, 2026-08-27).
					//weights/bias at this point ARE the just-improved values exportNet reads.
					exportNet();
				}
				if(ep%5==0 || ep==1 || ep==epochs){
					outstream.println(String.format(
						"epoch %d lr=%.5f trainMSE=%.6f valMSE=%.6f best=%.6f",
						ep, lrNow, trainMse/Math.max(1, active), validMse, bestValid));
				}
			}
		}finally{
			shutdownGradientWorkers();
		}

		if(bestWeights!=null){weights=bestWeights; bias=bestBias;}
	}

	/**
	 * Orders the training indices by priority and returns how many to train this epoch.
	 * <p>
	 * Priority follows ml.Sample.calcPivot(): error magnitude, minus a decay in the epoch
	 * the sample was last visited, sorted descending.  The decay term is what makes this
	 * safe — without it a sample the net happens to fit early would never be revisited and
	 * could drift arbitrarily far out of date.  The excess-error branch in calcPivot() is
	 * deliberately not ported: it distinguishes overshooting a positive label from
	 * undershooting it, which has no meaning for a continuous target.
	 * <p>
	 * Sorting is done on packed longs (sortable float key in the high half, index in the
	 * low half) to avoid boxing every index.  The selected prefix is then shuffled, so
	 * batches stay mixed rather than being ordered hardest-first inside the batch.
	 *
	 * @param order Training indices, reordered in place
	 * @param epoch Current epoch, 1-based
	 * @return Number of leading entries of order to train on
	 */
	private int prioritize(final int[] order, final int epoch){
		final long[] keys=new long[numTrain];
		for(int i=0; i<numTrain; i++){
			final int s=order[i];
			final float pivot=sampleError[s]-sampleEpoch[s]*EPOCH_MULT;
			//Negate for descending, then map to a signed-order-preserving int
			keys[i]=(((long)sortableInt(-pivot))<<32) | (s & 0xFFFFFFFFL);
		}
		Arrays.sort(keys);
		for(int i=0; i<numTrain; i++){order[i]=(int)(keys[i] & 0xFFFFFFFFL);}

		final int active=Math.min(setSize, numTrain);
		for(int i=active-1; i>0; i--){//shuffle only the selected prefix
			final int j=randy.nextInt(i+1);
			final int t=order[i]; order[i]=order[j]; order[j]=t;
		}
		for(int i=0; i<active; i++){sampleEpoch[order[i]]=epoch;}
		return active;
	}

	/** Maps a float onto an int whose signed ordering matches the float's ordering.
	 * @param f Value to map
	 * @return Order-preserving int */
	private static int sortableInt(final float f){
		final int i=Float.floatToIntBits(f);
		return (i>=0) ? i : (i ^ 0x7FFFFFFF);
	}

	/**
	 * Builds float mirrors of the model for the vectorized path.
	 * Each neuron gets its own contiguous float[] because simd.Vector.fma operates on whole
	 * arrays with no offset, and a transposed copy is kept per layer because the backward
	 * pass reads the weight matrix column-wise, which is neither contiguous nor vectorizable
	 * in row-major order.
	 */
	private void setupFloat(){
		weightsF=new float[layers][][]; gwF=new float[layers][][]; wT=new float[layers][][];
		mWF=new float[layers][][]; vWF=new float[layers][][];
		biasF=new float[layers][]; gbF=new float[layers][];
		mBF=new float[layers][]; vBF=new float[layers][];
		for(int l=0; l<layers; l++){
			final int rows=dims[l+1], cols=dims[l];
			weightsF[l]=new float[rows][cols]; gwF[l]=new float[rows][cols];
			mWF[l]=new float[rows][cols]; vWF[l]=new float[rows][cols];
			biasF[l]=new float[rows]; gbF[l]=new float[rows];
			mBF[l]=new float[rows]; vBF[l]=new float[rows];
			wT[l]=new float[cols][rows];
			final double[] w=weights[l], b=bias[l];
			for(int i=0, p=0; i<rows; i++, p+=cols){
				biasF[l][i]=(float)b[i];
				for(int j=0; j<cols; j++){weightsF[l][i][j]=(float)w[p+j];}
			}
		}
		actF=new float[dims.length][]; deltaF=new float[dims.length][];
		preActF=new float[dims.length][];
		for(int i=0; i<dims.length; i++){actF[i]=new float[dims[i]]; deltaF[i]=new float[dims[i]];
			preActF[i]=new float[dims[i]];}
		refreshTranspose();
		if(useSparse && edgeMask!=null){buildSparse();}
	}

	/** Allocates a fixed worker pool and one private scratch/gradient set per worker. */
	private void setupGradientWorkers(){
		gradientThreadCount=Math.min(threads, Math.min(batch, numTrain));
		if(gradientThreadCount<=1){return;}
		gradientWorkers=new ArrayList<GradWorker>(gradientThreadCount);
		gradientWorkers.add(new GradWorker(actF, preActF, deltaF, gwF, gbF));
		for(int i=1; i<gradientThreadCount; i++){
			gradientWorkers.add(new GradWorker(newFloatLayerBuffers(), newFloatLayerBuffers(),
				newFloatLayerBuffers(), newFloatWeightBuffers(), newFloatBiasBuffers()));
		}
		gradientPool=Executors.newFixedThreadPool(gradientThreadCount);
		outstream.println("threads="+gradientThreadCount
			+": deterministic per-worker SIMD gradients; Adam and reduction remain serial");
	}

	/** @return Fresh per-layer activation or delta buffers. */
	private float[][] newFloatLayerBuffers(){
		final float[][] x=new float[dims.length][];
		for(int i=0; i<dims.length; i++){x[i]=new float[dims[i]];}
		return x;
	}

	/** @return Fresh per-layer weight-gradient buffers. */
	private float[][][] newFloatWeightBuffers(){
		final float[][][] x=new float[layers][][];
		for(int l=0; l<layers; l++){x[l]=new float[dims[l+1]][dims[l]];}
		return x;
	}

	/** @return Fresh per-layer bias-gradient buffers. */
	private float[][] newFloatBiasBuffers(){
		final float[][] x=new float[layers][];
		for(int l=0; l<layers; l++){x[l]=new float[dims[l+1]];}
		return x;
	}

	/**
	 * Builds the index sets for sparse forward and backward passes, once.
	 * Only the VALUES change as training proceeds, so the indices are built here and only
	 * refreshed values are copied afterwards.
	 * <p>
	 * Index sets are built in ascending order deliberately. simd/Vector#002 documents that
	 * the sparse fma takes a dense shortcut ignoring bSet when a.length==b.length; with
	 * ascending indices, equal lengths mean bSet is the identity permutation, which is
	 * exactly the case where that shortcut is correct.
	 */
	private void buildSparse(){
		spWIdx=new int[layers][][]; spTIdx=new int[layers][][];
		spW=new float[layers][][]; spT=new float[layers][][];
		long keptFwd=0;
		for(int l=0; l<layers; l++){
			final int rows=dims[l+1], cols=dims[l];
			final boolean[] mask=edgeMask[l];
			spWIdx[l]=new int[rows][]; spW[l]=new float[rows][];
			for(int i=0, p=0; i<rows; i++, p+=cols){
				int n=0;
				for(int j=0; j<cols; j++){if(mask[p+j]){n++;}}
				final int[] idx=new int[n];
				for(int j=0, k=0; j<cols; j++){if(mask[p+j]){idx[k++]=j;}}
				spWIdx[l][i]=idx; spW[l][i]=new float[n];
				keptFwd+=n;
			}
			spTIdx[l]=new int[cols][]; spT[l]=new float[cols][];
			for(int j=0; j<cols; j++){
				int n=0;
				for(int i=0; i<rows; i++){if(mask[i*cols+j]){n++;}}
				final int[] idx=new int[n];
				for(int i=0, k=0; i<rows; i++){if(mask[i*cols+j]){idx[k++]=i;}}
				spTIdx[l][j]=idx; spT[l][j]=new float[n];
			}
		}
		refreshSparse();
		outstream.println("sparse kernels: "+keptFwd+" live edges; block="+edgeBlockSize
			+((edgeBlockSize&7)==0 ? " (vectorized)" : " (scalar gather; use edgeblock=8 to vectorize)"));
	}

	/** Copies current weight values into the sparse views; the index sets never change. */
	private void refreshSparse(){
		for(int l=0; l<layers; l++){
			final int rows=dims[l+1], cols=dims[l];
			final float[][] w=weightsF[l];
			for(int i=0; i<rows; i++){
				final int[] idx=spWIdx[l][i]; final float[] dst=spW[l][i]; final float[] wi=w[i];
				for(int k=0; k<idx.length; k++){dst[k]=wi[idx[k]];}
			}
			for(int j=0; j<cols; j++){
				final int[] idx=spTIdx[l][j]; final float[] dst=spT[l][j];
				for(int k=0; k<idx.length; k++){dst[k]=w[idx[k]][j];}
			}
		}
	}

	/** Rebuilds the transposed weight views after the weights change. */
	private void refreshTranspose(){
		for(int l=0; l<layers; l++){
			final int rows=dims[l+1], cols=dims[l];
			final float[][] w=weightsF[l], t=wT[l];
			for(int i=0; i<rows; i++){
				final float[] wi=w[i];
				for(int j=0; j<cols; j++){t[j][i]=wi[j];}
			}
		}
	}

	/**
	 * Vectorized forward and backward pass for one sample.
	 * @param sample Index of the sample
	 * @return Squared error for this sample
	 */
	private double accumulateGradientF(final int sample, final float[][] actLocal,
			final float[][] preActLocal, final float[][] deltaLocal,
			final float[][][] gwLocal, final float[][] gbLocal){
		final float[] vector=trainingData.inputs[sample];
		final float[] target=trainingData.targets[sample];
		final float[] a0=actLocal[0];
		for(int i=0; i<numInputs; i++){a0[i]=(float)((vector[i]-mean[i])/sd[i]);}

		for(int l=0; l<layers; l++){
			final int rows=dims[l+1];
			final float[][] w=weightsF[l];
			final float[] b=biasF[l], prev=actLocal[l], next=actLocal[l+1];
			if(spW==null){
				for(int i=0; i<rows; i++){
					final float z=b[i]+simd.Vector.fma(prev, w[i]);
					if(hiddenFunc!=null){preActLocal[l+1][i]=z;}
					next[i]=(l==layers-1) ? (float)fin(z, finalType)
						: (float)(hiddenFunc==null ? Math.tanh(z) : hiddenFunc[l+1][i].activate(z));
				}
			}else{//gather only the live edges instead of multiplying through the pruned zeros
				final float[][] sw=spW[l]; final int[][] si=spWIdx[l];
				for(int i=0; i<rows; i++){
					final float z=b[i]+simd.Vector.fma(sw[i], prev, si[i], edgeBlockSize, true);
					if(hiddenFunc!=null){preActLocal[l+1][i]=z;}
					next[i]=(l==layers-1) ? (float)fin(z, finalType)
						: (float)(hiddenFunc==null ? Math.tanh(z) : hiddenFunc[l+1][i].activate(z));
				}
			}
		}

		final float[] outs=actLocal[layers];
		double sqErr=0, absErr=0;
		for(int k=0; k<numOutputs; k++){
			final double outv=outs[k];
			final double err=outv-target[k];
			sqErr+=err*err;
			absErr+=Math.abs(err);
			deltaLocal[layers][k]=(float)(2*err*finDeriv(outv, finalType));
		}
		if(sampleError!=null){sampleError[sample]=(float)(absErr/numOutputs);}

		for(int l=layers-1; l>=0; l--){
			final int rows=dims[l+1], cols=dims[l];
			final float[] prev=actLocal[l], dNext=deltaLocal[l+1], gb=gbLocal[l];
			final float[][] gw=gwLocal[l];
			for(int i=0; i<rows; i++){
				final float d=dNext[i];
				gb[i]+=d;
				simd.Vector.addProduct(gw[i], prev, d);//gw[i] += prev*d
			}
			if(l>0){
				final float[] dCur=deltaLocal[l];
				if(spT==null){
					final float[][] t=wT[l];
					for(int j=0; j<cols; j++){
						dCur[j]=(float)(simd.Vector.fma(dNext, t[j])*(hiddenFunc==null
							? (1-prev[j]*prev[j])
							: hiddenFunc[l][j].derivativeXFX(preActLocal[l][j], prev[j])));
					}
				}else{
					final float[][] st=spT[l]; final int[][] si=spTIdx[l];
					for(int j=0; j<cols; j++){
						dCur[j]=(float)(simd.Vector.fma(st[j], dNext, si[j], edgeBlockSize, true)
							*(hiddenFunc==null ? (1-prev[j]*prev[j])
								: hiddenFunc[l][j].derivativeXFX(preActLocal[l][j], prev[j])));
					}
				}
			}
		}
		return sqErr;
	}

	/** Runs one batch on fixed worker partitions and reduces gradients in worker order. */
	private double accumulateParallelBatchF(final int[] order, final int start, final int end){
		final int activeWorkers=Math.min(gradientThreadCount, end-start);
		final int base=(end-start)/activeWorkers, extra=(end-start)%activeWorkers;
		int from=start;
		for(int i=0; i<activeWorkers; i++){
			final int to=from+base+(i<extra ? 1 : 0);
			gradientWorkers.get(i).prepare(order, from, to);
			from=to;
		}

		final List<Future<Void>> futures;
		try{
			futures=gradientPool.invokeAll(gradientWorkers.subList(0, activeWorkers));
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
			throw new RuntimeException("Interrupted while accumulating a training batch.", e);
		}
		for(Future<Void> future : futures){
			try{future.get();}
			catch(InterruptedException e){
				Thread.currentThread().interrupt();
				throw new RuntimeException("Interrupted while joining a training batch.", e);
			}catch(ExecutionException e){
				final Throwable cause=e.getCause();
				if(cause instanceof Error){throw (Error)cause;}
				if(cause instanceof RuntimeException){throw (RuntimeException)cause;}
				throw new RuntimeException("Gradient worker failed.", cause);
			}
		}

		double sqErr=0;
		for(int i=0; i<activeWorkers; i++){sqErr+=gradientWorkers.get(i).sqErr;}
		for(int i=1; i<activeWorkers; i++){addGradientsF(gradientWorkers.get(i));}
		return sqErr;
	}

	/** Adds one worker's private gradients into the master accumulators. */
	private void addGradientsF(final GradWorker worker){
		for(int l=0; l<layers; l++){
			for(int i=0; i<gwF[l].length; i++){simd.Vector.add(gwF[l][i], worker.gw[l][i]);}
			simd.Vector.add(gbF[l], worker.gb[l]);
		}
	}

	/** Stops the reusable pool on normal completion and every failure path. */
	private void shutdownGradientWorkers(){
		if(gradientPool!=null){
			gradientPool.shutdownNow();
			gradientPool=null;
		}
	}

	/** Float Adam update, then refresh the transpose the backward pass reads.
	 * @param bs Batch size
	 * @param step Global step count
	 * @param lrNow Current learning rate */
	private void applyAdamF(final int bs, final long step, final double lrNow){
		final float c1=(float)(1-Math.pow(BETA1, step)), c2=(float)(1-Math.pow(BETA2, step));
		final float b1=(float)BETA1, b2=(float)BETA2, eps=(float)EPS, lrF=(float)lrNow, wdF=(float)wd;
		for(int l=0; l<layers; l++){
			final int rows=dims[l+1], cols=dims[l];
			final boolean[] mask=(edgeMask==null) ? null : edgeMask[l];
			final float[][] w=weightsF[l], gw=gwF[l], mw=mWF[l], vw=vWF[l];
			final float[] b=biasF[l], gb=gbF[l], mb=mBF[l], vb=vBF[l];
			for(int i=0, p=0; i<rows; i++, p+=cols){
				final float[] wi=w[i], gi=gw[i], mi=mw[i], vi=vw[i];
				for(int j=0; j<cols; j++){
					if(mask!=null && !mask[p+j]){continue;}
					final float g=gi[j]/bs+wdF*wi[j];
					mi[j]=b1*mi[j]+(1-b1)*g;
					vi[j]=b2*vi[j]+(1-b2)*g*g;
					wi[j]-=lrF*(mi[j]/c1)/((float)Math.sqrt(vi[j]/c2)+eps);
				}
				final float g=gb[i]/bs;
				mb[i]=b1*mb[i]+(1-b1)*g;
				vb[i]=b2*vb[i]+(1-b2)*g*g;
				b[i]-=lrF*(mb[i]/c1)/((float)Math.sqrt(vb[i]/c2)+eps);
			}
		}
		//Only the values moved; sparse index sets are structural and never change.
		if(spW==null){refreshTranspose();}else{refreshSparse();}
	}

	/**
	 * Z-score weight standardization, ported from CellNet.normalizeDense().
	 * Statistics are taken over NON-ZERO weights only and zeros are left untouched, so a
	 * density mask survives normalization; the result is blended with the original by
	 * normFactor rather than replacing it, and magnitudes are floored to avoid underflow.
	 * The blend strength decays after each application, so normalization is strong early
	 * and fades as training settles.
	 */
	private void normalizeWeights(){
		for(int l=0; l<layers; l++){
			if(useSimd){normalizeLayerF(weightsF[l]);}
			else{normalizeLayer(weights[l]);}
		}
		if(useSimd){refreshTranspose();}
		normFactor*=normShrink;
	}

	/** @param w One layer's weights, flat row-major */
	private void normalizeLayer(final double[] w){
		double sum=0, sumSq=0; int n=0;
		for(double f : w){if(f!=0){sum+=f; sumSq+=f*f; n++;}}
		if(n<2){return;}
		final double mean=sum/n;
		double stdev=Math.sqrt(Math.max(0, sumSq/n-mean*mean));
		if(stdev<0.001){stdev=0.001;}
		final double invStdev=1/stdev, keep=1-normFactor;
		for(int i=0; i<w.length; i++){
			final double f=w[i];
			if(f!=0){
				final double target=(f-mean)*invStdev;
				final double v=keep*f+normFactor*target;
				w[i]=(Math.abs(v)<1e-12) ? 1e-12 : v;
			}
		}
	}

	/** @param layer One layer's weights, one row per neuron; statistics span the whole layer
	 * exactly as CellNet.normalizeDense does, not each neuron separately */
	private void normalizeLayerF(final float[][] layer){
		double sum=0, sumSq=0; int n=0;
		for(float[] w : layer){
			for(float f : w){if(f!=0){sum+=f; sumSq+=f*f; n++;}}
		}
		if(n<2){return;}
		final double mean=sum/n;
		double stdev=Math.sqrt(Math.max(0, sumSq/n-mean*mean));
		if(stdev<0.001){stdev=0.001;}
		final double invStdev=1/stdev, keep=1-normFactor;
		for(float[] w : layer){
			for(int i=0; i<w.length; i++){
				final float f=w[i];
				if(f!=0){
					final double target=(f-mean)*invStdev;
					final double v=keep*f+normFactor*target;
					w[i]=(float)((Math.abs(v)<1e-12) ? 1e-12 : v);
				}
			}
		}
	}

	/** Clears the float gradient accumulators. */
	private void zeroGradientsF(){
		zeroGradientsF(gwF, gbF);
	}

	/** Clears one worker's float gradient accumulators. */
	private void zeroGradientsF(final float[][][] gw, final float[][] gb){
		for(int l=0; l<layers; l++){
			for(float[] row : gw[l]){Arrays.fill(row, 0f);}
			Arrays.fill(gb[l], 0f);
		}
	}

	/** Copies the float model back into the double arrays, so validation, export and the
	 * round-trip check all keep running on one code path regardless of how training ran. */
	private void syncFloatToDouble(){
		for(int l=0; l<layers; l++){
			final int rows=dims[l+1], cols=dims[l];
			final double[] w=weights[l], b=bias[l];
			for(int i=0, p=0; i<rows; i++, p+=cols){
				b[i]=biasF[l][i];
				for(int j=0; j<cols; j++){w[p+j]=weightsF[l][i][j];}
			}
		}
	}

	/** Clears the gradient accumulators in place, rather than reallocating per batch. */
	private void zeroGradients(){
		for(int l=0; l<layers; l++){
			Arrays.fill(gW[l], 0);
			Arrays.fill(gB[l], 0);
		}
	}

	/**
	 * Forward and backward pass for one sample, accumulating into gW and gB.
	 * @param sample Index of the sample
	 * @return Squared error for this sample
	 */
	private double accumulateGradient(final int sample){
		final float[] vector=trainingData.inputs[sample];
		final float[] target=trainingData.targets[sample];
		final double[] a0=act[0];
		for(int i=0; i<numInputs; i++){a0[i]=(vector[i]-mean[i])/sd[i];}

		for(int l=0; l<layers; l++){
			final int rows=dims[l+1], cols=dims[l];
			final double[] w=weights[l], b=bias[l], prev=act[l], next=act[l+1];
			for(int i=0, p=0; i<rows; i++, p+=cols){
				double z=b[i];
				for(int j=0; j<cols; j++){z+=w[p+j]*prev[j];}
				if(hiddenFunc!=null){preAct[l+1][i]=z;}
				next[i]=(l==layers-1) ? fin(z, finalType)
					: (hiddenFunc==null ? Math.tanh(z) : hiddenFunc[l+1][i].activate(z));
			}
		}

		final double[] outs=act[layers];
		double sqErr=0, absErr=0;
		for(int k=0; k<numOutputs; k++){
			final double outv=outs[k];
			final double err=outv-target[k];
			sqErr+=err*err;
			absErr+=Math.abs(err);
			delta[layers][k]=2*err*finDeriv(outv, finalType);
		}
		if(sampleError!=null){sampleError[sample]=(float)(absErr/numOutputs);}

		for(int l=layers-1; l>=0; l--){
			final int rows=dims[l+1], cols=dims[l];
			final double[] w=weights[l], gw=gW[l], gb=gB[l];
			final double[] prev=act[l], dNext=delta[l+1];
			for(int i=0, p=0; i<rows; i++, p+=cols){
				final double d=dNext[i];
				gb[i]+=d;
				for(int j=0; j<cols; j++){gw[p+j]+=d*prev[j];}
			}
			if(l>0){
				final double[] dCur=delta[l];
				for(int j=0; j<cols; j++){
					double sum=0;
					for(int i=0; i<rows; i++){sum+=dNext[i]*w[i*cols+j];}
					dCur[j]=sum*(hiddenFunc==null ? (1-prev[j]*prev[j])
						: hiddenFunc[l][j].derivativeXFX(preAct[l][j], prev[j]));
				}
			}
		}
		return sqErr;
	}

	/**
	 * Applies one Adam update using the accumulated gradients.
	 * @param bs Batch size used for gradient averaging
	 * @param step Global step count, for bias correction
	 * @param lrNow Current learning rate
	 */
	private void applyAdam(final int bs, final long step, final double lrNow){
		final double c1=1-Math.pow(BETA1, step), c2=1-Math.pow(BETA2, step);
		for(int l=0; l<layers; l++){
			final int rows=dims[l+1], cols=dims[l];
			final double[] w=weights[l], b=bias[l];
			final double[] gw=gW[l], gb=gB[l];
			final double[] mw=mW[l], vw=vW[l], mb=mB[l], vb=vB[l];
			final boolean[] mask=(edgeMask==null) ? null : edgeMask[l];
			for(int i=0, p=0; i<rows; i++, p+=cols){
				for(int j=0; j<cols; j++){
					final int q=p+j;
					if(mask!=null && !mask[q]){continue;}//pruned edge stays exactly 0
					final double g=gw[q]/bs+wd*w[q];
					mw[q]=BETA1*mw[q]+(1-BETA1)*g;
					vw[q]=BETA2*vw[q]+(1-BETA2)*g*g;
					w[q]-=lrNow*(mw[q]/c1)/(Math.sqrt(vw[q]/c2)+EPS);
				}
				final double g=gb[i]/bs;
				mb[i]=BETA1*mb[i]+(1-BETA1)*g;
				vb[i]=BETA2*vb[i]+(1-BETA2)*g*g;
				b[i]-=lrNow*(mb[i]/c1)/(Math.sqrt(vb[i]/c2)+EPS);
			}
		}
	}

	/** @return Mean squared error over the validation split */
	private double validate(){
		if(numValid<1){return Double.MAX_VALUE;}
		double sum=0;
		if(validationData==null){
			for(int s=numTrain; s<numSamples; s++){
				final int sample=perm[s];
				final double[] p=predict(trainingData.inputs[sample], weights, bias);
				final float[] target=trainingData.targets[sample];
				for(int k=0; k<numOutputs; k++){
					final double e=p[k]-target[k];
					sum+=e*e;
				}
			}
		}else{
			for(int s=0; s<validationData.size(); s++){
				final double[] p=predict(validationData.inputs[s], weights, bias);
				final float[] target=validationData.targets[s];
				for(int k=0; k<numOutputs; k++){
					final double e=p[k]-target[k];
					sum+=e*e;
				}
			}
		}
		return sum/numValid;
	}

	/**
	 * Forward pass for one sample using reusable scratch buffers.
	 * @param vector Raw input vector
	 * @param w Weight arrays
	 * @param b Bias arrays
	 * @return Network output
	 */
	private double[] predict(final float[] vector, final double[][] w, final double[][] b){
		double[] cur=scratchA, next=scratchB;
		for(int i=0; i<numInputs; i++){cur[i]=(vector[i]-mean[i])/sd[i];}
		for(int l=0; l<layers; l++){
			final int rows=dims[l+1], cols=dims[l];
			final double[] wl=w[l], bl=b[l];
			for(int i=0, p=0; i<rows; i++, p+=cols){
				double z=bl[i];
				for(int j=0; j<cols; j++){z+=wl[p+j]*cur[j];}
				next[i]=(l==layers-1) ? fin(z, finalType)
					: (hiddenFunc==null ? Math.tanh(z) : hiddenFunc[l+1][i].activate(z));
			}
			final double[] t=cur; cur=next; next=t;
		}
		return cur;//caller reads the first numOutputs entries
	}

	/**
	 * Builds a CellNet from the trained weights and writes it in BBNet format,
	 * folding input standardization into the first layer so the saved net takes
	 * raw inputs.
	 */
	private void exportNet(){
		//net and its topology were built in buildNetAndMask(); only the values change here,
		//so a pruned edge trained as 0 is written into a slot the topology already agrees is 0.
		for(int l=0; l<layers; l++){
			final int rows=dims[l+1], cols=dims[l];
			final Cell[] layer=net.net[l+1];
			final double[] wl=weights[l], bl=bias[l];
			for(int i=0, p=0; i<rows; i++, p+=cols){
				final Cell c=layer[i];
				c.function=Function.getFunction((l==layers-1) ? finalFunction()
					: (hiddenType==null ? Function.TANH : hiddenType[l+1][i]));
				double b=bl[i];
				for(int j=0; j<cols; j++){
					double w=wl[p+j];
					if(l==0){//fold (x-mean)/sd into weights and bias
						w/=sd[j];
						b-=wl[p+j]*mean[j]/sd[j];
					}
					c.weights[j]=(float)w;
				}
				c.setBias((float)b, true);
			}
		}

		//Record the training command in the header as a #CL line, matching train.sh nets so the
		//saved .bbnet is self-documenting; CellNetParser reads #CL back into commands, and a netin=
		//continuation appends its command below the source net's, preserving the full history.
		//Guarded to run once: exportNet is now also called mid-training as a checkpoint (once per
		//improving epoch), and would otherwise append a duplicate #CL line on every checkpoint.
		if(!cmdLineRecorded){
			net.commands.add("#CL "+Shared.fullCommandline(false, true));
			cmdLineRecorded=true;
		}

		final ByteBuilder bb=net.toBytes();
		final ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();
		bsw.print(bb);
		errorState|=bsw.poisonAndWait();
		outstream.println("wrote "+out+" (valMSE="+String.format("%.6f", bestValid)+")");
	}

	/** @return The Function constant matching the selected output activation */
	private int finalFunction(){
		if(finalType==FINAL_LINEAR){return Function.LINEAR;}
		if(finalType==FINAL_RSLOG){return Function.RSLOG;}
		return Function.SIG;
	}

	/**
	 * Reloads the written net and compares it against the in-memory model, which
	 * catches export and serialization errors that training metrics cannot.
	 */
	private void roundTripCheck(){
		final CellNet check=CellNetParser.load(out);
		final int n=Math.min(1000, numSamples);
		double maxDiff=0;
		for(int s=0; s<n; s++){
			final float[] vector=trainingData.inputs[s];
			check.applyInput(vector);
			check.feedForward();
			final double[] b=predict(vector, weights, bias);
			for(int k=0; k<numOutputs; k++){
				maxDiff=Math.max(maxDiff, Math.abs(check.getOutput(k)-b[k]));
			}
		}
		final boolean ok=(maxDiff<1e-3);
		outstream.println(String.format(
			"round-trip check vs in-memory model: maxDiff=%.2e over %d samples %s",
			maxDiff, n, (ok ? "(OK)" : "(SUSPICIOUS!)")));
		errorState|=!ok;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Output-layer activation.
	 * @param z Pre-activation value
	 * @param finalType FINAL_LINEAR, FINAL_RSLOG, or FINAL_SIGMOID
	 * @return Activated value
	 */
	static double fin(double z, int finalType){
		if(finalType==FINAL_LINEAR){return z;}
		if(finalType==FINAL_RSLOG){return (z<0) ? -Math.log(1-z) : Math.log(1+z);}
		return 1.0/(1.0+Math.exp(-z));
	}

	/**
	 * Derivative of the output activation, expressed in terms of the ACTIVATED value.
	 * @param outv Activated output value
	 * @param finalType FINAL_LINEAR, FINAL_RSLOG, or FINAL_SIGMOID
	 * @return d(activation)/dz at that point
	 */
	static double finDeriv(double outv, int finalType){
		if(finalType==FINAL_LINEAR){return 1;}
		if(finalType==FINAL_RSLOG){return Math.exp(-Math.abs(outv));}
		return outv*(1-outv);
	}

	/** Deep copy of a ragged double array.
	 * @param x Source array
	 * @return An independent copy */
	private static double[][] deepCopy(double[][] x){
		final double[][] copy=new double[x.length][];
		for(int i=0; i<x.length; i++){copy[i]=x[i].clone();}
		return copy;
	}

	/** Row-oriented regression inputs and targets from one source file. */
	private static final class RegressionData {

		RegressionData(ArrayList<float[]> inputList, ArrayList<float[]> targetList){
			if(inputList.size()!=targetList.size()){
				throw new IllegalArgumentException("Input/target row count mismatch: "
					+inputList.size()+" != "+targetList.size());
			}
			inputs=inputList.toArray(new float[inputList.size()][]);
			targets=targetList.toArray(new float[targetList.size()][]);
		}

		int size(){return inputs.length;}

		final float[][] inputs;
		final float[][] targets;
	}

	/** Reusable task with private forward/backward scratch and gradient buffers. */
	private final class GradWorker implements Callable<Void> {

		GradWorker(float[][] act_, float[][] preAct_, float[][] delta_,
				float[][][] gw_, float[][] gb_){
			actLocal=act_; preActLocal=preAct_; deltaLocal=delta_; gw=gw_; gb=gb_;
		}

		void prepare(int[] order_, int start_, int end_){
			order=order_; start=start_; end=end_; sqErr=0;
			zeroGradientsF(gw, gb);
		}

		@Override
		public Void call() throws Exception{
			double sum=0;
			for(int s=start; s<end; s++){
				if((s&63)==0 && Thread.currentThread().isInterrupted()){
					throw new InterruptedException("Gradient worker interrupted.");
				}
				sum+=accumulateGradientF(order[s], actLocal, preActLocal, deltaLocal, gw, gb);
			}
			sqErr=sum;
			return null;
		}

		private final float[][] actLocal, preActLocal, deltaLocal;
		private final float[][][] gw;
		private final float[][] gb;
		private int[] order;
		private int start, end;
		private double sqErr;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String in=null;
	/** Optional external validation vectors; never included in normalization or training */
	private String valIn=null;
	private String out=null;
	/** Existing .bbnet to continue training from, instead of random initialization */
	private String netIn=null;
	/** Comma-separated hidden activations; one name is uniform, several are drawn per cell */
	private String hiddenSpec="tanh";
	/** Per-cell hidden activation, indexed [activationLayer][cell]; null means uniform tanh */
	private Function[][] hiddenFunc;
	/** Type code per cell, parallel to hiddenFunc, for export */
	private int[][] hiddenType;
	private int[] dims=null;

	private int epochs=60;
	private int batch=8192;
	private double lr=0.003;
	private double wd=1e-4;
	private long seed=1;
	private double vfraction=0.1;
	/** Whether vfraction was supplied explicitly, for rejecting ambiguous valin combinations */
	private boolean vfractionSet=false;
	private int finalType=FINAL_RSLOG;
	private boolean overwrite=true;
	/** Fraction of hidden-layer edges to keep; 1 is fully connected.  The output layer is always dense. */
	private float density=1f;
	/** Density override for the first hidden layer; 0 uses density */
	private float density1=0f;
	/** Block size for structured sparsity, passed through to CellNet */
	private int edgeBlockSize=1;
	/**
	 * Scan input lines for exponent notation, which Parse's fast path silently misreads.
	 * OFF by default: BBTools never emits scientific notation, and the scan costs ~8% of a
	 * load-dominated run, so it is not worth charging every run for a case that cannot occur.
	 * Turn it on for vector files produced by other tools, or use fjpd to force exact parsing.
	 */
	private boolean checkExponents=false;
	/**
	 * Train in float using the vectorized kernels in the simd package.
	 * OFF by default: SIMD reductions sum in a different order than scalar code, so results
	 * are close but not bit-identical (see simd/Vector.java:165), and this trainer's whole
	 * verification story rests on reproducing the scalar path exactly.
	 */
	private boolean useSimd=false;
	/** Maximum SIMD gradient workers; one preserves the historical accumulation path */
	private int threads=1;
	/** Round hidden-layer widths up to a whole number of SIMD lanes */
	private boolean padLayers=false;
	/**
	 * Use gathered sparse kernels instead of multiplying through the pruned zeros.
	 * OFF by default because it usually loses: the gather destroys the contiguous access
	 * SIMD depends on, and the overhead is comparable to the arithmetic it skips. Measured
	 * 26% SLOWER at density=0.75 on a narrow net, and only 9% faster at density=0.25 on a
	 * wide one despite skipping 71% of the edges. Worth trying only at aggressive pruning
	 * with wide layers.
	 */
	private boolean useSparse=false;
	/** Apply z-score weight standardization during training (CellNet's NORMALIZE_WEIGHTS is also off by default) */
	private boolean normalize=false;
	/** Blend strength for normalization; CellNet's default */
	private float normFactor=0.125f;
	/** Per-application decay of the blend strength; CellNet's default */
	private float normShrink=0.999f;
	/** Prioritize hard/stale samples instead of visiting everything uniformly */
	private boolean sortSamples=false;
	/** Training samples visited per epoch when sorting; 0 or more than available means all */
	private int setSize=0;

	/*--------------------------------------------------------------*/

	/** Training rows; also contains the internal validation rows when valin is absent */
	private RegressionData trainingData;
	/** Separate held-out rows when valin is supplied */
	private RegressionData validationData;
	private int numSamples;
	private int numInputs;
	/** Output width; 1 is the classic scalar case */
	private int numOutputs;
	private int numTrain;
	private int numValid;
	/** Shuffled training indices; the suffix is internal validation only when valin is absent */
	private int[] perm;

	private double[] mean;
	private double[] sd;

	private int layers;
	/** Per-layer weights, flat row-major, indexed [i*cols+j] */
	private double[][] weights;
	private double[][] bias;
	private double[][] mW, vW, mB, vB;
	private double[][] gW, gB;
	private double[][] bestWeights, bestBias;
	private double bestValid=Double.MAX_VALUE;
	/** Set once exportNet's #CL line has been recorded, so mid-training checkpoint
	 * writes (one per improving epoch) don't duplicate it on every call. */
	private boolean cmdLineRecorded=false;

	/** The net that will be exported; its topology also defines the edge mask */
	private CellNet net;
	/** Per-layer edge mask, flat row-major; null when fully connected */
	private boolean[][] edgeMask;

	private double[][] act;
	/** Pre-activation values; kept only when a hidden function needs them for its derivative */
	private double[][] preAct;
	private double[][] delta;
	private double[] scratchA, scratchB;

	private Random randy;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private final FileFormat ffin;
	private final FileFormat ffval;
	private final FileFormat ffout;

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	static final int FINAL_LINEAR=0, FINAL_RSLOG=1, FINAL_SIGMOID=2;

	/** Float mirrors used only by the vectorized path; each neuron's weights are contiguous */
	private float[][][] weightsF, gwF, mWF, vWF;
	/** Per-layer transposed weights, so the backward pass reads contiguous rows */
	private float[][][] wT;
	/** Sparse views used when a density mask is present: values plus the index sets they gather through */
	private float[][][] spW, spT;
	private int[][][] spWIdx, spTIdx;
	private float[][] biasF, gbF, mBF, vBF;
	private float[][] actF, deltaF;
	private float[][] preActF;
	private int gradientThreadCount=1;
	private ArrayList<GradWorker> gradientWorkers;
	private ExecutorService gradientPool;

	/** Most recent absolute error per sample; null unless sorting */
	private float[] sampleError;
	/** Epoch each sample was last visited; null unless sorting */
	private int[] sampleEpoch;

	/** Float lanes in a 256-bit vector; the pad8 rounding unit */
	private static final int LANE=8;

	private static final double BETA1=0.9, BETA2=0.999, EPS=1e-8;
	/** Epoch decay in the sort priority, matching ml.Sample.EPOCH_MULT */
	private static final float EPOCH_MULT=1/256f;

	private PrintStream outstream=System.err;
	public boolean errorState=false;

}
