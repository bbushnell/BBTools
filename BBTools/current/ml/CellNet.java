package ml;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Random;

import shared.Shared;
import shared.Tools;
import shared.Vector;
import structures.ByteBuilder;
import structures.FloatList;
import structures.IntList;

/**
 * Multi-layer feedforward neural network.
 * Supports dense and sparse connectivity, backpropagation training,
 * and configurable activation functions.
 * 
 * @author Brian Bushnell
 * @date October 25, 2013
 * @documentation Eru
 */
public class CellNet implements Cloneable, Comparable<CellNet> {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates an empty dummy neural network with null fields.
	 * Used as placeholder during construction or for testing purposes.
	 * All parameters are set to safe default values indicating an uninitialized state.
	 */
	public CellNet(){
		seed=0;
		dims=null;
		density=1;
		density1=0;
		edgeBlockSize=1;
		layers=0;
		values=null;
		weightsIn=weightsOut=null;
		eOverNet=null;
		list=null;
		net=null;
		finalLayer=null;
		transposeCounter=null;
		layerStart=null;
	}
	
	/**
	 * Creates multi-layer feedforward neural network with specified topology.
	 * Initializes network structure, allocates memory for activations and gradients,
	 * and sets up layer connectivity patterns based on density parameters.
	 * 
	 * Network uses dims[0] input neurons, dims[dims.length-1] output neurons,
	 * and dims.length-2 hidden layers. Connection density controls sparsity:
	 * density=1.0 creates fully connected layers, lower values create sparse connections.
	 * Edge block size enables structured sparsity for better memory locality.
	 * 
	 * @param dims_ Neurons per layer array (must have at least 2 layers)
	 * @param seed_ Random seed for weight initialization (negative uses system time)
	 * @param density_ Connection density for hidden layers (0.0-1.0, typically 0.1-1.0)
	 * @param density1_ First hidden layer density override (0 uses density_ value)
	 * @param edgeBlockSize_ Connectivity block size for structured sparsity (power of 2)
	 * @param commands_ Command line arguments used to create network (for reproducibility)
	 */
	public CellNet(int[] dims_, long seed_, float density_, float density1_, 
			int edgeBlockSize_, ArrayList<String> commands_) {
		seed=(seed_>=0 ? seed_ : new Random().nextLong()&Long.MAX_VALUE);
		commands=commands_; //Store creation parameters for reproducibility
		dims=dims_.clone(); //Deep copy to prevent external modification
		density=density_;
		density1=density1_;
		edgeBlockSize=edgeBlockSize_;
		layers=dims.length;
		values=makeFloatMatrix(dims); //Activation values per layer
		eOverNet=makeFloatMatrix(dims); //Error gradients per layer
		
		int cells=(int) shared.Vector.sum(dims); //Total neurons across all layers
		list=new ArrayList<Cell>(cells+1); //Global cell list (index 0 unused)
		transposeCounter=new int[cells+1]; //Sparse matrix transpose tracking
		net=makeNodes(dims, list, values, eOverNet); //Create network topology
		finalLayer=net[layers-1]; //Cache output layer for efficiency
		layerStart=new int[dims.length]; //Global cell ID offset per layer
		layerStart[0]=1; //Skip index 0 (reserved)
		for(int i=1; i<layerStart.length; i++){
			layerStart[i]=layerStart[i-1]+dims[i-1]; //Cumulative cell count
		}

//		weights=makeWeightsMatrix();
//		weightsOut=makeWeightsOutMatrix();
	}
	
	void makeWeightMatrices(){
//		assert(check());
		assert(weightsIn==null); //Ensure not already initialized
		weightsIn=makeWeightsInMatrix(); //Input weight matrix for forward pass
		weightsOut=makeWeightsOutMatrix(); //Output weight matrix for backprop
		edgesIn=makeEdgesInMatrix(); //Input connectivity for sparse networks
		edgesOut=makeEdgesOutMatrix(); //Output connectivity for transpose
//		assert(check());
	}
	
	private float[][][] makeWeightsInMatrix(){
		float[][][] x=new float[dims.length][][]; //[layer][neuron][weight]
		for(int layer=1; layer<x.length; layer++){ //Skip input layer
			Cell[] currentLayer=net[layer];
			float[][] y=x[layer]=new float[currentLayer.length][];
			for(int i=0; i<y.length; i++){
				y[i]=currentLayer[i].weights; //Direct reference for efficiency
			}
		}
		return x;
	}
	
	/**
	 * Creates input edge connectivity matrix for sparse networks.
	 * Structure: edgesIn[layer][neuron][input_index]
	 * Maps which neurons in previous layer connect to each neuron.
	 * 
	 * @return 3D array with direct references to cell input arrays
	 */
	private int[][][] makeEdgesInMatrix(){
		int[][][] x=new int[dims.length][][]; //[layer][neuron][input_id]
		for(int layer=1; layer<x.length; layer++){ //Skip input layer
			Cell[] currentLayer=net[layer];
			int[][] y=x[layer]=new int[currentLayer.length][];
			for(int i=0; i<y.length; i++){
				y[i]=currentLayer[i].inputs; //Sparse connectivity indices
			}
		}
		return x;
	}
	
	private float[][][] makeWeightsOutMatrix(){
		float[][][] x=new float[dims.length][][]; //[layer][neuron][output_weight]
		for(int layer=0; layer<x.length-1; layer++){ //All except final layer
			Cell[] current=net[layer];
			Cell[] next=net[layer+1];
			float[][] y=x[layer]=new float[current.length][];
			if(DENSE){
				for(int i=0; i<y.length; i++){
					y[i]=new float[next.length]; //Full connectivity
				}
			}else{
				for(int i=0; i<y.length; i++){
					assert(current[i].outputs!=null) : layer+", "+i+", "+current[i].id();
					y[i]=new float[current[i].outputs.length]; //Sparse connectivity
				}
			}
		}
		return x;
	}
	
	/**
	 * Creates output edge connectivity matrix for sparse networks.
	 * Structure: edgesOut[layer][neuron][output_index]
	 * Maps which neurons in next layer each neuron connects to.
	 * 
	 * @return 3D array with direct references to cell output arrays
	 */
	private int[][][] makeEdgesOutMatrix(){
		int[][][] x=new int[dims.length][][]; //[layer][neuron][output_id]
		for(int layer=0; layer<x.length-1; layer++){ //All except final layer
			Cell[] currentLayer=net[layer];
			int[][] y=x[layer]=new int[currentLayer.length][];
			for(int i=0; i<y.length; i++){
				y[i]=currentLayer[i].outputs; //Output connectivity indices
			}
		}
		return x;
	}
	
	public void randomize() {
		{
			Random randy=Shared.threadLocalRandom(seed); //Thread-safe RNG
			if(DENSE){
				makeEdgesDense(net, randy, density, density1, edgeBlockSize); //Full topology
//				check();
			}else{
				makeEdgesSparse(net, randy, density, density1, edgeBlockSize); //Sparse topology
//				check();
			}
			makeWeightMatrices(); //Optimize for training
//			check();
		}
		assert((Function.TYPE_RATES_CUM!=null)==(shared.Vector.sum(Function.TYPE_RATES)>0)) : Arrays.toString(Function.TYPE_RATES);
		if(Function.TYPE_RATES_CUM!=null || Cell.randomTypeRate>0){
			Random randy=Shared.threadLocalRandom(seed); //Same RNG for consistency
			if(Function.TYPE_RATES_CUM!=null){
				randomizeActivationA(net, list, randy, Function.TYPE_RATES_CUM); //Weighted selection
			}else{
				randomizeActivationB(net, list, randy, Cell.randomTypeRate); //Uniform random
			}
		}
	}
	
	/**
	 * Randomizes activation functions using cumulative distribution weights.
	 * Applies weighted random selection to hidden layer neurons only.
	 * Input and output layers retain their default activation functions.
	 * 
	 * @param net Network layer structure
	 * @param list Global cell list (unused but kept for signature compatibility)
	 * @param randy Random number generator (seeded for reproducibility)
	 * @param cumRate Cumulative probability distribution for function types
	 */
	static void randomizeActivationA(Cell[][] net, ArrayList<Cell> list, Random randy, float[] cumRate){
		assert(cumRate!=null && cumRate[cumRate.length-1]==1) : Arrays.toString(cumRate);
		
		for(int layerNum=1; layerNum<net.length-1; layerNum++){ //Only middle layers
			Cell[] layer=net[layerNum];
			for(Cell c : layer){
//				c.type=randomType(randy, cumRate);
				c.function=Function.randomFunction(randy); //Weighted random selection
			}
		}
	}
	
	/**
	 * Randomizes activation functions using uniform probability distribution.
	 * Each hidden layer neuron has 'rate' probability of changing its function.
	 * Ensures the new function differs from the current one to guarantee change.
	 * 
	 * @param net Network layer structure
	 * @param list Global cell list (unused but kept for signature compatibility)
	 * @param randy Random number generator (seeded for reproducibility)
	 * @param rate Probability (0.0-1.0) that each neuron changes function
	 */
	static void randomizeActivationB(Cell[][] net, ArrayList<Cell> list, Random randy, float rate){
		assert(rate>0); //Must have positive change probability
//		final int types=Cell.TYPES.length;
		final int types=Tools.min(Cell.MAX_TYPE+1, Function.TYPES.length); //Available function types
//		assert(false) : types;
//		System.err.println("types="+types);
		for(int layerNum=1; layerNum<net.length-1; layerNum++){ //Only middle layers
			Cell[] layer=net[layerNum];
			for(Cell c : layer){
				if(randy.nextFloat()<rate){ //Probabilistic change
//					int t=c.type;
//					while(t==c.type) {t=randy.nextInt(types);}
//					c.type=t;
					Function f=c.function;
					while(f==c.function){f=Function.getFunction(randy.nextInt(types));} //Ensure change
					c.function=f;
				}
			}
		}
	}
	
	/**
	 * Allocates 2D matrix for layer-wise neural network data storage.
	 * Creates one array per layer with size matching neuron count.
	 * Used for both activation values and error gradients.
	 * 
	 * @param dims Array containing neuron count per layer
	 * @return 2D matrix [layer][neuron] initialized to zero
	 */
	private static float[][] makeFloatMatrix(int[] dims){
		//TODO: Note - these can be made 1 longer with a constant of 1 to include bias
		float[][] matrix=new float[dims.length][]; //One array per layer
		for(int i=0; i<dims.length; i++){
			matrix[i]=new float[dims[i]]; //Size matches neuron count
		}
		return matrix;
	}
	
	private static Cell[][] makeNodes(int[] dims, ArrayList<Cell> list, float[][] values, float[][] eOverNext){
		final Cell[][] net=new Cell[dims.length][];
		assert(list.isEmpty()); //Must start with empty list
		list.add(null); //Reserve index 0
//		int prevWidth=-1, width=dims[0], nextWidth;
		int[] layerStart=new int[dims.length]; //Global ID offset per layer
		for(int i=0, sum=1; i<dims.length; i++){
			layerStart[i]=sum; //Cumulative neuron count
			sum+=dims[i];
		}
		for(int layerNum=0; layerNum<dims.length; layerNum++){
			final int prevLayerStart=(layerNum==0 ? -1 : layerStart[layerNum-1]);
			final int nextLayerStart=(layerNum>=dims.length-1 ? -1 : layerStart[layerNum+1]);
			final int width=dims[layerNum];
			final int type=(layerNum<dims.length-1) ? Cell.defaultActivationType : Cell.finalLayerType; //Hidden vs output
			Cell[] layer=net[layerNum]=new Cell[width];
			float[] lvals=values[layerNum]; //Shared activation array
			float[] eons=eOverNext[layerNum]; //Shared gradient array
			assert(lvals.length==width) : layerNum+", "+lvals.length+", "+width;
			for(int i=0; i<width; i++){
				Cell c=new Cell(list.size(), type, i, layerNum, dims.length-1, 
						prevLayerStart, nextLayerStart, width, lvals, eons); //Complete initialization
				layer[i]=c;
				list.add(c); //Add to global list
			}
		}
		return net;
	}
	
	/**
	 * Calculates effective connection density for a specific layer.
	 * First hidden layer can use special density1 parameter for different sparsity.
	 * Output layer is always fully connected (density=1.0).
	 * 
	 * @param density Default connection density for hidden layers
	 * @param density1 Special density for first hidden layer (0 uses default)
	 * @param layer Layer number being processed (1-based for hidden layers)
	 * @param layers Total number of layers in network
	 * @return Effective density (0.0-1.0) for this layer
	 */
	private static final float layerDensity(float density, float density1, int layer, int layers){
		return layer>=layers+1 ? 1f : (layer==1 && density1>0) ? density1 : density; //Output always dense
	}
	
	/**
	 * Creates dense network connectivity with probabilistic edge pruning.
	 * Allocates full weight matrices but randomly zeros edges based on density.
	 * Output layer remains fully connected regardless of density parameters.
	 * Uses structured sparsity via edgeBlockSize for memory locality.
	 * 
	 * @param net Network layer structure to initialize
	 * @param randy Random number generator for edge selection
	 * @param density Base connection probability for hidden layers
	 * @param density1 Special density for first hidden layer
	 * @param edgeBlockSize Block size for structured connectivity patterns
	 * @return Total number of non-zero edges created
	 */
	private static long makeEdgesDense(Cell[][] net, Random randy, 
			float density, float density1, int edgeBlockSize){
		long numEdges=0;
		for(int layerNum=1; layerNum<net.length; layerNum++){
			final float layerDensity=layerDensity(density, density1, layerNum, net.length);
			Cell[] layer=net[layerNum], prev=net[layerNum-1];
			final boolean finalLayer=(layerNum==net.length-1);
//			final int minLen=Tools.mid(1, 5, prev.length/3);
			for(Cell c : layer){
				assert(c.weights==null); //Ensure uninitialized
				c.weights=new float[prev.length]; //Full connectivity matrix
				c.deltas=new float[prev.length]; //Gradient accumulation
				c.setBias(randomBias(randy, 0.8f), true); //TODO: 0.8f should be a parameter
				final float[] weights=c.weights;
//				int nonzero=0;
//				for(int attempt=0; attempt<10 && nonzero<minLen; attempt++) {
//					nonzero=0;
//					for(int i=0; i<weights.length; i++) {
//						float w=randomWeight(randy, 0.5f);
//						w=(finalLayer || randy.nextFloat()<=density ? w : 0);
//						weights[i]=w;
//						numEdges+=(w==0 ? 0 : 1);
//						nonzero+=(w==0 ? 0 : 1);
//					}
//				}
				if(finalLayer){ //Output layer always fully connected
					for(int i=0; i<weights.length; i++){
						weights[i]=randomWeight(randy);
						numEdges++;
					}
				}else{ //Hidden layer uses density pruning
					BitSet bs=pickEdges(weights.length, edgeBlockSize, layerDensity, randy); //Structured sparsity
//					System.err.println(bs);
					for(int i=bs.nextSetBit(0); i>=0; i=bs.nextSetBit(i+1)){
						weights[i]=randomWeight(randy); //Initialize selected edges
						numEdges++;
					}
				}
//				assert(false) : Arrays.toString(c.weights);
			}
		}
		return numEdges;
	}
	
	/**
	 * Creates sparse network connectivity with explicit edge lists.
	 * Stores only non-zero connections to minimize memory usage.
	 * Output layer remains fully connected, hidden layers use density pruning.
	 * Must call makeOutputSets() afterward to build reverse connectivity.
	 * 
	 * @param net Network layer structure to initialize
	 * @param randy Random number generator for edge selection
	 * @param density Base connection probability for hidden layers
	 * @param density1 Special density for first hidden layer
	 * @param edgeBlockSize Block size for structured connectivity patterns
	 * @return Total number of edges created
	 */
	private static long makeEdgesSparse(Cell[][] net, Random randy, 
			float density, float density1, int edgeBlockSize){
		long numEdges=0;
		for(int layerNum=1; layerNum<net.length; layerNum++){
			final float layerDensity=layerDensity(density, density1, layerNum, net.length);
			Cell[] layer=net[layerNum], prev=net[layerNum-1];
			final boolean finalLayer=(layerNum==net.length-1);
//			final int minLen=Tools.mid(1, 5, prev.length/3);
			for(Cell c : layer){
				assert(c.weights==null); //Ensure uninitialized
				c.setBias(randomBias(randy, 0.8f), true); //TODO: 0.8f should be a parameter
				if(finalLayer){ //Output layer fully connected
					c.inputs=new int[prev.length]; //All previous neurons
					for(int i=0; i<c.inputs.length; i++){c.inputs[i]=i;} //Sequential indices
				}else{ //Hidden layer uses sparse connectivity
					BitSet bs=pickEdges(prev.length, edgeBlockSize, layerDensity, randy); //Select edges
//					System.err.println(bs);
					c.inputs=toArray(bs); //Convert to index array
				}
				c.weights=new float[c.inputs.length]; //Only store non-zero weights
				c.deltas=new float[c.inputs.length]; //Matching gradient array
				for(int i=0; i<c.weights.length; i++){
					c.weights[i]=randomWeight(randy); //Initialize all selected edges
				}
				numEdges+=(c.weights.length); //Count actual connections
			}
		}
		makeOutputSets(net); //Build reverse connectivity for backprop
		return numEdges;
	}
	
	/**
	 * Selects network edges using structured sparsity patterns.
	 * First determines target edge count based on density, then uses
	 * block-aligned selection for better memory locality. Ensures minimum
	 * connectivity to prevent disconnected neurons.
	 * 
	 * @param width Number of potential connections (previous layer size)
	 * @param edgeBlockSize Block alignment for structured sparsity
	 * @param density Connection probability (0.0-1.0)
	 * @param randy Random number generator for edge selection
	 * @return BitSet marking selected edge positions
	 */
	private static BitSet pickEdges(int width, int edgeBlockSize, float density, Random randy){
		int toMake=0;
		int min=Tools.mid(1, 5, width/3); //Minimum connectivity guarantee
		for(int i=0; i<width; i++){
			if(randy.nextFloat()<=density){toMake++;} //Probabilistic selection
		}
		toMake=Tools.max(toMake, min); //Enforce minimum
		int mod=toMake%edgeBlockSize;
		if(mod!=0){toMake=Tools.min(width, toMake-mod+edgeBlockSize);} //Block alignment
		BitSet bs=new BitSet(width);
		int range=(width-1)/edgeBlockSize; //Block range
		for(int made=0; made<toMake;){
			final int rand=randy.nextInt(range+1);
			assert(rand>=0 && rand<=range) : rand+", "+range+", "+edgeBlockSize+", "+made+", "+width;
			//TODO: BUG (Confirmed by Brian) - Causes intermittent crashes after ~12 hours of training
			//PROBLEM: Line 360 calculates 'rand' but never uses it
			//Line 363 calls randy.nextInt(range+1) AGAIN instead of using 'rand'
			//When second call returns exactly 'range', then range*edgeBlockSize can exceed array bounds
			//FIX: Replace 'randy.nextInt(range+1)' with 'rand' on line 363
			//IMPACT: 4 interns experiencing occasional training crashes
			int start=randy.nextInt(range+1)*edgeBlockSize; //Block start position
			assert(start>=0) : rand+", "+range+", "+edgeBlockSize+", "+made+", "+width+", "+start;
			if(start<width && !bs.get(start)){ //Valid unused block
				for(int i=0; i<edgeBlockSize && i+start<width; i++){
					int loc=i+start;
					assert(loc>=0) : rand+", "+range+", "+edgeBlockSize+", "+made+", "+width+", "+start+", "+loc+", "+i;
					bs.set(loc); //Mark edge as selected
					made++; //Count created edges
				}
			}
		}
//		System.err.println(bs);
		return bs;
	}
	
	static int[] toArray(BitSet bs){
		int[] array=new int[bs.cardinality()]; //Size equals number of set bits
		for(int i=bs.nextSetBit(0), j=0; i>=0; i=bs.nextSetBit(i+1), j++){
			array[j]=i; //Store bit position
		}
		return array;
	}
	
	static void makeOutputSets(Cell[][] net){
		final int cells; //Total neuron count
		{
			Cell[] lastLayer=net[net.length-1];
			Cell c=lastLayer[lastLayer.length-1];
			cells=c.id(); //Highest global ID
		}
		IntList[] map=new IntList[cells+1]; //Temporary output lists
		for(int lnum=0; lnum<net.length-1; lnum++){ //All except final layer
			for(Cell c : net[lnum]){
				assert((c.inputs!=null || c.layer==0) && c.outputs==null); //Valid state
				map[c.id()]=new IntList(); //Initialize output list
			}
		}
		for(int lnum=1; lnum<net.length; lnum++){ //All except input layer
			for(Cell c : net[lnum]){
				for(int i : c.inputs){ //Each input connection
					map[i+c.prevLayerStart].add(c.lpos); //Map reverse connection
				}
			}
		}
		for(int lnum=0; lnum<net.length-1; lnum++){ //All except final layer
			for(Cell c : net[lnum]){
				c.outputs=map[c.id()].toArray(); //Convert to final array
//				System.err.println("c"+c.id()+" outputs="+c.outputs.length);
			}
		}
	}
	
	private static float randomBias(Random randy, float probFlat){
		float weight=0;
		while(Math.abs(weight)<=Float.MIN_NORMAL){ //Ensure non-zero
			if(randy.nextFloat()<=probFlat){
				weight=(1-2*randy.nextFloat())*.25f; //Uniform distribution
			}else{
				weight=randy.nextFloat()*randy.nextFloat()*randy.nextFloat()*2*(randy.nextBoolean() ? 1 : -1); //Cubic distribution
			}
		}
		return weight;
	}
	
	private static float randomWeight(Random randy){
		float weight=0;
		while(Math.abs(weight)<=Float.MIN_NORMAL){ //Ensure non-zero
			final float f=randy.nextFloat();
			if(f<=PROB_FLAT){
				weight=(1-2*randy.nextFloat())*.25f; //Uniform component
			}else if(f<=PROB_FLAT+PROB_EXP){
				weight=(float)(Tools.exponential(randy, EXP_LAMDA)*(randy.nextBoolean() ? 1f : -1f)); //Exponential component
			}else{
				weight=randy.nextFloat()*randy.nextFloat()*randy.nextFloat()*2*(randy.nextBoolean() ? 1f : -1f); //Cubic component
			}
		}
		weight=Tools.mid(-RAND_WEIGHT_CAP, weight, RAND_WEIGHT_CAP); //Clamp extreme values
		return weight;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	
	/**
	 * Applies simulated annealing to all network weights.
	 * Adds controlled random noise to weights to escape local minima
	 * and improve generalization. Strength controls noise magnitude.
	 * 
	 * @param strength Annealing intensity (higher = more noise)
	 * @param randy Random number generator for noise generation
	 */
	void anneal(float strength, Random randy){
		//assert(check());
		for(Cell c : list){
			if(c!=null){c.anneal(strength, randy);} //Apply to each neuron
		}
		
		//assert(check());
	}
	
	/**
	 * Processes single training sample through complete forward/backward pass.
	 * Applies input, computes forward propagation, calculates error, and
	 * optionally performs backpropagation to update weights. Handles both
	 * dense and sparse network architectures automatically.
	 * 
	 * @param s Training sample with input, target, and metadata
	 * @param backProp Whether to perform backpropagation weight updates
	 * @param weightMult Learning rate multiplier for this sample
	 */
	public void processSample(Sample s, boolean backProp, float weightMult){
		//assert(check());
		weightMult*=s.weight; //Apply sample-specific weighting
		
		applyInput(s.in); //Load input values
		if(DENSE){
			feedForwardDense(); //Dense forward pass
		}else{
			feedForwardSparse(); //Sparse forward pass
		}

//		s.errorValue=0;
		for(int i=0; i<finalLayer.length; i++){
			s.result[i]=finalLayer[i].value(); //Extract output values
//			s.errorValue+=(s.result[i]-s.goal[i]);
		}
		
		s.calcError(weightMult); //Calculate prediction error
		
		if(!backProp){return;} //Skip training if inference only
		
//		double error=calcError(s.goal);
//		assert(error==s.errorMagnitude) : "\n"+s.result[0]+", "+s.goal[0]+", "+s.errorMagnitude+"\n"
//				+ finalLayer[0].value()+", "+s.goal[0]+", "+error+"\n";
		
//		s.errorMagnitude=(float)error;
//		s.weightedErrorMagnitude=(float)error*Cell.getWeightMultiplier(s.result[0], s.goal[0]);
		
//		errorSum+=s.errorMagnitude;
//		weightedErrorSum+=s.weightedErrorMagnitude;
		if(DENSE){
			backPropDense(s.goal, weightMult); //Dense backpropagation
		}else{
			backPropSparse(s.goal, weightMult); //Sparse backpropagation
		}
		//assert(check());
	}
	
	/**
	 * Apply input values to the network's input layer.
	 * Copies values into the first layer's value array for processing.
	 * 
	 * @param valuesIn Input values as FloatList (size must match input layer)
	 * @return This network for chaining
	 */
	public CellNet applyInput(FloatList valuesIn) {
		//assert(check());
		assert(valuesIn.size==dims[0]) : valuesIn.size+"; "+Arrays.toString(dims);
		Vector.copy(values[0], valuesIn.array); // Copy to input layer
		return this;
	}
	
	/**
	 * Apply input values to the network's input layer.
	 * 
	 * @param valuesIn Input values as float array (length must match input layer)
	 * @return This network for chaining
	 */
	public CellNet applyInput(float[] valuesIn) {
		//assert(check());
		assert(valuesIn.length==dims[0]) : valuesIn.length+"; "+Arrays.toString(dims);
		Vector.copy(values[0], valuesIn); // Copy to input layer
		return this;
	}
	
	/**
	 * Perform forward propagation through the network.
	 * Automatically selects dense or sparse implementation based on network type.
	 * 
	 * @return Output value from the final layer
	 */
	public float feedForward() {
		return (CellNet.DENSE ? feedForwardDense() : feedForwardSparse()); // Delegate to appropriate implementation
	}
	
	/**
	 * Forward propagation for sparse networks.
	 * Processes only active connections to improve performance.
	 * 
	 * @return Output value from the final layer
	 */
	public float feedForwardSparse(){
		//assert(check());
		for(int lnum=1; lnum<layers; lnum++){ // Process each layer
			final float[] valuesIn=values[lnum-1]; // Previous layer outputs
			Cell[] layer=net[lnum]; // Current layer cells
			for(int cnum=0; cnum<layer.length; cnum++) {
				Cell c=layer[cnum];
				//				//assert(c.check());
				c.summateSparse(valuesIn, edgeBlockSize); // Sparse summation
				//				//assert(c.check());
			}
		}
		return finalLayer[0].value(); // Return final output
	}
	
	public float feedForwardDense(){
		//assert(check());
		for(int lnum=1; lnum<layers; lnum++){ // Process each layer
			final float[] valuesIn=values[lnum-1]; // Previous layer outputs
			Cell[] layer=net[lnum]; // Current layer cells
			if(SPECIAL_FMA) { // Use fused multiply-add optimization
				Vector.feedForwardDense(layer, valuesIn);
			}else {
				for(int cnum=0; cnum<layer.length; cnum++) {
					Cell c=layer[cnum];
					//				//assert(c.check());
					c.summateDense(valuesIn); // Standard dense summation
					//				//assert(c.check());
				}
			}
		}
		return finalLayer[0].value(); // Return final output
	}
	
	/**
	 * Performs backpropagation for sparse networks using gradient descent.
	 * Computes gradients from output layer backwards through hidden layers,
	 * updating weight deltas for later batch application. Uses sparse connectivity
	 * to skip non-existent edges for computational efficiency.
	 * 
	 * @param truth Target output values for error computation
	 * @param weightMult Learning rate multiplier for gradient scaling
	 */
	public void backPropSparse(float[] truth, float weightMult){
		//assert(check());
		{ //Output layer processing
			final float[] valuesIn=values[values.length-2]; //Previous layer activations
			for(int i=0; i<finalLayer.length; i++){
				final Cell c=finalLayer[i];
				//Final layer is always dense
				c.updateEdgesFinalLayerDense(truth[i], valuesIn, weightMult); //Output layer gradients
			}
		}
		
		for(int lnum=layers-2; lnum>0; lnum--){ //Hidden layers (backwards)
			final float[] valuesIn=values[lnum-1]; //Input activations
			final float[] eOverNetNext=eOverNet[lnum+1]; //Next layer gradients
			Cell[] layer=net[lnum];
			float[][] weightsOutLnum=weightsOut[lnum]; //Forward weights for gradient propagation
			
			for(int i=0; i<layer.length; i++){
				Cell c=layer[i];
				//assert(c.check());
				c.updateEdgesHiddenLayerSparse(valuesIn, eOverNetNext, weightsOutLnum[i], edgeBlockSize); //Hidden layer gradients
				//assert(c.check());
			}
		}
		//assert(check());
	}
	
	/**
	 * Performs backpropagation for dense networks with optional optimizations.
	 * Uses vectorized operations when SPECIAL_FMA enabled for faster computation.
	 * Includes optional weight/bias normalization to prevent gradient explosion
	 * and improve training stability through standardization.
	 * 
	 * @param truth Target output values for error computation
	 * @param weightMult Learning rate multiplier for gradient scaling
	 */
	public void backPropDense(float[] truth, float weightMult){
		//assert(check());
		{ //Output layer processing
			final float[] valuesIn=values[values.length-2]; //Previous layer activations
			for(int i=0; i<finalLayer.length; i++){
				final Cell c=finalLayer[i];
				//assert(c.check());
				c.updateEdgesFinalLayerDense(truth[i], valuesIn, weightMult); //Output layer gradients
				//assert(c.check());
			}
		}
		
		for(int lnum=layers-2; lnum>0; lnum--){ //Hidden layers (backwards)
			final float[] valuesIn=values[lnum-1]; //Input activations
			final float[] eOverNetNext=eOverNet[lnum+1]; //Next layer gradients
			Cell[] layer=net[lnum];
			float[][] weightsOutLnum=weightsOut[lnum]; //Forward weights for gradient propagation
			
			if(SPECIAL_FMA){Vector.backPropFma(layer, eOverNetNext, weightsOutLnum);} //Vectorized optimization
			for(int i=0; i<layer.length; i++){
				Cell c=layer[i];
				//assert(c.check());
				c.updateEdgesHiddenLayerDense(valuesIn, eOverNetNext, weightsOutLnum[i]); //Hidden layer gradients
				//assert(c.check());
			}
		}
		//assert(check());
		if(normalization_factor>0.0001f){ //Apply normalization if enabled
			if(NORMALIZE_BIAS){normalizeBiasDense();} //Standardize bias values
			if(NORMALIZE_WEIGHTS){normalizeWeightsDense();} //Standardize weight values
			normalization_factor*=normalization_shrink_rate; //Decay normalization strength
		}
	}
	
	/**
	 * Normalizes bias values across hidden layers using z-score standardization.
	 * Computes mean and standard deviation for each layer, then applies
	 * weighted normalization to prevent extreme bias values that can
	 * destabilize training. Uses gradual adjustment to maintain stability.
	 */
	public void normalizeBiasDense(){
		FloatList list=new FloatList(); //Temporary storage for statistics
		for(int lnum=layers-2; lnum>0; lnum--){ //Hidden layers only
			list.clear();
			Cell[] layer=net[lnum];
			for(int i=0; i<layer.length; i++){
				Cell c=layer[i];
				list.add(c.bias); //Collect all bias values
			}
			float mean=(float)list.mean(); //Layer bias mean
			float stdev=list.stdev(); //Layer bias standard deviation
			if(stdev<0.001f){stdev=0.001f;} //Avoid division by zero
			float invStdev=1f/stdev; //Normalization factor
			for(int i=0; i<layer.length; i++){
				Cell c=layer[i];
				float f=c.bias;
				float target=(f-mean)*invStdev; //Z-score normalization
				f=(1f-normalization_factor)*f+normalization_factor*target; //Weighted update
				c.bias=f; //Apply normalized bias
			}
		}
	}
	
	/**
	 * Normalizes all network weights using layer-wise standardization.
	 * Applies z-score normalization to prevent gradient explosion and
	 * improve training convergence. Updates weight matrices and rebuilds
	 * transpose matrices for backpropagation consistency.
	 */
	public void normalizeWeightsDense(){
		FloatList list=new FloatList(); //Reusable statistics collector
		for(int layer=0; layer<weightsIn.length-1; layer++){
			float[][] weightsInL=weightsIn[layer+1]; //Current layer weights
			normalizeDense(weightsInL, list); //Apply normalization
		}
		transposeDense(); //Rebuild transpose matrices
	}
	
	/**
	 * Normalizes weight matrix using z-score standardization on non-zero values.
	 * Computes statistics over all non-zero weights to preserve sparsity,
	 * then applies weighted normalization to prevent extreme values.
	 * Ensures minimum weight magnitude to avoid underflow issues.
	 * 
	 * @param weights 2D weight matrix to normalize
	 * @param list Reusable FloatList for statistics computation
	 */
	public void normalizeDense(float[][] weights, FloatList list){
		if(weights==null){return;} //Skip null matrices
		list.clear();
		for(float[] array : weights){
			for(float f : array){
				if(f!=0){ //Only collect non-zero weights
					list.add(f);
				}
			}
		}
		float mean=(float)list.mean(); //Weight mean
		float stdev=list.stdev(); //Weight standard deviation
		if(stdev<0.001f){stdev=0.001f;} //Avoid division by zero
		float invStdev=1f/stdev; //Normalization factor
		for(float[] array : weights){
			for(int i=0; i<array.length; i++){
				float f=array[i];
				if(f!=0){ //Only normalize non-zero weights
					float target=(f-mean)*invStdev; //Z-score normalization
					f=(1f-normalization_factor)*f+normalization_factor*target; //Weighted update
					array[i]=(Math.abs(f)<1e-12f) ? 1e-12f : f; //Prevent underflow
				}
			}
		}
	}
	
	/**
	 * Updates transpose weight matrices for backpropagation.
	 * Synchronizes weightsOut matrices with current weightsIn values
	 * to ensure gradient computation uses latest weights. Required
	 * after weight updates or normalization.
	 */
	void transpose(){
		if(DENSE){transposeDense();} //Dense matrix transpose
		else{transposeSparse();} //Sparse matrix transpose
	}
	
	/**
	 * Transposes dense weight matrices for backward pass computation.
	 * Copies weights from input-oriented (weightsIn) to output-oriented
	 * (weightsOut) format needed for gradient propagation. Uses direct
	 * matrix indexing for full connectivity.
	 */
	void transposeDense(){
		for(int layer=0; layer<weightsIn.length-1; layer++){
			float[][] weightsInL=weightsIn[layer+1], weightsOutL=weightsOut[layer]; //Source and destination matrices
			for(int lnum=0; lnum<weightsOutL.length; lnum++){
				float[] weightsOutC=weightsOutL[lnum]; //Output neuron's weights
				for(int lnumOut=0; lnumOut<weightsInL.length; lnumOut++){
					weightsOutC[lnumOut]=weightsInL[lnumOut][lnum]; //Transpose weight value
				}
			}
		}
	}
	
	/**
	 * Transposes sparse weight matrices using connectivity mapping.
	 * Uses edge connectivity information to map weights from input-oriented
	 * to output-oriented format. Requires careful indexing due to sparse
	 * structure. Uses transposeCounter to track mapping progress.
	 */
	void transposeSparse(){
		Arrays.fill(transposeCounter, 0); //Reset mapping counters
		for(int layer=0; layer<weightsIn.length-1; layer++){
			float[][] weightsInL=weightsIn[layer+1], weightsOutL=weightsOut[layer]; //Source and destination
			int[][] edgesOutL=edgesOut[layer]; //Output connectivity map
			for(int lnum=0; lnum<weightsOutL.length; lnum++){
				float[] weightsOutC=weightsOutL[lnum]; //Output neuron's weights
				int[] outputs=edgesOutL[lnum]; //Connected output neurons
				assert(outputs.length==weightsOutC.length); //Consistency check
				for(int i=0; i<outputs.length; i++){
					int lnumOut=outputs[i]; //Output neuron local position
					int cidOut=lnumOut+layerStart[layer+1]; //Output neuron global ID
//					assert(false) : Arrays.toString(layerStart)+", "+layer+", "+lnumOut+", "+cidOut;
					int inputNum=transposeCounter[cidOut]; //Input index for this connection
					transposeCounter[cidOut]++; //Advance counter
					weightsOutC[i]=weightsInL[lnumOut][inputNum]; //Copy weight value
				}
			}
		}
		for(Cell c : list){ //Verify transpose correctness
			assert(c==null || c.inputs==null || transposeCounter[c.id()]==c.inputs.length);
		}
	}
	
	/**
	 * Apply accumulated weight changes after a training batch.
	 * Updates all network weights using the specified learning rate and sample count.
	 * 
	 * @param samples Number of samples processed in this batch
	 * @param alpha Learning rate for weight updates
	 */
	public void applyChanges(int samples, float alpha) {
		//assert(check());
		float invSamples=1f/samples; // Average gradients over batch
		assert(invSamples<=1) : samples+", "+invSamples;
		for(Cell c : list) {
			if(c!=null) {
				c.applyUpdates(invSamples, alpha); // Update weights and biases
			}
		}
		//assert(check());
		epochsTrained++; // Track training progress
		samplesTrained+=samples;
	}
	
	/**
	 * Get output value from a specific neuron in the final layer.
	 * 
	 * @param outnum Index of output neuron (0-based)
	 * @return Activation value of the specified output neuron
	 */
	public float getOutput(int outnum){
		return net[layers-1][outnum].value(); // Access final layer directly
	}
	
	/**
	 * Get all output values from the final layer.
	 * 
	 * @return Array containing activation values of all output neurons
	 */
	public float[] getOutput(){
		Cell[] outLayer=net[layers-1]; // Final layer reference
		float[] out=new float[outLayer.length];
		for(int i=0; i<outLayer.length; i++) {
			out[i]=outLayer[i].value(); // Copy each output value
		}
		return out;
	}
	
	/**
	 * Calculates total network error using legacy method.
	 * Computes sum of squared errors between predictions and targets.
	 * Deprecated in favor of Sample.calcError() for better error weighting.
	 * 
	 * @param truth Target output values
	 * @return Total squared error across all outputs
	 * @deprecated Use Sample.calcError() for weighted error computation
	 */
	@Deprecated
	public double calcError(float[] truth){
		double error=0;
		for(int i=0; i<finalLayer.length; i++){
			Cell c=finalLayer[i];
			float t=truth[i]; //Target value
			float e=c.calcError(t); //Cell-specific error
			error+=e; //Accumulate total error
		}
		return error;
	}

	
	/*--------------------------------------------------------------*/
	/*----------------            X            ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Generates detailed string representation of network structure.
	 * Shows layer-by-layer breakdown with neuron details for debugging.
	 * Includes activation values, weights, and connectivity information.
	 * 
	 * @return Multi-line string describing entire network state
	 */
	public String toString(){
		StringBuilder sb=new StringBuilder();
		for(int layernum=0; layernum<layers; layernum++){
			Cell[] layer=net[layernum];
			sb.append("\n* Layer "+layernum+", nodes="+layer.length+" *"); //Layer header
			for(Cell c : layer){
				sb.append("\n"+c.toString()); //Individual neuron details
			}
		}
		return sb.toString();
	}
	
	/**
	 * Generates network metadata header for serialization.
	 * Creates structured header with network parameters, training statistics,
	 * and performance metrics. Used for network file format output.
	 * 
	 * @return ByteBuilder containing formatted network header
	 */
	public ByteBuilder header(){
		ByteBuilder bb=new ByteBuilder();
		bb.append("##bbnet").nl(); //File format identifier
		bb.append("#version ").append(version).nl();
		if(CONCISE){bb.append("#concise").nl();} //Format options
		bb.append(OUT_DENSE ? "#dense" : OUT_SPARSE ? "#sparse" : DENSE ? "#dense" : "#sparse").nl();
		bb.append("#density ").append(density, 8, true).nl(); //Network parameters
		if(density1>0 && density1!=density){bb.append("#density1 ").append(density1, 8, true).nl();}
		bb.append("#blocksize ").append(edgeBlockSize).nl();
		bb.append("#seed ").append(seed).nl();
//		if(annealSeed!=seed) {bb.append("#annealseed ").append(annealSeed).nl();}
		bb.append("#layers ").append(layers).nl();
		if(epochsTrained>0){bb.append("#epochs ").append(epochsTrained).nl();} //Training history
		if(samplesTrained>0){bb.append("#samples ").append(samplesTrained).nl();}
		
		bb.append("#dims"); //Network architecture
		for(int d : dims){bb.space().append(d);}
		bb.nl();
		for(String s : commands){bb.append(s).nl();} //Creation commands
		if(lastStats!=null){bb.append("##stats ").append(lastStats).nl();}
		bb.append("##fpr ").append(fpRate, 6).nl(); //Performance metrics
		bb.append("##fnr ").append(fnRate, 6).nl();
		bb.append("##err ").append(errorRate, 6).nl();
		bb.append("##wer ").append(weightedErrorRate, 6).nl();
		bb.append("##ctf ").append(cutoff, 6).nl();
		return bb;
	}
	
	/**
	 * Serializes complete network to binary format with weights and structure.
	 * Outputs network in BBTools format with header, layer definitions,
	 * connectivity patterns, and weight values. Supports both dense and
	 * sparse representations with optional hexadecimal encoding.
	 * 
	 * @return ByteBuilder containing complete serialized network
	 */
	public ByteBuilder toBytes(){
		
		ByteBuilder bb=header(); //Include metadata header
		bb.append("#edges ").append(countEdges()).nl(); //Total edge count
		
		lastLinesWritten=18; //Track output size
		long edgeCount=0;
		assert(CONCISE); //Require compact format
		for(int lnum=1; lnum<layers; lnum++){ //Process each layer
			Cell[] layer=net[lnum];
			Cell[] prev=net[lnum-1];
			bb.append("##layer ").append(lnum).nl(); //Layer header
			for(Cell c : layer){ //Each neuron in layer
				if(DENSE){ //Dense network format
					if(OUT_SPARSE){ //Force sparse output
						lastLinesWritten+=2;
						
						if(OUT_HEX){ //Hexadecimal connectivity
							bb.append('H').append(c.id()).space();
							BitSet bs=new BitSet(c.weights.length);
							for(int i=0; i<c.weights.length; i++){
								if(c.weights[i]!=0){bs.set(i);} //Mark non-zero weights
							}
							toHex(bs, bb); //Encode as hex
						}else{ //Text connectivity
							bb.append('I').append(c.id()); //Input line
							for(int i=0; i<c.weights.length; i++){
								if(c.weights[i]!=0){
									bb.space().append(i); //Non-zero indices
								}
							}
						}
						bb.nl();

						bb.append('W').append(c.id()).space().append(c.typeString()); //Weight line
						bb.space().append(c.bias(), 6, true); //Bias value
						for(int i=0; i<c.weights.length; i++){
							if(c.weights[i]!=0){
								bb.space().append(c.weights[i], 6, true); //Weight values
								edgeCount++;
							}
						}
						bb.nl();
					}else{ //Dense output format
						lastLinesWritten++;
						bb.append('C').append(c.id()).space().append(c.typeString()); //Complete line
						bb.space().append(c.bias(), 6, true); //Bias
						for(int i=0; i<c.weights.length; i++){bb.space().append(c.weights[i], 6, true);} //All weights
						bb.nl();
						edgeCount+=c.weights.length;
					}
				}else{ //Sparse network format
					if(OUT_DENSE){ //Force dense output
						lastLinesWritten++;
						bb.append('C').append(c.id()).space().append(c.typeString());
						bb.space().append(c.bias(), 6, true);
//						for(int idx=0, nextInput=0; idx<c.inputs.length; idx++) {
//							for(int inum=c.inputs[idx]; nextInput<inum; nextInput++) {
//								bb.space().append(0);
//							}
//							bb.space().append(c.weights[idx], 6, true);
//							nextInput++;
//						}
						int inum=0;
						for(int idx=0; idx<c.inputs.length; inum++){ //Expand sparse to dense
							if(inum==c.inputs[idx]){
								bb.space().append(c.weights[idx], 6, true); //Weight value
								idx++;
							}else{
								bb.space().append(0); //Zero for missing connection
							}
						}
						for(; inum<dims[c.layer-1]; inum++){bb.space().append(0);} //Fill remaining
						bb.nl();
						edgeCount+=c.weights.length;
					}else{ //Native sparse output
						lastLinesWritten+=2;
						if(OUT_HEX){ //Hex connectivity
							bb.append('H').append(c.id()).space();
							toHex(c.inputs, bb); //Encode input indices
						}else{ //Text connectivity
							bb.append('I').append(c.id()); //Input line
							for(int i=0; i<c.inputs.length; i++){bb.space().append(c.inputs[i]);} //Input indices
						}
						bb.nl();

						bb.append('W').append(c.id()).space().append(c.typeString()); //Weight line
						bb.space().append(c.bias(), 6, true); //Bias
						for(int i=0; i<c.weights.length; i++){bb.space().append(c.weights[i], 6, true);} //Weights
						bb.nl();
						edgeCount+=c.weights.length;
					}
				}
			}
		}
		return bb;
	}
	
	/**
	 * Converts integer array to hexadecimal BitSet representation.
	 * Creates BitSet from integer indices then encodes as hex string.
	 * Used for compact sparse connectivity serialization.
	 * 
	 * @param set Array of integer indices to encode
	 * @param bb ByteBuilder to append hex encoding to
	 * @return ByteBuilder with appended hex representation
	 */
	static ByteBuilder toHex(int[] set, ByteBuilder bb){
		if(set==null || set.length==0){return bb.append(0);} //Empty set
		final int max=set[set.length-1]; //Highest index
		BitSet bs=new BitSet(max);
		for(int e : set){bs.set(e);} //Set bits for each index
		return toHex(bs, bb); //Convert to hex
	}
	
	/**
	 * Converts BitSet to hexadecimal string representation.
	 * Encodes BitSet bytes as hex digits for compact storage.
	 * Note: Current implementation only uses digits 0-9, not full hex.
	 * 
	 * @param bs BitSet to convert
	 * @param bb ByteBuilder to append hex string to
	 * @return ByteBuilder with appended hex representation
	 */
	static ByteBuilder toHex(BitSet bs, ByteBuilder bb){
		byte[] bytes=bs.toByteArray(); //Convert BitSet to bytes
		for(byte b : bytes){
			//TODO: POSSIBLE BUG - Only handles digits 0-9, not hex A-F. Should use full hex encoding.
			bb.append((byte)('0'+(b&15))); //Low nibble
			bb.append((byte)('0'+((b>>4)&15))); //High nibble
		}
		return bb;
	}
	
	/**
	 * Converts hexadecimal string back to integer array.
	 * Decodes hex-encoded BitSet back to original integer indices.
	 * Used for deserializing sparse connectivity patterns.
	 * 
	 * @param line Byte array containing hex-encoded data
	 * @param start Starting position in byte array
	 * @return Integer array of decoded indices
	 */
	static int[] fromHex(byte[] line, int start){
		byte[] bytes=new byte[(line.length-start)/2]; //Decode pairs
		for(int i=start, j=0; i<line.length; i+=2, j++){
			byte a=line[i], b=line[i+1]; //Hex digit pair
			//TODO: POSSIBLE BUG - No validation that a,b are valid hex chars. Invalid chars will produce wrong results.
			int x=((a-'0')&15)|(((b-'0')&15)<<4); //Reconstruct byte
			bytes[j]=(byte)x;
		}
		BitSet bs=BitSet.valueOf(bytes); //Convert to BitSet
		return toArray(bs); //Extract indices
	}
	
	/**
	 * Creates deep copy of network with identical architecture and state.
	 * Copies weights, biases, training statistics, and configuration parameters.
	 * Thread-safe operation for concurrent model evaluation and training.
	 * 
	 * @param copyDelta Whether to copy accumulated gradients (for training)
	 * @return New CellNet instance with copied state
	 */
	public synchronized CellNet copy(boolean copyDelta){
//		assert(check());
		CellNet copy=new CellNet(dims, seed, density, density1, edgeBlockSize, commands); //Create structure
		ArrayList<Cell> list2=copy.list;
		for(int i=1; i<list.size(); i++){ //Copy all neurons
			Cell c=list.get(i);
//			assert(c.check());
			Cell c2=list2.get(i);
			c2.setFrom(c, copyDelta); //Copy cell state
//			System.err.println("Checking "+c.id()+": copyDelta="+copyDelta+", "+c.weights+", "+c.outputs
//					+", "+c2.weights+", "+c2.outputs);
//			assert(c2.check());
		}
		copy.makeWeightMatrices(); //Build weight matrices
//		assert(copy.check());
		copy.commands=commands; //Copy metadata
		copy.errorRate=errorRate; //Copy performance metrics
		copy.weightedErrorRate=weightedErrorRate;
		copy.fpRate=fpRate;
		copy.fnRate=fnRate;
		copy.tpRate=tpRate;
		copy.tnRate=tnRate;
		copy.cutoff=cutoff; //Copy configuration
		copy.alpha=alpha;
		copy.annealStrength=annealStrength;
//		copy.annealSeed=annealSeed;
		copy.epoch=epoch; //Copy training progress
		copy.epochsTrained=epochsTrained;
		copy.lastStats=lastStats;
		copy.samplesTrained=samplesTrained;
		
//		assert(copy.check());
		return copy;
	}
	
	/**
	 * Copies state from another network of identical architecture.
	 * Updates weights, biases, training statistics, and configuration
	 * parameters. Used for model checkpointing and population-based training.
	 * 
	 * @param cn Source network to copy from
	 * @param copyWeight2 Whether to copy accumulated gradients (deltas)
	 * @return This network after copying
	 */
	public CellNet setFrom(CellNet cn, boolean copyWeight2){
		//assert(cn.check());
		//assert(check());
		ArrayList<Cell> list2=cn.list;
		for(int i=1; i<list.size(); i++){ //Copy all neurons
			Cell c=list.get(i);
			//assert(c.check());
			Cell c2=list2.get(i);
			c.setFrom(c2, copyWeight2); //Copy cell state
		}
		//assert(cn.check());
		//assert(check());
		commands=cn.commands; //Copy creation parameters
		errorRate=cn.errorRate; //Copy performance metrics
		weightedErrorRate=cn.weightedErrorRate;
		fpRate=cn.fpRate;
		fnRate=cn.fnRate;
		tpRate=cn.tpRate;
		tnRate=cn.tnRate;
		cutoff=cn.cutoff; //Copy classification threshold
		alpha=cn.alpha; //Copy learning parameters
		annealStrength=cn.annealStrength;
//		annealSeed=cn.annealSeed;
		epoch=cn.epoch; //Copy training progress
		epochsTrained=cn.epochsTrained;
		lastStats=cn.lastStats;
		samplesTrained=cn.samplesTrained;
		return this;
	}
	
	/**
	 * Validates internal network consistency and state.
	 * Checks all neurons for valid weights, connectivity, and data integrity.
	 * Used for debugging and ensuring network correctness.
	 * 
	 * @return true if network passes all consistency checks
	 */
	public boolean check(){
		for(int i=1; i<list.size(); i++){ //Check all neurons
			Cell c=list.get(i);
			boolean b=c.check(); //Individual cell validation
			if(!b){return false;} //Fail on first error
		}
		return true; //All checks passed
	}
	
	/**
	 * Accumulates error gradients from another network.
	 * Adds gradient values from parallel training threads or batches.
	 * Used in distributed training and gradient accumulation.
	 * 
	 * @param net2 Source network containing gradients to add
	 */
	public void addError(CellNet net2){
//		errorSum+=net2.errorSum;
//		weightedErrorSum+=net2.weightedErrorSum;
		for(int i=1; i<list.size(); i++){ //Add all neuron errors
			Cell c=list.get(i);
			Cell c2=net2.list.get(i);
			c.addError(c2); //Accumulate gradient
		}
	}
	
	/**
	 * Thread-safe gradient accumulation from another network.
	 * Synchronizes access to prevent race conditions during parallel training.
	 * Used for multi-threaded batch processing and gradient aggregation.
	 * 
	 * @param net2 Source network containing gradients to accumulate
	 */
	public void accumulate(CellNet net2){
		synchronized(this){ //Lock this network
			synchronized(net2){ //Lock source network
				for(int i=1; i<list.size(); i++){ //Accumulate all neurons
					Cell c=list.get(i);
					Cell c2=net2.list.get(i);
//					synchronized(c) {
						c.accumulate(c2); //Thread-safe accumulation
//					}
				}
			}
		}
	}
	
	/**
	 * Thread-safe clearing of temporary training data.
	 * Resets accumulated gradients and temporary variables.
	 * Called between training batches to prepare for next iteration.
	 */
	synchronized public void clear(){
//		errorSum=0;
//		weightedErrorSum=0;
		for(Cell c : list){ //Clear all neurons
			if(c!=null){c.clearTemp();} //Reset temporary data
		}
	}
	
	/**
	 * Counts total non-zero connections in the network.
	 * Iterates through all weight matrices to find active edges.
	 * Works for both dense and sparse networks by checking weight values.
	 * Used for network analysis and memory usage estimation.
	 * 
	 * @return Number of non-zero weighted connections
	 */
	public long countEdges(){
		long sum=0;
		for(float[][] y : weightsIn){ //Each layer's weight matrix
			if(y!=null){
				for(float[] z : y){ //Each neuron's weights
					if(z!=null){
						for(float f : z){ //Each individual weight
							sum+=(f==0 ? 0 : 1); //Count non-zero weights
						}
					}
				}
			}
		}
		return sum;
	}
	
	/**
	 * Compares networks based on performance metrics for ranking.
	 * Uses configurable comparison criteria (error rate, FPR, FNR, etc.)
	 * with fallback hierarchy for model selection and population sorting.
	 * Lower error rates are considered better (return 1).
	 * 
	 * @param b Network to compare against
	 * @return 1 if this network is better, -1 if worse, 0 if equal
	 */
	@Override
	public int compareTo(CellNet b){
		if(b==null){
			return 1; //Non-null is better than null
		}else if(compareCode==compareWER && weightedErrorRate!=b.weightedErrorRate){
			return weightedErrorRate<b.weightedErrorRate ? 1 : -1; //Lower WER is better
		}else if(compareCode==compareERR && errorRate!=b.errorRate){
			return errorRate<b.errorRate ? 1 : -1; //Lower error is better
		}else if(compareCode==compareFNR && fnRate!=b.fnRate){
			return fnRate<b.fnRate ? 1 : -1; //Lower FNR is better
		}else if(compareCode==compareFPR && fpRate!=b.fpRate){
			return fpRate<b.fpRate ? 1 : -1; //Lower FPR is better
		}
		
		//Fallback hierarchy for tie-breaking
		if(weightedErrorRate!=b.weightedErrorRate){
			return weightedErrorRate<b.weightedErrorRate ? 1 : -1;
		}else if(errorRate!=b.errorRate){
			return errorRate<b.errorRate ? 1 : -1;
		}else if(fnRate!=b.fnRate){
			return fnRate<b.fnRate ? 1 : -1;
		}else if(fpRate!=b.fpRate){
			return fpRate<b.fpRate ? 1 : -1;
		}
		
		return 0; //Networks are equivalent
	}
	
	/**
	 * Gets the primary comparison metric value for this network.
	 * Returns the error metric currently used for network ranking.
	 * Used for quick access to sorting criterion without full comparison.
	 * 
	 * @return Current primary error metric value
	 */
	final float pivot(){
		return compareCode==compareWER ? weightedErrorRate : 
			compareCode==compareERR ? errorRate : 
			compareCode==compareFNR ? fnRate :
			compareCode==compareFPR ? fpRate : weightedErrorRate; //Default to WER
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Get the number of layers in this network.
	 * 
	 * @return Total number of layers including input and output
	 */
	public int numLayers(){
		return dims.length;
	}
	
	/**
	 * Get the number of input neurons.
	 * 
	 * @return Size of the input layer
	 */
	public int numInputs(){
		return dims[0];
	}
	
	/**
	 * Get the number of output neurons.
	 * 
	 * @return Size of the output layer
	 */
	public int numOutputs(){
		return dims[dims.length-1];
	}
	
	/**
	 * Sets classification threshold for binary prediction decisions.
	 * Values above cutoff are classified as positive, below as negative.
	 * Critical parameter for balancing precision/recall tradeoffs.
	 * 
	 * @param f New cutoff threshold value (typically 0.0-1.0)
	 */
	public void setCutoff(float f){
//		System.err.println("Set cutoff from "+cutoff+" to "+f);
		cutoff=f;
	}
	
	/*--------------------------------------------------------------*/
	
	/**
	 * Checks equality based on training epoch number.
	 * Networks are considered equal if trained for same number of epochs.
	 * Used for sorting and comparison during model selection.
	 * 
	 * @param b Object to compare with
	 * @return true if both networks have same epoch count
	 */
	@Override
	public boolean equals(Object b){
		return epoch==((CellNet)b).epoch;
	}

	/**
	 * Generates hash code based on training epoch.
	 * Consistent with equals() method for proper hash table behavior.
	 * 
	 * @return Hash code derived from epoch number
	 */
	@Override
	public int hashCode(){return epoch;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Current classification error rate (0.0-1.0) */
	float errorRate=999;
	/** Weighted classification error rate accounting for class imbalance */
	float weightedErrorRate=999;
	/** False positive rate - incorrectly classified as positive */
	float fpRate=999;
	/** False negative rate - incorrectly classified as negative */
	float fnRate=999;
	/** True positive rate - correctly classified as positive */
	float tpRate=-999;
	/** True negative rate - correctly classified as negative */
	float tnRate=-999;
	/** Learning rate for weight updates */
	float alpha=-1;
	/** Simulated annealing strength for noise injection */
	float annealStrength=-1;
	/** Classification threshold for binary decisions */
	public float cutoff=-1;
	/** Current training epoch number */
	int epoch=-1;
	/** Display counter for printing progress */
	int count=1; //for printing
	
	/** Total epochs completed during training */
	long epochsTrained=0;
	/** Total samples processed during training */
	long samplesTrained=0;
	/** Last recorded training statistics string */
	String lastStats=null;
	
	/*--------------------------------------------------------------*/
	
	/** Random seed for reproducible weight initialization */
	final long seed;
	/** Number of layers in network (including input and output) */
	final int layers;
	/** Neuron count per layer [input, hidden1, hidden2, ..., output] */
	final int[] dims; //Stores widths of layers
	/** Connection density for hidden layers (0.0-1.0) */
	final float density;
	/** Special connection density for first hidden layer */
	final float density1;
	/** Block size for structured sparsity patterns */
	final int edgeBlockSize;
	/** 2D network structure [layer][neuron] */
	final Cell[][] net;
	/** Direct reference to output layer for efficiency */
	final Cell[] finalLayer;
	/** Global list of all neurons (index 0 reserved) */
	final ArrayList<Cell> list;
	
	/** Activation values per layer [layer][neuron] */
	final float[][] values;
	/** Error gradients per layer [layer][neuron] */
	final float[][] eOverNet;
	/** Input weight matrices [layer][neuron][input_weight] */
	float[][][] weightsIn;
	/** Output weight matrices for backpropagation [layer][neuron][output_weight] */
	float[][][] weightsOut;
	/** Input connectivity indices for sparse networks [layer][neuron][input_id] */
	int[][][] edgesIn;
	/** Output connectivity indices for sparse networks [layer][neuron][output_id] */
	int[][][] edgesOut;
	/** Temporary counters for sparse matrix transpose operations */
	final int[] transposeCounter; //??
	/** Starting global ID for each layer in the cell list */
	final int[] layerStart; //??
	
	/** 
	 * Command line arguments used to create this network.
	 * Enables network reproduction with identical parameters.
	 * Contains original training configuration and hyperparameters.
	 */
	public ArrayList<String> commands;
	
	/** Number of lines written during last serialization operation */
	long lastLinesWritten=0;
	
	/*--------------------------------------------------------------*/
	
	/** Enable vectorized fused multiply-add operations for ~20% speedup */
	public static final boolean SPECIAL_FMA=true; //~20% faster when true, but needs ml.Cell class
	/** Whether to apply bias normalization during training */
	public static boolean NORMALIZE_BIAS=false;
	/** Whether to apply weight normalization during training */
	public static boolean NORMALIZE_WEIGHTS=false;
	/** Strength of normalization applied (0.0-1.0) */
	public static float normalization_factor=0.125f;
	/** Rate of normalization strength decay per epoch */
	public static float normalization_shrink_rate=0.999f;
	/** Use compact serialization format */
	public static boolean CONCISE=true;
	/** Network type: true=dense connectivity, false=sparse */
	public static boolean DENSE=true;
	/** Output format: use hexadecimal encoding */
	public static boolean OUT_HEX=false;
	/** Output format: force dense representation */
	public static boolean OUT_DENSE=false;
	/** Output format: force sparse representation */
	public static boolean OUT_SPARSE=false;
	/** Enable verbose debug output */
	public static boolean verbose=false;
	/** Network file format version number */
	public static final int version=1;

	/** Probability of uniform distribution in weight initialization */
	public static float PROB_FLAT=0.3f;
	/** Probability of exponential distribution in weight initialization */
	public static float PROB_EXP=0.4f;
	/** Lambda parameter for exponential distribution */
	public static float EXP_LAMDA=5f;
	/** Maximum absolute weight value to prevent training instability */
	public static float RAND_WEIGHT_CAP=2.0f;
	
	/** Current comparison metric for network sorting (0=WER, 1=ERR, 2=FNR, 3=FPR) */
	static int compareCode=0;
	/** Comparison codes for different error metrics */
	final static int compareWER=0, compareERR=1, compareFNR=2, compareFPR=3;
	
	/*--------------------------------------------------------------*/
	
}
