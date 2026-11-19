package test.ml;

import ml.BSig;
import ml.BTanh;
import ml.BungeeAffine;
import ml.BungeeFormula;
import ml.Cell;
import ml.CellNet;
import ml.Function;

/**
 * Independent, immutable scalar inference oracle for ML tests.
 * The source must be quiescent while it is copied. Indexed rows are expanded
 * into dense rows; no source cells, mutable arrays or serving policies survive.
 * Arithmetic follows scalar dense inference, with a rounded float dot product
 * followed by a separate float bias addition. Sparse/SIMD reduction orders may
 * differ. This helper never delegates evaluation back to CellNet.
 * Adapted from the raw path of prot.MagQCForwardTemplate.
 * @author Brian Bushnell, Nilou, Yoimiya, Raiden
 */
strictfp final class ScalarForwardSnapshot {

	/** Copies one frozen model; serving interpretation and output caps are out of scope. */
	static ScalarForwardSnapshot snapshot(final CellNet source){
		require(source!=null, "A quiescent source network is required");
		return new ScalarForwardSnapshot(source);
	}

	/** Owns copies of every numeric array and shares only audited immutable functions. */
	private ScalarForwardSnapshot(final CellNet source){
		require(source.dims!=null && source.dims.length>=2, "Inference requires input and output layers");
		dims=source.dims.clone();
		require(source.net!=null && source.net.length==dims.length, "Cell layers must match dimensions");
		for(int layer=0; layer<dims.length; layer++){
			require(dims[layer]>0 && source.net[layer]!=null && source.net[layer].length==dims[layer],
				"Invalid cell layer at "+layer);
		}
		clip=source.inputClip();
		require(Float.isFinite(clip) && clip>=0, "Input clip must be finite and nonnegative");
		mean=source.inputMeanCopy();
		inverseStd=source.inputInverseStdCopy();
		require((mean==null)==(inverseStd==null), "Input normalization arrays must be paired");
		if(clip>0){
			require(mean!=null && dims.length>=3, "Clipping requires normalization and a hidden layer");
			for(final Cell cell:source.net[1]){
				require(cell!=null && (cell.function==Function.getFunction(Function.TANH)
					|| cell.function==Function.getFunction(Function.SIG)
					|| cell.function instanceof BTanh || cell.function instanceof BSig
					|| cell.function instanceof BungeeAffine || cell.function instanceof BungeeFormula
					|| standardClippedFunction(cell.function)),
					"Snapshot clipping must match CellNet.validateInputClip's activation policy");
			}
		}
		if(mean!=null){
			require(mean.length==inputs() && inverseStd.length==inputs(), "Input normalization width differs");
			for(int i=0; i<mean.length; i++){
				require(Float.isFinite(mean[i]) && Float.isFinite(inverseStd[i]) && inverseStd[i]>0,
					"Invalid normalization at input "+i);
			}
		}
		weights=new float[dims.length-1][][];
		biases=new float[dims.length-1][];
		functions=new Function[dims.length-1][];
		long total=0;
		for(int layer=1; layer<dims.length; layer++){
			final int width=dims[layer-1], count=dims[layer], at=layer-1;
			weights[at]=new float[count][];
			biases[at]=new float[count];
			functions[at]=new Function[count];
			for(int cell=0; cell<count; cell++){
				final Cell original=source.net[layer][cell];
				require(original!=null && original.weights!=null, "Missing inference row at "+layer+":"+cell);
				functions[at][cell]=auditedFunction(original.function);
				require(Float.isFinite(original.bias), "Nonfinite bias at "+layer+":"+cell);
				biases[at][cell]=original.bias;
				final float[] row=new float[width];
				if(original.inputs==null){
					require(original.weights.length==width, "Dense row width differs at "+layer+":"+cell);
					System.arraycopy(original.weights, 0, row, 0, width);
				}else{
					require(original.inputs.length==original.weights.length, "Indexed row length mismatch");
					final boolean[] seen=new boolean[width];
					for(int edge=0; edge<original.weights.length; edge++){
						final int index=original.inputs[edge];
						require(index>=0 && index<width && !seen[index], "Invalid or duplicate source index "+index);
						seen[index]=true;
						row[index]=original.weights[edge];
					}
				}
				for(final float value:row){require(Float.isFinite(value), "Nonfinite inference weight");}
				weights[at][cell]=row;
				total+=width;
			}
		}
		assert(total>0) : "Positive validated layer widths must produce a nonempty dense snapshot";
	}

	/** Rejects custom mutable functions without invoking their methods. */
	private static Function auditedFunction(final Function function){
		if(function instanceof BSig || function instanceof BTanh
			|| function instanceof BungeeAffine || function instanceof BungeeFormula){return function;}
		for(int type=Function.SIG; type<=Function.SOFTPLUS; type++){
			if(function==Function.getFunction(type)){return function;}
		}
		throw new IllegalArgumentException("Unsupported activation instance; explicit audit is required");
	}

	/** Only audited immutable standard singletons qualify for clipped input. */
	private static boolean standardClippedFunction(final Function function){
		if(function==Function.getFunction(Function.SWISH)){return true;}
		for(int type=Function.RELU; type<=Function.SOFTPLUS; type++){
			if(function==Function.getFunction(type)){return true;}
		}
		return false;
	}

	/** Allocates reusable evaluation scratch in the thread that will use it. */
	Worker newWorker(){return new Worker(this);}
	/** Number of caller-owned input values. */
	private int inputs(){return dims[0];}
	/** Number of raw model outputs. */
	private int outputs(){return dims[dims.length-1];}

	/** One thread's scratch; shared template parameters remain immutable. */
	static strictfp final class Worker {
		private Worker(final ScalarForwardSnapshot template_){
			template=template_;
			owner=Thread.currentThread();
			values=new float[template.dims.length][];
			for(int i=0; i<values.length; i++){values[i]=new float[template.dims[i]];}
			assert(values[0].length==template.inputs()) : "Worker scratch must follow the frozen input width";
		}

		/** Copies raw outputs to distinct caller storage, without per-evaluation allocation. */
		void score(final float[] input, final float[] output){
			require(owner==Thread.currentThread(), "Inference worker belongs to another thread");
			require(input!=null && input.length==template.inputs()
				&& output!=null && output.length==template.outputs(),
				"Input/output widths differ from the frozen model");
			require(input!=output, "Caller input and output buffers must be distinct");
			System.arraycopy(input, 0, values[0], 0, input.length);
			if(template.mean!=null){
				for(int i=0; i<input.length; i++){
					values[0][i]-=template.mean[i];
					values[0][i]*=template.inverseStd[i];
				}
				if(template.clip>0){
					for(int i=0; i<input.length; i++){
						final float value=values[0][i];
						// CellNet.normalizeInput rejects overflow before applying the bound.
						if(!Float.isFinite(value)){
							throw new IllegalArgumentException("Nonfinite normalized input before clipping at "+i);
						}
						values[0][i]=value>template.clip ? template.clip : value< -template.clip ? -template.clip : value;
					}
				}
			}
			for(int layer=0; layer<template.weights.length; layer++){
				final float[] previous=values[layer], next=values[layer+1];
				final float[][] matrix=template.weights[layer];
				for(int row=0; row<next.length; row++){
					final float[] weights=matrix[row];
					float dot=0;
					for(int column=0; column<previous.length; column++){dot+=previous[column]*weights[column];}
					// Scalar dense inference rounds the dot before adding the bias.
					float sum=template.biases[layer][row];
					sum+=dot;
					next[row]=(float)template.functions[layer][row].activate(sum);
				}
			}
			final float[] raw=values[values.length-1];
			System.arraycopy(raw, 0, output, 0, raw.length);
		}

		private final ScalarForwardSnapshot template;
		private final Thread owner;
		private final float[][] values;
	}

	/** Rejects invalid snapshots and API use even when assertions are disabled. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private final int[] dims;
	private final float[][][] weights;
	private final float[][] biases;
	private final Function[][] functions;
	private final float[] mean, inverseStd;
	private final float clip;
}
