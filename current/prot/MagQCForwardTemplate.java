package prot;

import ml.Cell;
import ml.CellNet;
import ml.Function;
import ml.BSig;
import ml.BTanh;
import ml.BungeeAffine;
import ml.BungeeFormula;

/**
 * Experimental immutable snapshot of scalar dense inference parameters.
 * The caller must hold the source quiescent during construction. Afterwards no
 * source CellNet, Cell, weight/index array or mutable metadata is retained.
 * Indexed rows are expanded exactly as CellNet.copyDenseForInference expands them.
 * Arithmetic follows Vector.feedForwardDense with Shared.SIMD=false; this is
 * not a bit-parity promise for the original sparse or SIMD reduction order.
 * No production caller selects this prototype. Compilation, actual-bundle parity
 * and measured resource/performance gates are required before integration.
 * @author Brian Bushnell, Nilou, Yoimiya
 */
public strictfp final class MagQCForwardTemplate {

	/** Snapshots raw inference, or the six-head ordinary-subnet serving contract. */
	public static MagQCForwardTemplate snapshot(final CellNet source, final boolean servedSubnet){
		return snapshot(source, servedSubnet, false);
	}

	/** Explicit bundle policy also applies to experimental replay; rawOutput remains raw. */
	public static MagQCForwardTemplate snapshot(final CellNet source, final boolean servedSubnet, final boolean capped){
		require(source!=null, "A quiescent source network is required");
		final MagQCPredictionContract interpretation=servedSubnet ? MagQCPredictionContract.parseSubnet(
			source.getTag("magqc_prediction_contract"), source.getTag("output_contract")) : null;
		require(!capped || servedSubnet || source.numOutputs()==10, "Size caps require ten model outputs");
		final MagQCSubnetOutputCaps caps=!capped ? null : servedSubnet ?
			MagQCSubnetOutputCaps.ordinary(interpretation) : MagQCSubnetOutputCaps.size();
		return new MagQCForwardTemplate(source, interpretation, caps);
	}

	/** Copies all mutable numeric state before publishing this immutable object. */
	private MagQCForwardTemplate(final CellNet source, final MagQCPredictionContract interpretation,
			final MagQCSubnetOutputCaps caps_){
		require(source.dims!=null && source.dims.length>=2, "Inference requires input and output layers");
		dims=source.dims.clone();
		require(source.net!=null && source.net.length==dims.length, "Cell layers must match dimensions");
		for(int layer=0; layer<dims.length; layer++){
			require(dims[layer]>0 && source.net[layer]!=null && source.net[layer].length==dims[layer],
				"Invalid cell layer at "+layer);
		}
		prediction=interpretation;
		caps=caps_;
		clip=source.inputClip();
		require(Float.isFinite(clip) && clip>=0, "Input clip must be finite and nonnegative");
		require(prediction==null || (outputs()==6 && prediction.rawOutputs()==outputs()),
			"Ordinary serving requires exactly six model outputs");
		mean=source.inputMeanCopy(); inverseStd=source.inputInverseStdCopy();
		require((mean==null)==(inverseStd==null), "Input normalization arrays must be paired");
		if(clip>0){
			require(mean!=null && dims.length>=3, "Clipping requires explicit normalization and a hidden layer");
			for(final Cell cell:source.net[1]){
				require(cell!=null && (cell.function==Function.getFunction(Function.TANH) || cell.function==Function.getFunction(Function.SIG)
					|| cell.function instanceof BTanh || cell.function instanceof BSig || cell.function instanceof BungeeAffine
					|| cell.function instanceof BungeeFormula || standardClippedFunction(cell.function)),
					"Snapshot clipping must match CellNet.validateInputClip's explicit activation policy");
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
			biases[at]=new float[count]; functions[at]=new Function[count];
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
						seen[index]=true; row[index]=original.weights[edge];
					}
				}
				for(final float value:row){require(Float.isFinite(value), "Nonfinite inference weight");}
				weights[at][cell]=row; total+=width;
			}
		}
		weightCount=total;
		assert(weightCount>0) : "Positive validated layer widths must produce a nonempty dense snapshot";
	}

	/** Audited singleton functions and final immutable bungee instances may be shared.
	 * Their activate methods use Functions/Math and final parameters; mutable
	 * Function.TYPE_RATES affects selection, not activation. A custom type/name is
	 * insufficient: reject it without calling its methods or retaining its state. */
	private static Function auditedFunction(final Function function){
		// These final classes own immutable parameters, so two models cannot interfere.
		if(function instanceof BSig || function instanceof BTanh || function instanceof BungeeAffine || function instanceof BungeeFormula){return function;}
		for(int type=Function.SIG; type<=Function.SOFTPLUS; type++){
			if(function==Function.getFunction(type)){return function;}
		}
		throw new IllegalArgumentException("Unsupported activation instance; explicit audit is required");
	}

	/** Only the audited immutable standard singletons qualify; never call custom functions. */
	private static boolean standardClippedFunction(final Function function){
		if(function==Function.getFunction(Function.SWISH)){return true;}
		for(int type=Function.RELU; type<=Function.SOFTPLUS; type++){
			if(function==Function.getFunction(type)){return true;}
		}
		return false;
	}

	/** Creates scratch in the executing thread, sharing only this template's immutable parameters. */
	public Worker newWorker(){return new Worker(this);}
	/** Number of caller-owned input values. */
	public int inputs(){return dims[0];}
	/** Number of model outputs, excluding private serving diagnostics. */
	public int outputs(){return dims[dims.length-1];}
	/** Defensive shape copy; no internal array is exposed. */
	public int[] dimensions(){return dims.clone();}
	/** Logical dense weight count, not a measured JVM heap or residency estimate. */
	public long weights(){return weightCount;}

	/** One thread's reusable scratch. No cells, weights, training or serialization API escapes. */
	public static strictfp final class Worker {
		/** Allocates only activation and optional report scratch, never parameter copies. */
		private Worker(final MagQCForwardTemplate template_){
			template=template_; owner=Thread.currentThread();
			values=new float[template.dims.length][];
			for(int i=0; i<values.length; i++){values[i]=new float[template.dims[i]];}
			reported=template.prediction!=null ? new float[template.prediction.outputNames().length] :
				template.caps==null ? null : new float[template.outputs()];
			assert(values[0].length==template.inputs()) : "Worker scratch must follow the frozen input width";
		}

		/** Copies served outputs to caller storage without allocating on the scoring path. */
		public void score(final float[] input, final float[] output){score(input, null, output);}

		/** Also copies pre-interpretation outputs when rawOutput is nonnull.
		 * Caller arrays must be distinct; inputs and prior output snapshots remain private
		 * to the caller. Create/use a worker in one thread; share templates across threads. */
		public void score(final float[] input, final float[] rawOutput, final float[] output){
			require(owner==Thread.currentThread(), "Inference worker belongs to another thread");
			require(input!=null && input.length==template.inputs() && output!=null && output.length==template.outputs(),
				"Input/output widths differ from the frozen model");
			require(input!=output && (rawOutput==null ||
				(rawOutput.length==output.length && rawOutput!=input && rawOutput!=output)), "Caller buffers must be distinct");
			System.arraycopy(input, 0, values[0], 0, input.length);
			if(template.mean!=null){
				for(int i=0; i<input.length; i++){
					values[0][i]-=template.mean[i];
					values[0][i]*=template.inverseStd[i];
				}
				if(template.clip>0){
					for(int i=0; i<input.length; i++){
						final float value=values[0][i];
						// CellNet.normalizeInput rejects overflow BEFORE applying the bound.
						if(!Float.isFinite(value)){throw new IllegalArgumentException("Nonfinite normalized input before clipping at "+i);}
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
					// Vector.feedForwardDense first rounds the dot, then float-adds the bias.
					float sum=template.biases[layer][row]; sum+=dot;
					next[row]=(float)template.functions[layer][row].activate(sum);
				}
			}
			final float[] raw=values[values.length-1];
			if(rawOutput!=null){System.arraycopy(raw, 0, rawOutput, 0, raw.length);}
			if(reported==null){System.arraycopy(raw, 0, output, 0, raw.length);}
			else{
				System.arraycopy(raw, 0, reported, 0, raw.length);
				if(template.prediction!=null){template.prediction.apply(reported);}
				if(template.caps!=null){template.caps.apply(reported);}
				if(template.prediction!=null){reported[2]*=100f; reported[3]*=100f;}
				System.arraycopy(reported, 0, output, 0, raw.length);
			}
		}

		private final MagQCForwardTemplate template;
		private final Thread owner;
		private final float[][] values;
		private final float[] reported;
	}

	/** Rejects invalid snapshots/API use even with assertions disabled. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private final int[] dims;
	private final float[][][] weights;
	private final float[][] biases;
	private final Function[][] functions;
	private final float[] mean, inverseStd;
	private final float clip;
	private final MagQCPredictionContract prediction;
	private final MagQCSubnetOutputCaps caps;
	private final long weightCount;
}
