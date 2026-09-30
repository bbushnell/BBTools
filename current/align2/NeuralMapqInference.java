package align2;

import java.io.IOException;

import ml.CellNet;
import ml.CellNetParser;
import stream.Read;

/** Single-end inference with frozen weights and mutable worker-local activations.
 * BBMapS loads separate loose/strict templates and copies each for every worker.
 * Calibration and length caps are immutable and shared; feature arrays, network
 * activations, and extraction scratch are private to each instance.
 * @author Collei
 */
public final class NeuralMapqInference{

	/** Loads a matched model/calibration/cap bundle from already resolved paths.
	 * BBMapS resolves '?' resources before calling this method. Width is checked;
	 * the caller selects the matching schema, loose/strict target, and regime. */
	public static NeuralMapqInference load(final String netPath, final String lutPath,
			final String capsPath, final int regime)throws IOException{
		final CellNet net=CellNetParser.load(netPath, false);
		if(net==null || net.numInputs()!=NeuralMapqFeatureTransform.WIDTH || net.numOutputs()!=1){
			throw new IllegalArgumentException("Neural MAPQ network differs from "+NeuralMapqFeatureTransform.SCHEMA_NAME);
		}
		return new NeuralMapqInference(net, NeuralMapqCalibrationLut.load(lutPath, 50), NeuralMapqLengthCaps.load(capsPath), regime);
	}

	private NeuralMapqInference(final CellNet net_, final NeuralMapqCalibrationLut lut_,
			final NeuralMapqLengthCaps caps_, final int regime_){
		net=net_; lut=lut_; caps=caps_; regime=regime_;
		raw=new float[NeuralMapqFeatureSchema.WIDTH];
		transformed=new float[NeuralMapqFeatureTransform.WIDTH];
		scratch=new NeuralMapqFeatureExtractor.Scratch();
	}

	/** Copies network state and allocates fresh feature/scratch storage for a worker. */
	public NeuralMapqInference copy(){return new NeuralMapqInference(net.copy(false), lut, caps, regime);}

	/** Returns calibrated Q, or -1 for an unsupported length/regime before extraction.
	 * Supported reads must be unpaired mapped primaries with bases, retained sites,
	 * and an existing traceback (FeatureExtractor.fillRuntime's contract).
	 * Mutates this instance's scratch, not the read or its cached MAPQ. */
	public int mapq(final Read read){
		final int cap=caps.cap(false, regime, read.length());
		if(cap<0){return -1;}
		NeuralMapqFeatureExtractor.fillRuntime(read, raw, scratch);
		NeuralMapqFeatureTransform.transform(raw, transformed);
		final float error=net.applyInput(transformed).feedForward();
		return lut.mapq(error, cap);
	}

	private final CellNet net;
	private final NeuralMapqCalibrationLut lut;
	private final NeuralMapqLengthCaps caps;
	private final int regime;
	private final float[] raw, transformed;
	private final NeuralMapqFeatureExtractor.Scratch scratch;
}
