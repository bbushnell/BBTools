package align2;

import java.io.IOException;

import ml.CellNet;
import ml.CellNetParser;
import stream.Read;

/** Paired inference with fixed weights and mutable worker-local activations.
 * Each worker owns its network, feature vectors, and extraction/transform scratch.
 * Immutable calibration and caps are shared by copies. Loose and strict targets
 * use separate model/calibration/cap bundles selected by BBMapS.
 * @author Collei
 */
public final class NeuralMapqPairedInference{

	/** Loads a matching bundle from paths already resolved by the caller. */
	public static NeuralMapqPairedInference load(final String netPath, final String lutPath,
			final String capsPath, final int regime)throws IOException{
		final CellNet net=CellNetParser.load(netPath, false);
		if(net==null || net.numInputs()!=NeuralMapqPairedFeatureTransformV1.WIDTH || net.numOutputs()!=1){
			throw new IllegalArgumentException("Paired neural MAPQ network differs from "+NeuralMapqPairedFeatureTransformV1.SCHEMA_NAME);
		}
		return new NeuralMapqPairedInference(net, NeuralMapqCalibrationLut.load(lutPath, 50), NeuralMapqLengthCaps.load(capsPath), regime);
	}

	private NeuralMapqPairedInference(final CellNet net_, final NeuralMapqCalibrationLut lut_,
			final NeuralMapqLengthCaps caps_, final int regime_){
		net=net_; lut=lut_; caps=caps_; regime=regime_;
		raw=new float[NeuralMapqPairedFeatureSchema.WIDTH];
		rawMate=new float[NeuralMapqPairedFeatureSchema.WIDTH];
		transformed=new float[NeuralMapqPairedFeatureTransformV1.WIDTH];
		rawScratch=new NeuralMapqPairedFeatureExtractor.Scratch();
		transformScratch=new NeuralMapqPairedFeatureTransformV1.Scratch();
	}

	/** Allocates independent network activations and scratch for another worker. */
	public NeuralMapqPairedInference copy(){return new NeuralMapqPairedInference(net.copy(false), lut, caps, regime);}

	/** Scores one mapped primary anchor; an unsupported length/regime returns -1
	 * before extraction. Eligible calls require reciprocal mates and the feature
	 * extractor's read/site/traceback state. Does not write either read's MAPQ. */
	public int mapq(final Read anchor, final Read mate, final int averagePairDistance,
			final boolean requireCorrectStrands, final boolean sameStrandPairs){
		final int cap=caps.cap(true, regime, anchor.length());
		if(cap<0){return -1;}
		NeuralMapqPairedFeatureExtractor.fillRuntime(anchor, mate, averagePairDistance,
				requireCorrectStrands, sameStrandPairs, raw, rawScratch);
		NeuralMapqPairedFeatureTransformV1.transform(raw, transformed, transformScratch);
		final float error=net.applyInput(transformed).feedForward();
		return lut.mapq(error, cap);
	}

	/** Writes per-end Q or -1 for unmapped/unsupported ends into the first two slots.
	 * Skips extraction if neither end is eligible. Otherwise both mapped mates'
	 * features are needed, including an unsupported-length mate of an eligible end.
	 * Extracts each raw block once and reuses transform/network scratch serially. */
	public void mapqs(final Read first, final Read second, final int averagePairDistance,
			final boolean requireCorrectStrands, final boolean sameStrandPairs, final int[] out){
		mapqs(first, second, averagePairDistance, requireCorrectStrands, sameStrandPairs, out, null, null);
	}

	/** Scores this model and an optional second model from one preparation per pair.
	 * Both models use the frozen V1 input schema, but retain independent caps,
	 * reference regimes, networks and calibration. All mutable state belongs to
	 * this worker. Neither model nor either output may be used concurrently. */
	void mapqs(final Read first, final Read second, final int averagePairDistance,
			final boolean requireCorrectStrands, final boolean sameStrandPairs, final int[] out,
			final NeuralMapqPairedInference other, final int[] otherOut){
		if(out==null || out.length<2){throw new IllegalArgumentException("Paired neural MAPQ output requires two slots");}
		if(other!=null && (otherOut==null || otherOut.length<2 || otherOut==out)){
			throw new IllegalArgumentException("Two MAPQ models need distinct two-slot outputs to preserve both predictions");
		}
		out[0]=out[1]=-1;
		if(other!=null){otherOut[0]=otherOut[1]=-1;}
		if(!first.mapped() && !second.mapped()){return;}
		final int firstCap=first.mapped() ? caps.cap(true, regime, first.length()) : -1;
		final int secondCap=second.mapped() ? caps.cap(true, regime, second.length()) : -1;
		final int otherFirstCap=other!=null && first.mapped() ? other.caps.cap(true, other.regime, first.length()) : -1;
		final int otherSecondCap=other!=null && second.mapped() ? other.caps.cap(true, other.regime, second.length()) : -1;
		if(firstCap<0 && secondCap<0 && otherFirstCap<0 && otherSecondCap<0){return;}
		NeuralMapqPairedFeatureExtractor.fillBothRuntime(first, second, averagePairDistance,
				requireCorrectStrands, sameStrandPairs, raw, rawMate, rawScratch);
		if(firstCap>=0 || otherFirstCap>=0){
			NeuralMapqPairedFeatureTransformV1.transform(raw, transformed, transformScratch);
			if(firstCap>=0){out[0]=lut.mapq(net.applyInput(transformed).feedForward(), firstCap);}
			if(otherFirstCap>=0){otherOut[0]=other.lut.mapq(other.net.applyInput(transformed).feedForward(), otherFirstCap);}
		}
		if(secondCap>=0 || otherSecondCap>=0){
			NeuralMapqPairedFeatureTransformV1.transform(rawMate, transformed, transformScratch);
			if(secondCap>=0){out[1]=lut.mapq(net.applyInput(transformed).feedForward(), secondCap);}
			if(otherSecondCap>=0){otherOut[1]=other.lut.mapq(other.net.applyInput(transformed).feedForward(), otherSecondCap);}
		}
	}

	private final CellNet net;
	private final NeuralMapqCalibrationLut lut;
	private final NeuralMapqLengthCaps caps;
	private final int regime;
	private final float[] raw, rawMate, transformed;
	private final NeuralMapqPairedFeatureExtractor.Scratch rawScratch;
	private final NeuralMapqPairedFeatureTransformV1.Scratch transformScratch;
}
