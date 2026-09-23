package align2;

import java.io.IOException;

import ml.CellNet;
import ml.CellNetParser;
import stream.Read;

/** One mapper thread's frozen paired neural-MAPQ inference state. */
public final class NeuralMapqPairedInference {

	public static NeuralMapqPairedInference load(final String netPath,final String lutPath,final String capsPath,final int regime)throws IOException{
		final CellNet net=CellNetParser.load(netPath,false);
		if(net==null||net.numInputs()!=NeuralMapqPairedFeatureTransformV1.WIDTH||net.numOutputs()!=1){
			throw new IllegalArgumentException("Paired neural MAPQ network differs from "+NeuralMapqPairedFeatureTransformV1.SCHEMA_NAME);
		}
		return new NeuralMapqPairedInference(net,NeuralMapqCalibrationLut.load(lutPath,50),NeuralMapqLengthCaps.load(capsPath),regime);
	}

	private NeuralMapqPairedInference(final CellNet net_,final NeuralMapqCalibrationLut lut_,final NeuralMapqLengthCaps caps_,final int regime_){
		net=net_;lut=lut_;caps=caps_;regime=regime_;
		raw=new float[NeuralMapqPairedFeatureSchema.WIDTH];
		rawMate=new float[NeuralMapqPairedFeatureSchema.WIDTH];
		transformed=new float[NeuralMapqPairedFeatureTransformV1.WIDTH];
		rawScratch=new NeuralMapqPairedFeatureExtractor.Scratch();
		transformScratch=new NeuralMapqPairedFeatureTransformV1.Scratch();
	}

	public NeuralMapqPairedInference copy(){return new NeuralMapqPairedInference(net.copy(false),lut,caps,regime);}

	public int mapq(final Read anchor,final Read mate,final int averagePairDistance,
			final boolean requireCorrectStrands,final boolean sameStrandPairs){
		final int cap=caps.cap(true,regime,anchor.length());if(cap<0){return -1;}
		NeuralMapqPairedFeatureExtractor.fill(anchor,mate,averagePairDistance,
				requireCorrectStrands,sameStrandPairs,raw,rawScratch);
		NeuralMapqPairedFeatureTransformV1.transform(raw,transformed,transformScratch);
		final float error=net.applyInput(transformed).feedForward();
		return lut.mapq(error,cap);
	}

	/** Scores both orientations after extracting each mapped mate's raw block once. */
	public void mapqs(final Read first,final Read second,final int averagePairDistance,
			final boolean requireCorrectStrands,final boolean sameStrandPairs,final int[] out){
		if(out==null||out.length<2){throw new IllegalArgumentException("Paired neural MAPQ output requires two slots");}
		out[0]=out[1]=-1;
		if(!first.mapped()&&!second.mapped()){return;}
		NeuralMapqPairedFeatureExtractor.fillBoth(first,second,averagePairDistance,
				requireCorrectStrands,sameStrandPairs,raw,rawMate,rawScratch);
		if(first.mapped()){
			final int cap=caps.cap(true,regime,first.length());
			if(cap>=0){NeuralMapqPairedFeatureTransformV1.transform(raw,transformed,transformScratch);out[0]=lut.mapq(net.applyInput(transformed).feedForward(),cap);}
		}
		if(second.mapped()){
			final int cap=caps.cap(true,regime,second.length());
			if(cap>=0){NeuralMapqPairedFeatureTransformV1.transform(rawMate,transformed,transformScratch);out[1]=lut.mapq(net.applyInput(transformed).feedForward(),cap);}
		}
	}

	private final CellNet net;private final NeuralMapqCalibrationLut lut;private final NeuralMapqLengthCaps caps;private final int regime;
	private final float[] raw,rawMate,transformed;
	private final NeuralMapqPairedFeatureExtractor.Scratch rawScratch;
	private final NeuralMapqPairedFeatureTransformV1.Scratch transformScratch;
}
