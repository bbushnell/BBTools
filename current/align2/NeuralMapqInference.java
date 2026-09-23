package align2;

import java.io.IOException;

import ml.CellNet;
import ml.CellNetParser;
import stream.Read;

/** One mapper-thread's frozen neural-MAPQ inference state. */
public final class NeuralMapqInference {

	public static NeuralMapqInference load(final String netPath,final String lutPath,final String capsPath,final int regime)throws IOException{
		final CellNet net=CellNetParser.load(netPath,false);
		if(net==null||net.numInputs()!=NeuralMapqFeatureTransform.WIDTH||net.numOutputs()!=1){
			throw new IllegalArgumentException("Neural MAPQ network differs from "+NeuralMapqFeatureTransform.SCHEMA_NAME);
		}
		return new NeuralMapqInference(net,NeuralMapqCalibrationLut.load(lutPath,50),NeuralMapqLengthCaps.load(capsPath),regime);
	}

	private NeuralMapqInference(final CellNet net_,final NeuralMapqCalibrationLut lut_,final NeuralMapqLengthCaps caps_,final int regime_){
		net=net_;lut=lut_;caps=caps_;regime=regime_;raw=new float[NeuralMapqFeatureSchema.WIDTH];transformed=new float[NeuralMapqFeatureTransform.WIDTH];scratch=new NeuralMapqFeatureExtractor.Scratch();
	}

	public NeuralMapqInference copy(){return new NeuralMapqInference(net.copy(false),lut,caps,regime);}

	public int mapq(final Read read){
		final int cap=caps.cap(false,regime,read.length());if(cap<0){return -1;}
		NeuralMapqFeatureExtractor.fillRuntime(read,raw,scratch);
		NeuralMapqFeatureTransform.transform(raw,transformed);
		final float error=net.applyInput(transformed).feedForward();
		return lut.mapq(error,cap);
	}

	private final CellNet net;private final NeuralMapqCalibrationLut lut;private final NeuralMapqLengthCaps caps;private final int regime;
	private final float[] raw,transformed;private final NeuralMapqFeatureExtractor.Scratch scratch;
}
