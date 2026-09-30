package assemble;

import parse.LineParser1;
import structures.ByteBuilder;

/** Exact candidate-order and fixed-decimal model-input parity tests. */
public final class LocalEditNeuralFeaturesTest {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private LocalEditNeuralFeaturesTest(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final int[] depthCases={0, 1, 2, 3, 15, 100, 1000000};
		final int[] depths=new int[LocalEditNeuralFeatures.CANDIDATE_COUNT];
		final double[] features=new double[LocalEditNeuralFeatures.FEATURE_COUNT];
		final float[] input=new float[LocalEditNeuralFeatures.FEATURE_COUNT];
		final ByteBuilder bb=new ByteBuilder(512);
		final LineParser1 parser=new LineParser1('\t');
		long vectors=0;
		for(int original=0; original<4; original++){
			for(int candidate=0; candidate<LocalEditNeuralFeatures.CANDIDATE_COUNT; candidate++){
				final int operation=LocalEditNeuralFeatures.operation(candidate);
				final int base=LocalEditNeuralFeatures.candidateBase(candidate, original);
				if(LocalEditNeuralFeatures.candidateIndex(original, operation, base)!=candidate){
					throw new AssertionError("Candidate order is not invertible.");
				}
				for(final int current:depthCases){
					for(int i=0; i<depths.length; i++){depths[i]=(i+1)*(current+2);}
					final int run=operation==LocalEditNeuralFeatures.DELETION ? 0 : candidate+1;
					if(!LocalEditNeuralFeatures.fill(features, current, depths, candidate, original, current<2 ? -1 : current+7,
							current+11, run, (candidate&1)!=0, 9973, 137+candidate)){
						throw new AssertionError("Complete feature evidence was rejected.");
					}
					if(!LocalEditNeuralFeatures.toModelInput(features, input)){
						throw new AssertionError("Finite feature vector was rejected.");
					}
					bb.clear(); LocalEditNeuralFeatures.append(bb, features, 1); parser.set(bb.toBytes());
					if(parser.terms()!=LocalEditNeuralFeatures.FEATURE_COUNT+1){throw new AssertionError("Wrong serialized dimensions.");}
					for(int i=0; i<LocalEditNeuralFeatures.FEATURE_COUNT; i++){
						final float parsed=parser.parseFloat(i);
						if(Float.floatToIntBits(parsed)!=Float.floatToIntBits(input[i])){
							throw new AssertionError("Model-input mismatch at feature "+i+": "+parsed+" != "+input[i]);
						}
					}
					vectors++;
				}
			}
		}
		if(LocalEditNeuralFeatures.fill(features, 1, depths, 0, 0, -1, -1, 1, false, 100, 50)){
			throw new AssertionError("Both missing flanks must abstain.");
		}
		System.out.println("LOCAL_EDIT_NEURAL_FEATURES_TEST_OK vectors="+vectors);
	}
}
