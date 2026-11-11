package test.assemble;

import assemble.LocalEditNeuralFeatures;
import parse.LineParser1;
import structures.ByteBuilder;

/**
 * Checks candidate-index round trips and 28-feature serialization/conversion parity.
 * Uses the production mapping in both directions, so invertibility does not by itself
 * establish the intended biological candidate order. Synthetic depth cases exercise
 * exact float-bit agreement between model-input conversion and serialized values;
 * they do not prove arbitrary-input quantization or model accuracy. No model is loaded.
 * Explicit failures remain active with JVM assertions disabled.
 */
public final class LocalEditNeuralFeaturesTest{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Static fixture entry point only. */
	private LocalEditNeuralFeaturesTest(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Checks four original bases, eight candidates and seven current depths: 224 vectors.
	 * One missing flank is exercised at depths zero/one; both missing flanks must abstain.
	 * Feature values come from production fill, not an independent formula oracle.
	 * Checks serialized dimensions and feature values, but not the appended label value.
	 * Arguments are ignored; the success line follows all checks. The 39-feature layout
	 * and model inference are covered separately by LocalEditReadDepthTest.
	 */
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
