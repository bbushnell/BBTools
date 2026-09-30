package assemble;

import ml.CellNet;
import ml.CellNetParser;
import ukmer.Kmer;

/** End-to-end frozen-model score parity over a literal complete count surface. */
public final class LocalEditNeuralGateTest {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private LocalEditNeuralGateTest(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		if(args.length!=3){throw new IllegalArgumentException("Expected model path, SHA80, and cutoff.");}
		final float cutoff=Float.parseFloat(args[2]);
		final CellNet gateNet=LocalEditNeuralGate.loadFrozen(args[0], args[1], cutoff);
		final HomopolymerIndelProposal.CountLookup lookup=new HomopolymerIndelProposal.CountLookup(){
			@Override public int count(final Kmer key){return 5;}
		};
		final LocalEditNeuralGate gate=new LocalEditNeuralGate(62, lookup, gateNet, cutoff);
		final byte[] bases=new byte[140];
		final byte[] alphabet={'A', 'C', 'G', 'T'};
		for(int i=0; i<bases.length; i++){bases[i]=alphabet[i&3];}
		final int position=70, original=LocalEditNeuralFeatures.baseIndex(bases[position]);
		final int candidate=LocalEditNeuralFeatures.candidateIndex(original, LocalEditNeuralFeatures.SUBSTITUTION, 0);
		final float observed=gate.score(bases, position, LocalSingleBaseEdit.Operation.SUBSTITUTION, (byte)'A');
		if(!Float.isFinite(observed) || gate.probeQueries!=8 || gate.windowQueries!=3 || gate.evaluations!=1){
			throw new AssertionError("Complete 3+1+4 evidence must produce one finite score with eleven lookups.");
		}
		final int[] depths={5, 5, 5, 5, 5, 5, 5, 5};
		final double[] features=new double[LocalEditNeuralFeatures.FEATURE_COUNT];
		final float[] input=new float[LocalEditNeuralFeatures.FEATURE_COUNT];
		if(!LocalEditNeuralFeatures.fill(features, 5, depths, candidate, original, 5, 5, 1, false, bases.length, position) ||
				!LocalEditNeuralFeatures.toModelInput(features, input)){
			throw new AssertionError("Literal complete feature vector was rejected.");
		}
		final CellNet direct=CellNetParser.load(args[0], false); direct.simdFF=true;
		final float expected=direct.applyInput(input).feedForwardDense();
		if(Float.floatToIntBits(observed)!=Float.floatToIntBits(expected)){
			throw new AssertionError("Runtime gate score differs from the frozen vector path: "+observed+" != "+expected);
		}
		bases[39]='N';
		if(!Float.isNaN(gate.score(bases, position, LocalSingleBaseEdit.Operation.SUBSTITUTION, (byte)'A'))){
			throw new AssertionError("Undefined central evidence must abstain.");
		}
		System.out.println("LOCAL_EDIT_NEURAL_GATE_TEST_OK score="+observed);
	}
}
