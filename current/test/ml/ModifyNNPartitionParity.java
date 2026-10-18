package test.ml;

import ml.Cell;
import ml.CellNet;
import ml.CellNetParser;
import ml.ModifyNNPartition;

import java.util.Arrays;

/**
 * Compares caller-supplied same-parent native and GPU exports, requiring identical
 * dimensions, input normalization, canonical partition ownership and stored FP32 parameters.
 * @author Nilou
 */
public final class ModifyNNPartitionParity{

	/** Loads native/GPU children and fails if topology, canonical ownership or parameter bits diverge. */
	public static void main(String[] args){
		if(args.length!=2){throw new IllegalArgumentException("Expected Java child and GPU child .bbnet paths");}
		final CellNet actual=CellNetParser.load(args[0], false), expected=CellNetParser.load(args[1], false);
		check(Arrays.equals(actual.dims, expected.dims), "Network dimensions differ");
		check(Arrays.equals(actual.inputMeanCopy(), expected.inputMeanCopy()), "Input means differ");
		check(Arrays.equals(actual.inputInverseStdCopy(), expected.inputInverseStdCopy()), "Input inverse deviations differ");
		check(ModifyNNPartition.grow(actual, actual.dims, 0).metadata().equals(
			ModifyNNPartition.grow(expected, expected.dims, 0).metadata()), "Partition ownership differs");
		long weights=0, biases=0;
		for(int layer=1; layer<actual.net.length; layer++){
			for(int h=0; h<actual.net[layer].length; h++){
				final Cell a=actual.net[layer][h], b=expected.net[layer][h];
				check(a.function==b.function, "Activation differs at layer="+layer+" row="+h);
				check(Float.floatToRawIntBits(a.bias)==Float.floatToRawIntBits(b.bias), "Bias bits differ at layer="+layer+" row="+h);
				check(Arrays.equals(a.inputs, b.inputs), "Sparse membership differs at layer="+layer+" row="+h);
				check(a.weights.length==b.weights.length, "Weight lengths differ");
				for(int i=0; i<a.weights.length; i++){
					check(Float.floatToRawIntBits(a.weights[i])==Float.floatToRawIntBits(b.weights[i]),
						"Weight bits differ at layer="+layer+" row="+h+" edge="+i+" actual="+a.weights[i]+" expected="+b.weights[i]);
					weights++;
				}
				biases++;
			}
		}
		System.out.println("PARTITION_PARAMETER_PARITY_PASS weights="+weights+" biases="+biases+" exact_fp32_bits=true");
	}

	/** Converts parity mismatches to diagnostic unchecked failures. */
	private static void check(boolean valid, String message){if(!valid){throw new IllegalStateException(message);}}
}
