package ml;

import java.util.ArrayList;
import java.util.Arrays;

/** Synthetic bypass, ownership and full-precision serving checks, without trained-model claims.
 * @author Brian Bushnell, Yoimiya */
public final class CellNetGraftTest {

	public static void main(String[] args){
		final CellNet base=linear(new int[]{2, 2, 2, 2});
		final CellNet privateNodes=linear(new int[]{2, 2});
		final CellNet heads=linear(new int[]{3, 2});
		heads.finalLayer[0].weights=new float[]{1, 1, 0};
		heads.finalLayer[1].weights=new float[]{1, 0, 1};
		for(Cell cell:heads.finalLayer){cell.function=Function.getFunction(Function.SIG);}
		final CellNetGraft graft=new CellNetGraft(base, privateNodes, heads, 1, 1, 1, "tanh");
		final float[] input={.25f, .5f};
		base.applyInput(input); base.feedForward();
		final float[] output={base.getOutput(0), base.getOutput(1), 0, 0};
		final float[] original=output.clone();
		graft.appendAfterBase(output, 2);
		check(bits(output[0])==bits(original[0]) && bits(output[1])==bits(original[1]), "Graft must preserve copied base outputs");
		check(bits(base.getOutput(0))==bits(original[0]) && bits(base.getOutput(1))==bits(original[1]), "Graft must not mutate base activations");
		near(output[2], 1/(1+Math.exp(-(.25+Math.tanh(.25)))));
		near(output[3], 1/(1+Math.exp(-(.25+Math.tanh(.5)))));
		check(Arrays.equals(input, new float[]{.25f, .5f}), "Inference must not mutate caller input");
		privateNodes.setInputNormalization(new float[]{0, 0}, new float[]{1, 1});
		heads.setInputNormalization(new float[]{0, 0, 0}, new float[]{1, 1, 1});
		final CellNetGraft explicit=new CellNetGraft(base, privateNodes, heads, 1, 1, 1, "tanh");
		final float[] identityOutput=output.clone(); explicit.appendAfterBase(identityOutput, 2);
		check(Arrays.equals(output, identityOutput), "Exported explicit identity normalization must preserve graft values");
		privateNodes.setInputNormalization(new float[]{.1f, 0}, new float[]{1, 1});
		boolean nonidentity=false;
		try{new CellNetGraft(base, privateNodes, heads, 1, 1, 1, "tanh");}
		catch(IllegalArgumentException expected){nonidentity=true;}
		check(nonidentity, "A component must not silently normalize the frozen activation a second time");
		privateNodes.setInputNormalization(null, null);
		heads.finalLayer[0].weights[2]=1;
		boolean failed=false;
		try{new CellNetGraft(base, privateNodes, heads, 1, 1, 1, "tanh");}
		catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Connections to another head's private nodes must be rejected");
		heads.finalLayer[0].weights[2]=0;
		final float precise=.25003f;
		base.applyInput(new float[]{precise, precise}); base.feedForward();
		final float[] full={base.getOutput(0), base.getOutput(1), 0, 0};
		graft.appendAfterBase(full, 2);
		near(full[2], 1/(1+Math.exp(-(precise+Math.tanh(precise)))));
		near(full[3], 1/(1+Math.exp(-(precise+Math.tanh(precise)))));
		check(bits(base.getOutput(0))==bits(precise), "The full-precision base output must remain untouched");
		System.out.println("CellNetGraftTest PASS: frozen base, shared bypass, private ownership, full-precision activations");
	}

	/** Explicit dense test workers avoid changing CellNet's global dispatch switch. */
	private static CellNet linear(final int[] dims){
		final CellNet net=new CellNet(dims, 1, 1, 1, 1, new ArrayList<String>()){
			@Override public float feedForward(){return feedForwardDense();}
		};
		for(int layer=1; layer<dims.length; layer++){
			for(int row=0; row<dims[layer]; row++){
				final Cell cell=net.net[layer][row];
				cell.function=Function.getFunction(Function.LINEAR);
				cell.weights=new float[dims[layer-1]]; cell.weights[row%dims[layer-1]]=1;
				cell.inputs=null; cell.bias=0;
			}
		}
		return net;
	}
	private static int bits(final float value){return Float.floatToRawIntBits(value);}
	private static void near(final double actual, final double expected){check(Math.abs(actual-expected)<1e-7, "Unexpected graft value: "+actual+" vs "+expected);}
	private static void check(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
}
