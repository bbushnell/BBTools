package ml;

import java.nio.charset.StandardCharsets;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import parse.PreParser;
import prot.MagQCNetBundle;

/**
 * Fixtures for versioned sparse I-line edge encodings.
 * @author Brian, Yelan
 */
public final class CellNetSparseEncodingTest {

	public static void main(String[] args) throws Exception{
		if(args.length>0){args=new PreParser(args, CellNetSparseEncodingTest.class, false).args;}
		boolean enabled=false; assert(enabled=true);
		if(!enabled || args.length>1 || (args.length==1 && !args[0].startsWith("bundle="))){
			throw new IllegalArgumentException("Run with -ea and optional bundle=<shipped .bbnets>");
		}
		checkNativeSparseDeltas();
		checkDenseAsSparseDeltas();
		checkMalformedDeltas();
		checkSparseTraining();
		if(args.length==1){checkBundle(args[0].substring(7));}
		System.out.println("CellNetSparseEncodingTest PASS deltaa48 roundtrip, dense-as-sparse zeros and malformed tokens");
	}

	/** Gapped output rows train like zero-masked dense rows; prefix rows retain the old positional arithmetic. */
	private static void checkSparseTraining(){
		final CellNet sparse=trainingNet(false), dense=parse(write(sparse, true, false, false));
		final CellNet prefix=trainingNet(true);
		CellNet.DENSE=false;
		final CellNet legacy=prefix.copy(false);
		final float[][] rows={{0.25f, -0.5f, 0.75f, 0.5f, -0.25f, 0.125f},
			{-0.5f, 0.25f, -0.125f, 0.75f, 0.5f, -0.25f}};
		for(int step=0; step<16; step++){
			final float[] input=rows[step%rows.length];
			final float[] target={step%2==0 ? 0.25f : 0.5f};
			final float a=predict(sparse, input)[0], b=predict(dense, input)[0];
			check(Math.abs(a-b)<=1e-6f, "Sparse/dense predictions diverged before update "+step);
			CellNet.DENSE=false;
			sparse.backPropSparse(target, 0);
			CellNet.DENSE=true;
			dense.backPropDense(target, 0);
			final Cell sc=sparse.finalLayer[0], dc=dense.finalLayer[0];
			for(int edge=0; edge<sc.inputs.length; edge++){
				check(Math.abs(sc.deltas[edge]-dc.deltas[sc.inputs[edge]])<=1e-6f,
					"Sparse gradient used the wrong source activation at edge "+edge);
			}
			CellNet.DENSE=false; sparse.applyChanges(1, 0.01f);
			CellNet.DENSE=true; dense.applyChanges(1, 0.01f);
			for(int edge=0; edge<sc.inputs.length; edge++){
				check(Math.abs(sc.weights[edge]-dc.weights[sc.inputs[edge]])<=1e-6f,
					"Sparse/dense trained weight mismatch at step "+step);
			}
			check(Math.abs(sc.bias-dc.bias)<=1e-6f, "Sparse/dense trained bias mismatch");
			for(int i=0; i<dc.weights.length; i++){
				if(i!=2 && i!=5){check(dc.weights[i]==0, "Dense equivalent gained an absent edge");}
			}

			checkPredict(prefix, legacy, input);
			prefix.backPropSparse(target, 0);
			final Cell old=legacy.finalLayer[0];
			final int[] indices=old.inputs;
			// The unchanged dense branch is the pre-fix valuesIn[i] arithmetic.
			// Use it only on prefix rows; restore topology before applying updates.
			try{old.inputs=null; old.updateEdgesFinalLayerDense(target[0], legacy.values[0], 0);}
			finally{old.inputs=indices;}
			check(Arrays.equals(prefix.finalLayer[0].deltas, old.deltas), "Prefix gradients changed from positional arithmetic");
			prefix.applyChanges(1, 0.01f); legacy.applyChanges(1, 0.01f);
			check(Arrays.equals(prefix.finalLayer[0].weights, old.weights) &&
				Float.floatToRawIntBits(prefix.finalLayer[0].bias)==Float.floatToRawIntBits(old.bias),
				"Prefix training must remain bit-identical to positional arithmetic");
		}
		checkPredict(sparse, dense, rows[0]);
		System.out.println("SPARSE_OUTPUT_TRAINING_PASS steps=16 indices=2,5 prefix_bit_identical=true");
	}

	/** A small linear output gives distinct source activations and preserves absent dense edges. */
	private static CellNet trainingNet(final boolean prefix){
		CellNet.DENSE=false;
		final CellNet result=new CellNet(new int[]{6, 1}, 7, 1f, 0f, 1, new ArrayList<String>());
		final Cell cell=result.finalLayer[0];
		cell.function=Function.getFunction(Function.LINEAR);
		cell.inputs=prefix ? new int[]{0, 1} : new int[]{2, 5};
		cell.weights=new float[]{0.0625f, -0.03125f};
		cell.deltas=new float[2];
		cell.bias=0.125f;
		CellNet.makeOutputSets(result.net);
		result.makeWeightMatrices();
		return result;
	}

	/** Native sparse rows keep weight order while edge IDs become per-row A48 deltas. */
	private static void checkNativeSparseDeltas(){
		final CellNet net=sparseLinear();
		final String text=write(net, false, true);
		check(text.contains("#indexencoding deltaa48\n"), "delta sparse files must carry a critical header");
		check(text.contains("\nI7 0 1 2\n"), "inputs 0,2,5 must delta-code as 0,1,2");
		check(!text.contains("\nI7 0 2 5\n"), "delta output must not leave literal sparse IDs");
		final CellNet parsed=parse(text);
		check(Arrays.equals(parsed.net[1][0].inputs, new int[]{0, 2, 5}), "delta row0 did not round-trip");
		check(Arrays.equals(parsed.net[1][1].inputs, new int[]{1, 4}), "delta row1 did not round-trip");
		checkPredict(net, parsed, new float[]{1, 2, 3, 4, 5, 6});
	}

	/** Dense-to-sparse output emits a zero delta for every consecutive live input. */
	private static void checkDenseAsSparseDeltas(){
		final CellNet dense=denseLinear();
		final String text=write(dense, true, true);
		check(text.contains("#indexencoding deltaa48\n"), "forced dense-to-sparse delta output must carry a critical header");
		check(text.contains("\nI5 0 0 0 0\n"), "dense-as-sparse IDs 0..3 must delta-code as all zeroes");
		final CellNet parsed=parse(text);
		check(Arrays.equals(parsed.net[1][0].inputs, new int[]{0, 1, 2, 3}), "dense-as-sparse inputs did not round-trip");
		checkPredict(dense, parsed, new float[]{-2, 0.5f, 3, 4});
	}

	/** Delta-A48 is a critical encoding: malformed headers and payloads must crash before misparsing. */
	private static void checkMalformedDeltas(){
		final String text=write(sparseLinear(), false, true);
		final String header="#indexencoding deltaa48\n";
		reject(text.replace("#indexencoding deltaa48", "#sparseencoding deltaa48"));
		reject(text.replace("#indexencoding deltaa48", "#indexencoding wat"));
		reject(text.replace(header, header+header));
		reject(text.replace("#sparse\n", "#dense\n"));
		reject(text.replace("\nI7 0 1 2\n", "\nI7 / 1 2\n"));
		reject(text.replace("\nI7 0 1 2\n", "\nI7 00 1 2\n"));
		reject(text.replace("\nI7 0 1 2\n", "\nI7 0000000 1 2\n"));
		reject(text.replace("\nI7 0 1 2\n", "\nI7 200000 1 2\n"));
		reject(text.replace("\nI7 0 1 2\n", "\nI7 6 1 2\n"));
	}

	/** Writes a net under selected global sparse output flags and restores previous process state. */
	private static String write(CellNet net, boolean outSparse, boolean outDelta){
		return write(net, false, outSparse, outDelta);
	}

	/** Selects one representation while restoring every process-global writer flag. */
	private static String write(CellNet net, boolean outDense, boolean outSparse, boolean outDelta){
		final boolean oldDense=CellNet.DENSE, oldOutDense=CellNet.OUT_DENSE;
		final boolean oldOutSparse=CellNet.OUT_SPARSE, oldOutHex=CellNet.OUT_HEX;
		final boolean oldOutDelta=CellNet.OUT_DELTA_A48;
		try{
			CellNet.DENSE=net.net[1][0].inputs==null;
			CellNet.OUT_DENSE=outDense;
			CellNet.OUT_SPARSE=outSparse;
			CellNet.OUT_HEX=false;
			CellNet.OUT_DELTA_A48=outDelta;
			return net.toBytes().toString();
		}finally{
			CellNet.DENSE=oldDense;
			CellNet.OUT_DENSE=oldOutDense;
			CellNet.OUT_SPARSE=oldOutSparse;
			CellNet.OUT_HEX=oldOutHex;
			CellNet.OUT_DELTA_A48=oldOutDelta;
		}
	}

	/** Reencodes every shipped subnet in three layouts and checks weights and raw-input predictions. */
	private static void checkBundle(final String path) throws Exception{
		final MagQCNetBundle bundle=MagQCNetBundle.loadMultiOutput(Paths.get(path));
		check(bundle.size()==131, "The shipped v1.2.1 fixture must contain131 subnets");
		long comparisons=0;
		double maximum=0;
		for(int index=0; index<bundle.size(); index++){
			final MagQCNetBundle.Subnet subnet=bundle.subnet(index);
			final CellNet original=parse(new String(subnet.netBytes, StandardCharsets.US_ASCII));
			for(int representation=0; representation<3; representation++){
				final String text=write(original, representation==0, representation!=0, representation==2);
				check(text.contains("#indexencoding deltaa48\n")== (representation==2),
					"Unexpected delta header for "+subnet.id+" representation="+representation);
				final CellNet decoded=parse(text);
				checkParameters(original, decoded, subnet.id);
				for(int row=0; row<3; row++){
					final float[] input=new float[original.dims[0]];
					for(int i=0; i<input.length; i++){input[i]=row==0 ? 0 : row==1 ? 1 : (i%23-11)/11f;}
					final float[] expected=predict(original, input), actual=predict(decoded, input);
					check(expected.length==actual.length, "Output width changed for "+subnet.id);
					for(int i=0; i<expected.length; i++){
						final double difference=Math.abs((double)expected[i]-actual[i]);
						check(Float.isFinite(expected[i]) && Float.isFinite(actual[i]) && difference<=1e-5,
							"Prediction changed for "+subnet.id+" layout="+representation+" row="+row+" output="+i+
							" expected="+expected[i]+" actual="+actual[i]);
						maximum=Math.max(maximum, difference);
						comparisons++;
					}
				}
			}
		}
		System.out.println("SHIPPED_BUNDLE_ENCODING_PASS subnets="+bundle.size()+
			" representations=3 output_comparisons="+comparisons+" max_abs="+maximum);
	}

	/** Checks decoded numerical parameters independently of dense versus sparse storage. */
	private static void checkParameters(final CellNet expected, final CellNet actual, final String id){
		check(Arrays.equals(expected.dims, actual.dims), "Layer dimensions changed for "+id);
		check(Arrays.equals(expected.inputMeanCopy(), actual.inputMeanCopy()) &&
			Arrays.equals(expected.inputInverseStdCopy(), actual.inputInverseStdCopy()), "Normalization changed for "+id);
		for(int layer=1; layer<expected.net.length; layer++){
			for(int sink=0; sink<expected.net[layer].length; sink++){
				final Cell a=expected.net[layer][sink], b=actual.net[layer][sink];
				check(a.function==b.function && Float.floatToRawIntBits(a.bias)==Float.floatToRawIntBits(b.bias),
					"Activation or bias changed for "+id+" layer="+layer+" sink="+sink);
				final float[] aw=expandedWeights(a, expected.dims[layer-1]);
				final float[] bw=expandedWeights(b, actual.dims[layer-1]);
				for(int i=0; i<aw.length; i++){
					check(aw[i]==bw[i], "Weight changed for "+id+" layer="+layer+" sink="+sink+" input="+i);
				}
			}
		}
	}

	/** Expands only a fixture row; absent sparse edges have numerical weight zero. */
	private static float[] expandedWeights(final Cell cell, final int width){
		if(cell.inputs==null){check(cell.weights.length==width, "Dense row width mismatch"); return cell.weights;}
		final float[] result=new float[width];
		check(cell.inputs.length==cell.weights.length, "Sparse index/weight width mismatch");
		for(int i=0; i<cell.inputs.length; i++){result[cell.inputs[i]]=cell.weights[i];}
		return result;
	}

	/** Creates a native sparse fixture with nonconsecutive inputs to exercise nonzero deltas. */
	private static CellNet sparseLinear(){
		return sparseLinear(false);
	}

	/** Builds each topology before its matrix aliases are initialized exactly once. */
	private static CellNet sparseLinear(final boolean prefix){
		CellNet.DENSE=false;
		final CellNet net=new CellNet(new int[]{6, 2}, 2, 1f, 0f, 1, new ArrayList<String>());
		net.net[1][0].function=Function.getFunction(Function.LINEAR);
		net.net[1][0].bias=1;
		net.net[1][0].inputs=prefix ? new int[]{0, 1, 2} : new int[]{0, 2, 5};
		net.net[1][0].weights=new float[]{0.5f, -0.25f, 0.125f};
		net.net[1][0].deltas=new float[3];
		net.net[1][1].function=Function.getFunction(Function.LINEAR);
		net.net[1][1].bias=-2;
		net.net[1][1].inputs=prefix ? new int[]{0, 1} : new int[]{1, 4};
		net.net[1][1].weights=new float[]{3, -4};
		net.net[1][1].deltas=new float[2];
		CellNet.makeOutputSets(net.net);
		net.makeWeightMatrices();
		return net;
	}

	/** Creates a dense fixture with all four inputs live. */
	private static CellNet denseLinear(){
		CellNet.DENSE=true;
		final CellNet net=new CellNet(new int[]{4, 1}, 3, 1f, 0f, 1, new ArrayList<String>());
		final Cell c=net.net[1][0];
		c.function=Function.getFunction(Function.LINEAR);
		c.bias=0.25f;
		c.weights=new float[]{1, -2, 3, -4};
		c.deltas=new float[c.weights.length];
		net.makeWeightMatrices();
		return net;
	}

	/** Checks exact parity on one deterministic input row. */
	private static void checkPredict(CellNet a, CellNet b, float[] input){
		final float[] pa=predict(a, input), pb=predict(b, input);
		check(pa.length==pb.length, "prediction widths differ");
		for(int i=0; i<pa.length; i++){
			check(Float.isFinite(pa[i]) && Float.isFinite(pb[i]), "Nonfinite prediction at output "+i);
			check(Float.floatToRawIntBits(pa[i])==Float.floatToRawIntBits(pb[i]), "prediction mismatch at output "+i);
		}
	}

	/** Runs a forward pass and returns its output. */
	private static float[] predict(CellNet net, float[] input){
		CellNet.DENSE=net.net[1][0].inputs==null;
		net.applyInput(input);
		net.feedForward();
		return net.getOutput();
	}

	/** Converts tiny fixture text to the production byte-line parser input. */
	private static ArrayList<byte[]> lines(String text){
		final ArrayList<byte[]> result=new ArrayList<byte[]>();
		for(String line:text.split("\n")){result.add(line.getBytes(StandardCharsets.US_ASCII));}
		return result;
	}

	private static CellNet parse(String text){return CellNetParser.loadFromLines(lines(text));}

	/** Malformed critical encoding must throw outside assertions-only validation. */
	private static void reject(String text){
		try{parse(text);}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Malformed delta-A48 sparse encoding accepted");
	}

	private static void check(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
