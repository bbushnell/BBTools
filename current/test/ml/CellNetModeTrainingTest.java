package test.ml;

import ml.Cell;
import ml.CellNet;
import ml.CellNetParser;
import ml.Sample;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;

/**
 * Selectable one-step -ea training fixture after dense/sparse parser loads in
 * the same JVM.
 * @author Brian Bushnell, Sucrose
 */
public final class CellNetModeTrainingTest{

	/** Runs one named ordering, or all four in sequence; a failure stops the invocation. */
	public static void main(String[] args){
		if(args.length==1){run(args[0]); return;}
		assert(args.length==0) : "This analytic fixture accepts only one optional case name";
		run("homogeneous_dense");
		run("homogeneous_sparse");
		run("dense_after_sparse");
		run("sparse_after_dense");
	}

	/** Selects the first network and the later parser load that should not affect it. */
	private static void run(String name){
		final boolean firstDense, secondDense;
		if(name.equals("homogeneous_dense")){
			firstDense=true;
			secondDense=true;
		}else if(name.equals("homogeneous_sparse")){
			firstDense=false;
			secondDense=false;
		}else if(name.equals("dense_after_sparse")){
			firstDense=true;
			secondDense=false;
		}else if(name.equals("sparse_after_dense")){
			firstDense=false;
			secondDense=true;
		}else{
			throw new IllegalArgumentException("Unknown fixture mode: "+name);
		}
		checkLoadOrder(firstDense, secondDense);
		System.out.println("CellNetModeTrainingTest PASS "+name);
	}

	/** Trains the first net after another load, which changed global mode before the per-network repair. */
	private static void checkLoadOrder(boolean firstDense, boolean secondDense){
		final CellNet first=parse(fixture(firstDense));
		final CellNet second=parse(fixture(secondDense));
		assert(first!=null && second!=null) : "Literal parser fixtures must construct both networks";
		checkTraining(first, firstDense);
	}

	/** Checks forward, gradients, one update, and delta clearing for the analytic row. */
	private static void checkTraining(CellNet net, boolean dense){
		final Sample sample=new Sample(new float[]{2f, 99f, -0.5f}, new float[]{8f}, 1f, 1);
		net.processSample(sample, true, 0f);
		checkClose(sample.result[0], 1f, "forward sample result");
		checkClose(net.getOutput(0), 1f, "forward network result");

		final Cell cell=net.net[1][0];
		checkClose(cell.bias, 7f, "pre-apply bias");
		checkClose(cell.weights[0], -2f, "pre-apply weight0");
		checkClose(cell.deltas[0], 14f, "pre-apply live delta0");
		if(dense){
			checkClose(cell.weights[1], 0f, "pre-apply absent dense weight1");
			checkClose(cell.deltas[1], 693f, "pre-apply absent dense delta1");
			checkClose(cell.weights[2], 4f, "pre-apply weight2");
			checkClose(cell.deltas[2], -3.5f, "pre-apply live delta2");
		}else{
			checkClose(cell.weights[1], 4f, "pre-apply weight2");
			checkClose(cell.deltas[1], -3.5f, "pre-apply live delta2");
		}

		net.applyChanges(1, 0.001f);
		checkClose(cell.bias, 7.007f, "post-apply bias");
		checkClose(cell.weights[0], -1.986f, "post-apply weight0");
		checkClose(cell.deltas[0], 0f, "post-apply delta0");
		if(dense){
			checkClose(cell.weights[1], 0f, "post-apply absent dense weight1");
			checkClose(cell.deltas[1], 0f, "post-apply absent dense delta1");
			checkClose(cell.weights[2], 3.9965f, "post-apply weight2");
			checkClose(cell.deltas[2], 0f, "post-apply delta2");
		}else{
			checkClose(cell.weights[1], 3.9965f, "post-apply weight2");
			checkClose(cell.deltas[1], 0f, "post-apply delta2");
		}
	}

	/** Two live edges decrease in magnitude, so default amplitude damping stays inactive. */
	private static String fixture(boolean dense){
		return "##bbnet\n#version 1\n#concise\n"+(dense ? "#dense\n" : "#sparse\n")+
			"#density 1\n#blocksize 1\n#seed 1\n#layers 2\n#dims 3 1\n#coding decimal\n"+
			(dense ? "#zeroabsent true\n#liveedges 2\nC4 LINEAR 7 -2 0 4\n" :
				"I4 0 2\nW4 LINEAR 7 -2 4\n");
	}

	/** Parses exactly the same text format used by normal CellNet file loading. */
	private static CellNet parse(String text){
		assert(text!=null) : "A parser fixture must supply its serialized network";
		final ArrayList<byte[]> lines=new ArrayList<byte[]>();
		for(String line : text.split("\n")){lines.add(line.getBytes(StandardCharsets.UTF_8));}
		return CellNetParser.loadFromLines(lines);
	}

	/** The 1e-6f tolerance covers float learning-rate rounding but not wrong inputs. */
	private static void checkClose(float observed, float expected, String label){
		if(!Float.isFinite(observed)){
			throw new AssertionError(label+" was nonfinite: "+observed);
		}
		if(Math.abs(observed-expected)>TOLERANCE){
			throw new AssertionError(label+" observed "+observed+" != expected "+expected);
		}
	}

	private static final float TOLERANCE=1e-6f;

}
