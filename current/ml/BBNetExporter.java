package ml;

import java.util.ArrayList;

/**
 * Builds a BBNet-format {@link CellNet} from one hidden layer of plain weight
 * arrays, replacing the allocated random weights and functions with
 * deterministic caller-supplied parameters.  Input standardization
 * {@code (x-mean)/sd} is folded into the hidden weights and biases, so the
 * resulting net consumes raw inputs.
 *
 * @author Amber
 * @date July 2026
 */
public class BBNetExporter{

	/**
	 * Exports a dense {@code D-H-1} network with TANH hidden cells and the
	 * requested output activation type.
	 *
	 * @param w1 rectangular hidden weight matrix, [H][D]
	 * @param b1 hidden biases, [H]
	 * @param w2 output weights, [H]
	 * @param b2 output bias
	 * @param mean per-feature means to fold into the hidden layer, [D]
	 * @param sd per-feature standard deviations to fold into the hidden layer, [D]
	 * @param finalType {@link Function} type constant for the output cell
	 * @return a fresh dense network whose first layer accepts unnormalized inputs
	 */
	public static CellNet export(double[][] w1, double[] b1, double[] w2,
			double b2, double[] mean, double[] sd, int finalType){
		final int H=w1.length, D=w1[0].length;
		assert(b1.length==H && w2.length==H) : "b1="+b1.length+" w2="+w2.length+" H="+H;
		assert(mean.length==D && sd.length==D) : "mean="+mean.length+" sd="+sd.length+" D="+D;
		final int[] dims={D, H, 1};
		final CellNet net=new CellNet(dims, 1, 1.0f, 0f, 1, new ArrayList<>());
		net.randomizeEdgesAndWeights(); //allocates dense edge topology + weight arrays
		final Cell[] hid=net.net[1];
		for(int h=0; h<H; h++){
			final Cell c=hid[h];
			c.function=Function.getFunction(Function.TANH);
			double b=b1[h];
			for(int d=0; d<D; d++){
				c.weights[d]=(float)(w1[h][d]/sd[d]);
				b-=w1[h][d]*mean[d]/sd[d];
			}
			c.setBias((float)b, true);
		}
		final Cell out=net.net[2][0];
		out.function=Function.getFunction(finalType);
		for(int h=0; h<H; h++){out.weights[h]=(float)w2[h];}
		out.setBias((float)b2, true);
		return net;
	}
}
