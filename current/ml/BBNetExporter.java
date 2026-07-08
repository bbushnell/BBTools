package ml;

import java.util.ArrayList;

/**
 * Builds a BBNet-format CellNet from plain weight arrays (one hidden tanh
 * layer, configurable output activation), folding input standardization
 * ((x-mean)/sd) into layer 1 so the exported net takes raw inputs.  Lives in
 * the ml package because CellNet/Cell internals are package-private.
 *
 * Used by cardinality.BigNetToBBNet to serialize ComplexityFarm.BigNet
 * window-complexity nets into standard BBTools ml format.
 *
 * @author Amber
 * @date July 2026
 */
public class BBNetExporter {

	/**
	 * @param w1 hidden weights [H][D]
	 * @param b1 hidden biases [H]
	 * @param w2 output weights [H]
	 * @param b2 output bias
	 * @param mean per-feature means [D] (folded into layer 1)
	 * @param sd per-feature sds [D] (folded into layer 1)
	 * @param finalType Function type constant for the output cell (e.g.
	 *        Function.RSLOG); hidden cells are TANH.
	 */
	public static CellNet export(double[][] w1, double[] b1, double[] w2,
			double b2, double[] mean, double[] sd, int finalType){
		final int H=w1.length, D=w1[0].length;
		assert(b1.length==H && w2.length==H) : "b1="+b1.length+" w2="+w2.length+" H="+H;
		assert(mean.length==D && sd.length==D) : "mean="+mean.length+" sd="+sd.length+" D="+D;
		final int[] dims={D, H, 1};
		final CellNet net=new CellNet(dims, 1, 1.0f, 0f, 1, new ArrayList<>());
		// randomize() consults activation-rate tables that ml.Trainer normally
		// initializes; every cell's function is assigned explicitly below, so
		// zero the rates (skips random selection, satisfies the assertion).
		java.util.Arrays.fill(Function.TYPE_RATES, 0f);
		net.randomize();   // allocates dense edge topology + weight arrays
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
